#!/usr/bin/env python
from __future__ import division, print_function

import magic
import os
import subprocess
import tempfile
from glob import glob

import dxpy


@dxpy.entry_point('main')
def main(tumor_bams=None, normal_bams=None, cn_reference=None,
         is_male_normal=True, baits=None, fasta=None, access=None,
         annotation=None, do_parallel=True):

    if cn_reference and any((normal_bams, baits, fasta, access, annotation)):
        raise dxpy.AppError("Reference profile (cn_reference) cannot be used "
                            "alongside normal_bams, baits, fasta, access "
                            "or annotation")
    if tumor_bams and not any((baits, cn_reference)):
        raise dxpy.AppError("Need cn_reference or baits to process tumor_bams")

    # Set R library path (pre-build packages are bundled here)
    # eg. "/home/dnanexus/R/x86_64-pc-linux-gnu-library/3.0/"
    os.environ["R_LIBS_USER"] = "/".join([os.environ["HOME"],
                                          "R/x86_64-pc-linux-gnu-library",
                                          "3.0"])
    # Install the Python dependencies, then the package itself
    sh("pip install -v --no-index --find-links=file:///wheelhouse -r /requirements.txt")
    sh("pip install -v --no-index --find-links=file:///wheelhouse --no-deps cnvkit")
    sh("chmod +x /usr/local/bin/cnvkit.py")

    print("Downloading file inputs to the local file system")
    cn_reference = download_link(cn_reference)
    baits = download_link(baits)
    fasta = download_link(fasta)
    access = download_link(access)
    annotation = download_link(annotation)
    if tumor_bams is not None:
        tumor_bams = map(download_link, tumor_bams)
    if normal_bams is not None:
        normal_bams = map(download_link, normal_bams)

    # If these input files are gzipped, decompress them
    fasta = maybe_gunzip(fasta, "ref", "fa")
    annotation = maybe_gunzip(annotation, "annot", "txt")

    out_fnames = run_cnvkit(tumor_bams, normal_bams, cn_reference,
                            is_male_normal, baits, fasta, access, annotation,
                            do_parallel)

    print("Uploading local file outputs to the DNAnexus platform")
    output = {}
    for filekey in ("cn_reference", "metrics", "genders", "scatter_pdf",
                    "diagram_pdf"):
        if filekey in out_fnames:
            output[filekey] = dxpy.dxlink(
                dxpy.upload_local_file(out_fnames[filekey]))
    for listkey in ("copy_ratios", "copy_segments", "gainloss", "breaks",
                    "nexus_cn"):
        if listkey in out_fnames:
            output[listkey] = [dxpy.dxlink(dxpy.upload_local_file(fname))
                               for fname in out_fnames[listkey]]
    return output


def run_cnvkit(tumor_bams, normal_bams, reference, is_male_normal, baits, fasta,
               access, annotation, do_parallel):
    """Run the CNVkit pipeline.

    Returns a dict of the generated file names.
    """
    yflag = "-y" if is_male_normal else ""

    print("Running the main CNVkit pipeline")
    command = ["cnvkit.py batch", yflag]
    if do_parallel:
        command.append("-p 0")
    if tumor_bams:
        command.extend(tumor_bams)
        command.extend(["--scatter", "--diagram"])
    if reference:
        command.extend(["-r", reference])
        print("Determining if the given reference profile is male or female")
        if shout("cnvkit.py gender", reference, "| cut -f 2"
                ).strip() == "Male":
            print("Looks like a male reference")
            yflag = "-y"
    else:
        reference = safe_fname("cnv-reference", "cnn")
        command.extend(["--output-reference", reference, "-n"])
        if normal_bams:
            command.extend(normal_bams)
        if baits:
            command.extend(["-t", baits, "--split", "--short-names"])
        if fasta:
            command.extend(["-f", fasta])
            if not access:
                access = safe_fname("access-10k", "bed")
                sh("cnvkit.py access -s 10000", fasta,
                   "-o", access)
        if access:
            command.extend(["-g", access])
        if annotation:
            command.extend(["--annotate", annotation])
    command.append(yflag)

    sh(*command)
    sh("ls -Altr")  # Show the generated files in the DNAnexus log

    # Collect the outputs
    if not tumor_bams:
        return {"cn_reference": reference}

    all_cnr = glob("*.cnr")
    all_cns = glob("*.cns")
    all_gainloss = []
    all_breaks = []
    all_nexus = []
    for acnr, acns in zip(all_cnr, all_cns):
        name = acnr.split('.')[0]

        gainloss = name + "-gainloss.csv"
        sh("cnvkit.py gainloss", acnr, "-s", acns, "-m 3 -t 0.3", yflag,
           "-o", gainloss)
        all_gainloss.append(gainloss)

        breaks = name + "-breaks.csv"
        sh("cnvkit.py breaks", acnr, acns, "-m 3", "-o", breaks)
        all_breaks.append(breaks)

        nexus = name + ".nexus"
        sh("cnvkit.py export nexus-basic", name + ".cnr", "-o", nexus)
        all_nexus.append(nexus)

    genders = safe_fname("gender", "csv")
    sh("cnvkit.py gender", yflag, "-o", genders, *all_cnr)

    metrics = safe_fname("metrics", "csv" if len(all_cnr) > 1 else "txt")
    sh("cnvkit.py metrics", " ".join(all_cnr), "-s", " ".join(all_cns),
       "-o", metrics)

    outputs = {
        "cn_reference": reference,
        "copy_ratios": all_cnr,
        "copy_segments": all_cns,
        "gainloss": all_gainloss,
        "breaks": all_breaks,
        "nexus_cn": all_nexus,
        "metrics": metrics,
        "genders": genders,
    }

    all_scatters = glob("*-scatter.pdf")
    if all_scatters:
        if len(all_scatters) == 1:
            scatter_pdf = all_scatters[0]
        else:
            scatter_pdf = safe_fname("cnv-scatters", "pdf")
            sh("pdfunite", " ".join(all_scatters), scatter_pdf)
        outputs["scatter_pdf"] = scatter_pdf

    all_diagrams = glob("*-diagram.pdf")
    if all_diagrams:
        if len(all_diagrams) == 1:
            diagram_pdf = all_diagrams[0]
        else:
            diagram_pdf = safe_fname("cnv-diagrams", "pdf")
            sh("pdfunite", " ".join(all_diagrams), diagram_pdf)
        outputs["diagram_pdf"] = diagram_pdf

    return outputs


# _____________________________________________________________________________

def download_link(dxlink):
    """Download a DNAnexus object reference to the original filename.

    If the object is actually None, return None.
    """
    if dxlink is not None:
        dxf = dxpy.DXFile(dxlink)
        fname = dxf.describe()['name']
        if ' ' in fname:
            fname_nospace = fname.replace(' ', '_')
            print("Filename", repr(fname), "contains spaces; renaming to",
                  repr(fname_nospace))
            fname = fname_nospace
        dxpy.download_dxfile(dxf.get_id(), fname)
        return fname


def maybe_gunzip(fname, base, ext):
    if fname and 'gzip' in magic.from_file(fname):
        newf = safe_fname(base, ext)
        sh("gunzip", fname, "-c >", newf)
        fname = newf
    return fname


def safe_fname(base, ext):
    """Ensure the chosen filename will not overwrite an existing file.

    If it would, insert some garbage in the middle of the chosen filename to
    avoid the naming conflict.
    """
    fname = base + "." + ext
    if os.path.exists(fname):
        fname = tempfile.mktemp(prefix=base + "-", suffix="." + ext, dir=".")
    return fname


def sh(*command):
    """Run a shell command."""
    cmd = " ".join(map(str, command))
    print("$", cmd)
    subprocess.check_call(cmd, shell=True)
    print()


def shout(*command):
    """Run a shell command and capture standard output."""
    cmd = " ".join(map(str, command))
    print("$", cmd)
    result = subprocess.check_output(cmd, shell=True)
    print()
    return result


# _____________________________________________________________________________

dxpy.run()
