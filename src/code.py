#!/usr/bin/env python
# cn-py 0.0.1
# Generated by dx-app-wizard.
#
# Basic execution pattern: Your app will run on a single machine from
# beginning to end.
#
# See https://wiki.dnanexus.com/Developer-Portal for documentation and
# tutorials on how to modify this file.
#
# DNAnexus Python Bindings (dxpy) documentation:
#   http://autodoc.dnanexus.com/bindings/python/current/
from __future__ import division, print_function

import os
import subprocess
import tempfile
from glob import glob

import dxpy


@dxpy.entry_point('main')
def main(tumor_bams=None, normal_bams=None, cn_reference=None, baits=None,
         fasta=None, access=None, annotation=None):

    if cn_reference and any((normal_bams, baits, fasta, access, annotation)):
        raise dxpy.AppError("Reference profile (cn_reference) cannot be used "
                            "alongside normal_bams, baits, fasta, access "
                            "or annotation")

    # Install R dependencies
    install_bioc()
    # Install CNVkit itself
    # (NB: doing this from dxapp.json tries and fails to rebuild matplotlib)
    sh("pip install --user cnvkit")
    os.environ["PATH"] += ":" + os.environ["HOME"] + "/.local/bin"

    # Initialize and download file inputs to the local file system
    cn_reference = download_link(cn_reference)
    baits = download_link(baits)
    fasta = download_link(fasta)
    access = download_link(access)
    annotation = download_link(annotation)
    if tumor_bams is not None:
        tumor_bams = map(download_link, tumor_bams)
    if normal_bams is not None:
        normal_bams = map(download_link, normal_bams)

    out_fnames = run_cnvkit(tumor_bams, normal_bams, cn_reference,
                            baits, fasta, access, annotation)

    # Upload file outputs from the local file system to the platform.
    def upload(key):
        return dxpy.dxlink(dxpy.upload_local_file(out_fnames[key]))

    def upload_list(key):
        return [dxpy.dxlink(dxpy.upload_local_file(fname))
                for fname in out_fnames[key]]

    return {
        "cn_reference": upload("reference"),
        "copy_ratios": upload_list("cnr"),
        "copy_segments": upload_list("cns"),
        "gainloss": upload_list("gainloss"),
        "breaks": upload_list("breaks"),
        "nexus_cn": upload_list("nexus"),
        "metrics": upload("metrics"),
        "genders": upload("genders"),
        "scatter_pdf": upload("scatter_pdf"),
        "diagram_pdf": upload("diagram_pdf"),
    }


def run_cnvkit(tumor_bams, normal_bams, reference, baits, fasta, access,
               annotation):
    """Run the CNVkit pipeline.

    Returns a dict of the generated file names.
    """
    print("Determining if the given reference profile is male or female")
    if shout("cnvkit.py gender", reference, "| cut -f 2").strip() == "Female":
        yflag = ""
    else:
        yflag = "-y"

    print("Running the main CNVkit pipeline")
    command = ["cnvkit.py batch -p 0", yflag]

    if tumor_bams:
        command.extend(tumor_bams)
        command.extend(["--scatter", "--diagram"])
    if reference:
        command.extend(["-r", reference])
    else:
        reference = safe_fname("cnv-reference", "cnn")
        command.extend("--output-reference", reference, "-n")
        if normal_bams:
            command.extend(normal_bams)
        if baits:
            command.extend(["-t", baits, "--split", "--short-names"])
        if fasta:
            command.extend(["-f", fasta])
            if not access:
                access = safe_fname("access-10k", "bed")
                sh("genome2access.py -s 10000", fasta, "-o", access)
        if access:
            command.extend(["-g", access])
        if annotation:
            command.extend(["--annotate", annotation])

    sh(*command)
    sh("ls -Altr")  # Show the generated files in the DNAnexus log

    # Collect the outputs
    all_cnr = glob("*.cnr")
    all_cns = glob("*.cns")
    all_gainloss = []
    all_breaks = []
    all_nexus = []
    for acnr, acns in zip(all_cnr, all_cns):
        name = acnr.split('.')[0]

        gainloss = name + "gainloss.csv"
        sh("cnvkit.py gainloss", acnr, "-s", acns, "-m 3 -t 0.3", yflag,
           "-o", gainloss)
        all_gainloss.append(gainloss)

        breaks = name + "-breaks.csv"
        sh("cnvkit.py breaks", acnr, acns, "-m 3", "-o", breaks)
        all_breaks.append(breaks)

        nexus = name + ".nexus"
        sh("cnvkit.py export nexus-basic", name + ".cnr", "-o", nexus)
        all_nexus.append(nexus)

    metrics = safe_fname("metrics", "csv")
    sh("cnvkit.py metrics", " ".join(all_cnr), "-s", " ".join(all_cns),
       "-o", metrics)

    genders = safe_fname("gender", "csv")
    sh("cnvkit.py gender", yflag, "-o", genders, *all_cnr)

    scatter_pdf = safe_fname("cnv-scatters", "pdf")
    sh("pdfunite", " ".join(glob("*-scatter.pdf")), scatter_pdf)

    diagram_pdf = safe_fname("cnv-diagrams", "pdf")
    sh("pdfunite", " ".join(glob("*-diagram.pdf")), diagram_pdf)

    return {
        "reference": reference,
        "cnr": all_cnr,
        "cns": all_cns,
        "gainloss": all_gainloss,
        "breaks": all_breaks,
        "nexus": all_nexus,
        "metrics": metrics,
        "genders": genders,
        "scatter_pdf": scatter_pdf,
        "diagram_pdf": diagram_pdf,
    }


def install_bioc():
    """In R, download and install dependencies including Bioconductor."""
    r_version = shout("R --version | head -n1 | cut -d' ' -f3 | cut -d. -f1-2"
                     ).strip()
    r_libs_user_dir = "/".join([os.environ["HOME"],
                                "R/x86_64-pc-linux-gnu-library",
                                r_version])
    # eg. "/home/dnanexus/R/x86_64-pc-linux-gnu-library/2.14/"
    os.environ["R_LIBS_USER"] = r_libs_user_dir
    if os.path.isdir(r_libs_user_dir):
        print("R_LIBS_USER dir exists:", r_libs_user_dir)
    else:
        print("Creating R_LIBS_USER dir:", r_libs_user_dir)
        os.makedirs(r_libs_user_dir)
    sh("Rscript -e '%s'"
       % "; ".join([
           '.libPaths(Sys.getenv("R_LIBS_USER"))',
           'source("http://bioconductor.org/biocLite.R")',
           'biocLite()',
           'biocLite("DNAcopy")',
           'install.packages("PSCBS", repos="http://cran.us.r-project.org")',
       ]))


# _____________________________________________________________________________

def download_link(dxlink):
    """Download a DNAnexus object reference to the original filename.

    If the object is actually None, return None.
    """
    if dxlink is not None:
        dxf = dxpy.DXFile(dxlink)
        fname = dxf.describe()['name']
        dxpy.download_dxfile(dxf.get_id(), fname)
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