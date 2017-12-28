#!/usr/bin/env python
from __future__ import division, print_function

import os
import shlex
import subprocess
import tempfile
from glob import glob

import magic
import psutil

import dxpy


@dxpy.entry_point('main')
def main(tumor_bams=None, normal_bams=None, vcfs=None, cn_reference=None,
         baits=None, fasta=None, annotation=None,
         method='hybrid', is_male_normal=True, drop_low_coverage=False,
         antitarget_avg_size=None, target_avg_size=None,
         purity=None, ploidy=None, do_parallel=True):
    if not tumor_bams and not normal_bams:
        raise dxpy.AppError("Must provide tumor_bams or normal_bams (or both)")
    if cn_reference and any((normal_bams, baits, fasta, annotation)):
        raise dxpy.AppError("Reference profile (cn_reference) cannot be used "
                            "alongside normal_bams, baits, fasta, "
                            "or annotation")
    if method != "wgs" and tumor_bams and not any((baits, cn_reference)):
        raise dxpy.AppError("Need cn_reference or baits to process tumor_bams")
    if tumor_bams:
        purities = validate_per_tumor(purity, len(tumor_bams), "purity values",
                                      lambda p: 0 < p <= 1)
        ploidies = validate_per_tumor(ploidy, len(tumor_bams), "ploidy values",
                                      lambda p: p > 0)
        vcfs = validate_per_tumor(vcfs, len(tumor_bams), "VCF files")
    else:
        purities = ploidies = None

    print("Downloading file inputs to the local file system")
    cn_reference = download_link(cn_reference)
    baits = download_link(baits)
    fasta = download_link(fasta)
    annotation = download_link(annotation)
    if tumor_bams is not None:
        tumor_bams = map(download_link, tumor_bams)
    if normal_bams is not None:
        normal_bams = map(download_link, normal_bams)
    if vcfs is not None:
        vcfs = map(download_link, vcfs)

    # If these input files are gzipped, decompress them
    fasta = maybe_gunzip(fasta, "ref", "fa")
    annotation = maybe_gunzip(annotation, "annot", "txt")

    out_fnames = run_cnvkit(tumor_bams, normal_bams, vcfs, cn_reference, baits,
                            fasta, annotation, method, is_male_normal,
                            drop_low_coverage, antitarget_avg_size,
                            target_avg_size, purities, ploidies, do_parallel)

    print("Uploading local file outputs to the DNAnexus platform")
    output = {}
    for filekey in ("cn_reference", "seg", "metrics", "sexes", "scatter_pdf",
                    "diagram_pdf"):
        if filekey in out_fnames:
            output[filekey] = dxpy.dxlink(
                dxpy.upload_local_file(out_fnames[filekey]))
    for listkey in ("copy_ratios", "copy_segments", "gainloss", "breaks"):
        if listkey in out_fnames:
            output[listkey] = [dxpy.dxlink(dxpy.upload_local_file(fname))
                               for fname in out_fnames[listkey]]
    return output


def validate_per_tumor(values, n_expected, title, criterion=None):
    """Ensure a per-tumor input array matches the number of tumor BAMs, etc.

    Also allow a value of None to skip checks & downstream processing, or a
    single value to apply to every tumor sample.

    Returns: list of values the same length as `n_expected`.
    """
    if values is None:
        out_vals = [None] * n_expected
    else:
        if criterion is not None and not all(map(criterion, values)):
            raise dxpy.AppError(
                """Tumor {} must all be between 0 and 1; got: {}"""
                .format(title, values))
        if len(values) == n_expected:
            out_vals = values
        elif len(values) == 1 and n_expected > 1:
            out_vals = [values] * n_expected
        else:
            raise dxpy.AppError(
                """Number of tumor {} specified ({}) does not
                match the number of tumor BAM files given ({})"""
                .format(title, len(values), n_expected))
    return out_vals


def run_cnvkit(tumor_bams, normal_bams, vcfs, reference, baits, fasta,
               annotation, method, is_male_normal, drop_low_coverage,
               antitarget_avg_size, target_avg_size, purities, ploidies,
               do_parallel):
    """Run the CNVkit pipeline.

    Returns a dict of the generated file names.
    """
    sh("mkdir -p /workdir")

    print("Running the main CNVkit pipeline")
    command = ["batch", "-m", method]
    if tumor_bams:
        command.extend(tumor_bams)
        command.extend(["--scatter", "--diagram"])
    if reference:
        # Use the given reference
        command.extend(["-r", reference])
    else:
        # Build a new reference
        reference = safe_fname("cnv-reference", "cnn")
        command.extend(["--output-reference", reference,
                        "-n"])
        if normal_bams:
            command.extend(normal_bams)
        if baits:
            command.extend(["-t", baits, "--short-names"])
        if fasta:
            command.extend(["-f", fasta])
        if annotation:
            command.extend(["--annotate", annotation])
        if antitarget_avg_size:
            command.extend(["--antitarget-avg-size", antitarget_avg_size])
        if target_avg_size:
            command.extend(["--target-avg-size", target_avg_size])
    if drop_low_coverage:
        command.append("--drop-low-coverage")
    if do_parallel:
        command.extend(["-p", str(psutil.cpu_count(logical=True))])
    yflag = "-y" if is_male_normal else ""
    command.append(yflag)

    cnvkit_docker(*command)
    sh("ls -Altr")  # Show the generated files in the DNAnexus log

    # Collect the outputs
    if not tumor_bams:
        return {"cn_reference": reference}

    all_cnr = glob("*.cnr")
    all_cns = glob("*.cns")
    if not len(all_cnr) == len(all_cns):
        raise dxpy.AppError(
            "Segmentation failed silently, try running again with "
            "do_parallel=false and see logs")

    all_gainloss = []
    all_breaks = []
    all_segmetrics = []
    all_call = []
    for acnr, acns, vcf, purity, ploidy in zip(all_cnr, all_cns, vcfs,
                                               purities, ploidies):
        name = acnr.split('.')[0]

        gainloss = name + ".gainloss.csv"
        cnvkit_docker("gainloss", acnr, "-s", acns, "-m 3 -t 0.3", yflag,
                      "-o", gainloss)
        all_gainloss.append(gainloss)

        breaks = name + ".breaks.csv"
        cnvkit_docker("breaks", acnr, acns, "-m 3", "-o", breaks)
        all_breaks.append(breaks)

        segmetrics = name + ".segmetrics.cns"
        cnvkit_docker("segmetrics", acnr, "-s", acns, "--ci -a .1",
                      "-o", segmetrics)
        all_segmetrics.append(segmetrics)

        call_fname = name + ".call.cns"
        call_opts = ["call", segmetrics, "-o", call_fname,
                     yflag, "--filter", "ci", "--center"]
        if purity or ploidy:
            call_opts.append("-m clonal")
            if purity:
                call_opts.extend(["--purity", purity])
            if ploidy:
                call_opts.extend(["--ploidy", ploidy])
        else:
            call_opts.append("-m threshold")
        if vcf:
            call_opts.extend(["--vcf", vcf])
        cnvkit_docker(*call_opts)
        all_call.append(call_fname)

    seg = safe_fname("cnv-segments", "seg")
    cnvkit_docker("export seg", " ".join(all_cns), "-o", seg)

    sexes = safe_fname("cnv-sex", "csv")
    cnvkit_docker("sex", yflag, "-o", sexes, *all_cnr)

    metrics = safe_fname("cnv-metrics", "csv")
    cnvkit_docker("metrics", " ".join(all_cnr), "-s", " ".join(all_cns), "-o", metrics)

    outputs = {
        "cn_reference": reference,
        "copy_ratios": all_cnr,
        "copy_segments": all_segmetrics,
        "gainloss": all_gainloss,
        "breaks": all_breaks,
        "seg": seg,
        "sexes": sexes,
        "metrics": metrics,
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
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as exc:
        check_files(command[1:])
        raise exc
    print()


def check_files(maybe_filenames):
    fnames = []
    for fname in maybe_filenames:
        if isinstance(fname, basestring):
            fnames.extend(shlex.split(fname))
    for fname in fnames:
        if '.' in fname:
            # It might be a filename
            if os.path.exists(fname):
                print("File:", os.path.abspath(fname))
            else:
                print("Not file:", os.path.abspath(fname))


def cnvkit_docker(*args):
    """Run a CNVkit sub-command."""
    docker_prefix = ["dx-docker", "run",
                     "-v", "/home/dnanexus:/workdir", "-w", "/workdir",
                     "etal/cnvkit:0.9.1", "cnvkit.py"]
    sh(*(docker_prefix + list(args)))


# _____________________________________________________________________________

dxpy.run()
