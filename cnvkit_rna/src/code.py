#!/usr/bin/env python
from __future__ import division, print_function

import os
import shlex
import subprocess
import tempfile
#from glob import glob
from os.path import basename, isfile
import re

import magic
import psutil

import dxpy


@dxpy.entry_point('main')
def main(case_counts=None, normal_counts=None, vcfs=None, cn_reference=None,
         baits=None, fasta=None, annotation=None,
         method='hybrid', haploid_x_reference=False, drop_low_coverage=False,
         exclude_access=None,
         antitarget_avg_size=None, target_avg_size=None,
         purity=None, ploidy=None, do_parallel=True):
    cnvkit_docker("version")


    sh("mkdir -p /workdir")

    output = {'copy_ratios': [], 'copy_segments': [],
              'scatters_png': [], 'genemetrics': [],
             }
             # Aggregate -- single file, no need to initialize list
             # seg, sexes, metrics, heatmap_pdf

    # Workflow:
    #   import-rna
    #   then, for each output .cnr, do segmentation &c.
    #   -> .cns, -scatter.png, -genemetrics.tsv
    #   -> .seg, .sex.tsv, -metrics.tsv, -heatmap.pdf

    # 1. import-rna

    # TODO put required files in resources/
    rna_cmd = ['import-rna',
            '--format', counts_format,
            '--gene-counts', '/',
            '--correlations', '/',
            ]
    if normal_counts:
        rna_cmd.extend(['--normal'] + normal_counts)
    if output_fname:
        rna_cmd.extend(['--output', output_fname])
    cnvkit_docker(*rna_cmd)

    # 2. segment, scatter, genemetrics

    sample_ids = []
    for sample_id in sample_ids:
        cnr_fname = sample_id + ".cnr"

        cns_fname = sample_id + ".cns"
        cnvkit_docker('segment', sample_id + ".cnr", "--output", cns_fname)
        if not isfile(cns_fname):
            raise dxpy.AppError(
                "Segmentation failed silently, try running again with "
                "do_parallel=false and see logs")

        genemetrics_fname = safe_fname(sample_id, "genemetrics.csv")
        cnvkit_docker('genemetrics', cnr_fname, '--segment', cns_fname,
                      '--min-probes', '3', '--threshold', '0.2',
                      '--output', genemetrics_fname)

        scatter_fname = safe_fname(sample_id, "-scatter.png")
        cnvkit_docker(['scatter', cnr_fname, '--segment', cns_fname,
                '--output', scatter_fname])

    sh("ls -Altr")  # Show the generated files in the DNAnexus log

    # 3. Aggregation

    seg_fname = safe_fname("rna-cnv-segments", "seg")
    cnvkit_docker('export', 'seg', ' '.join([sid + ".cns" for sid in sample_ids]),
                  '--output', seg_fname)

    return {'copy_ratios': upload_link(cnr_fname),
            'copy_segments': upload_link(cns_fname),
            'call_segments': upload_link(call_fname),
            'genemetrics': upload_link(genemetrics_fname),
            'breaks': upload_link(breaks_fname),
            'vcfs': upload_link(vcf_fname),
            'diagram': upload_link(sample_id + '-diagram.pdf'),
            'scatters': upload_link(sample_id + '-scatter.png'),
           }


@dxpy.entry_point('aggregate_outputs')
def aggregate_outputs(copy_ratios, copy_segments, haploid_x_reference):
    sh("mkdir -p /workdir")
    all_cnr = [download_link(cnr) for cnr in copy_ratios]
    all_cns = [download_link(cns) for cns in copy_segments]

    metrics_fname = safe_fname("cnv-metrics", "csv")
    cnvkit_docker("metrics", " ".join(all_cnr), "-s", " ".join(all_cns),
                  "--output", metrics_fname)

    seg_fname = safe_fname("cnv-segments", "seg")
    cnvkit_docker("export seg", " ".join(all_cns), "--output", seg_fname)

    if len(all_cns) > 1:
        heat_fname = safe_fname("cnv-heatmap", "pdf")
        cnvkit_docker("heatmap", "-d", " ".join(all_cns), "--output", heat_fname)
    else:
        heat_fname = None

    sex_fname = safe_fname("cnv-sex", "csv")
    sex_cmd = ["sex", "--output", sex_fname] + all_cnr
    if haploid_x_reference:
        sex_cmd.append("--haploid-x-reference")
    cnvkit_docker(*sex_cmd)

    outputs = {'metrics': upload_link(metrics_fname),
               'seg': upload_link(seg_fname),
               'sexes': upload_link(sex_fname)}
    if heat_fname:
        outputs['heatmap_pdf'] = upload_link(heat_fname)
    return outputs


@dxpy.entry_point('concat_pdfs')
def concat_pdfs(pdfs, name): # name="cnv-diagrams"
    if not pdfs:
        out_ref = None
    elif len(pdfs) == 1:
        # NB: make sure this operates on the list, not the link/dict
        out_ref = pdfs[0]
    else:
        # Download & concat
        pdf_fnames = [download_link(f) for f in pdfs]
        out_fname = safe_fname("cnv-" + name, "pdf")
        sh("pdfunite", " ".join(pdf_fnames), out_fname)
        out_ref = upload_link(out_fname)
    return {"pdf": out_ref}


def concat_pdfs_link(pdf_refs, name):
    job = dxpy.new_dxjob(fn_name='concat_pdfs',
            fn_input={'pdfs': pdf_refs, 'name': name})
    return job.get_output_ref('pdf')



# _____________________________________________________________________________

def download_link(dxlink):
    """Download a DNAnexus object reference to the original filename.

    If the object is actually None, return None.
    """
    if dxlink is not None:
        if isinstance(dxlink, dxpy.DXFile):
            dxf = dxlink
        else:
            dxf = dxpy.DXFile(dxlink)
        fname = dxf.describe()['name']
        if ' ' in fname:
            fname_nospace = fname.replace(' ', '_')
            print("Filename", repr(fname), "contains spaces; renaming to",
                  repr(fname_nospace))
            fname = fname_nospace
        dxpy.download_dxfile(dxf.get_id(), fname)
        return fname


def upload_link(fname):
    """Upload a local file to the DNAnexus platform and return a link."""
    if fname:
        return dxpy.dxlink(dxpy.upload_local_file(fname))


def fbase(fname):
    """Strip directory and extensions from a filename, just like CNVkit."""
    base = basename(fname)
    # Gzip extension usually follows another extension
    if base.endswith('.gz'):
        base = base[:-3]
    # Cases to drop more than just the last dot
    for ext in ('.recal.bam', '.deduplicated.realign.bam'):
        if base.endswith(ext):
            base = base[:-len(ext)]
            break
    else:
        # By default, just drop the last part of the extension
        base = base.rsplit('.', 1)[0]
    return base


def maybe_gunzip(fname, base, ext):
    """Unpack a gzipped file, if necessary."""
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
    """Log whether or not each input is an existing file."""
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
                     "etal/cnvkit:0.9.4", "cnvkit.py"]
    sh(*(docker_prefix + list(args)))


# _____________________________________________________________________________

dxpy.run()
