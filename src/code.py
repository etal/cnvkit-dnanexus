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
def main(case_bams=None, normal_bams=None, vcfs=None, cn_reference=None,
         baits=None, fasta=None, annotation=None,
         method='hybrid', haploid_x_reference=False, drop_low_coverage=False,
         exclude_access=None,
         antitarget_avg_size=None, target_avg_size=None,
         purity=None, ploidy=None, do_parallel=True):
    cnvkit_docker("version")

    # Validate inputs
    # (from cnvlib.commands._cmd_batch)
    if cn_reference:
        bad_flags = [flag
                     for is_used, flag in (
                         (normal_bams is not None,  'normal_bams'),
                         (fasta,                    'fasta'),
                         (baits,                    'baits'),
                         (annotation,               'annotation'),
                         (exclude_access,           'exclude_access'),
                         (target_avg_size,          'target_avg_size'),
                         (antitarget_avg_size,      'antitarget_avg_size'),
                         #(short_names,              'short_names'),
                         )
                     if is_used]
        if bad_flags:
            raise dxpy.AppError(
                    "If 'cn_reference' is given, options to construct a new "
                    "reference (%s) should not be used:" % ", ".join(bad_flags))
    elif method in ('hybrid', 'amplicon') and not baits:
        raise dxpy.AppError(
                "For the '%r' sequencing method, input 'baits' (at least) "
                "must be given to build a new reference if 'cn_reference' "
                "is not given." % baits)
    if case_bams:
        purities = validate_per_tumor(purity, len(case_bams), "purity values",
                                      lambda p: 0 < p <= 1)
        ploidies = validate_per_tumor(ploidy, len(case_bams), "ploidy values",
                                      lambda p: p > 0)
        vcfs = validate_per_tumor(vcfs, len(case_bams), "VCF files")
    else:
        purities = ploidies = None
    sh("mkdir -p /workdir")

    # If reference is not given, create one
    if not cn_reference:
        print("** About to call 'make_region_beds'")  # DBG
        targets, antitargets = make_region_beds(
                normal_bams, method, fasta, baits, annotation,
                exclude_access, antitarget_avg_size, target_avg_size)
        print("** Finished calling 'make_region_beds'")  # DBG
        normal_cvgs = []
        if normal_bams:
            # 'coverage' of each normal bam in a subjob
            for nbam in normal_bams:
                print("** About to launch 'run_coverage'")  # DBG
                job_cvg = dxpy.new_dxjob(fn_name='run_coverage',
                        fn_input={
                            'bam': nbam,
                            'targets': targets,
                            'antitargets': antitargets,
                            'do_parallel': do_parallel,
                            })
                normal_cvgs.append(job_cvg.get_output_ref('coverages'))
                print("** Got output ref from 'run_coverage'")  # DBG
        print("** About to launch 'run_reference'")  # DBG
        job_ref = dxpy.new_dxjob(fn_name='run_reference',
                fn_input={'coverages': normal_cvgs,
                          'fasta': fasta,
                          'targets': targets,
                          'antitargets': (antitargets
                              if method == 'hybrid' else None),
                          'haploid_x_reference': haploid_x_reference,
                          })
        cn_reference = job_ref.get_output_ref('cn_reference')
        print("** Got output ref from 'run_reference'")  # DBG
    output = {'cn_reference': cn_reference,
              'copy_ratios': [], 'copy_segments': [], 'call_segments': [],
              'genemetrics': [], 'breaks': [],
             }

    # Process each test/case/tumor individually using the given/built reference
    if case_bams:
        print("** About to process", len(case_bams), "'case_bams'")  # DBG
        diagram_pdfs = []
        scatter_pdfs = []
        for sample_bam, vcf, purity, ploidy in \
                zip(case_bams, vcfs, purities, ploidies):
            print("** About to launch 'run_sample'")  # DBG
            job_sample = dxpy.new_dxjob(fn_name='run_sample',
                    fn_input={
                        'sample_bam': sample_bam,
                        'method': method,
                        'cn_reference': cn_reference,
                        'vcf': vcf,
                        'purity': purity,
                        'ploidy': ploidy,
                        'drop_low_coverage': drop_low_coverage,
                        'haploid_x_reference': haploid_x_reference,
                        'do_parallel': do_parallel,
                        })
            for field in ('copy_ratios', 'copy_segments', 'call_segments',
                          'genemetrics', 'breaks', 'vcfs'):
                output[field].append(job_sample.get_output_ref(field))
            diagram_pdfs.append(job_sample.get_output_ref('diagram'))
            scatter_pdfs.append(job_sample.get_output_ref('scatter'))
            print("** Got outputs from 'run_sample'")  # DBG

        # Consolidate multi-sample outputs
        print("** About to launch 'aggregate_outputs'")  # DBG
        job_agg = dxpy.new_dxjob(fn_name='aggregate_outputs',
                fn_input={'copy_ratios': output['copy_ratios'],
                          'copy_segments': output['copy_segments'],
                          'haploid_x_reference': haploid_x_reference})
        for field in ('seg', 'heatmap_pdf', 'metrics', 'sexes'):
            output[field] = job_agg.get_output_ref(field)
        output['diagram_pdf'] = concat_pdfs_link(diagram_pdfs, 'diagrams')
        output['scatter_pdf'] = concat_pdfs_link(scatter_pdfs, 'scatters')
        print("** Got outputs from 'aggregate_outputs'")  # DBG

    print("** All done! Returning output:")
    from pprint import pprint
    pprint(output)
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


def make_region_beds(normal_bams, method, fasta, baits, annotation,
        exclude_access, antitarget_avg_size, target_avg_size):
    """Quickly calculate reasonable bin sizes from BAM read counts."""
    sh("mkdir -p /workdir")
    if method != 'amplicon':
        # Get 'access' from fasta
        fa_fname = maybe_gunzip(download_link(fasta), "ref", "fa")
        access_bed = safe_fname("access", "bed")
        cmd = ['access', fa_fname, '--output', access_bed]
        if exclude_access:
            if isinstance(exclude_access, (list, tuple)):
                for excl in exclude_access:
                    excl_fname = download_link(excl)
                    cmd.extend(['--exclude', excl_fname])
            else:
                raise dxpy.AppError("Unexpected type: %r" % exclude_access)
        if method == 'wgs':
            cmd.extend(['--min-gap-size', '100'])
        cnvkit_docker(*cmd)
    else:
        access_bed = None

    if annotation:
        annot_fname = maybe_gunzip(download_link(annotation), "annot", "txt")
    else:
        annot_fname = None

    if baits:
        bait_bed = download_link(baits)
    elif method == 'wgs':
        bait_bed = filter_bed_chroms(access_bed)
    else:
        bait_bed = None

    if target_avg_size or not normal_bams:
        # Use the given bin sizes with 'target' and 'antitarget'
        # (Or, if no samples, stick with default target size)
        tgt_bed = safe_fname('targets', 'bed')
        cmd_tgt = ['target', bait_bed, '--short-names', '--output', tgt_bed]
        if target_avg_size:
            cmd_tgt.extend(['--avg-size', target_avg_size])
        if annotation:
            cmd_tgt.extend(['--annotate', annot_fname])
        cnvkit_docker(*cmd_tgt)
        if method == 'hybrid':
            anti_bed = safe_fname('antitargets', 'bed')
            cmd_anti = ['antitarget', bait_bed, '--access', access_bed,
                        '--output', anti_bed]
            if antitarget_avg_size:
                cmd_anti.extend(['--avg-size', antitarget_avg_size])
            cnvkit_docker(*cmd_anti)
        else:
            anti_bed = None

    else:
        # Estimate bin sizes with 'autobin'
        # - wgs: targets ignored, access required
        #       - at this point bait_bed may be == access_bed
        # - hybrid: targets & access required
        # - amplicon: targets required, access ignored
        def midsize_dxfile(refs):
            dxfiles = [dxpy.DXFile(ref) for ref in refs]
            if len(dxfiles) == 1:
                return dxfiles[0]
            def dxsize(dxfile):
                """Get file size on the platform by API call."""
                return dxfile.describe()['size']
            return sorted(dxfiles, key=dxsize)[len(dxfiles) // 2 - 1]

        bam_fname = download_link(midsize_dxfile(normal_bams))
        cmd_autobin = ['autobin', bam_fname,
                       '--method', method, '--short-names']
        if annotation:
            cmd_autobin.extend(['--annotate', annot_fname])
        if baits:
            cmd_autobin.extend(['--targets', bait_bed])
            out_fname_base = fbase(bait_bed)
        if access_bed:
            cmd_autobin.extend(['--access', access_bed])
            if not baits:
                # For WGS
                #out_fname_base = fbase(access_bed)
                out_fname_base = fbase(bam_fname)
        tgt_bed = out_fname_base + '.target.bed'
        if method == 'hybrid':
            anti_bed = out_fname_base + '.antitarget.bed'
        else:
            anti_bed = None
        print("** About to run 'autobin' in docker")
        cnvkit_docker(*cmd_autobin)
        print("** Ran 'autobin' in docker")
    print("** About to upload links to", tgt_bed, "and", anti_bed)
    return upload_link(tgt_bed), upload_link(anti_bed)


def filter_bed_chroms(in_bed):
    """For WGS, don't survey alt contigs in the reference genome."""
    r_canonical = re.compile("(chr)?(\d+|[XYxy])\t")
    out_bed = fbase(in_bed) + ".filtered.bed"
    with open(in_bed) as infile:
        with open(out_bed, 'w') as outfile:
            for line in infile:
                if r_canonical.match(line):
                    outfile.write(line)
            did_match = (outfile.tell > 0)
    if did_match:
        print("Replacing", in_bed, "with", out_bed,
              "keeping only canonical contig names")
        return out_bed
    print("None of the contig names in", in_bed,
          "recognized as canonical; won't filter")
    return in_bed


@dxpy.entry_point('run_coverage')
def run_coverage(bam, targets, antitargets, do_parallel):
    """Calculate coverage in the given regions from BAM read depths."""
    sh("mkdir -p /workdir")
    cmd_base = ['coverage']
    if do_parallel:
        cmd_base.extend(['--processes', str(psutil.cpu_count(logical=True))])
    bam_fname = download_link(bam)
    target_fname = download_link(targets)
    sample_id = fbase(bam_fname)
    tcov_fname = safe_fname(sample_id, "targetcoverage.cnn")
    cnvkit_docker(*(cmd_base +
        [bam_fname, target_fname, '--output', tcov_fname]))
    coverages = [tcov_fname]
    if antitargets:
        antitarget_fname = download_link(antitargets)
        acov_fname = safe_fname(sample_id, "antitargetcoverage.cnn")
        cnvkit_docker(*(cmd_base +
            [bam_fname, antitarget_fname, '--output', acov_fname]))
        coverages.append(acov_fname)
    return {'coverages': [upload_link(f) for f in coverages]}


@dxpy.entry_point('run_reference')
def run_reference(coverages, fasta, targets, antitargets, haploid_x_reference):
    """Compile a coverage reference from the given files (normal samples)."""
    sh("mkdir -p /workdir")
    fa_fname = maybe_gunzip(download_link(fasta), "ref", "fa")
    out_fname = safe_fname("cnv-reference", "cnn")
    cmd = ['reference', '--fasta', fa_fname, '--output', out_fname]
    if haploid_x_reference:
        cmd.append('--haploid-x-reference')
    if coverages:
        cmd.extend(download_link(f) for f in flatten_dxarray(coverages))
    else:
        # Flat reference -- need targets, antitargets
        cmd.extend(['--targets', download_link(targets)])
        if antitargets:
            cmd.extend(['--antitargets', download_link(antitargets)])
        else:
            # XXX shim for <=0.9.1
            anti_fname = safe_fname("anti-MT", ".bed")
            sh("touch", anti_fname)
            cmd.extend(['--antitargets', anti_fname])
    cnvkit_docker(*cmd)
    return {'cn_reference': upload_link(out_fname)}


def flatten_dxarray(listish):
    """Ensure an array-of-refs is not an array-of-arrays.

    Specifically, a list of job outputs might resolve to a list-of-lists if each
    job output is itself a list -- but it's hard to know or handle that until
    the blocking job completes.

    Use this function within the subjob that takes an array as input.
    """
    if not isinstance(listish, list):
        return [listish]
    out = []
    for elem in listish:
        if isinstance(listish, list):
            out.extend(flatten_dxarray(elem))
        else:
            out.append(elem)
    return out


@dxpy.entry_point('run_sample')
def run_sample(sample_bam, method, cn_reference, vcf, purity, ploidy,
    drop_low_coverage, haploid_x_reference, do_parallel):
    """Run the CNVkit pipeline.

    Returns a dict of the generated file names.
    """
    sh("mkdir -p /workdir")

    ref_fname = download_link(cn_reference)
    bam_fname = download_link(sample_bam)
    sample_id = fbase(bam_fname)

    cmd = ['batch', bam_fname, '--method', method, '--reference', ref_fname,
           '--scatter', '--diagram']
    if do_parallel:
        cmd.extend(['--processes', str(psutil.cpu_count(logical=True))])
    shared_opts = []
    if drop_low_coverage:
        shared_opts.append('--drop-low-coverage')
    if haploid_x_reference:
        shared_opts.append('--haploid-x-reference')
    cmd.extend(shared_opts)

    cnvkit_docker(*cmd)
    sh("ls -Altr")  # Show the generated files in the DNAnexus log

    cnr_fname = sample_id + ".cnr"
    cns_fname = sample_id + ".cns"
    if not isfile(cnr_fname):
        raise dxpy.AppError(
            "CNVkit failed silently, try running again with "
            "do_parallel=false and see logs")
    if not isfile(cns_fname):
        raise dxpy.AppError(
            "Segmentation failed silently, try running again with "
            "do_parallel=false and see logs")

    # Post-processing

    sm_fname = safe_fname(sample_id, "segmetrics.cns")
    sm_cmd = ['segmetrics', cnr_fname, '-s', cns_fname, '--output', sm_fname,
              '--ci', '--alpha', '0.5',
              '--bootstrap', '50' if method == 'wgs' else '100']
    if drop_low_coverage:
        sm_cmd.append('--drop-low-coverage')
    cnvkit_docker(*sm_cmd)

    vcf_fname = sample_id + ".vcf"
    vcf_cmd = ['export', 'vcf', cns_fname, '--output', vcf_fname]
    if haploid_x_reference:
        vcf_cmd.append('--haploid-x-reference')
    if ploidy:
        vcf_cmd.extend(['--ploidy', ploidy])
    cnvkit_docker(*vcf_cmd)

    call_fname = safe_fname(sample_id, "call.cns")
    call_cmd = ['call', sm_fname, '--output', call_fname,
                '--center', '--filter', 'ci']
    #call_cmd.extend(shared_opts)
    if haploid_x_reference:
        call_cmd.append('--haploid-x-reference')
    if vcf:
        call_cmd.extend(['--vcf', download_link(vcf)])
    if purity or ploidy:
        call_cmd.extend(['--method', 'clonal'])
        if purity:
            call_cmd.extend(['--purity', purity])
        if ploidy:
            call_cmd.extend(['--ploidy', ploidy])
    else:
        call_cmd.extend(['--method', 'threshold'])
    cnvkit_docker(*call_cmd)

    genemetrics_fname = safe_fname(sample_id, "genemetrics.csv")
    gm_cmd = ['genemetrics', cnr_fname, '--segment', call_fname,
              '--min-probes', '3', '--threshold', '0.2',
              '--output', genemetrics_fname]
    gm_cmd.extend(shared_opts)
    cnvkit_docker(*gm_cmd)

    breaks_fname = safe_fname(sample_id, "breaks.csv")
    cnvkit_docker('breaks', cnr_fname, call_fname, '--min-probes', '3',
                  '--output', breaks_fname)

    return {'copy_ratios': upload_link(cnr_fname),
            'copy_segments': upload_link(cns_fname),
            'call_segments': upload_link(call_fname),
            'genemetrics': upload_link(genemetrics_fname),
            'breaks': upload_link(breaks_fname),
            'vcfs': upload_link(vcf_fname),
            'diagram': upload_link(sample_id + '-diagram.pdf'),
            'scatter': upload_link(sample_id + '-scatter.pdf'),
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
def concat_pdfs(pdfs, name): # name="cnv-diagrams"|"cnv-scatters"
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
                     "etal/cnvkit:0.9.1", "cnvkit.py"]
    sh(*(docker_prefix + list(args)))


# _____________________________________________________________________________

dxpy.run()
