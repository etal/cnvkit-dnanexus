<!-- dx-header -->
# CNVkit (DNAnexus Platform App)

Run the CNVkit pipeline to call copy number variants.

<!-- /dx-header -->

CNVkit is a software toolkit to infer and visualize copy number from high-throughput DNA sequencing data. It is designed for use with hybrid capture, including both whole-exome and custom target panels, and short-read sequencing platforms such as Illumina and Ion Torrent.

This app runs CNVkit's standard pipeline, equivalent to the "batch" command.

For details on CNVkit methods and best practices, see the online documentation:

* https://cnvkit.readthedocs.io/en/latest/


## Prerequisites

This app uses three inputs:

* **Sample BAM files**: Sequencing reads from each sample, mapped to the same reference genome. These can be cases (e.g. from tumor or constitutional tissue) or controls (usually from blood or saliva), or both.
* **Reference genome FASTA**: The same genome sequence that the BAMs were aligned to. Gzip compression is optional and will be handled automatically.
* **Capture region BED**: If targeted resequencing was performed, i.e. hybrid capture or targeted amplicon sequencing, then the capture region coordinates in UCSC BED format are specified here. For hybrid capture, use the captured/baited regions, not the targets. The capture kit vendor should be able to provide this file if you don't already have it.

The control samples are optional; CNVkit can perform control-free calling, but best results are usually obtained by using a pool of controls. The controls do not have to be matched to the same patients as the case samples (as in tumor-normal pairs), but controls should have been processed with the same lab preparation and sequencing protocol as the case samples, and preferably have no large copy number alterations (CNAs).


## Basic usage

To build a reference profile and then use it to call copy number, run this app with the inputs described above:

* `case_bams`: Case samples, e.g. diseased patient samples
* `normal_bams`: Normal samples from healthy individuals or cells (if available)
* `baits`: Capture regions, if targeted sequencing.
* `fasta`: Reference genome sequence

The app will first use the normal BAMs, capture regions, and genome sequence to construct a reference profile, which becomes the app's `cn_reference` output. This reference profile is reusable for subsequent analyses (see below).

Each of the case BAMs is then normalized to the reference profile, resulting in a table of normalized copy ratios (`copy_ratios` output). The copy ratios are segmented into discrete copy number regions (`copy_segments` output), which are then further processed to remove likely false positives, detect single-exon variants, and assign integer copy number (`call_segments` output).

Additional plots and summary tables are generated from the copy ratios and segments -- see the app input help text and CNVkit documentation for details on each.


## Build a reusable reference profile

Starting from scratch, use the `normal_bams`, `baits`, and `fasta` inputs as described above, but leave `case_bams` empty, to construct a new reference profile (`cn_reference` output).

Other relevant inputs used here:

* `method`: Target capture method, if used in sequencing. `hybrid` (hybridization capture) is the default; `amplicon` (PCR-based targeted amplicon capture) and `wgs` (whole genome sequencing -- no capture performed) are also supported.
* `annotation`: Gene annotation for the reference genome, e.g. UCSC refFlat.txt. Used to apply gene names to targets/baits, if the `baits` BED file does not have a fourth column of region names.
* `exclude_access`: Known problematic (e.g. unmappable) regions to exclude from off-target bins.
* `target_avg_size`, `antitarget_avg_size`: Average on- and off-target bin sizes to use. If any control samples are given, these values will be determined automatically from the control samples' sequencing coverage depth. Typically you do not need to specify these inputs.


## Call CNVs using a pre-built reference profile

Use the `cn_reference` output from the reference-building step (above) as the input here, alongside the test/tumor samples in the `case_bams` field.

Other relevant inputs used here:

* `snv_vcfs`: SNV calls for each case sample, including germline SNPs, in VCF format. Must be given in the same order as the case BAMs. See [VCF formatting guidelines](https://cnvkit.readthedocs.io/en/stable/fileformats.html#vcf).
* `ploidy`: Autosome ploidy. Number of copies of each chromosome expected in normal cells (excluding sex chromosomes and mitochondria). For humans, this is 2, the default.
* `purity`: Estimated tumor purity. Fraction, 0--1, of cells in each test sample that are lesional, rather than immune or normal-cell contamination. For germline samples, this is 1.0 and does not need to be specified.
