#!/bin/bash
# Run analysis on all TCRB paired samples
set -ex
dx select project-FFFkxQ00vzzYqffyGg9XyQXQ
dx run cnvkit_batch --watch -y \
    -imethod=wgs \
    -ifasta=file-B6ZY7VG2J35Vfvpkj8y0KZ01 \
    -itumor_bams=file-FFGKb6Q0vzzgkb2b0g3p8KGP \
    -itumor_bams=file-FFFxP4j0vzzQ96pFF6ZfpxXf \
    -itumor_bams=file-FFGKvZQ0vzzQ96pFF6Zk6456 \
    -inormal_bams=file-FFGK5jQ0vzzQVybZF6j27kK1 \
    -inormal_bams=file-FFGQY2j0vzzbbp1pF6B1FVjQ \
    -inormal_bams=file-FFGQF700vzzQVgb60gQBZ1f4 \
    -ido_parallel=false \
    --name="CNVkit $0" \
    --folder=tmp/
