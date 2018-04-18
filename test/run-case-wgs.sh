#!/bin/bash
# Run analysis on all TCRB paired samples
set -ex
dx select project-FFFkxQ00vzzYqffyGg9XyQXQ
dx run cnvkit_batch --watch -y \
    -imethod=wgs \
    -ifasta=file-B6ZY7VG2J35Vfvpkj8y0KZ01 \
    -itumor_bams=file-FFGKb6Q0vzzgkb2b0g3p8KGP \
    -ido_parallel=false \
    --name="CNVkit $0" \
    --folder=tmp/
