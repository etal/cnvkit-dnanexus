#!/bin/bash
# Run analysis on all TCRB paired samples
set -ex
dx select project-FFFkxQ00vzzYqffyGg9XyQXQ
dx run cnvkit_batch --watch -y \
    -imethod=hybrid \
    -ibaits=file-FFjBx4j0vzzv8QXVP877YYfk \
    -ifasta=file-B6ZY7VG2J35Vfvpkj8y0KZ01 \
    -icase_bams=file-FFGKb6Q0vzzgkb2b0g3p8KGP \
    -icase_bams=file-FFGKvZQ0vzzQ96pFF6Zk6456 \
    -inormal_bams=file-FFGK5jQ0vzzQVybZF6j27kK1 \
    --name="CNVkit $0" \
    --folder=tmp/
