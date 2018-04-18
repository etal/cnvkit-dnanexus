#!/bin/bash
# Build an amplicon flat reference (no samples, just targets)
set -ex
dx select project-FFFkxQ00vzzYqffyGg9XyQXQ
dx run cnvkit_batch --watch -y \
    -imethod=amplicon \
    -ibaits=file-FFjBx4j0vzzv8QXVP877YYfk \
    -ifasta=file-B6ZY7VG2J35Vfvpkj8y0KZ01 \
    --name="CNVkit $0" \
    --folder=tmp/
