#!/bin/bash
# Build a WGS flat reference (no samples, just the reference genome)
set -ex
dx select project-FFFkxQ00vzzYqffyGg9XyQXQ
dx run cnvkit_batch --watch -y \
    -imethod=wgs \
    -ifasta=file-B6ZY7VG2J35Vfvpkj8y0KZ01 \
    --name="CNVkit $0" \
    --folder=tmp/
