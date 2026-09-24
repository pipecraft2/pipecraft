#!/bin/bash

# Run MetaWorks v1.13.0 ESV workflow (snakefile_ESV).
# PipeCraft writes /input/metaworks_out/config_ESV.yaml and adapters.fasta
# before this script starts. The snakefile calls perl_scripts and bold.hmm
# by relative path, so the working directory stays the MetaWorks install.

set -euo pipefail

echo "=========================================="
echo "MetaWorks ESV"
echo "=========================================="
echo "Date: $(date)"

CONFIG="/input/metaworks_out/config_ESV.yaml"
if [[ ! -f "$CONFIG" ]]; then
    echo "ERROR: MetaWorks config was not written: $CONFIG" >&2
    exit 1
fi

INSTALL_DIR="/opt/tools/Metaworks1.13.0"
if [[ ! -d "$INSTALL_DIR" ]]; then
    echo "ERROR: MetaWorks install directory not found: $INSTALL_DIR" >&2
    exit 1
fi

cd "$INSTALL_DIR"
eval "$(conda shell.bash hook)"
# openjdk_activate.sh reads JAVA_HOME before conda defines it.
# set -u would abort on that line.
set +u
conda activate MetaWorks_v1.13.0
set -u

JOBS="${jobs:-1}"
echo "snakemake --jobs ${JOBS} --snakefile snakefile_ESV --configfile ${CONFIG}"

status=0
snakemake --jobs "${JOBS}" --snakefile snakefile_ESV --configfile "$CONFIG" || status=$?

if [[ -n "${HOST_UID:-}" && -n "${HOST_GID:-}" ]]; then
    chown -R "${HOST_UID}:${HOST_GID}" /input/metaworks_out || true
fi

if [[ "$status" -ne 0 ]]; then
    echo "ERROR: snakemake exited with status ${status}" >&2
    exit "$status"
fi

echo "workingDir=/input/metaworks_out"
echo "ESV table: /input/metaworks_out/ESV.table"
echo "Taxonomy: /input/metaworks_out/taxonomy.csv"
echo "Sequence FASTA and results.csv index are in the same folder."
