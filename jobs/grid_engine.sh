#!/bin/bash
#$ -N contact_analysis
#$ -cwd
#$ -S /bin/bash
# Queue, memory, walltime and CPU directives depend on your cluster.
set -euo pipefail
: "${CONTACT_PYTHON:?Set CONTACT_PYTHON to the absolute Python executable in your environment}"
: "${STAGE:?Set STAGE to 01_center or 02_calculate}"
case "$STAGE" in
    01_center|02_calculate) ;;
    *) echo "Invalid STAGE: $STAGE" >&2; exit 2 ;;
esac
set -- "$STAGE/run.py"
if [[ -n "${SIMULATION:-}" ]]; then set -- "$@" --simulation "$SIMULATION"; fi
exec "$CONTACT_PYTHON" -u "$@"
