#!/bin/bash
#$ -N contact_analysis
#$ -cwd
#$ -S /bin/bash
# Queue, memory, walltime and CPU directives depend on your cluster.
set -euo pipefail
: "${CONTACT_PYTHON:?Set CONTACT_PYTHON to the absolute Python executable in your environment}"
: "${STAGE:?Set STAGE to 01_center, 02_basic_calculate, 02_calculate, or any calculate folder/script}"
case "$STAGE" in
    01_center|02_basic_calculate|02_calculate|02_calculate/*) ;;
    *) echo "Invalid STAGE: $STAGE" >&2; exit 2 ;;
esac
if [[ -f "$STAGE" ]]; then
    set -- "$STAGE"
elif [[ -f "$STAGE/run.py" ]]; then
    set -- "$STAGE/run.py"
else
    echo "Cannot find script for STAGE: $STAGE" >&2; exit 2
fi
if [[ -n "${SIMULATION:-}" ]]; then set -- "$@" --simulation "$SIMULATION"; fi
if [[ -n "${TASK:-}" ]]; then set -- "$@" --task "$TASK"; fi
exec "$CONTACT_PYTHON" -u "$@"
