#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# Standard configuration
# ============================================================

MAKE_JOBS=8
MAX_EVENTS=100
BACKEND="cuda"

INPUT_FILE=""
RUN_TWICE=false

RUN_PYTHON=false
PYTHON_SCRIPT="scripts/autograph.py"
PYTHON_ARGS=()

# Opções comuns do alpaka
FROM_HITS=false
IS_PHASE2=false
IS_COLLIDER_ML=false
VALIDATION=false
TRACKER_TYPE=""

# Extra arguments for each execution
#
# Each --first-extra / --second-extra takes ONLY ONE argument
# passed directly to alpaka.
FIRST_EXTRA_ARGS=()
SECOND_EXTRA_ARGS=()

# ============================================================
# Auxiliary functions
# ============================================================

usage() {
    cat <<EOF
Usage:
    $0 [opções]

Description:
    Optionally executes a python script, recompiles alpaka
    executes code one or two times with the specified options.

General options:
    --input FILE
        Input txt file

    --max-events N
        Maximum number of events. Default: ${MAX_EVENTS}

    --make-jobs N
        NNumber of jobs used by make. Default: ${MAKE_JOBS}

    --no-cuda
        Do not passes --cuda to alpaka.

Reconstruction options:
    --from-hits
        Adds --fromHits.

    --phase2
        Adds --isPhase2.

    --collider-ml
        Adds --isColliderML.

    --validation
        Adds --validation.

    --tracker-type TYPE
        Adds --trackerType TYPE. Possible values: PixelOnly,
        PixelPlusShortStrips or AllTracker.

Script Python:
    --autograph
        Executes:
            python3 scripts/autograph.py --input FILE

    --python-script SCRIPT
        Executes a Python script specified by the user.

    --python-arg ARG
        Aditional argument for the Python script.
        Can be used multiple times.

Two executions:
    --run-twice
        Executes alpaka twice.

    --first-extra ARG
        Adds an argument only to the first execution.
        Can be used multiple times.

    --second-extra ARG
        Adds an argument only to the second execution.
        Can be used multiple times.

    --help
        Shows this message.

Examples:

1) Simple execution:

    $0 --input input.root --from-hits --phase2 --collider-ml --tracker-type PixelOnly --validation

2) Running autograph first:

    $0 \\
        --input input.root \\
        --autograph \\
        --from-hits \\
        --phase2 \\
        --collider-ml \\
        --validation

3) Two executions with different extra options:

    $0 \\
        --input input.root \\
        --autograph \\
        --from-hits \\
        --phase2 \\
        --run-twice \\
        --first-extra --someOption \\
        --second-extra --validation

4) Arbitrary python script:

    $0 \\
        --input input.root \\
        --python-script scripts/my_script.py \\
        --python-arg --foo \\
        --python-arg bar

EOF
}

die() {
    echo "Erro: $*" >&2
    exit 1
}

# ============================================================
# Arguments parsing
# ============================================================

while [[ $# -gt 0 ]]; do
    case "$1" in

        --input)
            [[ $# -ge 2 ]] || die "--input needs an argument."
            INPUT_FILE="$2"
            shift 2
            ;;

        --max-events)
            [[ $# -ge 2 ]] || die "--max-events needs an argument."
            MAX_EVENTS="$2"
            shift 2
            ;;

        --make-jobs)
            [[ $# -ge 2 ]] || die "--make-jobs needs an argument."
            MAKE_JOBS="$2"
            shift 2
            ;;

        --backend)
            [[ $# -ge 2 ]] || die "--backend needs an argument."
            BACKEND="$2"
            shift 2
            ;;

        --from-hits)
            FROM_HITS=true
            shift
            ;;

        --phase2)
            IS_PHASE2=true
            shift
            ;;

        --collider-ml)
            IS_COLLIDER_ML=true
            shift
            ;;

        --validation)
            VALIDATION=true
            shift
            ;;

        --tracker-type)
            [[ $# -ge 2 ]] || die "--tracker-type needs an argument."
            TRACKER_TYPE="$2"
            shift 2
            ;;

        --autograph)
            RUN_PYTHON=true
            PYTHON_SCRIPT="scripts/autograph.py"
            shift
            ;;

        --python-script)
            [[ $# -ge 2 ]] || die "--python-script needs an argument."
            RUN_PYTHON=true
            PYTHON_SCRIPT="$2"
            shift 2
            ;;

        --python-arg)
            [[ $# -ge 2 ]] || die "--python-arg needs an argument."
            PYTHON_ARGS+=("$2")
            shift 2
            ;;

        --run-twice)
            RUN_TWICE=true
            shift
            ;;

        --first-extra)
            [[ $# -ge 2 ]] || die "--first-extra needs an argument."
            FIRST_EXTRA_ARGS+=("$2")
            shift 2
            ;;

        --second-extra)
            [[ $# -ge 2 ]] || die "--second-extra needs an argument."
            SECOND_EXTRA_ARGS+=("$2")
            shift 2
            ;;

        --help|-h)
            usage
            exit 0
            ;;

        *)
            die "Unknown option: $1"
            ;;

    esac
done

# ============================================================
# Validations
# ============================================================

if [[ "$RUN_PYTHON" == true && -z "$INPUT_FILE" ]]; then
    die "--input is required when a Python script is executed."
fi

if [[ "$RUN_PYTHON" == true && ! -f "$PYTHON_SCRIPT" ]]; then
    die "Python script not found: $PYTHON_SCRIPT"
fi

if ! [[ "$MAX_EVENTS" =~ ^[0-9]+$ ]]; then
    die "--max-events must be a positive integer."
fi

if ! [[ "$MAKE_JOBS" =~ ^[0-9]+$ ]] || [[ "$MAKE_JOBS" -eq 0 ]]; then
    die "--make-jobs must be a positive integer."
fi

case "$BACKEND" in
    cuda|serial|rocm)
        ;;
    *)
        die "Invalid backend: $BACKEND. Use: cuda, serial or rocm."
        ;;
esac

if [[ -n "$TRACKER_TYPE" ]]; then
    case "$TRACKER_TYPE" in
        PixelOnly|PixelPlusShortStrips|AllTracker)
            ;;
        *)
            die "Invalid tracker: $TRACKER_TYPE. Use: PixelOnly, PixelPlusShortStrips or AllTracker."
            ;;
    esac
fi

# ============================================================
# Python script
# ============================================================

if [[ "$RUN_PYTHON" == true ]]; then
    echo
    echo "============================================================"
    echo "Executing python script: $PYTHON_SCRIPT"
    echo "============================================================"

    PYTHON_CMD=(
        python3
        "$PYTHON_SCRIPT"
        --input
        "$INPUT_FILE"
    )

    if [[ ${#PYTHON_ARGS[@]} -gt 0 ]]; then
        PYTHON_CMD+=("${PYTHON_ARGS[@]}")
    fi

    echo "+ ${PYTHON_CMD[*]}"
    "${PYTHON_CMD[@]}"
fi

# ============================================================
# Compilation
# ============================================================

echo
echo "============================================================"
echo "Compiling code with make -j $MAKE_JOBS"
echo "============================================================"

MAKE_CMD=(
    make
    -j
    "$MAKE_JOBS"
    alpaka
)

echo "+ ${MAKE_CMD[*]}"
"${MAKE_CMD[@]}"

# ============================================================
# Assembly of common reconstruction arguments
# ============================================================

ALPAKA_BASE_ARGS=(
    "--$BACKEND"
    --maxEvents "$MAX_EVENTS"
)

if [[ "$FROM_HITS" == true ]]; then
    ALPAKA_BASE_ARGS+=(--fromHits)
fi

if [[ "$IS_PHASE2" == true ]]; then
    ALPAKA_BASE_ARGS+=(--isPhase2)
fi

if [[ "$IS_COLLIDER_ML" == true ]]; then
    ALPAKA_BASE_ARGS+=(--isColliderML)
fi

if [[ "$VALIDATION" == true ]]; then
    ALPAKA_BASE_ARGS+=(--validation)
fi

if [[ -n "$TRACKER_TYPE" ]]; then
    ALPAKA_BASE_ARGS+=(--trackerType "$TRACKER_TYPE")
fi

# ============================================================
# First execution
# ============================================================

echo
echo "============================================================"
echo "Executing reconstruction - execution 1"
echo "============================================================"

FIRST_CMD=(
    ./alpaka
    "${ALPAKA_BASE_ARGS[@]}"
    "${FIRST_EXTRA_ARGS[@]}"
)

echo "+ ${FIRST_CMD[*]}"
"${FIRST_CMD[@]}"

# ============================================================
# Second execution, if requested
# ============================================================

if [[ "$RUN_TWICE" == true ]]; then

    echo
    echo "============================================================"
    echo "Executing reconstruction - execution 2"
    echo "============================================================"

    SECOND_CMD=(
        ./alpaka
        "${ALPAKA_BASE_ARGS[@]}"
        "${SECOND_EXTRA_ARGS[@]}"
    )

    echo "+ ${SECOND_CMD[*]}"
    "${SECOND_CMD[@]}"
fi

echo
echo "============================================================"
echo "Reconstruction completed successfully."
echo "============================================================"
