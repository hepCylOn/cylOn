#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# Configuração padrão
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

# Argumentos extras específicos de cada execução.
#
# Cada --first-extra / --second-extra representa UM argumento
# passado diretamente ao alpaka.
FIRST_EXTRA_ARGS=()
SECOND_EXTRA_ARGS=()

# ============================================================
# Funções auxiliares
# ============================================================

usage() {
    cat <<EOF
Uso:
    $0 [opções]

Descrição:
    Executa opcionalmente um script Python, recompila o alpaka
    e executa o programa uma ou duas vezes.

Opções gerais:
    --input FILE
        Arquivo de entrada.

    --max-events N
        Número máximo de eventos. Padrão: ${MAX_EVENTS}

    --make-jobs N
        Número de jobs usados pelo make. Padrão: ${MAKE_JOBS}

    --no-cuda
        Não passa --cuda para o alpaka.

Opções do alpaka:
    --from-hits
        Adiciona --fromHits.

    --phase2
        Adiciona --isPhase2.

    --collider-ml
        Adiciona --isColliderML.

    --validation
        Adiciona --validation.

Script Python:
    --autograph
        Executa:
            python3 scripts/autograph.py --input FILE

    --python-script SCRIPT
        Executa um script Python especificado pelo usuário.

    --python-arg ARG
        Argumento adicional para o script Python.
        Pode ser usado várias vezes.

Duas execuções:
    --run-twice
        Executa o alpaka duas vezes.

    --first-extra ARG
        Adiciona um argumento somente à primeira execução.
        Pode ser usado várias vezes.

    --second-extra ARG
        Adiciona um argumento somente à segunda execução.
        Pode ser usado várias vezes.

    --help
        Mostra esta mensagem.

Exemplos:

1) Execução simples:

    $0 --input input.root --from-hits --phase2 --collider-ml --validation

2) Executando o autograph antes:

    $0 \\
        --input input.root \\
        --autograph \\
        --from-hits \\
        --phase2 \\
        --collider-ml \\
        --validation

3) Duas execuções com opções extras diferentes:

    $0 \\
        --input input.root \\
        --autograph \\
        --from-hits \\
        --phase2 \\
        --run-twice \\
        --first-extra --someOption \\
        --second-extra --validation

4) Script Python arbitrário:

    $0 \\
        --input input.root \\
        --python-script scripts/meu_script.py \\
        --python-arg --foo \\
        --python-arg bar

EOF
}

die() {
    echo "Erro: $*" >&2
    exit 1
}

# ============================================================
# Parsing dos argumentos
# ============================================================

while [[ $# -gt 0 ]]; do
    case "$1" in

        --input)
            [[ $# -ge 2 ]] || die "--input requer um argumento."
            INPUT_FILE="$2"
            shift 2
            ;;

        --max-events)
            [[ $# -ge 2 ]] || die "--max-events requer um argumento."
            MAX_EVENTS="$2"
            shift 2
            ;;

        --make-jobs)
            [[ $# -ge 2 ]] || die "--make-jobs requer um argumento."
            MAKE_JOBS="$2"
            shift 2
            ;;

        --backend)
            [[ $# -ge 2 ]] || die "--backend requer um argumento."
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

        --autograph)
            RUN_PYTHON=true
            PYTHON_SCRIPT="scripts/autograph.py"
            shift
            ;;

        --python-script)
            [[ $# -ge 2 ]] || die "--python-script requer um argumento."
            RUN_PYTHON=true
            PYTHON_SCRIPT="$2"
            shift 2
            ;;

        --python-arg)
            [[ $# -ge 2 ]] || die "--python-arg requer um argumento."
            PYTHON_ARGS+=("$2")
            shift 2
            ;;

        --run-twice)
            RUN_TWICE=true
            shift
            ;;

        --first-extra)
            [[ $# -ge 2 ]] || die "--first-extra requer um argumento."
            FIRST_EXTRA_ARGS+=("$2")
            shift 2
            ;;

        --second-extra)
            [[ $# -ge 2 ]] || die "--second-extra requer um argumento."
            SECOND_EXTRA_ARGS+=("$2")
            shift 2
            ;;

        --help|-h)
            usage
            exit 0
            ;;

        *)
            die "Opção desconhecida: $1"
            ;;

    esac
done

# ============================================================
# Validações
# ============================================================

if [[ "$RUN_PYTHON" == true && -z "$INPUT_FILE" ]]; then
    die "--input é obrigatório quando um script Python é executado."
fi

if [[ "$RUN_PYTHON" == true && ! -f "$PYTHON_SCRIPT" ]]; then
    die "Script Python não encontrado: $PYTHON_SCRIPT"
fi

if ! [[ "$MAX_EVENTS" =~ ^[0-9]+$ ]]; then
    die "--max-events deve ser um número inteiro."
fi

if ! [[ "$MAKE_JOBS" =~ ^[0-9]+$ ]] || [[ "$MAKE_JOBS" -eq 0 ]]; then
    die "--make-jobs deve ser um inteiro maior que zero."
fi

case "$BACKEND" in
    cuda|serial|rocm)
        ;;
    *)
        die "Backend inválido: $BACKEND. Use: cuda, serial ou rocm."
        ;;
esac

# ============================================================
# Script Python
# ============================================================

if [[ "$RUN_PYTHON" == true ]]; then
    echo
    echo "============================================================"
    echo "Executando script Python"
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
# Compilação
# ============================================================

echo
echo "============================================================"
echo "Compilando alpaka"
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
# Montagem dos argumentos comuns do alpaka
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

# ============================================================
# Primeira execução
# ============================================================

echo
echo "============================================================"
echo "Executando alpaka - execução 1"
echo "============================================================"

FIRST_CMD=(
    ./alpaka
    "${ALPAKA_BASE_ARGS[@]}"
    "${FIRST_EXTRA_ARGS[@]}"
)

echo "+ ${FIRST_CMD[*]}"
"${FIRST_CMD[@]}"

# ============================================================
# Segunda execução, se solicitada
# ============================================================

if [[ "$RUN_TWICE" == true ]]; then

    echo
    echo "============================================================"
    echo "Executando alpaka - execução 2"
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
echo "Execução concluída com sucesso."
echo "============================================================"
