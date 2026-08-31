#!/usr/bin/env bash
# setup_env.sh — cria o ambiente conda dedicado ao pipeline de GCFs
set -euo pipefail

# 'conda' é uma função de shell definida via 'conda init', que scripts
# não-interativos não carregam; localizamos o conda.sh diretamente.
CONDA_SH=""
for candidate in \
    "${HOME}/miniconda3/etc/profile.d/conda.sh" \
    "${HOME}/anaconda3/etc/profile.d/conda.sh" \
    "/opt/conda/etc/profile.d/conda.sh" \
    "/opt/ohpc/pub/miniconda3/etc/profile.d/conda.sh"
do
    if [ -f "${candidate}" ]; then
        CONDA_SH="${candidate}"
        break
    fi
done

if [ -z "${CONDA_SH}" ]; then
    echo "ERRO: não encontrei conda.sh nos locais padrão." >&2
    echo "Rode 'which conda' num shell interativo, pegue o prefixo (sem /bin/conda)," >&2
    echo "e ajuste a lista de candidatos acima ou exporte CONDA_SH=/caminho/etc/profile.d/conda.sh" >&2
    exit 1
fi
source "${CONDA_SH}"

conda create -n gcf_clustering \
    --override-channels \
    -c conda-forge -c bioconda \
    --strict-channel-priority \
    python=3.11 \
    pandas numpy scipy scikit-learn tqdm \
    biopython \
    hmmer \
    -y

echo "Ambiente 'gcf_clustering' criado. Ative com: conda activate gcf_clustering"