# 1. Preparar IDs e baixar
awk 'NR > 1 {print $1}' metadata/gbk_hp_genomes.tsv > acessos.txt
datasets download genome accession --inputfile acessos.txt --include genome

# 2. Descompactar
unzip -q ncbi_dataset.zip

# 3. Organizar, compactar e mostrar progresso com verificação de existência
mkdir -p fastas

total=$(find ncbi_dataset/data -name "*.fna" | wc -l)
atual=0

find ncbi_dataset/data -name "*.fna" | while read f; do
    atual=$((atual + 1))
    acc=$(echo "$f" | cut -d'/' -f3)
    arquivo_final="fastas/${acc}.fasta.gz"
    
    # Verifica se o arquivo já existe na pasta de destino
    if [ -f "$arquivo_final" ]; then
        printf "\r[%d/%d] Pulando (já existe): %s          " "$atual" "$total" "$acc"
    else
        printf "\r[%d/%d] Compactando: %s...               " "$atual" "$total" "$acc"
        gzip -c "$f" > "$arquivo_final"
    fi
done

echo -e "\n\nProcesso concluído! Os arquivos estão em 'fastas/'"

# 4. Limpeza
rm -rf ncbi_dataset ncbi_dataset.zip acessos.txt
