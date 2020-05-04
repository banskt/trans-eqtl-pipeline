## How to run the cis-eGene characterization analysis
1. Create the list of gene ids with biotypes.
2. Map the biotype of the cis-mediating genes.
3. Get list of unique cis-eGenes from the GTEx results and also map their biotypes.
4. Calculate the fraction of each biotype.
5. Calculate the fraction of transcription factors
```
./create_gene_type_list.sh ../gencode_annotation_enrichment/gencode.v26.annotation.gtf.gz gencode.v26.gene_id_biotype_list.txt
mkdir gtex_v8
cd gtex_v8
wget https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL.tar
tar -xvf GTEx_Analysis_v8_eQTL.tar
../get_uniq_cis-eGenes.sh 
cut -f3 GTEx_Analysis_v8_uniq_cis-Genes.txt | sort | uniq -c | sort -nr > gtex_v8_cis-eGenes_counts.txt
cd ..
./map_biotype.sh CistransEQTL_targets.txt > cis-eGenes_all_tissues.txt
./calc_counts_biotypes.sh cis-eGenes_all_tissues.txt ../external/human_TF_annotation_gencode_v26.txt > cis-eGenes_uniq_all_tissues_biotypes.txt
```
