# mlst-phylogeny-pipeline

conda create -n raxmlng
conda install -c bioconda mafft -y
mafft --version
conda install -c bioconda gblocks -y
Gblocks
conda install bioconda::raxml-ng
conda install -c conda-forge parallel -y    
conda install -c conda-forge biopython

python mlst_ncbi_protein.py
ls mlst_extracted_proteins|wc
bash clean_name.sh  
python merge_all_mlst_protein.py


#https://biohpc.cornell.edu/doc/alignment_exercise3.html

  
mkdir -p aligned_gblock aln_merged_marker

ls merged_marker/*.faa | xargs -I {} echo "mafft --thread 1 --amino  --inputorder --quiet {}  > aln_{} ; Gblocks aln_{} -t=p -b5=h" > msa.batch

conda activate raxmlng

bash msa.batch 


#second way

ls merged_marker/*.faa | xargs -I {} bash -c '\
  base=$(basename {} .faa); \
  mafft --thread 1 --amino --inputorder --quiet {} > aligned_gblock/${base}_aln.fa; \
  Gblocks aligned_gblock/${base}_aln.fa -t=p -b5=h' \
> msa.batch

conda deactivate
mkdir MSAdir

mv aln_merged_marker/*gb MSAdir/

python concate_msa_cornell_order.py MSAdir

sed -i '' '/^>/ s/ /_/g' merged.fasta



#java -jar jmodeltest-2.1.10/jModelTest.jar -d merged.fasta -s 11 -f -g 4 -i -AIC -BIC

java -jar prottest-3.4.2/prottest-3.4.2.jar \
-i merged.fasta \
-all-distributions \
-F \
-AIC -BIC \
-threads 4


The best-fitting nucleotide substitution model (GTR+I+G) was selected using jModelTest v2.1.10 based on the BIC criterion. Phylogenetic analysis was performed using RAxML-NG under this model, with 100 bootstrap replicates to assess branch support.

conda activate raxmlng
raxml-ng --all \
--threads 6 \
--msa merged.fasta \
--model LG+I+G+F \
--prefix merged \
--bs-trees 1000 \
--seed 12345

raxml-ng --all \
--threads 6 \
--msa merged.fasta \
--model LG+G+F \
--prefix merged \
--bs-trees 1000 \
--seed 12345
