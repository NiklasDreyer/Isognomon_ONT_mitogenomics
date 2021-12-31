## Code for "The mitochondrial genome of *Isognomon nucleus* and mitophylogenomics of Pteriomorphia (Mollusca: Bivlavia: Autobranchia)" by Dreyer et al. 

## This document will outline how to perform:

- Data quality inspection
- Annotation and Assembly
- Phylogenetics (Alignment, Ambiguity assesment with GBlocks and Maximum-Likelihood inference with IQTREE)


### DATA QUALITY 

#inspect quality and read statistics for raw files

#install nanopack (https://github.com/wdecoster/nanopack)

```
conda create --name nanopack
conda activate nanopack
conda update conda
pip install nanopack --upgrade
NanoPlot -t 15 --format jpg --fastq <your_name>.fastq --title <your_name> --maxlength 50000 --N50 --plots dot kde -o <your_output_dir>
```

### QUALITY FILTERING AND ASSEMBLY

###  ANNOTATION 

#get annotated mtgenome sequences from polished contig with MitoS. Submit contig to http://mitos.bioinf.uni-leipzig.de/

#inspect ORFs and, if necessary, manually curate reads to stop codon. 

### PHYLOGENETICS 

#download annotated sequences (.faa and .fas files) from MitoS.

#make .txt file with NCBI accession numbers 

```
cat > <your_text_accesion_list>.txt
```

#download Pteriomorphia protein coding genes from GenBank with Batch Entrez (https://www.ncbi.nlm.nih.gov/sites/batchentrez)

#remove 'lcl|' from the download file with sed ('s/find/replace/')

```
sed 's/lcl|//' <your_file>
```

#organize sequences in folders for each gene (we used GeneiousPrime)

#align the sequences with MAFFT (https://github.com/GSLBiotech/mafft ; https://mafft.cbrc.jp/alignment/software/manual/manual.html#lbAI)

#automatic selction of appropriate alignment strategy with either L-INS-i, FFT-NS-i or FFT-NS-2

#use BLOSUM62 protin alignment matrix (--bl 62) 

```
conda create --name mafft
conda activate mafft
conda update conda
conda install -c bioconda/label/cf201901 mafft

mafft --auto --amino --maxiterate 100 --op 1.53 --ep 0.123 --lop -2.00 --lep 0.1 --lexp -0.1 --LOP -6.00 --LEXP 00 --bl 62 --jtt BLOSUM62 --tm BLOSUM62 --aamatrix BLOSUM62 --inputorder input_file > output_file
```
#export alignment files to separate folder remove everything after "_prot..." with sed before trimming alignment

```
sed 's/__prot_.*//' <your_input_file> > <your_output_file>
```

#import and trim alignments to shortest sequence (We used GeneiousPrime)

#concatenate (We used GeneiousPrime)

#asses ambiguities in trimmed alignment with Gblocks (https://github.com/atmaivancevic/Gblocks ; http://molevol.cmima.csic.es/castresana/Gblocks/Gblocks_documentation.html)

```
conda create --name gblocks 
conda activate gblocks 
conda update conda
conda install -c bioconda/label/cf201901 gblocks
```

#add "P1" to beginning of each ">" line with sed to allow Gblocks to run ambiguity assessment on protein sequences in command line (if you want to avoid this, run the alignment with prefered options through the webserver; http://molevol.cmima.csic.es/castresana/Gblocks_server.html)

```
sed 's/^/P1/' <your_input_file>
```

#manually add space after each sequence ID line and run Gblocks

```
gblocks <your_input_file> -t=p -b1=24 -b2=39 -b3=4 -b4=10 -b5=n -b6=y
```

#infer maximum likelihood phylogenetic trees with IQTREE (https://github.com/Cibiv/IQ-TREE ; http://www.iqtree.org/)

```
conda create --name iqtree
conda activate iqtree
conda update conda
conda install -c bioconda/label/cf201901 iqtree
```

#check installation

```
iqtree
```

#run maximum likelihood inference 

```
iqtree -s <your_input_file> --seqtype AA --prefix <your_prefix> --seed 321 --ninit 300 -B 30000 -bnni --alrt 2000 -m MFP -nt AUTO -mredo -redo
```

#visualize tree with ETE3 python toolkit with the Newick file (and optionally the alignment) as input file(s)
http://etetoolkit.org/treeview

#alternatively visualize tree with FigTree/iTOL (https://itol.embl.de) tree viewers (We used GeneiousPrime)
```
conda install -c bioconda figtree
figtree
```
