# isognomon_ONT_mitogenomics
This repository contains code and information related to our paper on the Isognomon nucleus mitogenome

#code for Dreyer et al. The mitochondrial genome of Isognomon nucleus and mitophylogenomics of Pteriomorphia (Mollusca: Bivlavia: Autobranchia)
#this document will outline
- Data quality inspection
- Annotation and Assembly
- Phylogenetics (Alignment, Ambiguity assesment with GBlocks and Maximum-Likelihood inference with IQTREE)

#concact: Niklas Dreyer, Biodiversity Research Center, Academia Sinica, Taiwan + Natural History Museum of Denmark, University of Copenhagen, Denmark

---- DATA QUALITY ----
#inspect quality and read statistics for raw files
#install nanopack (contains all scripts: https://github.com/wdecoster/nanopack)

conda create --name nanopack
conda activate nanopack
conda update conda
pip install nanopack --upgrade

```
NanoPlot -t 15 --format <your_format> --fastq <your_name>.fastq --title <your_name> --maxlength 50000 --N50 --plots dot kde -o <your_output_dir>
```

---- QUALITY FILTERING AND ASSEMBLY ----

---- ANNOTATION ----

#get annotated mtgenome sequences from polished contig1 with MitoS. Submit contig to http://mitos.bioinf.uni-leipzig.de/
#inspect ORFs and, if necessary, manually curate reads to stop codon. 

---- PHYLOGENETICS ----

#download annotated sequences (.faa and .fas files) from MitoS.

#download Pteriomorphia protein coding genes from GenBank
#accession numbers as follows

MW143047.1
MT419374.1
KU589290.1
KC153059.1
MT916741.1
MT916743.1
MT916745.1
KY270857.1
MK948426.1
EU715252.2
MH051332.1
MZ497416.1
EU266073.1
FJ841967.1
EU672833.1
NC_023384.1
MT985154.1
FJ595959.1
MT991018.1
KP100300.1
MF407676.1
FJ890850.1
AY823625.1
GU936626.1
HM015199.1
KC768038.1
MG766134.1
KM655841.1
JQ970425.1
MH330333.1
HM467838.1
KX669229.1
GQ452847.1
KU552127.1
KT992045.1
KU310914.1
KU310918.1
KU310920.1
MT026713.1
KP205428.1

#remove 'lcl|' from the download file with sed (model; sed 's/find/replace/' <your_file>)

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




















