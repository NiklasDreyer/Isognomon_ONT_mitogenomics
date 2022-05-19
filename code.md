## Code for "The mitochondrial genome of *Isognomon nucleus* and mitophylogenomics of Pteriomorphia (Mollusca: Bivlavia: Autobranchia)" by Dreyer et al. 

## This document will outline how to perform:

- Data quality inspection
- Annotation and Assembly
- Phylogenetics (Alignment, Ambiguity assesment with GBlocks and Maximum-Likelihood inference with IQTREE)


### INSPECT NANOPORE DATA QUALITY WITH NANOPACK

#inspect quality and read statistics for raw files with Nanopack (https://github.com/wdecoster/nanopack)
#create environment and activate

```
conda create --name nanopack
conda activate nanopack
conda update conda
pip install nanopack --upgrade
```

#run Nanoplot
```
NanoPlot -t 15 --format jpg --fastq <your_name>.fastq --title <your_name> --maxlength 50000 --N50 --plots dot kde -o <your_output_dir>
```

### QUALITY FILTERING, TRIMMING, ASSEMBLY AND ANNOTATION OF ILLUMINA READS
#make conda environment and activate

```
conda create --name trimassemble
conda activate trimassemble
conda install -c bioconda/label/cf201901 trimmomatic
conda install -c bioconda getorganelle
conda install -c bioconda/label/cf201901 seqtk
conda install -c bioconda/label/cf201901 fastqc
conda install -c bioconda/label/cf201901 bbmap ## install for removing phiX spike and human mt sequences 
```

#extract files with gunzip
```
gunzip <your_forward_reads>.fq.gz
gunzip <your_reverse_reads>.fq.gz
```

#run fastqc on extracted files, inspect html file
```
fastqc <your_forward_reads>.fq
fastqc <your_reverse_reads>.fq
```

#subsample with seqtk (optional, not necessarily recommended for mitochondrial genome assembly)
```
seqtk sample -s100 <your_file_R1>.fastq 100000 > <your_outputfile_sub_R1>.fq
seqtk sample -s100 <your_file_R2>.fastq 100000 > <your_outputfile_sub_R2>.fq
```

#trim reads
```
trimmomatic PE <your_forward_reads>.fq <your_reverse_reads>.fq output_forward_paired.fq.gz output_forward_unpaired.fq.gz output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36  

#note that this command runs trimmomatic in the conda environment. 
#If run outside the environment, run as java -jar trimmomatic-0.39.jar <commands ... >
```

#optional removal on phi-X and Homo sapiens contaminants and kmer filtering
#navigate to working directory and make text file, insert phiX sequences from the BBDuk folder and the mitogenome sequence of Homo sapiens from NCBI

```
cd /PATH/TO/YOUR/DIRECTORY
cat homo_phiadapter.txt

>PhiX_read1_adapter
AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTGAAA
>PhiX_read2_adapter
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAAA
>NC_012920.1_Homo_sapiens
GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGG # is continued (.........)
```

#run BBDuk on forward reads
```
bbduk.sh in=output_forward_paired.fq.gz out=unmatched.fq outm=matched.fq ref=phix_adapters.fa k=31 hdist=1 stats=stats.txt
```

#run BBDuk on forward reads
```
bbduk.sh in=output_reverse_paired.fq.gz out=unmatched_reversed.fq outm=matched_reversed.fq ref=phix_adapters.fa k=31 hdist=1 stats=stats_reversed.txt
```

#inspect fastqc post-trimming
```
fastqc output_forward_paired.fq.gz
fastqc output_reverse_paired.fq.gz
```

#assemble mitochondrial reads with GetOrganelle (https://github.com/Kinggerm/GetOrganelle)
```
conda create -n getorganelle python=3.6 
conda activate getorganelle 
conda install -c bioconda getorganelle
```

#download reference seeds from NCBI GenBank (we downloaded and concatenated coding sequences of Pteriomorphia mitogenomes)
#run get organelle with this expanded seed database

```
get_organelle_from_reads.py -1 output_forward_paired.fq.gz -2 output_reverse_paired.fq.gz -s <your_seed>.fasta -t 15 -F animal_mt -o <your_outdir> -R 15 
```

#submit output file to MitoS server (http://mitos.bioinf.uni-leipzig.de/index.py)

### PHYLOGENY ESTIMATION

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
#export alignment files to separate folder remove everything after "_prot_" with sed before trimming alignment.

```
sed 's/__prot_.*//' <your_input_file> > <your_output_file>
```

#We recommend downloading replicate mitogenome reads when/if possible from NCBI and inspect whether these are identical across samples after alignment.
#These can easily be removed from the alignment files with this script
#inearize the FASTA with awk, pipe to grep to filter for items of interest named in patterns.txt, then pipe to tr to delinearize.

```
awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' atp6_fix.fasta | grep -Ff pattern.txt - | tr "\t" "\n" > atp6_fix_pattern.fasta
```

#import alignments to GUI visualizing software (We used GeneiousPrime, MEGA can be downloaded for free here: https://www.megasoftware.net/)

#concatenate (We used GeneiousPrime)

#asses ambiguities in trimmed alignment with Gblocks (https://github.com/atmaivancevic/Gblocks ; http://molevol.cmima.csic.es/castresana/Gblocks/Gblocks_documentation.html; http://www.phylogeny.fr/one_task.cgi?task_type=gblocks&tab_index=2)

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
iqtree2 -h
```

#estimate phylogeny with ultrafast boostrap and SH-aLRT replicates, allow substitution model partitions to merge if this increases fit

```
iqtree2 -s <your_file>.fasta -p <your_partition>.nex -m MFP+MERGE --seqtype AA --prefix <your_prefix> -B 10000 -bnni --alrt 8000 -nt AUTO 
```

#visualize tree with ETE3 python toolkit with the Newick file (and optionally the alignment) as input file(s)
http://etetoolkit.org/treeview

#alternatively visualize tree with FigTree/iTOL (https://itol.embl.de) tree viewers (We used GeneiousPrime)
#navigate to the terminal and install figtree within environment, open FigTree from the terminal.
```
conda install -c bioconda figtree
figtree
```
