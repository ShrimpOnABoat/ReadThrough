# ReadThrough
A set of scripts identifying the differential readthrough expression in a ribosome profiling experiment

Work I did during my PhD for this article: 
Del Toro N, Lessard F, Bouchard J, Mobasheri N, Guillon J, Igelmann S, Tardif S, Buffard T, Bourdeau V, Brakier-Gingras L, Ferbeyre G. Cellular senescence limits translational readthrough. Biol Open. 2021 Dec 1;10(12):bio058688. doi: 10.1242/bio.058688. Epub 2021 Dec 2. PMID: 34676390; PMCID: PMC8649927.

Using the results from this ribosome profiling experiment : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE45833

We wanted to know how cellular senescence changes the expression profile of cells regarding translation readthrough, where ribosomes would keep translating through a STOP codon until the next one.

- riboseq : bash script to extract data from wig files
- wig files : tables associating the number of reads with the position for each chromosome (only the first 50 lines since these files are huge)
- UTR3, UTR5 files : for each gene, the first and last position of the UnTranslated Region, for the 3' or 5' extremity
- Short_UTR3.py : script computing, for each gene, the positions of the 3' and 5' UTR, the chromosome number, the strand, and mean and sum values from the riboseq
- KnownGene.txt : reference file from Human Genome 19 (GrCH37) with names, chromosomes, directions and positions
- kgXref.txt : reference file for each gene, with common names, NCBI refs and description
- ShortUTR_score.py : script creating a file with data for UTR only
- UTR3_RB_Norm : list of all the 3' UTR, their beginning and ending position, chromosome and gene names, as well as their direction
- UTR3_Shrt_codon : list of all the 3' UTR with the first STOP codon and the following 4 nucleotides
- CDS_file.py : script creating a list of all the coding sequences of the genome, and computing mean and sum values for each CDS
- CDS_hg19 : reference file for the coding sequences in the human genome
- CDS_KG : CDS beginning and ending positions, chromosome and gene names, and direction
- UTR_FirstNt.py : script finding matching STOP codons and following sequences for each 3' UTR 
- Comparaison readthrough.xlsx : result file
