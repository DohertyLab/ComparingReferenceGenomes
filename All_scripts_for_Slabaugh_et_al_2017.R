# Custom scripts for Slabaugh et al. 2017, "Analysis of differential gene expression and alternative splicing is significantly influenced by choice of mapping genome"
# Note: please do not source this code directly as many of the commands are not wrapped in functions

#####################################################################################
#   Parse and re-format GFF files to match downstream analysis input requirements   #
#####################################################################################

# Driver function= GFFtoGTF

# This function is used to parse a GFF file (Assumes no header by default)
# Columns: 1 - Chromosome, 2 - Genome Version, 3 - Feature Type ("chromosome","gene","mRNA","CDS","Protein", etc.), 4 - Start BP, 5 - End BP, 7 - Strand (+/-), 9 - Name 
ParseGFF<-function(filenm, sep="\t", header=F, ...)
{
  library(stringr)
  # load the file
  print(paste(" Reading GFF File: ",filenm))
  tgff<-read.csv(file = filenm, sep = sep, header = header, stringsAsFactors=F, ...)
  
  # construct an array
  # [[Chromosome]]
  #   [[FeatureType]]
  #       [Matrix]
  ret<-tapply(1:nrow(tgff), INDEX = tgff[,1], function(x) PGFF1(tgff[x,]))
  
  # convert chromosome names
  rnm<-names(ret)
  rnm<-gsub("Chr","",rnm)
  names(ret)<-rnm
  
  return(ret)
}

PGFF1<-function(chgff)
{
  chret<-tapply(1:nrow(chgff), INDEX = chgff[,3], function(y) chgff[y,])
}

# function to convert GFF to GTF in correct format
# feat=feature
# outfile=output file name of gtf file ("/path/to/filename.gtf")
GFFtoGTF<-function(GFFfile, sep="\t", outfile)
{
  require(stringr)
  GFFlist<-ParseGFF(filenm = GFFfile, sep=sep)
  GFFlist<-lapply(GFFlist, function(chr) lapply(chr, function(feat) DriveSplit(feat)))
  
  final<-deParseGFF(GFFin = GFFlist, filename = outfile)
  return(final)
}

# The splitV9 function splits up attributes in V9
# instr=input string
# kvsep=key value seperator
# type=exon or CDS
splitV9<-function(instr, type, sep2=";", kvsep="=")
{
  original<-instr
  instr<-str_replace_all(string = instr, pattern = "=", replacement = " ")
  instr<-gsub(x = instr, pattern = "(;|:)", replacement = "; ", fixed = T)
  # note: "pattern" for idstr and transid below is shown for the Nipponbare MSU genome, change "pattern" to match MH63 or ZS97 genomes accordingly
  idstr<-str_extract(string = instr, pattern = "ID LOC_Os[[:alpha:]]{0,2}[[:digit:]]{2}.{1}[[:digit:]]{5}[[^;]]*;")
  
  if(!is.na(idstr))
  {
    idstr<-str_split(string = idstr, pattern = "\\.")[[1]][1]
    idstr<-gsub(pattern = "ID", replacement = "gene_id", x = idstr)
  }
  
  retstr<-idstr
  if(type=="exon" | type=="CDS")
  {
    if (type=="exon")
    {
      nmstr<-str_extract(string = instr, pattern = "exon_[[:digit:]]{1,2}")
      nmstr<-gsub(pattern = "exon\\_", replacement = "", x = nmstr)
    }
    
    else
    {
      nmstr<-str_extract(string = instr, pattern = "cds_[[:digit:]]{1,2}")
      nmstr<-gsub(pattern = "cds\\_", replacement = "", x = nmstr)
    }
    
    transid<-str_extract(string = instr, pattern = "Parent LOC_Os[[:alpha:]]{0,2}[[:digit:]]{2}.{1}[[:digit:]]{5}\\.[[:digit:]]{1,2}")
    transid<-gsub(pattern = "Parent", replacement = "transcript_id", x = transid)
    retstr<-paste(retstr, "; ", transid, "; exon_number ", nmstr, sep = "")
  }
  
  return(retstr)
}

DriveSplit<-function(inDf)
{
  inDf$V9<-sapply(inDf$V9,function(x) splitV9(x,type=inDf$V3[1]))
  return(inDf)
}

# This function creates a flat file from GFF object
# input: GFF object (GFFin), created by ParseGFF
# input: filename (fname)
# output: flat file

deParseGFF<-function(GFFin, filename)
{
  GFFin<-lapply(GFFin, function(x) do.call(rbind, x))
  GFFin<-lapply(GFFin, function(x) x[order(x$V4),])
  retmat<-do.call(rbind, GFFin)
  write.table(x = retmat, file = filename, quote = F, sep = "\t", row.names = F, col.names = F, append = F)
  return(TRUE)
}

# remove NA's from Nipponbare_MSU gtf
FinalGTFMSU<-read.table(file = "~/Japonica_genome_files/MSU/FinalGTFMSU.gtf", sep = "\t", stringsAsFactors = F)
i<-which(FinalGTFMSU$V1=="ChrUn")
j<-which(FinalGTFMSU$V1=="ChrSy")
k<-c(i, j)
FinalGTFMSU2<-FinalGTFMSU[-k,]
write.table(x = FinalGTFMSU2, file="~/Japonica_genome_files/MSU/FinalGTFMSU_no_na.gtf", quote = F, sep = "\t", col.names = F, row.names = F)


##########################################################################
#   TopHat v2.1.0 Alignment of IR64 fastq files to three Oryza genomes   #
##########################################################################

# build bowtie index for each genome
bowtie2-build ~/Indica_genome_files/MH63_1.seq ~/MH_index/MH_index

# for MH63 genome
tophat --num-threads 12 --output-dir S10 --GTF ~/Indica_genome_files/FinalGTF_MH.gtf --library-type fr-firststrand ~/MH_index/MH_index /storage/data_1/jigar/IRRI_RNAseq/fastq_files/1-Dawn_C.1.fastq

# for ZS97 genome
tophat --num-threads 12 --output-dir ~/1-Dawn_C.1 --GTF ~/Indica_genome_files/FinalGTF_ZS.gtf --library-type fr-firststrand ~/ZS_index/ZS_index /storage/data_1/jigar/IRRI_RNAseq/fastq_files/Panicle_50_Control/1-Dawn_C.1.fastq 

# for MSU genome
tophat --num-threads 12 --output-dir ~/1-Dawn_C.1 --GTF ~/Japonica_genome_files/MSU/FinalGTFMSU_no_na.gtf --library-type fr-firststrand ~/MSU_MSU_index/MSU_MSU_index /storage/data_1/jigar/IRRI_RNAseq/fastq_files/Panicle_50_Control/1-Dawn_C.1.fastq

# use function "tophat_indica" to call tophat_indica_function
# !/bin/bash

for i in $(ls *.fastq)

do 

tophat --num-threads 12 --output-dir ~/${i%.fastq} --GTF ~/Indica_genome_files/FinalGTF_MH.gtf --library-type fr-firststrand ~/MH_index/MH_index ${i}

done

# change settings of function file
chmod 755 tophat_indica_function
# make bash profile:
nano ~/.bash_profile
# add to bash profile
PATH=$PATH: ~/Bash
export PATH
# save and exit nano
source ~/.bash_profile

tophat_indica_function

###############################################
#                 DESeq2                      #
###############################################
# use DESeq2 vingette and http://www.bioconductor.org/help/workflows/rnaseqGene/ as resources

library("GenomicAlignments")
library("DESeq2")
library("Rsamtools")
library("GenomicFeatures")
library("BiocParallel")

# generate SummarizedExperiment object to feed into DESeq2 package
# generate SampleTable
sampleTable = data.frame(
  row.names = c( "ZS_Dawn_C.1", "ZS_Dawn_C.2", "ZS_Dawn_C.3",
                 "ZS_Dawn_C.4", "ZS_Dusk_C.1", "ZS_Dusk_C.2", "ZS_Dusk_C.3", "ZS_Dusk_C.4" ),
  condition = c(rep("Dawn",4), rep("Dusk", 4)),
  samples = c( "ZS_Dawn_C.1", "ZS_Dawn_C.2", "ZS_Dawn_C.3",
               "ZS_Dawn_C.4", "ZS_Dusk_C.1", "ZS_Dusk_C.2", "ZS_Dusk_C.3", "ZS_Dusk_C.4"  ))

# generate path to bam files
filenames <- file.path("~/ZS_BAM", paste0(sampleTable$samples, ".bam"))
file.exists(filenames)

# provides an R interface to BAM files
bamfiles <- BamFileList(filenames)

# load gtf file
gtffile <- file.path("/storage/data_1/erins/Indica_genome_files","FinalGTF_ZS.gtf")
txdb <- makeTxDbFromGFF(gtffile, format="gtf", circ_seqs=character())

# produces a GRangesList of all the exons grouped by gene (Lawrence et al. 2013). Each element of the list is a GRanges object of the exons for a gene.
ebg <- exonsBy(txdb, by="gene")

# create the SummarizedExperiment object with counts
se <- summarizeOverlaps(features=ebg, reads=bamfiles,
                        mode="Union",
                        singleEnd=TRUE,
                        ignore.strand=TRUE)

colData(se) <- DataFrame(sampleTable)

# check the millions of fragments that uniquely aligned to the genes (the second argument of round tells how many decimal points to keep)
round( colSums(assay(se)) / 1e6, 1 )

# DESeq2
dds <- DESeqDataSet(se, design = ~ condition)

# remove the rows that have no or nearly no information about the amount of gene expression (additional filtering is performed later)
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)

# differential expression analysis
dds <- DESeq(dds)

# build results table
res <- results(dds)

# adjust FDR and LFC threshold
res.05LFC1 <- results(dds, lfcThreshold=1, alpha=.05)
sum(res.05$padj < 0.05, na.rm=TRUE)

# subset results table
res.05Sig <- subset(res.05, padj < 0.05)

# extract gene names
DEGenes.05<- row.names(res.05Sig)

# ggplot
geneCounts <- plotCounts(ddsMSU, gene="LOC_Os03g61160", intgroup=c("condition"), returnData=TRUE)
ggplot(geneCounts, aes(x=condition, y=count, color=condition)) +
  # scale_y_log10(breaks=c(0,10,100,1000), labels=(c(0,10,100,1000))) + 
  geom_point(position=position_jitter(width=.1,height=0), size=5) + 
  theme_bw() + #for white background panel
  theme(text = element_text(size=20), axis.text.x = element_text(angle=90, hjust=1)) + 
  coord_trans(y = "log10") +
  scale_color_manual(values=c("red", "blue"))

# MA plot
plotMA(res.05, ylim=c(-5,5))

# match DEG syntenic orthologs from each genome using the MSU nomeclauture as a referece from (see MCScanX script below)
# get list of DEGs only called differentially expressed for a specific genome
DEuniqueMH<- setdiff(x = DEMHorthologs, y = DEZSorthologs)
DEuniqueMH<- setdiff(x = DEuniqueMH, y = DEMSUorthologs)
DEuniqueMHnames<- match(x = DEuniqueMH, table = finalsynteny[,3])
DEuniqueMHnames<- finalsynteny[DEuniqueMHnames,]
DEuniqueMHnames<- DEuniqueMHnames[,1]

# to query counts for gene loci
counts(ddsMH)["MH03t0738500",]

# to query counts for gene loci with low counts
ddsprefilterZS<- DESeqDataSet(seZS, design = ~ condition)
counts(ddsprefilterZS)["ZS12t0519200",]

# make Venn diagram from syntenic lists
DEGeneListSynteny<- list(MH63=DEMHorthologs, ZS97=DEZSorthologs, MSU=DEMSUorthologs)
venn.diagram(DEGeneListSynteny, filename = "./vennDE.tif", imagetype = "tiff", fill = c("midnightblue","yellow","lightseagreen"), cex = 1.5, fontface = "bold", cat.cex = 1.5, cat.fontface = "bold", euler.d=TRUE, scaled=TRUE, na="remove" )

# Get gene_ids of genes only identified as DEGs in one genome
# SetLst is a named list of elements to compare
Plot3WayVenn<-function(SetLst)
{
  #this package requires "Vennerable"
  require(Vennerable)
  
  u<-unique(unlist(SetLst))
  c1<-rep(0,length(u))
  names(c1)<-u
  c2<-c1
  c3<-c1
  c1[SetLst[[1]]]<-1
  c2[SetLst[[2]]]<-1
  c3[SetLst[[3]]]<-1
  tlabs<-paste(c1,c2,c3,sep="")
  ttab<-table(tlabs)
  names(tlabs)<-names(c1)
  AllSets<-tapply(names(tlabs),tlabs,function(x) x)
  return(AllSets)
  ret<-Venn(Weight = ttab, SetNames = names(SetLst))
  return(ret)
  
}

a<- lapply(DEGeneListSynteny, function(x) which(is.na(x)))
DEGeneListSynteny_no_na<- lapply(1:3, function(x) DEGeneListSynteny[[x]][a[[x]]])
varDESeq2<- Plot3WayVenn(DEGeneListSynteny_no_na)
DEMHonly<- varDESeq2[[4]]
DEZSonly<- varDESeq2[[2]]
DEMSUonly<- varDESeq1[[1]]

# get adjusted pvals for all sig DEGs for density plot
DEMHpvals<- which(res.05LFC1MH[,6]<=0.05)
DEMHpvals<- res.05LFC1MH[DEMHpvals,6]

# compare adjusted pvals of unique DEGs between MH and MSU
DEMSUuniquevsMH<- setdiff(DEGeneListSyntenyNEW_no_na[[3]], DEGeneListSyntenyNEW_no_na[[1]])
a<- match(DEMHuniquevsMSU, finalsyntenyNEW$V3)
a<- finalsyntenyNEW$V2[a]
DEMHpvals_unique_vMSU<- res.05LFC1MH[a,]
DEMHpvals_unique_vMSU<- DEMHpvals_unique_vMSU[,6]

temp_res.05MSU<- as.data.frame(res.05LFC1MSU)
a<- temp_res.05MSU[DEMHuniquevsMSU,]
DEMHpvals_unique_vMSU_MSU_pvals<- a[,6]

DEMHunique_pvals_list<- list(MH63=DEMHpvals_unique_vMSU, MSU=DEMHpvals_unique_vMSU_MSU_pvals)

a<- temp_res.05MSU[DEMSUuniquevsMH,]
DEMSUuniquepvals_vMH<- a[,6]
a<- match(DEMSUuniquevsMH, finalsyntenyNEW$V3)
a<- finalsyntenyNEW$V2[a] 
temp_res.05MH<- as.data.frame(res.05LFC1MH)
DEMSUuniquepvals_vMH_MH_pvals<- temp_res.05MH[a,]
DEMSUuniquepvals_vMH_MH_pvals<- DEMSUuniquepvals_vMH_MH_pvals[,6]

DEMSUunique_pvals_list<- list(MH63=DEMSUuniquepvals_vMH_MH_pvals, MSU=DEMSUuniquepvals_vMH)  

# make histograms of adjusted pvals of unique DEGs between MH and MSU
a<- DEMHpvals_unique_vMSU
b<- DEMHpvals_unique_vMSU_MSU_pvals
ah<- hist(a, breaks = seq(0,1,0.025))
bh<- hist(b, alpha.f = 0.6), breaks= seq(0,1,0.025)) 
plot(ah, col=adjustcolor(col="midnightblue", alpha.f = 0.6))
plot(bh, add=T, col=adjustcolor(col="lightseagreen", alpha.f = 0.6))

d<- DEMSUuniquepvals_vMH
e<- DEMSUuniquepvals_vMH_MH_pvals
dh<- hist(d, breaks = seq(0,1,0.025))
eh<- hist(e, breaks= seq(0,1,0.025)) 
plot(dh, col=adjustcolor(col="lightseagreen", alpha.f = 0.6))
plot(eh, add=T, col=adjustcolor(col="midnightblue", alpha.f = 0.6))

# compare DEGs with > 4-fold change
temp_res.05MSUsig<- as.data.frame(res.05SigLFC1MSU)
MSUpos4<- temp_res.05MSUsig[which(temp_res.05MSUsig[,2]>4),]
MSUneg4<- temp_res.05MSUsig[which(temp_res.05MSUsig[,2]< -4),]
MSUnamesFC4<- c(row.names(MSUpos4), row.names(MSUneg4))
a<- match(MSUnamesFC4, finalsyntenyNEW$V3)
MSUnamesFC4<- finalsyntenyNEW$V3[a]
a<- which(!is.na(MSUnamesFC4))
MSUnamesFC4<- MSUnamesFC4[a]

temp_res.05MHsig<- as.data.frame(res.05SigLFC1MH)
MHpos4<- temp_res.05MHsig[which(temp_res.05MHsig[,2]>4),]
MHneg4<- temp_res.05MHsig[which(temp_res.05MHsig[,2]< -4),]
MHnamesFC4<- c(row.names(MHpos4), row.names(MHneg4))
a<- match(MHnamesFC4, finalsyntenyNEW$V2)
MHnamesFC4_MSUnames<- finalsyntenyNEW$V3[a]
a<- which(!is.na(MHnamesFC4_MSUnames))
MHnamesFC4_MSUnames<- MHnamesFC4_MSUnames[a]

temp_res.05ZSsig<- as.data.frame(res.05SigLFC1ZS)
ZSpos4<- temp_res.05ZSsig[which(temp_res.05ZSsig[,2]>4),]
ZSneg4<- temp_res.05ZSsig[which(temp_res.05ZSsig[,2]< -4),]
ZSnamesFC4<- c(row.names(ZSpos4), row.names(ZSneg4)) 
a<- match(ZSnamesFC4, finalsyntenyNEW$V7)
ZSnamesFC4_MSUnames<- finalsyntenyNEW$V3[a]
a<- which(!is.na(ZSnamesFC4_MSUnames))
ZSnamesFC4_MSUnames<- ZSnamesFC4_MSUnames[a]

venn.diagram(list(MH63=MHnamesFC4_MSUnames, ZS97=ZSnamesFC4_MSUnames, MSU=MSUnamesFC4), filename = "./vennFC4All.tif", imagetype = "tiff", fill = c("midnightblue","yellow","lightseagreen"), cex = 1.5, fontface = "bold", cat.cex = 1.5, cat.fontface = "bold", euler.d=TRUE, scaled=TRUE )


####################################################
#                 JunctionSeq                      #
####################################################

# manually installed JunctionSeq package 
# CRAN package dependencies:
install.packages("statmod");
install.packages("plotrix");
install.packages("stringr");
install.packages("Rcpp");
install.packages("RcppArmadillo");
install.packages("locfit");
install.packages("Hmisc");
# Bioconductor dependencies:
source("http://bioconductor.org/biocLite.R");
biocLite();
biocLite("Biobase");
biocLite("BiocGenerics");
biocLite("BiocParallel");
biocLite("GenomicRanges");
biocLite("IRanges");
biocLite("S4Vectors");
biocLite("genefilter");
biocLite("geneplotter");
biocLite("SummarizedExperiment");
biocLite("DESeq2");
install.packages("http://hartleys.github.io/JunctionSeq/install/JunctionSeq_LATEST.tar.gz",
                 repos = NULL,
                 type = "source");

# 4.1 Generating raw counts via QoRT for all bam files being used (in terminal)
# for MH63
java -jar ~/JunctionSeq/QoRTs.jar QC --stranded  --singleEnded --minMAPQ 50 ~/Indica_BAM/S46_Dawn_C.1.bam ~/Indica_genome_files/FinalGFF_MH.gtf ~/JunctionSeq/MH63/MH63_Dawn_C.1
# for ZS97
java -jar ~/JunctionSeq/QoRTs.jar QC --stranded  --singleEnded --minMAPQ 50 ~/ZS_BAM/ZS_Dawn_C.1.bam ~/Indica_genome_files/FinalGTF_ZS.gtf ~/JunctionSeq/ZS97/ZS_Dawn_C.1
# for Nipponbare MSU
java -jar ~/JunctionSeq/QoRTs.jar QC --stranded  --singleEnded --minMAPQ 50 ~/MSU_BAM/1-Dawn_C.1.bam ~/Japonica_genome_files/MSU/FinalGTFMSU_no_na.gtf ~/JunctionSeq/MSU/MSU_Dawn_C.1


# Make sample table and decoder file
sampleTable = data.frame(
row.names = c( "MSU_Dawn_C.1", "MSU_Dawn_C.2", "MSU_Dawn_C.3",
"MSU_Dawn_C.4", "MSU_Dusk_C.1", "MSU_Dusk_C.2", "MSU_Dusk_C.3", "MSU_Dusk_C.4" ),
condition = c(rep("Dawn",4), rep("Dusk", 4)))

decoder<-sampleTable
decoder<-cbind(decoder,c("MSU_Dawn_C.1","MSU_Dawn_C.2","MSU_Dawn_C.3","MSU_Dawn_C.4", "MSU_Dusk_C.1", "MSU_Dusk_C.2","MSU_Dusk_C.3", "MSU_Dusk_C.4"))
colnames(decoder)<-c("condition", "sample.ID")

countFiles <- paste0("/storage/data_1/erins/JunctionSeq/MSU/",
                     decoder$sample.ID,
                     "/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt.gz")

# 4.4 Novel splice junctions
write.table(decoder, file="decoder.txt", sep = "\t", row.names = F, quote=F)

# In terminal
java -jar /storage/data_1/erins/JunctionSeq/QoRTs.jar mergeNovelSplices --minCount 6 --stranded ~/JunctionSeq/MSU ~/JunctionSeq/MSU/decoder.txt ~/Japonica_gene_files/MSU/FinalGTFMSU_no_na.gtf ~/JunctionSeq/MSU

# 5.1 for novel junctions
jscs <- runJunctionSeqAnalyses(sample.files = countFiles,
                               sample.names = decoder$sample.ID,
                               condition=factor(decoder$condition),
                               flat.gff.file = "withNovel.forJunctionSeq.gff.gz",
                               nCores = 1,
                               analysis.type = "junctionsAndExons")

# 5.1 Exporting size factors (optional)
writeSizeFactors(jscs, file = "sizeFactors.txt")

# 5.3 Extracting test results
writeCompleteResults(jscs,
                     outfile.prefix="./test",
                     save.jscs = TRUE)

# 6 Visualization and interpretation
buildAllPlots(jscs=jscs,
              outfile.prefix = "./plots/",
              use.plotting.device = "png",
              FDR.threshold = 0.05)

# plot only one gene at a time
buildAllPlotsForGene(geneID = "LOC_Os12g39630", jscs = jscsMSU, colorRed.FDR.threshold = 0.05)

# match JunctionSeq sig 0.05 genes from each genome to MSU orthologs (see MCScanX script below)
AllDataMH<- jscsMH@featureData@data
AllDataMH0.05<- AllDataMH[which(AllDataMH$padjust< 0.05),]
JSeqMHnames<- AllDataMH0.05$geneID
JSeqMHorthologs<- match(x = JSeqMHnames, table = finalsyntenyNEW$V2)
JSeqMHorthologs<- finalsyntenyNEW[JSeqMHorthologs,]
JSeqMHorthologs<- JSeqMHorthologs$V3

########################################
#               MCScanX                #
########################################

# in terminal
unzip MCscanX.zip
cd MCScanX
make

# generate xyz.blast file for each genome (in terminal)
# blastall -i query_file -d database –p blastp –e 1e-10 –b 5 –v 5 –m 8 –o xyz.blast

blastall -i ~/Reciprocal_best_blast/blastdb/MH63.pacbio.pep.LocusLevel.fa -d ~/Reciprocal_best_blast/blastdb/ZS97_pep_Locus.fa -p blastp -e 1e-10 -b 5 -v 5 -m 8 -o MH63.blast
blastall -i ~/Reciprocal_best_blast/blastdb/MH63.pacbio.pep.LocusLevel.fa -d ~/Reciprocal_best_blast/blastdb/MSU_pep_Locus.fa -p blastp -e 1e-10 -b 5 -v 5 -m 8 -o MH63_MSU.blast
blastall -i ~/Reciprocal_best_blast/blastdb/ZS97_pep_Locus.fa -d ~/Reciprocal_best_blast/blastdb/MSU_pep_Locus.fa -p blastp -e 1e-10 -b 5 -v 5 -m 8 -o ZS97_MSU.blast

# re-format gtf files to hae columns: sp#, gene, start_position, stop_position
MH63gtf<- ParseGFF(filenm = "~/Indica_genome_files/FinalGTF_MH.gtf")
MH63genes<- GetChromLocs(tGFF = MH63gtf, GeneFeat = "gene")
n<- MH63genes$V1
n<- gsub(pattern = "Chr", replacement = "MH", x = n, fixed = T)
MH63genes$V1<- n

m<- MH63genes$V9
m<- gsub(pattern = "g", replacement = "t", x = m, fixed = T)
MH63genes$V9<- m

MH63newgtf<- MH63genes[,c(1, 9, 4, 5)]
write.table(MH63newgtf, file="MH63newgff.gff", sep = "\t", row.names = F, col.names = F, quote = F)

# extracts chromosome locations of a list of genes from a GFF object (ParseGFF())
GetChromLocs<-function(tGFF, qGenes, GeneFeat="gene")
{
  AllCords<-NULL
  for(i in (1:length(tGFF)))
  {
    AllCords<-rbind(AllCords,tGFF[[i]][[GeneFeat]])
  }
  rmna<-which(is.na(AllCords[,9]))
  if(length(rmna)>0)
    AllCords<-AllCords[-rmna,]
  AllCords<-AllCords[qGenes,]
  
  #Sort The Cords by Chr/Start
  AllCords<-AllCords[order(AllCords[,1],AllCords[,4]),]
  return(AllCords)
}

# concatonate files (in terminal)
cat MH63_MSU.blast ZS97_MSU.blast MH63.blast > all.blast
cat ~/Indica_genome_files/MH63newgff.gff ~/Indica_genome_files/ZS97newgff.gff ~/Japonica_genome_files/MSU/MSUnewgff.gff > all.gtf

# put all.blast and all.gtf file in same folder (in terminal)
./MCScanX/MCScanX ./MCScanX_input/all

# match orthologs across all three genomes
# read in collinearity file
MCS<-read.table(file = "~/MCScanX_input/all.collinearity", header = F, comment.char = "#", sep = "\t")

# compile data.frame of genes matched in all three genomes
DriveGroupSynt<-function(MCS)
{
  MHhits<-grep(pattern = "MH", x = MCS$V2)
  MHhits<- MCS[MHhits,]
  MSUhits<- grep(pattern = "LOC", x = MCS$V2)
  MSUhits<- MCS[MSUhits,]
  testMHZS<- grep(pattern = "ZS", x = MHhits$V3)
  testMHZS<- MHhits[testMHZS,]
  testMHMSU<- grep(pattern = "LOC", x = MHhits$V3)
  testMHMSU<- MHhits[testMHMSU,]
  testMSUZS<- grep(pattern = "ZS", x = MSUhits$V3)
  testMSUZS<- MSUhits[testMSUZS,]
  
  SyntLst<-lapply(1:nrow(testMHMSU),function(x) GroupSynt(testMHMSU[x,],qDF1 = testMHZS,qDF2 = testMSUZS))
  Nulls<-which(sapply(SyntLst,function(x) is.null(x)))
  SyntLst<-SyntLst[-Nulls]
  SyntMat<-do.call(rbind,SyntLst)
  
  Sums<-t(apply(SyntMat,1,function(x) as.numeric(c(x[4],x[8],x[12]))))
  Sums<-apply(Sums,1,sum)
  
  SyntMat<-cbind(SyntMat,Sums)
  
  MHGroups<-tapply(1:nrow(SyntMat),SyntMat$V2,function(x) SyntMat[x,])
  SyntMat<-do.call(rbind,lapply(MHGroups,function(x) x[which(x$Sums==min(x$Sums))[1],]))
  
  MSUGroups<-tapply(1:nrow(SyntMat),SyntMat$V3,function(x) SyntMat[x,])
  SyntMat<-do.call(rbind,lapply(MSUGroups,function(x) x[which(x$Sums==min(x$Sums))[1],]))
  
  ZSGroups<-tapply(1:nrow(SyntMat),SyntMat$V7,function(x) SyntMat[x,])
  SyntMat<-do.call(rbind,lapply(ZSGroups,function(x) x[which(x$Sums==min(x$Sums))[1],]))
  
  return(SyntMat)
}

# InDFrow comes from MH->MSU
# qDF1 is MH->ZS
# qDF2 is MUS->ZS
GroupSynt<-function(InDFrow,qDF1,qDF2)
{
  inds1<-grep(InDFrow$V2[1],qDF1$V2)
  inds2<-grep(InDFrow$V3[1],qDF2$V2)
  if(length(inds1)==0 | length(inds2)==0)
    return(NULL)
  
  indsGrid<-expand.grid(inds1,inds2)
  geneGrid<-apply(indsGrid,1,function(x) cbind(InDFrow,qDF1[x[1],],qDF2[x[2],]))
  geneGrid<-do.call(rbind,geneGrid)
  colnames(geneGrid)<-paste("V",1:ncol(geneGrid),sep="")
  
  keeps<-apply(geneGrid,1,function(x) x["V7"]==x["V11"])
  if(!any(keeps))
    return(NULL)
  
  geneGrid<-geneGrid[which(keeps)[1],]
  
  return(geneGrid)
}

finalsyntenyNEW<- DriveGroupSynt(MCS)

############################################
#     Identificiation of exonic SNPs       #
############################################

# SNP caller GATK function
# !/bin/sh
# ref.fa is reference genome used for aligment for bam files ***Chr names must be identical in alignment files**
echo "Please enter SAM/BAM file"
read fruit

java -jar picard.jar ReorderSam I=$fruit O=${fruit}.bam R=ref.fa CREATE_INDEX=TRUE TMP_DIR=~/Path/to/temp 

java -jar picard.jar AddOrReplaceReadGroups I=${fruit}.bam O=${fruit}2.bam SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample TMP_DIR=~/Path/to/temp  

java -jar picard.jar MarkDuplicates I=${fruit}2.bam O=${fruit}3.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=${fruit}output.metrics TMP_DIR=~/Path/to/temp

java -jar GenomeAnalysisTK.jar -T SplitNCigarReads -R ref.fa -I ${fruit}3.bam -o ${fruit}4.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS  

java -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R ref.fa -I ${fruit}4.bam -dontUseSoftClippedBases -stand_call_conf 10.0 -o ${fruit}.vcf 

java -jar GenomeAnalysisTK.jar -T VariantFiltration -R ref.fa -V ${fruit}.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o ${fruit}2.vcf 

#Packages

library(dplyr)
# This function identifies the number of exonic SNPs from the longest annotated transcript for each gene ***GTF format***
SNPperTx2<-function(TorF_SNP_trnascript_vector, gff_of_trnascripts){
  
  t1=NULL
  for(i in 1:length(gff_of_trnascripts[,5])){
    txsum=sum(TorF_SNP_trnascript_vector[gff_of_trnascripts[i,5]:gff_of_trnascripts[i,4]])
    t1=c(t1,txsum)
  }
  return(t1)
}

MH63=read.vcfR(file= "MHmerge.bam2.vcf")
ZS97=read.vcfR(file="ZSmerge.bam2.vcf")
MSU=read.vcfR(file= "MSUmerge.bam2.vcf")

MH63mat=MH63@fix
ZS97mat=ZS97@fix
MSUmat=MSU@fix

for(i in unique(MH63mat[,1])){
  print(length(which(MH63mat[,1]==i)))
}

for(i in unique(ZS97mat[,1])){
  print(length(which(ZS97mat[,1]==i)))
}

for(i in unique(MSUmat[,1])){
  print(length(which(MSUmat[,1]==i)))
}


MH_GTF=read.table(file = "FinalGTF_MH.gtf", sep = "\t", stringsAsFactors = F)
ZS_GTF=read.table(file = "FinalGTF_ZS.gtf", sep = "\t", stringsAsFactors = F)
MSU_GTF=read.table(file = "FinalGTFMSU_no_na.gtf", sep = "\t", stringsAsFactors = F)

# length 12 becuase of 12 rice chr
MH63snp=vector("list", 12)
ZS97snp=vector("list", 12)
MSUsnp=vector("list", 12)

# get chr names
names(MH63snp)<-unique(MH63mat[,1])[1:12]
names(ZS97snp)<-unique(ZS97mat[,1])[1:12]
names(MSUsnp)<-unique(MSUmat[,1])[1:12]

# seperate vcf to chr list
MH63snp=lapply(1:12,function(x) as.numeric(MH63mat[which(MH63mat[,1]==names(MH63snp)[x]),2]))
ZS97snp=lapply(1:12,function(x) as.numeric(ZS97mat[which(ZS97mat[,1]==names(ZS97snp)[x]),2]))
MSUsnp=lapply(1:12,function(x) as.numeric(MSUmat[which(MSUmat[,1]==names(MSUsnp)[x]),2]))

# name snp list
names(MH63snp)<-unique(MH63mat[,1])[1:12]
names(ZS97snp)<-unique(ZS97mat[,1])[1:12]
names(MSUsnp)<-unique(MSUmat[,1])[1:12]

# create for snp position
MH_tx=vector("list", 12)
ZS_tx=vector("list", 12)
MSU_tx=vector("list", 12)

# name snp position list
MH_tx=lapply(1:12,function(x) (MH_GTF[which(MH_GTF[,1]==names(MH63snp)[x]),]))
ZS_tx=lapply(1:12,function(x) (ZS_GTF[which(ZS_GTF[,1]==names(ZS97snp)[x]),]))
MSU_tx=lapply(1:12,function(x) (MSU_GTF[which(MSU_GTF[,1]==names(MSUsnp)[x]),]))

names(MH_tx)<-unique(MH63mat[,1])[1:12]
names(ZS_tx)<-unique(ZS97mat[,1])[1:12]
names(MSU_tx)<-unique(MSUmat[,1])[1:12]

# get gene name
MH_tx=lapply(MH_tx, function(x) x[which(x[,3]=="gene"),])
ZS_tx=lapply(ZS_tx, function(x) x[which(x[,3]=="gene"),])
MSU_tx=lapply(MSU_tx, function(x) x[which(x[,3]=="gene"),])

# create vector from 1 to length of last nt postion that is exonic
TorF_MH=lapply(1:12,function(x) 1:max(MH_tx[[x]][,5]) %in% MH63snp[[x]])
TorF_ZS=lapply(1:12,function(x) 1:max(ZS_tx[[x]][,5]) %in% ZS97snp[[x]])
TorF_MSU=lapply(1:12,function(x) 1:max(MSU_tx[[x]][,5]) %in% MSUsnp[[x]])

# call nt position within vector as a SNP
MH_SNPperTx=lapply(1:12,function(x) SNPperTx2(TorF_SNP_trnascript_vector = TorF_MH[[x]], gff_of_trnascripts = MH_tx[[x]]))
ZS_SNPperTx=lapply(1:12,function(x) SNPperTx2(TorF_SNP_trnascript_vector = TorF_ZS[[x]], gff_of_trnascripts = ZS_tx[[x]]))
MSU_SNPperTx=lapply(1:12,function(x) SNPperTx2(TorF_SNP_trnascript_vector = TorF_MSU[[x]], gff_of_trnascripts = MSU_tx[[x]]))

# get exon names of trnascipt with the most exons
MH_SNPperTx=lapply(1:12,function(y) setNames(MH_SNPperTx[[y]],sapply(strsplit(MH_tx[[y]][,9], split = " "), function(x) x[2])))
ZS_SNPperTx=lapply(1:12,function(y) setNames(ZS_SNPperTx[[y]],sapply(strsplit(ZS_tx[[y]][,9], split = " "), function(x) x[2])))
MSU_SNPperTx=lapply(1:12,function(y) setNames(MSU_SNPperTx[[y]],substring(sapply(strsplit(MSU_tx[[y]][,9], split = " "), function(x) x[2]),first = 1,last = 14)))

MH_SNPperTx=lapply(1:12,function(y) setNames(MH_SNPperTx[[y]],paste(substring(sapply(strsplit(MH_tx[[y]][,9], split = " "), function(x) x[4]),first = 1,last = 15),sapply(strsplit(MH_tx[[y]][,9], split = " "), function(x) x[8]),sep = ":E0")))
ZS_SNPperTx=lapply(1:12,function(y) setNames(ZS_SNPperTx[[y]],paste(substring(sapply(strsplit(ZS_tx[[y]][,9], split = " "), function(x) x[4]),first = 1,last = 15),sapply(strsplit(ZS_tx[[y]][,9], split = " "), function(x) x[8]),sep = ":E0")))
MSU_SNPperTx=lapply(1:12,function(y) setNames(MSU_SNPperTx[[y]],paste(substring(sapply(strsplit(MSU_tx[[y]][,9], split = " "), function(x) x[2]),first = 1,last = 16),sapply(strsplit(MSU_tx[[y]][,9], split = " "), function(x) x[6]),sep = ":E00")))

MH_SNPperTx=lapply(1:12,function(x) setNames(MH_SNPperTx[[x]],paste(substring(names(MH_SNPperTx[[x]]),first = 1, last=12), substring(names(MH_SNPperTx[[x]]),first = 16, last=20), sep = "")))
ZS_SNPperTx=lapply(1:12,function(x) setNames(ZS_SNPperTx[[x]],paste(substring(names(ZS_SNPperTx[[x]]),first = 1, last=12), substring(names(ZS_SNPperTx[[x]]),first = 16, last=20), sep = "")))
MSU_SNPperTx=lapply(1:12,function(x) setNames(MSU_SNPperTx[[x]],paste(substring(names(MSU_SNPperTx[[x]]),first = 1, last=14), substring(names(MSU_SNPperTx[[x]]),first = 16, last=20), sep = "")))

MH_SNPperTx=lapply(1:12,function(x) MH_SNPperTx[[x]][unique(names(MH_SNPperTx[[x]]))])
ZS_SNPperTx=lapply(1:12,function(x) ZS_SNPperTx[[x]][unique(names(ZS_SNPperTx[[x]]))])
MSU_SNPperTx=lapply(1:12,function(x) MSU_SNPperTx[[x]][unique(names(MSU_SNPperTx[[x]]))])

# calculate SNPs per exon
MH_SNPexon<- unlist(MH_SNPperTx)
MH_SNPexon<- as.matrix(MH_SNPexon)
MH_SNPexonNames<- gsub(pattern = "-[[:digit:]]{2}", replacement = "", x = row.names(MH_SNPexon))
row.names(MH_SNPexon)<- MH_SNPexonNames
colnames(MH_SNPexon)<- c("Exon_SNPs")

ZS_SNPexon<- unlist(ZS_SNPperTx)
ZS_SNPexon<- as.matrix(ZS_SNPexon)
colnames(ZS_SNPexon)<- c("Exon_SNPs")

MSU_SNPexon<- unlist(MSU_SNPperTx)
MSU_SNPexon<- as.matrix(MSU_SNPexon)
colnames(MSU_SNPexon)<- c("Exon_SNPs")

# re-format to add only gene_id in additional column that can be matched with finalsyntenyNEW 
# for MH SNPs
MH_SNPexonDF<- matrix(ncol = 2, nrow = 246338)
MH_SNPexonDF[,1]<- MH_SNPexon[,1]
row.names(MH_SNPexonDF)<- row.names(MH_SNPexon)
MH_SNPexonDF<- as.data.frame(MH_SNPexonDF)
MH_SNPexonDF$V2<- row.names(MH_SNPexon)
MH_SNPexonDF$V2<- gsub(pattern = ":[[:alpha:]]{1}[[:digit:]]{3}", replacement = "", x = MH_SNPexonDF$V2)

# for MSU SNPs
MSU_SNPexonDF<- matrix(ncol = 2, nrow = 214173)
MSU_SNPexonDF[,1]<- MSU_SNPexon[,1]
row.names(MSU_SNPexonDF)<- row.names(MSU_SNPexon)
MSU_SNPexonDF<- as.data.frame(MSU_SNPexonDF)
MSU_SNPexonDF$V2<- row.names(MSU_SNPexon)
MSU_SNPexonDF$V2<- gsub(pattern = ":[[:alpha:]]{1}[[:digit:]]{3}", replacement = "", x = MSU_SNPexonDF$V2)

# for ZS SNPs
ZS_SNPexonDF<- matrix(ncol = 2, nrow = 232892)
ZS_SNPexonDF[,1]<- ZS_SNPexon[,1]
row.names(ZS_SNPexonDF)<- row.names(ZS_SNPexon)
ZS_SNPexonDF<- as.data.frame(ZS_SNPexonDF)
ZS_SNPexonDF$V2<- row.names(ZS_SNPexon)
ZS_SNPexonDF$V2<- gsub(pattern = ":[[:alpha:]]{1}[[:digit:]]{3}", replacement = "", x = ZS_SNPexonDF$V2)

# match syntenic orthologs
MH_SNPexonDFSynNames<- intersect(x = finalsyntenyNEW$V2, y = MH_SNPexonDF$V2)
MHSNPind<- which(MH_SNPexonDF$V2 %in% MH_SNPexonDFSynNames)
MH_SNPexonDFSyn<- MH_SNPexonDF[MHSNPind,]

ZS_SNPexonDFSynNames<- intersect(x = finalsyntenyNEW$V7, y = ZS_SNPexonDF$V2)
ZSSNPind<- which(ZS_SNPexonDF$V2 %in% ZS_SNPexonDFSynNames)
ZS_SNPexonDFSyn<- ZS_SNPexonDF[ZSSNPind,]

MSU_SNPexonDFSynNames<- intersect(x = finalsyntenyNEW$V3, y = MSU_SNPexonDF$V2)
MSUSNPind<- which(MSU_SNPexonDF$V2 %in% MSU_SNPexonDFSynNames)
MSU_SNPexonDFSyn<- MSU_SNPexonDF[MSUSNPind,]

# This function sums the number of exonic SNPs per gene
# use objects MH_SNPexonDFSyn, ZS_SNPexonDFSyn, and MSU_SNPexonDFSyn to compare exonic SNPs per gene for syntenic genes
# inputs:
# data.frame of gene_exon_ids, SNPs/exon, gene_ids that match format of gtf and synteny files (e.g. MH_SNPexonDFSyn)
# list of syntenic gene_ids (e.g. MH_SNPexonDFSynNames)

ExonSNPCount<- function(df, genenames)
{
  exonsum<- matrix(ncol = 1, nrow = length(genenames))
  row.names(exonsum)<- genenames
  colnames(exonsum)<- c("exonic_SNPs_per_gene")
  for(i in genenames)
  {
    geneSNPind<- which(df$V2 %in% i)
    geneSNP<- df[geneSNPind,]
    geneSNPcount<- sum(geneSNP$V1)

    exonsum[i,1]<- geneSNPcount
  }
  return(exonsum)
}

# calculate exonic SNPs/gene
MH_SNPexonSynSums<- ExonSNPCount(df = MH_SNPexonDFSyn, genenames = MH_SNPexonDFSynNames)
ZS_SNPexonSynSums<- ExonSNPCount(df = ZS_SNPexonDFSyn, genenames = ZS_SNPexonDFSynNames)
MSU_SNPexonSynSums<- ExonSNPCount(df = MSU_SNPexonDFSyn, genenames = MSU_SNPexonDFSynNames)

# assemble matrix of exonic SNPs
exonicSNPs<- matrix(ncol = 3, nrow = 21145)
exonicSNPs[,1]<- MH_SNPexonSynSums[,1]
exonicSNPs[,2]<- ZS_SNPexonSynSums[,1]
exonicSNPs[,3]<- MSU_SNPexonSynSums[,1]
row.names(exonicSNPs)<- row.names(MH_SNPexonSynSums)
colnames(exonicSNPs)<- c("MH63", "ZS97", "MSU")

# note: see below for values of j and naAll
exonicSNPs<- exonicSNPs[-j,]
exonicSNPs<- exonicSNPstest[-naAll,]

##################################################
#         Calculate counts per gene              #
##################################################

# counts per detectable gene in dds object from DESeq2
countsMHrows<- counts(ddsMH)
countsMHrows<- rowSums(countsMHrows)
countsMHsyn<- match(finalsyntenyNEW$V2, names(countsMHrows))
countsMHsyn<- countsMHrows[countsMHsyn]

countsZSAll<- counts(ddsZS)
countsZSAllrows<- rowSums(countsZSAll)
countsZSsyn<- match(finalsyntenyNEW$V7, names(countsZSAllrows))
countsZSsyn<- countsZSAllrows[countsZSsyn]

countsMSUrows<- counts(ddsMSU)
countsMSUrows<- rowSums(countsMSUrows)
countsMSUsyn<- match(finalsyntenyNEW$V3, names(countsMSUrows))
countsMSUsyn<- countsMSUrows[countsMSUsyn]

countsmat<- matrix(data = countsMHsyn, ncol = 3, nrow = 21145)
countsmat[,2]<- countsZSsyn
countsmat[,3]<- countsMSUsyn
row.names(countsmat)<- names(countsMSUsyn)
j<- which(is.na(countsmat))
countsmat_no_na<- countsmat[-j]
nacount<- which(apply(countsmat_no_na, 1, function(x) any(is.na(x)))
naAll<- union(nacount, exonicSNPs)
nacount<- countsmat_no_na[-naAll,]

###################################################
#   Calculate gene length and exonic SNP density  #
###################################################

# exonic SNP density= exonic SNPs per gene/gene length
# determine gene length
# for MH63
MH63genelength<- matrix(nrow = 57174, ncol = 2)
MH63genelength[,1]<- MH63newgtf$V4
MH63genelength[,2]<- MH63newgtf$V5
row.names(MH63genelength)<- MH63newgtf$V9
b<- sapply(1:nrow(MH63genelength), function(x) abs(MH63genelength[x,2]-MH63genelength[x,1]))
MH63genelength<- cbind(MH63genelength, b)
colnames(MH63genelength)<- c("gene_start", "gene_stop", "gene_length")

# for ZS97
ZS97genelength<- matrix(nrow = 54831, ncol = 2)
ZS97genelength[,1]<- ZS97newgtf$V4
ZS97genelength[,2]<- ZS97newgtf$V5
row.names(ZS97genelength)<- ZS97newgtf$V9
d<- sapply(1:nrow(ZS97genelength), function(x) abs(ZS97genelength[x,2]-ZS97genelength[x,1]))
ZS97genelength<- cbind(ZS97genelength, d)
colnames(ZS97genelength)<- c("gene_start", "gene_stop", "gene_length")

# for MSU
MSUgenelength<- matrix(nrow = 55801, ncol = 2)
MSUgenelength[,1]<- MSUnewgtf$V4
MSUgenelength[,2]<- MSUnewgtf$V5
row.names(MSUgenelength)<- MSUnewgtf$V9
e<- sapply(1:nrow(MSUgenelength), function(x) abs(MSUgenelength[x,2]-MSUgenelength[x,1]))
MSUgenelength<- cbind(MSUgenelength, e)
colnames(MSUgenelength)<- c("gene_start", "gene_stop", "gene_length")

# gene length for only syntenic orthologs
MH63genelengthsyn<- match(x = finalsyntenyNEW$V2, row.names(MH63genelength))
MH63genelengthsyn<- MH63genelength[MH63genelengthsyn,]

ZS97genelengthsyn<- match(x = finalsyntenyNEW$V7, table = row.names(ZS97genelength))
ZS97genelengthsyn<- ZS97genelength[ZS97genelengthsyn,]

MSUgenelengthsyn<- match(x = finalsyntenyNEW$V3, table = row.names(MSUgenelength))
MSUgenelengthsyn<- MSUgenelength[MSUgenelengthsyn,]

genelength<- matrix(ncol = 3, nrow = 21145)
genelength[,1]<- MH63genelengthsyn[,3]
genelength[,2]<- ZS97genelengthsyn[,3]
genelength[,3]<- MSUgenelengthsyn[,3]
row.names(genelength)<- row.names(MSUgenelengthsyn)
colnames(genelength)<- c("MH63", "ZS97", "MSU")
genelength_no_na<- genelength[-j,]
nagenelength<- genelength_no_na[-naAll,]
genelengthsubset<- interset(row.names(nacount), row.names(nagenelength))
genelengthsubset<- nagenelength[genelengthsubset,]

# exonic SNP density
exonicSNPdensity<- data.matrix(exonicSNPs/genelengthsubset)

#########################################################
#       Calculate exon length and exon number           #
#########################################################

MH63gtf_tx<- ParseGFF(filenm = "~/Indica_genome_files/FinalGTF_MH.gtf") #change to 5 for transcript_id in AccNms<-sapply(strsplit(x = nms, split = "\\ |\\;"),function(x) x[5])
MH63exons<- GetChromLocs(tGFF = MH63gtf_tx, GeneFeat = "exon")
ZS97gtf_tx<- ParseGFF(filenm = "~/Indica_genome_files/FinalGTF_ZS.gtf")
ZS97exons<- GetChromLocs(tGFF = ZS97gtf_tx, GeneFeat = "exon")
MSUgtf_tx<- ParseGFF(filenm = "~/Japonica_genome_files/MSU/FinalGTFMSU_no_na.gtf")
MSUexons<- GetChromLocs(tGFF = MSUgtf_tx, GeneFeat = "exon")

# for MH63
exonlength<- abs(MH63exons$V5-MH63exons$V4)
exonlist<- tapply(exonlength,INDEX = MH63exons[,9], function(x) x)
geneidvec<- MH63newgtf$V9

exonsums<- lapply(geneidvec, function(x) exonlist[grep(pattern = x, names(exonlist))])
names(exonsums)<- geneidvec
exonsum2<- lapply(exonsums, function(x) lapply(x, function(y) sum(y)))

exonunlist<- unlist(exonsum2)
exonunlist<- sapply(geneidvec, function(x) exonunlist[grep(pattern = x, names(exonunlist))])
exonmax<- lapply(exonunlist, function(x) max(x))
MH63exonlength<- as.matrix(exonmax)
colnames(MH63exonlength)<- c("exon_length")

a<- sapply(exonsums, function(x) sapply(x, function(y) length(y)))
MHexonnumber<- sapply(a, function(x) max(x))

# for ZS97
exonlength<- abs(ZS97exons$V5-ZS97exons$V4)
exonlist<- tapply(exonlength,INDEX = ZS97exons[,9], function(x) x)
geneidvec<- ZS97newgtf$V9

exonsums<- lapply(geneidvec, function(x) exonlist[grep(pattern = x, names(exonlist))])
names(exonsums)<- geneidvec
exonsum2<- lapply(exonsums, function(x) lapply(x, function(y) sum(y)))

exonunlist<- unlist(exonsum2)
exonunlist<- sapply(geneidvec, function(x) exonunlist[grep(pattern = x, names(exonunlist))])
exonmax<- lapply(exonunlist, function(x) max(x))
ZS97exonlength<- as.matrix(exonmax)
colnames(ZS97exonlength)<- c("exon_length")

# for MSU
exonlength<- abs(MSUexons$V5-MSUexons$V4)
exonlist<- tapply(exonlength,INDEX = MSUexons[,9], function(x) x)
geneidvec<- MSUnewgtf$V9

exonsums<- lapply(geneidvec, function(x) exonlist[grep(pattern = x, names(exonlist))])
names(exonsums)<- geneidvec
exonsum2<- lapply(exonsums, function(x) lapply(x, function(y) sum(y)))

exonunlist<- unlist(exonsum2)
exonunlist<- lapply(geneidvec, function(x) exonunlist[grep(pattern = x, names(exonunlist))])
exonmax<- lapply(exonunlist, function(x) max(x))
MSUexonlength2<- as.matrix(exonmax)
colnames(MSUexonlength)<- c("exon_length")

a<- sapply(exonsums, function(x) sapply(x, function(y) length(y)))
MSUexonnumber<- sapply(a, function(x) max(x))

# assemble exonlength matrix
a<- sapply(MH63exonlength, function(x) sapply(x, function(y) y))
a<- unlist(a)
a<- as.matrix(a)
test<-matrix(nrow = 57174, ncol=1)
test[,1]<-a
row.names(test)<- MH63genes$V9
b<- match(x = finalsyntenyNEW$V2, table = row.names(test))
MH63exonsyn<- test[b,]
MH63exonsyn<- as.matrix(MH63exonsyn)

a<- sapply(ZS97exonlength, function(x) sapply(x, function(y) y))
a<- unlist(a)
a<- as.matrix(a)
test<-matrix(nrow = 54831, ncol=1)
test[,1]<-a
row.names(test)<- ZS97genes$V9
b<- match(x = finalsyntenyNEW$V7, table = row.names(test))
ZS97exonsyn<- test[b,]
ZS97exonsyn<- as.matrix(ZS97exonsyn)

a<- sapply(MSUexonlength, function(x) sapply(x, function(y) y))
a<- unlist(a)
a<- as.matrix(a)
test<-matrix(nrow = 55801, ncol=1)
test[,1]<-a
row.names(test)<- MSUgenes$V9
b<- match(x = finalsyntenyNEW$V3, table = row.names(test))
MSUexonsyn<- test[b,]
MSUexonsyn<- as.matrix(MSUexonsyn)

exonlength<- matrix(ncol = 3, nrow = 21145)
exonlength[,1]<- MH63exonsyn[,1]
exonlength[,2]<- ZS97exonsyn[,1]
exonlength[,3]<- MSUexonsyn[,1]
row.names(exonlength)<- row.names(MSUexonsyn) 
colnames(exonlength)<- c("MH63", "ZS97", "MSU")

# syntenic genes
exonlengthsubset<- intersect(row.names(nacount), row.names(exonlength))
exonlengthsubset<- exonlength[exonlengthsubset,]

# assemble exonnumber matrix
MSUexonnumbersyn<- intersect(finalsyntenyNEW$V3, names(MSUexonnumber))
MSUexonnumbersyn<- MSUexonnumber[MSUexonnumbersyn]
MHexonnumbersyn<- intersect(finalsyntenyNEW$V2, names(MHexonnumber))
MHexonnumbersyn<- MHexonnumber[MHexonnumbersyn]

MHMSUexonnumber<- matrix(nrow=21145, ncol=2)
MHMSUexonnumber[,1]<- MHexonnumbersyn
MHMSUexonnumber[,2]<- MSUexonnumbersyn
row.names(MHMSUexonnumber)<- names(MSUexonnumbersyn)
colnames(MHMSUexonnumber)<- c("MH_exon", "MSU_exon")

###################################################
#             Linear regression model             #
###################################################

# assemble log2 values of MH63:MSU ratios of response and explanatory variables to use for linear model
# counts
counts_ratio_log<- nacount[,1]/nacount[,3]
counts_ratio_log<- log2(counts_ratio_log)
a<- names(counts_ratio_log)

# exon length
exlength_ratio_log<- exonlength[,1]/exonlength[,3]
exlength_ratio_log<- log2(exlength_ratio_log)
exlength_ratio_log<- exlength_ratio_log[a]

# gene length
genelength_ratio_log<- genelengthsubset[,1]/genelengthsubset[,3]
genelength_ratio_log<- log2(genelength_ratio_log)
genelength_ratio_log<- genelength_ratio_log[a]

# exon number
MHMSUexonnumberratiolog<- MHMSUexonnumber[,1]/MHMSUexonnumber[,2]
MHMSUexonnumberratiolog<- log2(MHMSUexonnumberratiolog)
exnumber_ratio_log<- MHMSUexonnumberratiolog[a]

# exonic SNP density
snp_ratio<- exonicSNPsdensity
snp_ratio<- snp_ratio+0.0000000001
snp_ratio<- snp_ratio[,1]/snp_ratio[,3]
snp_ratio_log<- log2(snp_ratio)
snp_ratio_log<- snp_ratio_log[a]

# assemble data frame
df<- cbind(counts_ratio_log, exlength_ratio_log, genelength_ratio_log, exnumber_ratio_log, snp_ratio_log)
df<- as.data.frame(df)

lmfinal<- lm(df$counts_ratio_log ~ df$exlength_ratio_log + df$genelength_ratio_log + df$exnumber_ratio_log + df$snp_ratio_log)
summary(lmfinal)

Call:
  lm(formula = df$counts_ratio_log ~ df$exlength_ratio_log + 
       df$genelength_ratio_log + df$exnumber_ratio_log + 
       df$snp_ratio_log)

Residuals:
  Min       1Q   Median       3Q      Max 
-13.7960  -0.1939   0.0435   0.1959  11.8888 

Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
(Intercept)                   -0.1302786  0.0101702 -12.810  < 2e-16 ***
  df$exlength_ratio_log        0.4512914  0.0169183  26.675  < 2e-16 ***
  df$genelength_ratio_log      0.1658273  0.0197290   8.405  < 2e-16 ***
  df$exnumber_ratio_log        0.2233650  0.0205788  10.854  < 2e-16 ***
  df$snp_ratio_log             0.0015790  0.0006088   2.594  0.00951 ** 
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 1.075 on 17034 degrees of freedom
Multiple R-squared:  0.1595,	Adjusted R-squared:  0.1593 
F-statistic: 808.1 on 4 and 17034 DF,  p-value: < 2.2e-16

anovalm<- anova(lmfinal)
Analysis of Variance Table

Response: genematall$counts_ratio_log
                          Df  Sum Sq Mean Sq   F value    Pr(>F)    
df$exlength_ratio_log     1  3192.2  3192.2 2760.3421 < 2.2e-16 ***
df$genelength_ratio_log   1   402.6   402.6  348.1157 < 2.2e-16 ***
df$exnumber_ratio_log     1   135.5   135.5  117.1687 < 2.2e-16 ***
df$snp_ratio_log          1     7.8     7.8    6.7267  0.009506 ** 
Residuals                     17034 19698.9     1.2                        
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

