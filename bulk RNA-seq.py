cd ~/reference/genome/buffalo/hisat2_index
hisat2_extract_exons.py ../GCF_003121395.1_ASM312139v1_genomic.gtf > exons_buffalo.txt
hisat2_extract_splice_sites.py ../GCF_003121395.1_ASM312139v1_genomic.gtf >ss_buffalo.txt
nohup hisat2-build -p 20 --ss ss_buffalo.txt --exon exons_buffalo.txt ../GCF_003121395.1_ASM312139v1_genomic.fna  buffalo &
cd ~/buffalo_bulk
mkdir {raw,clean,align}
mv *.gz ./raw
cd raw 
mkdir rawqc
cd rawqc
ls ../*gz|xargs fastqc -t 10 -o  ./  
multiqc .
cd ~/buffalo_bulk/raw
ls *_1.fq.gz >1
ls *_2.fq.gz >2
cat 1
cat 2
paste 1 1 2 > config.raw
cat config.raw  |while read id;
do echo $id
arr=($id)
fq2=${arr[2]}
fq1=${arr[1]}
sample=${arr[0]}
nohup  trim_galore -q 25 --phred33 --length 36 -e 0.1 --stringency 4 --paired -o  ../clean  $fq1   $fq2  & 
done 
cd ~/buffalo_bulk/clean
mkdir cleanqc
cd cleanqc
ls ../*gz|xargs fastqc -t 10 -o  ./  
multiqc .
ls *_1_val_1.fq.gz >0
index=~/reference/genome/buffalo/hisat2_index
inputdir=~/buffalo_bulk/clean
outdir=~/buffalo_bulk/align
cat 0 | while read id
do
  echo "hisat2 -p 10 -x ${index} -1 ${inputdir}/${id}_1_val_1.fq.gz -2 ${inputdir}/${id}_2_val_2.fq.gz 2>${id}.log  | samtools sort -@ 5 -o ${outdir}/${id}.Hisat_aln.sorted.bam -  && samtools index ${outdir}/${id}.Hisat_aln.sorted.bam ${outdir}/${id}.Hisat_aln.sorted.bam.bai"
done >Hisat.sh

nohup sh Hisat.sh >Hisat.log &
gtf="~/reference/genome/genome/buffalo/GCF_003121395.1_ASM312139v1_genomic.gtf"   
ls -lh $gtf 
ls *bam |while read id;do
nohup qualimap rnaseq --java-mem-size=20G -gtf $gtf -bam $id  -pe -oc ${id%%.*}  & 
done 
gtf="~/reference/genome/genome/buffalo/GCF_003121395.1_ASM312139v1_genomic.gtf"   
nohup featureCounts -T 5 -p -t exon -g gene_id  -a $gtf \
-o  all.id.txt  *.bam  1>counts.id.log 2>&1 &
nohup featureCounts -T 5 -p -t exon -g transcript_id  -a $gtf \
-o  all.transcript_id.name.txt  *.bam  1>counts.transcript_id.log 2>&1 &
nohup featureCounts -f -O -s 2 -p -T 5 \
-F GTF -a  $gtf  \
-o  all.exon.txt  *.bam  1>counts.exon.log 2>&1 &
multiqc all.id.txt.summary

##Running in R
##转录组下游分析
library(RNAseqStat)
library(edgeR)
library(airway)
rm(list=ls())
Sys.setenv(LANGUAGE = "en")
library("parathyroidSE")
setwd("~/buffalo_bulk/coount")
rawcount <- read.table("all.id.txt", header=TRUE, row.names=1,sep='\t')
group_list <- c(rep("ILC", 3), rep("ALC", 3))
library(clusterProfiler)
library(AnnotationHub)
ah <- AnnotationHub()
ah[ah$species=="Bubalus bubalis"]
ah[ah$species=="Bubalus bubalis" & ah$rdataclass=="OrgDb"]
bbub_orgdb <- ah[["AH72312"]]
saveDb(bbub_orgdb, file="bbub.orgdb")
table(rownames(rawcount)%in%ids$ENTREZID)
expr <- rawcount[rownames(rawcount)%in%ids$ENTREZID,]
ids <- ids[match(rownames(expr),ids$ENTREZID),]
rownames(count) <- ids$SYMBOL
expr<-rawcount
dup_gene<-data.frame(gene=rownames(expr),mean=apply(expr,1,mean))
expr <-expr[order(dup_gene$gene,dup_gene$mean,decreasing = T),]
expr <- expr[!duplicated(rownames(expr)),]
keep <- rowSums(expr > 0) >= floor(0.75 * ncol(expr))
expr <- expr[keep,]
count_data <- as.data.frame(expr)
dir.create("./deg_result_onestep")
dir="./deg_result_onestep/"
?runAll
runAll(count_data = count_data, group_list = group_list, 
       OrgDb = 'bbub_orgdb', dir = dir)
#1. pre_check
dir.create("./deg_result")
dir="./deg_result/"
?pre_check
pre_check(counts_data = count_data, group_list = group_list, dir = dir)
#2. deg_run
?deg_run
deg_results <- deg_run(count_data, group_list, 
                       test_group = "ALC", 
                       control_group = "ILC",
                       dir = dir)
#3. enhance_volcano
cut_FC = 2
cut_P = 0.05
OrgDb <- bbub_orgdb
names(deg_results@deg_df_DESeq2)
x = "log2FoldChange"
y = "pvalue"
?enhance_volcano
res <- enhance_volcano(deg_results@deg_df_DESeq2,x = x, y = y,
                       label = c("Down","Stable","Up"), label_ns = "Stable",
                       palette = c('#66c2a5', 'grey50', '#fc8d62'),
                       cut_FC = cut_FC,cut_P = cut_P,top = 10, size = 2.0,
                       expand = c(0.25,0.25),
                       highlight = NULL)
ggplot2::ggsave(res,filename = file.path(dir,"volcano.pdf"),
       width = 1600,height = 1600,units = "px",limitsize = FALSE)
#4. enrichGO_run
enrichGO_run(deg_results@deg_df_DESeq2, x = "log2FoldChange", y = "pvalue",
             OrgDb = OrgDb, dir = dir,prefix = "3-EnrichGO-DESeq2",
             cut_FC = cut_FC,cut_P = cut_P)
enrichKEGG_run(deg_results@deg_df_DESeq2, x = "log2FoldChange", 
      y = "pvalue", OrgDb = OrgDb, dir = dir, prefix = "4-EnrichKEGG-DESeq2",
      cut_FC = cut_FC,cut_P = cut_P)
enrichgesKEGG_run(deg_results@deg_df_DESeq2, x = "log2FoldChange", 
    OrgDb = OrgDb, dir = dir, prefix = "5-EnrichgseKEGG-DESeq2")  
