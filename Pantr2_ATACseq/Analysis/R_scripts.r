library(DiffBind)
library(ChIPQC)
library(ChIPpeakAnno)
library(tidyverse)
library(EnsDb.Mmusculus.v79)
library(Repitools)
library(motifmatchr)
library(MotifDb)
library(JASPAR2020)
library(TFBSTools)
library(chromVAR)
library(BSgenome.Mmusculus.UCSC.mm10)
library(pheatmap)
library(memes)
library(GenomicRanges)
library(clusterProfiler)
BiocParallel::register(BiocParallel::MulticoreParam(12), default = T)

#read in the sample sheet
samples<-read.csv("../Assets/samples.csv")

#create the dba object
dat<-dba(sampleSheet = samples)
dat.qc<-ChIPQC(samples, annotation = "mm10")

pdf("../Plots/coverage_estimate.pdf")
plotCoverageHist(dat.qc, facetBy = "Condition") + ggpubr::theme_pubr()
dev.off()

pdf("../Plots/cross-coverage.pdf", width = 20)
plotCC(dat.qc, facetBy = "Condition") + ggpubr::theme_pubr()
dev.off()

plotPrincomp(dat.qc)

#changed summits to 100 for a 201bp window around the summit that is recommended for ATACseq
dat<-dba.blacklist(dat, blacklist=DBA_BLACKLIST_MM10, greylist=FALSE)
dat<-dba.count(dat, bParallel = T, summits = 100, bUseSummarizeOverlaps = T)
dat<-dba.normalize(dat)

dat<-dba.contrast(dat, reorderMeta = list(Condition = "WT"), minMembers = 2)
dat<-dba.analyze(dat, bParallel = T)

pdf("../Plots/volcano_plot.pdf")
dba.plotVolcano(dat, bUsePval = T)
dev.off()

report<-dba.report(dat, bUsePval = T)
write.csv(report, "../Tables/diff_peaks.csv")
#for better volcano plot
all.plot<-dba.report(dat, th =1, bUsePval = T)
annoData<-toGRanges(EnsDb.Mmusculus.v79)
all.plot<-annoPeaks(all.plot, annoData, bindingType = "startSite", bindingRegion = c(-15000, 3000))
names(all.plot)<-NULL
outs<-as.data.frame(all.plot)
outs$label<-ifelse(outs$p.value < 0.003, outs$gene_name, NA)

pdf("../Plots/better_volcano_plot.pdf")
ggplot(outs, aes(x=Fold, y=-log10(p.value), color = -log10(p.value) > -log10(0.05), label = label)) + 
  geom_point() + 
  ggpubr::theme_pubr(legend = "none") + 
  xlab("log2(Fold Change)") + 
  ylab("-log10(p-value)") +
  scale_color_manual(values = c("black", "red")) +
  ggtitle("WT vs Pantr2 KO [847 p<-0.05]") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggrepel::geom_text_repel()
dev.off()

#annotating peak features for differential peaks
annoData<-toGRanges(EnsDb.Mmusculus.v79)
report<-annoPeaks(report, annoData, bindingType = "startSite", bindingRegion = c(-15000, 3000))
names(report)<-NULL
View(as.data.frame(report))
dim(as.data.frame(report))
write_csv(as.data.frame(report), "../Tables/diffpeaks_annotated.csv")
export.bed(object = report, con = "../Tables/diffpeaks.bed", format = "bed")

#motif enrichment analysis for differential peaks
seq<-report %>% get_sequence(BSgenome.Mmusculus.UCSC.mm10)
options(meme_db = "~/meme/bin/")
enr_mo<-runAme(seq, database = "../Assets/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt")
write.csv(enr_mo, "../Tables/diffpeaks_enriched_motifs.csv")
enr_mo<-read.csv("../Tables/diffpeaks_enriched_motifs.csv")

#motif enrichment using control samples (this results in no enriched peaks vs the control)
all_seq<-list(dat$peaks[[1]] %>% get_sequence(BSgenome.Mmusculus.UCSC.mm10), dat$peaks[[2]] %>% get_sequence(BSgenome.Mmusculus.UCSC.mm10), dat$peaks[[3]] %>% get_sequence(BSgenome.Mmusculus.UCSC.mm10), dat$peaks[[4]] %>% get_sequence(BSgenome.Mmusculus.UCSC.mm10))
names(all_seq)<-c("AW1", "AX3", "AW4", "AW2")
all_enr_mo<-runAme(all_seq, control = c("AW4", "AW2"))
do.call(rbind, all_enr_mo)
plot_ame_heatmap(all_enr_mo[1:20,], id = motif_alt_id)

#Gene Ontology analysis for differential peaks list
ugh<-as.data.frame(report)
genes<-unique(ugh$gene_name)
genes<-bitr(geneID = genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")$ENTREZID
go_bp<-enrichGO(genes,OrgDb = "org.Mm.eg.db", keyType = "ENTREZID", ont = "BP", readable = TRUE, pvalueCutoff  = 0.05)
pdf("../Plots/GO_diffpeaks.pdf")
dotplot(go_bp)
dev.off()
result_gobp_peaks<-as.data.frame(go_bp)

genes<-str_to_title(enr_mo$motif_alt_id)
genes<-bitr(geneID = genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")$ENTREZID
go_bp<-enrichGO(genes,OrgDb = "org.Mm.eg.db", keyType = "ENTREZID", ont = "BP", readable = TRUE, pvalueCutoff  = 0.05)
pdf("../Plots/GO_diffmotifs.pdf")
dotplot(go_bp)
dev.off()
result_gobp_motifs<-as.data.frame(go_bp)

write.csv(result_gobp_motifs, "../Tables/GO_enriched_motifs.csv")
write.csv(result_gobp_peaks, "../Tables/GO_diff_peaks.csv")

tmp<-as.data.frame(go_bp)
tmp<-tmp[c(4,9,11,15,19,46,49,58,61,64,69,80,82,101, 100, 134, 135, 138, 180, 203, 274, 333, 356, 360, 430, 434, 523, 532, 556),]
tmp$GeneRatio<-sapply(tmp$GeneRatio, function(x) eval(parse(text=x)))
tmp<-tmp[order(tmp$GeneRatio),]
tmp$Description<-as.factor(tmp$Description)
reorder(tmp$Description, tmp$GeneRatio)
pdf("../Plots/better_go_motifs.pdf", width = 8)
ggplot(tmp, # you can replace the numbers to the row number of pathway of your interest
       aes(x = GeneRatio, y = reorder(Description, GeneRatio))) + 
  geom_point(aes(size = Count, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.05), low="red") +
  ylab(NULL) +
  ggtitle("GO pathway enrichment")
dev.off()
