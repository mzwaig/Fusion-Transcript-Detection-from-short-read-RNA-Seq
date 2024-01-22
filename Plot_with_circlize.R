#setwd("")
library(dplyr)
library(tidyr)
library(ggplot2)
require(data.table) ## 1.9.4+
library(openxlsx)
library(circlize)
library(RColorBrewer)

datasets = read.delim(file= paste0("datasets.tsv"))
AllSamples_AllCalls = read.delim(file="Fusion_filtering_Output/AllSamples_AllCallers.noFiltering.tsv")
#Pull out genes of interest 
GOI = c("BRAF", "RAF1") 
GOI.fusion.calls= AllSamples_AllCalls %>% filter(gene1 %in% GOI | gene2 %in% GOI)



# Load GTF to add gene names  ####
print("Load gtf")
# gene_span_interval=1000
# gtf <- fread(paste0('refdata-b37-2.1.0_gene_annotations.gtf'), header = FALSE, sep = '\t')
# colnames(gtf) = c("chromosome","source", "feature_type", "start", "end","score", "strand", "phase", "info")
# 
# gtf = gtf %>% filter(feature_type=="gene") %>% mutate(info=gsub("; ", ";", info)) %>% mutate(info=gsub('"', "", info)) %>% 
#   separate(info, c("X1","gene_id", "X2","gene_name","X3", "gene_source", "X4", "gene_biotype"),
#            sep = "([\\ \\;])") %>% select(-contains("X"),-source)
# 
# gtf2 = gtf %>% mutate(start.exp= start - gene_span_interval, end.exp= end + gene_span_interval) %>% 
#   select(-c(feature_type,score,strand,phase,gene_source))
# fwrite(gtf2, "refdata-b37-2.1.0_gene_annotations.Filt.Genes.gtf")
gtf2 <- fread(paste0('refdata-b37-2.1.0_gene_annotations.Filt.Genes.gtf'), header = FALSE, sep = '\t')

#### Plotting with circlize ####
GOI.fusion.calls= AllSamples_AllCalls %>% filter(gene1 %in% GOI | gene2 %in% GOI) %>% 
  filter(gene2 %in% GOI) %>% select(FusionName, chr_1, break1, gene1, chr_2, break2, gene2) %>% distinct()
col = brewer.pal(length(unique(GOI.fusion.calls$FusionName)),"Paired")
col = cbind.data.frame(FusionName = unique(GOI.fusion.calls$FusionName), col)
col = col[length(unique(GOI.fusion.calls$FusionName)),] #for when there's <3 SV
GOI.fusion.calls = merge(GOI.fusion.calls, col, by = "FusionName")
GOI.fusion.calls = GOI.fusion.calls %>% mutate(break1 = as.numeric(break1), break1.1 = as.numeric(break1),
                                               break2 = as.numeric(break2), break2.1 = as.numeric(break2),
                                               gene1 = as.character(gene1), gene2 = as.character(gene2),
                                               chr_1 = as.character(chr_1), chr_2 = as.character(chr_2))

fusion.chr = unique(as.vector(cbind(GOI.fusion.calls$chr_1, GOI.fusion.calls$chr_2)))
fusion.gene = character()
fusion.gene = c(fusion.gene, na.omit(unique(GOI.fusion.calls$gene1)), na.omit(unique(GOI.fusion.calls$gene2)))
# do this so gene names are only plotted once with multiple breaks
bed = gtf2 %>% filter(gene_name %in% fusion.gene) %>% 
  mutate(chr = paste0("chr", chromosome)) %>% 
  filter(chr %in% fusion.chr) %>%
  select(chr, start, end, gene = gene_name) 
bed = as.data.frame(bed)
bed = transform(bed, start = as.numeric(as.character(start)), end = as.numeric(as.character(end)))


chr_cyto = fread("refdata-b37-2.1.0_cytoBand.txt") #use this to open files
colnames(chr_cyto)=c("chr", "start", "end", "cytoband", "stain")
chr_cyto = as.data.frame(chr_cyto)
chr_cyto["omosome"] = sub("chr", "", chr_cyto$chr)
chr_cyto = chr_cyto %>% filter(chr %in% fusion.chr) %>% mutate(omosome = as.numeric(omosome)) %>% arrange(omosome, start)

png(file=paste0("Fusion_filtering_Output/GOI_circular-diagram.png"), width=3000, height=2000, res=300)
circos.par(gap.after = 2, start.degree = 90)
circos.genomicInitialize(chr_cyto, plotType = NULL)
circos.genomicLabels(bed, labels.column = 4, side = "outside", cex = 1)
circos.track(ylim = c(0, 0.01), track.height = 0.05,
             bg.col = brewer.pal(length(fusion.chr),"Set2"), 
             panel.fun = function(x, y) circos.genomicAxis())
circos.genomicIdeogram(chr_cyto, track.height = convert_height(4, "mm"))

region1 = GOI.fusion.calls %>% select(chr_1, break1, break1.1, col) 
colnames(region1)[1:3] = c("chr", "start", "end")
region2 = GOI.fusion.calls %>% select(chr_2, break2, break2.1, col)
colnames(region2)[1:3] = c("chr", "start", "end")
if (nrow(region1)!=0) {
  for (i in 1:dim(region1)[1]) {
    circos.link(region1[i,1], region1[i,2], region2[i,1], region2[i,2], col = region1[i,4], lwd = 2)
  }
}
legend(1,0.15, title = "Fusions", legend = col$FusionName, bty = "n",
       col = col$col, lty = 1, lwd = 3, pt.cex = 2)
circos.clear()
dev.off()
#####