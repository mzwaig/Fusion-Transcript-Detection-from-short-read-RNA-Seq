#setwd("")
library(dplyr)
library(tidyr)
library(ggplot2)
require(data.table) ## 1.9.4+
library(openxlsx)
library(circlize)
library(RColorBrewer)

datasets = read.delim(file= paste0("datasets.tsv"))

load.calls <- function(x){
  # load STAR-Fusion ####
  print(paste0(x$Sample, " - STAR-Fusion"))
  star.fusion = read.delim(file=paste0(x$STAR.Fusion))
  colnames(star.fusion)[1]="FusionName"
  star.fusion = star.fusion %>% 
    separate(col = LeftGene, into = c("LeftGene", "LeftGene_ID"), sep = "\\^") %>% 
    separate(col = RightGene, into = c("RightGene", "RightGene_ID"), sep = "\\^") %>% 
    separate(col = LeftBreakpoint, into = c("LeftCHROM", "LeftBreak","LeftOrientation"), sep = "\\:") %>% 
    separate(col = RightBreakpoint, into = c("RightCHROM", "RightBreak","RightOrientation"), sep = "\\:") %>%
    mutate(Sample=paste0(x$Sample), Caller="STAR.Fusion")
  # Give files common headers ####
  common.star.fusion = star.fusion %>% select(c(Sample = Sample, FusionName=FusionName,  Caller=Caller, chr_1=LeftCHROM, break1=LeftBreak, gene1=LeftGene,
                                                strand1=LeftOrientation, chr_2=RightCHROM, break2=RightBreak, gene2=RightGene,
                                                strand2=RightOrientation, SplitReads=JunctionReadCount,
                                                DiscordantMates=SpanningFragCount, FusionType=PROT_FUSION_TYPE,
                                                Feature1=CDS_LEFT_ID, Feature2=CDS_RIGHT_ID))
    # load Arriba ####
  print(paste0(x$Sample, " - Arriba"))
  arriba = read.delim(file=paste0(x$Arriba))
  colnames(arriba)[1]="gene1"
  arriba = arriba %>% 
    separate(col = breakpoint1, into = c("Chrom1", "Breakpoint1"), sep = "\\:") %>% 
    separate(col = breakpoint2, into = c("Chrom2", "Breakpoint2"), sep = "\\:") %>% 
    mutate(Sample=paste0(x$Sample), Chrom1=paste0("chr",Chrom1), Chrom2=paste0("chr",Chrom2),
           FusionName=paste0(gene1,"--",gene2), split_reads=as.numeric(split_reads1)+as.numeric(split_reads2), Caller= "Arriba")
  # Give files common headers ####
  common.arriba = arriba %>% select(c(Sample=Sample, FusionName=FusionName, Caller=Caller,chr_1=Chrom1, break1=Breakpoint1, gene1=gene1,
                                      strand1=strand1.gene.fusion., chr_2=Chrom2, break2=Breakpoint2, gene2=gene2,
                                      strand2=strand2.gene.fusion., SplitReads=split_reads,
                                      DiscordantMates=discordant_mates, FusionType=type, Feature1=site1,
                                      Feature2=site2))
  
  # load InFusion ####
  print(paste0(x$Sample, " - InFusion"))
  infusion =  read.delim(file=paste0(x$InFusion)) #pre-filtered for only gene-gene fusions
  colnames(infusion)[1]="ID"
  infusion["break_on_exon"] = gsub("yes", "ExonBreak", infusion$break_on_exon)
  infusion = infusion %>% mutate(Sample=paste0(x$Sample), ref1=paste0("chr",ref1), ref2=paste0("chr",ref2), FusionName=paste0(gene_1,"--",gene_2))
  infusion_stringent =  read.delim(paste0(x$InFusion_stringent)) #pre-filtered for only gene-gene fusions
  colnames(infusion_stringent)[1]="ID"
  infusion_stringent = infusion_stringent %>% mutate(Sample=paste0(x$Sample))
  infusion_stringent = unique(infusion_stringent$ID)
  infusion = infusion %>% mutate(stringent= ifelse(Sample %in% infusion_stringent, "stringent", NA), Caller="InFusion")
  # Give files common headers ####
  common.infusion = infusion %>% select(c(Sample = Sample, FusionName=FusionName, Caller=Caller,chr_1=ref1, break1=break_pos1, gene1=gene_1,
                                          strand1=gene_1_strand, chr_2=ref2, break2=break_pos2, gene2=gene_2,
                                          strand2=gene_2_strand, SplitReads=num_split, DiscordantMates=num_paired, 
                                          FusionType=break_on_exon, Feature1=feature_1, Feature2=feature_2))
  #Create masterlist of all fusion calls by all callers
  masterlist = as.data.frame(rbind(common.star.fusion, common.arriba, common.infusion))
  masterlist = masterlist %>% select(Sample, Caller, FusionName, everything())
  return(masterlist)
}


AllSamples_AllCalls=NULL
for(i in 1:nrow(datasets)) {
  data = datasets[i,]
  print(data$Sample)
  dir.create("Fusion_filtering_Output", showWarnings = FALSE) #stops warnings if folder already exists
  #Call load.calls function ####
  masterlist = load.calls(data)
  fwrite(masterlist, paste0("Fusion_filtering_Output/", data$Sample,"_AllCallers.tsv"), sep="\t") #Creates an output per sample
  
  #Create masterlist across all samples
  AllSamples_AllCalls = as.data.frame(rbind(AllSamples_AllCalls,masterlist)) #all sample to cross-sample masterlist
  }
fwrite(AllSamples_AllCalls, "Fusion_filtering_Output/AllSamples_AllCallers.noFiltering.tsv", sep="\t") #Create ouput across all samples

#Look for fusions called more than once
recurrent_fusions = AllSamples_AllCalls %>% select(Sample, Caller, FusionName) %>% distinct() %>%
  group_by(Sample, FusionName) %>% summarise(count=n()) %>% filter(count>=2)
recurrent_fusions = as.data.frame(recurrent_fusions)

#Generate summaries 
breakpoint_summary = AllSamples_AllCalls %>% group_by(chr_1, break1, gene1, chr_2, break2, gene2) %>% summarise(count=n()) 
breakpoint_summary = as.data.frame(breakpoint_summary)
summary_by_sample = AllSamples_AllCalls %>% group_by(Sample, chr_1, break1, gene1, chr_2, break2, gene2) %>% 
  summarise(count=n(),meanSplitReads=round(mean(SplitReads),digits=1), meanDiscordantMates = round(mean(DiscordantMates),digits=1))
summary_by_fusion = AllSamples_AllCalls %>% group_by(FusionName, chr_1, break1, gene1, chr_2, break2, gene2) %>% summarise(count=n()) 

#Pull out genes of interest 
GOI = c("BRAF", "RAF1") 
GOI.fusion.calls= AllSamples_AllCalls %>% filter(gene1 %in% GOI | gene2 %in% GOI)


