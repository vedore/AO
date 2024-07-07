library(readxl)
diff_genes <- read_excel("GSE39210.top.table2.xlsx")

#only logFold abs > 3.5

diffgenes_up <- diff_genes[diff_genes$logFC > 3.5,]
diffgenes_down <- diff_genes[diff_genes$logFC > - 3.5,]
diffgenes_down <- diffgenes_down[diffgenes_down$logFC < 0,]

#only p < 0.05

diffgenes_up <- diffgenes_up[diffgenes_up$adj.P.Val < 0.05,]
diffgenes_down <- diffgenes_down[diffgenes_down$adj.P.Val <  0.05,]

#remove null symbols

diffgenes_up <- diffgenes_up[!is.na(diffgenes_up$Gene.symbol) & !grepl("^\\s*$", diffgenes_up$Gene.symbol), ]
diffgenes_down <- diffgenes_down[!is.na(diffgenes_down$Gene.symbol) & !grepl("^\\s*$", diffgenes_down$Gene.symbol), ]

#split and remove repeated reads of the same gene (in different probes)
diffgenesup=unlist(strsplit(diffgenes_up$Gene.symbol,split="///",fixed=T))
diffgenesup=unique(diffgenesup)

diffgenesdown=unlist(strsplit(diffgenes_down$Gene.symbol,split="///",fixed=T))
diffgenesdown=unique(diffgenesdown)

#define universe

univgenes=unique(allgenes)

##up

library(clusterProfiler)
library(org.Hs.eg.db)

ego_up <- enrichGO(gene          = diffgenesup,
                universe      = univgenes,
                OrgDb         = org.Hs.eg.db,
                keyType = "SYMBOL",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.2,
                qvalueCutoff  = 1,
                minGSSize = 10,
                maxGSSize = 500,
                readable      = TRUE)

library(DOSE)

ego_up@result$FoldEnrichment=parse_ratio(ego_up@result$GeneRatio)/parse_ratio(ego_up@result$BgRatio)

#to check the significant results in a table
ego_updf=as.data.frame(ego_up)

#to remove enriched terms with small number of associated genes
todrop=ego_up@result$ID[ego_up@result$Count<4]
ego_up=dropGO(ego_up,term=todrop)

ego_updf=as.data.frame(ego_up)

#to view results in a plot

dotplot(ego_up, x="FoldEnrichment",color="p.adjust",size="Count",showCategory=20,orderBy="Count")


##down

library(clusterProfiler)
library(org.Hs.eg.db)

ego_down <- enrichGO(gene          = diffgenesdown,
                   universe      = univgenes,
                   OrgDb         = org.Hs.eg.db,
                   keyType = "SYMBOL",
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 1,
                   minGSSize = 10,
                   maxGSSize = 500,
                   readable      = TRUE)

library(DOSE)

ego_down@result$FoldEnrichment=parse_ratio(ego_down@result$GeneRatio)/parse_ratio(ego_down@result$BgRatio)

#to check the significant results in a table
ego_downdf=as.data.frame(ego_down)

#to remove enriched terms with small number of associated genes
todrop=ego_down@result$ID[ego_down@result$Count<4]
ego_down=dropGO(ego_down,term=todrop)

ego_downdf=as.data.frame(ego_down)

#to view results in a plot

dotplot(ego_down, x="FoldEnrichment",color="p.adjust",size="Count",showCategory=20,orderBy="Count")

