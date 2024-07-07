#install necessary packages (only need to run once)
install.packages("BiocManager")
install.packages("readxl")
install.packages("readr")
BiocManager::install('clusterProfiler')
BiocManager::install('DOSE')
BiocManager::install('org.Hs.eg.db')


#load table with list of differential expressed genes

library(readxl)
de_genes <- read_excel("")

#load table with all genes detected to be expressed in the experiment
#this set of genes will the universe for the over-representation analysis 

library(readr)
detected_genes <- read_csv("detected_genes.csv")

univgenes=unique(detected_genes$symbol)

#perform over-representation analysis of GO BP terms

library(clusterProfiler)
library(org.Hs.eg.db)

ego <- enrichGO(gene          = de_genes$`Gene ID`,
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

ego@result$FoldEnrichment=parse_ratio(ego@result$GeneRatio)/parse_ratio(ego@result$BgRatio)


#to check the significant results in a table
egodf=as.data.frame(ego)

#to remove enriched terms with small number of associated genes
todrop=ego@result$ID[ego@result$Count<4]
ego=dropGO(ego,term=todrop)

egodf=as.data.frame(ego)

#to view results in a plot

dotplot(ego, x="FoldEnrichment",color="p.adjust",size="Count",showCategory=20,orderBy="Count")


#to export the results table

write.table(egodf,"egodf.txt",row.names=F)


#perform over-representation analysis of KEGG pathways

#this analysis requires using a different gene id
gene.df <- bitr(gene= de_genes$`Gene ID`, fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)

univ.df <- bitr(gene= univgenes, fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)

kk <- enrichKEGG(gene          = gene.df$ENTREZID,
                 organism     = 'hsa',
                 pAdjustMethod = "BH",
                 universe= univ.df$ENTREZID,
                 pvalueCutoff = 0.05)

kk@result$FoldEnrichment=parse_ratio(kk@result$GeneRatio)/parse_ratio(kk@result$BgRatio)

dotplot(kk, x="FoldEnrichment",color="p.adjust",size="Count",showCategory=20,orderBy="Count")


kkdf=as.data.frame(kk)

#view results on kegg map
browseKEGG(kk, 'hsa04919')


