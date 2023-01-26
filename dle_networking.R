# Title     : TODO
# Objective : TODO
# Created by: coltongarelli
# Created on: 7/15/20
library(clusterProfiler)
library(Biobase)
library(org.Hs.eg.db)




e_ctcl <- read.csv('/Users/coltongarelli/Desktop/GSE113113_stage1_DE.csv')
filt <- c(
    'gene', 'logFC', 'adj.P.Val')
l_ctcl <- read.csv('/Users/coltongarelli/Desktop/GSE113113_stage2_DE.csv')
dle <- read.csv('/Users/coltongarelli/Desktop/GSE95474 DLE diffex collapsed probes.csv')
e_ctcl <- e_ctcl[ ,filt]
l_ctcl <- l_ctcl[ ,filt]
dle <- dle[ ,filt]
dle <- dle[abs(dle$logFC) > 1.0,]
dle <- dle[dle$adj.P.Val < 0.05,]

e_ctcl <- e_ctcl[abs(e_ctcl$logFC) > 1.0,]
e_ctcl <- e_ctcl[e_ctcl$adj.P.Val < 0.05,]

l_ctcl <- l_ctcl[abs(l_ctcl$logFC) > 1.0,]
l_ctcl <- l_ctcl[l_ctcl$adj.P.Val < 0.05,]




dle_simp <- dle[c('gene', 'logFC')]
ec_simp <- e_ctcl[c('gene', 'logFC')]
lc_simp <- l_ctcl[c('gene', 'logFC')]






rownames(dle_simp) <- dle_simp$gene
# dle_simp$gene <- NULL
dle_simp <- dle_simp[order(dle_simp$logFC, decreasing=TRUE),]
# dle_simp

rownames(ec_simp) <- ec_simp$gene
# ec_simp$gene <- NULL
ec_simp <- ec_simp[order(ec_simp$logFC, decreasing=TRUE),]

rownames(lc_simp) <- lc_simp$gene
# lc_simp$gene <- NULL
lc_simp <- lc_simp[order(lc_simp$logFC, decreasing=TRUE),]




# IDs
dle_symbols <- clusterProfiler::bitr(rownames(dle_simp), fromType='ALIAS', toType='ENTREZID', OrgDb='org.Hs.eg.db')
dle_symbols$gene <- dle_symbols$ALIAS
dle_symbols$ALIAS <- NULL
dle_entrez <- merge(dle_simp, dle_symbols, by.y='gene')




options(max.print=1000)


dle_agg <- aggregate(dle_simp$logFC,
          by = list(gene = dle_simp$gene),
          FUN = mean,
          na.rm = TRUE)
# dle_agg <- dle_agg[order(dle_agg$x, decreasing=TRUE),]
# dle_agg
# dle_entrez
# rownames(dle_agg) <- dle_agg$entrez
# dle_agg$entrez <- NULL
rownames(dle_simp)<-dle_simp$gene
dle_simp$gene <- NULL
dle_vec <- unlist(dle_simp)
names(dle_vec) <- rownames(dle_simp)
# dle_agg
dle_gsego <- clusterProfiler::gseGO(geneList= dle_vec,
              OrgDb        = 'org.Hs.eg.db',
              ont          = "BP",
              nPerm        = 1000,
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE,
              keyType      = 'SYMBOL')










dle_gsego


# pdf(file='/Users/coltongarelli/Desktop/dle_go_graph_GSE95474.pdf', height=15, width=15)
# dle_gog <- clusterProfiler::plotGOgraph(dle_gsego[100:nrow(dle_gsego),], firstSigNodes=100)
# print(dle_gog)
# dev.off()




pdf(file='/Users/coltongarelli/Desktop/dle_circ_cnet_fc_GSE95474.pdf', height=15, width=15)
dle_gog <- clusterProfiler::cnetplot(dle_gsego, foldChange=dle_vec, showCategory=10, colorEdge=TRUE)
print(dle_gog)
dev.off()




# install.packages("devtools")
devtools::install_github("YuLab-SMU/clusterProfiler.dplyr", repos='http://cran.us.r-project.org')
library(dplyr)
y <- dplyr::arrange(dle_gsego, abs(NES)) %>%
        group_by(sign(NES)) %>%
        slice(1:5)
clusterProfiler::cnetplot(y)





dle_kegg <- clusterProfiler::compareCluster(rownames(dle_simp), fun="enrichKEGG", organism='hsa', pvalueCutoff=0.05)
clusterProfiler::emapplot(dle_kegg)




dle_symbols <- clusterProfiler::bitr(rownames(dle_simp), fromType='ALIAS', toType='KEGGID', OrgDb='org.Hs.eg.db')

# clusterProfiler::emapplot(dle_kegg)





library(DOSE)
pdf(file='/Users/coltongarelli/Desktop/dle_emap_fc_GSE95474.pdf', height=15, width=15)
dle_gog <- clusterProfiler::emapplot(dle_gsego[100:,])
print(dle_gog)
dev.off()




ec_gsego <- clusterProfiler::gseGO(geneList= ec_simp,
              OrgDb        = 'org.Hs.eg.db',
              ont          = "BP",
              nPerm        = 1000,
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)




lc_gsego <- clusterProfiler::gseGO(geneList= lc_simp,
              OrgDb        = 'org.Hs.eg.db',
              ont          = "BP",
              nPerm        = 1000,
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)




