devtools::install_version('BiocManager', version = '3.10')

library(BiocManager)
BiocManager::install('oligo')
BiocManager::install('affy')
BiocManager::install('sva')
BiocManager::install('limma')
BiocManager::install('RSQLite')
# BiocManager::install('annotate')
# BiocManager::install('biomaRt')
library('hugene21stprobeset.db')
library('hugene21sttranscriptcluster.db')
library(oligo)
library(limma)
library(sva)
wp <-"/Users/coltongarelli/Desktop/DLE paper data/data/human/files/raw files/GSE81071_RAW/"

# library('biomaRt')

# Load probeset information
fd<-read.AnnotatedDataFrame("/Users/coltongarelli/Desktop/DLE paper data/data/human/files/raw files/GSE81071_RAW/GPL19983_hugene21st_Hs_ENTREZG_desc.txt.gz",header=TRUE ,sep="\t")
# Load the data
data <- oligo::read.celfiles(filenames=affy::list.celfiles('/Users/coltongarelli/Desktop/DLE paper data/data/human/files/raw files/GSE81071_RAW/', full.names = TRUE))
oligo::availProbeInfo(data)
# add featureData from file to data object 


# find probe info @ https://www.bioconductor.org/packages/release/data/annotation/
# id_map <- data.frame(affy(hugene21stprobeset.db, keys = data@featureData@varMetadata, keytype = 'PROBEID', column = 'SYMBOL'))
# fData(data)$altnames <- out
# Define the groups (disease vs healthy or in this case two diseases vs healthy)
names <- ifelse(grepl('.*?_DLE', sampleNames(data)), 'DLE',
                ifelse(grepl('.*?_SCLE', sampleNames(data)), 'SCLE', 'Control'))
pData(data)$groups <- names
# Make unique sample names from file names
samples <- sub('.*?_', '', sampleNames(data))
samples <- sub('.*?_', '', samples)
samples <- sub('*.CEL.gz', '', samples)
sampleNames(data) <- samples

# Translate Affymetrix probes into entrez/gene names (310819 total na)
out <- AnnotationDbi::select(hugene21stprobeset.db, keys =oligo::probeNames(data),c('SYMBOL', 'GENENAME', 'ENTREZID'))
# remove any annotations that return an NA )
annotations <- out[!is.na(out['SYMBOL']),]
annotationsVec = c(annotations$ENTREZID)
names(annotationsVec) = annotations$SYMBOL
annotationsVec = as.list(annotationsVec)
d = data[rownames(data@featureData) %in% annotations$SYMBOL,]
featd <-combine(fd,  AnnotatedDataFrame(out))
data@featureData <- fd

# Robust microarray analysis to correct background and normalize
data.rma <- oligo::rma(data)
pheno <- pData(data)

# Make model matrix for batch correction
pheno <- pData(data)
batch1 = c(rep(1, 56))
batch2 = c(rep(2, 47))
batches <- c(batch1, batch2)
gs <- as.factor(pData(data)$groups)
bs <- as.factor(batches)
mod <- model.matrix(~as.factor(pData(data)$groups), data = pheno)
mod0 = model.matrix(~1, data=pheno)
affy::mva.pairs(exprs(data.rma)[,1:4],log.it = FALSE,plot.method="smoothScatter")

# Correct for batch effects
# https://www.coursera.org/lecture/statistical-genomics/batch-effects-in-r-part-a-8-18-zsDd4
data.normed_batched <- sva::ComBat(exprs(data.rma), batch = batches, mod = mod)
affy::mva.pairs(exprs(data.rma)[,1:4],log.it = FALSE,plot.method="smoothScatter")
affy::mva.pairs(data.normed_batched[,1:4],log.it = FALSE,plot.method="smoothScatter")

# Generate a linear model used to determine true DEGs
# https://medium.com/biosyntax/using-limma-to-find-differentially-expressed-genes-1c12c905e3b1
fit <- limma::lmFit(data.normed_batched, mod)
bay <- limma::eBayes(fit) 
data.de <- limma::topTable(bay, adjust.method = 'BH', sort.by = 'none', number = Inf)




# color=c('green','green','green','red','red','red','blue','blue','blue')
# data.PC = prcomp(t(data.matrix),scale.=TRUE)
# plot(data.PC$x[1:2],col=color)

