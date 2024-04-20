

####################
#### step0: 
## Prepare htseq-count output files of S2 mapping in directory "S2/", each file were named as "sample_rep[1-3].unstranded.txt (for example: CtrlAct0h_rep1.unstranded.txt)".
## Prepare htseq-count output files of mm10 mapping in directory "GRCm38/", each file were named as "sample_rep[1-3].unstranded.txt (for example: CtrlAct0h_rep1.unstranded.txt)".

####################
#### step1:
## DESeq2 normalization
library(DESeq2)
library(apeglm)
library(tidyverse)

#Store "S2/" and "GRCm38/" files in the current directory
setwd("C:/Users/Lenovo/Desktop/github_JiachenBu/test_data")  

#S2 counts
dat.merged <- c()
for(file in dir(path="S2/", pattern="unstranded.txt", full.names=TRUE)) {
     x <- read.table(file, header=FALSE, sep="\t")
     colnames(x) <- c("gene_id", gsub(".*/|.txt|.unstranded.txt|S2_no24h/", "", file))
     if(is.null(dat.merged)) dat.merged <- x else dat.merged <- merge(x, dat.merged, by="gene_id")
}
spikeInCounts <- dat.merged[,!(colnames(dat.merged) %in% c("gene_id","gene_symbol", "exonLength"))]
rownames(spikeInCounts) <- dat.merged[,"gene_id"]

samp_labels <- colnames(spikeInCounts)
label <- colnames(spikeInCounts)
sample <- data.frame(condition=as.factor(sub("_rep[1-3]", "", label)), label=as.factor(label))
rownames(sample) <- label

if(!all(rownames(sample) == colnames(spikeInCounts))) {
    stop("Must keep consistent between the row names of samples and the column names of rawCount matrix")
 }
design="~condition"
design <- as.formula(design)
deseq.obj <- DESeqDataSetFromMatrix(countData = spikeInCounts, colData = sample, design = design)
preFilterFactor=0
deseq.obj <- deseq.obj[rowSums(counts(deseq.obj)) >= preFilterFactor,]   #preFilter             
deseq.obj <- DESeq(deseq.obj)  
sizeFactor <-sizeFactors(deseq.obj)  #sizeFactors of spike-in
sizeFactor

#GRCm38 counts 
dat.merged2 <- c()
for(file in dir(path="GRCm38/", pattern="unstranded.txt", full.names=TRUE)) {
     x <- read.table(file, header=FALSE, sep="\t")
     colnames(x) <- c("gene_id", gsub(".*/|.txt|.unstranded.txt|GRCm38_no24h/", "", file))
     if(is.null(dat.merged2)) dat.merged2 <- x else dat.merged2 <- merge(x, dat.merged2, by="gene_id")
}
rawCounts <- dat.merged2[,!(colnames(dat.merged2) %in% c("gene_id","gene_symbol", "exonLength"))]
rownames(rawCounts) <- dat.merged2[,"gene_id"]
rawCounts <- round(rawCounts,0)


#### spike-in normalized DE analysis
deseq.obj2 <- DESeqDataSetFromMatrix(countData = rawCounts, colData = sample, design = design)
preFilterFactor=0
deseq.obj2 <- deseq.obj2[rowSums(counts(deseq.obj2)) >= preFilterFactor,]
sizeFactors(deseq.obj2) <- sizeFactor 

#### spike-in normalized RPM
normalized_counts <- counts(deseq.obj2, normalized=TRUE) 
normalized_count <- data.frame(gene_id=rownames(normalized_counts), normalized_counts)


### output
annotation <-read.table("Z:/YJ/ES_RNAseq/Loc/gencode.vM20.annotation.txt", h=T)   ##annotation file 
colnames(annotation) <- c("gene_id", "gene_symbol", "gene_type")
result <- merge(annotation, normalized_count, by="gene_id")
write.table(result, file="DESeq2.S2normalized.normalized_counts.txt", quote=F, sep="\t", row.names=F)


####################
#### step2:
## half-life calculation

mse <- function(residual) {
	sum(residual**2)
}

nls.fit <- function(abundance, time) {
	data <- data.frame(abundance=abundance, time=time)
	fit <- nls(abundance ~ exp(-kdecay*time), data = data, start=list(kdecay=0.01), algorithm='plinear')
	kdecay <- coef(fit)['kdecay']
	t.halflife <- log(2) / kdecay
	names(t.halflife) <- 't.halflife'
	list(fit=fit, kdecay=kdecay, t.halflife=t.halflife, mse=mse(summary(fit)$residuals))
}

nls.plot <- function(data, fit, geneSymbol=NA) {
	t <- seq(from=min(data$time), to=max(data$time), len=20)
	expr <- try(sprintf('y=%.2f*exp(-%.2f*x), halflife=%.2fh', coef(fit)[2], coef(fit)[1], log(2)/coef(fit)[1]))
	if(is.na(geneSymbol)) {
		plot(data$time, data$abundance, xlab='Time(h)', ylab='Gene expression remaining')
		lines(t, predict(fit, list(time=t)), col='red', lty='dashed')
		text(3, 1, expr)
	} else {
		plot(data$time, data$abundance, xlab='Time(h)', ylab='Gene expression remaining')
		lines(t, predict(fit, list(time=t)), col='red', lty='dashed')
		title(geneSymbol)
		text(3, 1, expr)
	}
}

nls.genePlot <- function(geneName, geneDict=NA) {
	if(is.na(geneDict)) {
		geneDict <- 1:nrow(x)
		as.matrix(x[,'gene_symbol']) -> names(geneDict)
	}
	index <- geneDict[geneName]
	nsfit <- nls.fit(gene.expr[index,]/gene.expr[index,1], time)
	nls.plot(data.frame(abundance=gene.expr[index,]/gene.expr[index,1], time=time), nsfit$fit, geneSymbol=geneName)
}
#nls.genePlot('ABHD15')

x <- result

#preFilter of input matrix 
y <- x[,(colnames(x) %in% c("CtrlAct0h_rep1","CtrlAct0h_rep2", "Martre1OEAct0h_rep1", "Martre1OEAct0h_rep2"))]
preFilterFactor=200
sub_x <- x[rowSums(y) > preFilterFactor,]  

#Calculate Ctrl and OE respectively
pattern <- 'CtrlAct'
x.subset <- sub_x[,grepl(pattern, colnames(sub_x))]

subpattern <- c()
for(name in colnames(x.subset)) subpattern <- c(subpattern, unlist(strsplit(name, '_'))[1])
subpattern <- factor(subpattern)
gene.expr <- c()
for(sp in levels(subpattern)) {
	xt <- x.subset[,grepl(sp, colnames(x.subset))]
	gene.expr  <- cbind(gene.expr, rowMeans(xt))
}
colnames(gene.expr) <- levels(subpattern)
time <- c()
for(n in colnames(gene.expr)) time <- c(time, as.double(unlist(strsplit(n, paste0('^',pattern,'|h$')))[2]))

result <- c()
for(i in 1:nrow(x.subset)) {
	if(sd(gene.expr[i,]) <= 0) {
		result <- c(result, list(list(NA)))
	} else {
		msg <- try(list(append(list(gene_symbol=as.matrix(sub_x[i,'gene_id'])[1,1]), nls.fit(gene.expr[i,]/gene.expr[i,1], time))), silent=T)
		if(is(msg)[1] != 'try-error') result <- c(result, msg)
	}
}

output <- c()
for(i in 1:length(result)) {
	output <- rbind(output, try(
		c(gene_symbol=result[[i]]$gene_symbol, result[[i]]$t.halflife, coef(result[[i]]$fit), result[[i]]$mse)
	))
}
colnames(output) <- c('gene_id', 'halflife', 'kdecay', 'constant', 'mse')
write.table(output,'ES_RNAseq.normalized_counts.halflife.Ctrl.txt', sep="\t",quote=F, row.names=F ) 


pattern <- 'Martre1OEAct'
x.subset <- sub_x[,grepl(pattern, colnames(sub_x))]
subpattern <- c()
for(name in colnames(x.subset)) subpattern <- c(subpattern, unlist(strsplit(name, '_'))[1])
subpattern <- factor(subpattern)
gene.expr <- c()
for(sp in levels(subpattern)) {
	xt <- x.subset[,grepl(sp, colnames(x.subset))]
	gene.expr  <- cbind(gene.expr, rowMeans(xt))
}
colnames(gene.expr) <- levels(subpattern)
time <- c()
for(n in colnames(gene.expr)) time <- c(time, as.double(unlist(strsplit(n, paste0('^',pattern,'|h$')))[2]))

result <- c()
for(i in 1:nrow(x.subset)) {
	if(sd(gene.expr[i,]) <= 0) {
		result <- c(result, list(list(NA)))
	} else {
		msg <- try(list(append(list(gene_symbol=as.matrix(sub_x[i,'gene_id'])[1,1]), nls.fit(gene.expr[i,]/gene.expr[i,1], time))), silent=T)
		if(is(msg)[1] != 'try-error') result <- c(result, msg)
	}
}

output <- c()
for(i in 1:length(result)) {
	output <- rbind(output, try(
		c(gene_symbol=result[[i]]$gene_symbol, result[[i]]$t.halflife, coef(result[[i]]$fit), result[[i]]$mse)
	))
}
colnames(output) <- c('gene_id', 'halflife', 'kdecay', 'constant', 'mse')
write.table(output,'ES_RNAseq.normalized_counts.halflife.Martre1OE.txt', sep="\t",quote=F, row.names=F )


####################
#### step3:
## plot 

library(dplyr)

# Merge data 
Ctrl <- read.table("ES_RNAseq.normalized_counts.halflife.Ctrl.txt",h=T)
Ctrl_f <- filter(Ctrl,Ctrl$mse < 0.05 & Ctrl$halflife > 0)

OE <- read.table("ES_RNAseq.normalized_counts.halflife.Martre1OE.txt",h=T)
OE_f <- filter(OE,OE$mse < 0.05 & OE$halflife > 0)

res <- merge(Ctrl_f, OE_f, by="gene_id")
colnames(res)<- c('gene_id', 'Ctrl.halflife', 'Ctrl.kdecay', 'Ctrl.constant', 'Ctrl.mse', 'Martre1OE.halflife', 'Martre1OE.kdecay', 'Martre1OE.constant', 'Martre1OE.mse')
annotation <- read.table("gencode.vM20.annotation.txt",h=T)
colnames(annotation) <- c("gene_id", "gene_symbol", "gene_type")
result <- merge(annotation, res, by="gene_id")

write.table(result, "ES_RNAseq.Martre1OEvsCtrl.normalized_counts.halflife.result.mse0.05.txt", quote=F, sep="\t", row.names=F) 

# box_plot
library(ggplot2)
library(ggpubr)

HL_result <- result
HL <- HL_result[, c("Ctrl.halflife", "Martre1OE.halflife") ]
x <- log2(HL+1)
colnames(x) <- c("Ctrl", "Martre1OE")

df <- c()
for(cn in colnames(x)) {
  strs <- unlist(strsplit(cn, '_rep'))
	df <- rbind(
	  df,
		data.frame(gene = rownames(x), sample = strs[1], replicate = paste0('rep', strs[2]), expression = x[,cn])
	)
}


p <- ggplot(df, aes(x = sample, y = expression)) + 
geom_boxplot(outlier.colour = 'transparent', fill = "lightblue", colour="black", linetype="dashed") + 
stat_boxplot(aes(ymin=..lower..,ymax=..upper..), fill = "lightblue", outlier.colour = 'transparent') +
stat_boxplot(geom = "errorbar",aes(ymin=..ymax..),width=0.2, color="black") +
stat_boxplot(geom = "errorbar",aes(ymax=..ymin..),width=0.2, color="black") +
theme_classic()  + 
labs(x = NULL, y= "Log2 halflife(h)") + 
theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
theme(axis.text = element_text(colour = "black")) +
theme(axis.title.y = element_text(size = rel(1), angle = 90))
p + stat_compare_means() + ylim(0, 8)
ggsave("ES_RNAseq.Martre1OEvsCtrl.normalized_counts.halflife.result.mse0.05.boxplot.pdf", width = 2.5, height = 6)

median(x$Ctrl)
median(x$Martre1OE)

