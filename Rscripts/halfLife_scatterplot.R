library(dplyr)
setwd("C:/Users/Lenovo/Desktop/github_JiachenBu/output")


HL_result <- read.table('ES_RNAseq.Martre1OEvsCtrl.normalized_counts.halflife.result.mse0.05.txt', h=T, row.names=1)
HL <- HL_result[, c("Ctrl.halflife", "Martre1OE.halflife") ]
x <- log2(HL+1)
colnames(x) <- c("Ctrl", "LocOE")
nrow(x)  #11,829

#######################
library(MASS)
library(ggplot2)
library(viridis)

theme_set(theme_bw(base_size = 16))
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

set.seed(1)
x$density <- get_density(x$Ctrl,x$LocOE, n = 100)
ggplot(x) + geom_point(aes(Ctrl, LocOE, color = density)) + scale_color_viridis() + xlim(0,10) + ylim(0,10) + geom_abline(intercept = 0, slope=1, linetype="dashed", color="red", size = 1.5) + 
labs(x = "mRNA t1/2 Ctrl(h, log2)", y= "mRNA t1/2 Martre1 OE(h, log2)") +
annotate("text", x=7.5, y=1, label= "n = 11,829")
ggsave("ES_RNAseq.Martre1OEvsCtrl.normalized_counts.halflife.result.mse0.05.heatscatterplot.pdf", width=5.5, height=4)
