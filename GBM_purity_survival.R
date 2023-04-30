
#' ----------------
#'  load genes    |
#' ----------------
gen1 <- c('SDC1','SDC4','CD44','ITGAV','ITGB8')

#' ---------------------------------------------------- 
#'  TCGA survival 
#' ---------------------------------------------------- 
#' --- RNAseq data ---
#' @param rnaseq2 RNAseq data with matched clinical profile
cur.dir <- '.'
rnaseq <- read.delim(file='TCGA_GBM_HG-U133A',stringsAsFactors=F,row.names=1,check.names=F)

#' --- phenotype data ---
#' @import data.surv2 clinical data with matched RNAseq profile
surv <- read.delim(file='TCGA_GBM_clinicalMatrix.txt',stringsAsFactors=F,row.names=1,check.names=F)

ps1 <- which(colnames(surv)=='OS.time')
ps2 <- which(colnames(surv)=='OS')
ps.o <- grep('PFI',colnames(surv))
surv1 <- subset(surv,select=c(ps1,ps2,ps.o))
colnames(surv1) <- c('time','event','pfiE','PFI')

#' ------- tumor purity -------
puri <- read.csv(file='TCGA_tumor_purity.csv',stringsAsFactors=F,row.names=1,check.names=F)
puri1 <- puri[puri[,1]=='GBM',1:6]
tem <- sapply(rownames(puri1),function(i){
    t <- strsplit(i,'-')[[1]]
    t1 <- paste0(strsplit(t[4],'')[[1]][1:2],collapse='')
    t2 <- paste(t[1],t[2],t[3],sep='-')
    return (paste(t2,t1,sep='-'))})
puri1$sample <- tem
puri2 <- puri1[-which(is.na(puri1$ABSOLUTE)),]

#' patient samples in common 
coms <- intersect(intersect(colnames(rnaseq),rownames(surv1)),
                  puri2$sample)
df.rnaseq <- rnaseq[,coms]; df.surv <- surv1[coms,]
all(identical(colnames(df.rnaseq),rownames(df.surv)))

puri3 <- puri2[match(coms,puri2$sample),]
all(identical(colnames(df.rnaseq),puri3$sample))

#' calcuate new rnaseq considering tumor purity
new.rna <- df.rnaseq * puri3$ABSOLUTE

#' ------------------------------
#'  survival analysis in GBM    |
#' ------------------------------
gbm.df <- dfs
gbm.rna <- gbm.df[[1]]; gbm.surv <- gbm.df[[2]]

library(matrixStats)
gbm.rna1 <- colMedians(as.matrix(gbm.rna[rownames(gbm.rna)%in%gen1,]))
    
a <- as.numeric(gbm.rna1)
ps1 <- which(a>quantile(a,0.7))
ps2 <- which(a<quantile(a,0.3))
gbm.surv1 <- gbm.surv[c(ps1,ps2),]
gbm.surv1$level <- 'low'
gbm.surv1$level[1:length(ps1)] <- 'high'

#' ------- os analysis ----------
gbm.surv.1 <- gbm.surv1[which(gbm.surv1$time< 365*5),]
gbm.surv.1$time <- gbm.surv.1$time/365

library(survival)
sdf <- survdiff(Surv(time,event) ~ level, data = gbm.surv.1)
print (sdf)

library(ggsci); col1 <- c(pal_material('red')(10)[8], pal_material('blue')(10)[8])
library(ggfortify)
gbm.fit1 <- survfit(Surv(time,event) ~ level, data = gbm.surv.1)

library(survminer)
pdf('FiveYrs_OS_purity.pdf',height=4,width=4)
print (ggsurvplot(gbm.fit1,data = gbm.surv.1, 
           size = 1, palette = col1,
           pval = TRUE,
           ggtheme = theme_bw(
               base_size=12)) + theme_survminer(
     rect = element_rect(fill = "transparent") 
   ))
graphics.off()

#'  PFS curves 
gbm.surv.2 <- gbm.surv1[which(gbm.surv1$PFI< 365*5),]
gbm.surv.2$PFI <- gbm.surv.2$PFI/365

gbm.fit2 <- survfit(Surv(PFI,pfiE) ~ level, data = gbm.surv.2)
sdf <- survdiff(Surv(PFI,pfiE) ~ level, data = gbm.surv.2)
print (sdf)

library(survminer)
pdf('FiveYrs_PFS_purity.pdf',height=4,width=4)
print (ggsurvplot(gbm.fit2,data = gbm.surv.2, 
           size = 1, palette = col1,
           pval = TRUE,
           ggtheme = theme_bw(
               base_size=12)) + theme_survminer(
     rect = element_rect(fill = "transparent") 
   ))
graphics.off()
