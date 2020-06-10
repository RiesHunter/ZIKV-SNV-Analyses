#### Session prep ####
## Clear Global Environment
rm(list = ls())

## Install packages and load libraries as required
if(!require(tidyverse)){
  install.packages("tidyverse",dependencies = T)
  library(tidyverse)
}
if(!require(ggplot2)){
  install.packages("ggplot2",dependencies = T)
  library(ggplot2)
}
if(!require(gridExtra)){
  install.packages("gridExtra",dependencies = T)
  library(gridExtra)
}
if(!require(grid)){
  install.packages("grid",dependencies = T)
  library(grid)
}

#### Mouse Passage 10 ####
## Create list of all .vcf files for Passage Number 10 
setwd("~/Documents/Research/1-Friedrich/1-ZIKV-SNVs/Kasen-Tutorial/Datasets/SCP_refaligned_vcf")
mo.myvcf10 <- dir(pattern="SCP10_")    
mo.n10 <- length(mo.myvcf10)             
mo.mylist10 <- vector("list",mo.n10)     

## Read all the tables in myvcf and apply to mylist10
for (i in 1:mo.n10) {mo.mylist10[[i]] <- read.table(mo.myvcf10[i])} 

## Split INFO into 4 columns, respectively
for (i in 1:mo.n10) {mo.mylist10[[i]] <- mo.mylist10[[i]] %>% 
  separate("V8",c("V8_1","V8_2"),sep=";",convert=FALSE) %>% 
  separate("V8_1",c("V8_1_1","V8_1_2"),sep="=",convert=FALSE) %>% 
  separate("V8_2",c("V8_2_1","V8_2_2"),sep="=",convert=FALSE)
}

## Remove unnecessary columns and rename remaining  
for (i in 1:mo.n10) {mo.mylist10[[i]] <- mo.mylist10[[i]][,!(colnames(mo.mylist10[[i]]) %in% c("V1","V3","V6","V7","V8_1_1","V8_1_2","V8_2_1"))]}
for (i in 1:mo.n10) {mo.mylist10[[i]] <- dplyr::rename(mo.mylist10[[i]],CHROM=V2,POS=V2,REF=V4,ALT=V5,AF=V8_2_2)}

## Create "data" data.frame for all mylist10
for (i in 1:mo.n10) {mo.mylist10[[i]]$AF <- as.numeric(mo.mylist10[[i]]$AF)}
for (i in 1:mo.n10) {mo.mylist10[[i]]$AF <- mo.mylist10[[i]]$AF*100} #As % frequency
for (i in 1:mo.n10) {mo.mylist10[[i]]$SUB <- paste(mo.mylist10[[i]]$POS,mo.mylist10[[i]]$REF,">",mo.mylist10[[i]]$ALT)}
for (i in 1:mo.n10) {mo.mylist10[[i]]$NUC <- paste(mo.mylist10[[i]]$POS,mo.mylist10[[i]]$ALT)}
for (i in 1:mo.n10) {mo.mylist10[[i]]$SUB <- gsub(" ", "",mo.mylist10[[i]]$SUB,fixed=TRUE)}
for (i in 1:mo.n10) {mo.mylist10[[i]] <- mo.mylist10[[i]][,!(colnames(mo.mylist10[[i]]) %in% c("REF","column_label"))]}

## Filter out all barcodes with AF <= .3%
for (i in 1:mo.n10) {mo.mylist10[[i]] <- filter(mo.mylist10[[i]],AF>.3)}

## Give each list its own data frame
mo.data10.1 <- mo.mylist10[[1]] %>% 
  add_column(PN=10) %>% 
  add_column(REP=1)
mo.data10.2 <- mo.mylist10[[2]] %>% 
  add_column(PN=10) %>% 
  add_column(REP=2)
mo.data10.3 <- mo.mylist10[[3]] %>% 
  add_column(PN=10) %>% 
  add_column(REP=3)
mo.data10.4 <- mo.mylist10[[4]] %>% 
  add_column(PN=10) %>% 
  add_column(REP=4)
mo.data10.5 <- mo.mylist10[[5]] %>% 
  add_column(PN=10) %>% 
  add_column(REP=5)
mo.data10 <- rbind(mo.data10.1,mo.data10.2,mo.data10.3,mo.data10.4,mo.data10.5)

mo.data10.factor <- mo.data10
mo.data10.factor$REP <- as.factor(mo.data10$REP)
levels(mo.data10.factor$REP) <- c("Lineage A", "Lineage B", "Lineage C", "Lineage D", "Lineage E")
mo.data10.factor$PN <- as.factor(mo.data10$PN)
levels(mo.data10.factor$PN) <- c("Passage 10")
mo.data10.factor$REP <- as.factor(mo.data10.factor$REP)

#### Mosquito Passage 10 ####
## Create list of all .vcf files for Passage Number 10 
setwd("~/Documents/Research/1-Friedrich/1-ZIKV-SNVs/Kasen-Tutorial/Datasets/MP_refaligned_vcf")
mq.myvcf10 <- dir(pattern="MP10_")    
mq.n10 <- length(mq.myvcf10)             
mq.mylist10 <- vector("list",mq.n10)     

## Read all the tables in myvcf and apply to mylist10
for (i in 1:mq.n10) {mq.mylist10[[i]] <- read.table(mq.myvcf10[i])} 

## Split INFO into 4 columns, respectively
for (i in 1:mq.n10) {mq.mylist10[[i]] <- mq.mylist10[[i]] %>% 
  separate("V8",c("V8_1","V8_2"),sep=";",convert=FALSE) %>% 
  separate("V8_1",c("V8_1_1","V8_1_2"),sep="=",convert=FALSE) %>% 
  separate("V8_2",c("V8_2_1","V8_2_2"),sep="=",convert=FALSE)
}

## Remqve unnecessary columns and rename remaining  
for (i in 1:mq.n10) {mq.mylist10[[i]] <- mq.mylist10[[i]][,!(colnames(mq.mylist10[[i]]) %in% c("V1","V3","V6","V7","V8_1_1","V8_1_2","V8_2_1"))]}
for (i in 1:mq.n10) {mq.mylist10[[i]] <- dplyr::rename(mq.mylist10[[i]],CHROM=V2,POS=V2,REF=V4,ALT=V5,AF=V8_2_2)}

## Create "data" data.frame for all mylist10
for (i in 1:mq.n10) {mq.mylist10[[i]]$AF <- as.numeric(mq.mylist10[[i]]$AF)}
for (i in 1:mq.n10) {mq.mylist10[[i]]$AF <- mq.mylist10[[i]]$AF*100} #As % frequency
for (i in 1:mq.n10) {mq.mylist10[[i]]$SUB <- paste(mq.mylist10[[i]]$POS,mq.mylist10[[i]]$REF,">",mq.mylist10[[i]]$ALT)}
for (i in 1:mq.n10) {mq.mylist10[[i]]$NUC <- paste(mq.mylist10[[i]]$POS,mq.mylist10[[i]]$ALT)}
for (i in 1:mq.n10) {mq.mylist10[[i]]$SUB <- gsub(" ", "",mq.mylist10[[i]]$SUB,fixed=TRUE)}
for (i in 1:mq.n10) {mq.mylist10[[i]] <- mq.mylist10[[i]][,!(colnames(mq.mylist10[[i]]) %in% c("REF","column_label"))]}

## Filter out all barcodes with AF <= .3%
for (i in 1:mq.n10) {mq.mylist10[[i]] <- filter(mq.mylist10[[i]],AF>.3)}

## Give each list its own data frame
mq.data10.1 <- mq.mylist10[[1]] %>% 
  add_column(PN=10) %>% 
  add_column(REP=1)
mq.data10.2 <- mq.mylist10[[2]] %>% 
  add_column(PN=10) %>% 
  add_column(REP=2)
mq.data10.3 <- mq.mylist10[[3]] %>% 
  add_column(PN=10) %>% 
  add_column(REP=3)
mq.data10.4 <- mq.mylist10[[4]] %>% 
  add_column(PN=10) %>% 
  add_column(REP=4)
mq.data10.5 <- mq.mylist10[[5]] %>% 
  add_column(PN=10) %>% 
  add_column(REP=5)
mq.data10 <- rbind(mq.data10.1,mq.data10.2,mq.data10.3,mq.data10.4,mq.data10.5)

mq.data10.factor <- mq.data10
mq.data10.factor$REP <- as.factor(mq.data10$REP)
levels(mq.data10.factor$REP) <- c("Lineage A", "Lineage B", "Lineage C", "Lineage D", "Lineage E")
mq.data10.factor$PN <- as.factor(mq.data10$PN)
levels(mq.data10.factor$PN) <- c("Passage 10")
mq.data10.factor$REP <- as.factor(mq.data10.factor$REP)

#### Alternate Passage 10 ####
## Create list of all .vcf files for Passage Number 10 
setwd("~/Documents/Research/1-Friedrich/1-ZIKV-SNVs/Kasen-Tutorial/Datasets/SCP_refaligned_vcf")
alt.myvcf10 <- dir(pattern="SCP10_")    
alt.n10 <- length(alt.myvcf10)             
alt.mylist10 <- vector("list",alt.n10)     

## Read all the tables in myvcf and apply to mylist10
for (i in 1:alt.n10) {alt.mylist10[[i]] <- read.table(alt.myvcf10[i])} 

## Split INFO into 4 columns, respectively
for (i in 1:alt.n10) {alt.mylist10[[i]] <- alt.mylist10[[i]] %>% 
  separate("V8",c("V8_1","V8_2"),sep=";",convert=FALSE) %>% 
  separate("V8_1",c("V8_1_1","V8_1_2"),sep="=",convert=FALSE) %>% 
  separate("V8_2",c("V8_2_1","V8_2_2"),sep="=",convert=FALSE)
}

## Realtve unnecessary columns and rename remaining  
for (i in 1:alt.n10) {alt.mylist10[[i]] <- alt.mylist10[[i]][,!(colnames(alt.mylist10[[i]]) %in% c("V1","V3","V6","V7","V8_1_1","V8_1_2","V8_2_1"))]}
for (i in 1:alt.n10) {alt.mylist10[[i]] <- dplyr::rename(alt.mylist10[[i]],CHROM=V2,POS=V2,REF=V4,ALT=V5,AF=V8_2_2)}

## Create "data" data.frame for all mylist10
for (i in 1:alt.n10) {alt.mylist10[[i]]$AF <- as.numeric(alt.mylist10[[i]]$AF)}
for (i in 1:alt.n10) {alt.mylist10[[i]]$AF <- alt.mylist10[[i]]$AF*100} #As % frequency
for (i in 1:alt.n10) {alt.mylist10[[i]]$SUB <- paste(alt.mylist10[[i]]$POS,alt.mylist10[[i]]$REF,">",alt.mylist10[[i]]$ALT)}
for (i in 1:alt.n10) {alt.mylist10[[i]]$NUC <- paste(alt.mylist10[[i]]$POS,alt.mylist10[[i]]$ALT)}
for (i in 1:alt.n10) {alt.mylist10[[i]]$SUB <- gsub(" ", "",alt.mylist10[[i]]$SUB,fixed=TRUE)}
for (i in 1:alt.n10) {alt.mylist10[[i]] <- alt.mylist10[[i]][,!(colnames(alt.mylist10[[i]]) %in% c("REF","column_label"))]}

## Filter out all barcodes with AF <= .3%
for (i in 1:alt.n10) {alt.mylist10[[i]] <- filter(alt.mylist10[[i]],AF>.3)}

## Give each list its own data frame
alt.data10.1 <- alt.mylist10[[1]] %>% 
  add_column(PN=10) %>% 
  add_column(REP=1)
alt.data10.2 <- alt.mylist10[[2]] %>% 
  add_column(PN=10) %>% 
  add_column(REP=2)
alt.data10.3 <- alt.mylist10[[3]] %>% 
  add_column(PN=10) %>% 
  add_column(REP=3)
alt.data10.4 <- alt.mylist10[[4]] %>% 
  add_column(PN=10) %>% 
  add_column(REP=4)
alt.data10.5 <- alt.mylist10[[5]] %>% 
  add_column(PN=10) %>% 
  add_column(REP=5)
alt.data10 <- rbind(alt.data10.1,alt.data10.2,alt.data10.3,alt.data10.4,alt.data10.5)

alt.data10.factor <- alt.data10
alt.data10.factor$REP <- as.factor(alt.data10$REP)
levels(alt.data10.factor$REP) <- c("Lineage A", "Lineage B", "Lineage C", "Lineage D", "Lineage E")
alt.data10.factor$PN <- as.factor(alt.data10$PN)
levels(alt.data10.factor$PN) <- c("Passage 10")
alt.data10.factor$REP <- as.factor(alt.data10.factor$REP)

#### Plots ####
TEN.mo.facet <- ggplot(mo.data10.factor,aes(POS,AF)) + 
  geom_point(size=.75,color="steelblue3") + 
  scale_y_continuous(expand=c(0,0),limits=(c(0,100))) + 
  scale_x_continuous(expand=c(0,0),limits=(c(0,10400)),breaks=c(106,473,977,2489,3545,4223,4613,6464,6845,6914,7667,10378,10675)) + 
  labs(y="SNV Frequency (%)",title="Mice - Passage 10") + 
  #5'UTR	1	106
  #C	107	472
  #prM	473	976
  #E	977	2488
  #NS1	2489	3544
  #NS2A	3545	4222
  #NS2B	4223	4612
  #NS3	4613	6463
  #NS4A	6464	6844
  #2K	6845	6913
  #NS4B	6914	7666
#NS5	7667	10378
#3'UTR	10379	10675
facet_grid(cols=vars(REP)) + 
  theme(panel.background=element_rect(fill="transparent",color=NA),panel.spacing.x = unit(2.5,"lines"), 
        axis.line.x=element_line(color="black"),axis.line.y=element_line(color="black"),
        axis.ticks=element_line(color="black"),
        axis.title.x=element_text(color="white"),axis.title.y=element_text(color="black",size=10),
        axis.text.x=element_text(color="white",angle = 90, hjust = 1,size=8),
        strip.background=element_rect(color="white",fill="white"),
        axis.text.y=element_text(color="black"),panel.spacing = unit(1.5,"lines"))

TEN.mq.facet <- ggplot(mq.data10.factor,aes(POS,AF)) + 
  geom_point(size=.75,color="steelblue3") + 
  scale_y_continuous(expand=c(0,0),limits=(c(0,100))) + 
  scale_x_continuous(expand=c(0,0),limits=(c(0,10400)),breaks=c(106,473,977,2489,3545,4223,4613,6464,6845,6914,7667,10378,10675)) + 
  labs(y="SNV Frequency (%)",title="Mosquito - Passage 10") + 
  facet_grid(cols=vars(REP)) + 
  theme(panel.background=element_rect(fill="transparent",color=NA),panel.spacing.x = unit(2.5,"lines"), 
        axis.line.x=element_line(color="black"),axis.line.y=element_line(color="black"),
        axis.ticks=element_line(color="black"),
        axis.title.x=element_text(color="white"),axis.title.y=element_text(color="black",size=10),
        axis.text.x=element_text(color="white",angle = 90, hjust = 1,size=8),
        strip.background=element_rect(color="white",fill="white"),
        axis.text.y=element_text(color="black"),panel.spacing = unit(1.5,"lines"))

TEN.alt.facet <- ggplot(alt.data10.factor,aes(POS,AF)) + 
  geom_point(size=.75,color="steelblue3") + 
  scale_y_continuous(expand=c(0,0),limits=(c(0,100))) + 
  scale_x_continuous(expand=c(0,0),limits=(c(0,10400)),breaks=c(106,473,977,2489,3545,4223,4613,6464,6845,6914,7667,10378,10675)) + 
  labs(y="SNV Frequency (%)",title="Placeholder - Passage 10") + 
  facet_grid(cols=vars(REP)) + 
  theme(panel.background=element_rect(fill="transparent",color=NA),panel.spacing.x = unit(2.5,"lines"), 
        axis.line.x=element_line(color="black"),axis.line.y=element_line(color="black"),
        axis.ticks=element_line(color="black"),
        axis.title.x=element_text(color="white"),axis.title.y=element_text(color="black",size=10),
        axis.text.x=element_text(color="white",angle = 90, hjust = 1,size=8),
        strip.background=element_rect(color="white",fill="white"),
        axis.text.y=element_text(color="black"),panel.spacing = unit(1.5,"lines"))

#### Figure Layout ####
grid.arrange(TEN.mo.facet, TEN.mq.facet, TEN.alt.facet, nrow = 3) # 800x700

#### Mouse Passage CC ####
## Create list of all .vcf files for Passage Number CC 
setwd("~/Documents/Research/1-Friedrich/1-ZIKV-SNVs/Kasen-Tutorial/Datasets/SCP_refaligned_vcf")
mo.myvcfCC <- dir(pattern="SCPCC_")    
mo.nCC <- length(mo.myvcfCC)             
mo.mylistCC <- vector("list",mo.nCC)     

## Read all the tables in myvcf and apply to mylistCC
for (i in 1:mo.nCC) {mo.mylistCC[[i]] <- read.table(mo.myvcfCC[i])} 

## Split INFO into 4 columns, respectively
for (i in 1:mo.nCC) {mo.mylistCC[[i]] <- mo.mylistCC[[i]] %>% 
  separate("V8",c("V8_1","V8_2"),sep=";",convert=FALSE) %>% 
  separate("V8_1",c("V8_1_1","V8_1_2"),sep="=",convert=FALSE) %>% 
  separate("V8_2",c("V8_2_1","V8_2_2"),sep="=",convert=FALSE)
}

## Remove unnecessary columns and rename remaining  
for (i in 1:mo.nCC) {mo.mylistCC[[i]] <- mo.mylistCC[[i]][,!(colnames(mo.mylistCC[[i]]) %in% c("V1","V3","V6","V7","V8_1_1","V8_1_2","V8_2_1"))]}
for (i in 1:mo.nCC) {mo.mylistCC[[i]] <- dplyr::rename(mo.mylistCC[[i]],CHROM=V2,POS=V2,REF=V4,ALT=V5,AF=V8_2_2)}

## Create "data" data.frame for all mylistCC
for (i in 1:mo.nCC) {mo.mylistCC[[i]]$AF <- as.numeric(mo.mylistCC[[i]]$AF)}
for (i in 1:mo.nCC) {mo.mylistCC[[i]]$AF <- mo.mylistCC[[i]]$AF*100} #As % frequency
for (i in 1:mo.nCC) {mo.mylistCC[[i]]$SUB <- paste(mo.mylistCC[[i]]$POS,mo.mylistCC[[i]]$REF,">",mo.mylistCC[[i]]$ALT)}
for (i in 1:mo.nCC) {mo.mylistCC[[i]]$NUC <- paste(mo.mylistCC[[i]]$POS,mo.mylistCC[[i]]$ALT)}
for (i in 1:mo.nCC) {mo.mylistCC[[i]]$SUB <- gsub(" ", "",mo.mylistCC[[i]]$SUB,fixed=TRUE)}
for (i in 1:mo.nCC) {mo.mylistCC[[i]] <- mo.mylistCC[[i]][,!(colnames(mo.mylistCC[[i]]) %in% c("REF","column_label"))]}

## Filter out all barcodes with AF <= .3%
for (i in 1:mo.nCC) {mo.mylistCC[[i]] <- filter(mo.mylistCC[[i]],AF>.3)}

## Give each list its own data frame
mo.dataCC.1 <- mo.mylistCC[[1]] %>% 
  add_column(PN="CC") %>% 
  add_column(REP=1)
mo.dataCC.2 <- mo.mylistCC[[2]] %>% 
  add_column(PN="CC") %>% 
  add_column(REP=2)
mo.dataCC.3 <- mo.mylistCC[[3]] %>% 
  add_column(PN="CC") %>% 
  add_column(REP=3)
mo.dataCC.4 <- mo.mylistCC[[4]] %>% 
  add_column(PN="CC") %>% 
  add_column(REP=4)
mo.dataCC.5 <- mo.mylistCC[[5]] %>% 
  add_column(PN="CC") %>% 
  add_column(REP=5)
mo.dataCC <- rbind(mo.dataCC.1,mo.dataCC.2,mo.dataCC.3,mo.dataCC.4,mo.dataCC.5)

mo.dataCC.factor <- mo.dataCC
mo.dataCC.factor$REP <- as.factor(mo.dataCC$REP)
levels(mo.dataCC.factor$REP) <- c("Lineage A", "Lineage B", "Lineage C", "Lineage D", "Lineage E")
mo.dataCC.factor$PN <- as.factor(mo.dataCC$PN)
levels(mo.dataCC.factor$PN) <- c("Passage CC")
mo.dataCC.factor$REP <- as.factor(mo.dataCC.factor$REP)

#### Mosquito Passage CC ####
## Create list of all .vcf files for Passage Number CC 
setwd("~/Documents/Research/1-Friedrich/1-ZIKV-SNVs/Kasen-Tutorial/Datasets/MP_refaligned_vcf")
mq.myvcfCC <- dir(pattern="MPCC_")    
mq.nCC <- length(mq.myvcfCC)             
mq.mylistCC <- vector("list",mq.nCC)     

## Read all the tables in myvcf and apply to mylistCC
for (i in 1:mq.nCC) {mq.mylistCC[[i]] <- read.table(mq.myvcfCC[i])} 

## Split INFO into 4 columns, respectively
for (i in 1:mq.nCC) {mq.mylistCC[[i]] <- mq.mylistCC[[i]] %>% 
  separate("V8",c("V8_1","V8_2"),sep=";",convert=FALSE) %>% 
  separate("V8_1",c("V8_1_1","V8_1_2"),sep="=",convert=FALSE) %>% 
  separate("V8_2",c("V8_2_1","V8_2_2"),sep="=",convert=FALSE)
}

## Remqve unnecessary columns and rename remaining  
for (i in 1:mq.nCC) {mq.mylistCC[[i]] <- mq.mylistCC[[i]][,!(colnames(mq.mylistCC[[i]]) %in% c("V1","V3","V6","V7","V8_1_1","V8_1_2","V8_2_1"))]}
for (i in 1:mq.nCC) {mq.mylistCC[[i]] <- dplyr::rename(mq.mylistCC[[i]],CHROM=V2,POS=V2,REF=V4,ALT=V5,AF=V8_2_2)}

## Create "data" data.frame for all mylistCC
for (i in 1:mq.nCC) {mq.mylistCC[[i]]$AF <- as.numeric(mq.mylistCC[[i]]$AF)}
for (i in 1:mq.nCC) {mq.mylistCC[[i]]$AF <- mq.mylistCC[[i]]$AF*100} #As % frequency
for (i in 1:mq.nCC) {mq.mylistCC[[i]]$SUB <- paste(mq.mylistCC[[i]]$POS,mq.mylistCC[[i]]$REF,">",mq.mylistCC[[i]]$ALT)}
for (i in 1:mq.nCC) {mq.mylistCC[[i]]$NUC <- paste(mq.mylistCC[[i]]$POS,mq.mylistCC[[i]]$ALT)}
for (i in 1:mq.nCC) {mq.mylistCC[[i]]$SUB <- gsub(" ", "",mq.mylistCC[[i]]$SUB,fixed=TRUE)}
for (i in 1:mq.nCC) {mq.mylistCC[[i]] <- mq.mylistCC[[i]][,!(colnames(mq.mylistCC[[i]]) %in% c("REF","column_label"))]}

## Filter out all barcodes with AF <= .3%
for (i in 1:mq.nCC) {mq.mylistCC[[i]] <- filter(mq.mylistCC[[i]],AF>.3)}

## Give each list its own data frame
mq.dataCC.1 <- mq.mylistCC[[1]] %>% 
  add_column(PN="CC") %>% 
  add_column(REP=1)
mq.dataCC.2 <- mq.mylistCC[[2]] %>% 
  add_column(PN="CC") %>% 
  add_column(REP=2)
mq.dataCC.3 <- mq.mylistCC[[3]] %>%                #PLACEHOLDER VALUES!!!
  add_column(PN="CC") %>% 
  add_column(REP=3)
mq.dataCC.4 <- mq.mylistCC[[3]] %>% 
  add_column(PN="CC") %>% 
  add_column(REP=4)
mq.dataCC.5 <- mq.mylistCC[[4]] %>%
  add_column(PN="CC") %>% 
  add_column(REP=5)
mq.dataCC <- rbind(mq.dataCC.1,mq.dataCC.2,mq.dataCC.3,mq.dataCC.4,mq.dataCC.5)

mq.dataCC.factor <- mq.dataCC
mq.dataCC.factor$REP <- as.factor(mq.dataCC$REP)
levels(mq.dataCC.factor$REP) <- c("Lineage A", "Lineage B", "Placeholder", "Lineage D", "Lineage E")
mq.dataCC.factor$PN <- as.factor(mq.dataCC$PN)
levels(mq.dataCC.factor$PN) <- c("Passage CC")
mq.dataCC.factor$REP <- as.factor(mq.dataCC.factor$REP)

#### Alternate Passage CC ####
## Create list of all .vcf files for Passage Number CC 
setwd("~/Documents/Research/1-Friedrich/1-ZIKV-SNVs/Kasen-Tutorial/Datasets/SCP_refaligned_vcf")
alt.myvcfCC <- dir(pattern="SCPCC_")    
alt.nCC <- length(alt.myvcfCC)             
alt.mylistCC <- vector("list",alt.nCC)     

## Read all the tables in myvcf and apply to mylistCC
for (i in 1:alt.nCC) {alt.mylistCC[[i]] <- read.table(alt.myvcfCC[i])} 

## Split INFO into 4 columns, respectively
for (i in 1:alt.nCC) {alt.mylistCC[[i]] <- alt.mylistCC[[i]] %>% 
  separate("V8",c("V8_1","V8_2"),sep=";",convert=FALSE) %>% 
  separate("V8_1",c("V8_1_1","V8_1_2"),sep="=",convert=FALSE) %>% 
  separate("V8_2",c("V8_2_1","V8_2_2"),sep="=",convert=FALSE)
}

## Realtve unnecessary columns and rename remaining  
for (i in 1:alt.nCC) {alt.mylistCC[[i]] <- alt.mylistCC[[i]][,!(colnames(alt.mylistCC[[i]]) %in% c("V1","V3","V6","V7","V8_1_1","V8_1_2","V8_2_1"))]}
for (i in 1:alt.nCC) {alt.mylistCC[[i]] <- dplyr::rename(alt.mylistCC[[i]],CHROM=V2,POS=V2,REF=V4,ALT=V5,AF=V8_2_2)}

## Create "data" data.frame for all mylistCC
for (i in 1:alt.nCC) {alt.mylistCC[[i]]$AF <- as.numeric(alt.mylistCC[[i]]$AF)}
for (i in 1:alt.nCC) {alt.mylistCC[[i]]$AF <- alt.mylistCC[[i]]$AF*100} #As % frequency
for (i in 1:alt.nCC) {alt.mylistCC[[i]]$SUB <- paste(alt.mylistCC[[i]]$POS,alt.mylistCC[[i]]$REF,">",alt.mylistCC[[i]]$ALT)}
for (i in 1:alt.nCC) {alt.mylistCC[[i]]$NUC <- paste(alt.mylistCC[[i]]$POS,alt.mylistCC[[i]]$ALT)}
for (i in 1:alt.nCC) {alt.mylistCC[[i]]$SUB <- gsub(" ", "",alt.mylistCC[[i]]$SUB,fixed=TRUE)}
for (i in 1:alt.nCC) {alt.mylistCC[[i]] <- alt.mylistCC[[i]][,!(colnames(alt.mylistCC[[i]]) %in% c("REF","column_label"))]}

## Filter out all barcodes with AF <= .3%
for (i in 1:alt.nCC) {alt.mylistCC[[i]] <- filter(alt.mylistCC[[i]],AF>.3)}

## Give each list its own data frame
alt.dataCC.1 <- alt.mylistCC[[1]] %>% 
  add_column(PN="CC") %>% 
  add_column(REP=1)
alt.dataCC.2 <- alt.mylistCC[[2]] %>% 
  add_column(PN="CC") %>% 
  add_column(REP=2)
alt.dataCC.3 <- alt.mylistCC[[3]] %>% 
  add_column(PN="CC") %>% 
  add_column(REP=3)
alt.dataCC.4 <- alt.mylistCC[[4]] %>% 
  add_column(PN="CC") %>% 
  add_column(REP=4)
alt.dataCC.5 <- alt.mylistCC[[5]] %>% 
  add_column(PN="CC") %>% 
  add_column(REP=5)
alt.dataCC <- rbind(alt.dataCC.1,alt.dataCC.2,alt.dataCC.3,alt.dataCC.4,alt.dataCC.5)

alt.dataCC.factor <- alt.dataCC
alt.dataCC.factor$REP <- as.factor(alt.dataCC$REP)
levels(alt.dataCC.factor$REP) <- c("Lineage A", "Lineage B", "Lineage C", "Lineage D", "Lineage E")
alt.dataCC.factor$PN <- as.factor(alt.dataCC$PN)
levels(alt.dataCC.factor$PN) <- c("Passage CC")
alt.dataCC.factor$REP <- as.factor(alt.dataCC.factor$REP)

#### Plots ####
CC.mo.facet <- ggplot(mo.dataCC.factor,aes(POS,AF)) + 
  geom_point(size=.75,color="steelblue3") + 
  scale_y_continuous(expand=c(0,0),limits=(c(0,100))) + 
  scale_x_continuous(expand=c(0,0),limits=(c(0,10400)),breaks=c(106,473,977,2489,3545,4223,4613,6464,6845,6914,7667,10378,10675)) + 
  labs(y="SNV Frequency (%)",title="Mice - Passage CC") + 
  facet_grid(cols=vars(REP)) + 
  theme(panel.background=element_rect(fill="transparent",color=NA),panel.spacing.x = unit(2.5,"lines"), 
        axis.line.x=element_line(color="black"),axis.line.y=element_line(color="black"),
        axis.ticks=element_line(color="black"),
        axis.title.x=element_text(color="white"),axis.title.y=element_text(color="black",size=10),
        axis.text.x=element_text(color="white",angle = 90, hjust = 1,size=8),
        strip.background=element_rect(color="white",fill="white"),
        axis.text.y=element_text(color="black"),panel.spacing = unit(1.5,"lines"))

CC.mq.facet <- ggplot(mq.dataCC.factor,aes(POS,AF)) + 
  geom_point(size=.75,color="steelblue3") + 
  scale_y_continuous(expand=c(0,0),limits=(c(0,100))) + 
  scale_x_continuous(expand=c(0,0),limits=(c(0,10400)),breaks=c(106,473,977,2489,3545,4223,4613,6464,6845,6914,7667,10378,10675)) + 
  labs(y="SNV Frequency (%)",title="Mosquito - Passage CC") + 
  facet_grid(cols=vars(REP)) + 
  theme(panel.background=element_rect(fill="transparent",color=NA),panel.spacing.x = unit(2.5,"lines"), 
        axis.line.x=element_line(color="black"),axis.line.y=element_line(color="black"),
        axis.ticks=element_line(color="black"),
        axis.title.x=element_text(color="white"),axis.title.y=element_text(color="black",size=10),
        axis.text.x=element_text(color="white",angle = 90, hjust = 1,size=8),
        strip.background=element_rect(color="white",fill="white"),
        axis.text.y=element_text(color="black"),panel.spacing = unit(1.5,"lines"))

CC.alt.facet <- ggplot(alt.dataCC.factor,aes(POS,AF)) + 
  geom_point(size=.75,color="steelblue3") + 
  scale_y_continuous(expand=c(0,0),limits=(c(0,100))) + 
  scale_x_continuous(expand=c(0,0),limits=(c(0,10400)),breaks=c(106,473,977,2489,3545,4223,4613,6464,6845,6914,7667,10378,10675)) + 
  labs(y="SNV Frequency (%)",title="Placeholder - Passage CC") + 
  facet_grid(cols=vars(REP)) + 
  theme(panel.background=element_rect(fill="transparent",color=NA),panel.spacing.x = unit(2.5,"lines"), 
        axis.line.x=element_line(color="black"),axis.line.y=element_line(color="black"),
        axis.ticks=element_line(color="black"),
        axis.title.x=element_text(color="white"),axis.title.y=element_text(color="black",size=10),
        axis.text.x=element_text(color="white",angle = 90, hjust = 1,size=8),
        strip.background=element_rect(color="white",fill="white"),
        axis.text.y=element_text(color="black"),panel.spacing = unit(1.5,"lines"))

#### Figure Layout ####
grid.arrange(CC.mo.facet, CC.mq.facet, CC.alt.facet, nrow = 3)

# Save figures as 700x700 in RStudio