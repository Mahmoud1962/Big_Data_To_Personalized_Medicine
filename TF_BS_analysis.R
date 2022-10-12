if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")

install.packages("bedr")

## [1] "/Users/irishelenejonkers/Downloads"
# first load the genes that you wish to do enrichment for. In this case, we use the top 250 differentially expressed genes that are upregulated after stimulation. 
DE_IFNg <- read.table("C:/Users/jespe/OneDrive/Documenten/RijksUniversiteit Groningen/Master BMS/From Big data to personalized medicine/Practical/Analysis/Final data analysis/IFNg_vs_Control_Sign.txt", header=T)
# make the column called Row.names the actual rownames and remove the now duplicate column
rownames(DE_IFNg) <- DE_IFNg$Row.names
DE_IFNg <- DE_IFNg[,-1]

# add additional information from biomaRt
library(biomaRt)

# define which biomart data you would like to use. Make sure you use the reference genome compatible with your initial analysis! (in our case GRCh37 aka hg19)
ensembl  <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
# select the attributes you want
genemap <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol","gene_biotype","chromosome_name","start_position","end_position","strand"),
                 filters = "ensembl_gene_id",
                 values = rownames(DE_IFNg),
                 mart = ensembl )

# Now merge the the DE data frame with the newly acquired ensembl data that is compiled in genemap. This can be done based on column names, where '0' represents row names.

DE_IFNg_extended <- merge(x = DE_IFNg, y = genemap, by.x = 0, by.y = "ensembl_gene_id", all.x = TRUE, all.y = FALSE)
head(DE_IFNg_extended)

# get promoters of each differentially expressed gene 
# first, create an empty matrix called "prom""
prom_IFNg=matrix(ncol=3, nrow=nrow(DE_IFNg_extended))
# by strand, add 350 bp upstream and 150 bp downstream to TSS of each gene of interest in a forward loop. Please note that the genes on the antisense strand will have the TSS at the end_position!
for (i in 1:nrow(DE_IFNg_extended)){
  if (DE_IFNg_extended$strand[i] == "1") 
  {prom_IFNg[i,1] <- as.character(DE_IFNg_extended$Row.names[i])
  prom_IFNg[i,2] <- DE_IFNg_extended$start_position[i]-350
  prom_IFNg[i,3] <- DE_IFNg_extended$start_position[i]+150
  }
  else {
    prom_IFNg[i,1] <- as.character(DE_IFNg_extended$Row.names[i])
    prom_IFNg[i,2] <- DE_IFNg_extended$end_position[i]-150
    prom_IFNg[i,3] <- DE_IFNg_extended$end_position[i]+350}   
}

colnames(prom_IFNg)=c("Row.names","prom_IFNg_start","prom_IFNg_end")

# make bedtools compatible bed format (chr, prom_start, prom_end, strand, gene.id, gene.name)
comb <- merge(DE_IFNg_extended, prom_IFNg, by.x="Row.names", by.y="Row.names")
DE_IFNg_ALL <- comb[,c("chromosome_name","prom_IFNg_start", "prom_IFNg_end","Row.names","hgnc_symbol","strand")]


# convert the file into something that the bedtools package can deal with
# So add new column names:
colnames(DE_IFNg_ALL) <- c("chr","start", "end", "ENSG", "hgnc", "strand")
# bedtools only recognizes chromosome names with chr in front of it. 
DE_IFNg_ALL$chr <- paste("chr", DE_IFNg_ALL$chr, sep="")
# And finally, the columns need to have certain characteristics for them to be read correctly
DE_IFNg_ALL$chr <- as.character(as.factor(DE_IFNg_ALL$chr))
DE_IFNg_ALL$start <- as.numeric(as.character(DE_IFNg_ALL$start))
DE_IFNg_ALL$end <- as.numeric(as.character(DE_IFNg_ALL$end))

# first load the genes that you want to use as controls. These are a random set of non-differentially expressed genes
CTR_IFNg <- read.table("C:/Users/jespe/OneDrive/Documenten/RijksUniversiteit Groningen/Master BMS/From Big data to personalized medicine/Practical/Analysis/Final data analysis/IFNg_vs_Control_Control_file1.txt", header=T)
# make the column called Row.names the actual rownames and remove the now duplicate column
rownames(CTR_IFNg) <- CTR_IFNg$Row.names
CTR_IFNg <- CTR_IFNg[,-c(1)]

# add additional information from biomaRt
library(biomaRt)

# define which biomart data you would like to use. Make sure you use the reference genome compatible with your initial analysis! (in our case GRCh37 aka hg19)
ensembl  <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
# select the attributes you want
genemap <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol","gene_biotype","chromosome_name","start_position","end_position","strand"),
                 filters = "ensembl_gene_id",
                 values = rownames(CTR_IFNg),
                 mart = ensembl )

# Now merge the the CTR_IFNg data frame with the newly acquired ensembl data that is compiled in genemap. This can be done based on column names, where '0' represents row names.

CTR_IFNg_extended <- merge(x = CTR_IFNg, y = genemap, by.x = 0, by.y = "ensembl_gene_id", all.x = TRUE, all.y = FALSE)
head(CTR_IFNg_extended)

prom_IFNg_Ctrl=matrix(ncol=3, nrow=nrow(CTR_IFNg_extended))
# by strand, add 350 bp upstream and 150 bp downstream to TSS of each gene of interest in a forward loop. Please note that the genes on the antisense strand will have the TSS at the end_position!
for (i in 1:nrow(CTR_IFNg_extended)){
  if (CTR_IFNg_extended$strand[i] == "1") 
  {prom_IFNg_Ctrl[i,1] <- as.character(CTR_IFNg_extended$Row.names[i])
  prom_IFNg_Ctrl[i,2] <- CTR_IFNg_extended$start_position[i]-350
  prom_IFNg_Ctrl[i,3] <- CTR_IFNg_extended$start_position[i]+150
  }
  else {
    prom_IFNg_Ctrl[i,1] <- as.character(CTR_IFNg_extended$Row.names[i])
    prom_IFNg_Ctrl[i,2] <- CTR_IFNg_extended$end_position[i]-150
    prom_IFNg_Ctrl[i,3] <- CTR_IFNg_extended$end_position[i]+350}  
}

colnames(prom_IFNg_Ctrl)=c("Row.names","prom_IFNg_Ctrl_start","prom_IFNg_Ctrl_end")

# make bedtools compatible bed format (chr, prom_IFNg_Ctrl_start, prom_IFNg_Ctrl_end, strand, gene.id, gene.name)
comb_IFNg <- merge(CTR_IFNg_extended, prom_IFNg_Ctrl, by.x="Row.names", by.y="Row.names")
CTR_IFNg150 <- comb_IFNg[,c("chromosome_name","prom_IFNg_Ctrl_start", "prom_IFNg_Ctrl_end","Row.names","hgnc_symbol","strand")]


# convert the file into something that the bedtools package can deal with
# So add new column names:
colnames(CTR_IFNg150) <- c("chr","start", "end", "ENSG", "hgnc", "strand")
# bedtools only recognizes chromosome names with chr in front of it. 
CTR_IFNg150$chr <- paste("chr", CTR_IFNg150$chr, sep="")
# And finally, the columns need to have certain characteristics for them to be read correctly
CTR_IFNg150$chr <- as.character(as.factor(CTR_IFNg150$chr))
CTR_IFNg150$start <- as.numeric(as.character(CTR_IFNg150$start))
CTR_IFNg150$end <- as.numeric(as.character(CTR_IFNg150$end))  


library(bedr)


DE_IFNg_Sign = subset(DE_IFNg_ALL, select = -c(ENSG,hgnc,strand))
CTR_IFNg150_Sign = subset(CTR_IFNg150, select = -c(ENSG,hgnc,strand))
sort_IFNg_DE150 <- bedr.sort.region(DE_IFNg_Sign)
ctr_IFNg_DE150 <- bedr.sort.region(CTR_IFNg150_Sign)


write.table(sort_IFNg_DE150, file = "IFNg_sort_Sign.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

write.table(ctr_IFNg_DE150, file = "IFNg_sort_Ctrl_150.txt", sep = "\t",
            row.names = TRUE, col.names = NA)





## [1] "/Users/irishelenejonkers/Downloads"
# first load the genes that you wish to do enrichment for. In this case, we use the top 250 differentially expressed genes that are upregulated after stimulation. 
DE_IL17 <- read.table("C:/Users/jespe/OneDrive/Documenten/RijksUniversiteit Groningen/Master BMS/From Big data to personalized medicine/Practical/Analysis/Final data analysis/IL17_vs_Control_Sign.txt", header=T)
# make the column called Row.names the actual rownames and remove the now duplicate column
rownames(DE_IL17) <- DE_IL17$Row.names
DE_IL17 <- DE_IL17[,-1]

# add additional information from biomaRt
library(biomaRt)

# define which biomart data you would like to use. Make sure you use the reference genome compatible with your initial analysis! (in our case GRCh37 aka hg19)
ensembl  <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
# select the attributes you want
genemap <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol","gene_biotype","chromosome_name","start_position","end_position","strand"),
                 filters = "ensembl_gene_id",
                 values = rownames(DE_IL17),
                 mart = ensembl )

# Now merge the the DE data frame with the newly acquired ensembl data that is compiled in genemap. This can be done based on column names, where '0' represents row names.

DE_IL17_extended <- merge(x = DE_IL17, y = genemap, by.x = 0, by.y = "ensembl_gene_id", all.x = TRUE, all.y = FALSE)
head(DE_IL17_extended)

# get prom_IL17oters of each differentially expressed gene 
# first, create an empty matrix called "prom_IL17""
prom_IL17=matrix(ncol=3, nrow=nrow(DE_IL17_extended))
# by strand, add 350 bp upstream and 150 bp downstream to TSS of each gene of interest in a forward loop. Please note that the genes on the antisense strand will have the TSS at the end_position!
for (i in 1:nrow(DE_IL17_extended)){
  if (DE_IL17_extended$strand[i] == "1") 
  {prom_IL17[i,1] <- as.character(DE_IL17_extended$Row.names[i])
  prom_IL17[i,2] <- DE_IL17_extended$start_position[i]-350
  prom_IL17[i,3] <- DE_IL17_extended$start_position[i]+150
  }
  else {
    prom_IL17[i,1] <- as.character(DE_IL17_extended$Row.names[i])
    prom_IL17[i,2] <- DE_IL17_extended$end_position[i]-150
    prom_IL17[i,3] <- DE_IL17_extended$end_position[i]+350}   
}

colnames(prom_IL17)=c("Row.names","prom_IL17_start","prom_IL17_end")

# make bedtools compatible bed format (chr, prom_IL17_start, prom_IL17_end, strand, gene.id, gene.name)
comb_IL17 <- merge(DE_IL17_extended, prom_IL17, by.x="Row.names", by.y="Row.names")
DE_IL17_ALL <- comb_IL17[,c("chromosome_name","prom_IL17_start", "prom_IL17_end","Row.names","hgnc_symbol","strand")]


# convert the file into something that the bedtools package can deal with
# So add new column names:
colnames(DE_IL17_ALL) <- c("chr","start", "end", "ENSG", "hgnc", "strand")
# bedtools only recognizes chromosome names with chr in front of it. 
DE_IL17_ALL$chr <- paste("chr", DE_IL17_ALL$chr, sep="")
# And finally, the columns need to have certain characteristics for them to be read correctly
DE_IL17_ALL$chr <- as.character(as.factor(DE_IL17_ALL$chr))
DE_IL17_ALL$start <- as.numeric(as.character(DE_IL17_ALL$start))
DE_IL17_ALL$end <- as.numeric(as.character(DE_IL17_ALL$end))

# first load the genes that you want to use as controls. These are a random set of non-differentially expressed genes
CTR_IL17 <- read.table("C:/Users/jespe/OneDrive/Documenten/RijksUniversiteit Groningen/Master BMS/From Big data to personalized medicine/Practical/Analysis/Final data analysis/IL17_vs_Control_Control_file1.txt", header=T)
# make the column called Row.names the actual rownames and remove the now duplicate column
rownames(CTR_IL17) <- CTR_IL17$Row.names
CTR_IL17 <- CTR_IL17[,-c(1)]

# add additional information from biomaRt
library(biomaRt)

# define which biomart data you would like to use. Make sure you use the reference genome compatible with your initial analysis! (in our case GRCh37 aka hg19)
ensembl  <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
# select the attributes you want
genemap <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol","gene_biotype","chromosome_name","start_position","end_position","strand"),
                 filters = "ensembl_gene_id",
                 values = rownames(CTR_IL17),
                 mart = ensembl )

# Now merge the the CTR_IL17 data frame with the newly acquired ensembl data that is compiled in genemap. This can be done based on column names, where '0' represents row names.

CTR_IL17_extended <- merge(x = CTR_IL17, y = genemap, by.x = 0, by.y = "ensembl_gene_id", all.x = TRUE, all.y = FALSE)
head(CTR_IL17_extended)

prom_IL17_ctrl=matrix(ncol=3, nrow=nrow(CTR_IL17_extended))
# by strand, add 350 bp upstream and 150 bp downstream to TSS of each gene of interest in a forward loop. Please note that the genes on the antisense strand will have the TSS at the end_position!
for (i in 1:nrow(CTR_IL17_extended)){
  if (CTR_IL17_extended$strand[i] == "1") 
  {prom_IL17_ctrl[i,1] <- as.character(CTR_IL17_extended$Row.names[i])
  prom_IL17_ctrl[i,2] <- CTR_IL17_extended$start_position[i]-350
  prom_IL17_ctrl[i,3] <- CTR_IL17_extended$start_position[i]+150
  }
  else {
    prom_IL17_ctrl[i,1] <- as.character(CTR_IL17_extended$Row.names[i])
    prom_IL17_ctrl[i,2] <- CTR_IL17_extended$end_position[i]-150
    prom_IL17_ctrl[i,3] <- CTR_IL17_extended$end_position[i]+350}  
}

colnames(prom_IL17_ctrl)=c("Row.names","prom_IL17_ctrl_start","prom_IL17_ctrl_end")

# make bedtools compatible bed format (chr, prom_IL17_ctrl_start, prom_IL17_ctrl_end, strand, gene.id, gene.name)
comb <- merge(CTR_IL17_extended, prom_IL17_ctrl, by.x="Row.names", by.y="Row.names")
CTR_IL17150 <- comb[,c("chromosome_name","prom_IL17_ctrl_start", "prom_IL17_ctrl_end","Row.names","hgnc_symbol","strand")]


# convert the file into something that the bedtools package can deal with
# So add new column names:
colnames(CTR_IL17150) <- c("chr","start", "end", "ENSG", "hgnc", "strand")
# bedtools only recognizes chromosome names with chr in front of it. 
CTR_IL17150$chr <- paste("chr", CTR_IL17150$chr, sep="")
# And finally, the columns need to have certain characteristics for them to be read correctly
CTR_IL17150$chr <- as.character(as.factor(CTR_IL17150$chr))
CTR_IL17150$start <- as.numeric(as.character(CTR_IL17150$start))
CTR_IL17150$end <- as.numeric(as.character(CTR_IL17150$end))  


library(bedr)


DE_IL17_Sign = subset(DE_IL17_ALL, select = -c(ENSG,hgnc,strand))
CTR_IL17150_Sign = subset(CTR_IL17150, select = -c(ENSG,hgnc,strand))
sort_IL17_DESIGN <- bedr.sort.region(DE_IL17_Sign)
ctr_IL17_DE40 <- bedr.sort.region(CTR_IL17150_Sign)



write.table(sort_IL17_DESIGN, file = "IL17_sort_Sign.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

write.table(ctr_IL17_DE40, file = "IL17_sort_Ctrl_40.txt", sep = "\t",
            row.names = TRUE, col.names = NA)
