#Code for BiSSE (implemented via HiSSE package) w/sexual system traits 
#This is for sexual system data for which there is only two states 

#set working directory to location of download 
setwd("~/scripts") 

#load relevant packages
library(devtools)
library(dplyr)
library(hisse)
library(phytools)
library(readxl)
library(TreeTools)
library(diversitree)

#install name-checking tool if needed
#devtools::install_github("ecoinfor/U.Taxonstand")

#calling in "outersect" function for later - see: https://tonybreyal.wordpress.com/2011/11/29/outersect-the-opposite-of-rs-intersect-function/
#this is basically just the opposite of the 'intersect' function 
outersect <- function(x, y, ...) {
  big.vec <- c(x, y, ...)
  duplicates <- big.vec[duplicated(big.vec)]
  setdiff(big.vec, unique(duplicates))
}

#calling in data for sexual system (mdh) + name-checking 
#from Wang et al. (2021): https://datadryad.org/dataset/doi:10.5061/dryad.s1rn8pk7n
mono_dio_herm <- read_excel("mono.dio.herm.xlsx")

#reading in reference for correcting names from World Flora Online 
dat1 <- read_excel("Plants_WFO_database_part1.xlsx")
dat2 <- read_excel("Plants_WFO_database_part2.xlsx")
dat3 <- read_excel("Plants_WFO_database_part3.xlsx")

#rbinding to single datatable
database <- rbind(dat1, dat2, dat3)

#run name check
res.mdh.all <- nameMatch(spList=mono_dio_herm, spSource=database, author=FALSE, max.distance=1, genusPairs=NULL, Append=FALSE)

#reading in large multi-family tree - Smith and Brown (2018), available: 
#https://github.com/FePhyFoFum/big_seed_plant_trees/releases
big.tree <- read.newick("GBOTB.tre")

#creates a character vector of tip labs
tips<- TipLabels(big.tree) 
tips.df <- as.data.frame(tips)
tips.df <- rename(.data = tips.df, Name = "tips")

#set object with family of choice - filtering for family of choice, set here
fam.name <- "Lamiaceae"
size.fam.data <- filter(res.mdh.all, Family == fam.name)

##set fam size for family for later + filter 
#numbers from Christenhusz & Byng (2016):
#https://doi.org/10.11646/phytotaxa.261.3.1)
fam.number <- 7530

name.trait <- rename(.data = size.fam.data, Name = Accepted_SPNAME)
size.fam.trait <- merge(mono_dio_herm, name.trait, by = "Name")

#frequency of trait across family (based on whole dataset)
all.fam.freq <- table(size.fam.trait$SexSyst)

#getting distinct values for family 
dist.fam <- distinct(.data=size.fam.data, Accepted_SPNAME, .keep_all = TRUE) 

#fixing formatting with underscore and making dataframe 
string.fam <- dist.fam$Accepted_SPNAME
string.fam.under <- gsub(" ", "_", string.fam)
df.fam.under <- as.data.frame(string.fam.under)  

#renaming columns to match as "Name"
fam.mdh <- rename(df.fam.under, Name = "string.fam.under")

#overlapping pt. 1 - tree + mdh 
overlap.fam.bigtree <- merge(tips.df, fam.mdh, by = "Name") 

#did underscore not in R for some reason (?) -- we will match this for traits
mdh.underscore <- read.csv("mdh2.csv")

#overlap 2 - matching for traits 
match <- merge(mdh.underscore, overlap.fam.bigtree, by = "Name")
table(match$SexSyst)
sex.syst <- table(match$SexSyst)

#recoding traits 
match$SexSyst[match$SexSyst == "dioecy"] <- 0
match$SexSyst[match$SexSyst == "hermaphroditism"] <- 1

#formating for tree 
A1 <- tips.df$Name #tip labels from tree
A2 <- match$Name #species names for trait data

#get number of non-overlaps
x <- setdiff(A1, A2) 

intersect.phy <- intersect(A2, A1) 
outersect.phy <- outersect(A2, A1) 

tip <- outersect.phy
tip2 <- gsub(" ", "_", tip)
tree.dropped <- drop.tip(big.tree, tip2) 
tree.dropped.u <- force.ultrametric(tree.dropped)

#formatting for data 
match.df <- data.frame(match, row.names = 1)
final.match <- as.integer(match.df$SexSyst)
names(final.match) <- rownames(match.df)

#sample fraction calculations 
#must use fractions from whole mdh dataset for this family for "real" calcs 
#object we need is "all.fam.freqs" 

ms.tot = table(match.df$SexSyst) #generating table of traits
dio.tot = ms.tot[1] #assigning from table - dio 
mono.tot = ms.tot[2] #assigning from table - herm

dio.prop <- all.fam.freq[1]/(all.fam.freq[1] + all.fam.freq[2]) #prop real 
mono.prop <- all.fam.freq[2]/(all.fam.freq[1] + all.fam.freq[2]) #prop real

n.theory.dio <- dio.prop * fam.number #prop theoretical 
n.theory.mono <- mono.prop * fam.number #prop theoretical 

rho.d <- dio.tot/n.theory.dio
rho.m <- mono.tot/n.theory.mono

rho <- c(rho.d, rho.m)

#HiSSE + BiSSE implimentation in HiSSE package

#creating function to back-transform - from Revell and Harmon (2022) book 
repar.bd <- function(object, k=2){ 
  pars<-object$solution
  tt<-pars[grep("turnover", names(pars))][1:k]
  ee <-pars[grep("eps", names(pars))][1:k]
  lambda <- tt/(1+ee) 
  mu <- tt-lambda
  nn <-sapply(strsplit(names(tt), "turnover"), 
              function(x) x[2])
  matrix(c(lambda,mu), k,2, dimnames=list(nn,
                                          c("lambda", "mu")))
}

#formatting matched data 
df.match <- as.data.frame(final.match)
d <- cbind(rownames(df.match), data.frame(df.match, row.names=NULL))
names(d)[names(d) == 'rownames(df.match)'] <- 'Genus.species'
df.match.rename <- rename(d, x = final.match)

#create HiSSE design matrix for BiSSE 
rates.bisse <- TransMatMakerHiSSE(hidden.traits = 0)

#matrix for HiSSE 
rates.hisse <- TransMatMakerHiSSE(hidden.traits = 1) 

#fit BiSSE with HiSSE
bisse <- hisse(tree.dropped.u, df.match.rename, turnover = c(1,2), 
                    eps=c(1,2), hidden.states = FALSE, 
                    trans.rate = rates.bisse, f=rho)

#CID model w/HiSSE 
cid <- hisse(tree.dropped.u, df.match.rename, turnover = c(1,1), 
                  eps=c(1,1), hidden.states = FALSE, 
                  trans.rate = rates.bisse, f=rho)

#CID-2 model 
rates.cid2 <- rates.hisse
rates.cid2[!is.na(rates.cid2)] <- 1 

cid.2 <- hisse(tree.dropped.u, df.match.rename, f=rho, 
                   turnover = c(1,1,2,2), eps = c(1,1,2,2),
                   hidden.states = TRUE,
                   trans.rate = rates.cid2)

#CID-4
rates.cid4 <- TransMatMakerHiSSE(hidden.traits = 3)
rates.cid4[!is.na(rates.cid4)] <- 1 

cid.4 <- hisse(tree.dropped.u, df.match.rename, f=rho, 
                   turnover = c(1,1,2,2,3,3,4,4), eps = c(1,1,2,2,3,3,4,4),
                   hidden.states = TRUE,
                   trans.rate = rates.cid4)

#Full HiSSE

hisse.full <- hisse(tree.dropped.u, df.match.rename, f=rho, 
                       turnover = c(1,2,3,4,5,6,7,8), eps = c(1,2,3,4,5,6,7,8),
                       hidden.states = TRUE,
                       trans.rate = rates.cid4)
