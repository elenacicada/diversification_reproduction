#Workflow for BiSSE + HiSSE only 

#set working directory to location of download 
setwd("~/scripts") 

library(hisse)
library(phytools)
library(TreeTools)
library(dplyr)

#remove all objects except those that take a long time to create
rm(list=ls()[! ls() %in% c("res.mdh.all","mono_dio_herm")]) 

#calling in "outersect" function for later - see: https://tonybreyal.wordpress.com/2011/11/29/outersect-the-opposite-of-rs-intersect-function/
#this is basically just the opposite of the 'intersect' function 
outersect <- function(x, y, ...) {
  big.vec <- c(x, y, ...)
  duplicates <- big.vec[duplicated(big.vec)]
  setdiff(big.vec, unique(duplicates))
}

#reading in large multi-family tree - Smith and Brown (2018), available: 
#https://github.com/FePhyFoFum/big_seed_plant_trees/releases
big.tree <- read.newick("GBOTB.tre")

#getting a character vector of tip labels
tips<- TipLabels(big.tree) 
tips.df <- as.data.frame(tips)
tips.df <- rename(.data = tips.df, Name = "tips")

#trait data - from Meyer et al.2025
ms.dat.fin <- read.csv("ms.data.ch2.csv") 
msdat.rename <- rename(ms.dat.fin, Name = "genus_species")

#filter for specific family -- change this manually to family of interest
fam.name <- "Solanaceae"
size.fam.ms <- filter(msdat.rename, Family == fam.name)

#overlap tips of tree with mating system data
overlap.fam.bigtree.ms <- merge(tips.df, size.fam.ms, by = "Name") 

#set fam size for family for later + filter 
#numbers from Christenhusz & Byng (2016):
#https://doi.org/10.11646/phytotaxa.261.3.1)
fam.number <- 2600 

#overlapping pt. 1 - tree + mating system
overlap.fam.bigtree.ms <- merge(tips.df, size.fam.ms, by = "Name") 

#tabulate 
table(overlap.fam.bigtree.ms$mating_system)
real.ms <- table(overlap.fam.bigtree.ms$mating_system)

#re-coding traits 
overlap.fam.bigtree.ms$mating_system[overlap.fam.bigtree.ms$mating_system == "sc"] <- 1
overlap.fam.bigtree.ms$mating_system[overlap.fam.bigtree.ms$mating_system == "si"] <- 0

#formatting for tree 
A1 <- tips.df$Name #tip labels from tree
A2 <- overlap.fam.bigtree.ms$Name #species names for trait data

x <- setdiff(A1, A2) #get number of non-overlaps

intersect.phy <- intersect(A2, A1) 
outersect.phy <- outersect(A2, A1) 

tip <- outersect.phy
tip2 <- gsub(" ", "_", tip)
tree.dropped <- drop.tip(big.tree, tip2) 
tree.dropped.u <- force.ultrametric(tree.dropped)

#formatting for data 
match.df <- data.frame(overlap.fam.bigtree.ms, row.names = 1)
final.match <- as.integer(match.df$mating_system)
names(final.match) <- rownames(match.df)

#getting sampling props (using family size assigned earlier)
n = length(final.match) #getting sample size 
ms.tot = table(size.fam.ms$mating_system) #generating table of traits
sc.tot = ms.tot[1] #assigning from table 
si.tot = ms.tot[2]

sc.prop <- sc.tot/(sc.tot + si.tot)
si.prop <- si.tot/(sc.tot + si.tot)

n.theory.sc <- sc.prop * fam.number
n.theory.si <- si.prop * fam.number

sc.real = real.ms[1]
si.real = real.ms[2]

rho.si <- si.real/n.theory.si
rho.sc <- sc.real/n.theory.sc

rho <- c(rho.si, rho.sc)

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

#create HiSSE design matrix for BiSSE-like model
rates.bisse <- TransMatMakerHiSSE(hidden.traits = 0)

#matrix for HiSSE 
rates.hisse <- TransMatMakerHiSSE(hidden.traits = 1) 

#fit bisse with HiSSE
p <- starting.point.bisse(tree.dropped.u) 
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
rates.cid2

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

##CID-4 with fixed extinction
rates.cid4fe <- TransMatMakerHiSSE(hidden.traits = 3)
rates.cid4fe[!is.na(rates.cid4fe)] <- 1 

cid.4fe <- hisse(tree.dropped.u, df.match.rename, f=rho, 
                   turnover = c(1,1,2,2,3,3,4,4), eps = c(1,1,1,1,1,1,1,1),
                   hidden.states = TRUE,
                   trans.rate = rates.cid4fe)

#Full HiSSE
hisse.full <- hisse(tree.dropped.u, df.match.rename, f=rho, 
                       turnover = c(1,2,3,4,5,6,7,8), eps = c(1,2,3,4,5,6,7,8),
                       hidden.states = TRUE,
                       trans.rate = rates.cid4)

#Full HiSSE w/extinction set to be uniform 

hisse.setE <- hisse(tree.dropped.u, df.match.rename, f=rho, 
                       turnover = c(1,2,3,4,5,6,7,8), eps = c(1,1,1,1,1,1,1,1),
                       hidden.states = TRUE,
                       trans.rate = rates.cid4)
