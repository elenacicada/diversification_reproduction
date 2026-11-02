#Workflow for MuSSE + multi-state HiSSE (MuHiSSE)

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

#creating function to back-transform - from Liam Revell's book
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

#calling in "outersect" function for later 
outersect <- function(x, y, ...) {
  big.vec <- c(x, y, ...)
  duplicates <- big.vec[duplicated(big.vec)]
  setdiff(big.vec, unique(duplicates))
}

#calling in data for sexual system (mdh) + name-checking 
#originally from Wang et al. (2021) - link:
#https://onlinelibrary.wiley.com/doi/full/10.1111/ele.13815
mono_dio_herm <- read_excel("mono.dio.herm.xlsx")

#reading in reference for correcting names from World Flora Online 
#from https://www.sciencedirect.com/science/article/pii/S2468265922000944 
dat1 <- read_excel("Plants_WFO_database_part1.xlsx")
dat2 <- read_excel("Plants_WFO_database_part2.xlsx")
dat3 <- read_excel("Plants_WFO_database_part3.xlsx")

#rbinding to single datatable
database <- rbind(dat1, dat2, dat3)

#name checking procedure 
res.mdh.all <- nameMatch(spList=mono_dio_herm, 
                         spSource=database, 
                         author=FALSE, 
                         max.distance=1, 
                         genusPairs=NULL, 
                         Append=FALSE)

#reading in large multi-family tree + extracting tip labels
big.tree <- read.newick("GBOTB.tre")
tips<- TipLabels(big.tree) #creates a character vector of tip labs
tips.df <- as.data.frame(tips)
tips.df <- rename(.data = tips.df, Name = "tips")

#filtering for family of choice, set here
fam.name <- "Poaceae"
size.fam.data <- filter(res.mdh.all, Family == fam.name)

#add number of species in chosen family number of species
#we used Christenhusz & Byng (2016) to define family sizes:
#https://doi.org/10.11646/phytotaxa.261.3.1)
fam.number <- 12000 #set manually 
name.trait <- rename(.data = size.fam.data, Name = Accepted_SPNAME)
size.fam.trait <- merge(mono_dio_herm, name.trait, by = "Name")

#make table + check 
all.fam.freq <- table(size.fam.trait$SexSyst)
all.fam.freq

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

#what we have now does not have traits attached - re-match here 
traits.og <- mono_dio_herm
traits.og $Name <- gsub(" ", "_", traits.og$Name)
mdh.underscore <- merge(traits.og, overlap.fam.bigtree, by = "Name")

#overlap 2 - matching for traits 
match <- merge(mdh.underscore, overlap.fam.bigtree, by = "Name")
table(match$SexSyst)
sex.syst <- table(match$SexSyst)

#recoding traits 
match$SexSyst[match$SexSyst == "dioecy"] <- 1
match$SexSyst[match$SexSyst == "monoecy"] <- 2
match$SexSyst[match$SexSyst == "hermaphroditism"] <- 3

#formating for tree 
A1 <- tips.df$Name #tip labels from tree
A2 <- match$Name #species names for trait data

x <- setdiff(A1, A2) #number of non-overlaps

intersect.phy <- intersect(A2, A1) 
outersect.phy <- outersect(A2, A1) 

tip <- outersect.phy
tip2 <- gsub(" ", "_", tip)
tree.dropped <- drop.tip(big.tree, tip2) 
tree.dropped.u <- force.ultrametric(tree.dropped)

#formatting for data 
match.df <- data.frame(match, row.names = 1)
final.match <- as.double(match.df$SexSyst)
names(final.match) <- rownames(match.df)

#tree should be done
#need generalized sampling fraction calc 

#number of dio, mono, and herm for our present matched set + total n 
n = length(final.match)
dio.tot.fam = sex.syst[1]
mono.tot.fam = sex.syst[3]
herm.tot.fam = sex.syst[2]
sum.tot.data = (dio.tot.fam + mono.tot.fam + herm.tot.fam)

#number of dio, mono, and herm for the broader datset (mdh) + total n 
dio.tot.real <- all.fam.freq[1] 
mono.tot.real <- all.fam.freq [3]
herm.tot.real <- all.fam.freq [2] 
sum.tot.real = (herm.tot.real + mono.tot.real + dio.tot.real)

#calculation proportions of our data represented by each trait 
dio.prop.data <- dio.tot.fam/(sum.tot.data)
mono.prop.data <- mono.tot.fam/(sum.tot.data)
herm.prop.data <- herm.tot.fam/(sum.tot.data)

#calculation proportions from TOTAL data set (i.e., our best existing estimates)
dio.prop.real <- dio.tot.real/(sum.tot.real)
mono.prop.real <- mono.tot.real/(sum.tot.real)
herm.prop.real <- herm.tot.real/(sum.tot.real)

#theoretical proportions
#i.e., what our dataset would contain if we had full coverage 
n.theory.dio <- dio.prop.real * fam.number
n.theory.mono <- mono.prop.real * fam.number
n.theory.herm <- herm.prop.real * fam.number

#checking total for theoretical number
fam.tot.check <- n.theory.dio + n.theory.herm + n.theory.mono

rho.dio <- dio.tot.fam/n.theory.dio
rho.mono <- mono.tot.fam/n.theory.mono
rho.herm <- herm.tot.fam/n.theory.herm

tot.fam.check <- dio.tot.fam + mono.tot.fam + herm.tot.fam

#setting sampling fraction object below 
rho <- setNames(c(rho.dio, rho.mono, rho.herm), 1:3)

#initial setup for MuHiSSE - based on: 
#https://cran.r-project.org/web/packages/hisse/vignettes/MuHiSSE-vignette.pdf

#setting a 4th rate w/zero for MuHiSSE
rho2 <- setNames(c(rho.dio, rho.mono, rho.herm, 0), 1:4) 

#re-namimg objects to match sample code 
phy <- tree.dropped.u
tip.state <- final.match

#making all transition rate matrices 
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=0)
print(trans.rate)

#dropping 4th state ("11")
trans.rate.mod <- ParDrop(trans.rate, c(4,6,7,8))
print(trans.rate.mod)

#making state formatted correctly - mostly directly from tutorial 

states <- data.frame(tip.state, tip.state,
                     row.names=names(tip.state))

states.trans <- states
for(i in 1:Ntip(phy)){
  if(states[i,1] == 1){
    states.trans[i,1] = 0
    states.trans[i,2] = 0
  }
  if(states[i,1] == 2){
    states.trans[i,1] = 0
    states.trans[i,2] = 1
  }
  if(states[i,1] == 3){
    states.trans[i,1] = 1
    states.trans[i,2] = 0
  }
}

#changing rownames to their own column 
states.trans2 <- cbind(rownames(states.trans), 
                       data.frame(states.trans, row.names=NULL))

f = rho2

#running CID model (a minimal, a true "dull" null)
turnover <- c(1,1,1,0)
extinction.fraction <- c(1,1,1,0)

cid.null <- MuHiSSE(phy=phy, data=states.trans2, f=f, turnover=turnover,
                    eps=extinction.fraction, hidden.states=FALSE,
                    trans.rate=trans.rate.mod)

#running MuSSE equiv 
turnover <- c(1,2,3,0)
extinction.fraction <- c(1,1,1,0)

musse.equv <- MuHiSSE(phy=phy, data=states.trans2, f=f, turnover=turnover,
                     eps=extinction.fraction, hidden.states=FALSE,
                     trans.rate=trans.rate.mod)

#full MuHiSSE model w/hidden state but no extinction variation 

turnover <- c(1,2,3,0, 4,5,6,0)
extinction.fraction <- c(1,1,1,0, 1,1,1,0)
f = rho2
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=1)
trans.rate.mod <- ParDrop(trans.rate, c(4,6,7,8,12,14,15,16))

full.multi <- MuHiSSE(phy=phy, data=states.trans2, f=f, turnover=turnover,
                      eps=extinction.fraction, hidden.states=TRUE,
                      trans.rate=trans.rate.mod)

##full MuHiSSE model w/hidden state AND extinction variation 

turnover <- c(1,2,3,0, 4,5,6,0)
extinction.fraction <- c(1,2,3,0, 4,5,6,0)
f = rho2
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=1)
trans.rate.mod <- ParDrop(trans.rate, c(4,6,7,8,12,14,15,16))

full.multi.vary <- MuHiSSE(phy=phy, data=states.trans2, f=f, turnover=turnover,
                           eps=extinction.fraction, hidden.states=TRUE,
                           trans.rate=trans.rate.mod)

##Full MuHiSSE w/4 hidden states 

trans.rate.4.hidden <- TransMatMakerMuHiSSE(hidden.traits = 3)
View(trans.rate.4.hidden)

#full MuHiSSE model w/hidden state but no extinction variation 

turnover <- c(1,2,3,0, 4,5,6,0)
extinction.fraction <- c(1,1,1,0, 1,1,1,0)
f = rho2
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=1)
trans.rate.mod <- ParDrop(trans.rate, c(4,6,7,8,12,14,15,16))

full.multi <- MuHiSSE(phy=phy, data=states.trans2, f=f, turnover=turnover,
                      eps=extinction.fraction, hidden.states=TRUE,
                      trans.rate=trans.rate.mod)

#now doing character-independent models

#CID-2, which would be MuSSE + hidden state as a null 

#making cid2 model and dropping 4th state 
trans.rate.cid2 <- TransMatMakerMuHiSSE(hidden.traits=1)
trans.rate.mod.cid2 <- ParDrop(trans.rate.cid2, c(4,6,7,8,
                                                  12,14,15,16))

#setting all transition rates to be equal for cid2
trans.rate.mod.cid2[!is.na(trans.rate.mod.cid2) & !trans.rate.mod.cid2 == 0 ] = 1 

turnover = c(1,1,1,0, 2,2,2,0)
extinction.fraction = c(1,1,1,0, 1,1,1,0)

cid2 <- MuHiSSE(phy=phy, data=states.trans2, f=f, turnover=turnover,
                eps=extinction.fraction, hidden.states=TRUE,
                trans.rate=trans.rate.mod.cid2)

#CID-4 equiv for null to hidden state

#making cid4
trans.rates.cid4 <- TransMatMakerMuHiSSE(hidden.traits=3)
trans.rate.mod.cid4 <- ParDrop(trans.rates.cid4, c(4,6,7,8,
                                                   12,14,15,16, 
                                                   20,22,23,24,
                                                   28,30,31,32))

#setting all transition rates to be equal for cid4
trans.rate.mod.cid4[!is.na(trans.rate.mod.cid4) & !trans.rate.mod.cid4 == 0 ] = 1 


turnover = c(1,1,1,0, 2,2,2,0, 3,3,3,0, 4,4,4,0)
extinction.fraction = c(1,1,1,0, 1,1,1,0, 1,1,1,0, 1,1,1,0)

cid4 <- MuHiSSE(phy=phy, data=states.trans2, f=f, turnover=turnover,
                eps=extinction.fraction, hidden.states=TRUE,
                trans.rate = trans.rate.mod.cid4)
