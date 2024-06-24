####STRUCTURE run first to identify potential pure and hybrid individuals in the dataset

####

####INTROGRESS analysis to confirm hybrid ID of individuals, estimate hybrid index,
####construct triangle plots and investigate clines

## load the introgress library 
library(introgress)

## load helper libraries
library(adegenet)  #readstructure command to genind 
library(dplyr)  # data manipulation functions
# If you don't already have devtools installed
#install.packages("devtools")
#devtools::install_github("btmartin721/ClineHelpR")
library("ClineHelpR") #genind2introgress function to create input object

## set the working directory in R to an existing directory where you will work
## and where the downloaded data files have been saved
setwd("C:/Users/ataylor66/OneDrive - University of Central Oklahoma/Desktop/R_projects/Impound_study_R/Impound_study_R")

###importing structure file with all data (P1 & P2 refs + admixed wild pop)
###column 3 has pop names (e.g., refSHB, refALB, wildpop); column 2 has site numbers saved as "other (use @other in genind)
x <- read.structure("Lanier_inputs_for_R_Introgress2.stru", n.ind = 204, n.loc = 16, onerowperind = TRUE,
                    col.lab = 1, col.pop = 3, col.others = 2, NA.char = "-9", ask = TRUE, quiet = TRUE)

###examine genind object for accuracy and completeness
is.genind(x)
summary(x)
###view as a dataframe to check import 
xdf <- genind2df(x, oneColPerAll = FALSE, sep="/")


## read in marker information (this is not needed b/c the above function also provides this info)
loci.data<- read.csv(file="locidata_lanier.txt", header=TRUE)

###"If pop.id = TRUE and ind.id = TRUE the first row of admix.gen should give the population identification (i.e. sampling locality) 
###of each individual and the second row should provide a unique individual identification; genotype information would then begin on row three.
xdf <- cbind(ind = rownames(xdf), xdf)
rownames(xdf) <- 1:nrow(xdf)
xdf %>% relocate(pop, .before = ind)

####P1 and P2 don't need the ind and pop info, just genotypes.
####Then the tables must be transposed before passed thru prepare.data 
p1.data <- filter(xdf, pop=="refCHB")
p1.data <- select(p1.data, -c(ind, pop))
p1.data <- t(p1.data)
p2.data <- filter(xdf, pop=="refALB")
p2.data <- select(p2.data, -c(ind, pop))
p2.data <- t(p2.data)
admix.data <- filter(xdf, pop=="wildpop")
admix.data <- t(admix.data)


introgress_input <- prepare.data(admix.gen=admix.data, loci.data=loci.data, parental1=p1.data, parental2=p2.data, pop.id=TRUE, ind.id=TRUE, fixed=FALSE, sep.rows=FALSE, sep.columns=FALSE)


## estimate hybrid index values and save the results
##see https://www.rdocumentation.org/packages/introgress/versions/1.2.3/topics/est.h
## max likelihood est of the proportion of genome inherited from pop 2 (h is pt est) and 95% CI's
hi.index <- est.h(introgress.data=introgress_input, loci.data=loci.data)
hi.index
write.csv(hi.index, "CHBxALB_hi_index.csv")

## estimate inter-specific heterozygosity
##This function calculates an admixed individual's interspecific heterozygosity 
##(i.e. the proportion of the individual's genome with alleles inherited from both parental populations) based on allele counts. 
##This function should only be used for co-dominant markers

in.het<-calc.intersp.het(introgress.data=introgress_input)
in.het
write.csv(in.het, "CHBxALB_in_het.csv")

##triangle plot
##This function plots interspecific heterozygosity as a function 
###of hybrid index for individuals from an admixed population. 
###Individuals that are the progeny of at least one parent from one of the pure parental 
###populations should have maximal heterozygosity for the observed hybrid index. 
###The plot has lines that correspond to these theoretical maximum values. 
###Individuals that fall on the maximal line are likely F1s or backcross progeny. 
###Evidence for individuals of this type will be more likely if the data set 
###consists of loci with no alleles in common between parental species (delta=1), 
###whereas shared alleles will lead to ambiguity in inferring ancestry. 
###Hybrid index estimates can be obtained from the est.h function and interspecific 
###heterozygosity estimates can be obtained from the calc.intersp.het function.

triPlot <- triangle.plot(hi.index=hi.index, int.het=in.het, pdf=TRUE, out.file="CHBxALB_tri.pdf")
triplot <- triangle.plot(hi.index=hi.index, int.het=in.het, pdf=FALSE)


## make plot to visualize patterns of introgression
## this saves the plot to a pdf in the current directory for R
## Xlab.h is "population 2 ancestry"
mk.image(introgress.data=introgress_input, loci.data=loci.data, 
         hi.index=hi.index, ylab.image="Individuals",
         xlab.h="ALB ancestry", pdf=TRUE, out.file="CHBxALB_image.pdf")

## conduct the genomic clines analysis and save the results to a
## data object (list) called gen.out
## this uses the permutation procedure with 1000 permutations
cline.out <- genomic.clines(introgress.data=introgress_input,
                            hi.index=hi.index, loci.data=loci.data,
                            sig.test=TRUE, method="parametric", het.cor=TRUE)

## make plots to visualize the genomic clines
clineplot <- clines.plot(cline.data=cline.out, rplots=3, cplots=3, pdf=FALSE)
clineplot
## this saves the plots to a pdf in the current directory for R
clines.plot(cline.data=cline.out, rplots=3, cplots=3, pdf=TRUE, out.file="CHBxALB_clines.pdf")




