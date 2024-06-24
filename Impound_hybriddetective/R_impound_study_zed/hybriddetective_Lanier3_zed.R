###Andrew Taylor, 5/17/2022 based on example codes - WORKS ON MACBOOK ONLY
###Trying to run with zeds for known parental types simulated 

###see:https://github.com/bwringe/hybriddetective/blob/master/README.Rmd
###See: https://github.com/bwringe/hybriddetective/blob/master/hybriddetective_example_RScript.pdf 
###Trying hybriddetective from their example code to ensure all works on local machine
###Installed NewHybrids, PLINK, and PGDspider to local machine first in "software" folder


library(genepopedit)
library(parallelnewhybrid)
library(genepopedit)
library(hybriddetective)


####PREPARE DATA FOR INPUT####

###set working directory, hold it to add folders to later...
setwd("/Users/andrewtaylor/Desktop/R_Projects/R_impound_study/R_impound_study_zed/Lanier3/")
path.hold <-getwd()


####SIMULATED HYBRIDS FROM REF POPS####
## Here I am simulating genotypes for 6 hybrid categories using the top loci discovered in the PLINK validation samples for each parent dataset. 
## The freqbasedsim_AlleleSample() function generates a number of individuals in each hybrid category based on allele frequencies in references

## Genotypes were simulated 3 separate times (S1, S2, S3), and for each set of simulations (unique genotype datasest), 
##I performed 3 replicate hybrid generation simulations (R1, R2, R3). 
##This way, I can assess deviations across simulations (S) AND within simulations (R).

##Make sure to read in the file created from getTopLoc function
##*outputName*| An optional character vector to be applied as the name of the output. 
####The default is NULL, in which case the output name is constructed from the name of the input, 
##with the *suffix _SiRj_NH* added where *i* is the number of simulations corresponding to the output, 
##and *j* is the number of replicates of the *ith* simulation. 
##NH refers to the fact that the output is in NewHybrids format 
##default sample size is 200 individuals per sim

freqbasedsim_GTFreq(GenePopData ="/Users/andrewtaylor/Desktop/R_Projects/R_impound_study/R_impound_study_zed/Lanier3/Lanier3_ref_for_hybriddetective_copy.txt", 
                        NumSims = 3, 
                        NumReps = 3,
                        sample.sizePure = 50,
                        sample.sizeF1 = 50,
                        sample.sizeF2 = 50,
                        sample.sizeBC = 50,
                        pop.groups = c("SHB", "REB"))

####NOTE:  Cut and past the simulated files out of the inputs and back into the Lanier3 folder!


#####ZEDS OPTION#####
###Optional step to add Zeds (known hybrid class designations) from .csv file to simulated individuals
###Can help NewHybrids more accurately model the expected genotypes of potential mixtures of two parentals
#z indicates that it is known beforehand the hybrid class of an individual
#0 indicates pure for pop 1, 1 indicates pure for pop 2
#s indicates the individual is to be used for calc of allele freqs only
#For more info, see NewHybrids User Guide

####Adding Zeds file here to match the reference pops created above
###Creating dataframe to match the number generated above of each pure pop (in this case, 50 each)
Individual <- c(1)
Zscore <- c("z0", "z0", "z0", "z0", "z0", "z0", "z0", "z0", "z0", "z0","z0", "z0", "z0", "z0", "z0","z0", "z0", "z0", "z0", "z0",
            "z0", "z0", "z0", "z0", "z0","z0", "z0", "z0", "z0", "z0","z0", "z0", "z0", "z0", "z0","z0", "z0", "z0", "z0", "z0",
            "z0", "z0", "z0", "z0", "z0","z0", "z0", "z0", "z0", "z0",
            "z1", "z1", "z1", "z1", "z1","z1", "z1", "z1", "z1", "z1","z1", "z1", "z1", "z1", "z1","z1", "z1", "z1", "z1", "z1",
            "z1", "z1", "z1", "z1", "z1","z1", "z1", "z1", "z1", "z1","z1", "z1", "z1", "z1", "z1","z1", "z1", "z1", "z1", "z1",
            "z1", "z1", "z1", "z1", "z1","z1", "z1", "z1", "z1", "z1")
Zeds_simulated <- data.frame(Individual, Zscore)
Zeds_simulated$Individual<-1:nrow(Zeds_simulated)

## Save the Zeds data to the working directory as a file called "SimZeds.txt"
write.table(x = Zeds_simulated, file = "/Users/andrewtaylor/Desktop/R_Projects/R_impound_study/R_impound_study_zed/Lanier3_inputs/SimZeds.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.csv(Zeds_simulated, file = "/Users/andrewtaylor/Desktop/R_Projects/R_impound_study/R_impound_study_zed/Lanier3_inputs/SimZeds.csv", row.names = FALSE, quote = FALSE)

###Applies the saved zed file across the simulations to be run by newhybrids
nh_Zcore(GetstheZdir = "/Users/andrewtaylor/Desktop/R_Projects/R_impound_study/R_impound_study_zed/Lanier3/",
multiapplyZvec = "/Users/andrewtaylor/Desktop/R_Projects/R_impound_study/R_impound_study_zed/Lanier3_inputs/SimZeds.csv")

####RUN SIMULATED HYBRIDS######
###Run Simulated hybrids with parallelnewhybrids to assign simulated genotypes to hybrid class 
###(Bayesian Gibbs sampling algorithm). 
###NewHybrids was installed on my local machine (MAC OS) to conduct the analysis
###Be sure to record final burnin and sweeps for methods

###NOTE: had to install newhybrids from Git, then XCode app was needed to compile and run on local machine successfully

parallelnh_OSX("/Users/andrewtaylor/Desktop/R_Projects/R_impound_study/R_impound_study_zed/Lanier3/", 
               where.NH = "/Users/andrewtaylor/newhybrids/", 
              burnin = 500000, 
              sweeps = 1000000)

## precheckR will check for NewHybrids results in all folders within the folder specified by "PreDir". i.e. all results in the folders within the "NH.Results" folder created by parallelnh_xx from the package parallelnewhybrid
##NOTE: before running this, remove from NH.Results folder 3 subfolders "..Panel.txt", "SimPurePops.txt", and "SimZeds.txt" _Results
nh_preCheckR(PreDir="/Users/andrewtaylor/Desktop/R_Projects/R_impound_study/R_impound_study_zed/Lanier3/NH.Results/", propCutOff = 0.5, PofZCutOff = 0.1)
## run with the default cut-off values; these could have been left blank.
### preCheckR will now flag any results where more than 50% (propCutOFF) of EITHER Pure1 or Pure2 individuals are given a probability of being an F2 greater than 10% (PofZCutOff)


#####EVALUATE MARKER POWER ON SIMULATED DATA######
##Plot the results so far..."structure plot" style
##order should be mostly pure 1, pure 2, F1, F2, BC1, BC2...otherwise it may show failure to converge
nh_multiplotR(NHResults="/Users/andrewtaylor/Desktop/R_Projects/R_impound_study/R_impound_study_zed/Lanier3/NH.Results/")

####Evaluate the hybrid-class specific efficacy of the panel across a range of
####posterior probability of assignment thresholds

###Note: if this command is used with only one panel size, it will issue a "geom_path" warning that can be ignored
###Creates 31 plots "Plot_1" thru Plot_31 (see ReadME), saves them all as tables, too with .csv files
hybridPowerComp(dir="/Users/andrewtaylor/Desktop/R_Projects/R_impound_study/R_impound_study_zed/Lanier3/NH.Results/", samplesize = c(50,50,50,50,50,50))


###NOTE: Getting weird error Error in dimnames(x) <- dn : length of 'dimnames' [2] not equal to array extent
####If I remove one of the three simulations (all three reps) it works, seems as if a coding error.  Maybe need to deconstruct this code to make it work.
####For now, I will just use 2 simulations, 3 replicates each, to calculate and avoid this error.
####It is unresolved bug in GitHub.  :(



###########COMBINE WITH UNKNOWN GENOTYPES#########
#####Combine the experimental (unknown) genotypes with the simulated pure genotypes
nh_analysis_generateR(
  ReferencePopsData = "/Users/andrewtaylor/Desktop/R_Projects/R_impound_study/R_impound_study_zed/Lanier3/Lanier3_ref_for_hybriddetective_copy_S2R2_NH.txt",
  UnknownIndivs = "/Users/andrewtaylor/Desktop/R_Projects/R_impound_study/R_impound_study_zed/Lanier3_inputs/Lanier3_wild_for_hybriddetective.txt",
  output.name = "NH_Lanier3_assignments_final_output.txt")

###Manually Take the resulting two files and save them in a subfolder of NH.Results called "Combined"

###Make the zscore file manually from the individuals file
###Applies the saved zed file across the simulations to be run by newhybrids
nh_Zcore(GetstheZdir = "/Users/andrewtaylor/Desktop/R_Projects/R_impound_study/R_impound_study_zed/Lanier3/Combined/",
         multiapplyZvec = "/Users/andrewtaylor/Desktop/R_Projects/R_impound_study/R_impound_study_zed/Lanier3_inputs/SimZeds.csv")


###Run NewHybrids on the combined data
parallelnh_OSX("/Users/andrewtaylor/Desktop/R_Projects/R_impound_study/R_impound_study_zed/Lanier3/Combined/", 
               where.NH = "/Users/andrewtaylor/newhybrids/", 
               burnin = 500000, 
               sweeps = 1000000)

#check for convergence
nh_preCheckR(PreDir="/Users/andrewtaylor/Desktop/R_Projects/R_impound_study/R_impound_study_zed/Lanier3/Combined/NH.Results/", propCutOff = 0.5, PofZCutOff = 0.1)

###Visualize results
###First shows pure 1 (blue) and pure 2 (red) from simulated pops, then unknown individuals
nh_plotR("/Users/andrewtaylor/Desktop/R_Projects/R_impound_study/R_impound_study_zed/Lanier3/Combined/NH.Results/NH_Lanier3_assignments_final_output.txt_Results/NH_Lanier3_assignments_final_output.txt_PofZ.txt")

###Plotting both ways so plot can be saved with colors that match original simulations for the system
nh_plotR("/Users/andrewtaylor/Desktop/R_Projects/R_impound_study/R_impound_study_zed/Lanier3/Combined/NH.Results/NH_Lanier3_assignments_final_output.txt_Results/NH_Lanier3_assignments_final_output.txt_PofZ.txt", ReversePure = 1)
nh_plotR("/Users/andrewtaylor/Desktop/R_Projects/R_impound_study/R_impound_study_zed/Lanier3/Combined/NH.Results/NH_Lanier3_assignments_final_output.txt_Results/NH_Lanier3_assignments_final_output.txt_PofZ.txt", ReversePure = 2)
#Manually save as desired, e.g., jpeg 300 dpi and pdf

