###Code to visualize outputs from Structure, NewHybrids, and Introgress in integrative ways
##Primarily using packages dplyr, ggplot2, patchwork, and fishualize palettes
#Written A. Taylor 8/30/23 to 5/21/24, referenced examples from various linked authors in code 


#####IMPORT DATA#####
##### Import the combined data to dataframe & clean
dat <- read.csv("/Users/andrewtaylor/Desktop/R_Projects/R_impound_IntrogressNPlots/InputData/Lanier_combined_2024_05_21.csv")

#Remove refs to keep only wild caught fish
dat <- subset(dat, SampType == 'wildpop')

##IF we want to remove the specimens that were had mixed ancestry, not included in NewHybrids 
dat2 <- subset(dat, No.Newhyb. != "N")

####Create new dataframes to hold the bounding lines for triangle plot
# Create a dataframe
line1 <- data.frame(x = seq(0, 0.5, by = 0.01), y = seq(0, 1, by = 0.02))
line2 <- data.frame(x = seq(0.5, 1, by = 0.01), y = seq(1, 0, by = -0.02))

#####Importing NewHybrids P of Z's for structure-like plots later
NHdat <- read.csv("/Users/andrewtaylor/Desktop/R_Projects/R_impound_IntrogressNPlots/InputData/Lanier_NHPofZ.csv", na.strings = "")

##Also removing the specimens that were had mixed ancestry, not included in NewHybrids 
NHdat <- subset(NHdat, Final_NH_assign != "MULTIPLE")

#####STRUCTURE PLOT DATA PREP#####

#With help from Arregoitia code here: https://luisdva.github.io/rstats/model-cluster-plots/
#With help from French code too: https://connor-french.github.io/intro-pop-structure-r/


#First need to mutate data so that individuals have multiple rows, each with STRU val
#Select columns to keep
struplotdat <- subset(dat, select = c(Samp, Pop, class_nh_50, FLB_LMB_stru, ALB_stru, REB_stru, CHB_stru, SHB_stru))

#Sort by NEO high to low within each population first (for plotting viz)
library(dplyr)
struplotdat <- dplyr::arrange(struplotdat, Pop, SHB_stru)

#see idea here for idea on how to order: https://youtu.be/yZSWFND3ca0
#Make a samp_order attribute as a column here
struplotdat$samp_ord <- row.names(struplotdat)
#checking structure of data
str(struplotdat)
#changing samp_ord from character to numeric
struplotdat$samp_ord <- as.numeric(struplotdat$samp_ord)


#reshaping data for plotting (each row now has one structure prob per cluster)
library(tidyr)
struplotdat <- tidyr::pivot_longer(struplotdat, cols = c(FLB_LMB_stru, ALB_stru, REB_stru, CHB_stru, SHB_stru), names_to = "cluster", values_to = "q")


#NH P of Z's for later use in plotting
NHplotdat <-dplyr::right_join(dat2,NHdat, by = "Samp")
NHplotdat <- dplyr::arrange(NHplotdat, Pop)
NHplotdat$samp_ord <- as.numeric(row.names(NHplotdat))

####NOTE: THIS PORTION IS DIFFERENT FROM TENKILLER BC OF MULTIPLE HEWHYBRIDS COMPARISONS
####Different PofZ df subset for each cross-comparison of parental taxa
#renaming each class so they progress logically to match color palettes
NHplotdat1 <- tidyr::pivot_longer(NHplotdat, cols = c(X1_P1_ALB, X1_P2_SHB, X1_F1, X1_F2, X1_BCALB, X1_BCSHB), names_to = "class", values_to = "PofZ")
NHplotdat1 <- subset(NHplotdat1, Newhyb.Run1..y != "NA")
NHplotdat1$class[NHplotdat1$class=="X1_P1_ALB"] <- "1.ALB"
NHplotdat1$class[NHplotdat1$class=="X1_BCALB"] <- "2.BCALB"
NHplotdat1$class[NHplotdat1$class=="X1_F1"] <- "3.F1"
NHplotdat1$class[NHplotdat1$class=="X1_F2"] <- "4.F2"
NHplotdat1$class[NHplotdat1$class=="X1_BCSHB"] <- "5.BCSHB"
NHplotdat1$class[NHplotdat1$class=="X1_P2_SHB"] <- "6.SHB"

NHplotdat2 <- tidyr::pivot_longer(NHplotdat, cols = c(X2_P1_ALB, X2_P2_CHB, X2_F1, X2_F2, X2_BCALB, X2_BCCHB), names_to = "class", values_to = "PofZ")
NHplotdat2 <- subset(NHplotdat2, Newhyb.Run2..y != "NA")
NHplotdat2$class[NHplotdat2$class=="X2_P1_ALB"] <- "1.ALB"
NHplotdat2$class[NHplotdat2$class=="X2_BCALB"] <- "2.BCALB"
NHplotdat2$class[NHplotdat2$class=="X2_F1"] <- "3.F1"
NHplotdat2$class[NHplotdat2$class=="X2_F2"] <- "4.F2"
NHplotdat2$class[NHplotdat2$class=="X2_BCCHB"] <- "5.BCCHB"
NHplotdat2$class[NHplotdat2$class=="X2_P2_CHB"] <- "6.CHB"

NHplotdat3 <- tidyr::pivot_longer(NHplotdat, cols = c(X3_P1_SHB, X3_P2_REB, X3_F1, X3_F2, X3_BCSHB, X3_BCREB), names_to = "class", values_to = "PofZ")
NHplotdat3 <- subset(NHplotdat3, Newhyb.Run3..y != "NA")
NHplotdat3$class[NHplotdat3$class=="X3_P1_SHB"] <- "1.SHB"
NHplotdat3$class[NHplotdat3$class=="X3_BCSHB"] <- "2.BCSHB"
NHplotdat3$class[NHplotdat3$class=="X3_F1"] <- "3.F1"
NHplotdat3$class[NHplotdat3$class=="X3_F2"] <- "4.F2"
NHplotdat3$class[NHplotdat3$class=="X3_BCREB"] <- "5.BCREB"
NHplotdat3$class[NHplotdat3$class=="X3_P2_REB"] <- "6.REB"

#####COLOR PALETTES & PLOTTING MADNESS####
library(ggplot2)
library(fishualize)


# Could use RColorBrewer or Fishualize, for example
# Create a divergent set for the clusters in Structure assignments
# in this case, ordered ALPHABETICALLY BY DEFAULT - ALB, CHB, LMB_FLB, REB, SHB  
# pulling from pal1 and a random complementary color for SPB
#library(RColorBrewer)
#q_palette <- brewer.pal(n=5, "BrBG")
#spectral q_palette <- c("#D7191C" "#FDAE61" "#FFFFBF" "#ABDDA4" "#2B83BA")

# Create a divergent color palette for hybrid classes
#library(fishualize)
#fish(n = 5, option = "Hypsypops_rubicundus")
pal_q <- c("#0C59FEFF", "#FEC700FF", "#FB7200FF", "#FC0F00FF", "#15E0FAFF")


###Palettes for NH comparison sets 1,2, & 3
#using this website and "lightgray" in hex #D3D3D3 to find intermediates for BC's
#https://meyerweb.com/eric/tools/color-blend/#D3D3D3:0C59FE:1:hex 
#pal1 = ALB, BCALB, F1, F2, BCSHB, SHB
pal_1 <-c("#0C59FEFF","#7096E9", "darkgray", "black", "#74DAE7", "#15E0FAFF")
pal_1.1 <-c("#0C59FEFF","#7096E9", "#119DFC", "black", "#74DAE7", "#15E0FAFF")
pal_1.2 <-c("#0C59FEFF","#7096E9", "darkgray", "black", "magenta", "#15E0FAFF")
#THEN, subset ordered as appears for triplots (SHB, BC SHB, ALB)
pal_1.2_tri <-c("#15E0FAFF", "magenta", "#0C59FEFF")

#pal2 = ALB, BCALB, F1, F2, BCCHB, CHB
pal_2 <-c("#0C59FEFF","#7096E9", "darkgray", "black", "#D4B855", "#FEC700FF")
pal_2.1<-c("#0C59FEFF","#7096E9", "#84667F", "black", "#D4B855", "#FEC700FF")
#THEN, subset ordered as appears for triplots ("CHB", "BC CHB x ALB","F1 ALB x CHB", "F2 ALB x CHB", "ALB")
pal_2_tri <- c("#FEC700FF","#D4B855", "darkgray", "black", "#0C59FEFF")

#pal3 = SHB, BCSHB, F1, F2, BCREB, REB
pal_3 <-c("#15E0FAFF","#74DAE7", "#9933CC", "black", "#E8716A", "#FC0F00FF")
pal_3.1 <-c("#15E0FAFF","#74DAE7", "#9933CC", "black", "#E8716A", "#FC0F00FF")
#THEN, subset ordered as appears for triplots ("F1 SHB x REB", "SHB")
pal_3_tri <-c("#15E0FAFF", "#9933CC")
###Need palette for all possible types in one
#run below to pull color options
#fish(n=8, option = "Hypsypops_rubicundus")
#Final plot has classes in this order: ALB, BC CHBxALB, BC SHBxALB, CHB, F1 ALBxCHB, F1 SHBxREB, F2 ALBxCHB, SHB
pal_all_ordered <-c("#0C59FEFF", "#D4B855", "magenta", "#FEC700FF", "darkgray", "#9933CC","black", "#15E0FAFF")


#####STRUCTURE PLOTS######

# Stacked percent Structure plot

Lanier_structure_by_site <- 
  ggplot(struplotdat, aes(factor(samp_ord), q, fill = factor(cluster))) +
  facet_grid(~Pop, scales = "free", switch = "x", space = "free") +
  labs(x = "Individuals by Site", y = "Proportional Assignment (q)", title = "A") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  #*changes to theme for plotting appearance*
  geom_col(color=NA, linewidth=0) +
  theme_minimal(base_size=20) + 
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
      #Make title font bold
    plot.title = element_text(face = "bold", size = 20),
    axis.title.x = element_text(color = "black", size = 18),
    axis.title.y = element_text(color = "black", size = 18),
    axis.text = element_text(color = "black", size = 16),
    strip.text.x = element_text(size=16)
  ) +
  scale_fill_manual(values = pal_q, labels = c("ALB", "CHB", "FLB/LMB", "REB", "SHB")) +
  labs(fill = "Taxon") 
Lanier_structure_by_site


#Now lets see facet grid groups by nh_class

#Need dataframe with "N/A" for those not assigned at 0.5 threshold in NewHybrids
struplotdatv2 <- struplotdat
struplotdatv2["class_nh_50"][struplotdatv2["class_nh_50"] == ""] <-"N/A"
struplotdatv2
struplotdat

Lanier_structure_by_class <- 
  ggplot(struplotdatv2, aes(factor(samp_ord), q, fill = factor(cluster))) +
  facet_grid(~class_nh_50, scales = "free", switch = "x", space = "free") +
  theme_minimal() + labs(x = "Individuals by Hybrid Class", y = "Proportional Assignment (q)", title = "B") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  #*changes to theme for plotting appearance*
  geom_col(color=NA, linewidth=0) +
  theme_minimal(base_size=20) + 
  theme(
    panel.spacing.x = unit(1, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    #Make title font bold
    plot.title = element_text(face = "bold", size = 20),
    axis.title.x = element_text(color = "black", size = 18),
    axis.title.y = element_text(color = "black", size = 18),
    axis.text = element_text(color = "black", size = 16),
    #Change angle of facet labels
    strip.text.x = element_text(angle=45, size=12),
    #don't cut off BCSMB label!!!
    strip.clip = "off"
  ) +
  scale_fill_manual(values = pal_q, labels = c("ALB", "CHB", "FLB/LMB", "REB", "SHB")) +
  labs(fill = "Taxon") 
Lanier_structure_by_class

  
#Plotting NewHybrids P of Z's to place in line with Structure assignments for comparison
#Bring in those data

Lanier_classes_by_site1 <-
  ggplot(NHplotdat1, aes(factor(samp_ord), PofZ, fill = factor(class))) +
  geom_col(color = "NA", size = 0) +
  facet_grid(~Pop, scales = "free", switch = "x", space = "free") +
  theme_minimal() + labs(x = "Individuals by Site", y = "P of Z", title = "A") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  #*changes to theme for plotting appearance*
  geom_col(color=NA, linewidth=0) +
  theme_minimal(base_size=16) + 
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    #Make title font bold
    plot.title = element_text(face = "bold", size = 20),
    axis.title.x = element_text(color = "black", size = 18),
    axis.title.y = element_text(color = "black", size = 18),
    axis.text = element_text(color = "black", size = 12)
  ) +
  scale_fill_manual(values = pal_1.2, labels = c("ALB", "BCALB", "F1", "F2", "BCSHB", "SHB")) +
  labs(fill = "Hybrid Class") 
Lanier_classes_by_site1

Lanier_classes_by_site2 <-
  ggplot(NHplotdat2, aes(factor(samp_ord), PofZ, fill = factor(class))) +
  geom_col(color = "NA", size = 0) +
  facet_grid(~Pop, scales = "free", switch = "x", space = "free") +
  theme_minimal() + labs(x = "Individuals by Site", y = "P of Z", title = "B") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  #*changes to theme for plotting appearance*
  geom_col(color=NA, linewidth=0) +
  theme_minimal(base_size=16) + 
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    #Make title font bold
    plot.title = element_text(face = "bold", size = 20),
    axis.title.x = element_text(color = "black", size = 18),
    axis.title.y = element_text(color = "black", size = 18),
    axis.text = element_text(color = "black", size = 12)
  ) +
  scale_fill_manual(values = pal_2, labels = c("ALB", "BCALB", "F1", "F2", "BCCHB", "CHB")) +
  labs(fill = "Hybrid Class") 
Lanier_classes_by_site2

Lanier_classes_by_site3 <-
  ggplot(NHplotdat3, aes(factor(samp_ord), PofZ, fill = factor(class))) +
  geom_col(color = "NA", size = 0) +
  facet_grid(~Pop, scales = "free", switch = "x", space = "free") +
  theme_minimal() + labs(x = "Individuals by Site", y = "P of Z", title = "C") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  #*changes to theme for plotting appearance*
  geom_col(color=NA, linewidth=0) +
  theme_minimal(base_size=16) + 
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    #Make title font bold
    plot.title = element_text(face = "bold", size = 20),
    axis.title.x = element_text(color = "black", size = 18),
    axis.title.y = element_text(color = "black", size = 18),
    axis.text = element_text(color = "black", size = 12)
  ) +
  scale_fill_manual(values = pal_3, labels = c("SHB", "BCSHB", "F1", "F2", "BCREB", "REB")) +
  labs(fill = "Hybrid Class") 
Lanier_classes_by_site3


#NOTE:  Need to fix above color palette, maybe with override : https://aosmith.rbind.io/2020/07/09/ggplot2-override-aes/
#Then make other two plots match colors


#####TRIANGLE PLOT SHOWING NEWHYBRIDS CLASSES######


library(dplyr)
library(forcats)

#Making df subsets for each newhybrids comparison grouping of two
dat_comp1 <- subset(dat2, Newhyb.Run1. == '1')
dat_comp2 <- subset(dat2, Newhyb.Run2. == '2')
dat_comp3 <- subset(dat2, Newhyb.Run3. == '3')


dat_comp1 %>% mutate(name = fct_relevel(Newhybs1_assign, "SHB", "BC SHB x ALB", "ALB"))
dat_comp2 %>% mutate(name = fct_relevel(Newhybs2_assign, "CHB", "BC CHB x ALB", "F1 ALB x CHB", "F2 ALB x CHB", "ALB"))
dat_comp3 %>% mutate(name = fct_relevel(Newhybs3_assign, "SHB", "F1 SHB x REB"))


# Create a ggplot object

library(ggh4x)
# Create a ggplot object

triplot1 <- ggplot() +
  
  # Add a scatterplot layer of Introgress Data that is color-coded by NewHybrids class assignments
  #See this example code: https://stackoverflow.com/questions/70919700/ggplot-how-to-assign-both-color-and-shape-for-one-factor-and-also-shape-for-an
  geom_point(data = dat_comp1, size=3, alpha=0.75, aes(x = X1_h_introgress, y = X1_in_het_introgress, color = fct_relevel(Newhybs1_assign, "SHB", "BC SHB x ALB", "ALB"), shape = fct_relevel(Newhybs1_assign, "SHB", "BC SHB x ALB", "ALB"), group = fct_relevel(Newhybs1_assign, "SHB", "BC SHB x ALB", "ALB"))) +
  #scale_shape_manual(values=c(16, 15, 17, 15, 16), labels=c("NEO","BCNEO", "F2", "BCSMB", "SMB"), drop = FALSE) +
  scale_shape_manual(values=c(16, 15, 16), labels=c("SHB", "BCSHB", "ALB"), drop = FALSE) +
  scale_color_manual(values=pal_1.2_tri, labels=c("SHB", "BCSHB", "ALB"), drop = FALSE) +  
  scale_fill_manual(values=pal_1.2_tri, labels=c("SHB", "BCSHB", "ALB"), drop = FALSE) +
  
  #c("black","black","black","black","black")
  # Plot line1 & line 2 bounds
  geom_line(data = line1, aes(x = x, y = y), color = "darkgray", alpha = 0.75, size = 1) +
  geom_line(data = line2, aes(x = x, y = y), color = "darkgray", alpha = 0.75, size = 1) +
  
  
  ggtitle("A") +
  xlab("Hybrid Index") +
  ylab("Interspecific Heterozygosity") +
  theme(panel.background = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
  theme(axis.ticks = element_line(color = "black")) +
  theme(legend.key = element_rect(fill = "NA", color="NA")) +
  theme(legend.position = "bottom") +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size=12)) + # mod to get full legend to fit in final plotted area
  theme(legend.spacing.x = unit(1, "pt")) +
  #Keeping plots square to preserve visual interpretation of triangle
  theme(aspect.ratio = 1) +
  #*changes to theme for plotting appearance*
  theme(base_size=16) + 
  theme(
    #Make title font bold
    plot.title = element_text(face = "bold", size = 20),
    axis.title.x = element_text(color = "black", size = 18),
    axis.title.y = element_text(color = "black", size = 18),
    axis.text = element_text(color = "black", size = 12)
  ) +
  ggh4x::force_panelsizes(rows = unit(3, "in"),
                          cols = unit(3.5, "in"))


triplot1

triplot2 <- ggplot() +
  geom_point(data = dat_comp2, size=3, alpha=0.75, aes(x = X2_h_introgress, y = X2_in_het_introgress, color = fct_relevel(Newhybs2_assign, "CHB", "BC CHB x ALB", "F1 ALB x CHB", "F2 ALB x CHB", "ALB"), shape = fct_relevel(Newhybs2_assign, "CHB", "BC CHB x ALB", "F1 ALB x CHB", "F2 ALB x CHB", "ALB"), group = fct_relevel(Newhybs2_assign, "CHB", "BC CHB x ALB", "F1 ALB x CHB", "F2 ALB x CHB", "ALB"))) +
  scale_shape_manual(values=c(16, 15, 17, 17, 16), labels=c("CHB", "BCCHB", "F1", "F2", "ALB"), drop = FALSE) +
  scale_color_manual(values=pal_2_tri, labels=c("CHB", "BCCHB", "F1", "F2", "ALB"), drop = FALSE) +  
  scale_fill_manual(values=pal_2_tri, labels=c("CHB", "BCCHB", "F1", "F2", "ALB"), drop = FALSE) +
  
  # Plot line1 & line 2 bounds
  geom_line(data = line1, aes(x = x, y = y), color = "darkgray", alpha = 0.75, size = 1) +
  geom_line(data = line2, aes(x = x, y = y), color = "darkgray", alpha = 0.75, size = 1) +
  
  ggtitle("B") +
  xlab("Hybrid Index") +
  ylab("Interspecific Heterozygosity") +
  theme(panel.background = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
  theme(axis.ticks = element_line(color = "black")) +
  theme(legend.key = element_rect(fill = "NA", color="NA")) +
  theme(legend.position = "bottom") +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size=11)) + # mod to get full legend to fit in final plotted area
  theme(legend.spacing.x = unit(1, "pt")) +
  #Keeping plots square to preserve visual interpretation of triangle
  theme(aspect.ratio = 1) +
  #*changes to theme for plotting appearance*
  theme(base_size=16) + 
  theme(
    #Make title font bold
    plot.title = element_text(face = "bold", size = 20),
    axis.title.x = element_text(color = "black", size = 18),
    axis.title.y = element_text(color = "black", size = 18),
    axis.text = element_text(color = "black", size = 12)
  ) +
  ggh4x::force_panelsizes(rows = unit(3, "in"),
                          cols = unit(3.5, "in"))
  
triplot2


triplot3 <- ggplot() +
  geom_point(data = dat_comp3, size=3, alpha=0.75, aes(x = X3_h_introgress, y = X3_in_het_introgress, color = fct_relevel(Newhybs3_assign, "SHB", "F1 SHB x REB"), shape = fct_relevel(Newhybs3_assign, "SHB", "F1 SHB x REB"), group = fct_relevel(Newhybs3_assign, "SHB", "F1 SHB x REB"))) +
  scale_shape_manual(values=c(16, 17), labels=c("SHB", "F1"), drop = FALSE) +
  scale_color_manual(values=pal_3_tri, labels=c("SHB", "F1"), drop = FALSE) +  
  scale_fill_manual(values=pal_3_tri, labels=c("SHB", "F1"), drop = FALSE) +
  
  # Plot line1 & line 2 bounds
  geom_line(data = line1, aes(x = x, y = y), color = "darkgray", alpha = 0.75, size = 1) +
  geom_line(data = line2, aes(x = x, y = y), color = "darkgray", alpha = 0.75, size = 1) +
  
  ggtitle("C") +
  xlab("Hybrid Index") +
  ylab("Interspecific Heterozygosity") +
  theme(panel.background = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
  theme(axis.ticks = element_line(color = "black")) +
  theme(legend.key = element_rect(fill = "transparent", color="NA")) +
  theme(legend.position = "bottom") +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size=12)) + # mod to get full legend to fit in final plotted area
  theme(legend.spacing.x = unit(1, "pt")) +
  #Keeping plots square to preserve visual interpretation of triangle
  theme(aspect.ratio = 1) +
  #*changes to theme for plotting appearance*
  theme(base_size=16) + 
  theme(
    #Make title font bold
    plot.title = element_text(face = "bold", size = 20),
    axis.title.x = element_text(color = "black", size = 18),
    axis.title.y = element_text(color = "black", size = 18),
    axis.text = element_text(color = "black", size = 12)
  ) +
  ggh4x::force_panelsizes(rows = unit(3, "in"),
                          cols = unit(3.5, "in"))

triplot3


#####PIE CHARTS BY SITE with HYBRID CLASSES######
#Showing Proportion of individuals belonging to HYBRID CLASSES
# Calculate the proportion of individuals per site and factor level
library(dplyr)

#Calculating relative frequency aka proportions of individuals of each hybrid class by site
dat2$class_nh_50 <- as.factor(dat2$class_nh_50)
dat2$Pop <- as.factor (dat2$Pop)
proportions <- dat2 %>%
  group_by(Pop, class_nh_50) %>%
  summarise(n_class = n()) %>%
  mutate(freq = n_class / sum(n_class))

# Print the new data frame
print(proportions)

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size = 16),
    legend.text = element_text(size=14)
  )

pie_charts <- list()
for (Pop in unique(proportions$Pop)) {
  data_subset <- proportions[proportions$Pop == Pop,]
  pie_charts[[Pop]] <- ggplot(data_subset, aes(x = "", y = freq, fill = class_nh_50)) +
    geom_bar(stat="identity", width=0.05) +
    coord_polar("y", start=0) +
    labs(title = paste0("Site ", Pop)) +
    scale_fill_manual(values = pal_all_ordered, drop = FALSE) +
    blank_theme +
    theme(axis.text.x=element_blank(), legend.title=element_blank()) +
   #black borders between colors
    geom_col(color = "black")
}


  
#####PIE CHARTS BY SITE with STRUCTURE PROPORTIONS######
#Showing Population level proportional assignment to taxonomic clusters w/ Structure outputs
# Calculate the proportion of alleles across individuals by site 
library(dplyr)

#Calculating mean allelic proportions of all individuals per site ("Q")
stru_pie_dat <- NHplotdat %>%
  group_by(Pop) %>%
  summarise_at(.vars = vars(FLB_LMB_stru, ALB_stru, REB_stru, CHB_stru, SHB_stru),
             .funs = c(mean="mean"))
stru_pie_dat

#Make data into long format for ease of plotting
stru_pie_dat <- tidyr::pivot_longer(stru_pie_dat, cols = c(FLB_LMB_stru_mean, ALB_stru_mean, REB_stru_mean, CHB_stru_mean, SHB_stru_mean),names_to = "cluster", values_to = "Q")

#For plotting, keeping same theme "blank_theme" as in previous pie chart
#But should adopt same cluster colors as used previously, aka q_pallete
Stru_pie_charts <- list()
for (Pop in unique(stru_pie_dat$Pop)) {
  data_subset2 <- stru_pie_dat[stru_pie_dat$Pop == Pop,]
  Stru_pie_charts[[Pop]] <- ggplot(data_subset2, aes(x = "", y = Q, fill = cluster)) +
    geom_bar(stat="identity", width=0.05) +
    coord_polar("y", start=0) +
    labs(title = paste0("Site ", Pop)) +
    scale_fill_manual(values = pal_q, labels = c("ALB", "CHB", "FLB/LMB", "REB", "SHB"), drop = FALSE) +
    blank_theme +
    theme(axis.text.x=element_blank(), legend.title=element_blank()) +
    #black borders between colors
    geom_col(color = "black")
}

#####Combine panels into Single Plots#####
library(patchwork)
library(gridExtra)

#Triangle plots showing introgress results paired with NH classes and Structure q values
#NOTE: Patchwork throws error because of sizing to keep aspect ratios intact, use gridExtra

TRI2 <-gridExtra::grid.arrange(triplot1, triplot2, triplot3, ncol=1)
TRI2

ggsave("/Users/andrewtaylor/Desktop/R_Projects/R_impound_IntrogressNPlots/Lanier_triangleplots.pdf", 
       plot = TRI2,
       width = 6,
       height = 14, 
       units=c("in"))


#Structure-like bar plots showing assignments of individuals
PLOTZ2 <-gridExtra::grid.arrange(Lanier_structure_by_site, Lanier_structure_by_class, ncol=1)
PLOTZ2

ggsave("/Users/andrewtaylor/Desktop/R_Projects/R_impound_IntrogressNPlots/Lanier_Structure_plots.pdf", 
       plot = PLOTZ2,
       width = 11,
       height = 8.5, 
       units=c("in"))


ClassPlot2 <-gridExtra::grid.arrange(Lanier_classes_by_site1, Lanier_classes_by_site2, Lanier_classes_by_site3, ncol=1)
ClassPlot2

ggsave("/Users/andrewtaylor/Desktop/R_Projects/R_impound_IntrogressNPlots/Lanier_Hybrid_classes_wStru_plots.pdf", 
       plot = ClassPlot2,
       width = 11,
       height = 8.5, 
       units=c("in"))

#Pies showing Hybrid class frequency\
Pies_class_organized2 <-gridExtra::grid.arrange(Lanier_classes_by_site1, Lanier_classes_by_site2, Lanier_classes_by_site3, ncol=1)
ClassPlot2

library(patchwork)

Pies_class_organized <- pie_charts[[4]] + pie_charts[[8]] +  
                        pie_charts[[3]] + pie_charts[[7]] + 
                        pie_charts[[2]] + pie_charts[[6]] + 
                        pie_charts[[1]] + pie_charts[[5]] + 
  plot_layout(nrow=4, ncol=2) + guides(color="collect") + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')

Pies_class_organized

pdf(file = "/Users/andrewtaylor/Desktop/R_Projects/R_impound_IntrogressNPlots/Lanier_Pies_classes_organized.pdf",
    width = 8.5,
    height = 11) 
extrafont::loadfonts(device="pdf")
pdfFonts('ArialMT')
Pies_class_organized
dev.off()

Stru_pies_organized <- Stru_pie_charts[[4]] + Stru_pie_charts[[8]] +  
  Stru_pie_charts[[3]] + Stru_pie_charts[[7]] + 
  Stru_pie_charts[[2]] + Stru_pie_charts[[6]] + 
  Stru_pie_charts[[1]] + Stru_pie_charts[[5]] + 
  plot_layout(nrow=4, ncol=2) + guides(color="collect") + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')

Stru_pies_organized

pdf(file = "/Users/andrewtaylor/Desktop/R_Projects/R_impound_IntrogressNPlots/Lanier_Pies_stru_organized.pdf",
    width = 8.5,
    height = 11) 
extrafont::loadfonts(device="pdf")
pdfFonts('ArialMT')
Stru_pies_organized
dev.off()



dev.off()
