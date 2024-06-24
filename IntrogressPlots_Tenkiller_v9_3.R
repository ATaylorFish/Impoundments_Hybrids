###Code to visualize outputs from Structure, NewHybrids, and Introgress in integrative ways
##Primarily using packages dplyr, ggplot2, patchwork, and fishualize palettes
#Written A. Taylor 8/15/23 - 5/22/24, help from various linked authors


#####IMPORT DATA#####
##### Import the combined data to dataframe & clean
dat <- read.csv("/Users/andrewtaylor/Desktop/R_Projects/R_impound_IntrogressNPlots/InputData/Tenkiller_combined.csv")

#Remove refs to keep only wild caught fish
dat <- subset(dat, SampType == 'wildpop')

##IF we want to remove the specimens that NewHybrids failed to classify for plotting
dat2 <- subset(dat, class_nh_50 != "")

####Create new dataframes to hold the bounding lines for triangle plot
# Create a dataframe
line1 <- data.frame(x = seq(0, 0.5, by = 0.01), y = seq(0, 1, by = 0.02))
line2 <- data.frame(x = seq(0.5, 1, by = 0.01), y = seq(1, 0, by = -0.02))

#####Importing NewHybrids P of Z's for structure-like plots later
NHdat <- read.csv("/Users/andrewtaylor/Desktop/R_Projects/R_impound_IntrogressNPlots/InputData/Tenkiller_NHPofZ.csv", na.strings = "")




#####COLOR PALETTES & PLOTTING MADNESS####
library(ggplot2)


# Could use RColorBrewer or Fishualize, for example

# Create a divergent color palette for hybrid classes
#In this case ordered "NEO","BCNEO","F2", "BCSMB", "SMB"
library(fishualize)
pal1 <- fish(n = 5, option = "Lepomis_megalotis")
pal1

# Create a divergent color palette for hybrid classes PofZ's
#In this case ordered 'NEO', 'BCNEO', 'F1', 'F2', 'BCSMB', 'SMB'
nh_palette <- fish(n = 6, option = "Lepomis_megalotis")

# Create a divergent set for the clusters in Structure assignments
# in this case, ordered NEO, SMB, SPB
# pulling from pal1 and a random complementary color for SPB
q_palette <- c("#56C0BFFF", "#5E1A0AFF","pink")

#####TRIANGLE PLOT SHOWING NEWHYBRIDS CLASSES######


library(dplyr)
library(forcats)

dat3 <- dat2 %>%
  mutate(name = fct_relevel(class_nh_50, "NEO", "BCNEO", "F2", "BCSMB", "SMB"))

# Create a ggplot object
library(ggh4x)

triplot <- ggplot() +
  
  # Add a scatterplot layer of Introgress Data that is color-coded by NewHybrids class assignments
  #See this example code: https://stackoverflow.com/questions/70919700/ggplot-how-to-assign-both-color-and-shape-for-one-factor-and-also-shape-for-an
  geom_point(data = dat3, size=3, alpha=0.75, aes(x = h_introgress, y = in_het_introgress, color = fct_relevel(class_nh_50, "NEO", "BCNEO", "F2", "BCSMB", "SMB"), shape = fct_relevel(class_nh_50, "NEO", "BCNEO", "F2", "BCSMB", "SMB"), group = fct_relevel(class_nh_50, "NEO", "BCNEO", "F2", "BCSMB", "SMB"))) +
  scale_shape_manual(values=c(16, 15, 17, 15, 16), labels=c("NEO","BCNEO", "F1/F2", "BCSMB", "SMB"), drop = FALSE) +
  scale_color_manual(values=pal1, labels=c("NEO","BCNEO", "F1/F2", "BCSMB", "SMB"), drop = FALSE) +  
  scale_fill_manual(values=pal1, labels=c("NEO","BCNEO", "F1/F2", "BCSMB", "SMB"), drop = FALSE) +
  
  #c("black","black","black","black","black")
  # Plot line1 & line 2 bounds
  geom_line(data = line1, aes(x = x, y = y), color = "darkgray", alpha = 0.75, size = 1) +
  geom_line(data = line2, aes(x = x, y = y), color = "darkgray", alpha = 0.75, size = 1) +
  
  #Point appearance, adding shape and size for visual display
  
  
  # Add a title - don't need this for publication figures
  #ggtitle("Triangle Plot with HewHybrids Classes") +
  
  # Add a x-axis label
  ggtitle("A") +
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

triplot

###NOTES

#Figure out how to size points and add black borders to shapes
#Interpret all f2 as f1 or f2 since power tests showed weaker ID here.


#####TRIANGLE PLOT SHOWING STRUCTURE PROPORTIONS#####
# Create a ggplot object
triplot_stru <- ggplot() +
  
  # Add a scatterplot layer of Introgress Data that is a gradient showing proportion NEO from Structure output 
  geom_point(data = dat2, size=3, alpha=0.75, aes(x = h_introgress, y = in_het_introgress, color = SMB_stru)) +
  
  scale_colour_gradient2(
    low = "#56C0BFFF",
    mid = "#F6D25CFF",
    high = "#5E1A0AFF",
    midpoint = 0.50,
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "colour") +
  
  # Plot line1 & line 2 bounds
  geom_line(data = line1, aes(x = x, y = y), color = "darkgray", alpha = 0.75, size = 1) +
  geom_line(data = line2, aes(x = x, y = y), color = "darkgray", alpha = 0.75, size = 1) +
  
  #Add transparency to points for overlap visual display
  geom_point(alpha=0.5) +
  
  # Add a title - don't need this for publication figures
  #ggtitle("Triangle Plot with HewHybrids Classes") +
  
  # Add a x-axis label
  ggtitle("B") +
  xlab("Hybrid Index") +
  ylab("Interspecific Heterozygosity") +
  theme(panel.background = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
  theme(axis.ticks = element_line(color = "black")) +
  theme(legend.key = element_rect(fill = "NA", color="NA")) +
  theme(legend.position = "bottom") +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size=10)) + # mod to get full legend to fit in final plotted area
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
  
triplot_stru



#####STRUCTURE PLOT DATA PREP#####

#With help from Arregoitia code here: https://luisdva.github.io/rstats/model-cluster-plots/
#With help from French code too: https://connor-french.github.io/intro-pop-structure-r/


#First need to mutate data so that individuals have multiple rows, each with STRU val
#Select columns to keep
struplotdat <- subset(dat, select = c(Samp, Pop, class_nh_50, SPB_stru, SMB_stru, NEO_stru))

#Sort by NEO high to low within each population first (for plotting viz)
library(dplyr)
struplotdat <- dplyr::arrange(struplotdat, Pop, NEO_stru)

#see idea here for idea on how to order: https://youtu.be/yZSWFND3ca0
#Make a samp_order attribute as a column here
struplotdat$samp_ord <- row.names(struplotdat)
#checking structure of data
str(struplotdat)
#changing samp_ord from character to numeric
struplotdat$samp_ord <- as.numeric(struplotdat$samp_ord)


#reshaping data for plotting (each row now has one structure prob per cluster)
library(tidyr)
struplotdat <- tidyr::pivot_longer(struplotdat, cols = c(SPB_stru, SMB_stru, NEO_stru), names_to = "cluster", values_to = "q")

#NH P of Z's for later use in plotting
NHplotdat <-dplyr::right_join(dat,NHdat, by = "Samp")
NHplotdat <- dplyr::arrange(NHplotdat, Pop, NEO_stru)
NHplotdat$samp_ord <- as.numeric(row.names(NHplotdat))
NHplotdat <- tidyr::pivot_longer(NHplotdat, cols = c(NEO, SMB, F1, F2, BCNEO, BCSMB), names_to = "class", values_to = "PofZ")

#####STRUCTURE PLOTS######

# Stacked percent Structure plot

Tenkiller_structure_by_site <- 
  ggplot(struplotdat, aes(factor(samp_ord), q, fill = factor(cluster))) +
  geom_col(color = "NA", size = 0) +
  facet_grid(~Pop, scales = "free", switch = "x", space = "free") +
  theme_minimal() + labs(x = "Individuals by Site", y = "Proportional Assignment (q)", title = "A") +
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
    axis.title.x = element_text(color = "black", size = 16),
    axis.title.y = element_text(color = "black", size = 14),
    axis.text = element_text(color = "black", size = 12)
  ) +
  scale_fill_manual(values = q_palette, labels = c("NEO", "SMB", "SPB")) +
  labs(fill = "Taxon") 
Tenkiller_structure_by_site


#Now lets see facet grid groups by nh_class

#Need dataframe with "N/A" for those not assigned at 0.5 threshold in NewHybrids
struplotdatv2 <- struplotdat
struplotdatv2["class_nh_50"][struplotdatv2["class_nh_50"] == ""] <-"N/A"
struplotdatv2
struplotdat

Tenkiller_structure_by_class <- 
  ggplot(struplotdatv2, aes(factor(samp_ord), q, fill = factor(cluster))) +
  geom_col(color = "NA", size = 0) +
  facet_grid(~class_nh_50, scales = "free", switch = "x", space = "free") +
  theme_minimal() + labs(x = "Individuals by Hybrid Class", y = "Proportional Assignment (q)", title = "C") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  #*changes to theme for plotting appearance*
  geom_col(color=NA, linewidth=0) +
  theme_minimal(base_size=16) + 
  theme(
    panel.spacing.x = unit(0.2, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    #Make title font bold
    plot.title = element_text(face = "bold", size = 20),
    axis.title.x = element_text(color = "black", size = 16),
    axis.title.y = element_text(color = "black", size = 14),
    axis.text = element_text(color = "black", size = 12),
    #Change angle of facet labels
    strip.text.x = element_text(angle=45),
    #don't cut off BCSMB label!!!
    strip.clip = "off"
  ) +
  scale_fill_manual(values = q_palette, labels = c("NEO", "SMB", "SPB")) +
  labs(fill = "Taxon") 
Tenkiller_structure_by_class

  
#Plotting NewHybrids P of Z's to place in line with Structure assignments for comparison
#Bring in those data

Tenkiller_classes_by_site <-
  ggplot(NHplotdat, aes(factor(samp_ord), PofZ, fill = factor(class))) +
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
    axis.title.x = element_text(color = "black", size = 16),
    axis.title.y = element_text(color = "black", size = 16),
    axis.text = element_text(color = "black", size = 12)
  ) +
  #scale_fill_manual(values = nh_palette, labels = c("BCNEO", "BCSMB", "F1", "F2", "NEO", "SMB")) +
  scale_fill_manual(values = nh_palette, breaks=c('NEO', 'BCNEO', 'F1', 'F2', 'BCSMB', 'SMB')) +
  labs(fill = "Hybrid Class") 
Tenkiller_classes_by_site



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
    scale_fill_manual(values = nh_palette, breaks=c('NEO', 'BCNEO', 'F1', 'F2', 'BCSMB', 'SMB'), drop = FALSE) +
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
  summarise_at(.vars = vars(SPB_stru, SMB_stru, NEO_stru),
             .funs = c(mean="mean"))
stru_pie_dat
colnames(stru_pie_dat)[2] ="SPB"
colnames(stru_pie_dat)[3] ="SMB"
colnames(stru_pie_dat)[4] ="NEO"
#Make data into long format for ease of plotting
stru_pie_dat <- tidyr::pivot_longer(stru_pie_dat, cols = c(SPB, SMB, NEO), names_to = "cluster", values_to = "Q")

#For plotting, keeping same theme "blank_theme" as in previous pie chart
#But should adopt same cluster colors as used previously, aka q_pallete
Stru_pie_charts <- list()
for (Pop in unique(stru_pie_dat$Pop)) {
  data_subset2 <- stru_pie_dat[stru_pie_dat$Pop == Pop,]
  Stru_pie_charts[[Pop]] <- ggplot(data_subset2, aes(x = "", y = Q, fill = cluster)) +
    geom_bar(stat="identity", width=0.05) +
    coord_polar("y", start=0) +
    labs(title = paste0("Site ", Pop)) +
    scale_fill_manual(values = q_palette, labels = c("NEO", "SMB", "SPB"), drop = FALSE) +
    blank_theme +
    theme(axis.text.x=element_blank(), legend.title=element_blank()) +
    #black borders between colors
    geom_col(color = "black")
}

#####Combine panels into Single Plots#####
library(patchwork)
library(gridExtra)

#Triangle plots showing introgress results paired with NH classes and Structure q values
TRI2 <- gridExtra::grid.arrange(triplot, triplot_stru, ncol=1)
TRI2

ggsave("/Users/andrewtaylor/Desktop/R_Projects/R_impound_IntrogressNPlots/Tenkiller_triangleplots.pdf", 
       plot = TRI2,
       width = 6,
       height = 10, 
       units=c("in"))


#Structure-like bar plots showing assignments of individuals
PLOTZ<-Tenkiller_structure_by_site + Tenkiller_classes_by_site + Tenkiller_structure_by_class + plot_layout(ncol = 1)
PLOTZ  

pdf(file = "/Users/andrewtaylor/Desktop/R_Projects/R_impound_IntrogressNPlots/Tenkiller_Structure_plots.pdf",
  width = 11,
  height = 8.5) 
extrafont::loadfonts(device="pdf")
pdfFonts('ArialMT')
PLOTZ
dev.off()

#Pies showing Hybrid class frequency
library(patchwork)

Pies_class_organized <- pie_charts[[4]] + pie_charts[[8]] + pie_charts[[12]] + 
  pie_charts[[3]] + pie_charts[[7]] + pie_charts[[11]] + 
  pie_charts[[2]] + pie_charts[[6]] + pie_charts[[10]] + 
  pie_charts[[1]] + pie_charts[[5]] + pie_charts[[9]] + 
  plot_layout(nrow=4, ncol=3) + guides(color="collect") + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')

Pies_class_organized

pdf(file = "/Users/andrewtaylor/Desktop/R_Projects/R_impound_IntrogressNPlots/Tenkiller_Pies_classes_organized.pdf",
    width = 8.5,
    height = 11) 
extrafont::loadfonts(device="pdf")
pdfFonts('ArialMT')
Pies_class_organized
dev.off()

Stru_pies_organized <- Stru_pie_charts[[4]] + Stru_pie_charts[[8]] + Stru_pie_charts[[12]] + 
  Stru_pie_charts[[3]] + Stru_pie_charts[[7]] + Stru_pie_charts[[11]] + 
  Stru_pie_charts[[2]] + Stru_pie_charts[[6]] + Stru_pie_charts[[10]] + 
  Stru_pie_charts[[1]] + Stru_pie_charts[[5]] + Stru_pie_charts[[9]] + 
  plot_layout(nrow=4, ncol=3) + guides(color="collect") + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')

Stru_pies_organized

pdf(file = "/Users/andrewtaylor/Desktop/R_Projects/R_impound_IntrogressNPlots/Tenkiller_Pies_stru_organized.pdf",
    width = 8.5,
    height = 11) 
extrafont::loadfonts(device="pdf")
pdfFonts('ArialMT')
Stru_pies_organized
dev.off()

