######################################
#
# BUSCO summary figure
# @version 4.0.0
# @since BUSCO 2.0.0
# 
# Copyright (c) 2015-2024, Evgeny Zdobnov (ez@ezlab.org)
# Licensed under the MIT license. See LICENSE.md file.
#
######################################

# Load the required libraries
library(ggplot2)
library("grid")

#install.packages("extrafont")
library(extrafont)


# Import system fonts into R
#font_import()

# Register fonts for use in the PDF device
#loadfonts(device = "pdf")

# Now try saving the plot again using the specified font
#my_family <- "Arial"  # or any other system-installed font

# !!! CONFIGURE YOUR PLOT HERE !!! 
# Output
#my_output <- paste("./busco_summary/","busco_figure.png",sep="/")
setwd("/home/yzliu/eDNA/faststorage/yzliu/DK_proj/data/busco/hifi_asm")
my_output <- paste("./busco_summary/","busco_figure.pdf",sep="/")  
my_width <- 20
my_height <- 15
my_unit <- "cm"

# Colors
my_colors <- c("#56B4E9", "#3492C7", "#F0E442", "#F04442")
# Bar height ratio
my_bar_height <- 0.75

# Legend
my_title <- "BUSCO Assessment Results"

# Font
my_family <- "sans"
#my_family <- "Arial"
my_size_ratio <- 1

# !!! SEE YOUR DATA HERE !!! 
# Your data as generated by python, remove or add more
#my_species <- c('BUSCO_Bombus_veteranus_pl_2', 'BUSCO_Bombus_veteranus_pl_2', 'BUSCO_Bombus_veteranus_pl_2', 'BUSCO_Bombus_veteranus_pl_2', 'BUSCO_Bombus_veteranus_pl_0', 'BUSCO_Bombus_veteranus_pl_0', 'BUSCO_Bombus_veteranus_pl_0', 'BUSCO_Bombus_veteranus_pl_0', 'BUSCO_Bombus_veteranus_pl_1', 'BUSCO_Bombus_veteranus_pl_1', 'BUSCO_Bombus_veteranus_pl_1', 'BUSCO_Bombus_veteranus_pl_1')
#my_species <- c(expression("BUSCO_",italic('Bombus_veteranus'),"_pl_2"), 
                #expression("BUSCO_",italic('Bombus_veteranus'),"_pl_2"),
                #expression("BUSCO_",italic('Bombus_veteranus'),"_pl_2"),
                #expression("BUSCO_",italic('Bombus_veteranus'),"_pl_2"),
                #expression("BUSCO_",italic('Bombus_veteranus'),"_pl_1"),
                #expression("BUSCO_",italic('Bombus_veteranus'),"_pl_1"),
                #expression("BUSCO_",italic('Bombus_veteranus'),"_pl_1"),
                #expression("BUSCO_",italic('Bombus_veteranus'),"_pl_1"),
                #expression("BUSCO_",italic('Bombus_veteranus'),"_pl_0"),
                #expression("BUSCO_",italic('Bombus_veteranus'),"_pl_0"),
                #expression("BUSCO_",italic('Bombus_veteranus'),"_pl_0"),
                #expression("BUSCO_",italic('Bombus_veteranus'),"_pl_0"))
# above codes not working properly

my_species <- factor(my_species)
my_species <- factor(my_species,levels(my_species)[c(length(levels(my_species)):1)]) # reorder your species here just by changing the values in the vector :
my_percentage <- c(97.7, 0.6, 0.7, 1.0, 16.2, 82.3, 0.7, 0.8, 95.5, 2.8, 0.7, 1.0)
my_values <- c(5856, 38, 41, 56, 968, 4932, 40, 51, 5724, 170, 40, 57)

######################################
######################################
######################################
# Code to produce the graph
labsize = 1
if (length(levels(my_species)) > 10){
 labsize = 0.66
}
print("Plotting the figure ...")
category <- c(rep(c("S","D","F","M"),c(1)))
category <-factor(category)
category = factor(category,levels(category)[c(4,1,2,3)])
df = data.frame(my_species,my_percentage,my_values,category)

figure <- ggplot() + 
  
  geom_bar(aes(y = my_percentage, x = my_species, fill = category), position = position_stack(reverse = TRUE), data = df, stat="identity", width=my_bar_height) + 
  coord_flip() + 
  theme_gray(base_size = 8) + 
  scale_y_continuous(labels = c("0","20","40","60","80","100"), breaks = c(0,20,40,60,80,100)) + 
  scale_fill_manual(values = my_colors,labels =c(" Complete (C) and single-copy (S)  ",
                                                 " Complete (C) and duplicated (D)",
                                                 " Fragmented (F)  ",
                                                 " Missing (M)")) +   
  ggtitle(my_title) + 
  xlab("Purging level of assemblies") + 
  ylab("\n%BUSCOs") + 

# Italicizing species names in y-axis
  scale_x_discrete(labels = c(expression(italic('Bombus veteranus') ~ 'pl_2'), 
                              expression(italic('Bombus veteranus') ~ 'pl_1'), 
                              expression(italic('Bombus veteranus') ~ 'pl_0'))) + 

  theme(plot.title = element_text(family=my_family, hjust=0.5, colour = "black", size = rel(2.2)*my_size_ratio, face = "bold")) + 
  theme(legend.position="top",legend.title = element_blank()) + 
  theme(legend.text = element_text(family=my_family, size = rel(1.2)*my_size_ratio)) + 
  theme(panel.background = element_rect(color="#FFFFFF", fill="white")) + 
  theme(panel.grid.minor = element_blank()) + 
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.y = element_text(family=my_family, colour = "black", size = rel(1.66)*my_size_ratio)) + 
  theme(axis.text.x = element_text(family=my_family, colour = "black", size = rel(1.66)*my_size_ratio)) + 
  theme(axis.line = element_line(size=1*my_size_ratio, colour = "black")) + 
  theme(axis.ticks.length = unit(.85, "cm")) + 
  theme(axis.ticks.y = element_line(colour="white", size = 0)) + 
  theme(axis.ticks.x = element_line(colour="#222222")) + 
  theme(axis.ticks.length = unit(0.4, "cm")) + 
  theme(axis.title = element_text(family=my_family, size=rel(1.5)*my_size_ratio)) + 
  
  guides(fill = guide_legend(override.aes = list(colour = NULL))) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))
  
  for(i in rev(c(1:length(levels(my_species))))){
    detailed_values <- my_values[my_species==my_species[my_species==levels(my_species)[i]]]
    total_buscos <- sum(detailed_values)
    figure <- figure + 
    annotate("text", label=paste("C:", detailed_values[1] + detailed_values[2], " [S:", detailed_values[1], ", D:", detailed_values[2], "], F:", detailed_values[3], ", M:", detailed_values[4], ", n:", total_buscos, sep=""), 
             y=3, x = i, size = labsize*4*my_size_ratio, colour = "black", hjust=0, family=my_family)
  }
  
ggsave(figure, file=my_output, width = my_width, height = my_height, unit = my_unit)
print("Done")
