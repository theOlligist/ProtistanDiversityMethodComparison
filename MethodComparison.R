rm(list = ls())
source("~/Dropbox/R_general_working_directory/Dictionary_of_Colors.R")

setwd("~/Dropbox/SMP_R_working_directory/MetaTranscriptomics/Method_comparison/")
load("~/Dropbox/SMP_R_working_directory/MetaTranscriptomics/Method_comparison/MethodComparisson.RData")
library(tidyverse);library(jishonoiro);library(patchwork);library(vegan);library(reshape2);library(edgeR)
one41= c("#FF5200", "#0D2B52","#A6FF47")
pal1 = c("#8c6510", "#0D2B52","#FF7399")
color_pal = c("black", "light grey", "steel blue")
# Direct enumeration ------------------------------------------------------

#load the data
count_table = read.csv("~/Dropbox/SMP_R_working_directory/Physical_chemical/SMP2018-Cell_counts.csv", header = T)
#count_table %>% view()

##Create figure S1 - Whole community, genus-level 
Tidy_df = count_table %>% 
  #filter(Group == "Diatom") %>% 
  pivot_longer(cols = starts_with("Day"), names_to = "Day", values_to = "value") %>% 
  mutate(Day = str_replace_all(Day, "Day", ""),
         Day = factor(Day, levels = c("1","2","3","4",
                                      "5","6","7","8",
                                      "9","10","11","12",
                                      "13","14","15")),
         Group = factor(Group, levels = c("All.others","Dinoflagellate", "Diatom")),
         value = ifelse(value > 0, value, NA), # I did this to get the zero values to not show a point at all.
         tax_color = case_when(Group == "Diatom" ~"#56AA69", 
                               Group == "Dinoflagellate" ~"#A90636",
                               str_detect(Group, "other") ~"#FFBF6E")) %>% 
  mutate(Genus = factor(Genus, levels = c("Lauderia", "Lioloma", 
                                          "Thalassionema","Pleurosigma", "Coscinodiscus", 
                                          "Odontella", "Entomoneis",   
                                          "Leptocylindrus",  "Hemiaulus", 
                                          "Chaetoceros", "Navicula", "Cylindrotheca", "Thalassiosira","Guinardia", "Pseudo-nitzschia",
                                          "Protoperidinium", "Polykrikos", "Dinophysis", "Oxytoxum",
                                            "Alexandrium", "Gymnodinium", 
                                           "Cochlodinium", "Ceratium", "Lingulodinium", "Gonyaulax", "Scrippsiella","Prorocentrum", 
                                              
                                               
                                          "Euglenoid", "Silicoflagellate","Other ciliate", "Tintinnid", "Oligotrich")))  

#Genus level resolution: supplimental figure
Tidy_df %>% 
  mutate(Group = factor(Group, levels = c("All.others", "Dinoflagellate", "Diatom"))) %>% 
  ggplot(., aes(x = Day, y = Genus, size = value, fill = Group)) +
  geom_point(aes(color = Group), shape = 21) +
  scale_fill_manual(name = NULL,
                     breaks = c( "All.others","Dinoflagellate","Diatom"),
                     values = pal1) +
  scale_color_manual(name = NULL,
                    breaks = c( "All.others","Dinoflagellate","Diatom"),
                    values = pal1) +
  scale_size_area(name = "Cells / mL", 
                  max_size = 5) +
  labs(y = "",
       x = "Day") +
  theme_minimal() +
  theme(panel.grid.major.y = element_line(size = .4),
        panel.grid.major.x = element_line(size = .4)) 

#Save the plot 
ggsave(plot = last_plot(), filename = "~/Dropbox/Project-data/SMPier/Transcriptomics/Figures/Method_comparison/enumeration_wholecommunty_point-SMP2018.pdf", width = 8, height = 6)


#Low resolution stacked barplot: Figure 1
Tidy_df %>% 
  ggplot(., aes(x = Day, y = value, fill = Group)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(breaks = c( "Dinoflagellate", "Diatom","All.others"),
                    values = c(  "#8c6510", "#FF7399", "#0D2B52")) +
  labs(y = "Cells / ml",
       x = "Day") +
  theme_minimal() +
  theme(panel.grid.major.x = element_line(size = 0.4),
        panel.grid.major.y = element_line(size = 0.4))

#Save the plot
ggsave(plot = last_plot(), filename = "~/Dropbox/Project-data/SMPier/Transcriptomics/Figures/Method_comparison/enumeration_bars-SMP2018.pdf", width = 8, height = 6)

# Alternate view of figure 1 - area plot
Fig1.enum = Tidy_df %>% 
  group_by(Day, Group) %>% 
  summarise(value = sum(value, na.rm = TRUE)) %>% 
  ggplot(., aes(x = as.numeric(Day), y = value, fill = Group)) +
  geom_area(stat = "identity", position = "stack", aes(fill = Group)) +
  scale_fill_manual(name = "",
                    breaks = c("All.others","Dinoflagellate", "Diatom"),
                    values = c("black", "light gray", "steel blue")) +
  scale_x_continuous(breaks = c(1:15),
                     labels = c(1:15)) +
  labs(y = "Cells / ml",
       x = "") +
  theme_minimal() +
  theme(panel.grid.major.x = element_line(size = 0.4),
        panel.grid.major.y = element_line(size = 0.4),
        #axis.text.x = element_blank(),
        axis.title.y = element_text(size = 13))

#Alternate view of figure 1 - Lines
Tidy_df %>% 
  group_by(Group, Day) %>% 
  summarise(value = sum(value, na.rm = TRUE)) %>% 
  ungroup() %>% 
  ggplot(., aes(x = Day, y = value, color = Group, group = Group)) +
  geom_line() +
  scale_color_manual(breaks = c( "Dinoflagellate", "Diatom","All.others"),
                    values = c(  "#8c6510", "#FF7399", "#0D2B52")) +
  labs(y = "Cells / ml",
       x = "Day") +
  theme_minimal() +
  theme(panel.grid.major.x = element_line(size = 0.4),
        panel.grid.major.y = element_line(size = 0.4))

#diversity measures
Tidy_df %>% 
  select(Genus, Day, value) %>% 
  dcast(formula = Genus~Day, value.var = .$value)
  
### Stats
StatsTable = Tidy_df %>% 
  group_by(Day, Group) %>% 
  summarise(value = sum(value, na.rm = TRUE)) %>% 
  group_by(Group) %>% 
  summarise(stdev = sd(value),
            var = var(value),
            avg = mean(value),
            dif = max(value)-min(value))

### Sandbox
#daily diatom abundance
diversity_table = count_table %>% select(starts_with("Day"))
DE_shannon<-diversity(diversity_table,index="shannon",2)
#Calculate inverse simpson index
DE_invsimp<-diversity(diversity_table,index="invsimpson",2)
#Calculate the number of ASVs
DE_sp_count <-colSums(diversity_table>0) 


# Amplicon sequencing -----------------------------------------------------

load("~/Dropbox/SMP_R_working_directory/Table_2018.rda")
write.csv(Table_2018, file = "~/Dropbox/Project-data/SMPier/Transcriptomics/Figures/Method_comparison/Tables/Norm_ASVtable.csv")
Table_2018 %>% 
  mutate(low_res = case_when(str_detect(Tax, regex("diatom", ignore_case = T)) ~"Diatom", 
                             str_detect(Tax, regex("dino", ignore_case = T)) ~"Dinoflagellate", 
                             !str_detect(Tax, regex("(diatom)|(dino)", ignore_case = T)) ~"All.others")) %>% 
  pivot_longer(cols = contains("Day"), names_to = "Day", values_to = "value") %>% 
  mutate(Day = str_replace_all(Day, "(Day)|(\\.2018)", ""),
         Day = factor(Day, levels = c("1","2","3","4",
                                      "5","6","7","8",
                                      "9","10","11","12",
                                      "13","14","15")),
         low_res = factor(low_res, levels = c("Dinoflagellate", "Diatom", "All.others")),
         frac = value/sum(value, na.rm = T)) %>% 
  select(Level7, Day, value, frac, Level8, low_res) %>% 
  ggplot(., aes(x = Day, y = value, fill = low_res)) + 
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(breaks = c("All.others", "Dinoflagellate", "Diatom"),
                    values = c(  "#8c6510","#FF7399","#0D2B52" )) +
  labs(y = "Relative abundance",
       x = "Day") +
  theme_minimal() +
  theme(panel.grid.major.x = element_line(size = 0.4),
        panel.grid.major.y = element_line(size = 0.4))

ggsave(plot = last_plot(), filename = "~/Dropbox/Project-data/SMPier/Transcriptomics/Figures/Method_comparison/metabarcode_bars-SMP2018.pdf", width = 8, height = 6)

#Alternate view figure 2/3 - area
Fig1.ASV = Table_2018 %>%
  mutate(low_res = case_when(str_detect(Tax, regex("diatom", ignore_case = T)) ~"Diatom", 
                             str_detect(Tax, regex("dino", ignore_case = T)) ~"Dinoflagellate", 
                             !str_detect(Tax, regex("(diatom)|(dino)", ignore_case = T)) ~"All.others")) %>% 
  pivot_longer(cols = contains("Day"), names_to = "Day", values_to = "value") %>% 
  mutate(Day = str_replace_all(Day, "(Day)|(\\.2018)", ""),
         Day = factor(Day, levels = c("1","2","3","4",
                                      "5","6","7","8",
                                      "9","10","11","12",
                                      "13","14","15")),
         low_res = factor(low_res, levels = c("All.others", "Dinoflagellate", "Diatom")),
         frac = value/sum(value, na.rm = T)) %>% 
  select(Level7, Day, value, frac, Level8, low_res) %>%  
  group_by(Day, low_res) %>% 
  summarise(value = sum(value, na.rm = TRUE)) %>% 
  #filter(low_res == "Diatom") %>% 
  ggplot(., aes(x = as.numeric(Day), y = value, fill = low_res, group = low_res)) + 
  geom_area(stat = "identity", position = "fill", aes(fill = low_res), show.legend = F) +
  scale_fill_manual(breaks = c("All.others","Dinoflagellate", "Diatom"),
                    values = c(  "black", "light gray", "steel blue" )) +
  scale_x_continuous(breaks = c(1:15),
                     labels = c(1:15)) +
  labs(y = "Relative abundance",
       x = "") +
  theme_minimal() +
  theme(panel.grid.major.x = element_line(size = 0.4),
        panel.grid.major.y = element_line(size = 0.4),
        #axis.text.x = element_blank(),
        axis.title.y = element_text(size = 13))

#alternate view figure 1 - lines
Table_2018 %>%
  mutate(low_res = case_when(str_detect(Tax, regex("diatom", ignore_case = T)) ~"Diatom", 
                             str_detect(Tax, regex("dino", ignore_case = T)) ~"Dinoflagellate", 
                             !str_detect(Tax, regex("(diatom)|(dino)", ignore_case = T)) ~"All.others")) %>% 
  pivot_longer(cols = contains("Day"), names_to = "Day", values_to = "value") %>% 
  mutate(Day = str_replace_all(Day, "(Day)|(\\.2018)", ""),
         Day = factor(Day, levels = c("1","2","3","4",
                                      "5","6","7","8",
                                      "9","10","11","12",
                                      "13","14","15")),
         low_res = factor(low_res, levels = c("Dinoflagellate", "Diatom", "All.others")),
         frac = value/sum(value, na.rm = T)) %>% 
  select(Level7, Day, value, frac, Level8, low_res) %>% 
  group_by(low_res, Day) %>% 
  summarise(value = sum(value, na.rm = TRUE)) %>% 
  ggplot(., aes(x = Day, y = value, color = low_res, group = low_res)) + 
  geom_line() +
  scale_color_manual(breaks = c( "Dinoflagellate", "Diatom","All.others"),
                    values = c(  "#8c6510","#FF7399","#0D2B52" )) +
  labs(y = "Relative abundance",
       x = "Day") +
  theme_minimal() +
  theme(panel.grid.major.x = element_line(size = 0.4),
        panel.grid.major.y = element_line(size = 0.4))

# alternative view figure 2-3 - combine lines of diatoms

# Calculate diversity indices
## For diversity, I am going to aggregate at approxamately genus level (Level 7)
Sp_table = Table_2018 %>% 
  #filter(str_detect(Tax, regex("diatom", ignore_case = T))) %>% 
  pivot_longer(cols = contains("Day"), names_to = "Day", values_to = "value") %>% 
  select(Level7, Day, value) %>% 
  mutate(Day = str_replace_all(Day, "(Day)|(\\.2018)", ""),
         Day = factor(Day, levels = c("1","2","3","4",
                                      "5","6","7","8",
                                      "9","10","11","12",
                                      "13","14","15"))) %>%
  group_by(Level7, Day) %>% 
  summarise(count = sum(value, na.rm = T)) %>% 
  dcast(Level7~Day) %>% 
  column_to_rownames("Level7")

shannon<-diversity(Sp_table,index="shannon",2)
#Calculate inverse simpson index
invsimp<-diversity(Sp_table,index="invsimpson",2)
#Calculate the number of ASVs
ASV.sp_count <-colSums(Sp_table>0) #to evaluate species richness

## Some stats
stat_table_ASV = Table_2018 %>%
  mutate(low_res = case_when(str_detect(Tax, regex("diatom", ignore_case = T)) ~"Diatom", 
                             str_detect(Tax, regex("dino", ignore_case = T)) ~"Dinoflagellate", 
                             !str_detect(Tax, regex("(diatom)|(dino)", ignore_case = T)) ~"All.others")) %>% 
  pivot_longer(cols = contains("Day"), names_to = "Day", values_to = "value") %>% 
  mutate(Day = str_replace_all(Day, "(Day)|(\\.2018)", ""),
         Day = factor(Day, levels = c("1","2","3","4",
                                      "5","6","7","8",
                                      "9","10","11","12",
                                      "13","14","15")),
         low_res = factor(low_res, levels = c("Dinoflagellate", "Diatom", "All.others")),
         frac = value/sum(value, na.rm = T)) %>% 
  select(Level7, Day, value, frac, Level8, low_res) %>% 
  group_by(low_res) %>% 
  summarise(stdev = sd(value, na.rm = TRUE),
            var = var(value, na.rm = T),
            avg = mean(value, na.rm = TRUE),
            dif = max(value)- min(value))

#I imagine i will have to put these in a tidy dataframe to plot with the other methods.
data.frame()

## Ordination species level community
ASV_relabund = decostand(Table_2018[2:16], MARGIN = 2,  method = "total")
ASV_clust = hclust(dist(t(ASV_relabund)), method = "average")
plot(ASV_clust)

# Metatranscriptome sequencing: rRNA Fuhrman --------------------------------------------
load("~/Dropbox/SMP_R_working_directory/MetaTranscriptomics/Fuhrman_seqs/agg_rRNATax.rda")
#barplot
data.agg %>% 
  filter(!str_detect(Taxa, regex("(unknown)|(unas)", ignore_case = TRUE))) %>% 
  mutate(low_res = case_when(str_detect(Taxa, regex("diatom", ignore_case = T)) ~"Diatom", 
                             str_detect(Taxa, regex("dino", ignore_case = T)) ~"Dinoflagellate", 
                             !str_detect(Taxa, regex("(diatom)|(dino)", ignore_case = T)) ~"All.others"),
         low_res = factor(low_res, levels = c("Dinoflagellate", "Diatom", "All.others")),
         Day = str_replace_all(Day, "Day", ""),
         Day = factor(Day, levels = c("1","2","3","4",
                                      "5","6","7","8",
                                      "9","10","11","12",
                                      "13","14","15"))) %>% 
  ggplot(., aes(x = Day, y = x, fill = low_res)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(name = NULL,
                    values = c("#e31a1c", "#a1d99b", "#2d004b")) +
  labs(y = "Relative abundance",
       x = "Day") +
  theme_minimal() +
  theme(panel.grid.major.x = element_line(size = 0.4),
        panel.grid.major.y = element_line(size = 0.4))



tax_order=c("Ciliates","Alveolates-Other","Dinophyceae","Syndiniales","Amoebozoa","Excavates","Rhizaria-other","Radiolaria","Cercozoa","Foramnifera","Archaeplastids","Cryptophytes",    "Unclassified-Ochrophyte","Diatoms","Pelagophytes","Stramenopiles-Other","MAST",    "Hacrobia-other","Haptophytes","Telonemia","Choanoflagellida","Fungi","Opisthokonts-Other","Other/unknown")
tax_color = c("#7f0000","#b30000","#e31a1c","#ec7014",   "#8c510a","#dfc27d",    "#08306b","#08519c","#6baed6","#7fcdbb",     "#e7298a","#ffff99",      "#e5f5e0","#a1d99b","#41ab5d","#006d2c","#00441b",      "#8c510a","#f6e8c3","coral1",   "#d8daeb","#b2abd2","#542788","#2d004b")
names(tax_color)<-tax_order

# Recast as wide format in order to TMM normalize across samples.
data.agg_wide = data.agg %>% 
  dcast(Taxa~Day, value.var = "x") 
data.agg_wide[is.na(data.agg_wide)] = 0
rownames(data.agg_wide) = data.agg_wide$Taxa

#### Normalize using edgeR ###
ListDGE = DGEList(data.agg_wide[,2:15])
#ListDGE
ListDGE = calcNormFactors(ListDGE, method = "TMM") 
#ListDGE
TMMNorm_ASV_table = cpm(ListDGE)
#head(TMMNorm_ASV_table)
TMMNorm_ASV_table = as.data.frame(TMMNorm_ASV_table)
TMMNorm_ASV_table$Taxa = row.names(TMMNorm_ASV_table) #slap the ASV IDs back on the table for joining with the other sample info
write.csv(TMMNorm_ASV_table, file = "~/Dropbox/Project-data/SMPier/Transcriptomics/Figures/Method_comparison/Tables/Norm_rRNA.csv")

rRNA_for_plot = TMMNorm_ASV_table %>% 
  mutate("Day5" = NA) %>% 
  pivot_longer(cols = contains("Day")) %>% 
  rename(Day = name) %>% 
  mutate(low_res = case_when(str_detect(Taxa, regex("diatom", ignore_case = T)) ~"Diatom", 
                             str_detect(Taxa, regex("dino", ignore_case = T)) ~"Dinoflagellate", 
                             !str_detect(Taxa, regex("(diatom)|(dino)", ignore_case = T)) ~"All.others"),
         low_res = factor(low_res, levels = c("All.others","Dinoflagellate", "Diatom")),
         Day = str_replace_all(Day, "Day", ""),
         Day = factor(Day, levels = c("1","2","3","4",
                                      "5","6","7","8",
                                      "9","10","11","12",
                                      "13","14","15")),
         Taxa = factor(Taxa, levels = tax_order),
         material = "rRNA") %>% 
  filter(!str_detect(Taxa, regex("(unknown)|(unas)", ignore_case = TRUE)))

#alternative view - area plot
Fig1.rRNA = rRNA_for_plot %>% 
  group_by(Day, low_res) %>% 
  summarise(value = sum(value, na.rm = TRUE)) %>% 
  #filter(low_res == "Diatom") %>% 
  ggplot(., aes(y = value, fill = low_res, order = low_res)) +
  geom_area(na.rm = T, position = "fill", stat = "identity", aes(x = as.numeric(Day)), show.legend = F)+
  scale_x_continuous(breaks = c(1:15),
                     labels = c(1:15)) +
  scale_fill_manual(name = NULL,
                    values = c("black", "light gray", "steel blue")) +
  labs(y = "Relative abundance",
       x = "") +
  theme_minimal() +
  theme(panel.grid.major.x = element_line(size = 0.4),
        panel.grid.major.y = element_line(size = 0.4),
        #axis.text.x = element_blank(),
        axis.title.y = element_text(size = 13))

### Alternate view - diatoms
rRNA_diatoms = rRNA_for_plot %>% 
  #filter(low_res == "Diatom") %>% 
  ggplot(., aes(x = Day, y = value, fill = low_res)) +
  geom_bar(stat = "identity", position = "stack", aes(fill = low_res)) +
  scale_fill_manual(name = NULL,
                    values = c("#e31a1c", "#a1d99b", "#2d004b")) +
  labs(y = "Relative abundance",
       x = "Day") +
  theme_minimal() +
  theme(panel.grid.major.x = element_line(size = 0.4),
        panel.grid.major.y = element_line(size = 0.4))

ggsave(plot = last_plot(), filename = "~/Dropbox/Project-data/SMPier/Transcriptomics/Figures/Method_comparison/metat-rRNA_bars-SMP2018.pdf", width = 8, height = 6)

## Ordination for outliers
rRNA_relabund = vegan::decostand(TMMNorm_ASV_table[-15], MARGIN = 2, method = "total")
rRNAclust = hclust(dist(t(relabund)), method = "average")
plot(rRNAclust)

## Some stats
rRNA_for_plot %>% 
  group_by(low_res) %>% 
  summarise(stdev = sd(value, na.rm = T),
            var = var(value, na.rm = T),
            avg = mean(value, na.rm = T),
            diff = max(value, na.rm = T) - min(value, na.rm = T))

# Metatranscriptome sequencing: mRNA Fuhrman --------------------------------------------
load("~/Dropbox/SMP_R_working_directory/MetaTranscriptomics/Fuhrman_seqs/Normed_avg_annotated-Fuhrman.RData")
df_wtax = df_wtax %>% 
  select(-Day) %>% 
  mutate(Year = case_when(Year == "Apr_16_CMT" ~"Day1", 
                          Year == "Apr_17_CMT" ~"Day2", 
                          Year == "Apr_18_CMT" ~"Day3",
                          Year == "Apr_19_CMT" ~"Day4",
                          Year == "Apr_21_CMT" ~"Day6",
                          Year == "Apr_22_CMT" ~"Day7",
                          Year == "Apr_23_CMT" ~"Day8",
                          Year == "Apr_24_CMT" ~"Day9",
                          Year == "Apr_25_CMT" ~"Day10",
                          Year == "Apr_26_CMT" ~"Day11",
                          Year == "Apr_27_CMT" ~"Day12",
                          Year == "Apr_28_CMT" ~"Day13",
                          Year == "Apr_29_CMT" ~"Day14",
                          Year == "Apr_30_CMT" ~"Day15")) %>% 
  rename("Day" = Year)


df_wtax_noNA = df_wtax %>% 
  filter(!str_detect(Taxonomy, regex("Not assigned", ignore_case = TRUE)))

## Re-compile taxonomy breakdown for visualization purposes
compile_tax<-function(df){
  df$Nextlevel<-df$Phylum
  df$Nextlevel[df$Phylum==""]="Other"
  df$Nextlevel[df$Class=="Foraminifera"]="Foraminifera"
  df$Nextlevel[df$Class=="Acantharia"]="Acantharia"
  df$Nextlevel[df$Class=="Polycystinea"]="Polycystinea" 
  df$Nextlevel[df$Class=="Syndinians"]="Syndiniales" 
  df$Nextlevel[df$Class=="Bacillariophyceae"]="Diatom" 
  df$Nextlevel[df$Class=="Pelagophyceae"]="Pelagophyceae" 
  df$tax_compiled<-df$Supergroup
  df$tax_compiled[df$Supergroup == "Alveolate"]<-"Other Alveolate"
  df$tax_compiled[df$Nextlevel=="Ciliate"]="Ciliate"
  df$tax_compiled[df$Nextlevel=="Dinoflagellate"]="Dinoflagellate"
  df$tax_compiled[df$Nextlevel=="Syndiniales"]="Syndiniales"
  #
  df$tax_compiled[df$Supergroup == "Archaeplastida"]<-"Other Archaeplastida"
  df$tax_compiled[df$Nextlevel=="Chlorophyta"]="Chlorophyta"
  #
  df$tax_compiled[df$Supergroup=="Rhizaria"]<-"Other Rhizaria"
  #df$tax_compiled[df$Nextlevel=="Cercozoa"]="Cercozoa"
  #df$tax_compiled[df$Nextlevel=="Retaria"]="Retaria"
  #df$tax_compiled[df$Nextlevel=="Acantharia"]="Acantharia"
  #df$tax_compiled[df$Nextlevel=="Foraminifera"]="Foraminifera"
  #df$tax_compiled[df$Nextlevel=="Polycystinea"]="Polycystinea"
  #
  df$tax_compiled[df$Supergroup=="Stramenopile"]<-"Other Stramenopile"
  df$tax_compiled[df$Nextlevel=="MAST"]="MAST"
  df$tax_compiled[df$Nextlevel=="Diatom"]="Diatom"
  df$tax_compiled[df$Nextlevel=="Pelagophyceae"]="Pelagophytes"
  #
  other<-c("Amoebozoa", "Cryptista", "Discoba")
  df$tax_compiled[df$Supergroup %in% other]="Other"
  return(df)
}
# add on more levels of taxonomy
df_tax<-compile_tax(df_wtax_noNA)
mRNA_for_plot = df_tax %>% 
  mutate(Phylum = ifelse(Phylum == "", "XXX", Phylum),
         Class = ifelse(Class == "", "XXX", Class),
         Order = ifelse(Order == "", "XXX", Order),
         Family = ifelse(Family == "", "XXX", Family),
         Genus = ifelse(Genus == "", "XXX", Genus),
         Species = ifelse(Species == "", "XXX", Species)) %>% 
  filter(Genus != "XXX") %>% 
  select(Nextlevel, Genus, Day, mean_TPM, tax_compiled) %>% 
  mutate(low_res = case_when(str_detect(tax_compiled, regex("diatom", ignore_case = T)) ~"Diatom", 
                             str_detect(tax_compiled, regex("dino", ignore_case = T)) ~"Dinoflagellate", 
                             !str_detect(tax_compiled, regex("(diatom)|(dino)", ignore_case = T)) ~"All.others"),
         Day = str_replace_all(Day, "Day", ""),
         Day = factor(Day, levels = c("1","2","3","4",
                                      "5","6","7","8",
                                      "9","10","11","12",
                                      "13","14","15")),
         low_res = factor(low_res, levels = c("All.others", "Dinoflagellate", "Diatom")),
         material = "mRNA")
write.csv(mRNA_for_plot, file = "~/Dropbox/Project-data/SMPier/Transcriptomics/Figures/Method_comparison/Tables/Norm_mRNA.csv")

# alternative view - area plot
Fig1.mRNA = mRNA_for_plot %>% 
  group_by(low_res, Day) %>% 
  summarise(TPM = sum(mean_TPM, na.rm = T)) %>% 
  dcast(low_res~Day) %>% mutate(Day5 = NA) %>% 
  pivot_longer(cols = -low_res, names_to = "Day", values_to = "TPM") %>% 
  mutate(Day = str_replace_all(Day, "Day", ""),
         Day = factor(Day, levels = c("1","2","3","4",
                                      "5","6","7","8",
                                      "9","10","11","12",
                                      "13","14","15"))) %>% 
  #filter(low_res == "Diatom") %>% 
  ggplot(., aes(y = TPM, fill = low_res, order = low_res)) +
  #geom_bar(stat = "identity", position = "stack") +
  #geom_line(aes(group = 1)) + geom_point() +
  geom_area(position = "fill", stat = "identity", aes(x = as.numeric(Day)), show.legend = F)+
  scale_x_continuous(breaks = c(1:15),
                     labels = c(1:15)) +
  scale_fill_manual(name = NULL,
                    values = c("black", "light gray", "steel blue")) +
  labs(y = "Relative abundance",
       x = "Day") +
  theme_minimal() +
  theme(panel.grid.major.x = element_line(size = 0.4),
        panel.grid.major.y = element_line(size = 0.4))

# alternate view - figure 1
 mRNA_for_plot %>% 
   group_by(low_res, Day) %>% 
   summarise(TPM = sum(mean_TPM, na.rm = T)) %>% 
  #filter(low_res == "Diatom") %>% 
  ggplot(., aes(x = Day, y = TPM, fill = low_res)) +
  geom_bar(stat = "identity", position = "stack", color = "white") +
  scale_fill_manual(name = NULL,
                    values = c("light gray", "black", "steel blue")) +
  labs(y = "Relative abundance",
       x = "Day") +
  theme_minimal() +
  theme(panel.grid.major.x = element_line(size = 0.4),
        panel.grid.major.y = element_line(size = 0.4))

ggsave(plot = last_plot(), filename = "~/Dropbox/Project-data/SMPier/Transcriptomics/Figures/Method_comparison/metat-mRNA_bars-SMP2018.pdf", width = 8, height = 6)

Fig1.enum + Fig1.ASV + Fig1.rRNA + Fig1.mRNA + patchwork::plot_layout(ncol = 1)

ggsave(plot = last_plot(), filename = "~/Dropbox/Project-data/SMPier/Transcriptomics/Figures/Method_comparison/Fig1-CombinedDynamics.pdf", width = 8, height = 6)


# Differential analysis of diatom reads from mRNA -------------------------

#Load in the count table
Diatom_raw_df = read.delim("~/Dropbox/SMP_R_working_directory/MetaTranscriptomics/2018_assembly/raw-count-data-metaT-SMP2018.txt") %>%
  filter(str_detect(Taxonomy, regex("bacillariophyceae", ignore_case = T))) %>%  #filter for diatom reads using string detect
  select(-Taxonomy)
view(Diatom_raw_df %>% arrange(KO))
#Note that we can get the whole community by removing the filter statement.
#Count the number of diatom reads in the filtered dataset
summary(colSums(Diatom_raw_df[-1]))
#And count the number of KO terms. I plan to use the KO IDs as the rownames;
#therefore I'll will aggregate them for each sample. Flattening the tax diversity to FX diversity.
Diatom_raw_df %>%
  distinct(KO) %>% nrow()
# I can expect a max o ~5622
#aggregate
agg_df = Diatom_raw_df %>%
  pivot_longer(cols = starts_with("SMP"), names_to = "sample", values_to = "value") %>%
  mutate(sample = str_replace(sample, "SMPier\\.", ""),
         sample = str_replace(sample, regex("\\.\\d{2}\\.\\d{2}\\.2018_S\\d+"), ""),
         sample = factor(sample, levels = c("Day1.A", "Day1.B", "Day1.C", 
                                            "Day3.A", "Day3.B", "Day3.C", 
                                            "Day5.A", "Day5.B", "Day5.C", 
                                            "Day9.A", "Day9.B", "Day9.C", 
                                            "Day11.A", "Day11.B", "Day11.C"))) %>% 
  group_by(KO, sample) %>%
  summarise(count = round(sum(value, na.rm = TRUE))) %>% 
  dcast(KO~sample)
rownames(agg_df) = agg_df$KO

# At this point i have a count table containing only diatom reads where KO counts have been aggregated using a sum. Maybe this shoudl be the mean.
agg_df[1:10, 1:10]
# Normalize data with edgeR
# create the dge object
sample_list = c("Day1","Day3","Day5","Day9","Day11")
#metadata = factor(c("day1","day1","day1","day3","day3","day3","day5","day5","day5","day9","day9","day9","day11","day11","day11"), levels = c("day1","day3","day5","day9","day11"))
metadata = factor(c(rep("day1",3),
                    rep("day3",3),
                    rep("day5",3),
                    rep("day9",3),
                    rep("day11",3)),
                  levels = str_to_lower(sample_list)) # this is important for purposefully setting the sample order.
# make a vector for the sample list


# Use this information to create the DGEList object
dge_obj = DGEList(counts = as.data.frame.matrix(agg_df[-1]), # Counts is the actual columns of the dataframe
                  genes = agg_df[1], # This is the name of the column(similar to genes)
                  group = metadata)

# Filter out low count rows using a conditional filter--that is, to keep rows that meet a certian condition.
keep = filterByExpr(dge_obj)
dge_obj = dge_obj[keep, ,keep.lib.sizes=FALSE]
# Calculate the normalizing factors using the TMM process (actually optional)
dge_obj = calcNormFactors(dge_obj, method = "TMM")
# setup the design matrix without an intercept as day 1!
design = model.matrix(~0+group, data = dge_obj$samples)

# Below is the design with day 1 as an intercept (incorrect for what I'm doing)
#design_w.intercept = model.matrix(~group, data = dge_obj1$samples)

#Set the column names for the design to match the
colnames(design) = levels(dge_obj$samples$group)



# using makeContrasts() to specify the pariwise comparisons
# comparing day 12 to all days prior (1 & 4); day12 to day18; day18 to day20.
#conts = makeContrasts(day1+day4-day12, day12-day18, day18-day20, levels = str_to_lower(sample_list)) #for dinoflagellate bloom
#I think the below follows the correct convention; above is incorrect
#conts = makeContrasts(day3-day1, day5-day3, day9-day5, day11-day9, levels = str_to_lower(sample_list)) #for pairwise series contrasts
conts = makeContrasts(day3-day1, day5-day1, day9-day1, day11-day1, levels = str_to_lower(sample_list)) #for contrasting all samples to prebloom state

# generate the design matrix
#design = model.matrix(~metadata)
#olnames(design) = str_replace_all(colnames(design), regex("metadata", ignore_case = T), "")

# estimate common and tagwise(pairwise) dispersion accrding the study design matrix
disp = estimateGLMCommonDisp(dge_obj, design = design)
disp = estimateGLMTagwiseDisp(disp, design = design)

# Determine differentially expressed genes
fit = glmQLFit(disp, design = design)
DEGs = glmQLFTest(fit, contrast = conts)
DEG_table = DEGs$table %>% data.frame()

# adjust p values due to so many comparissons to decrease chance of type 1 error.
DEG_table$P.adjusted = p.adjust(DEG_table$PValue, method = "fdr")
DEG_table$KO = rownames(DEG_table)
DEG_table = DEG_table %>% filter(P.adjusted < 0.01)


#Before moving on, what is the total number of differentially expressed genes were there on the different Days?
DEG_table %>% nrow

# The point of this analysis is to investigate diatom gene expression during this bloom so lets filter out the diatom reads


Diatom_Degs_wcat = DEG_table %>% 
  #separate(ID, into = c("tax","KO"), sep="_") %>%
  left_join(., read.csv("~/Dropbox/SMP_R_working_directory/MetaTranscriptomics/K0_geneIDs_fulllist_08192019.csv", header = T)) %>%
  left_join(., read.delim("~/Dropbox/SMP_R_working_directory/MetaTranscriptomics/Custom_KO_list_30-04-2020.txt", header = TRUE, sep = "\t"))  

save(DDegs_wcat, file = "~/Dropbox/SMP_R_working_directory/MetaTranscriptomics/2018_assembly/DEG_table-assembly2018.rda")
load("~/Dropbox/SMP_R_working_directory/MetaTranscriptomics/2018_assembly/DEG_table-assembly2018.rda")

# Arrange data for plotting.
##use the below for the full contrast list
# Or everything contrasted against day 1
##use the below for the full contrast list
DEG_table = DEG_table %>% 
  dplyr::rename("Day3" = logFC.day3...day1,
                "Day5" = logFC.day5...day1,
                "Day9" = logFC.day9...day1,
                "Day11" = logFC.day11...day1) %>%
  pivot_longer(cols = c("Day3", "Day5", "Day9", "Day11"), names_to = "contrasts", values_to = "fold_chng") %>%
  select(-logCPM, -F, -PValue) %>%
  select(KO, contrasts, fold_chng, everything())
#how many different KOs are differentially expressed
Diatom_Degs_wcat %>% filter(contrasts == "Day3") %>%
  count(KO) %>% arrange(desc(n)) %>% pull(n) %>% sum()

#Diatom_Degs_wcat %>% filter(contrasts == "Day3") %>% count(KO) %>% arrange(desc(n)) %>% pull(n) %>% sum() == Diatom_Degs_wcat %>% filter(contrasts == "Day5") %>% count(KO) %>% arrange(desc(n)) %>% pull(n) %>% sum() == Diatom_Degs_wcat %>% filter(contrasts == "Day9") %>% count(KO) %>% arrange(desc(n)) %>% pull(n) %>% sum() == Diatom_Degs_wcat %>% filter(contrasts == "Day11") %>% count(KO) %>% arrange(desc(n)) %>% pull(n) %>% sum()
# 1,565 of them are differentially expressed

## how many differentially expressed genes are there with an logfc greatr than |1|?
# i didn't select a higher logFC for two reasons:
# -We dont know how much of a logFC is enough to produce a significant cellular reaction.
# -Better chance of finding genes once cross referenced with the WGCNA modules.

#How many DEGs with logFCs greater than |1| on each day?
DEG_table$direction = "minor"
DEG_table$direction[DEG_table$fold_chng > 1] = "upreg"
DEG_table$direction[DEG_table$fold_chng < -1] = "downreg"

#major_FCs is a table displaying the number and fraction of up and down > and < |1| 
major_FCs = Diatom_Degs_wcat %>%
  group_by(contrasts) %>%
  count(direction) %>%
  pivot_wider(names_from = direction, values_from = n) %>%
  mutate(totalDE = sum(downreg,upreg, na.rm = TRUE),
         fracDE = totalDE / sum(downreg, minor, upreg, na.rm = T),
         fracUp = upreg / sum(downreg, minor, upreg, na.rm = T),
         fracDn = downreg / sum(downreg, minor, upreg, na.rm = T),
         fracM = minor / sum(downreg, minor, upreg, na.rm = T)) %>% ungroup()
write.csv(major_FCs, file = "~/Dropbox/SMP_R_working_directory/MetaTranscriptomics/2018_assembly/logFC_table.csv", row.names = F)
#use permutations of the below to ask specific questions about the proportion of up down and muted DEGs
major_FCs %>% filter(contrasts != "Day9") %>% select(fracDn) %>% pull() %>% mean(., na.rm = T)

# plot this information
DEG_table %>%
  #filter(KO %in% D5_subexpr$KO) %>%
  mutate(hilight = abs(fold_chng) > abs(1),
         contrasts = factor(contrasts, levels = c("Day3","Day5","Day9","Day11"))) %>%
  #mutate(hilight = ifelse(abs(fold_chng) > abs(7), "super", hilight)) %>%
  filter(hilight) %>% 
  ggplot(., aes(x = contrasts, y = fold_chng)) +
  #geom_jitter(width = .2, size = .8, shape = 21, color = "black", alpha = 0.9, show.legend = F) +
  geom_boxplot(alpha = .4, outlier.size = .9, lwd = .4,  show.legend = T, outlier.colour = "white", outlier.shape = "triangle") +
  stat_summary(fun = "mean", shape = 23, size = .4, fill = "black") +
  labs(y = "Log Fold Change",
       x = "") +
  #scale_fill_manual(values =  jisho_picker("american"), breaks = c("FALSE","TRUE"), labels = c("< logFC |1|", "> logFC |1|"), name = "") +
  #scale_color_manual(values = jisho_picker("american")) +
  theme_minimal()
#or as a boxplot

ggsave(plot = last_plot(), filename = "~/Dropbox/Project-data/SMPier/Transcriptomics/Figures/Method_comparison/DEGS-figure5_wjitter.pdf", width = 8, height = 6)
ggsave(plot = last_plot(), filename = "~/Dropbox/Project-data/SMPier/Transcriptomics/Figures/Method_comparison/DEGS-figure5_filtdboxplotWmean.pdf", width = 8, height = 6)

# Correlations of abundances: dot plot ------------------------------------
#ASV table treatment
ASV_table = Table_2018 %>% rownames_to_column("toss") %>% select(-toss) %>% 
  mutate(Group = case_when(str_detect(Level4, regex("dino", ignore_case = T)) ~ "Dinoflagellate", 
                           str_detect(Level4, regex("bacillariophyta", ignore_case = T)) ~ "Diatom",
                           !str_detect(Level4, regex("(dino)|(bacillario)", ignore_case = TRUE)) ~ "All.others")) %>% 
  select(starts_with("Day"), Level7, Group)


ASV_FracTable = ASV_table %>% 
  select(starts_with("Day")) %>% 
  decostand(., MARGIN = 2, method = "total") %>% signif(., 2) %>% 
  mutate(ID = rownames(.)) %>% 
  left_join(., ASV_table %>% mutate(ID = rownames(.)) %>% select(Level7, Group, ID), by = "ID") 

ASV_FracTable_long = ASV_FracTable %>% 
  pivot_longer(cols = starts_with("Day"), names_to = "Day", values_to = "count") %>% 
  mutate(Day = str_replace(Day, "\\.2018","")) %>% 
  group_by(Group, Day) %>% 
  summarise(ASV.sum_frac = sum(count, na.rm = T)) %>% 
  unite(ID, Group, Day, sep = "-") 
# I think here If i want to use the mean and sd, I'll need to first aggregate at the Level7(genus). 
# For now I'll just just teh sum to represent the fraction of diatoms, dinos, and others.

#Enumeration table treatment
enum_FracTable = count_table %>% 
  select(starts_with("Day")) %>% 
  decostand(., MARGIN = 2, method = "total") %>% signif(., 2) %>% 
  mutate(ID = rownames(.)) %>% #so that i can kdeep things straight and rejoin with the names
  left_join(., count_table %>% mutate(ID = rownames(.)) %>% select(Group, Genus, ID), by = "ID")

enum_FracTable_long = enum_FracTable %>% 
  pivot_longer(cols = starts_with("Day"), names_to = "Day", values_to = "count") %>% 
  group_by(Group, Day) %>% 
  summarise(E.sum_frac = sum(count, na.rm = T)) %>% 
  unite(ID, Group, Day, sep = "-") 
#This table will be left joined with the others using the ID column.
enum_FracTable_long %>% filter(str_detect(ID, "Day14")) #sanity check

#rRNA table treatment
data.agg_wide = data.agg %>% 
  mutate(Day = factor(Day, levels = c("Day1", "Day2", "Day3", "Day4", "Day6", "Day7", "Day8", "Day9", "Day10", "Day11", "Day12", "Day13", "Day14", "Day15"))) %>% 
  dcast(Taxa~Day) 
  

data.agg_long = data.agg_wide %>% 
  select(-Taxa) %>% 
  vegan::decostand(., MARGIN = 2, method = "total", na.rm = T) %>% 
  mutate(ID = rownames(.)) %>% 
  left_join(., data.agg_wide %>% mutate(ID = rownames(.)) %>% select(Taxa, ID), by = "ID") %>% 
  mutate(Group = case_when(str_detect(Taxa, regex("dino", ignore_case = T)) ~ "Dinoflagellate",
                           str_detect(Taxa, regex("Diatom", ignore_case = T)) ~ "Diatom",
                           !str_detect(Taxa, regex("(dino)|(bacillario)", ignore_case = TRUE)) ~ "All.others")) %>% 
  pivot_longer(cols = starts_with("Day"), names_to = "Day", values_to = "count") %>% 
  filter(!str_detect(Taxa, regex("(unknown)|(unas)", ignore_case = TRUE))) %>% 
  group_by(Group, Day) %>% 
  summarise(rRNA.sum_frac = sum(count, na.rm = T)) %>% 
  unite(col = ID, Group, Day, sep = "-")


TMMNorm_ASV_table %>% 
  #mutate("Day5" = NA) %>% 
  pivot_longer(cols = contains("Day"), names_to = "Day", values_to = "value") %>% 
  mutate(low_res = case_when(str_detect(Taxa, regex("diatom", ignore_case = T)) ~"Diatom", 
                             str_detect(Taxa, regex("dino", ignore_case = T)) ~"Dinoflagellate", 
                             !str_detect(Taxa, regex("(diatom)|(dino)", ignore_case = T)) ~"All.others"),
         low_res = factor(low_res, levels = c("All.others","Dinoflagellate", "Diatom")),
         Day = factor(Day, levels = c("Day1", "Day2", "Day3", "Day4", 
                                      "Day6", "Day7", "Day8", "Day9", "Day10",
                                      "Day11", "Day12", "Day13", "Day14", "Day15"))) %>% 
  filter(!str_detect(Taxa, regex("(unknown)|(unas)", ignore_case = TRUE))) %>% 
  mutate(frac = value/sum(value, na.rm = TRUE)) %>% 
  group_by(low_res, Day) %>% 
  summarise(rRNA.sum_frac = sum(frac, na.rm = TRUE)) %>% filter(Day == "Day15") %>% pull(rRNA.sum_frac) %>% sum()

rRNA_long = TMMNorm_ASV_table %>% 
  #mutate("Day5" = NA) %>% 
  pivot_longer(cols = contains("Day"), names_to = "Day", values_to = "value") %>% 
  mutate(low_res = case_when(str_detect(Taxa, regex("diatom", ignore_case = T)) ~"Diatom", 
                             str_detect(Taxa, regex("dino", ignore_case = T)) ~"Dinoflagellate", 
                             !str_detect(Taxa, regex("(diatom)|(dino)", ignore_case = T)) ~"All.others"),
         low_res = factor(low_res, levels = c("All.others","Dinoflagellate", "Diatom")),
         Day = factor(Day, levels = c("Day1", "Day2", "Day3", "Day4", 
                                      "Day6", "Day7", "Day8", "Day9", "Day10",
                                      "Day11", "Day12", "Day13", "Day14", "Day15"))) %>% 
  filter(!str_detect(Taxa, regex("(unknown)|(unas)", ignore_case = TRUE))) %>% 
  group_by(Day, low_res) %>% 
  summarise(agg = sum(value)) %>% ungroup %>% 
  select(low_res,Day, agg) %>%  # filter(low_res == "All.others") %>% head()
  dcast(low_res~Day) %>% 
  column_to_rownames("low_res") %>% 
  vegan::decostand(., MARGIN = 2, method = "total") %>% 
  rownames_to_column("Group") %>% 
  pivot_longer(cols = contains("Day"), names_to = "Day", values_to = "rRNA.sum_frac") %>% 
  unite(ID, Group, Day, sep = "-")
  
   
  mutate(frac = value/sum(value, na.rm = TRUE)) %>% 
  group_by(low_res, Day) %>% 
  summarise(rRNA.sum_frac = sum(frac, na.rm = TRUE)) %>% filter(Day == "Day15") %>% pull(rRNA.sum_frac) %>% sum()

#mRNA table treatment (undone)
mRNA_FracTable_long = df_tax %>% 
    mutate(Phylum = ifelse(Phylum == "", "XXX", Phylum),
           Class = ifelse(Class == "", "XXX", Class),
           Order = ifelse(Order == "", "XXX", Order),
           Family = ifelse(Family == "", "XXX", Family),
           Genus = ifelse(Genus == "", "XXX", Genus),
           Species = ifelse(Species == "", "XXX", Species)) %>% 
    filter(Genus != "XXX") %>% 
    select(Nextlevel, Genus, Day, mean_TPM, tax_compiled) %>% 
    mutate(low_res = case_when(str_detect(tax_compiled, regex("diatom", ignore_case = T)) ~"Diatom", 
                               str_detect(tax_compiled, regex("dino", ignore_case = T)) ~"Dinoflagellate", 
                               !str_detect(tax_compiled, regex("(diatom)|(dino)", ignore_case = T)) ~"All.others"),
           Day = factor(Day, levels = c("Day1","Day2","Day3","Day4","Day6",
                                        "Day7","Day8","Day9","Day10",
                                        "Day11","Day12","Day13","Day14","Day15"))) %>% 
  group_by(low_res, Day) %>% 
  summarise(count = sum(mean_TPM)) %>% 
  ungroup() %>% 
  dcast(low_res~Day) %>% 
  column_to_rownames("low_res") %>%
  vegan::decostand(MARGIN = 2, method = "total") %>% 
  rownames_to_column("tax") %>% 
  pivot_longer(cols = starts_with("Day"), names_to = "Day", values_to = "mRNA.sum_frac") %>% 
  unite(ID, tax, Day, sep = "-")
    
           

##Combine some tables
for_plot = rRNA_long %>% 
  left_join(., enum_FracTable_long) %>% 
  left_join(., ASV_FracTable_long) %>% 
  left_join(., mRNA_FracTable_long) %>% 
  separate(ID, into = c("Tax", "Day"), sep = "-", remove = F)
  #filter(str_detect(ID, "Day15$"))
  filter(str_detect(ID, regex("diatom", ignore_case = T))) 

#Plot Enumeration vs rRNA
E_vs_rRNA.plot = for_plot %>% 
  #filter(str_detect(ID, "Day10")) %>% 
  ggplot(., mapping = aes(x = rRNA.sum_frac, y = E.sum_frac, fill = Tax)) +
  geom_abline(alpha = 0.3, linetype = "dashed") +
  #geom_smooth(method = "lm",se = F, size = .2, colour = "red", linetype = "dashed", alpha = .5) +
  geom_point(alpha = 0.8, size = 2, color = "black", aes(shape = Tax), show.legend = F) +
  coord_cartesian(xlim = c(0,.75),
                  ylim = c(0,1)) +
  scale_x_continuous(breaks = c(seq(0,1,.1))) +
  scale_y_continuous(breaks = c(seq(0,1,.1))) +
  scale_fill_manual(name = NULL,
                    breaks = c("All.others","Dinoflagellate", "Diatom"),
                    values = c("black", "light grey", "steel blue")) +
  scale_shape_manual(name = NULL,
                     breaks = c("All.others","Dinoflagellate", "Diatom"),
                     values = c(21, 22, 24)) +
  labs(y = "Proportion (microscopy)",
       x = "Proportion (rRNA)") +
  theme_minimal() +
  theme(panel.grid.minor = element_line(linetype = "dashed"),
        panel.grid.major = element_line(size = .3))

#Enumeration vs ASV
E_vs_ASV.plot = for_plot %>% 
  #filter(str_detect(ID, "Day10")) %>% 
  ggplot(., mapping = aes(x = ASV.sum_frac, y = E.sum_frac, fill = Tax)) +
  geom_abline(alpha = 0.3, linetype = "dashed") +
  #geom_smooth(method = "lm",se = F, size = .2, colour = "red", linetype = "dashed", alpha = .5) +
  geom_point(alpha = 0.8, size = 2, color = "black", aes(shape = Tax), show.legend = F) +
  coord_cartesian(xlim = c(0,.75),
                  ylim = c(0,1)) +
  scale_x_continuous(breaks = c(seq(0,1,.1))) +
  scale_y_continuous(breaks = c(seq(0,1,.1))) +
  scale_fill_manual(name = NULL,
                     breaks = c("All.others","Dinoflagellate", "Diatom"),
                     values = c("black", "light grey", "steel blue")) +
  scale_shape_manual(name = NULL,
                     breaks = c("All.others","Dinoflagellate", "Diatom"),
                     values = c(21, 22, 24)) +
  labs(y = "Proportion (microscopy)",
       x = "Proportion (ASV)") +
  theme_minimal() +
  theme(panel.grid.minor = element_line(linetype = "dashed"),
        panel.grid.major = element_line(size = .3))

#ASV vs rRNA
rRNA_v_ASV.plot = for_plot %>% 
  #filter(str_detect(ID, "Day10")) %>% 
  ggplot(., mapping = aes(y = rRNA.sum_frac, x = ASV.sum_frac, fill = Tax)) +
  geom_abline(alpha = 0.3, linetype = "dashed") +
  #geom_smooth(method = "lm",se = F, size = .2, colour = "red", linetype = "dashed", alpha = .5) +
  geom_point(alpha = 0.8, size = 2, show.legend = F, color = "black", aes(shape = Tax)) +
  coord_cartesian(xlim = c(0,.75),
                  ylim = c(0,.75)) +
  scale_x_continuous(breaks = c(seq(0,1,.1))) +
  scale_y_continuous(breaks = c(seq(0,1,.1))) +
  scale_fill_manual(name = NULL,
                     breaks = c("All.others","Dinoflagellate", "Diatom"),
                     values = color_pal) +
  scale_shape_manual(name = NULL,
                     breaks = c("All.others","Dinoflagellate", "Diatom"),
                     values = c(21, 22, 24)) +
  labs(y = "Proportion (rRNA)",
       x = "Proportion (ASV)") +
  theme_minimal() +
  theme(panel.grid.minor = element_line(linetype = "dashed"),
        panel.grid.major = element_line(size = .3))

#rRNA vs mRNA plot
rRNA_v_mRNA.plot = for_plot %>% 
  #filter(str_detect(ID, "Day10")) %>% 
  ggplot(., mapping = aes(y = rRNA.sum_frac, x = mRNA.sum_frac, fill = Tax)) +
  #geom_smooth(method = "lm",se = F, size = .2, colour = "red", linetype = "dashed", alpha = .5) +
  geom_abline(alpha = 0.3, linetype = "dashed") +
  #geom_smooth(method = "lm",se = F, size = .2, colour = "red", linetype = "dashed", alpha = .5) +
  geom_point(size = 2, show.legend = T, color = "black", aes(shape = Tax)) +
  coord_cartesian(xlim = c(0,.75),
                  ylim = c(0,.75)) +
  scale_x_continuous(breaks = c(seq(0,1,.1))) +
  scale_y_continuous(breaks = c(seq(0,1,.1))) +
  scale_fill_manual(name = NULL,
                    breaks = c("All.others","Dinoflagellate", "Diatom"),
                    values = color_pal) +
  scale_shape_manual(name = NULL,
                     breaks = c("All.others","Dinoflagellate", "Diatom"),
                     values = c(21, 22, 24)) +
  labs(y = "Proportion (rRNA)",
       x = "Proportion (mRNA)") +
  theme_minimal() +
  theme(panel.grid.minor = element_line(linetype = "dashed"),
        panel.grid.major = element_line(size = .3))

# stitch
E_vs_rRNA.plot + E_vs_ASV.plot + rRNA_v_ASV.plot + rRNA_v_mRNA.plot + patchwork::plot_layout(ncol = 2)

ggsave(plot = last_plot(), filename = "~/Dropbox/Project-data/SMPier/Transcriptomics/Figures/Method_comparison/correlativeDots_alltax-SMP2018.pdf", width = 10, height = 8)

one55= c("#00D973", "#FA2B00","#000831")
american = c("darkgrey","#A90636", "steelblue")


# correlations ------------------------------------------------------------
corrs = cor(for_plot[4:7])
cor.matrix(for_plot[4], for_plot[5], method = "pearson")
corrplot::corrplot(corrs, method = "shade" )


WGCNA::cor(for_plot[4:7], use = "p") %>% #pairwise.complete.obs
#calculate p values for the correlations.
WGCNA::corPvalueStudent(., 9) 


## are the differences significant? Anova vs Kruskal
#1. are the values normal
shapiro.test(for_plot$E.sum_frac)
hist(for_plot$E.sum_frac)
#not normally distributed: go with kruskal
penguins

for_plot %>% 
  pivot_longer(cols = c(rRNA.sum_frac, E.sum_frac, ASV.sum_frac, mRNA.sum_frac), names_to = "material", values_to = "values") %>% 
  #pull(values) %>% shapiro.test() #data are normally distributed.
  ggplot(aes(sample = values, group = material)) + geom_qq() + geom_qq_line() + facet_wrap(~material)
# data is not normally distributed go with the non-parametric (kruskal) test of variation between and within samples
  #kruskal.test(values ~ material, data = .) %>% 
  #aov(values ~ material, data = .) %>% the data not normally distributed
#  summary()


# Ordination -------------------------------------------------------------------
#rRNA samples


# format wide table for cluster dendrogram
wide_table = for_plot %>% 
  pivot_longer(cols = c(rRNA.sum_frac, E.sum_frac, ASV.sum_frac, mRNA.sum_frac), names_to = "material", values_to = "values") %>% 
  unite(Sample, Day, material, sep = "-") %>% 
  dcast(Tax~Sample)

#calculate distances and clusters
clust = hclust(dist(t(wide_table[-1])), method = "average")
#plot dendrogram
plot(clust)

## NMDS
NMDS.Bray = metaMDS(t(wide_table[-1]), distance = "bray")
NMDS.bray_pts<-NMDS.Bray$points[1:nrow(NMDS.Bray$points),] %>% data.frame()
NMDS.bray_pts$Sample<-row.names(NMDS.bray_pts)
NMDS.bray_pts$Sample=factor(NMDS.bray_pts$Sample, levels = NMDS.bray_pts$Sample)
NMDS.bray_pts %>% 
  mutate(Material = case_when(str_detect(Sample, "ASV") ~ "ASV",
                              str_detect(Sample, "E") ~ "Micro",
                              str_detect(Sample, "mRNA") ~ "mRNA",
                              str_detect(Sample, "rRNA") ~ "rRNA")) %>% 
  ggplot(., aes(x = MDS1, y = MDS2, fill = Material)) +
  geom_point(size = 3, alpha = 0.8, aes(shape = Material)) +
  scale_shape_manual(name = "",
                     values = c(21, 22, 23, 24)) +
  scale_fill_manual(name = "",
                    values = c("light grey","#A90636","steel blue","black")) +
  theme_minimal()
  
ggsave(plot = last_plot(), filename = "~/Dropbox/Project-data/SMPier/Transcriptomics/Figures/Method_comparison/NMDS_materialsSMP2019.pdf", width = 6, height = 5) 
  
###### Working on getting some basic stats from these distributions to address the question of quantititatively demonstrating that there is a difference between the distributions presented by the differenyt methods.
# I think the question is a matter of differences in standard deviation / variance of the proportions. 
# The first step that I need to correct in the below is to convert these abundances from reads to percentages.

### Stats-microscopy
StatsTable = Tidy_df %>% 
  group_by(Day, Group) %>% 
  summarise(value = sum(value, na.rm = TRUE)) %>% 
  group_by(Group) %>% 
  summarise(stdev = sd(value),
            var = var(value),
            avg = mean(value),
            dif = max(value)-min(value))

## Some stats - ASV
stat_table_ASV = Table_2018 %>%
  mutate(low_res = case_when(str_detect(Tax, regex("diatom", ignore_case = T)) ~"Diatom", 
                             str_detect(Tax, regex("dino", ignore_case = T)) ~"Dinoflagellate", 
                             !str_detect(Tax, regex("(diatom)|(dino)", ignore_case = T)) ~"All.others")) %>% 
  pivot_longer(cols = contains("Day"), names_to = "Day", values_to = "value") %>% 
  mutate(Day = str_replace_all(Day, "(Day)|(\\.2018)", ""),
         Day = factor(Day, levels = c("1","2","3","4",
                                      "5","6","7","8",
                                      "9","10","11","12",
                                      "13","14","15")),
         low_res = factor(low_res, levels = c("Dinoflagellate", "Diatom", "All.others")),
         frac = value/sum(value, na.rm = T)) %>% 
  select(Level7, Day, value, frac, Level8, low_res) %>% 
  group_by(low_res) %>% 
  summarise(stdev = sd(value, na.rm = TRUE),
            var = var(value, na.rm = T),
            avg = mean(value, na.rm = TRUE),
            dif = max(value)- min(value))
###
## Some stats-rRNA
rRNA_for_plot %>% 
  group_by(low_res) %>% 
  summarise(stdev = sd(value, na.rm = T),
            var = var(value, na.rm = T),
            avg = mean(value, na.rm = T),
            diff = max(value, na.rm = T) - min(value, na.rm = T))

