library(tidyverse)
library(dplyr)
library(ggpubr)
library(reshape2)
library(viridis) 
library(ggrepel)

#CTD
#Filter out disease gene
CTD_chem_cgixns<-CTD_chem_cgixns_1698339792344[which(CTD_chem_cgixns_1698339792344$GeneSymbol %in% unique_values),]

PFDA_gene<-CTD_chem_cgixns%>%
  mutate(ChemicalName = case_when(
    ChemicalName=="perfluorodecanoic acid" ~"PFDA",
    ChemicalName=="perfluorohexanesulfonic acid" ~"PFHxS",
    ChemicalName=="perfluorooctanoic acid" ~"PFOA",
    ChemicalName=="perfluorooctane sulfonic acid" ~"PFOS"
  ))%>%
  filter(ChemicalName=="PFDA")

PFDA_gene_unique<-as.data.frame(unique(PFDA_gene$GeneSymbol))

write.csv(PFDA_gene_unique,"PFDA_gene_unique.csv")



PFHxS_gene<-CTD_chem_cgixns%>%
  mutate(ChemicalName = case_when(
    ChemicalName=="perfluorodecanoic acid" ~"PFDA",
    ChemicalName=="perfluorohexanesulfonic acid" ~"PFHxS",
    ChemicalName=="perfluorooctanoic acid" ~"PFOA",
    ChemicalName=="perfluorooctane sulfonic acid" ~"PFOS"
  ))%>%
  filter(ChemicalName=="PFHxS")

PFHxS_gene_unique<-as.data.frame(unique(PFHxS_gene$GeneSymbol))

write.csv(PFHxS_gene_unique,"PFHxS_gene_unique.csv")



PFOA_gene<-CTD_chem_cgixns%>%
  mutate(ChemicalName = case_when(
    ChemicalName=="perfluorodecanoic acid" ~"PFDA",
    ChemicalName=="perfluorohexanesulfonic acid" ~"PFHxS",
    ChemicalName=="perfluorooctanoic acid" ~"PFOA",
    ChemicalName=="perfluorooctane sulfonic acid" ~"PFOS"
  ))%>%
  filter(ChemicalName=="PFOA")

PFOA_gene_unique<-as.data.frame(unique(PFOA_gene$GeneSymbol))

write.csv(PFOA_gene_unique,"PFOA_gene_unique.csv")


PFOS_gene<-CTD_chem_cgixns%>%
  mutate(ChemicalName = case_when(
    ChemicalName=="perfluorodecanoic acid" ~"PFDA",
    ChemicalName=="perfluorohexanesulfonic acid" ~"PFHxS",
    ChemicalName=="perfluorooctanoic acid" ~"PFOA",
    ChemicalName=="perfluorooctane sulfonic acid" ~"PFOS"
  ))%>%
  filter(ChemicalName=="PFOS")

PFOS_gene_unique<-as.data.frame(unique(PFOS_gene$GeneSymbol))

write.csv(PFOS_gene_unique,"PFOS_gene_unique.csv")


write.csv(PFOS_gene_unique,"PFOS_gene_unique.csv")
write.csv(PFHxS_gene,"PFHxS_gene.csv")
write.csv(PFOA_gene,"PFOA_gene.csv")
write.csv(PFOS_gene,"PFOS_gene.csv")


#GO terms
CTD_chem_go<-CTD_chem_go_enriched_1698339852934%>%
  mutate(ChemicalName = case_when(
    ChemicalName=="perfluorodecanoic acid" ~"PFDA",
    ChemicalName=="perfluorohexanesulfonic acid" ~"PFHxS",
    ChemicalName=="perfluorooctanoic acid" ~"PFOA",
    ChemicalName=="perfluorooctane sulfonic acid" ~"PFOS"
  ))

count_CTD_chem_go<-CTD_chem_go%>%
  group_by(ChemicalName)%>%
  summarize(count = n())

common_CTD_chem_go <- CTD_chem_go %>%
  group_by(ChemicalName) %>%
  summarize(Common_GoTermName = list(unique(GoTermName)))

common_CTD_chem_go <- reduce(common_CTD_chem_go$Common_GoTermName, intersect)
print(common_CTD_chem_go)

write.csv(CTD_chem_go,"output/CTD_chem_go.csv")

#Sleep only
filtered_data <- CTD_chem_go %>%
  filter(grepl("circadian|sleep", GoTermName, ignore.case = TRUE))
write.csv(filtered_data,"output/CTD_chem_go_sleep.csv")
filtered_data$qscore<- -log(filtered_data$CorrectedPValue, base=10)




#Tox21
my_colors <- c("#E69F00", "#3182bd", "#CC79A7", "#33a02c")
p1 <- ggplot(data = activity_npod_subset_df_melt_2, aes(x = Cell.Type, y = value, color = id, label = Assay.Target)) +
  geom_point(size=2.5) +
  geom_text_repel(size = 4) +
  labs(x = "Cell Type", y = "Activity", color = "Chemical") +
  theme_minimal() +
  theme(text= element_text(size=14))+
  scale_color_manual(values = my_colors)
p1

