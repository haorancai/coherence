library(tidyverse)
library(broman)
library(wesanderson)
library(patchwork)
library(ggrepel)
library(nord)
library(ggpubr)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#load("../Data/CCF_lag_2.rdata")
load("../Data/CCF_lag_1.rdata")


CCF <- CCF_lag_1 %>% rename(Class = CLASS) %>% 
  mutate(Class = replace(Class, Class == 'CONTROL', 'Control')) %>% 
  mutate(Class = replace(Class, Class == 'HEAT', 'Heat')) %>% 
  mutate(Class = replace(Class, Class == 'DROUGHT', 'Drought'))


CCF<- CCF %>% 
  rename(regulator = X1, target = X2)


#get TF family

TF_family <- read_delim("gene_family.txt", 
                        "\t", escape_double = FALSE, trim_ws = TRUE)
TF_family <- TF_family %>% as_tibble() %>% 
  rename(gene = `Protein ID`) %>% 
  mutate(gene = str_sub(gene, 1,14)) %>% dplyr::select(gene, Family) 
TF_family <- TF_family %>% as_tibble() %>% unique()

#Get drought related TF 

geneKeyword_table <- read_delim("geneKeyword.table.txt", 
                                "\t", escape_double = FALSE, trim_ws = TRUE)
geneKeyword_table <- geneKeyword_table %>% select(-Title) %>% 
  unique() 

#Get gene symble
OryzabaseGeneListEn_20210114010048 <- read_delim("gene_table.txt", 
                                                 "\t", escape_double = FALSE, trim_ws = TRUE)
gene <- OryzabaseGeneListEn_20210114010048 %>% select(`MSU ID`, `CGSNL Gene Symbol`, `Gene symbol synonym(s)`, 'Gene Ontology','Trait Ontology', 'Plant Ontology') 
colnames(gene) <- c('MSU_ID', 'gene_symbol', 'synonym', 'Gene_Ontology', 'Trait_Ontology', 'Plant_Ontology')

gene <- gene %>% filter(is.na(MSU_ID)== FALSE)
gene %>% mutate(MSU_ID = str_sub(MSU_ID, 1,14))

#get differential expression in rice gene expression time series
DE_time_point <- read_csv("DE_time_point.csv", 
                          skip = 3)
DE <- DE_time_point %>% mutate(Heat = HEAT_UP - HEAT_DOWN, Drought = DROUGHT_UP - DROUGHT_DOWN) %>% 
  select(gene, Heat, Drought) 
de_drought <- DE %>% select(gene,Drought)
de_heat <- DE %>% select(gene, Heat)
# main+fig_5 -------------------------------------------------------------------




my_comparisons <- list( c("HSF_family", "non_HSF_family"))


dummy<-  CCF %>% as_tibble() %>% 
  filter(regulator != target) %>% 
  mutate_at("max_ccf", as.numeric) %>% 
  mutate(max_ccf = abs(max_ccf)) %>% 
  filter(Class == 'Heat' | Class =='Control') %>% 
  group_by(regulator,target,Class) %>% 
  summarise( max_ccf = mean(max_ccf,na.rm = TRUE)) %>% 
  ungroup() %>% 
  pivot_wider(c(regulator,target), values_from = max_ccf, names_from = Class) %>% 
  left_join(TF_family, by = c('regulator' ='gene')) %>% 
  mutate(Category = ifelse(Family == "HSF","HSF_family", "non_HSF_family")) %>% 
  select(regulator,target, Control,Heat, Category) %>% 
  #mutate(diff = ifelse(abs(HEAT)>0.54, abs(HEAT) - 0.54, 0)) %>% 
  mutate(diff = abs(Heat) - abs(Control)) %>% 
  #arrange(desc(diff)) %>% view()
  group_by(regulator, Category) %>% 
  summarise(diff_mean = mean(diff, na.rm = TRUE)) %>% 
  filter(is.na(Category) == FALSE) 




fig_5_1 <-
  dummy %>% 
  ungroup() %>% 
  mutate(
    Category = ifelse(Category == "HSF_family", "HSF family", "non-HSF family")
  ) %>% 
  ggstatsplot::ggbetweenstats(
    x = Category,
    y = diff_mean,
    plot.type = 'violin',
    title = "Heat",
    results.subtitle = FALSE,
    pairwise.comparisons = TRUE,
    package = 'nord',
    palette = 'algoma_forest',
    centrality.label.args = list(size  = 5,nudge_x = 0.4, segment.linetype = 4)
  ) +
  stat_compare_means(label.y = 0.90, label.x = 0.63,size = 6) +
  theme_bw()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.title = element_text(size = 18), 
        axis.text = element_text(size = 12))+
  theme(plot.title = element_text(hjust = 0.5, size = 20))+
  labs(
    x = ' ',
    y = 'Differential MCC'
  )



dummy <- CCF %>% as_tibble() %>% 
  filter(regulator != target) %>% 
  mutate_at("max_ccf", as.numeric) %>% 
  mutate(max_ccf = abs(max_ccf)) %>% 
  filter(Class == 'Drought' | Class =='Control') %>% 
  group_by(regulator,target,Class) %>% 
  summarise( max_ccf = mean(max_ccf,na.rm = TRUE)) %>% 
  ungroup() %>% 
  pivot_wider(c(regulator,target), values_from = max_ccf, names_from = Class) %>% 
  left_join(geneKeyword_table %>% group_by(MSU) %>%  summarize(Keyword=paste(Keyword,collapse=",")), by = c('regulator' ='MSU')) %>% 
  mutate(Category = ifelse((str_detect(Keyword, "ABA")  | str_detect(Keyword, "drought")),"Known_drought_related_TF", "Unknown_TF")) %>% 
  select(regulator,target, Control,Drought, Category) %>% 
  mutate(diff = abs(Drought)- abs(Control) ) %>% 
  group_by(regulator, Category) %>% 
  summarise(diff_mean = mean(diff,na.rm = TRUE)) %>% 
  mutate(Category = ifelse((is.na(Category) == TRUE | Category == 'Unknown_TF'), 'Uncharacterized_TF', Category))





fig_5_2 <- dummy %>% 
  ungroup() %>% 
  mutate(Category = ifelse(Category == 'Known_drought_related_TF',"Known drought-related TF",  "Uncharacterized_TF")) %>% 
  ggstatsplot::ggbetweenstats(
    x = Category,
    y = diff_mean,
    plot.type = 'violin',
    title = "Drought",
    results.subtitle = FALSE,
    pairwise.comparisons = TRUE,
    package = 'nord',
    palette = 'algoma_forest',
    centrality.label.args = list(size  = 5,nudge_x = 0.4, segment.linetype = 4)
  ) +
  stat_compare_means(label.y = 0.90, label.x = 0.60,size = 6) +
  theme_bw()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 12))+
  theme(plot.title = element_text(hjust = 0.5, size = 20))+
  labs(
    x = ' ',
    y = 'Differential MCC'
  )






fig_5_3 <- CCF %>% as_tibble() %>% 
  mutate_at("max_ccf", as.numeric) %>% 
  filter(Class == 'Heat' | Class =='Control') %>% 
  group_by(regulator,target,Class) %>% 
  summarise(max_ccf = mean(max_ccf,na.rm = TRUE)) %>% 
  ungroup() %>% 
  pivot_wider(c(regulator,target), values_from = max_ccf, names_from = Class) %>% 
  mutate(diff = abs(Heat) - abs(Control)) %>% 
  group_by(regulator) %>% 
  summarise(diff_mean = mean(diff)) %>% 
  left_join(de_heat, by = c('regulator' ='gene')) %>% 
  mutate(Heat = ifelse(is.na(Heat), 0, Heat) ) %>% 
  left_join(gene %>% mutate(MSU_ID = str_sub(MSU_ID, 1,14)), by = c('regulator' = 'MSU_ID')) %>% 
  left_join(TF_family, by = c('regulator' ='gene')) %>% 
  mutate(label = ifelse(diff_mean > 0.34, gene_symbol,'')) %>% 
  mutate(label = ifelse(label == '_', synonym, label)) %>% 
  mutate(Category = ifelse(Family == "HSF","HSF_family", "non_HSF_family")) %>% 
  arrange(desc(diff_mean)) %>% 
  arrange(desc(Category))%>% 
  filter(is.na(Category) == FALSE) %>% 
  ggplot(aes(Heat, diff_mean, label = label)) + 
  geom_text_repel(box.padding = 0.5, max.overlaps = Inf)+
  geom_point(aes(color = Category),alpha = .7, size = 2)+ 
  ylab('Differential MCC aganst Control')+
  xlab('Number of time points \n with differential expression')+
  labs(color = ' ', shape = '')+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(legend.position = c(0.70, 0.100),
        legend.text = element_text(size = 15))+ 
  theme(axis.title= element_text(size = 18), 
        axis.text = element_text(size = 12))+
  ggtitle('Heat')+
  theme(plot.title = element_text(hjust = 0.5, size = 20))+
  theme(legend.title= element_blank(), legend.background = element_rect(linetype = 'solid', color = 'black'))+
  scale_color_manual(values = c("#BF616A", "#868686CC"),labels = c("HSF Family", "Non-HSF Family"))+
  theme(aspect.ratio = 1/1)




fig_5_4 <- CCF %>% as_tibble() %>% 
  mutate_at("max_ccf", as.numeric) %>% 
  filter(Class == 'Drought' | Class =='Control') %>% 
  group_by(regulator,target,Class) %>% 
  summarise(max_ccf = mean(max_ccf,na.rm = TRUE)) %>% 
  ungroup() %>% 
  pivot_wider(c(regulator,target), values_from = max_ccf, names_from = Class) %>% 
  mutate(diff = abs(Drought) - abs(Control)) %>% 
  group_by(regulator) %>% 
  summarise(diff_mean = mean(diff)) %>% 
  left_join(de_drought, by = c('regulator' ='gene')) %>% 
  mutate(Drought = ifelse(is.na(Drought), 0, Drought) ) %>% 
  left_join(gene %>% mutate(MSU_ID = str_sub(MSU_ID, 1,14)), by = c('regulator' = 'MSU_ID')) %>% 
  arrange(desc(diff_mean)) %>% 
  mutate(rank = row_number()) %>% 
  mutate(label = ifelse(rank < 14, gene_symbol,'')) %>% 
  mutate(label = ifelse(label == '_', synonym, label)) %>% 
  dplyr::rename(Drought_response = Drought)




fig_5_4 <- fig_5_4 %>% 
  ggplot(aes(x = Drought_response, y = diff_mean, label = label)) + 
  geom_point(alpha = .5, size = 2)+
  geom_text_repel(box.padding = 0.5, max.overlaps = 20)+
  ylab('Differential MCC against Control')+
  xlab('Number of time points \n with differential expression')+
  theme_bw()+
  ggtitle('Drought response')+
  theme(plot.title = element_text(hjust = 0.5, size = 20))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(legend.position = c(0.85, 0.90),
        legend.text = element_text(size = 15))+ 
  theme(axis.title= element_text(size = 18), 
        axis.text = element_text(size = 12))+
  theme(legend.title= element_blank())+
  scale_color_manual(values = c("#EFC000CC", "#868686CC"),labels = c("HSF Family", "Non-HSF Family"))+
  theme(aspect.ratio = 1/1)




main_fig_5<- (fig_5_1 | fig_5_2) / (fig_5_3 + fig_5_4) + plot_annotation(tag_levels = 'A')
ggsave(main_fig_5, filename = 'main_fig_5.pdf',height = 15 , width = 15)














 

# main_fig_5_suppl --------------------------------------------------------


fig_5_suppl<- CCF %>% as_tibble() %>% 
  mutate_at('max_ccf', as.numeric) %>% 
  filter(regulator != target) %>% 
  group_by(regulator,target, Class) %>%
  summarise(mean = mean(max_ccf, na.rm = TRUE)) %>%
  ungroup() %>%
  inner_join(TF_family, by = c('regulator' ='gene')) %>% 
  ggplot(aes(x = Family, y = mean, fill = Class, color = Class))+ 
  geom_boxplot( size = 0.7)+
  scale_fill_manual(values = c("#81A1C1", "#A3BE8C",'#B48EAD')) +
  scale_color_manual(values = c("#81A1C1", "#A3BE8C",'#B48EAD')) +
  xlab('')+
  ylab('Mean MCC')+ 
  theme_light()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.title.y = element_text(size = 18))+
  theme(legend.title= element_blank(), 
        legend.text = element_text(size = 18),
        legend.position = 'top') +
  theme(axis.text.x = element_text( color="black", size=15, angle=45))




ggsave(fig_5_suppl, filename = 'fig_5_suppl.pdf', height = 7, width = 14)












