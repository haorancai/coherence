
library(tidygraph)  
library(tidyverse)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
RNA_seq_summary_Experiment_ID <- read_table2("../Data/RNA-seq summary_Experiment_ID.txt", skip=  6) 
RNA_seq_summary_Experiment_ID <- RNA_seq_summary_Experiment_ID %>% filter(is.na(MINUTES) == FALSE)         



network_prior <- read_delim("../Data/network_prior.txt", 
                            "\t", escape_double = FALSE, col_names = FALSE, 
                            trim_ws = TRUE)




network_prior<- network_prior %>% rename(from = X1, to = X2)


core_net<- network_prior %>% 
  filter(to %in% from)



# Get TF activity ---------------------------------------------------------


load('../Data/rice_tidy_RNA_seq.rdata')




subset <- data %>% 
  filter( CLASS == 'DROUGHT' ) %>% 
  select(-key, -CULTIVAR) %>% 
  nest(value = c(value))



df <- RNA_seq_summary_Experiment_ID %>%
  select(CLASS, MINUTES) %>% 
  filter(CLASS == 'DROUGHT' ) %>% 
  filter(MINUTES < 140) %>% 
  unique() %>% 
  mutate(prior = map(CLASS, ~network_prior)) %>%
  unnest(prior) %>%
  left_join(subset, by = c('from' = 'gene', 'CLASS' = 'CLASS', 'MINUTES' = 'MINUTES')) %>%
  dplyr::rename(regulator_value = value) %>% 
  left_join(subset, by = c('to' = 'gene', 'CLASS' = 'CLASS', 'MINUTES' = 'MINUTES')) %>%
  dplyr::rename(target_value = value) 



df <- df %>%
  ungroup() %>% 
  rowwise() %>% 
  mutate(cor = cor(regulator_value %>% 
                     pull(),
                   target_value %>% 
                     pull(), use = 'everything',method = 'pearson')) 



drought <- df 




subset <- data %>% 
  filter( CLASS == 'CONTROL' ) %>% 
  select(-key, -CULTIVAR) %>% 
  nest(value = c(value))



df <- RNA_seq_summary_Experiment_ID %>%
  select(CLASS, MINUTES) %>% 
  filter(CLASS == 'CONTROL' ) %>% 
  filter(MINUTES < 140) %>% 
  unique() %>% 
  mutate(prior = map(CLASS, ~network_prior)) %>%
  unnest(prior) %>%
  left_join(subset, by = c('from' = 'gene', 'CLASS' = 'CLASS', 'MINUTES' = 'MINUTES')) %>%
  dplyr::rename(regulator_value = value) %>% 
  left_join(subset, by = c('to' = 'gene', 'CLASS' = 'CLASS', 'MINUTES' = 'MINUTES')) %>%
  dplyr::rename(target_value = value) 



df <- df %>%
  ungroup() %>% 
  rowwise() %>% 
  mutate(cor = cor(regulator_value %>% 
                     pull(),
                   target_value %>% 
                     pull(), use = 'everything',method = 'pearson')) 



control <- df 



# Get paired t test comparison of activity --------------------------------





pair_t <- control %>% 
  mutate(cor = ifelse(is.na(cor) == TRUE, 0, cor)) %>% 
  mutate(cor = abs(cor)) %>%  
  group_by(MINUTES, from) %>% 
  summarise(mean_control = mean(cor)) %>% 
  left_join(drought %>% 
              mutate(cor = ifelse(is.na(cor) == TRUE, 0, cor)) %>% 
              mutate(cor = abs(cor)) %>%  
              group_by(MINUTES, from) %>% 
              summarise(mean_drought = mean(cor)), by = c('MINUTES' = 'MINUTES', "from" = 'from')) %>% 
  ungroup() %>% 
  select(-MINUTES) %>% 
  nest(control = mean_control, drought = mean_drought) %>% 
  rowwise() %>% 
  mutate(pair_t = list(t.test(control$mean_control, drought$mean_drought, paired = TRUE, alternative = 'less')$p.value))





# Bottom up hierarchy -----------------------------------------------------



level_1 <- subset(core_net, !(to %in% from)) %>% 
  select(to) %>% 
  unique()


level_2 <- subset(core_net, to %in% level_1$to) %>% 
  select(from) %>% 
  unique() %>% 
  filter( !(from %in% level_1$to) ) 


level_3 <- subset(core_net, to %in% level_2$from) %>% 
  select(from) %>% 
  unique() %>% 
  filter(!(from %in% level_2$from | from %in% level_1$to))

level_4 <- subset(core_net, to %in% level_3$from) %>% 
  select(from) %>% 
  unique() %>% 
  filter(!(from %in% level_2$from | from %in% level_1$to | from %in% level_3$from))



TF_hierarchy_bottom_up <- bind_rows(bind_rows(bind_rows(level_1 %>% rename(name = to) %>% mutate(level = 1), 
                                                        level_2 %>% rename(name = from) %>% mutate(level = 2)),
                                              level_3 %>% rename(name = from) %>% mutate(level = 3)),
                                    level_4 %>% rename(name = from) %>% mutate(level = 4))




# Stem input --------------------------------------------------------------



stem_activity<-  drought %>% 
  mutate(cor = ifelse(is.na(cor) == TRUE, 0, cor)) %>% 
  mutate(cor = abs(cor)) %>%  
  group_by(MINUTES, X1) %>% 
  summarise(mean = mean(cor), n = n()) %>% 
  filter(n > 2) %>% 
  left_join(TF_hierarchy, by = c('X1' = 'name')) %>% 
  ungroup() %>% 
  left_join(
    pair_t %>% mutate(p.value = as.numeric(pair_t))) %>% 
  filter(p.value < 0.05) %>% 
  select(X1, MINUTES,mean) %>% 
  mutate(Hour = MINUTES / 60) %>% 
  mutate(Hour = paste0(Hour,"h")) %>% 
  pivot_wider(X1, names_from = 'Hour', values_from = 'mean')


write.table(stem_activity %>% rename(Gene = X1), file = 'stem_activity.txt',row.names = FALSE,sep = '\t',quote = FALSE)  




# height_based hierarchy --------------------------------------------------



# as_tbl_graph(core_net)  %>% 
#   mutate(out_degree = centrality_degree(mode = 'out', loops = TRUE)) %>% 
#   mutate(in_degree = centrality_degree(mode = 'in', loops = TRUE)) %>% 
#   mutate(height = (out_degree - in_degree )/ (out_degree + in_degree)) %>% 
#   as_tibble() %>% 
#   left_join(TF_hierarchy_bottom_up) %>% 
#   mutate(level = factor(level)) %>% 
# ggplot(aes(level, height))+geom_boxplot()+theme_bw()
# 
# 
# 
# 
# as_tbl_graph(network_prior)  %>% 
#   mutate(out_degree = centrality_degree(mode = 'out', loops = TRUE)) %>% 
#   mutate(in_degree = centrality_degree(mode = 'in', loops = TRUE)) %>% 
#   mutate(height = (out_degree - in_degree )/ (out_degree + in_degree)) %>% 
#   as_tibble() %>% 
#   left_join(TF_hierarchy_bottom_up) %>% 
#   mutate(level = factor(level)) %>% 
#   ggplot(aes(level, height))+geom_boxplot()+theme_bw()

library(tidygraph)
TF_hierarchy <- 
  as_tbl_graph(core_net)  %>% 
  mutate(out_degree = centrality_degree(mode = 'out', loops = TRUE)) %>% 
  mutate(in_degree = centrality_degree(mode = 'in', loops = TRUE)) %>% 
  mutate(height = (out_degree - in_degree )/ (out_degree + in_degree)) %>% 
  as_tibble() %>% 
  left_join(TF_hierarchy_bottom_up) %>% 
  mutate(level = factor(level)) 







TF_hierarchy <- TF_hierarchy %>% 
  mutate(height = case_when(height == 1~ 4,
                            height <1 & height >0~ 3,
                            height < 0.000001 & height >-1 ~ 2,
                            height == -1 ~ 1))  





# gene expression level ---------------------------------------------------

fig_5_2 <- drought %>% 
  select(-cor,-target_value,-to) %>% 
  unique() %>% 
  rowwise() %>% 
  mutate(mean_count = mean(regulator_value$value)) %>% 
  select(-regulator_value) %>% 
  left_join(
    pair_t %>% mutate(p.value = as.numeric(pair_t))) %>% 
  filter(p.value < 0.05) %>%
  left_join(TF_hierarchy,by = c('from' = 'name')) %>% 
  filter(is.na(level) == FALSE) %>% 
  ungroup() %>% 
  group_by(level, MINUTES) %>% 
  summarise(mean = mean(mean_count), n = n()) %>% 
  mutate(Num = as.factor(paste(n, 'TFs'))) %>% 
  mutate(Level = as.factor(level)) %>% 
  ggplot(aes(MINUTES, mean, color = Level))+ 
  geom_line(size = 1.5)+
  theme_light()+
  theme(legend.title = element_blank(), legend.position = 'top')+
  xlab('Time (min)')+
  ylab('Mean Expression Level')+
  geom_point(size = 2.5)+
  scale_color_manual(values= c('#8FBCBB','#EBCB8B','#81A1C1','#D08770')) +
  theme(
    axis.title.x = element_text(color="black", size=15),
    axis.title.y = element_text(color="black", size=15),
    legend.text = element_text(color = 'black', size = 12)
  )+
  gghighlight::gghighlight(label_key = Num , label_params = list(size = 8))+
  theme(
    axis.text.y=element_blank(),
    axis.text.x = element_text(size = 12),
    axis.ticks.y=element_blank())+
  theme(aspect.ratio=1/1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())



# activity at drought -----------------------------------------------------


fig_5_3 <- 
drought %>% 
  mutate(cor = ifelse(is.na(cor) == TRUE, 0, cor)) %>% 
  mutate(cor = abs(cor)) %>%  
  group_by(MINUTES, from) %>% 
  summarise(mean = mean(cor), n = n()) %>% 
  filter(n > 2) %>% 
  left_join(TF_hierarchy ,by = c('from' = 'name')) %>% 
  filter(is.na(level) == FALSE) %>% 
  ungroup() %>% 
  left_join(
    pair_t %>% mutate(p.value = as.numeric(pair_t))) %>% 
  filter(p.value < 0.05) %>%
  group_by(MINUTES, level) %>% 
  summarise(mean_activity = mean(mean), n = n()) %>%
  mutate(Num = as.factor(paste(n, 'TFs'))) %>% 
  mutate(Level = as.factor(level)) %>% 
  ggplot(aes(MINUTES, mean_activity, color = Level))+ 
  geom_line(size = 1.5)+
  geom_point(size = 2)+
  theme_light()+
  scale_color_manual(values= c('#8FBCBB','#EBCB8B','#81A1C1','#D08770'), label = c('1 (Bottom)','2','3','4 (Top)')) +
  ylim(0.27,0.67)+
  theme(legend.title = element_blank(), legend.position = 'top')+
  xlab('Time (min)')+
  ylab('TF Activity')+
  geom_point(size = 2)+
  theme(
    axis.title.x = element_text(color="black", size=15),
    axis.title.y = element_text(color="black", size=15),
    legend.text = element_text(color = 'black', size = 20)
  )+
  theme(
    axis.text.y=element_blank(),
    axis.text.x = element_text(size = 12),
    axis.ticks.y=element_blank())+
  annotate(
    geom = "curve", x = 90, y = 0.55, xend = 100, yend = 0.42, 
    curvature = .4, arrow = arrow(length = unit(3, "mm"))
  ) +
  annotate(geom = "text", x = 103, y = 0.42, label = "Second Wave", hjust = "left",size = 6)+
  annotate(
    geom = "curve", x = 60, y = 0.59, xend = 50, yend = 0.62, 
    curvature = .4, arrow = arrow(length = unit(3, "mm"))
  ) +
  annotate(geom = "text", x = 51, y = 0.63, label = "First Wave", hjust = "left",size = 6)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())







# Network_visualization ---------------------------------------------------


library(ggraph)
library(igraph)
library(colormap)
mycolor <- colormap(colormap=colormaps$viridis, nshades=max(TF_hierarchy_bottom_up$level))
mycolor <- sample(mycolor, length(mycolor))





graph <- graph_from_data_frame(core_net,directed = TRUE, vertices = TF_hierarchy %>% arrange(desc(level)))


fig_5_1 <- ggraph(graph, layout="linear") + 
  geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
  geom_node_point(aes(size=0.01, color=as.factor(level), fill=level), alpha=0.5) +
  scale_size_continuous(range=c(0.5,8)) +
  scale_color_manual(values= c('#8FBCBB','#EBCB8B','#81A1C1','#D08770')) +
  geom_node_text(aes(label= ''), angle=65, hjust=1, nudge_y = -1.1, size=2) +
  theme_void() +
  theme(
    legend.position="none",
    plot.margin=unit(c(0,0,0.4,0), "null"),
    panel.spacing=unit(c(0,0,3.4,0), "null")
  ) +
  expand_limits(x = c(-1.2, 1.2), y = c(-5.6, 1.2)) +
  theme(aspect.ratio = 2/5)







# Clustering of activity --------------------------------------------------







# continuously increase ---------------------------------------------------


fig_5_4 <- drought %>% 
  mutate(cor = ifelse(is.na(cor) == TRUE, 0, cor)) %>% 
  mutate(cor = abs(cor)) %>%  
  group_by(MINUTES, from) %>% 
  summarise(mean = mean(cor), n = n()) %>% 
  filter(n > 2) %>% 
  filter(from %in% c('LOC_Os01g09640', 
                     'LOC_Os03g12860',
                     'LOC_Os04g48070',
                     "LOC_Os05g34050",
                     'LOC_Os09g38340',
                     'LOC_Os12g03050',
                     'LOC_Os09g38340',
                     'LOC_Os03g20550',
                     'LOC_OS01G66120'
  )) %>% 
  left_join(TF_hierarchy, by = c('from'= 'name')) %>% 
  filter(is.na(level) == FALSE) %>% 
  ggplot(aes(MINUTES, mean, group = from, color = level))+ 
  geom_line(size = 1.5)+
  geom_point(size = 2)+
  theme_light()+
  ylim(0,0.9)+
  theme(legend.title = element_blank(), legend.position = 'none')+
  xlab('Time (min)')+
  ylab('TF Activity')+
  geom_point(size = 2)+
  scale_color_manual(values= c('#EBCB8B','#81A1C1','#D08770')) +
  theme(
    axis.title.x = element_text(color="black", size=15),
    axis.title.y = element_text(color="black", size=15),
    legend.text = element_text(color = 'black', size = 15)
  )+
  theme(
    axis.text.y=element_blank(),
    axis.text.x = element_text(size = 12),
    axis.ticks.y=element_blank())+
  theme(aspect.ratio = 1/1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
  
  


# 1.00	LOC_OS01G09640	ID_0	0.00	0.35	0.17	0.41	0.21	0.57	0.79	0.84	1.16
# 1.00	LOC_OS02G43330	ID_17	0.00	0.52	0.76	0.90	0.80	1.14	1.20	1.12	1.27
# 1.00	LOC_OS04G48070	ID_42	0.00	0.37	0.32	0.93	1.05	1.08	1.03	0.78	1.03
# 1.00	LOC_OS03G12860	ID_28	0.00	-0.73	-0.35	0.25	0.05	0.29	0.37	0.44	0.44
# 1.00	LOC_OS03G12350	ID_26	0.00	0.62	0.20	0.73	0.74	0.57	1.20	0.94	1.24
# 1.00	LOC_OS05G34050	ID_47	0.00	0.04	0.48	0.87	0.76	0.94	0.91	0.85	1.01
# 1.00	LOC_OS07G44950	ID_59	0.00	0.37	0.70	1.01	0.80	0.92	1.18	0.94	1.19
# 1.00	LOC_OS09G38340	ID_79	0.00	-0.12	0.46	0.57	0.91	1.09	1.11	1.00	1.14
# 1.00	LOC_OS12G03050	ID_89	0.00	-0.24	0.35	0.12	0.53	0.72	0.49	0.53	0.79
# 1.00	LOC_OS03G20550	ID_29	0.00	0.50	0.91	1.03	0.86	0.95	0.82	0.94	1.04
# 1.00	LOC_OS01G66120	ID_10	0.00	0.51	0.76	1.07	0.76	0.90	0.91	0.94	1.00
# 1.00	LOC_OS01G09640	ID_0	0.00	0.35	0.17	0.41	0.21	0.57	0.79	0.84	1.16
# 1.00	LOC_OS03G12350	ID_26	0.00	0.62	0.20	0.73	0.74	0.57	1.20	0.94	1.24
# 1.00	LOC_OS09G38340	ID_79	0.00	-0.12	0.46	0.57	0.91	1.09	1.11	1.00	1.14
# 1.00	LOC_OS08G43550	ID_70	0.00	-0.36	-0.55	-0.20	-0.07	-0.18	0.42	-0.02	0.54














# first wave driver -------------------------------------------------------

fig_5_6 <- drought %>% 
  mutate(cor = ifelse(is.na(cor) == TRUE, 0, cor)) %>% 
  mutate(cor = abs(cor)) %>%  
  group_by(MINUTES, from) %>% 
  summarise(mean = mean(cor), n = n()) %>% 
  filter(n > 2) %>% 
  filter(from %in% c('LOC_Os10g01470', 
                     'LOC_Os02g32590',
                     'LOC_Os07g48596',
                     'LOC_Os04g48070'
  )) %>% 
  left_join(TF_hierarchy, by = c('from'= 'name')) %>% 
  filter(is.na(level) == FALSE) %>% 
  mutate(from = ifelse(from == 'LOC_Os02g32590', 'HSFA3', 'ROC4')) %>% 
  ggplot(aes(MINUTES, mean, color = level, group = from))+ 
  geom_line(size = 1.5)+
  geom_point(size = 2)+
  theme_light()+
  ylim(0,0.8)+
  theme(legend.title = element_blank(), legend.position = 'top')+
  xlab('Time (min)')+
  ylab('TF Activity')+
  geom_point(size = 2)+
  scale_color_manual(values= c('#8FBCBB','#D08770')) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.text = element_text(color = 'black', size = 15)
  )+
  gghighlight::gghighlight(label_key = from, label_params = list(size = 7))+#, MINUTES %in% 60:90)+
  theme(
    axis.text.y=element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())






# 1.00	LOC_OS10G01470	ID_80	0.00	0.02	-0.63	0.31	0.21	0.56	0.50	0.39	0.49
# 1.00	LOC_OS02G32590	ID_14	0.00	-0.37	0.40	0.56	0.69	0.44	0.50	0.28	0.45
# 1.00	LOC_OS02G47660	ID_20	0.00	0.33	0.48	1.04	0.86	1.12	0.97	0.78	0.80
# 1.00	LOC_OS04G44670	ID_39	0.00	-0.24	0.68	1.17	1.05	1.24	0.69	1.06	0.95
# 1.00	LOC_OS07G48596	ID_62	0.00	-0.20	-0.46	0.59	0.29	0.37	0.40	0.39	0.51
# 



c <- igraph::all_simple_paths(graph, from = 'LOC_Os04g48070', to = 'LOC_Os02g32590', cutoff = 50)


# second wave driver ------------------------------------------------------



fig_5_7 <- drought %>% 
  mutate(cor = ifelse(is.na(cor) == TRUE, 0, cor)) %>% 
  mutate(cor = abs(cor)) %>%  
  group_by(MINUTES, from) %>% 
  summarise(mean = mean(cor), n = n()) %>% 
  filter(n > 2) %>% 
  filter(from %in% c('LOC_Os03g12350', 
                     'LOC_Os08g43550'
  )) %>% 
  left_join(TF_hierarchy, by = c('from'= 'name')) %>% 
  mutate(from = ifelse(from == 'LOC_Os03g12350', 'PR21', 'GL1A')) %>% 
  ggplot(aes(MINUTES, mean, color = level, group = from))+ 
  geom_line(size = 1.5)+
  geom_point(size = 2)+
  theme_light()+
  ylim(0,0.8)+
  theme(legend.title = element_blank(), legend.position = 'top')+
  xlab('Time (min)')+
  ylab('TF Activity')+
  geom_point(size = 2)+
  scale_color_manual(values= c('#8FBCBB','#D08770')) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.text = element_text(color = 'black', size = 15)
  )+
  gghighlight::gghighlight(label_key = from, label_params = list(size = 7))+#, MINUTES %in% 60:90)+
  theme(
    axis.text.y=element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())



# 1.00	LOC_OS03G12350	ID_26	0.00	0.62	0.20	0.73	0.74	0.57	1.20	0.94	1.24
# 1.00	LOC_OS08G43550	ID_70	0.00	-0.36	-0.55	-0.20	-0.07	-0.18	0.42	-0.02	0.54


#d <-distances(graph)



d <- igraph::all_simple_paths(graph, from = 'LOC_Os08g43550', to = 'LOC_Os03g12350', cutoff = 6)




drought %>%
  ungroup() %>% 
  filter(from == 'LOC_Os08g43550', to == 'LOC_Os09g32510')



d %>% as_tibble() %>% 
  mutate(from = rownames(d)) %>% 
  pivot_longer(-from,names_to = 'to',values_to = 'distance') %>% 
  filter(from == 'LOC_Os08g43550', to == 'LOC_Os03g12350')
  
  



# all rounder driver ------------------------------------------------------


fig_5_5 <- drought %>% 
  mutate(cor = ifelse(is.na(cor) == TRUE, 0, cor)) %>% 
  mutate(cor = abs(cor)) %>%  
  group_by(MINUTES, from) %>% 
  summarise(mean = mean(cor), n = n()) %>% 
  filter(n > 2) %>% 
  filter(from %in% c('LOC_Os05g49420', 
                     'LOC_Os01g46970',
                     'LOC_Os02g47660',
                     'LOC_Os07g44950',
                     'LOC_Os04g44670',
                     'LOC_Os09g29460',
                     'LOC_Os02g43330', 
                     'LOC_Os08g36790')) %>% 
  #'LOC_Os10g01470',  
  left_join(TF_hierarchy, by = c('from'= 'name')) %>% 
  filter(is.na(level) == FALSE) %>% 
  ggplot(aes(MINUTES, mean, group = from, color = level))+ 
  geom_line(size = 1.5)+
  geom_point(size = 2)+
  theme_light()+
  ylim(0,0.9)+
  theme(legend.position = 'none')+
  xlab('Time (min)')+
  ylab('TF Activity')+
  geom_point(size = 2)+
  scale_color_manual(values= c('#8FBCBB','#EBCB8B','#81A1C1','#D08770')) +
  theme(
    axis.title.x = element_text(color="black", size=15),
    axis.title.y = element_text(color="black", size=15),
    legend.text = element_text(color = 'black', size = 15)
  )+
  theme(
    axis.text.y=element_blank(),
    axis.text.x = element_text(size = 12),
    axis.ticks.y=element_blank())+
  theme(aspect.ratio = 1/1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank())

# 1.00	LOC_OS05G49420	ID_51	0.00	1.26	1.70	1.96	1.06	1.88	1.98	1.28	1.82
# 1.00	LOC_OS01G46970	ID_2	0.00	1.04	1.32	1.61	0.66	1.43	1.56	0.85	1.50

# 1.00	LOC_OS09G29460	ID_76	0.00	0.24	-0.18	0.46	0.26	0.83	1.03	0.47	0.60
# 1.00	LOC_OS08G36790	ID_66	0.00	-0.19	0.24	0.81	0.56	0.99	0.90	0.70	0.48





# aggregating different group of TFs --------------------------------------



OryzabaseGeneListEn_20210114010048 <- read_delim("gene_table.txt", 
                                                 "\t", escape_double = FALSE, trim_ws = TRUE)
gene <- OryzabaseGeneListEn_20210114010048 %>% select(`MSU ID`, `CGSNL Gene Symbol`, `Gene symbol synonym(s)`, 'Gene Ontology','Trait Ontology', 'Plant Ontology') 
colnames(gene) <- c('MSU_ID', 'gene_symbol', 'synonym', 'Gene_Ontology', 'Trait_Ontology', 'Plant_Ontology')

gene <- gene %>% filter(is.na(MSU_ID)== FALSE)

gene <- gene %>% mutate(MSU_ID = str_sub(MSU_ID, 1,14))

aggregate <- bind_rows(bind_rows(bind_rows(drought %>% 
  mutate(cor = ifelse(is.na(cor) == TRUE, 0, cor)) %>% 
  mutate(cor = abs(cor)) %>%  
  group_by(MINUTES, from) %>% 
  summarise(mean = mean(cor), n = n()) %>% 
  filter(n > 2) %>% 
  filter(from %in% c('LOC_Os05g49420', 
                     'LOC_Os01g46970',
                     'LOC_Os02g47660',
                     'LOC_Os07g44950',
                     'LOC_Os04g44670',
                     'LOC_Os09g29460',
                     'LOC_Os02g43330', 
                     'LOC_Os08g36790')) %>% 
  ungroup() %>% 
  select(from) %>% 
  unique() %>% 
  mutate(group = 'Both_two_waves'),
  drought %>% 
    mutate(cor = ifelse(is.na(cor) == TRUE, 0, cor)) %>% 
    mutate(cor = abs(cor)) %>%  
    group_by(MINUTES, from) %>% 
    summarise(mean = mean(cor), n = n()) %>% 
    filter(n > 2) %>% 
    filter(from %in% c('LOC_Os03g12350', 
                       'LOC_Os08g43550'
    )) %>% ungroup() %>% 
    select(from) %>% 
    unique() %>% 
    mutate(group = 'Second_wave_driver')
  
  ),
  drought %>% 
    mutate(cor = ifelse(is.na(cor) == TRUE, 0, cor)) %>% 
    mutate(cor = abs(cor)) %>%  
    group_by(MINUTES, from) %>% 
    summarise(mean = mean(cor), n = n()) %>% 
    filter(n > 2) %>% 
    filter(from %in% c('LOC_Os10g01470', 
                       'LOC_Os02g32590',
                       'LOC_Os07g48596',
                       'LOC_Os04g48070'
    )) %>% ungroup() %>% 
    select(from) %>% 
    unique() %>% 
    mutate(group = 'First_wave_driver')
  
  ), drought %>% 
    mutate(cor = ifelse(is.na(cor) == TRUE, 0, cor)) %>% 
    mutate(cor = abs(cor)) %>%  
    group_by(MINUTES, from) %>% 
    summarise(mean = mean(cor), n = n()) %>% 
    filter(n > 2) %>% 
    filter(from %in% c('LOC_Os01g09640', 
                       'LOC_Os03g12860',
                       'LOC_Os04g48070',
                       "LOC_Os05g34050",
                       'LOC_Os09g38340',
                       'LOC_Os12g03050',
                       'LOC_Os09g38340',
                       'LOC_Os03g20550',
                       'LOC_OS01G66120'
    )) %>% ungroup() %>% 
    select(from) %>% 
    unique() %>% 
    mutate(group = 'Continous_increasing')
  ) %>%
left_join(TF_hierarchy, by = c('from'= 'name')) %>% 
  left_join(gene, by = c('from' = 'MSU_ID'))
  


aggregate <- aggregate %>% rename(MSU_ID = from)

write.csv(aggregate, file = 'GO_hierarchy_drought_response.csv')







drought %>% 
  mutate(cor = ifelse(is.na(cor) == TRUE, 0, cor)) %>% 
  mutate(cor = abs(cor)) %>%  
  group_by(MINUTES, from) %>% 
  summarise(mean = mean(cor), n = n()) %>% 
  filter(n > 2) %>% 
  filter(from %in% c('LOC_Os10g01470', 
                     'LOC_Os02g32590',
                     'LOC_Os07g48596',
                     'LOC_Os04g48070'
  )) %>% 
  left_join(TF_hierarchy, by = c('from'= 'name')) 

# Fig  --------------------------------------------------------------------






library(patchwork)





p2 <- fig_5_3 + inset_element(fig_5_7, left = 0.68, bottom = 0.0, right = 1, top = 0.35)+inset_element(fig_5_6, left = 0, bottom = 0.65, right = 0.31, top = 1.00)


layout <- '
AAAAAA
AAAAAA
BBBBDD
BBBBDD
BBBBEE
BBBBEE
CCCCCC
CCCCCC
CCCCCC
CCCCCC'
ggsave(wrap_plots(A = fig_5_1,
                  B = fig_5_2,   
                  C = p2, 
                  D = fig_5_4, 
                  E = fig_5_5, 
                  design = layout, 
                  guides = 'auto',
                  tag_level = 'new'), filename = 'main_fig_5.pdf', height = 20)








#theme(plot.tag = element_text(size = 8))












# Suppl figures for main fig 5 --------------------------------------------







#with non-responsive TF





fig_5_suppl_1 <- drought %>% 
  select(-cor,-target_value,-to) %>% 
  unique() %>% 
  rowwise() %>% 
  mutate(mean_count = mean(regulator_value$value)) %>% 
  select(-regulator_value) %>% 
  left_join(TF_hierarchy,by = c('from' = 'name')) %>% 
  filter(is.na(level) == FALSE) %>% 
  ungroup() %>% 
  group_by(level, MINUTES) %>% 
  summarise(mean = mean(mean_count), n = n()) %>% 
  mutate(Num = as.factor(paste(n, 'TFs'))) %>% 
  mutate(Level = as.factor(level)) %>% 
  ggplot(aes(MINUTES, mean, color = Level))+ 
  geom_line(size = 1.5)+
  theme_light()+
  theme(legend.title = element_blank(), legend.position = 'top')+
  xlab('Time (min)')+
  ylab('Mean Expression Level')+
  geom_point(size = 2.5)+
  scale_color_manual(values= c('#8FBCBB','#EBCB8B','#81A1C1','#D08770')) +
  theme(
    axis.title.x = element_text(color="black", size=15),
    axis.title.y = element_text(color="black", size=15),
    legend.text = element_text(color = 'black', size = 12)
  )+
  gghighlight::gghighlight(label_key = Num , label_params = list(size = 8))+
  theme(
    axis.text.y=element_blank(),
    axis.text.x = element_text(size = 12),
    axis.ticks.y=element_blank())+
  theme(aspect.ratio=1/1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())






fig_5_suppl_2 <- 
  drought %>% 
  mutate(cor = ifelse(is.na(cor) == TRUE, 0, cor)) %>% 
  mutate(cor = abs(cor)) %>%  
  group_by(MINUTES, from) %>% 
  summarise(mean = mean(cor), n = n()) %>% 
  filter(n > 2) %>% 
  left_join(TF_hierarchy ,by = c('from' = 'name')) %>% 
  filter(is.na(level) == FALSE) %>% 
  ungroup() %>% 
  group_by(MINUTES, level) %>% 
  summarise(mean_activity = mean(mean), n = n()) %>%
  mutate(Num = as.factor(paste(n, 'TFs'))) %>% 
  mutate(Level = as.factor(level)) %>% 
  ggplot(aes(MINUTES, mean_activity, color = Level))+ 
  geom_line(size = 1.5)+
  geom_point(size = 2)+
  theme_light()+
  scale_color_manual(values= c('#8FBCBB','#EBCB8B','#81A1C1','#D08770'), label = c('1 (Bottom)','2','3','4 (Top)')) +
  ylim(0.25,0.5)+
  theme(legend.title = element_blank(), legend.position = 'top')+
  xlab('Time (min)')+
  ylab('TF Activity')+
  geom_point(size = 2)+
  theme(
    axis.title.x = element_text(color="black", size=15),
    axis.title.y = element_text(color="black", size=15),
    legend.text = element_text(color = 'black', size = 20)
  )+
  theme(
    axis.text.y=element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.ticks.y=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())




fig_5_suppl_1 + fig_5_suppl_2 +plot_annotation(tag_levels = 'A')+plot_layout(guides = 'collect')&theme(legend.position = 'bottom')
ggsave(filename = 'main_fig_5_suppl_1.pdf')







# Another way to reconstruct hierarchy ------------------------------------





#remove non-responsive TF









fig_5_suppl_3 <- drought %>% 
  select(-cor,-target_value,-to) %>% 
  unique() %>% 
  rowwise() %>% 
  mutate(mean_count = mean(regulator_value$value)) %>% 
  select(-regulator_value) %>% 
  left_join(
    pair_t %>% mutate(p.value = as.numeric(pair_t))) %>%
  filter(p.value < 0.05) %>%
  left_join(TF_hierarchy,by = c('from' = 'name')) %>% 
  filter(is.na(height) == FALSE) %>% 
  ungroup() %>% 
  group_by(height, MINUTES) %>% 
  summarise(mean = mean(mean_count), n = n()) %>% 
  mutate(Num = as.factor(paste(n, 'TFs'))) %>% 
  mutate(Level = as.factor(height)) %>% 
  ggplot(aes(MINUTES, mean, color = Level))+ 
  geom_line(size = 1.5)+
  theme_light()+
  theme(legend.title = element_blank(), legend.position = 'top')+
  xlab('Time (min)')+
  ylab('Mean Expression Level')+
  geom_point(size = 2.5)+
  scale_color_manual(values= c('#8FBCBB','#EBCB8B','#81A1C1','#D08770')) +
  theme(
    axis.title.x = element_text(color="black", size=15),
    axis.title.y = element_text(color="black", size=15),
    legend.text = element_text(color = 'black', size = 12)
  )+
  gghighlight::gghighlight(label_key = Num , label_params = list(size = 4))+
  theme(
    axis.text.y=element_blank(),
    axis.text.x = element_text(size = 12),
    axis.ticks.y=element_blank())+
  theme(aspect.ratio=1/1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())






fig_5_suppl_4 <- 
  drought %>% 
  mutate(cor = ifelse(is.na(cor) == TRUE, 0, cor)) %>% 
  mutate(cor = abs(cor)) %>%  
  group_by(MINUTES, from) %>% 
  summarise(mean = mean(cor), n = n()) %>% 
  filter(n > 2) %>% 
  left_join(TF_hierarchy ,by = c('from' = 'name')) %>% 
  filter(is.na(height) == FALSE) %>% 
  ungroup() %>% 
  left_join(
    pair_t %>% mutate(p.value = as.numeric(pair_t))) %>%
  filter(p.value < 0.05) %>%
  group_by(MINUTES, height) %>% 
  summarise(mean_activity = mean(mean), n = n()) %>%
  mutate(Num = as.factor(paste(n, 'TFs'))) %>% 
  mutate(Level = as.factor(height)) %>% 
  ggplot(aes(MINUTES, mean_activity, color = Level))+ 
  geom_line(size = 1.5)+
  geom_point(size = 2)+
  theme_light()+
  scale_color_manual(values= c('#8FBCBB','#EBCB8B','#81A1C1','#D08770'), label = c('1 (Bottom)','2','3','4 (Top)')) +
  ylim(0.25,0.6)+
  theme(legend.title = element_blank(), legend.position = 'top')+
  xlab('Time (min)')+
  ylab('TF Activity')+
  geom_point(size = 2)+
  theme(
    axis.title.x = element_text(color="black", size=15),
    axis.title.y = element_text(color="black", size=15),
    legend.text = element_text(color = 'black', size = 20)
  )+
  theme(
    axis.text.y=element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.ticks.y=element_blank())+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank())+
theme(aspect.ratio=1/1)










#with non-responsive TF





fig_5_suppl_5 <- drought %>% 
  select(-cor,-target_value,-to) %>% 
  unique() %>% 
  rowwise() %>% 
  mutate(mean_count = mean(regulator_value$value)) %>% 
  select(-regulator_value) %>% 
  left_join(TF_hierarchy,by = c('from' = 'name')) %>% 
  filter(is.na(height) == FALSE) %>% 
  ungroup() %>% 
  group_by(height, MINUTES) %>% 
  summarise(mean = mean(mean_count), n = n()) %>% 
  mutate(Num = as.factor(paste(n, 'TFs'))) %>% 
  mutate(Level = as.factor(height)) %>% 
  ggplot(aes(MINUTES, mean, color = Level))+ 
  geom_line(size = 1.5)+
  theme_light()+
  theme(legend.title = element_blank(), legend.position = 'top')+
  xlab('Time (min)')+
  ylab('Mean Expression Level')+
  geom_point(size = 2.5)+
  scale_color_manual(values= c('#8FBCBB','#EBCB8B','#81A1C1','#D08770')) +
  theme(
    axis.title.x = element_text(color="black", size=15),
    axis.title.y = element_text(color="black", size=15),
    legend.text = element_text(color = 'black', size = 12)
  )+
  gghighlight::gghighlight(label_key = Num , label_params = list(size = 4))+
  theme(
    axis.text.y=element_blank(),
    axis.text.x = element_text(size = 12),
    axis.ticks.y=element_blank())+
  theme(aspect.ratio=1/1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())






fig_5_suppl_6 <- 
  drought %>% 
  mutate(cor = ifelse(is.na(cor) == TRUE, 0, cor)) %>% 
  mutate(cor = abs(cor)) %>%  
  group_by(MINUTES, from) %>% 
  summarise(mean = mean(cor), n = n()) %>% 
  filter(n > 2) %>% 
  left_join(TF_hierarchy ,by = c('from' = 'name')) %>% 
  filter(is.na(height) == FALSE) %>% 
  ungroup() %>% 
  group_by(MINUTES, height) %>% 
  summarise(mean_activity = mean(mean), n = n()) %>%
  mutate(Num = as.factor(paste(n, 'TFs'))) %>% 
  mutate(Level = as.factor(height)) %>% 
  ggplot(aes(MINUTES, mean_activity, color = Level))+ 
  geom_line(size = 1.5)+
  geom_point(size = 2)+
  theme_light()+
  scale_color_manual(values= c('#8FBCBB','#EBCB8B','#81A1C1','#D08770'), label = c('1 (Bottom)','2','3','4 (Top)')) +
  ylim(0.25,0.6)+
  theme(legend.title = element_blank(), legend.position = 'top')+
  xlab('Time (min)')+
  ylab('TF Activity')+
  geom_point(size = 2)+
  theme(
    axis.title.x = element_text(color="black", size=15),
    axis.title.y = element_text(color="black", size=15),
    legend.text = element_text(color = 'black', size = 20)
  )+
  theme(
    axis.text.y=element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.ticks.y=element_blank())+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank())+
theme(aspect.ratio=1/1)




(fig_5_suppl_3 | fig_5_suppl_4) /(fig_5_suppl_5 | fig_5_suppl_6)+plot_annotation(tag_levels = 'A') +plot_layout(guides = 'collect')&theme(legend.position = 'bottom')
ggsave(filename = 'main_fig_5_suppl_2.pdf')
















