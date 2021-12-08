library(tidyverse)
library(patchwork)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
RNA_seq_summary_Experiment_ID <- read_table2("../Data/RNA-seq summary_Experiment_ID.txt", skip=  6) 
RNA_seq_summary_Experiment_ID <- RNA_seq_summary_Experiment_ID %>% filter(is.na(MINUTES) == FALSE)         
network_prior <- read_delim("../Data/network_prior.txt", 
                            "\t", escape_double = FALSE, col_names = FALSE, 
                            trim_ws = TRUE)



# TEMP robust test (Generate CCF_lag_1 and CCF_lag_2)--------------------------------------------------------------------


max_abs <- function(x) {x[which.max(abs(x))]}
data <- data %>% ungroup() %>%
  filter(CLASS == "CONTROL"  | CLASS == 'DROUGHT'| CLASS == 'HEAT')

data <- data %>%
  select(-key) %>%
  nest(value = c(value))

data <- data %>%
    mutate(value = map(value, ~mean(as.vector(.x$value))))

data <- data %>%
  select(-MINUTES) %>%
  nest(value = c(value))
 
save(data, file = '../Data/rice_tidy_RNA_seq_avg_rep.rdata')



#control 1-9 lag = 1
load('../Data/rice_tidy_RNA_seq_avg_rep.rdata')
control <- data %>% filter(CLASS == "CONTROL") %>% 
mutate(value = map(value, ~ .x[1:9,]))
  

  
df <- RNA_seq_summary_Experiment_ID %>% 
  dplyr::select(CLASS,CULTIVAR) %>% 
  filter(CLASS == 'CONTROL') %>% 
  unique() %>% 
  mutate(prior = map(CLASS, ~network_prior)) %>% 
  unnest(prior) %>% 
  left_join(control, by = c('X1' = 'gene', 'CULTIVAR' = "CULTIVAR",'CLASS' = 'CLASS')) %>% 
  rename(regulator_value = value) %>% 
  left_join(control, by = c('X2' = 'gene',  'CULTIVAR' = "CULTIVAR",'CLASS' = 'CLASS')) %>% 
  rename(target_value = value) #%>% 



df <- df  %>%
    mutate(CCF = map2(regulator_value, target_value,
                      ~ ccf(unlist(.x), unlist(.y), lag = 1, correlation = TRUE,pl = FALSE)))
  
  
  
  CCF <- df %>% dplyr::select(-regulator_value, -target_value) %>%
    cbind(as.data.frame(do.call(rbind, df$CCF))) %>%
    dplyr::select(-CCF)
  
  
  CCF <- CCF %>% select(CLASS,CULTIVAR,X1,X2,acf) %>%
    cbind(as.data.frame(do.call(rbind, CCF$acf))) %>%
    select(-acf)
  
 
  
  
  CCF <- CCF %>% as_tibble() %>%
    filter(is.na(V1) == FALSE) %>%
    rowwise() %>%
    mutate(max_ccf = max_abs(c(V1,V2,V3))) %>%
    mutate_at("max_ccf", as.numeric)
  
  CCF <- CCF %>% as_tibble() %>% 
    group_by(CLASS,X1,X2) %>% 
    summarize(max_ccf = mean(max_ccf))
  

  CCF_control <- CCF


#Drought 1-9 lag = 1
 
 df <- RNA_seq_summary_Experiment_ID %>% 
   dplyr::select(CLASS,CULTIVAR) %>% 
   filter(CLASS == 'DROUGHT') %>% 
   unique() %>% 
   #filter(CLASS == "CONTROL" | CLASS == 'DROUGHT') %>% 
   mutate(prior = map(CLASS, ~network_prior)) %>% 
   unnest(prior) %>% 
   left_join(data, by = c('X1' = 'gene', 'CULTIVAR' = "CULTIVAR",'CLASS' = 'CLASS')) %>% 
   rename(regulator_value = value) %>% 
   left_join(data, by = c('X2' = 'gene',  'CULTIVAR' = "CULTIVAR",'CLASS' = 'CLASS')) %>% 
   rename(target_value = value) #%>% 
 
 
 
 df <- df  %>%
   mutate(CCF = map2(regulator_value, target_value,
                     ~ ccf(unlist(.x), unlist(.y), lag = 1, correlation = TRUE,pl = FALSE)))
 
 
 
 CCF <- df %>% dplyr::select(-regulator_value, -target_value) %>%
   cbind(as.data.frame(do.call(rbind, df$CCF))) %>%
   dplyr::select(-CCF)
 
 
 CCF <- CCF %>% select(CLASS,CULTIVAR,X1,X2,acf) %>%
   cbind(as.data.frame(do.call(rbind, CCF$acf))) %>%
   select(-acf)
 
 
 
 
 CCF <- CCF %>% as_tibble() %>%
   filter(is.na(V1) == FALSE) %>%
   rowwise() %>%
   mutate(max_ccf = max_abs(c(V1,V2,V3))) %>%
   mutate_at("max_ccf", as.numeric)
 
 CCF <- CCF %>% as_tibble() %>% 
   group_by(CLASS,X1,X2) %>% 
   summarize(max_ccf = mean(max_ccf))
 
 
 CCF_drought <- CCF



# heat 1-9 lag = 1


 
 heat <- data %>% filter(CLASS == "HEAT") %>% 
   mutate(value = map(value, ~ .x[1:9,]))
 
 
 df <- RNA_seq_summary_Experiment_ID %>% 
   dplyr::select(CLASS,CULTIVAR) %>% 
   filter(CLASS == 'HEAT') %>% 
   unique() %>% 
   mutate(prior = map(CLASS, ~network_prior)) %>% 
   unnest(prior) %>% 
   left_join(heat, by = c('X1' = 'gene', 'CULTIVAR' = "CULTIVAR",'CLASS' = 'CLASS')) %>% 
   rename(regulator_value = value) %>% 
   left_join(heat, by = c('X2' = 'gene',  'CULTIVAR' = "CULTIVAR",'CLASS' = 'CLASS')) %>% 
   rename(target_value = value) 
 
 
 
 df <- df  %>%
   mutate(CCF = map2(regulator_value, target_value,
                     ~ ccf(unlist(.x), unlist(.y), lag = 1, correlation = TRUE,pl = FALSE)))
 
 
 
 CCF <- df %>% dplyr::select(-regulator_value, -target_value) %>%
   cbind(as.data.frame(do.call(rbind, df$CCF))) %>%
   dplyr::select(-CCF)
 
 
 CCF <- CCF %>% select(CLASS,CULTIVAR,X1,X2,acf) %>%
   cbind(as.data.frame(do.call(rbind, CCF$acf))) %>%
   select(-acf)
 
 
 
 
 CCF <- CCF %>% as_tibble() %>%
   filter(is.na(V1) == FALSE) %>%
   rowwise() %>%
   mutate(max_ccf = max_abs(c(V1,V2,V3))) %>%
   mutate_at("max_ccf", as.numeric)
 
 CCF <- CCF %>% as_tibble() %>% 
   group_by(CLASS,X1,X2) %>% 
   summarize(max_ccf = mean(max_ccf))
 
 
 CCF_heat <- CCF
 
 
 CCF %>%
   mutate_at("max_ccf", as.numeric) %>%
   ggplot(aes(max_ccf,fill = CLASS))+ geom_density(alpha = .35)+
   theme_bw()

 
 
 
 bind_rows(bind_rows(CCF_control,CCF_drought), CCF_heat) %>% 
   mutate_at("max_ccf", as.numeric) %>%
   ggplot(aes(max_ccf,fill = CLASS))+ geom_density(alpha = .35)+
   theme_bw()
   

 
 CCF_lag_1 <- bind_rows(bind_rows(CCF_control,CCF_drought), CCF_heat)

 
 save(CCF_lag_1, file = '../Data/CCF_lag_1.rdata')

 
 #control 1-9 lag = 2
 
 load('../Data/rice_tidy_RNA_seq_avg_rep.rdata')
 
 
 control <- data %>% filter(CLASS == "CONTROL") %>% 
   mutate(value = map(value, ~ .x[1:9,]))
 
 
 
 df <- RNA_seq_summary_Experiment_ID %>% 
   dplyr::select(CLASS,CULTIVAR) %>% 
   filter(CLASS == 'CONTROL') %>% 
   unique() %>% 
   #filter(CLASS == "CONTROL" | CLASS == 'DROUGHT') %>% 
   mutate(prior = map(CLASS, ~network_prior)) %>% 
   unnest(prior) %>% 
   left_join(control, by = c('X1' = 'gene', 'CULTIVAR' = "CULTIVAR",'CLASS' = 'CLASS')) %>% 
   rename(regulator_value = value) %>% 
   left_join(control, by = c('X2' = 'gene',  'CULTIVAR' = "CULTIVAR",'CLASS' = 'CLASS')) %>% 
   rename(target_value = value) #%>% 
 
 
 
 df <- df  %>%
   mutate(CCF = map2(regulator_value, target_value,
                     ~ ccf(unlist(.x), unlist(.y), lag = 2, correlation = TRUE,pl = FALSE)))
 
 
 
 CCF <- df %>% dplyr::select(-regulator_value, -target_value) %>%
   #sample_n(100) %>%
   cbind(as.data.frame(do.call(rbind, df$CCF))) %>%
   dplyr::select(-CCF)
 
 
 CCF <- CCF %>% select(CLASS,CULTIVAR,X1,X2,acf) %>%
   cbind(as.data.frame(do.call(rbind, CCF$acf))) %>%
   select(-acf)
 
 
 
 
 CCF <- CCF %>% as_tibble() %>%
   filter(is.na(V1) == FALSE) %>%
   rowwise() %>%
   mutate(max_ccf = max_abs(c(V1,V2,V3,V4,V5))) %>%
   mutate_at("max_ccf", as.numeric)
 
 CCF <- CCF %>% as_tibble() %>% 
   group_by(CLASS,X1,X2) %>% 
   summarize(max_ccf = mean(max_ccf))
 
 
 CCF_control <- CCF

 
 
# drought lag = 2
 
 
 
 
 
 
 df <- RNA_seq_summary_Experiment_ID %>% 
   dplyr::select(CLASS,CULTIVAR) %>% 
   filter(CLASS == 'DROUGHT') %>% 
   unique() %>% 
   mutate(prior = map(CLASS, ~network_prior)) %>% 
   unnest(prior) %>% 
   left_join(data, by = c('X1' = 'gene', 'CULTIVAR' = "CULTIVAR",'CLASS' = 'CLASS')) %>% 
   rename(regulator_value = value) %>% 
   left_join(data, by = c('X2' = 'gene',  'CULTIVAR' = "CULTIVAR",'CLASS' = 'CLASS')) %>% 
   rename(target_value = value) 
 
 
 
 df <- df  %>%
   mutate(CCF = map2(regulator_value, target_value,
                     ~ ccf(unlist(.x), unlist(.y), lag = 2, correlation = TRUE,pl = FALSE)))
 
 
 
 CCF <- df %>% dplyr::select(-regulator_value, -target_value) %>%
   cbind(as.data.frame(do.call(rbind, df$CCF))) %>%
   dplyr::select(-CCF)
 
 
 CCF <- CCF %>% select(CLASS,CULTIVAR,X1,X2,acf) %>%
   cbind(as.data.frame(do.call(rbind, CCF$acf))) %>%
   select(-acf)
 
 
 CCF <- CCF %>% as_tibble() %>%
   filter(is.na(V1) == FALSE) %>%
   rowwise() %>%
   mutate(max_ccf = max_abs(c(V1,V2,V3,V4,V5))) %>%
   mutate_at("max_ccf", as.numeric)
 
 CCF <- CCF %>% as_tibble() %>% 
   group_by(CLASS,X1,X2) %>% 
   summarize(max_ccf = mean(max_ccf))
 
 
 CCF_drought <- CCF
 

 
 

 #heat 1-9 lag = 2
 
 
 
 heat <- data %>% filter(CLASS == "HEAT") %>% 
   mutate(value = map(value, ~ .x[1:9,]))
 
 
 
 
 df <- RNA_seq_summary_Experiment_ID %>% 
   dplyr::select(CLASS,CULTIVAR) %>% 
   filter(CLASS == 'HEAT') %>% 
   unique() %>% 
   mutate(prior = map(CLASS, ~network_prior)) %>% 
   unnest(prior) %>% 
   left_join(heat, by = c('X1' = 'gene', 'CULTIVAR' = "CULTIVAR",'CLASS' = 'CLASS')) %>% 
   rename(regulator_value = value) %>% 
   left_join(heat, by = c('X2' = 'gene',  'CULTIVAR' = "CULTIVAR",'CLASS' = 'CLASS')) %>% 
   rename(target_value = value) #%>% 
 
 
 
 df <- df  %>%
   mutate(CCF = map2(regulator_value, target_value,
                     ~ ccf(unlist(.x), unlist(.y), lag = 2, correlation = TRUE,pl = FALSE)))
 
 
 
 CCF <- df %>% dplyr::select(-regulator_value, -target_value) %>%
   cbind(as.data.frame(do.call(rbind, df$CCF))) %>%
   dplyr::select(-CCF)
 
 
 CCF <- CCF %>% select(CLASS,CULTIVAR,X1,X2,acf) %>%
   cbind(as.data.frame(do.call(rbind, CCF$acf))) %>%
   select(-acf)
 
 
 
 
 CCF <- CCF %>% as_tibble() %>%
   filter(is.na(V1) == FALSE) %>%
   rowwise() %>%
   mutate(max_ccf = max_abs(c(V1,V2,V3,V4,V5))) %>%
   mutate_at("max_ccf", as.numeric)
 
 CCF <- CCF %>% as_tibble() %>% 
   group_by(CLASS,X1,X2) %>% 
   summarize(max_ccf = mean(max_ccf))
 
 
 CCF_heat <- CCF
 
 

 
 
 
 bind_rows(bind_rows(CCF_control,CCF_drought), CCF_heat) %>% 
   mutate_at("max_ccf", as.numeric) %>%
   ggplot(aes(max_ccf,fill = CLASS))+ geom_density(alpha = .35)+
   theme_bw()
 
 
 
 CCF_lag_2 <- bind_rows(bind_rows(CCF_control,CCF_drought), CCF_heat)
 
 
 save(CCF_lag_2, file = '../Data/CCF_lag_2.rdata')

 
 



# Generate pop_cor for heat and drought -----------------------------------
 load('../Data/rice_tidy_RNA_seq.rdata')
 df_drought <- tibble()
 
 for (i in c(1:9)){
   data_filter <-   data %>% ungroup() %>% 
     filter(CLASS == "CONTROL" | CLASS == 'DROUGHT') %>% 
     filter(MINUTES == i*15) %>% 
     select(-MINUTES,-key, -CULTIVAR) %>% 
     nest(value = c(value))
   
   time = i * 15
   
   df <- RNA_seq_summary_Experiment_ID %>%
     dplyr::select(CLASS) %>%
     unique() %>%
     slice(1,4) %>%
     mutate(prior = map(CLASS, ~network_prior)) %>%
     unnest(prior) %>%
     left_join(data_filter, by = c('X1' = 'gene', 'CLASS' = 'CLASS')) %>%
     rename(regulator_value = value) %>%
     left_join(data_filter, by = c('X2' = 'gene', 'CLASS' = 'CLASS')) %>%
     rename(target_value = value)
   
   
   
   df <- df %>%
     ungroup() %>% 
     rowwise() %>% 
     mutate(cor = cor(regulator_value %>% 
                        pull(),
                      target_value %>% 
                        pull(), use = 'everything',method = 'pearson')) 
   
   
   
   
   df <- df %>% select(-regulator_value, -target_value) %>% 
     mutate(Time = 15*i)
   
   
   df_drought <- bind_rows(df_drought,df)
   
 }
 
 
 save(df_drought, file = '../Data/df_drought_pop.rdata')
 
 
 
 df_heat <- tibble()
 for (i in c(1:9)){
   
   data_filter <-   data %>% ungroup() %>% #sample_n(10000) %>% 
     filter( CLASS == "CONTROL" | CLASS == 'HEAT') %>% 
     filter(MINUTES == i*15) %>% 
     select(-MINUTES,-key, -CULTIVAR) %>% 
     nest(value = c(value))
   
   time = i * 15
   
   df <- RNA_seq_summary_Experiment_ID %>%
     dplyr::select(CLASS) %>%
     unique() %>%
     slice(1,2) %>%
     mutate(prior = map(CLASS, ~network_prior)) %>%
     unnest(prior) %>%
     left_join(data_filter, by = c('X1' = 'gene', 'CLASS' = 'CLASS')) %>%
     rename(regulator_value = value) %>%
     left_join(data_filter, by = c('X2' = 'gene', 'CLASS' = 'CLASS')) %>%
     rename(target_value = value) #%>%
   
   
   
   df <- df %>%
     ungroup() %>% 
     rowwise() %>% 
     mutate(cor = cor(regulator_value %>% 
                        #unnest(regulator_value_wet) %>% 
                        pull(),
                      target_value %>% 
                        #unnest(target_value_wet) %>% 
                        pull(), use = 'everything',method = 'pearson')) 
   
   
   df <- df %>% select(-regulator_value, -target_value) %>% 
     mutate(Time = 15*i)
   
   
   df_heat <- bind_rows(df_heat,df)
   
 }
 
 save(df_heat, file='../Data/df_heat_pop.rdata')
 
# Main_figure_1 -----------------------------------------------------------

 
 
 CCF %>% 
   ggplot(aes(max_ccf, color = Class, fill = Class))+
   theme_bw()+
   theme(legend.title = element_blank(), legend.position = 'top')+
   xlab('MCC')+
   ylab('Density')+
   scale_fill_manual(values=c("#81A1C1", "#A3BE8C",'#B48EAD'),name = ' ',labels = c("Control", 'Drought',"Heat"))+
   scale_color_manual(values=c( "#81A1C1", "#A3BE8C",'#B48EAD'),name = ' ',labels = c("Control",'Drought', "Heat"))+
   theme(
     #plot.title = element_text(color="black", size=5, face="bold"),
     axis.title.x = element_text(color="black", size=12, face="bold"),
     axis.title.y = element_text(color="black", size=12, face="bold"),
     legend.text = element_text(color = 'black', size = 12, face = 'bold')
   )+
   theme(
     axis.text.y=element_blank(),
     axis.ticks.y=element_blank())+
   theme(aspect.ratio=1/1)+
   geom_density(alpha = .4)
 
 
 ggsave(filename = 'main_fig_1.pdf')
 
 
 


# Main_figure_2 -----------------------------------------------------------

 load("../Data/CCF_lag_1.rdata")#load max cross correlation for each pair
 
 CCF <- CCF_lag_1 %>% rename(Class = CLASS) %>% 
   mutate(Class = replace(Class, Class == 'CONTROL', 'Control')) %>% 
   mutate(Class = replace(Class, Class == 'HEAT', 'Heat')) %>% 
   mutate(Class = replace(Class, Class == 'DROUGHT', 'Drought'))
 
 
 CCF<- CCF %>% 
   rename(regulator = X1, target = X2)
 
 
 fun.1  <- function(x) 5*x
 fun.2 <- function(x) x/5 
 
 
 fig_2_1 <- 
   CCF %>% as_tibble() %>% 
   mutate_at("max_ccf", as.numeric) %>% 
   filter(Class == 'Heat' | Class =='Control') %>% 
   group_by(regulator,target,Class) %>% 
   summarise( max_ccf = mean(max_ccf,na.rm = TRUE)) %>% 
   ungroup() %>% 
   mutate(max_ccf = abs(max_ccf)) %>% 
   pivot_wider(c(regulator,target), values_from = max_ccf, names_from = Class) %>% 
   mutate(Category = case_when((abs(Heat) / abs(Control) > 5) ~ 'Heat' ,
                               (abs(Control) / abs(Heat) > 5) ~ 'Control')) %>%
   mutate(Category = ifelse((abs(Control) > 0.69 | abs(Heat) > 0.69), Category, NA)) %>% 
   geom_point(alpha = 0.9, aes(color = Category),size =  0.6)+
   scale_colour_brewer(palette = "Pastel1", na.value = 'grey',
                       labels = c("Regulatory decoherence", "Regulatory coherence",''),
                       guide = guide_legend(override.aes = list(size = 5,alpha = 1) ))+
   geom_function(fun = fun.1)+
   geom_function(fun = fun.2)+
   annotate("label", label = "Heat / Control > 5", x = 0.4, y = 0.6, size = 5) +
   annotate("label", label = "Control / Heat > 5", x = 0.7, y = 0.4, size = 5) +
   annotate(
     geom = "curve", x = 0.3, y = 0.63, xend = 0.15, yend = 0.75, 
     curvature = 0.6, arrow = arrow(length = unit(2, "mm"))
   ) +
   annotate(
     geom = "curve", x = 0.78, y = 0.35, xend = 0.75, yend = 0.15, 
     curvature = 0.5, arrow = arrow(length = unit(2, "mm"))
   ) +
   ylim(0, 1)+
   theme_light()+
   theme(panel.grid.major = element_blank(),
         panel.grid.minor = element_blank())+
   xlab('MCC (Control)')+
   ylab('MCC (Heat)')+
   theme(axis.title = element_text(size = 15),
         axis.text = element_text(size = 12),
         legend.title = element_blank(),
         legend.position = 'top',
         legend.text = element_text(size = 15))+
   theme(legend.background = element_rect(fill="grey",
                                          size=0.5, linetype="solid", 
                                          colour ="white"))+
   theme(aspect.ratio=1/1)
 
 
 fig_2_2 <- 
   CCF %>% as_tibble() %>% 
   mutate_at("max_ccf", as.numeric) %>% 
   filter(Class == 'Drought' | Class =='Control') %>% 
   group_by(regulator,target,Class) %>% 
   summarise( max_ccf = mean(max_ccf,na.rm = TRUE)) %>% 
   ungroup() %>% 
   mutate(max_ccf = abs(max_ccf)) %>% 
   pivot_wider(c(regulator,target), values_from = max_ccf, names_from = Class) %>% 
   mutate(Category = case_when((abs(Drought) / abs(Control) > 5) ~ 'Drought' ,
                               (abs(Control) / abs(Drought) > 5) ~ 'Control')) %>%
   mutate(Category = ifelse((abs(Control) > 0.69 | abs(Drought) > 0.69), Category, NA)) %>% 
   ggplot(aes(Control, Drought))+ 
   geom_point(alpha = 0.9, aes(color = Category),size =  0.6)+
   scale_colour_brewer(palette = "Pastel1", na.value = 'grey',labels = c("Regulatory decoherence", "Regulatory coherence", "NA"))+
   geom_function(fun = fun.1)+
   geom_function(fun = fun.2)+
   annotate("label", label = "Drought / Control > 5", x = 0.4, y = 0.6, size = 5) +
   annotate("label", label = "Control / Drought > 5", x = 0.7, y = 0.4, size = 5) +
   annotate(
     geom = "curve", x = 0.3, y = 0.63, xend = 0.15, yend = 0.75, 
     curvature = 0.6, arrow = arrow(length = unit(2, "mm"))
   ) +
   annotate(
     geom = "curve", x = 0.78, y = 0.35, xend = 0.75, yend = 0.15, 
     curvature = 0.5, arrow = arrow(length = unit(2, "mm"))
   ) +
   ylim(0, 1)+
   theme_light()+
   theme(legend.position = 'none')+
   theme(legend.key.size = unit(0.35, "cm"))+
   theme(panel.grid.major = element_blank(),
         panel.grid.minor = element_blank())+
   xlab('MCC (Control)')+
   ylab('MCC (Drought)')+
   theme(legend.position = 'none',
         legend.title = element_blank(), 
         legend.text = element_text(size = 15)) +
   theme(axis.title = element_text(size = 15),
         axis.text = element_text(size = 12))+
   theme(aspect.ratio=1/1)
 
 
 
 #Main_figure_1_pop
 #drought
 setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
 RNA_seq_summary_Experiment_ID <- read_table2("../Data/RNA-seq summary_Experiment_ID.txt", skip=  6) 
 RNA_seq_summary_Experiment_ID <- RNA_seq_summary_Experiment_ID %>% filter(is.na(MINUTES) == FALSE)         
 network_prior <- read_delim("../Data/network_prior.txt", 
                             "\t", escape_double = FALSE, col_names = FALSE, 
                             trim_ws = TRUE)
 network_prior<- network_prior %>% dplyr::rename(from = X1, to = X2)
 load('../Data/rice_tidy_RNA_seq.rdata')
 
 
 i = 9
 
 data_filter <-   data %>% ungroup() %>% 
   filter(CLASS == "CONTROL" | CLASS == 'DROUGHT') %>% 
   filter(MINUTES == i*15) %>% 
   select(-MINUTES,-key, -CULTIVAR) %>% 
   nest(value = c(value))
 
 time = i * 15
 
 df <- RNA_seq_summary_Experiment_ID %>%
   dplyr::select(CLASS) %>%
   unique() %>%
   dplyr::slice(1,4) %>%
   mutate(prior = map(CLASS, ~network_prior)) %>%
   unnest(prior) %>%
   left_join(data_filter, by = c('from' = 'gene', 'CLASS' = 'CLASS')) %>%
   rename(regulator_value = value) %>%
   left_join(data_filter, by = c('to' = 'gene', 'CLASS' = 'CLASS')) %>%
   rename(target_value = value) 
 
 
 
 df <- df %>%
   ungroup() %>% 
   rowwise() %>% 
   mutate(cor = cor(regulator_value %>% 
                      pull(),
                    target_value %>% 
                      pull(), use = 'everything',method = 'pearson')) 
 
 
 
 
 fun.1  <- function(x) 5*x
 fun.2 <- function(x) x/5 
 fig_2_4 <- df %>% mutate(cor = ifelse(is.na(cor) == FALSE, cor, 0)) %>% 
   pivot_wider(c(from,to), names_from = 'CLASS', values_from = 'cor') %>% 
   mutate(Category = case_when((abs(DROUGHT) / abs(CONTROL) > 5) ~ 'Drought' ,
                               (abs(CONTROL) / abs(DROUGHT) > 5) ~ 'Control')) %>%
   mutate(Category = ifelse((abs(CONTROL) > 0.69 | abs(DROUGHT) > 0.69), Category, NA)) %>%  count(Category)
 ggplot(aes(abs(CONTROL), abs(DROUGHT)))+ 
   geom_point(alpha = 0.9, aes(color = Category),size =  0.6)+
   scale_colour_brewer(palette = "Pastel1", na.value = 'grey',labels = c("Regulatory decoherence", "Regulatory coherence", "NA"))+
   geom_function(fun = fun.1)+
   geom_function(fun = fun.2)+
   annotate("label", label = "Drought / Control > 5", x = 0.4, y = 0.6, size = 5) +
   annotate("label", label = "Control / Drought > 5", x = 0.7, y = 0.4, size = 5) +
   annotate(
     geom = "curve", x = 0.3, y = 0.63, xend = 0.15, yend = 0.75, 
     curvature = 0.6, arrow = arrow(length = unit(2, "mm"))
   ) +
   annotate(
     geom = "curve", x = 0.78, y = 0.35, xend = 0.75, yend = 0.15, 
     curvature = 0.5, arrow = arrow(length = unit(2, "mm"))
   ) +
   ylim(0, 1)+
   theme_light()+
   theme(panel.grid.major = element_blank(),
         panel.grid.minor = element_blank())+
   xlab('PCC (Control) at t = 135 min')+
   ylab('PCC (Drought) at t = 135 min')+
   theme(legend.position = 'none',
         legend.title = element_blank(),
         legend.text = element_text(size = 15)) +
   theme(axis.title = element_text(size = 15),
         axis.text = element_text(size = 12))+
   theme(aspect.ratio=1/1)
 
 
 
 
 
 
 #Heat
 

 i = 9
 data_filter <-   data %>% ungroup() %>% 
   filter(CLASS == "CONTROL" | CLASS == 'HEAT') %>% 
   filter(MINUTES == i*15) %>% 
   select(-MINUTES,-key, -CULTIVAR) %>% 
   nest(value = c(value))
 
 time = i * 15
 
 
 
 df <- RNA_seq_summary_Experiment_ID %>%
   dplyr::select(CLASS) %>%
   unique() %>%
   slice(1,2) %>%
   mutate(prior = map(CLASS, ~network_prior)) %>%
   unnest(prior) %>%
   left_join(data_filter, by = c('from' = 'gene', 'CLASS' = 'CLASS')) %>%
   rename(regulator_value = value) %>%
   left_join(data_filter, by = c('to' = 'gene', 'CLASS' = 'CLASS')) %>%
   rename(target_value = value) #%>%
 
 
 
 df <- df %>%
   ungroup() %>% 
   rowwise() %>% 
   mutate(cor = cor(regulator_value %>% 
                      pull(),
                    target_value %>% 
                      pull(), use = 'everything',method = 'pearson')) 
 
 fun.1  <- function(x) 5*x
 fun.2 <- function(x) x/5 
 fig_2_3 <- df %>% mutate(cor = ifelse(is.na(cor) == FALSE, cor, 0)) %>% 
   pivot_wider(c(from,to), names_from = 'CLASS', values_from = 'cor') %>% 
   mutate(Category = case_when((abs(HEAT) / abs(CONTROL) > 5) ~ 'Heat' ,
                               (abs(CONTROL) / abs(HEAT) > 5) ~ 'Control')) %>%
   mutate(Category = ifelse((abs(CONTROL) > 0.69 | abs(HEAT) > 0.69), Category, NA)) %>% 
   ggplot(aes(abs(CONTROL), abs(HEAT)))+ #geom_point(alpha = 0.3,color = "#999999", size = 2.5)+
   geom_point(alpha = 0.9, aes(color = Category),size =  0.6)+
   scale_colour_brewer(palette = "Pastel1", na.value = 'grey',labels = c("Regulatory decoherence", "Regulatory coherence", "NA"))+
   geom_function(fun = fun.1)+
   geom_function(fun = fun.2)+
   annotate("label", label = "Heat / Control > 5", x = 0.4, y = 0.6, size = 5) +
   annotate("label", label = "Control / Heat > 5", x = 0.7, y = 0.4, size = 5) +
   annotate(
     geom = "curve", x = 0.3, y = 0.63, xend = 0.15, yend = 0.75, 
     curvature = 0.6, arrow = arrow(length = unit(2, "mm"))
   ) +
   annotate(
     geom = "curve", x = 0.78, y = 0.35, xend = 0.75, yend = 0.15, 
     curvature = 0.5, arrow = arrow(length = unit(2, "mm"))
   ) +
   ylim(0, 1)+
   theme_light()+
   theme(panel.grid.major = element_blank(),
         panel.grid.minor = element_blank())+
   xlab('PCC (Control) at t = 135 min')+
   ylab('PCC (Heat) at t = 135 min')+
   theme(axis.title = element_text(size = 15),
         axis.text = element_text(size = 12))+
   theme(legend.position = 'none',
         legend.title = element_blank(),
         legend.text = element_text(size = 12))+
   theme(aspect.ratio=1/1)
 
 
 main_fig_2 <- (fig_2_1 + fig_2_2) / guide_area()/(fig_2_3 + fig_2_4)+ plot_annotation(tag_levels = 'A',)+plot_layout(guides = 'collect', heights = c(8,1,8))& theme(plot.tag = element_text(size = 16))
 ggsave(main_fig_1, filename = 'main_fig_2.pdf', height = 12, width = 12)
 

# Suppl_figure_for_main_figure 2 ------------------------------------------
 
 #main_fig_2_suppl_1
 CCF %>% filter(Class == "Control") %>% 
   ungroup() %>% 
   mutate(max_ccf = abs(max_ccf)) %>% 
   summarise(x = quantile(max_ccf, c(0.5, 0.90, 0.95, 0.99, 0.999)), q = c(0.5, 0.90, 0.95, 0.99,0.999))
 

 
 fig_2_suppl_1_1 <- CCF %>% as_tibble() %>% 
   mutate_at("max_ccf", as.numeric) %>% 
   filter(Class == 'Heat' | Class =='Control') %>% 
   group_by(regulator,target,Class) %>% 
   summarise( max_ccf = mean(max_ccf,na.rm = TRUE)) %>% 
   ungroup() %>% 
   pivot_wider(c(regulator,target), values_from = max_ccf, names_from = Class) %>% 
   mutate(Category = case_when((abs(Heat) > 0.566& abs(Control) < 0.566) ~ 'Heat' ,
                               (abs(Heat) < 0.566& abs(Control) > 0.566) ~ 'Control')) %>%
   ggplot(aes(Control, Heat))+ 
   geom_point(alpha = 0.9, aes(color = Category),size = 0.6)+
   scale_colour_brewer(palette = "Pastel1", na.value = 'grey',labels = c("Regulatory decoherence", "Regulatory coherence", "NA"))+
   theme_bw()+
   annotate("label", label = "Regulatory coherence", x = 0.1, y = 1.1, size = 4) +
   annotate("label", label = "Regulatory decoherence", x = 0.1, y = -1.05, size = 4) +
   annotate(
     geom = "curve", x = 0.1, y = 1.0, xend = 0.15, yend = 0.75, 
     curvature = -0.5, arrow = arrow(length = unit(2, "mm"))
   ) +
   annotate(
     geom = "curve", x = 0.1, y = -0.92, xend = 0.70, yend = -0.35, 
     curvature = 0.5, arrow = arrow(length = unit(2, "mm"))
   ) +
   theme(legend.position  = 'none')+
   theme(panel.grid.major = element_blank(),
         panel.grid.minor = element_blank())+
   xlab('MCC (Control)')+
   ylab('MCC(Heat)')+
   theme(#legend.position = c(0.88, 0.1),
     legend.title = element_blank(),
     legend.text = element_text(size = 10)) +
   theme(axis.title = element_text(size = 15),
         axis.text = element_text(size = 7))+
   theme(aspect.ratio=1/1)
 
 fun.1  <- function(x) 5*x
 fun.2 <- function(x) x/5 
 fig_2_suppl_1_2 <- 
   CCF %>% as_tibble() %>% 
   mutate_at("max_ccf", as.numeric) %>% 
   filter(Class == 'Heat' | Class =='Control') %>% 
   group_by(regulator,target,Class) %>% 
   summarise( max_ccf = mean(max_ccf,na.rm = TRUE)) %>% 
   ungroup() %>% 
   mutate(max_ccf = abs(max_ccf)) %>% 
   pivot_wider(c(regulator,target), values_from = max_ccf, names_from = Class) %>% 
   mutate(Category = case_when((abs(Heat) / abs(Control) > 5) ~ 'Heat' ,
                               (abs(Control) / abs(Heat) > 5) ~ 'Control')) %>%
   mutate(Category = ifelse((abs(Control) > 0.566 | abs(Heat) > 0.566), Category, NA)) %>% 
   ggplot(aes(Control, Heat))+ 
   geom_point(alpha = 0.9, aes(color = Category),size =  0.6)+
   scale_colour_brewer(palette = "Pastel1", na.value = 'grey',labels = c("Regulatory decoherence", "Regulatory coherence", "NA"))+
   geom_function(fun = fun.1)+
   geom_function(fun = fun.2)+
   annotate("label", label = "Heat / Control > 5", x = 0.4, y = 0.6, size = 4) +
   annotate("label", label = "Control / Heat > 5", x = 0.7, y = 0.4, size = 4) +
   annotate(
     geom = "curve", x = 0.3, y = 0.63, xend = 0.15, yend = 0.75, 
     curvature = 0.6, arrow = arrow(length = unit(2, "mm"))
   ) +
   annotate(
     geom = "curve", x = 0.78, y = 0.35, xend = 0.75, yend = 0.15, 
     curvature = 0.5, arrow = arrow(length = unit(2, "mm"))
   ) +
   ylim(0, 1)+
   theme_bw()+
   theme(panel.grid.major = element_blank(),
         panel.grid.minor = element_blank())+
   xlab('|MCC (Control)|')+
   ylab('|MCC (Heat)|')+
   theme(legend.position = c(0.70, 0.9),
         legend.title = element_blank(),
         legend.text = element_text(size = 10)) +
   theme(axis.title = element_text(size = 15),
         axis.text = element_text(size = 7))+
   theme(legend.position = 'none')+
   theme(aspect.ratio=1/1)

 fig_2_suppl_1_3 <- CCF %>% as_tibble() %>% 
   mutate_at("max_ccf", as.numeric) %>% 
   filter(Class == 'Drought' | Class =='Control') %>% 
   group_by(regulator,target,Class) %>% 
   summarise( max_ccf = mean(max_ccf,na.rm = TRUE)) %>% 
   ungroup() %>% 
   pivot_wider(c(regulator,target), values_from = max_ccf, names_from = Class) %>% 
   mutate(Category = case_when((abs(Drought) > 0.566& abs(Control) < 0.566) ~ 'Drought' ,
                               (abs(Drought) < 0.566& abs(Control) > 0.566) ~ 'Control')) %>%
   ggplot(aes(Control, Drought))+
   geom_point(alpha = 0.9, aes(color = Category),size =  0.6)+
   scale_colour_brewer(palette = "Pastel1", na.value = 'grey',labels = c("Regulatory decoherence", "Regulatory coherence", "NA"))+
   theme_bw()+
   annotate("label", label = "Regulatory coherence", x = 0.3, y = 1.1, size = 4) +
   annotate("label", label = "Regulatory decoherence", x = 0.3, y = -1.05, size = 4) +
   annotate(
     geom = "curve", x = 0.3, y = 1.0, xend = 0.15, yend = 0.75, 
     curvature = -0.5, arrow = arrow(length = unit(2, "mm"))
   ) +
   annotate(
     geom = "curve", x = 0.3, y = -0.92, xend = 0.70, yend = -0.35, 
     curvature = 0.5, arrow = arrow(length = unit(2, "mm"))
   ) +
   theme(legend.position  = 'none')+
   theme(panel.grid.major = element_blank(),
         panel.grid.minor = element_blank())+
   xlab('MCC (Control)')+
   ylab('MCC (Drought)')+
   theme(
     legend.title = element_blank(),
     legend.text = element_text(size = 10)) +
   theme(axis.title = element_text(size = 15),
         axis.text = element_text(size = 7))+
   theme(aspect.ratio=1/1)

 
 

 fun.1  <- function(x) 5*x
 fun.2 <- function(x) x/5 
 fig_2_suppl_1_4 <- 
   CCF %>% as_tibble() %>% 
   mutate_at("max_ccf", as.numeric) %>% 
   filter(Class == 'Drought' | Class =='Control') %>% 
   group_by(regulator,target,Class) %>% 
   summarise( max_ccf = mean(max_ccf,na.rm = TRUE)) %>% 
   ungroup() %>% 
   mutate(max_ccf = abs(max_ccf)) %>% 
   pivot_wider(c(regulator,target), values_from = max_ccf, names_from = Class) %>% 
   mutate(Category = case_when((abs(Drought) / abs(Control) > 5) ~ 'Heat' ,
                               (abs(Control) / abs(Drought) > 5) ~ 'Control')) %>%
   mutate(Category = ifelse((abs(Control) > 0.566 | abs(Drought) > 0.566), Category, NA)) %>% 
   #filter(is.na(Category) == FALSE)
   #filter(abs(HEAT) < 0.7| abs(CONTROL) < 0.7) %>%
   #mutate(Category = ifelse(is.na(Category) == TRUE, 'NA',Category)) %>% 
   #filter(Category == 'Control' | Category == 'Heat') %>%
   ggplot(aes(Control, Drought))+ #geom_point(alpha = 0.3,color = "#999999", size = 2.5)+
   geom_point(alpha = 0.9, aes(color = Category),size =  0.6)+
   #stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white")
   scale_colour_brewer(palette = "Pastel1", na.value = 'grey',labels = c("Regulatory decoherence", "Regulatory coherence", "NA"))+
   geom_function(fun = fun.1)+
   geom_function(fun = fun.2)+
   annotate("label", label = "Drought / Control > 5", x = 0.4, y = 0.6, size = 4) +
   annotate("label", label = "Control / Drought > 5", x = 0.7, y = 0.4, size = 4) +
   annotate(
     geom = "curve", x = 0.3, y = 0.63, xend = 0.15, yend = 0.75, 
     curvature = 0.6, arrow = arrow(length = unit(2, "mm"))
   ) +
   annotate(
     geom = "curve", x = 0.78, y = 0.35, xend = 0.75, yend = 0.15, 
     curvature = 0.5, arrow = arrow(length = unit(2, "mm"))
   ) +
   ylim(0, 1)+
   theme_bw()+
   theme(legend.position = 'none')+
   theme(legend.key.size = unit(0.35, "cm"))+
   theme(panel.grid.major = element_blank(),
         panel.grid.minor = element_blank())+
   xlab('|MCC (Control)|')+
   ylab('|MCC (Drought)|')+
   theme(legend.position = c(0.70, 0.9),
         legend.title = element_blank(),
         legend.text = element_text(size = 10)) +
   theme(axis.title = element_text(size = 15),
         axis.text = element_text(size = 7))+
   theme(legend.position = 'none')+
   theme(aspect.ratio=1/1)
 
 
main_fig_2_suppl_1 <- (fig_2_suppl_1_1 + fig_2_suppl_1_3) / (fig_2_suppl_1_2 + fig_2_suppl_1_4)+ plot_annotation(tag_levels = 'A')
ggsave(main_fig_2_suppl_1, filename = 'main_fig_2_suppl_1.pdf')
 
 
 
#main_fig_2_suppl_2
 
fig_2_suppl_2_1 <- CCF %>% as_tibble() %>% 
   mutate_at("max_ccf", as.numeric) %>% 
   filter(Class == 'Heat' | Class =='Control') %>% 
   group_by(regulator,target,Class) %>% 
   summarise( max_ccf = mean(max_ccf,na.rm = TRUE)) %>% 
   ungroup() %>% 
   pivot_wider(c(regulator,target), values_from = max_ccf, names_from = Class) %>% 
   mutate(Category = case_when((abs(Heat) > 0.867& abs(Control) < 0.867) ~ 'Heat' ,
                               (abs(Heat) < 0.867& abs(Control) > 0.867) ~ 'Control')) %>%
   geom_point(alpha = 0.9, aes(color = Category),size = 0.6)+
   scale_colour_brewer(palette = "Pastel1", na.value = 'grey',labels = c("Regulatory decoherence", "Regulatory coherence", "NA"))+
   theme_bw()+
   annotate("label", label = "Regulatory coherence", x = 0.1, y = 1.1, size = 4) +
   annotate("label", label = "Regulatory decoherence", x = 0.1, y = -1.05, size = 4) +
   annotate(
     geom = "curve", x = 0.1, y = 1.0, xend = 0.15, yend = 0.85, 
     curvature = -0.5, arrow = arrow(length = unit(2, "mm"))
   ) +
   annotate(
     geom = "curve", x = 0.1, y = -0.92, xend = 0.820, yend = -0.35, 
     curvature = 0.5, arrow = arrow(length = unit(2, "mm"))
   ) +
   theme(legend.position  = 'none')+
   theme(panel.grid.major = element_blank(),
         panel.grid.minor = element_blank())+
   xlab('MCC (Control)')+
   ylab('MCC (Heat)')+
   theme(
     legend.title = element_blank(),
     legend.text = element_text(size = 10)) +
   theme(axis.title = element_text(size = 15),
         axis.text = element_text(size = 7))+
   theme(aspect.ratio=1/1)

 
 fun.1  <- function(x) 5*x
 fun.2 <- function(x) x/5 
 fig_2_suppl_2_2 <- 
   CCF %>% as_tibble() %>% 
   mutate_at("max_ccf", as.numeric) %>% 
   filter(Class == 'Heat' | Class =='Control') %>% 
   group_by(regulator,target,Class) %>% 
   summarise( max_ccf = mean(max_ccf,na.rm = TRUE)) %>% 
   ungroup() %>% 
   mutate(max_ccf = abs(max_ccf)) %>% 
   pivot_wider(c(regulator,target), values_from = max_ccf, names_from = Class) %>% 
   mutate(Category = case_when((abs(Heat) / abs(Control) > 5) ~ 'Heat' ,
                               (abs(Control) / abs(Heat) > 5) ~ 'Control')) %>%
   mutate(Category = ifelse((abs(Control) > 0.867 | abs(Heat) > 0.867), Category, NA)) %>% 
   ggplot(aes(Control, Heat))+ 
   geom_point(alpha = 0.9, aes(color = Category),size =  0.6)+
   scale_colour_brewer(palette = "Pastel1", na.value = 'grey',labels = c("Regulatory decoherence", "Regulatory coherence", "NA"))+
   geom_function(fun = fun.1)+
   geom_function(fun = fun.2)+
   annotate("label", label = "Heat / Control > 5", x = 0.4, y = 0.6, size = 4) +
   annotate("label", label = "Control / Heat > 5", x = 0.7, y = 0.4, size = 4) +
   annotate(
     geom = "curve", x = 0.3, y = 0.63, xend = 0.15, yend = 0.75, 
     curvature = 0.6, arrow = arrow(length = unit(2, "mm"))
   ) +
   annotate(
     geom = "curve", x = 0.78, y = 0.35, xend = 0.75, yend = 0.15, 
     curvature = 0.5, arrow = arrow(length = unit(2, "mm"))
   ) +
   ylim(0, 1)+
   theme_bw()+
   theme(panel.grid.major = element_blank(),
         panel.grid.minor = element_blank())+
   xlab('|MCC (Control)|')+
   ylab('|MCC (Heat)|')+
   theme(legend.position = c(0.70, 0.9),
         legend.title = element_blank(),
         legend.text = element_text(size = 10)) +
   theme(axis.title = element_text(size = 15),
         axis.text = element_text(size = 7))+
   theme(legend.position = 'none')+
   theme(aspect.ratio=1/1)
 

 fig_2_suppl_2_3 <- CCF %>% as_tibble() %>% 
   mutate_at("max_ccf", as.numeric) %>% 
   filter(Class == 'Drought' | Class =='Control') %>% 
   group_by(regulator,target,Class) %>% 
   summarise( max_ccf = mean(max_ccf,na.rm = TRUE)) %>% 
   ungroup() %>% 
   pivot_wider(c(regulator,target), values_from = max_ccf, names_from = Class) %>% 
   mutate(Category = case_when((abs(Drought) > 0.867& abs(Control) < 0.867) ~ 'Drought' ,
                               (abs(Drought) < 0.867& abs(Control) > 0.867) ~ 'Control')) %>%
   geom_point(alpha = 0.9, aes(color = Category),size =  0.6)+
   scale_colour_brewer(palette = "Pastel1", na.value = 'grey',labels = c("Regulatory decoherence", "Regulatory coherence", "NA"))+
   theme_bw()+
   annotate("label", label = "Regulatory coherence", x = 0.3, y = 1.1, size = 4) +
   annotate("label", label = "Regulatory decoherence", x = 0.3, y = -1.05, size = 4) +
   annotate(
     geom = "curve", x = 0.3, y = 1.0, xend = 0.15, yend = 0.85, 
     curvature = -0.5, arrow = arrow(length = unit(2, "mm"))
   ) +
   annotate(
     geom = "curve", x = 0.3, y = -0.92, xend = 0.820, yend = -0.35, 
     curvature = 0.5, arrow = arrow(length = unit(2, "mm"))
   ) +
   theme(legend.position  = 'none')+
   theme(panel.grid.major = element_blank(),
         panel.grid.minor = element_blank())+
   xlab('MCC (Control)')+
   ylab('MCC (Drought)')+
   theme(
     legend.title = element_blank(),
     legend.text = element_text(size = 10)) +
   theme(axis.title = element_text(size = 15),
         axis.text = element_text(size = 7))+
   theme(aspect.ratio=1/1)

 
 
 
 
 
 fun.1  <- function(x) 5*x
 fun.2 <- function(x) x/5 
 fig_2_suppl_2_4 <- 
   CCF %>% as_tibble() %>% 
   mutate_at("max_ccf", as.numeric) %>% 
   filter(Class == 'Drought' | Class =='Control') %>% 
   group_by(regulator,target,Class) %>% 
   summarise( max_ccf = mean(max_ccf,na.rm = TRUE)) %>% 
   ungroup() %>% 
   mutate(max_ccf = abs(max_ccf)) %>% 
   pivot_wider(c(regulator,target), values_from = max_ccf, names_from = Class) %>% 
   mutate(Category = case_when((abs(Drought) / abs(Control) > 5) ~ 'Heat' ,
                               (abs(Control) / abs(Drought) > 5) ~ 'Control')) %>%
   mutate(Category = ifelse((abs(Control) > 0.867 | abs(Drought) > 0.867), Category, NA)) %>% 
   ggplot(aes(Control, Drought))+
   geom_point(alpha = 0.9, aes(color = Category),size =  0.6)+
   scale_colour_brewer(palette = "Pastel1", na.value = 'grey',labels = c("Regulatory decoherence", "Regulatory coherence", "NA"))+
   geom_function(fun = fun.1)+
   geom_function(fun = fun.2)+
   annotate("label", label = "Drought / Control > 5", x = 0.4, y = 0.6, size = 4) +
   annotate("label", label = "Control / Drought > 5", x = 0.7, y = 0.4, size = 4) +
   annotate(
     geom = "curve", x = 0.3, y = 0.63, xend = 0.15, yend = 0.75, 
     curvature = 0.6, arrow = arrow(length = unit(2, "mm"))
   ) +
   annotate(
     geom = "curve", x = 0.78, y = 0.35, xend = 0.75, yend = 0.15, 
     curvature = 0.5, arrow = arrow(length = unit(2, "mm"))
   ) +
   ylim(0, 1)+
   theme_bw()+
   theme(legend.position = 'none')+
   theme(legend.key.size = unit(0.35, "cm"))+
   theme(panel.grid.major = element_blank(),
         panel.grid.minor = element_blank())+
   xlab('|MCC (Control)|')+
   ylab('|MCC (Drought)|')+
   theme(legend.position = c(0.70, 0.9),
         legend.title = element_blank(),
         legend.text = element_text(size = 10)) +
   theme(axis.title = element_text(size = 15),
         axis.text = element_text(size = 7))+
   theme(legend.position = 'none')+
   theme(aspect.ratio=1/1)
 
 
 fig_2_suppl_2 <- (fig_2_suppl_2_1 + fig_2_suppl_2_3) / (fig_2_suppl_2_2 + fig_2_suppl_2_4)+ plot_annotation(tag_levels = 'A')
 ggsave(main_fig_2_suppl_2, filename = 'main_fig_2_suppl_2.pdf') 
 
 
 
 
 
 
 # #main_fig_2_suppl_3 with lag time = 2
 
 load("../Data/CCF_lag_2.rdata") #load max cross correlation for each pair with lag = 2
 
 CCF <- CCF_lag_2 %>% rename(Class = CLASS) %>% 
   mutate(Class = replace(Class, Class == 'CONTROL', 'Control')) %>% 
   mutate(Class = replace(Class, Class == 'HEAT', 'Heat')) %>% 
   mutate(Class = replace(Class, Class == 'DROUGHT', 'Drought'))
 
 
 CCF<- CCF %>% 
   rename(regulator = X1, target = X2)
 
 
 fig_2_suppl_3_1 <- CCF %>% as_tibble() %>% 
   mutate_at("max_ccf", as.numeric) %>% 
   filter(Class == 'Heat' | Class =='Control') %>% 
   group_by(regulator,target,Class) %>% 
   summarise( max_ccf = mean(max_ccf,na.rm = TRUE)) %>% 
   ungroup() %>% 
   pivot_wider(c(regulator,target), values_from = max_ccf, names_from = Class) %>% 
   mutate(Category = case_when((abs(Heat) > 0.69& abs(Control) < 0.69) ~ 'Heat' ,
                               (abs(Heat) < 0.69& abs(Control) > 0.69) ~ 'Control')) %>%
   ggplot(aes(Control, Heat))+
   geom_point(alpha = 0.9, aes(color = Category),size = 0.6)+
   scale_colour_brewer(palette = "Pastel1", na.value = 'grey',labels = c("Regulatory decoherence", "Regulatory coherence", "NA"))+
   theme_bw()+
   annotate("label", label = "Regulatory coherence", x = 0.1, y = 1.1, size = 4) +
   annotate("label", label = "Regulatory decoherence", x = 0.1, y = -1.05, size = 4) +
   annotate(
     geom = "curve", x = 0.1, y = 1.0, xend = 0.15, yend = 0.85, 
     curvature = -0.5, arrow = arrow(length = unit(2, "mm"))
   ) +
   annotate(
     geom = "curve", x = 0.1, y = -0.92, xend = 0.820, yend = -0.35, 
     curvature = 0.5, arrow = arrow(length = unit(2, "mm"))
   ) +
   theme(legend.position  = 'none')+
   theme(panel.grid.major = element_blank(),
         panel.grid.minor = element_blank())+
   xlab('MCC (Control)')+
   ylab('MCC (Heat)')+
   theme(
     legend.title = element_blank(),
     legend.text = element_text(size = 10)) +
   theme(axis.title = element_text(size = 15),
         axis.text = element_text(size = 7))+
   theme(aspect.ratio=1/1)

 
 fun.1  <- function(x) 5*x
 fun.2 <- function(x) x/5 
 fig_2_suppl_3_2 <- 
   CCF %>% as_tibble() %>% 
   mutate_at("max_ccf", as.numeric) %>% 
   filter(Class == 'Heat' | Class =='Control') %>% 
   group_by(regulator,target,Class) %>% 
   summarise( max_ccf = mean(max_ccf,na.rm = TRUE)) %>% 
   ungroup() %>% 
   mutate(max_ccf = abs(max_ccf)) %>% 
   pivot_wider(c(regulator,target), values_from = max_ccf, names_from = Class) %>% 
   mutate(Category = case_when((abs(Heat) / abs(Control) > 5) ~ 'Heat' ,
                               (abs(Control) / abs(Heat) > 5) ~ 'Control')) %>%
   mutate(Category = ifelse((abs(Control) > 0.69 | abs(Heat) > 0.69), Category, NA)) %>% 
   ggplot(aes(Control, Heat))+ 
   geom_point(alpha = 0.9, aes(color = Category),size =  0.6)+
   scale_colour_brewer(palette = "Pastel1", na.value = 'grey',labels = c("Regulatory decoherence", "Regulatory coherence", "NA"))+
   geom_function(fun = fun.1)+
   geom_function(fun = fun.2)+
   annotate("label", label = "Heat / Control > 5", x = 0.4, y = 0.6, size = 4) +
   annotate("label", label = "Control / Heat > 5", x = 0.7, y = 0.4, size = 4) +
   annotate(
     geom = "curve", x = 0.3, y = 0.63, xend = 0.15, yend = 0.75, 
     curvature = 0.6, arrow = arrow(length = unit(2, "mm"))
   ) +
   annotate(
     geom = "curve", x = 0.78, y = 0.35, xend = 0.75, yend = 0.15, 
     curvature = 0.5, arrow = arrow(length = unit(2, "mm"))
   ) +
   ylim(0, 1)+
   theme_bw()+
   theme(panel.grid.major = element_blank(),
         panel.grid.minor = element_blank())+
   xlab('|MCC (Control)|')+
   ylab('|MCC (Heat)|')+
   theme(legend.position = c(0.70, 0.9),
         legend.title = element_blank(),
         legend.text = element_text(size = 10)) +
   theme(axis.title = element_text(size = 15),
         axis.text = element_text(size = 7))+
   theme(legend.position = 'none')+
   theme(aspect.ratio=1/1)
 

 fig_2_suppl_3_3 <- CCF %>% as_tibble() %>% 
   mutate_at("max_ccf", as.numeric) %>% 
   filter(Class == 'Drought' | Class =='Control') %>% 
   group_by(regulator,target,Class) %>% 
   summarise( max_ccf = mean(max_ccf,na.rm = TRUE)) %>% 
   ungroup() %>% 
   pivot_wider(c(regulator,target), values_from = max_ccf, names_from = Class) %>% 
   mutate(Category = case_when((abs(Drought) > 0.69& abs(Control) < 0.69) ~ 'Drought' ,
                               (abs(Drought) < 0.69& abs(Control) > 0.69) ~ 'Control')) %>%
   ggplot(aes(Control, Drought))+ 
   geom_point(alpha = 0.9, aes(color = Category),size =  0.6)+
   scale_colour_brewer(palette = "Pastel1", na.value = 'grey',labels = c("Regulatory decoherence", "Regulatory coherence", "NA"))+
   theme_bw()+
   annotate("label", label = "Regulatory coherence", x = 0.3, y = 1.1, size = 4) +
   annotate("label", label = "Regulatory decoherence", x = 0.3, y = -1.05, size = 4) +
   annotate(
     geom = "curve", x = 0.3, y = 1.0, xend = 0.15, yend = 0.85, 
     curvature = -0.5, arrow = arrow(length = unit(2, "mm"))
   ) +
   annotate(
     geom = "curve", x = 0.3, y = -0.92, xend = 0.820, yend = -0.35, 
     curvature = 0.5, arrow = arrow(length = unit(2, "mm"))
   ) +
   theme(legend.position  = 'none')+
   theme(panel.grid.major = element_blank(),
         panel.grid.minor = element_blank())+
   xlab('MCC (Control)')+
   ylab('MCC (Drought)')+
   theme(
     legend.title = element_blank(),
     legend.text = element_text(size = 10)) +
   theme(axis.title = element_text(size = 15),
         axis.text = element_text(size = 7))+
   theme(aspect.ratio=1/1)

 
 
 
 
 
 
 fun.1  <- function(x) 5*x
 fun.2 <- function(x) x/5 
 fig_2_suppl_3_4 <- 
   CCF %>% as_tibble() %>% 
   mutate_at("max_ccf", as.numeric) %>% 
   filter(Class == 'Drought' | Class =='Control') %>% 
   group_by(regulator,target,Class) %>% 
   summarise( max_ccf = mean(max_ccf,na.rm = TRUE)) %>% 
   ungroup() %>% 
   mutate(max_ccf = abs(max_ccf)) %>% 
   pivot_wider(c(regulator,target), values_from = max_ccf, names_from = Class) %>% 
   mutate(Category = case_when((abs(Drought) / abs(Control) > 5) ~ 'Heat' ,
                               (abs(Control) / abs(Drought) > 5) ~ 'Control')) %>%
   mutate(Category = ifelse((abs(Control) > 0.69 | abs(Drought) > 0.69), Category, NA)) %>% 
   ggplot(aes(Control, Drought))+ 
   geom_point(alpha = 0.9, aes(color = Category),size =  0.6)+
   scale_colour_brewer(palette = "Pastel1", na.value = 'grey',labels = c("Regulatory decoherence", "Regulatory coherence", "NA"))+
   geom_function(fun = fun.2)+
   annotate("label", label = "Drought / Control > 5", x = 0.4, y = 0.6, size = 4) +
   annotate("label", label = "Control / Drought > 5", x = 0.7, y = 0.4, size = 4) +
   annotate(
     geom = "curve", x = 0.3, y = 0.63, xend = 0.15, yend = 0.75, 
     curvature = 0.6, arrow = arrow(length = unit(2, "mm"))
   ) +
   annotate(
     geom = "curve", x = 0.78, y = 0.35, xend = 0.75, yend = 0.15, 
     curvature = 0.5, arrow = arrow(length = unit(2, "mm"))
   ) +
   ylim(0, 1)+
   theme_bw()+
   theme(legend.position = 'none')+
   theme(legend.key.size = unit(0.35, "cm"))+
   theme(panel.grid.major = element_blank(),
         panel.grid.minor = element_blank())+
   xlab('|MCC (Control)|')+
   ylab('|MCC (Drought)|')+
   theme(legend.position = c(0.70, 0.9),
         legend.title = element_blank(),
         legend.text = element_text(size = 10)) +
   theme(axis.title = element_text(size = 15),
         axis.text = element_text(size = 7))+
   theme(legend.position = 'none')+
   theme(aspect.ratio=1/1)
 
 
 main_fig_2_suppl_3 <- (fig_2_suppl_3_1 + fig_2_suppl_3_3) / (fig_2_suppl_3_2 + fig_2_suppl_3_4)+ plot_annotation(tag_levels = 'A')
 ggsave(main_fig_2_suppl_3, filename = 'main_fig_2_suppl_3.pdf') 
 
 
 
 
 
 # #main_fig_2_suppl_4 (brachy)
 library(DESeq2)
 library(statmod)
 library(corrr)
 library(gdata)
 normalized_counts_with_chloroplast <- read_csv("../Data/normalized_counts_with_chloroplast.csv")
 samplelabels <- read_csv("../Data/samplelabels.csv")
 samplelabels <- samplelabels %>% select(-X1)
 
 counts <- normalized_counts_with_chloroplast %>% select(samplelabels %>% select(`rownames(matrixbd)`) %>% 
                                                           pull())
 
 counts <- counts %>% mutate(across(everything(), round, 0))
 rownames(counts) <- normalized_counts_with_chloroplast$X1
 # counts_ave5 <-  filter(counts, rowMeans(counts[,1:192])>=5)
 # rownames(counts_ave5) <- normalized_counts_with_chloroplast$X1
 # 
 # 
 # 
 # dds <- DESeqDataSetFromMatrix(countData = counts_ave5, #counts
 #                               colData = samplelabels %>% select(Accession, Treatment),
 #                               design= ~ Treatment)
 # dds <- DESeq(dds)
 # res <- results(dds)
 
 

 counts_bd21 <- counts[,samplelabels %>% filter(Accession == 'BD21') %>% pull('rownames(matrixbd)')]
 rownames(counts_bd21) <- normalized_counts_with_chloroplast$X1
 counts_bd21 <-  filter(counts_bd21, rowMeans(counts[,1:192])>=5)
 dds_bd21 <- DESeqDataSetFromMatrix(countData = counts_bd21, #counts
                                    colData = samplelabels %>% filter(Accession == 'BD21')  %>% select(Treatment),
                                    design= ~ Treatment)
 dds_bd21 <- DESeq(dds_bd21)
 res_bd21 <- results(dds_bd21)
 DE_gene_bd21 <- 
   tibble(gene = rownames(res_bd21), log2Foldchange = res_bd21$log2FoldChange, padj = res_bd21$padj) %>% 
   filter(padj < 0.0001) %>% 
   pull(gene)

 counts_bd21_DE <-  filter(counts_bd21, rownames(counts_bd21) %in% DE_gene_bd21)

 
 drought <- counts_bd21_DE[,samplelabels %>% filter(Treatment == 'drought'& Accession == 'BD21') %>% pull('rownames(matrixbd)')]
 rownames(drought) <- rownames(counts_bd21_DE)
 normal <- counts_bd21_DE[,samplelabels %>% filter(Treatment == 'normal'& Accession == 'BD21') %>% pull('rownames(matrixbd)')]
 rownames(normal) <- rownames(counts_bd21_DE)
 cor_drought <- correlate(t(drought))
 cor_normal  <- correlate(t(normal))
 
 
 main_fig_2_suppl_4_1 <- bind_rows(as.vector(as.matrix(cor_drought[,-1])) %>% as_tibble() %>% 
                                     filter(is.na(value)==FALSE) %>% mutate(class = 'Drought'),
                                   as.vector(as.matrix(cor_normal[,-1])) %>% as_tibble() %>% 
                                     filter(is.na(value)==FALSE) %>% mutate(class = 'Control')) %>% 
   ggplot(aes(value, color = class, fill = class))+
   theme_bw()+
   theme(legend.title = element_blank(), legend.position = 'top')+
   xlab('PCC')+
   ylab('Density')+
   scale_fill_manual(values=c("#81A1C1", "#A3BE8C"),name = ' ',labels = c("Control", 'Drought'))+
   scale_color_manual(values=c( "#81A1C1", "#A3BE8C"),name = ' ',labels = c("Control",'Drought'))+
   theme(
     #plot.title = element_text(color="black", size=5, face="bold"),
     axis.title.x = element_text(color="black", size=12, face="bold"),
     axis.title.y = element_text(color="black", size=12, face="bold"),
     legend.text = element_text(color = 'black', size = 12, face = 'bold')
   )+
   theme(
     axis.text.y=element_blank(),
     axis.ticks.y=element_blank())+
   theme(aspect.ratio=1/1)+
   geom_density(alpha = .4)
 
 
 
 
 
 
 
 
 

 counts_bd31 <- counts[,samplelabels %>% filter(Accession == 'BD3-1') %>% pull('rownames(matrixbd)')]
 rownames(counts_bd31) <- normalized_counts_with_chloroplast$X1
 counts_bd31 <-  filter(counts_bd31, rowMeans(counts[,1:192])>=5)
 
 dds_bd31 <- DESeqDataSetFromMatrix(countData = counts_bd31, #counts
                                    colData = samplelabels %>% filter(Accession == 'BD3-1')  %>% select(Treatment),
                                    design= ~ Treatment)
 dds_bd31 <- DESeq(dds_bd31)
 res_bd31 <- results(dds_bd31)
 DE_gene_bd31 <- 
   tibble(gene = rownames(res_bd31), log2Foldchange = res_bd31$log2FoldChange, padj = res_bd31$padj) %>% 
   filter(padj < 0.0001) %>% 
   pull(gene)
 
 
 counts_bd31_DE <-  filter(counts_bd31, rownames(counts_bd31) %in% DE_gene_bd31)
 
 
 drought <- counts_bd31_DE[,samplelabels %>% filter(Treatment == 'drought'& Accession == 'BD3-1') %>% pull('rownames(matrixbd)')]
 rownames(drought) <- rownames(counts_bd31_DE)
 normal <- counts_bd31_DE[,samplelabels %>% filter(Treatment == 'normal'& Accession == 'BD3-1') %>% pull('rownames(matrixbd)')]
 rownames(normal) <- rownames(counts_bd31_DE)
 cor_drought <- correlate(t(drought))
 cor_normal  <- correlate(t(normal))
 
 
 main_fig_2_suppl_4_2 <- bind_rows(as.vector(as.matrix(cor_drought[,-1])) %>% as_tibble() %>% 
                                     filter(is.na(value)==FALSE) %>% mutate(class = 'Drought'),
                                   as.vector(as.matrix(cor_normal[,-1])) %>% as_tibble() %>% 
                                     filter(is.na(value)==FALSE) %>% mutate(class = 'Control')) %>% 
   ggplot(aes(value, color = class, fill = class))+
   theme_bw()+
   theme(legend.title = element_blank(), legend.position = 'top')+
   xlab('PCC')+
   ylab('Density')+
   scale_fill_manual(values=c("#81A1C1", "#A3BE8C"),name = ' ',labels = c("Control", 'Drought'))+
   scale_color_manual(values=c( "#81A1C1", "#A3BE8C"),name = ' ',labels = c("Control",'Drought'))+
   theme(
     #plot.title = element_text(color="black", size=5, face="bold"),
     axis.title.x = element_text(color="black", size=12, face="bold"),
     axis.title.y = element_text(color="black", size=12, face="bold"),
     legend.text = element_text(color = 'black', size = 12, face = 'bold')
   )+
   theme(
     axis.text.y=element_blank(), 
     axis.ticks.y=element_blank())+
   theme(aspect.ratio=1/1)+
   geom_density(alpha = .4)
 
 
 
 
 
 main_fig_2_suppl_4 <-   main_fig_2_suppl_4_1 + main_fig_2_suppl_4_2 + plot_annotation(tag_levels = 'A')+plot_layout(guides = 'collect')&theme(legend.position = 'bottom')
 ggsave(filename = 'main_fig_2_suppl_4.pdf')   
 
 
 

# Main_figure_3 ------------------------------------
 
 load("../Data/CCF_lag_1.rdata")#load max cross correlation for each pair
 load("../Data/df_heat_pop.rdata") #load population correlation for each pair under heat condition
 CCF <- CCF_lag_1 %>% rename(Class = CLASS) %>% 
   mutate(Class = replace(Class, Class == 'CONTROL', 'Control')) %>% 
   mutate(Class = replace(Class, Class == 'HEAT', 'Heat')) %>% 
   mutate(Class = replace(Class, Class == 'DROUGHT', 'Drought'))
 CCF<- CCF %>% 
   rename(regulator = X1, target = X2)
 
 
 #get TF family
 TF_family <- read_delim("../Data/gene_family.txt", 
                         "\t", escape_double = FALSE, trim_ws = TRUE)
 TF_family <- TF_family %>% as_tibble() %>% 
   rename(gene = `Protein ID`) %>% 
   mutate(gene = str_sub(gene, 1,14)) %>% dplyr::select(gene, Family) 
 TF_family <- TF_family %>% as_tibble() %>% unique()
 
 
 
 fig_3_3 <-  df_heat %>% 
   filter(Time < 120) %>% 
   left_join(TF_family,  by = c('X1' = 'gene')) %>% 
   filter(Family == 'HSF') %>% 
   mutate(cor = ifelse(is.na(cor) == TRUE, 0, cor)) %>% 
   ggplot(aes(cor, factor(Time) , color = CLASS, fill = CLASS))+
   geom_boxplot(alpha = .8,notch = TRUE,varwidth = TRUE, size = 1.2)+
   theme_light()+
   theme(legend.title = element_blank())+
   theme(legend.position="bottom")+
   xlab('PCC')+
   ylab('Time (min)')+
   scale_fill_manual(values=c("#81A1C1", '#B48EAD'),name = ' ',labels = c("Control", "Heat"))+
   scale_color_manual(values=c( "#81A1C1", '#B48EAD'),name = ' ',labels = c("Control", "Heat"))+
   theme(
     axis.title = element_text(color="black", size=15),
     legend.text = element_text(color = 'black', size = 15)
   )+
   theme(aspect.ratio=1.5/1)+
   theme(plot.title = element_text(size = 16, face = 'bold', vjust = 0.7),
         axis.text=element_text(size = 12))
 
 
 
 
 
 fig_3_4 <- CCF %>% as_tibble() %>% 
   filter(Class == 'Heat' | Class == 'Control') %>% 
   left_join(TF_family,  by = c('regulator' = 'gene')) %>% 
   filter(Family == 'HSF') %>% 
   ggplot(aes(x = max_ccf, color = Class, fill = Class))+
   geom_boxplot(alpha = .8,notch = TRUE,varwidth = TRUE, size = 1.2)+
   theme_light()+
   theme(legend.title = element_blank())+
   theme(legend.position="none")+
   xlab('MCC')+
   scale_fill_manual(values=c("#81A1C1", '#B48EAD'),name = ' ',labels = c("Control", "Heat"))+
   scale_color_manual(values=c( "#81A1C1", '#B48EAD'),name = ' ',labels = c("Control", "Heat"))+
   theme(
     axis.title = element_text(color="black", size=15),
     legend.text = element_text(color = 'black', size = 15)
   )+
   theme(aspect.ratio=1/1)+
   theme(axis.title.y=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks.y=element_blank())+
   theme(axis.text.x=element_text(size = 12))
 
 
 layout <- "
AA##
BB##"
 
 fig_3 <-  (fig_3_3) / ( fig_3_4)+  plot_annotation(tag_levels = 'A')+plot_layout(widths = c(1, 2), heights  = c(1.5,1))& theme(plot.tag = element_text(size = 16))
 ggsave(fig_3, filename = 'main_fig_3.pdf', height = 12, width = 5.2)
 
 
 
 
# Suppl_figure_for_main_figure_3 --------------------------------------------------------
 

 
 load("../Data/df_drought_pop.rdata")
 
 
 plot_drought <- list()
 
 for (i in 1:9) {
   
   
   plot_drought[[i]] <-  df_drought %>%
     filter(Time == i * 15) %>% 
     ggplot(aes(cor,  color = CLASS))+geom_density(alpha = .8, size = 1)+
     theme_light()+labs(title = paste0(" Time after treatment(min): ",time), x = 'PCC', y = 'Density')+
     guides(fill = FALSE, size = FALSE)+
     theme(legend.title = element_blank())+
     scale_color_discrete(name = ' ',labels = c("Control", "Drought"))+
     scale_color_manual(values=c( "#81A1C1", "#A3BE8C"))+
     theme(
       plot.title = element_text(color="black", size=5, face="bold"),
       axis.title.x = element_text(color="black", size=5, face="bold"),
       axis.title.y = element_text(color="black", size=5, face="bold"),
       legend.text = element_text(color = 'black', size = 5, face = 'bold')
     )
   
   
 }
 
 
 
 
 #https://stackoverflow.com/questions/65291723/merging-two-y-axes-titles-in-patchwork
 fig_3_supp_1 <- ((plot_drought[[1]] + plot_drought[[2]] + plot_drought[[3]]+ plot_layout(guides = 'collect', ncol = 3)) / (plot_drought[[4]] + plot_drought[[5]] + plot_drought[[6]]+ plot_layout(guides = 'collect', ncol = 3)) / (plot_drought[[7]] + plot_drought[[8]] + plot_drought[[9]]+ plot_layout(guides = 'collect', ncol = 3))) 
 ggsave(fig_3_supp_1, filename = 'drought_pop_cor_over_time.pdf', height = 10, width = 11.5)
 
 
 
 
 #heat
 plot_heat <- list()
 for (i in c(1:16)){
   
   data_filter <-   data %>% ungroup() %>% 
     filter( CLASS == "CONTROL" | CLASS == 'HEAT') %>% 
     filter(MINUTES == i*15) %>% 
     select(-MINUTES,-key, -CULTIVAR) %>% 
     nest(value = c(value))
   
   time = i * 15
   
   df <- RNA_seq_summary_Experiment_ID %>%
     dplyr::select(CLASS) %>%
     unique() %>%
     slice(1,2) %>%
     mutate(prior = map(CLASS, ~network_prior)) %>%
     unnest(prior) %>%
     left_join(data_filter, by = c('X1' = 'gene', 'CLASS' = 'CLASS')) %>%
     rename(regulator_value = value) %>%
     left_join(data_filter, by = c('X2' = 'gene', 'CLASS' = 'CLASS')) %>%
     rename(target_value = value) #%>%
   
   
   
   df <- df %>%
     ungroup() %>% 
     rowwise() %>% 
     mutate(cor = cor(regulator_value %>% 
                        pull(),
                      target_value %>% 
                        pull(), use = 'everything',method = 'pearson')) 
   
   
   
   plot_heat[[i]] <-  df %>% 
     ggplot(aes(cor,  color = CLASS))+geom_density(alpha = .8, size = 1)+
     theme_light()+labs(title = paste0("  Time after treatment (min): ",time), x = 'PCC', y = 'Density')+
     guides(fill = FALSE, size = FALSE)+
     theme(legend.title = element_blank())+
     scale_color_discrete(name = ' ',labels = c("Control", "Heat"))+
     scale_color_manual(values=c( "#81A1C1", '#B48EAD'))+
     theme(
       plot.title = element_text(color="black", size=10, face="bold"),
       axis.title.x = element_text(color="black", size=10, face="bold"),
       axis.title.y = element_text(color="black", size=10, face="bold"),
       legend.text = element_text(color = 'black', size = 7, face = 'bold')
     )
   
 }
 
 
 fig_3_supp_2 <- ((plot_heat[[1]] + plot_heat[[2]] + plot_heat[[3]]+ plot_heat[[4]]+plot_layout(guides = 'collect', ncol = 4)) / (plot_heat[[5]] + plot_heat[[6]] + plot_heat[[7]]+plot_heat[[8]]+ plot_layout(guides = 'collect', ncol = 4)) / (plot_heat[[9]] + plot_heat[[10]] + plot_heat[[11]]+ plot_heat[[12]]+ plot_layout(guides = 'collect', ncol = 4)) + (plot_heat[[13]] + plot_heat[[14]] + plot_heat[[15]] + plot_heat[[16]]+ plot_layout(guides = 'collect', ncol = 4))) 
 ggsave(fig_3_supp_2, filename = 'heat_pop_cor_over_time.pdf', height = 11, width = 11)
 
 
 
 
 
 
 
 
 
 
 load('../Data/df_drought_pop.rdata')
 
 
 
 fig_3_suppl_3 <- df_drought %>% 
   mutate(cor = ifelse(is.na(cor) == TRUE, 0, cor)) %>% 
   ggplot(aes(factor(Time), cor, color = CLASS, fill = CLASS))+
   geom_boxplot(alpha = .8,notch = TRUE,varwidth = TRUE,size = 1.2)+
   theme_light()+
   theme(legend.title = element_blank())+
   theme(legend.position="top")+
   ylab('PCC')+
   xlab('Time (min)')+
   scale_fill_manual(values=c("#81A1C1", "#A3BE8C"),name = ' ',labels = c("Control", "Drought"))+
   scale_color_manual(values=c( "#81A1C1", "#A3BE8C"),name = ' ',labels = c("Control", "Drought"))+
   theme(
     axis.title.x = element_text(color="black", size=12, face="bold"),
     axis.title.y = element_text(color="black", size=12, face="bold"),
     legend.text = element_text(color = 'black', size = 12, face = 'bold')
   )+
   theme(
     axis.text.x=element_text(size = 10),
     axis.text.y=element_text(size = 10))
 
 
 
 
 
 ##Heat
 
 
 load('../Data/df_heat_pop.rdata')
 
 
 fig_3_suppl_4<-df_heat %>% 
   mutate(cor = ifelse(is.na(cor) == TRUE, 0, cor)) %>% 
   ggplot(aes(factor(Time), cor, color = CLASS,fill = CLASS))+
   geom_boxplot(alpha = .8,notch = TRUE,varwidth = TRUE, size = 1.2)+
   theme_light()+
   theme(legend.title = element_blank())+
   theme(legend.position="top")+
   ylab('PCC')+
   xlab('Time (min)')+
   scale_fill_manual(values=c("#81A1C1", '#B48EAD'),name = ' ',labels = c("Control", "Heat"))+
   scale_color_manual(values=c( "#81A1C1", '#B48EAD'),name = ' ',labels = c("Control", "Heat"))+
   theme(
     #plot.title = element_text(color="black", size=5, face="bold"),
     axis.title.x = element_text(color="black", size=12, face="bold"),
     axis.title.y = element_text(color="black", size=12, face="bold"),
     legend.text = element_text(color = 'black', size = 12, face = 'bold')
   )+
   theme(
     axis.text.x=element_text(size = 10),
     axis.text.y=element_text(size = 10))
 
 
 
 
 fig_3_suppl_3 / fig_3_suppl_4 + plot_annotation(tag_levels = 'A')
 ggsave('main_fig_3_suppl_3.pdf')
 