#This Code is Meant to Normalize to Control Treatment (ie no rif is set to 1)
#Enter data on excel sheet and change name of file on line 6
#To alter the names of columns based on different treatments, change treatment and gene vectors on lines 59 and 60
#To change names on graphs, change cols vector on line 94

library("tidyverse")
raw_cps <- read_csv("raw_data_qpcr.csv")
raw_long <- gather(raw_cps,
                             key = "replicate",
                             value = "cp_value",
                             InputRep1:IPRep3)
#this gets summary stats for our three technical replicates
average_cps_inputs <- filter(raw_long, grepl("Input", replicate)) %>% 
                        group_by(Sample) %>% 
                        summarize(Inputaverage=mean(cp_value), Inputstandard_deviation=sd(cp_value))

average_cps_IPs <- filter(raw_long, grepl("IP", replicate)) %>% 
  group_by(Sample) %>% 
  summarize(IPaverage=mean(cp_value), IPstandard_deviation=sd(cp_value))


#combine averages with description of treatments
treatments <- subset(raw_cps, select = -(InputRep1:IPRep3))
summary_cps <- left_join(average_cps_IPs, treatments, by = "Sample")
summary_cps <- left_join(average_cps_inputs, summary_cps, by = "Sample")

#find delta cp
summary_cps$delta_cps <- (summary_cps$IPaverage - summary_cps$Inputaverage)

#average delta cp
ave_delta_cp <- summary_cps %>% group_by(Treatment) %>%
                            summarize(ave_dcp = mean(delta_cps), sd_delta_cp = sd(delta_cps))
treatments <- subset(treatments, select = -(Sample))
treatments <- treatments[order(treatments$Treatment), ]
treatments <- treatments[!duplicated(treatments), ]
ave_delta_cp <- left_join(treatments, ave_delta_cp, by = "Treatment")

#delta delta cp
gene1 <- filter(ave_delta_cp, ControlGene == "n")
gene1 <- gene1[order(gene1$ControlTreatment), c(1,4:5)]
gene1_ddcp_trtd <- -(diff(gene1$ave_dcp))
gene1_ddcp_control <- 0
gene2 <- filter(ave_delta_cp, ControlGene == "y")
gene2 <- gene2[order(gene2$ControlTreatment), c(1,4:5)]
gene2_ddcp_trtd <- -(diff(gene2$ave_dcp))
gene2_ddcp_control <- 0

#s
gene1_sd <- as.vector(gene1$sd_delta_cp)
gene1_s_trtd <- ((gene1_sd[1]^2) + (gene1_sd[2]^2))^(1/2)
gene1_s_control <- ((gene1_sd[1]^2) + (gene1_sd[1]^2))^(1/2)
gene2_sd <- as.vector(gene2$sd_delta_cp)
gene2_s_trtd <- ((gene2_sd[1]^2) + (gene2_sd[2]^2))^(1/2)
gene2_s_control <- ((gene2_sd[1]^2) + (gene2_sd[1]^2))^(1/2)

#put into df
ddcp <- data.frame(gene1_ddcp_trtd, gene1_ddcp_control, gene2_ddcp_trtd, gene2_ddcp_control)
s <- data.frame(gene1_s_trtd, gene1_s_control, gene2_s_trtd, gene2_s_control)
#label with your treatments and controls
treatment_vector <- c('Rifampicin', 'NoRifampicin', 'Rifampicin', 'NoRifampicin')
gene_vector <- c('ProteinCoding', 'ProteinCoding', 'NonProteinCoding', 'NonProteinCoding')
treatment_number <- c('1', '2', '3', '4')
ddcp_long <- gather(ddcp,
                key = "type",
                value = "ddcp",
                gene1_ddcp_trtd:gene2_ddcp_control)
s_long <- gather(s,
                    key = "type",
                    value = "sd",
                    gene1_s_trtd:gene2_s_control)
ddcp_long <- cbind(ddcp_long, treatment=treatment_vector)
ddcp_long <- cbind(ddcp_long, gene=gene_vector)
ddcp_long <- cbind(ddcp_long, treat_num=treatment_number)
s_long <- cbind(s_long, treatment=treatment_vector)
s_long <- cbind(s_long, gene=gene_vector)
s_long <- cbind(s_long, treat_num=treatment_number)

#normalized for primer efficiency
ddcp_long$norm_ddcp <- cbind(1.8^(-ddcp_long$ddcp))


#combine_columns
norm_ddcp_and_s <- left_join(ddcp_long, s_long, by = "treat_num")

#plus and minus s
norm_ddcp_and_s$ddcp_plus_s <- (norm_ddcp_and_s$ddcp + norm_ddcp_and_s$sd)
norm_ddcp_and_s$ddcp_minus_s <- (norm_ddcp_and_s$ddcp - norm_ddcp_and_s$sd)

#normalize plus and minus s
norm_ddcp_and_s$norm_ddcp_plus_s <- cbind(1.8^(-norm_ddcp_and_s$ddcp_plus_s))
norm_ddcp_and_s$norm_ddcp_minus_s <- cbind(1.8^(-norm_ddcp_and_s$ddcp_minus_s))

#plot
cols <- c("NoRifampicin" ="blue", "Rifampicin" = "red")
ggplot(norm_ddcp_and_s, aes(gene.x, norm_ddcp, fill=treatment.x)) + 
  geom_bar(stat="identity", position=position_dodge(width=0.8), width=0.5) +
  geom_errorbar(aes(ymin=norm_ddcp_minus_s, ymax=norm_ddcp_plus_s), position = position_dodge(0.8), width = 0.1) +
  scale_fill_manual(values=cols, name="Treatment") + 
  labs(x="Gene", y="Protein Fold Enrichment")

  