#This Code is Meant to Normalize to Control Gene (ie non-protein coding gene is set to 1)
#Enter data on excel sheet and change name of file on line 7
#To alter the names of columns based on different treatments, change treatment vector and gene vector on line 61 and 62
#Change labels on line 94

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
treated <- filter(ave_delta_cp, ControlTreatment == "n")
treated <- treated[order(treated$ControlGene), c(2,4:5)]
treated_ddcp_gene1 <- -(diff(treated$ave_dcp))
treated_ddcp_gene2 <- 0
control <- filter(ave_delta_cp, ControlTreatment == "y")
control <- control[order(control$ControlGene), c(2,4:5)]
control_ddcp_gene1 <- -(diff(control$ave_dcp))
control_ddcp_gene2 <- 0

#s
treated_sd <- as.vector(treated$sd_delta_cp)
treated_s_gene1 <- ((treated_sd[1]^2) + (treated_sd[2]^2))^(1/2)
treated_s_gene2 <- ((treated_sd[1]^2) + (treated_sd[1]^2))^(1/2)
control_sd <- as.vector(control$sd_delta_cp)
control_s_gene1 <- ((control_sd[1]^2) + (control_sd[2]^2))^(1/2)
control_s_gene2 <- ((control_sd[1]^2) + (control_sd[1]^2))^(1/2)

#put into df
ddcp <- data.frame(treated_ddcp_gene1, treated_ddcp_gene2, control_ddcp_gene1, control_ddcp_gene2)
s <- data.frame(treated_s_gene1, treated_s_gene2, control_s_gene1, control_s_gene2)

#label with your treatments and controls
treatment_vector <- c('Rifampicin', 'Rifampicin', 'NoRifampicin', 'NoRifampicin')
gene_vector <- c('ProteinCoding', 'NonProteinCoding', 'ProteinCoding', 'NonProteinCoding')
treatment_number <- c('1', '2', '3', '4')
ddcp_long <- gather(ddcp,
                    key = "type",
                    value = "ddcp",
                    treated_ddcp_gene1:control_ddcp_gene2)
s_long <- gather(s,
                 key = "type",
                 value = "sd",
                 treated_s_gene1:control_s_gene2)
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


