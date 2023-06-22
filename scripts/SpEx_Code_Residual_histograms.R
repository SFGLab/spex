#----------------------------------
setwd("/mnt/raid/abhishek/spex/expression/")
#---------------------------------

suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))
suppressMessages(library(tidyverse))
suppressMessages(library(VennDiagram))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggplot2))
suppressMessages(library(data.table))
suppressMessages(library(mixtools))

theme_DP <- theme_bw() +
  theme(panel.grid.major=element_line(color="white",linetype="dashed"),
        panel.grid.minor=element_blank(),
        axis.title=element_text(size=rel(1.25)),
        axis.text=element_text(size=rel(1.1)),
        legend.position= c(.98, .98),
        legend.justification=c("right","top"),
        legend.text=element_text(size=rel(0.75)))

#-----------------------------------------------------------------------#
#Expecto
exp = (read.csv(file = "data/GM12878_EXPECTO.csv", header = T, sep = ",", row.names = 2)[,-c(1)])
colnames(exp) = c("Predicted_Value","Real_Value")
exp$Residuals = exp$Real_Value - exp$Predicted_Value
setDT(exp, keep.rownames = "Genes")[]
exp <- exp[order(exp$Genes),]
exp = exp[,-c(2:3)]
exp_1 = exp %>%  remove_rownames() %>% column_to_rownames(var = 'Genes')

#CTCF
ctcf = (read.csv(file = "data/GM12878_CTCF.csv", header = T, sep = ",", row.names = 2)[,-c(1)])
colnames(ctcf) = c("Predicted_Value","Real_Value")
ctcf$Residuals = ctcf$Real_Value - ctcf$Predicted_Value
setDT(ctcf, keep.rownames = "Genes")[]
ctcf <- ctcf[order(ctcf$Genes),]
r_ctcf = ctcf[,-c(2:3)]
r1_ctcf = r_ctcf %>%  remove_rownames() %>% column_to_rownames(var = 'Genes')

#Cohesin
cohesin = (read.csv(file = "data/GM12878_COHESIN.csv", header = T, sep = ",", row.names = 2)[,-c(1)])
colnames(cohesin) = c("Predicted_Value","Real_Value")
cohesin$Residuals = cohesin$Real_Value - cohesin$Predicted_Value
setDT(cohesin, keep.rownames = "Genes")[]
cohesin <- cohesin[order(cohesin$Genes),]
r_cohesin = cohesin[,-c(2:3)]
r1_cohesin = r_cohesin %>%  remove_rownames() %>% column_to_rownames(var = 'Genes')

#RNAP2
rnap2 = (read.csv(file = "data/GM12878_RNAPOL2.csv", header = T, sep = ",", row.names = 2)[,-c(1)])
colnames(rnap2) = c("Predicted_Value","Real_Value")
rnap2$Residuals = rnap2$Real_Value - rnap2$Predicted_Value
setDT(rnap2, keep.rownames = "Genes")[]
rnap2 <- rnap2[order(rnap2$Genes),]
r_rnap2 = rnap2[,-c(2:3)]
r1_rnap2 = r_rnap2 %>%  remove_rownames() %>% column_to_rownames(var = 'Genes')

#SpEx
spex = cbind(r_ctcf,r_cohesin,r_rnap2)
spex = spex[,-c(3,5)]
colnames(spex) = c("Genes", "CTCF","Cohesin","RNAP2")
spex= spex %>%  remove_rownames() %>% column_to_rownames(var = 'Genes')
spex$spex <- apply(spex, 1, function(x) x[which.min(abs(x))])
spex = setDT(spex, keep.rownames = "Genes")[]
spex = spex[,-c(2:4)]
colnames(spex) = c("Genes","Residuals")
r1_spex = spex %>%  remove_rownames() %>% column_to_rownames(var = 'Genes')

#-----------------------------------------------------------------------#
#Cutoff
out <- normalmixEM(spex$Residuals, k=2)
summary(out)
lower_spex <- (out$mu[1] - (out$sigma[1]/2))
upper_spex <- (out$mu[2] + (out$sigma[2]/2))

h_spex <- hist(r1_spex$Residuals, breaks = 25, plot = FALSE)
df_spex <- data.frame(x = h_spex$mids, y = h_spex$density, fill = h_spex$mids > lower_spex & h_spex$mids < upper_spex)

h_ctcf <- hist(r1_ctcf$Residuals, breaks = 25, plot = FALSE)
df_ctcf <- data.frame(x = h_ctcf$mids, y = h_ctcf$density, fill = h_ctcf$mids > lower_spex & h_ctcf$mids < upper_spex)

h_cohesin <- hist(r1_cohesin$Residuals, breaks = 25, plot = FALSE)
df_cohesin <- data.frame(x = h_cohesin$mids, y = h_cohesin$density, fill = h_cohesin$mids > lower_spex & h_cohesin$mids < upper_spex)

h_rnap2 <- hist(r1_rnap2$Residuals, breaks = 25, plot = FALSE)
df_rnap2 <- data.frame(x = h_rnap2$mids, y = h_rnap2$density, fill = h_rnap2$mids > lower_spex & h_rnap2$mids < upper_spex)

h_exp <- hist(exp$Residuals, breaks = 25, plot = FALSE)
df_exp <- data.frame(x = h_exp$mids, y = h_exp$density, fill = h_exp$mids > lower_spex & h_exp$mids < upper_spex)

#-----------------------------------------------------------------------#

spex_cut.off <- r1_spex$Residuals >= lower_spex & r1_spex$Residuals <= upper_spex
spex_sig <- r1_spex[spex_cut.off,]

ctcf_cut.off <- ctcf$Residuals >= lower_spex & ctcf$Residuals <= upper_spex
ctcf_sig <- ctcf[ctcf_cut.off, ]
ctcf_sig_gene = ctcf_sig$Genes

cohesin_cut.off <- cohesin$Residuals >= lower_spex & cohesin$Residuals <= upper_spex
cohesin_sig <- cohesin[cohesin_cut.off, ]
cohesin_sig_gene = cohesin_sig$Genes

rnap2_cut.off <- rnap2$Residuals >= lower_spex & rnap2$Residuals <= upper_spex
rnap2_sig <- rnap2[rnap2_cut.off, ]
rnap2_sig_gene = rnap2_sig$Genes

exp_cut.off <- exp$Residuals >= lower_spex & exp$Residuals <= upper_spex
exp_sig <- exp[exp_cut.off, ]
exp_sig_gene = exp_sig$Genes

#-----------------------------------------------------------------------#

p_spex = ggplot(df_spex) +
  geom_col(aes(x, y, fill = fill), width = 1, color = 'black') +
  geom_density(data = data.frame(x = spex$Residuals), 
               aes(x = x, color = 'Actual density'),
               key_glyph = 'path') +
  geom_function(fun = function(x) {
    dnorm(x, mean = mean(spex$Residuals), sd = sd(spex$Residuals)) },
    aes(color = 'Theoretical density')) +
  scale_fill_manual(values = c(`TRUE` = '#FF9933', 'FALSE' = 'gray'), 
                    name = 'Within Cutoff') +
  scale_color_manual(values = c('black', 'brown'), name = 'Density lines') +
  scale_x_continuous(breaks = pretty(spex$Residuals, n = 10)) +
  labs(x = 'Residuals', y = 'Density') +
  ggtitle("SpEx - Best") +
  theme_DP

p_ctcf = ggplot(df_ctcf) +
  geom_col(aes(x, y, fill = fill), width = 1, color = 'black') +
  geom_density(data = data.frame(x = r1_ctcf$Residuals), 
               aes(x = x, color = 'Actual density'),
               key_glyph = 'path') +
  geom_function(fun = function(x) {
    dnorm(x, mean = mean(r1_ctcf$Residuals), sd = sd(r1_ctcf$Residuals)) },
    aes(color = 'Theoretical density')) +
  scale_fill_manual(values = c(`TRUE` = '#5C83AC', 'FALSE' = 'gray'), 
                    name = 'Within Cutoff') +
  scale_color_manual(values = c('black', 'brown'), name = 'Density lines') +
  scale_x_continuous(breaks = pretty(r1_ctcf$Residuals, n = 10)) +
  labs(x = 'Residuals', y = 'Density') +
  ggtitle("SpEx - CTCF") +
  theme_DP

p_cohesin = ggplot(df_cohesin) +
  geom_col(aes(x, y, fill = fill), width = 1, color = 'black') +
  geom_density(data = data.frame(x = r1_cohesin$Residuals), 
               aes(x = x, color = 'Actual density'),
               key_glyph = 'path') +
  geom_function(fun = function(x) {
    dnorm(x, mean = mean(r1_cohesin$Residuals), sd = sd(r1_cohesin$Residuals)) },
    aes(color = 'Theoretical density')) +
  scale_fill_manual(values = c(`TRUE` = '#5C8000', 'FALSE' = 'gray'), 
                    name = 'Within Cutoff') +
  scale_color_manual(values = c('black', 'brown'), name = 'Density lines') +
  scale_x_continuous(breaks = pretty(r1_cohesin$Residuals, n = 10)) +
  labs(x = 'Residuals', y = 'Density') +
  ggtitle("SpEx - Cohesin") +
  theme_DP

p_rnap2 = ggplot(df_rnap2) +
  geom_col(aes(x, y, fill = fill), width = 1, color = 'black') +
  geom_density(data = data.frame(x = r1_rnap2$Residuals), 
               aes(x = x, color = 'Actual density'),
               key_glyph = 'path') +
  geom_function(fun = function(x) {
    dnorm(x, mean = mean(r1_rnap2$Residuals), sd = sd(r1_rnap2$Residuals)) },
    aes(color = 'Theoretical density')) +
  scale_fill_manual(values = c(`TRUE` = '#FF9933', 'FALSE' = 'gray'), 
                    name = 'Within Cutoff') +
  scale_color_manual(values = c('black', 'brown'), name = 'Density lines') +
  scale_x_continuous(breaks = pretty(r1_rnap2$Residuals, n = 10)) +
  labs(x = 'Residuals', y = 'Density') +
  ggtitle("SpEx - RNA POL II") +
  theme_DP

p_exp = ggplot(df_exp) +
  geom_col(aes(x, y, fill = fill), width = 1, color = 'black') +
  geom_density(data = data.frame(x = exp$Residuals), 
               aes(x = x, color = 'Actual density'),
               key_glyph = 'path') +
  geom_function(fun = function(x) {
    dnorm(x, mean = mean(exp$Residuals), sd = sd(exp$Residuals)) },
    aes(color = 'Theoretical density')) +
  scale_fill_manual(values = c(`TRUE` = '#D56565', 'FALSE' = 'gray'), 
                    name = 'Within Cutoff') +
  scale_color_manual(values = c('black', 'brown'), name = 'Density lines') +
  scale_x_continuous(breaks = pretty(exp$Residuals, n = 10)) +
  labs(x = 'Residuals', y = 'Density') +
  ggtitle("Baseline") +
  theme_DP

#-----------------------------------------------------------------------#

spex$Methord = "SpEx - Best"
spex= spex %>%  remove_rownames() %>% column_to_rownames(var = 'Genes')

exp_1$Methord = "Baseline"
dist_exp = rbind(spex, exp_1)

r1_ctcf$Methord = "SpEx - CTCF"
dist_ctcf = rbind(r1_ctcf, exp_1)

r1_cohesin$Methord = "SpEx - Cohesin"
dist_cohesin = rbind(r1_cohesin, exp_1)

r1_rnap2$Methord = "SpEx - RNA POL II"
dist_rnap2 = rbind(r1_rnap2, exp_1)

p_dist_exp = ggplot(dist_exp, aes(x=Residuals,fill = Methord)) +
  geom_density(alpha = 0.25) + ggtitle("Baseline vs SpEx - Best")  + 
  labs(x = "Residuals", y = "Density") + theme_DP

plot_dist_ctcf = ggplot(dist_ctcf, aes(x=Residuals,fill = Methord)) +
  geom_density(alpha = 0.25) + ggtitle("Baseline vs SpEx - CTCF")  + 
  labs(x = "Residuals", y = "Density") + theme_DP

plot_dist_cohesin = ggplot(dist_cohesin, aes(x=Residuals,fill = Methord)) + 
  geom_density(alpha = 0.25)+ ggtitle("Baseline vs SpEx - Cohesin") + 
  labs(x = "Residuals", y = "Density") + theme_DP

plot_dist_rnap2 = ggplot(dist_rnap2, aes(x=Residuals,fill = Methord)) +
  geom_density(alpha = 0.25) + ggtitle("Baseline vs SpEx - RNA POL II") +  
  labs(x = "Residuals", y = "Density") + theme_DP

#-----------------------------------------------------------------------#

patchwork = p_dist_exp / (p_exp | p_spex) / (plot_dist_ctcf | plot_dist_cohesin | plot_dist_rnap2) /(p_ctcf | p_cohesin | p_rnap2)
patchwork + plot_annotation(title = "Distribution of Residual Value on GM12878 cell line",subtitle = "[Residual Value | Cut-Off = -1.39 / +2.10 of SpEx]")  & theme(plot.title = element_text(face = 'bold',size = 15), plot.subtitle = element_text(face = 'bold',size = 12,colour = "#660000"))

#-----------------------------------------------------------------------#
