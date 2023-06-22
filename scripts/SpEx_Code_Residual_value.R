#----------------------------------
setwd("/mnt/raid/abhishek/spex/expression/")
#----------------------------------

suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))
suppressMessages(library(tidyverse))
suppressMessages(library(VennDiagram))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggplot2))
suppressMessages(library(data.table))
suppressMessages(library(mixtools))
#----------------------------------

theme_DHO <- theme_bw() +
  theme(panel.grid.major=element_line(color="white",linetype="dashed"),
        panel.grid.minor=element_blank(),
        axis.title=element_text(size=rel(1.25)),
        axis.text=element_text(size=rel(1.1)),
        legend.position= "none",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.justification=c(-0.05,1.02),
        legend.title=element_blank(),
        legend.text=element_text(size=rel(1.1)))
theme_DP <- theme_bw() +
  theme(panel.grid.major=element_line(color="white",linetype="dashed"),
        panel.grid.minor=element_blank(),
        axis.title=element_text(size=rel(1.25)),
        axis.text=element_text(size=rel(1.1)),
        legend.position= c(.98, .98),
        legend.justification=c("right","top"),
        legend.title=element_blank(),
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

#Q-Q Plot
qq_exp = ggplot(exp, aes(sample = Residuals)) + ggtitle("Baseline") + xlab("Theoretical Quantiles") + ylab("Sample Quantiles") + stat_qq() + stat_qq_line(col = "red") + theme_DHO
qq_CTCF = ggplot(r_ctcf, aes(sample = Residuals)) + ggtitle("SpEx - CTCF") + xlab("Theoretical Quantiles") + ylab("Sample Quantiles") + stat_qq() + stat_qq_line(col = "red") + theme_DHO
qq_Cohesin = ggplot(r_cohesin, aes(sample = Residuals)) + ggtitle("SpEx - Cohesin") + xlab("Theoretical Quantiles") + ylab("Sample Quantiles") + stat_qq() + stat_qq_line(col = "red") + theme_DHO
qq_RNAP2 = ggplot(r_rnap2, aes(sample = Residuals)) + ggtitle("SpEx - RNA POL II") + xlab("Theoretical Quantiles") + ylab("Sample Quantiles") + stat_qq() + stat_qq_line(col = "red") + theme_DHO
patchwork3 <- (qq_exp | qq_CTCF) / (qq_Cohesin | qq_RNAP2)
patchwork3 + plot_annotation(title = "Quantile-Quantile Plots")  & theme(plot.title = element_text(face = 'bold',size = 15)) & theme(plot.subtitle = element_text(face = 'bold',size = 10,color = "#660000"))

#Cutoff
out <- normalmixEM(spex$Residuals, k=2)
summary(out)
lower_spex <- (out$mu[1] - (out$sigma[1]/2))
upper_spex <- (out$mu[2] + (out$sigma[2]/2))
#-----------------------------------------------------------------------#

spex_cut.off <- spex$Residuals >= lower_spex & spex$Residuals <= upper_spex
spex_sig <- spex[spex_cut.off,]
spex_sig_gene = spex_sig$Genes
spex$code = ifelse(spex$Genes %in% spex_sig_gene, TRUE, FALSE)
exp_spex = exp
exp_spex$code = ifelse(exp_spex$Genes %in% spex_sig_gene, TRUE, FALSE)

ctcf_cut.off <- ctcf$Residuals >= lower_spex & ctcf$Residuals <= upper_spex
ctcf_sig <- ctcf[ctcf_cut.off, ]
ctcf_sig_gene = ctcf_sig$Genes
ctcf$code = ifelse(ctcf$Genes %in% ctcf_sig_gene, TRUE, FALSE)
exp_ctcf = exp
exp_ctcf$code = ifelse(exp$Genes %in% ctcf_sig_gene, TRUE, FALSE)

cohesin_cut.off <- cohesin$Residuals >= lower_spex & cohesin$Residuals <= upper_spex
cohesin_sig <- cohesin[cohesin_cut.off, ]
cohesin_sig_gene = cohesin_sig$Genes
cohesin$code = ifelse(cohesin$Genes %in% cohesin_sig_gene, TRUE, FALSE)
exp_cohesin = exp
exp_cohesin$code = ifelse(exp_cohesin$Genes %in% cohesin_sig_gene, TRUE, FALSE)

rnap2_cut.off <- rnap2$Residuals >= lower_spex & rnap2$Residuals <= upper_spex
rnap2_sig <- rnap2[rnap2_cut.off, ]
rnap2_sig_gene = rnap2_sig$Genes
rnap2$code = ifelse(rnap2$Genes %in% rnap2_sig_gene, TRUE, FALSE)
exp_rnap2 = exp
exp_rnap2$code = ifelse(exp_rnap2$Genes %in% rnap2_sig_gene, TRUE, FALSE)

exp_cut.off <- exp$Residuals >= lower_spex & exp$Residuals <= upper_spex
exp_sig <- exp[exp_cut.off, ]
exp_sig_gene = exp_sig$Genes
#-----------------------------------------------------------------------#

plot_ctcf = ggplot(ctcf, aes(x = Genes, y = Residuals)) +
  geom_point(aes(colour = code),ctcf) +
  geom_smooth(method = "lm", se = T) +
  labs(x = "Genes", y = "Residuals") +
  ylim(-10,10) +
  ggtitle("SpEx - CTCF") +
  geom_hline(yintercept=0,linetype="dashed") +
  theme_DHO +
  scale_color_manual(name = "Residuals", 
                     values = c("lightgrey", "#00cccc"))

plot_exp_1 = ggplot(exp_ctcf, aes(x = Genes, y = Residuals)) +
  geom_point(aes(colour = code),exp_ctcf) +
  geom_smooth(method = "lm", se = T) +
  labs(x = "Genes", y = "Residuals") +
  ylim(-10,10) +
  ggtitle("Baseline") +
  geom_hline(yintercept=0,linetype="dashed") +
  theme_DHO +
  scale_color_manual(name = "Residuals", 
                     values = c("lightgrey", "darkblue"))

plot_cohesin = ggplot(cohesin, aes(x = Genes, y = Residuals)) +
  geom_point(aes(colour = code),cohesin) +
  geom_smooth(method = "lm", se = T) +
  labs(x = "Genes", y = "Residuals") +
  ylim(-10,10) +
  ggtitle("SpEx - Cohesin") +
  geom_hline(yintercept=0,linetype="dashed") +
  theme_DHO +
  scale_color_manual(name = "Residuals", 
                     values = c("lightgrey", "#00CC66"))


plot_exp_2 = ggplot(exp_cohesin, aes(x = Genes, y = Residuals)) +
  geom_point(aes(colour = code),exp_cohesin) +
  geom_smooth(method = "lm", se = T) +
  labs(x = "Genes", y = "Residuals") +
  ylim(-10,10) +
  ggtitle("Baseline") +
  geom_hline(yintercept=0,linetype="dashed") +
  theme_DHO +
  scale_color_manual(name = "Residuals", 
                     values = c("lightgrey", "darkgreen"))

plot_rnap2 = ggplot(rnap2, aes(x = Genes, y = Residuals)) +
  geom_point(aes(colour = code),rnap2) +
  geom_smooth(method = "lm", se = T) +
  labs(x = "Genes", y = "Residuals") +
  ylim(-10,10) +
  ggtitle("SpEx - RNA POL II") +
  geom_hline(yintercept=0,linetype="dashed") +
  theme_DHO +
  scale_color_manual(name = "Residuals", 
                     values = c("lightgrey", "#DB7093"))

plot_exp_3 = ggplot(exp_rnap2, aes(x = Genes, y = Residuals)) +
  geom_point(aes(colour = code),exp_rnap2) +
  geom_smooth(method = "lm", se = T) +
  labs(x = "Genes", y = "Residuals") +
  ylim(-10,10) +
  ggtitle("Baseline") +
  geom_hline(yintercept=0,linetype="dashed") +
  theme_DHO +
  scale_color_manual(name = "Residuals", 
                     values = c("lightgrey", "#990000"))

plot_spex = ggplot(spex, aes(x = Genes, y = Residuals)) +
  geom_point(aes(colour = code),spex) +
  geom_smooth(method = "lm", se = T) +
  labs(x = "Genes", y = "Residuals") +
  ylim(-10,10) +
  ggtitle("SpEx - Best") +
  geom_hline(yintercept=0,linetype="dashed") +
  theme_DHO +
  scale_color_manual(name = "Residuals", 
                     values = c("lightgrey", "#ff8000"))

plot_exp_sx = ggplot(exp_spex, aes(x = Genes, y = Residuals)) +
  geom_point(aes(colour = code),exp_spex) +
  geom_smooth(method = "lm", se = T) +
  labs(x = "Genes", y = "Residuals") +
  ylim(-10,10) +
  ggtitle("Baseline") +
  geom_hline(yintercept=0,linetype="dashed") +
  theme_DHO +
  scale_color_manual(name = "Residuals", 
                     values = c("lightgrey", "#663300"))
#-----------------------------------------------------------------------#

r1_spex$Methord = "SpEx - Best"

exp_1$Methord = "Baseline"
dist_spex = rbind(r1_spex, exp_1)

r1_ctcf$Methord = "SpEx - CTCF"
dist_ctcf = rbind(r1_ctcf, exp_1)

r1_cohesin$Methord = "SpEx - Cohesin"
dist_cohesin = rbind(r1_cohesin, exp_1)

r1_rnap2$Methord = "SpEx - RNA POL II"
dist_rnap2 = rbind(r1_rnap2, exp_1)

plot_dist_spex = ggplot(dist_spex, aes(x=Residuals,fill = Methord)) + 
  geom_density(alpha = 0.25) + ggtitle("Baseline vs SpEx - Best") +  
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

patchwork2 <- (plot_dist_ctcf | plot_ctcf | plot_exp_1) / (plot_dist_cohesin | plot_cohesin | plot_exp_2) / (plot_dist_rnap2 | plot_rnap2 | plot_exp_3) / (plot_dist_spex | plot_spex | plot_exp_sx)
patchwork2 + theme(plot.title = element_text(face = 'bold',size = 15)) & theme(plot.subtitle = element_text(face = 'bold',size = 10,color = "#660000"))
dev.off(dev.list()["RStudioGD"])
#-----------------------------------------------------------------------#

set1 <- as.vector(unlist(spex_sig$Genes))
set2 <- as.vector(unlist(exp_sig$Genes))
set3 <- as.vector(unlist(ctcf_sig$Genes))
set4 <- as.vector(unlist(cohesin_sig$Genes))
set5 <- as.vector(unlist(rnap2_sig$Genes))

Venn1 = venn.diagram(
  x = list(set2, set3, set4, set5),
  category.names = c("Baseline" , "SpEx - CTCF", "SpEx - Cohesin" , "SpEx - RNA POL II"),
  filename = NULL,
  output =T,
  imagetype="png" ,
  compression = "lzw",
  lwd = 1,
  col=c('#D56565',"#5C83AC", '#FF9933', "#5C8000"),
  fill = c(alpha('#D56565',0.25),alpha('#5C83AC',0.25), alpha('#FF9933',0.25),alpha('#5C8000',0.25)),
  cex = 1.5,
  cat.cex = 1.75,
  fontfamily = "Arial",
  sub.fontfamily = "Arial",
  cat.fontfamily ="Arial",
  fontface = "bold",
  height = 700 ,
  width = 700 ,
  resolution = 150,
  main.cex = 6,
  sub.col = "#660000",
  main.fontface = 15,
  sub.fontface = 15,
  disable.logging = T)
grid::grid.draw(Venn1)
#ggsave(Venn1, file="#revision_3.svg", device = "svg", width = 12,height = 12,dpi = 1200,limitsize = F)
dev.off(dev.list()["RStudioGD"])

Venn2 = venn.diagram(
  x = list(set1, set2),
  category.names = c("SpEx - Best" , "Baseline"),
  filename = NULL,
  output =T,
  compression = "lzw",
  lwd = 0.50,
  col=c('#FF9933', "#5C8000"),
  fill = c(alpha('#FF9933',0.25),alpha('#5C8000',0.25)),
  cex = 1.5,
  cat.cex = 1.80,
  catface ="bold",
  fontfamily = "Arial",
  sub.fontfamily = "Arial",
  cat.fontfamily ="Arial",
  fontface = "bold",
  height = 700,
  width = 700 ,
  resolution = 100,
  main.cex = 2,
  sub.col = "#660000",
  main.fontface = 2,
  sub.fontface = 2,
  scaled =3,
  ext.text = 0,
  disable.logging = T)
grid::grid.draw(Venn2)
#ggsave(Venn2, file="#revision_4.svg", device = "svg", width = 12,height = 12,dpi = 1200,limitsize = F)
dev.off(dev.list()["RStudioGD"])
#-----------------------------------------------------------------------#