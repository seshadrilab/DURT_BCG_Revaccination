date <- 20190314

library(openCyto)
library(ggcyto)

gs_b1_controls <- load_gs(here::here("out/GatingSets/20190314_GatingSet_Batch1_Controls"))
gs_b1 <- load_gs(here::here("out/GatingSets/20180908_GatingSet_AllBatches"))
gs131 <- subset(gs_b1, `SAMPLE ID` == "131 V13")
# use Patient 131 at Week 3

pData(parameters(getData(gs131[[1]])))[,c(1,2)]
# name          desc
# $P1      Time          Time
# $P2     FSC-A         FSC-A
# $P3     FSC-H         FSC-H
# $P4     SSC-A         SSC-A
# $P5     SSC-H         SSC-H
# $P6  <V450-A>           MR1
# $P7  <V510-A> Unloaded_CD1b
# $P8  <V610-A>       TRAV1-2
# $P9  <V710-A>          CCR7
# $P10 <V780-A>     CD14/CD19
# $P11 <B515-A>           L/D
# $P12 <B710-A>          CD8b
# $P13 <U395-A>           CD3
# $P14 <U730-A>        CD45RA
# $P15 <R660-A>       MA_R660
# $P16 <R780-A>           CD4
# $P17 <G575-A>       aGalCer
# $P18 <G610-A>       MA_G610
# $P19 <G780-A>            GD
getNodes(gs131)
# "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/aGalCer+"
# "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/GammaDelta"
# "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/MR1+"

attributes(getGate(gs131, "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/aGalCer+")[[1]])$boundaries # vs FSC-A
attributes(getGate(gs131, "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/GammaDelta")[[1]])$boundaries # vs CD3
attributes(getGate(gs131, "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/MR1+")[[1]])$boundaries # vs TRAV1-2

# Fig 2B row 1 (subject samples)

png(filename = sprintf(here::here("out/20190314_PlotsForPoster/%s_Fig2B_pt1.png"), date), width=5, height=5, units = "in", res = 600)
ggcyto(gs131, aes(x = CD3, y = GD), subset = "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes", filter = marginalFilter) +
  coord_cartesian(xlim = c(1700,3600), ylim = c(0,3800)) + 
  ggcyto::labs_cyto("marker") +
  geom_hex(bins=240) +
  geom_gate("/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/GammaDelta") + 
  ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"),
                 legend.position="none", axis.text=ggplot2::element_blank(), axis.ticks=ggplot2::element_blank(),
                 axis.title=element_text(size=26), plot.title = element_text(size = 30, hjust = 0.5, face = "plain"),
                 strip.background = element_blank(), strip.text = element_blank(),
                 line = element_blank(),
                 panel.border = element_rect(fill = NA, color = "black"),
                 panel.spacing = unit(0,"line")) +
  ggplot2::labs(title="Subject Sample", y = bquote("Pan-" ~ gamma * delta ~ " (PE Cy7 G780)"))
dev.off()

png(filename = sprintf(here::here("out/20190314_PlotsForPoster/%s_Fig2B_pt2.png"), date), width=5, height=5, units = "in", res = 600)
ggcyto(gs131, aes(x = `TRAV1-2`, y = MR1), subset = "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes", filter = marginalFilter) +
  coord_cartesian(xlim = c(0,3300), ylim = c(0,3300)) + 
  ggcyto::labs_cyto("marker") +
  geom_hex(bins=240) +
  geom_gate("/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/MR1+") + 
  ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"),
                 legend.position="none", axis.text=ggplot2::element_blank(), axis.ticks=ggplot2::element_blank(),
                 axis.title=element_text(size=28), plot.title = element_text(size = 30, hjust = 0.5, face = "plain"),
                 strip.background = element_blank(), strip.text = element_blank(),
                 line = element_blank(),
                 panel.border = element_rect(fill = NA, color = "black"),
                 panel.spacing = unit(0,"line")) +
  ggplot2::labs(title="Subject Sample", y = "MR1 5-OP-RU (V450)")
dev.off()

png(filename = sprintf(here::here("out/20190314_PlotsForPoster/%s_Fig2B_pt3.png"), date), width=5, height=5, units = "in", res = 600)
ggcyto(gs131, aes(x = `FSC-A`, y = aGalCer), subset = "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes", filter = marginalFilter) +
  ggcyto_par_set(limits = "instrument") +
  ggcyto::labs_cyto("marker") +
  geom_hex(bins=240) +
  geom_gate("/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/aGalCer+") + 
  ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"),
                 legend.position="none", axis.text=ggplot2::element_blank(), axis.ticks=ggplot2::element_blank(),
                 axis.title=element_text(size=28), plot.title = element_text(size = 30, hjust = 0.5, face = "plain"),
                 strip.background = element_blank(), strip.text = element_blank(),
                 line = element_blank(),
                 panel.border = element_rect(fill = NA, color = "black"),
                 panel.spacing = unit(0,"line")) +
  ggplot2::labs(title="Subject Sample", y = "aGalCer CD1d (PE G575)")
dev.off()

##############################

# Fig 2B row 2 (Control samples)
pData(gs_b1_controls)

pData(parameters(getData(subset(gs_b1_controls, `SAMPLE ID` == "MR1 FMO")[[1]])))[,c(1,2)] # here MR1 channel is "FMO"
png(filename = sprintf(here::here("out/20190314_PlotsForPoster/%s_Fig2B_pt5.png"), date), width=5, height=5, units = "in", res = 600)
ggcyto(subset(gs_b1_controls, `SAMPLE ID` == "MR1 FMO"),
       aes(x = `TRAV1-2`, y = FMO), subset = "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes", filter = marginalFilter) +
  coord_cartesian(xlim = c(0,3300), ylim = c(0,3300)) + 
  ggcyto::labs_cyto("marker") +
  geom_hex(bins=240) +
  geom_gate("/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/MR1+") + 
  ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"),
                 legend.position="none", axis.text=ggplot2::element_blank(), axis.ticks=ggplot2::element_blank(),
                 axis.title=element_text(size=28), plot.title = element_text(size = 30, hjust = 0.5, face = "plain"),
                 strip.background = element_blank(), strip.text = element_blank(),
                 line = element_blank(),
                 panel.border = element_rect(fill = NA, color = "black"),
                 panel.spacing = unit(0,"line")) +
  ggplot2::labs(title="MR1 5-OP-RU\nFMO (V450)", y = "MR1 5-OP-RU (V450)")
dev.off()

png(filename = sprintf(here::here("out/20190314_PlotsForPoster/%s_Fig2B_pt6.png"), date), width=5, height=5, units = "in", res = 600)
ggcyto(subset(gs_b1_controls, `SAMPLE ID` == "aGalCer FMO"),
       aes(x = `FSC-A`, y = FMO), subset = "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes", filter = marginalFilter) +
  ggcyto_par_set(limits = "instrument") +
  ggcyto::labs_cyto("marker") +
  geom_hex(bins=240) +
  geom_gate("/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/aGalCer+") + 
  ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"),
                 legend.position="none", axis.text=ggplot2::element_blank(), axis.ticks=ggplot2::element_blank(),
                 axis.title=element_text(size=28), plot.title = element_text(size = 30, hjust = 0.5, face = "plain"),
                 strip.background = element_blank(), strip.text = element_blank(),
                 line = element_blank(),
                 panel.border = element_rect(fill = NA, color = "black"),
                 panel.spacing = unit(0,"line")) +
  ggplot2::labs(title="aGalCer CD1d\nFMO (PE G575)", y = "aGalCer CD1d (PE G575)")
dev.off()
