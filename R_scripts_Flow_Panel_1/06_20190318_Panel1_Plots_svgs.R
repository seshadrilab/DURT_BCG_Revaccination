# Plan:
# Use an individual from Batch 2 (PTID 122) as an example for making the flow plots showing the gating strategy. (Fig 2A)
# Then for Fig 4A and 4B make flow plots showing GMM gating and staining (PTID 122 has strong 3 weeks post-BCG GMM staining)
# Then for Fig 4C make GMM line plot
# Fig 5A make line plot of GEM T cells
# Fig 5B t-SNE showing GEM Memory at pre and 3 wks

date <- 20190314

library(openCyto)
library(ggcyto)
library(svglite)

gs_dir <- here::here("out/GatingSets/20180907_GatingSets_Combined_AutoGate")
gs <- load_gs(gs_dir)
gs122 <- subset(gs, Individual == "122" & TimePoint_v2 == "3 weeks")

pData(parameters(getData(gs122[[1]])))[,c(1,2)]
#          name          desc
# $P1      Time          Time
# $P2     FSC-A         FSC-A
# $P3     FSC-H         FSC-H
# $P4     SSC-A         SSC-A
# $P5     SSC-H         SSC-H
# $P6  <B710-A>          CD8b
# $P7  <B515-A>           L/D
# $P8  <G780-A>        CD45RA
# $P9  <G610-A>           CD3
# $P10 <G575-A> Unloaded CD1b
# $P11 <R780-A>           CD4
# $P12 <R710-A>          CD56
# $P13 <R660-A>      GMM_R660
# $P14 <V780-A>     CD14/CD19
# $P15 <V710-A>          CCR7
# $P16 <V655-A>      GMM_V655
# $P17 <V610-A>        HLA-DR
# $P18 <V510-A>       TRAV1-2
# $P19 <V450-A>          CD38

###############################

# Figure 2A

getNodes(gs122) # "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/GMM_R660+"
attributes(getGate(gs122, "/Live-CD3+")[[1]])$boundaries
svglite(file=sprintf(here::here("out/20190314_PlotsForPoster/%s_Fig2a_pt1.svg"), date), width=4, height=5)
# png(filename = sprintf(here::here("out/20190314_PlotsForPoster/%s_Fig2A_pt1.svg"), date), width=4, height=5)
ggcyto(gs122, aes(x = CD3, y = `L/D`), subset = "root", filter = marginalFilter) +
  ggcyto_par_set(limits = "instrument") +
  ggcyto::labs_cyto("marker") +
  geom_hex(bins=240) +
  geom_gate("/Live-CD3+") + 
  ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"),
                 legend.position="none", axis.text=ggplot2::element_blank(), axis.ticks=ggplot2::element_blank(),
                 axis.title=element_text(size=36), plot.title = element_text(size = 48, hjust = 0.5, face = "plain"),
                 strip.background = element_blank(), strip.text = element_blank(),
                 line = element_blank(),
                 panel.border = element_rect(fill = NA, color = "black"),
                 panel.spacing = unit(0,"line")) +
  ggplot2::labs(title="\nT Cells", y = "Live/Dead")
dev.off()

svglite(file=sprintf(here::here("out/20190314_PlotsForPoster/%s_Fig2A_pt2.svg"), date), width=4, height=5,)# units = "in", res = 600)
ggcyto(gs122, aes(x = CD3, y = `CD14/CD19`), subset = "/Live-CD3+", filter = marginalFilter) +
  ggcyto_par_set(limits = "instrument") +
  ggcyto::labs_cyto("marker") +
  geom_hex(bins=240) +
  geom_gate("/Live-CD3+/CD14-CD19-") + 
  ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"),
                 legend.position="none", axis.text=ggplot2::element_blank(), axis.ticks=ggplot2::element_blank(),
                 axis.title=element_text(size=36), plot.title = element_text(size = 35, hjust = 0.5, face = "plain", margin = margin(t = 16, b = 16)),
                 strip.background = element_blank(), strip.text = element_blank(),
                 line = element_blank(),
                 panel.border = element_rect(fill = NA, color = "black"),
                 panel.spacing = unit(0,"line")) +
  ggplot2::labs(title="Monocyte/B Cell\nExclusion")
dev.off()

svglite(file=sprintf(here::here("out/20190314_PlotsForPoster/%s_Fig2A_pt3.svg"), date), width=4, height=5)
ggcyto(gs122, aes(x = `FSC-A`, y = `FSC-H`), subset = "/Live-CD3+/CD14-CD19-", filter = marginalFilter) +
  ggcyto_par_set(limits = "instrument") +
  ggcyto::labs_cyto("marker") +
  geom_hex(bins=240) +
  geom_gate("/Live-CD3+/CD14-CD19-/S") + 
  ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"),
                 legend.position="none", axis.text=ggplot2::element_blank(), axis.ticks=ggplot2::element_blank(),
                 axis.title=element_text(size=36), plot.title = element_text(size = 48, hjust = 0.5, face = "plain"),
                 strip.background = element_blank(), strip.text = element_blank(),
                 line = element_blank(),
                 panel.border = element_rect(fill = NA, color = "black"),
                 panel.spacing = unit(0,"line")) +
  ggplot2::labs(title="\nSinglet")
dev.off()

svglite(file=sprintf(here::here("out/20190314_PlotsForPoster/%s_Fig2A_pt4.svg"), date), width=4, height=5)
ggcyto(gs122, aes(x = `FSC-A`, y = CCR7), subset = "/Live-CD3+/CD14-CD19-/S", filter = marginalFilter) +
  ggcyto_par_set(limits = "instrument") +
  ggcyto::labs_cyto("marker") +
  geom_hex(bins=240) +
  geom_gate("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper") + 
  ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"),
                 legend.position="none", axis.text=ggplot2::element_blank(), axis.ticks=ggplot2::element_blank(),
                 axis.title=element_text(size=36), plot.title = element_text(size = 48, hjust = 0.5, face = "plain"),
                 strip.background = element_blank(), strip.text = element_blank(),
                 line = element_blank(),
                 panel.border = element_rect(fill = NA, color = "black"),
                 panel.spacing = unit(0,"line")) +
  ggplot2::labs(title="\nKeeper")
dev.off()

svglite(file=sprintf(here::here("out/20190314_PlotsForPoster/%s_Fig2A_pt5.svg"), date), width=4, height=5)
ggcyto(gs122, aes(x = `FSC-A`, y = CD38), subset = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper", filter = marginalFilter) +
  ggcyto_par_set(limits = "instrument") +
  ggcyto::labs_cyto("marker") +
  geom_hex(bins=240) +
  geom_gate("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper") + 
  ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"),
                 legend.position="none", axis.text=ggplot2::element_blank(), axis.ticks=ggplot2::element_blank(),
                 axis.title=element_text(size=36), plot.title = element_text(size = 48, hjust = 0.5, face = "plain"),
                 strip.background = element_blank(), strip.text = element_blank(),
                 line = element_blank(),
                 panel.border = element_rect(fill = NA, color = "black"),
                 panel.spacing = unit(0,"line")) +
  ggplot2::labs(title="\nKeeper")
dev.off()

svglite(file=sprintf(here::here("out/20190314_PlotsForPoster/%s_Fig2A_pt6.svg"), date), width=4, height=5)
ggcyto(gs122, aes(x = `FSC-A`, y = `SSC-A`), subset = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper", filter = marginalFilter) +
  ggcyto_par_set(limits = "instrument") +
  ggcyto::labs_cyto("marker") +
  geom_hex(bins=240) +
  geom_gate("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes") + 
  ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"),
                 legend.position="none", axis.text=ggplot2::element_blank(), axis.ticks=ggplot2::element_blank(),
                 axis.title=element_text(size=36), plot.title = element_text(size = 48, hjust = 0.5, face = "plain"),
                 strip.background = element_blank(), strip.text = element_blank(),
                 line = element_blank(),
                 panel.border = element_rect(fill = NA, color = "black"),
                 panel.spacing = unit(0,"line")) +
  ggplot2::labs(title="\nLymphocyte")
dev.off()

###########################################

# Figure 4A

gsFMO <- load_gs(here::here("out/GatingSets/20180906_GatingSet_Batch2_FMO_AutoGate"))
pData(parameters(getData(gsFMO[[1]])))[,c(1,2)]
# name          desc
# $P1      Time          Time
# $P2     FSC-A         FSC-A
# $P3     FSC-H         FSC-H
# $P4     SSC-A         SSC-A
# $P5     SSC-H         SSC-H
# $P6  <B710-A>          CD8b
# $P7  <B515-A>           L/D
# $P8  <G780-A>        CD45RA
# $P9  <G610-A>           CD3
# $P10 <G575-A> Unloaded CD1b
# $P11 <R780-A>           CD4
# $P12 <R710-A>          CD56
# $P13 <R660-A>      GMM_R660
# $P14 <V780-A>     CD14/CD19
# $P15 <V710-A>          CCR7
# $P16 <V655-A>      GMM_V655
# $P17 <V610-A>        HLA-DR
# $P18 <V510-A>       TRAV1-2
# $P19 <V450-A>          CD38

pData(gsFMO)$facet_labels <- factor(c("No Tetramer", "Mock Loaded\nCD1b FMO", "GMM CD1b FMO\n(BV650)", "GMM CD1b FMO\n(APC R660)", "Full Stain"),
                                    levels = c(c("No Tetramer", "Mock Loaded\nCD1b FMO", "GMM CD1b FMO\n(BV650)", "GMM CD1b FMO\n(APC R660)", "Full Stain")))

svglite(file=sprintf(here::here("out/20190314_PlotsForPoster/%s_Fig4A_pt1.svg"), date), width=8, height=5)
ggcyto(subset(gsFMO, `SAMPLE ID` %in% c("GMM CD1b BV650 FMO", "Full Stain")),
       mapping = aes(x = "Unloaded CD1b", y = "GMM_V655"),
       subset = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes",
       filter = marginalFilter) +
  geom_hex(bins=240) +
  ggcyto::labs_cyto("marker") +
  geom_gate(c("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/GMM_V655+")) +
  facet_grid(. ~ facet_labels) +
  coord_cartesian(xlim = c(-300,1800), ylim = c(-200,3100)) + 
  ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"),
                 legend.position="none", axis.text=ggplot2::element_blank(), axis.ticks=ggplot2::element_blank(),
                 axis.title=element_text(size=33),
                 strip.text = element_text(size=27), strip.background = element_rect(fill="lavender"),
                 line = element_blank(),
                 panel.border = element_rect(fill = NA, color = "black"),
                 panel.spacing = unit(0,"line"),
                 plot.title = element_blank()) +
  labs(y = "GMM CD1b\n(BV650 V655)", x = "Mock Loaded CD1b")
dev.off()

unloaded_gates_to_draw <- sapply(c("28341.fcs_215406", "28335.fcs_212355"), function(sn) { getGate(subset(gsFMO, `SAMPLE ID` %in% c("Full Stain")),
                                                                                                   "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/Unloaded_CD1b+")[[1]] })
R660_gates_to_draw <- sapply(c("28341.fcs_215406", "28339.fcs_215217"), function(sn) { getGate(subset(gsFMO, `SAMPLE ID` %in% c("Full Stain")),
                                                                                               "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/GMM_R660+")[[1]] })
svglite(file=sprintf(here::here("out/20190314_PlotsForPoster/%s_Fig4A_pt2.svg"), date), width=12, height=5)
ggcyto(subset(gsFMO, `SAMPLE ID` %in% c("Unloaded CD1b FMO", "Full Stain", "GMM CD1b APC FMO")), #
       mapping = aes(x = "Unloaded CD1b", y = "GMM_R660"),
       subset = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes",
       filter = marginalFilter) +
  geom_hex(bins=240) +
  ggcyto::labs_cyto("marker") +
  geom_gate(unloaded_gates_to_draw) + geom_gate(R660_gates_to_draw) +
  # geom_gate(c("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/Unloaded_CD1b+", "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/GMM_R660+")) +
  facet_grid(. ~ facet_labels) +
  coord_cartesian(xlim = c(-300,1800), ylim = c(-200,3100)) + 
  ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"),
                 legend.position="none", axis.text=ggplot2::element_blank(), axis.ticks=ggplot2::element_blank(),
                 axis.title=element_text(size=33),
                 strip.text = element_text(size=27), strip.background = element_rect(fill="lavender"),
                 line = element_blank(),
                 panel.border = element_rect(fill = NA, color = "black"),
                 panel.spacing = unit(0,"line"),
                 plot.title = element_blank()) +
  labs(y = "GMM CD1b\n(APC R660)", x = "Mock Loaded CD1b")
dev.off()

svglite(file=sprintf(here::here("out/20190314_PlotsForPoster/%s_Fig4A_pt3.svg"), date), width=5, height=5.05)
ggcyto(subset(gsFMO, `SAMPLE ID` %in% c("Full Stain")),
       mapping = aes(x = "GMM_R660", y = "GMM_V655"),
       subset = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes",
       filter = marginalFilter) +
  geom_hex(bins=240) +
  ggcyto::labs_cyto("marker") +
  geom_gate(c("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/GMM_V655+", "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/GMM_R660+")) +
  facet_grid(. ~ facet_labels) +
  coord_cartesian(xlim = c(-200,3100), ylim = c(-200,3100)) + 
  ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"),
                 legend.position="none", axis.text=ggplot2::element_blank(), axis.ticks=ggplot2::element_blank(),
                 axis.title=element_text(size=33),
                 strip.text = element_text(size=27), strip.background = element_rect(fill="lavender"),
                 line = element_blank(),
                 panel.border = element_rect(fill = NA, color = "black"),
                 panel.spacing = unit(0,"line"),
                 plot.title = element_blank()) +
  labs(y = "GMM CD1b\n(BV650 V655)", x = "GMM CD1b\n(APC R660)")
dev.off()

###########################################

# Figure 4B

svglite(file=sprintf(here::here("out/20190314_PlotsForPoster/%s_Fig4B_pt1_v2.svg"), date), width=5, height=5)
ggcyto(gs122,
       mapping = aes(x = "GMM_R660", y = "GMM_V655"),
       subset = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes",
       filter = marginalFilter) +
  geom_hex(bins=300) +
  ggcyto::labs_cyto("marker") +
  geom_gate(c("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/GMM_V655+", "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/GMM_R660+")) +
  coord_cartesian(xlim = c(-200,3100), ylim = c(-200,3100)) + 
  ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"),
                 legend.position="none", axis.text=ggplot2::element_blank(), axis.ticks=ggplot2::element_blank(),
                 axis.title=element_text(size=28), plot.title = element_text(size = 30, hjust = 0.5, face = "plain"),
                 strip.background = element_blank(), strip.text = element_blank(),
                 line = element_blank(),
                 panel.border = element_rect(fill = NA, color = "black"),
                 panel.spacing = unit(0,"line")) +
  labs(y = "GMM CD1b (BV650 V655)", x = "GMM CD1b (APC R660)", title="Subject Sample") +
  geom_overlay("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/GMM_clean", size = 0.4, color = "red")
dev.off()

recompute(gs122, "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes")
svglite(file=sprintf(here::here("out/20190314_PlotsForPoster/%s_Fig4B_pt2_v2.svg"), date), width=5, height=5)
ggcyto(gs122, mapping = aes(x = "GMM_R660", y = "TRAV1-2"),
       subset = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes",
       filter = marginalFilter) +
  geom_hex(bins=240) +
  ggcyto::labs_cyto("marker") +
  coord_cartesian(xlim = c(-200,3100), ylim = c(-500,3000)) + 
  geom_gate("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/TRAV1-2+") + 
  ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"),
                 legend.position="none", axis.text=ggplot2::element_blank(), axis.ticks=ggplot2::element_blank(),
                 axis.title=element_text(size=28), plot.title = element_text(size = 30, hjust = 0.5, face = "plain"),
                 strip.background = element_blank(), strip.text = element_blank(),
                 line = element_blank(),
                 panel.border = element_rect(fill = NA, color = "black"),
                 panel.spacing = unit(0,"line")) +
  labs(title = "TRAV1-2 and GEM Cells", x = "GMM CD1b (APC R660)") +
  geom_overlay("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/GEM_GMM", size = 0.4, color = "red")
dev.off()
getPopStats(gs122, subpopulations = c("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/GEM_GMM", "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD3_clean/GEM_GMM"))
# name                                                                    Population      Parent Count ParentCount
# 1: Concatenate_122 V13.FCS_2418687           /Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/GEM_GMM Lymphocytes    68      770786
# 2: Concatenate_122 V13.FCS_2418687 /Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD3_clean/GEM_GMM   CD3_clean    68      768976