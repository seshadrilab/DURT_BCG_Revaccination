library(here)
library(openCyto)
library(ggcyto)
library(lubridate)
library(CytoML)

batch <- "Batch 1"

pre_inh_metadata_file <- here::here("data/151127_PreINH_Manifest.xlsx")
bcg_vacc_metadata_file <- here::here("data/150518_BCGVacc_Manifest.xlsx")
pre_inh_metadata <- readxl::read_excel(pre_inh_metadata_file, sheet = 1, col_names = T, n_max = 20)
pre_inh_metadata$Timepoint <- "Pre"
bcg_vacc_metadata <- readxl::read_excel(bcg_vacc_metadata_file, sheet = 1, col_names = T)
all_ext_metadata <- rbind(pre_inh_metadata[, c("Sample ID", "Sample Date", "Timepoint")],
                          bcg_vacc_metadata[, c("Sample ID", "Sample Date", "Timepoint")])
all_ext_metadata$`Sample ID` <- as.character(all_ext_metadata$`Sample ID`)
colnames(all_ext_metadata)[which(colnames(all_ext_metadata) == "Sample Date")] <- "SampleDate"

# workspace_path_b1 <- here::here("data/20180711_SATVI-TBRU_Memory_B1/20180711_SATVI-TBRU_Mem_B1_dl20180906.xml")
workspace_path_b1 <- here::here("data/20180711_SATVI-TBRU_Memory_B1/20180711_SATVI-TBRU_Mem_B1.xml")
fcs_path_b1 <- here::here("data/20180711_SATVI-TBRU_Memory_B1/20180711_SATVI-TBRU_Revacc_Memory_Batch1")
gs_out_dir_b1_krystle <- here::here("out/GatingSets/20180906_GatingSet_Batch1_Krystle")
gs_out_dir_b1_fmo_krystle <- here::here("out/GatingSets/20180906_GatingSet_Batch1_FMO_Krystle")

gs_out_dir_b1_autogate <- here::here("out/GatingSets/20180906_GatingSet_Batch1_AutoGate")
gs_out_dir_b1_fmo_autogate <- here::here("out/GatingSets/20180906_GatingSet_Batch1_FMO_AutoGate")

gating_plots_dir_b1 <- here::here("out/GatingPlots/20180905_Batch1_Revisited")

ws_b1 <- open_flowjo_xml(workspace_path_b1)
tail(plyr::join(getSampleGroups(ws_b1), getSamples(ws_b1), by = "sampleID"))
unlist(getKeywords(ws_b1, "Concatenate_52 V08.FCS")[grep("\\$|P[0-9]", names(getKeywords(ws_b1, "Concatenate_52 V08.FCS")), invert = T, value = T)])
keywords2import=c("SAMPLE ID", "$DATE", "EXPORT TIME", "WELL ID", "EXPERIMENT NAME")

gs_b1 <- parseWorkspace(ws_b1, name = "Concatenated Samples",
                        keywords = keywords2import,
                        path = fcs_path_b1)
pData(gs_b1)$filename <- sapply(rownames(pData(gs_b1)), function(x) {
  gsub("/.+/", "", description(getData(gs_b1[[x]]))$FILENAME)
}, USE.NAMES = F)
pData(gs_b1)$rowname <- rownames(pData(gs_b1))
pData(gs_b1)$Individual <- sapply(strsplit(pData(gs_b1)$`SAMPLE ID`, " "), function(x) { x[[1]] })
pData(gs_b1)$TimePoint <- sapply(strsplit(pData(gs_b1)$`SAMPLE ID`, " "), function(x) { x[[2]] })
gs_b1_new <- dplyr::left_join(pData(gs_b1), all_ext_metadata, by = c("Individual" = "Sample ID", "TimePoint" = "Timepoint"))
rownames(gs_b1_new) <- gs_b1_new$rowname
pData(gs_b1) <- gs_b1_new
pData(gs_b1)$TimePoint_v2 <- factor(plyr::revalue(pData(gs_b1)$TimePoint, c("Pre" = "Pre-INH",
                                                                            "V08" = "Pre-BCG",
                                                                            "V13" = "3 weeks",
                                                                            "V28" = "1 year")),
                                    levels = c("Pre-INH", "Pre-BCG", "3 weeks", "1 year"))
# Add a Day column which sets the Pre-BCG SampleDate as Day 0 and calculates all the other days based on that. Then I can possibly use a linear x-axis.
pre_bcg_indices <- which(pData(gs_b1)$TimePoint_v2 == "Pre-BCG")
pre_bcg_dates <- pData(gs_b1)[pre_bcg_indices, "SampleDate"]
names(pre_bcg_dates) <- pData(gs_b1)[pre_bcg_indices, "Individual"]
pData(gs_b1)$Day <- purrr::map2_dbl(pData(gs_b1)$Individual, pData(gs_b1)$SampleDate,
                                    function(ind, current_date) {
                                      ymd(current_date) - ymd(pre_bcg_dates[[ind]])
                                    })
pData(gs_b1)

# Patient Samples: Set marker names where needed
pData(parameters(getData(gs_b1[[1]])))[,c(1,2)]
gs_b1_markernames <- c("Time", "FSC-A", "FSC-H", "SSC-A", "SSC-H", "CD8b", "L/D", "CD45RA", "CD3", 
                       "Unloaded CD1b", "CD4", "CD56", "GMM_R660", "CD14/CD19", "CCR7", "GMM_V655", 
                       "HLA-DR", "TRAV1-2", "CD38")
names(gs_b1_markernames) <- c("Time", "FSC-A", "FSC-H", "SSC-A", "SSC-H", "<B710-A>", "<B515-A>", 
                              "<G780-A>", "<G610-A>", "<G575-A>", "<R780-A>", "<R710-A>", "<R660-A>", 
                              "<V780-A>", "<V710-A>", "<V655-A>", "<V610-A>", "<V510-A>", "<V450-A>")
markernames(gs_b1) <- gs_b1_markernames
markernames(gs_b1)
pData(parameters(getData(gs_b1[[1]])))[,c(1,2)]

# Save Krystle's version of the gating as a Gating Set for now.
save_gs(gs_b1, gs_out_dir_b1_krystle)

fmo_b1 <- flowjo_to_gatingset(ws_b1, name = "FMO Controls",
                         keywords = keywords2import,
                         path = fcs_path_b1)
pData(fmo_b1)$filename <- sapply(rownames(pData(fmo_b1)), function(x) {
  gsub("/.+/", "", description(getData(fmo_b1[[x]]))$FILENAME)
}, USE.NAMES = F)
pData(fmo_b1)$rowname <- rownames(pData(fmo_b1))
pData(fmo_b1)

# 20201023 remove extra nodes from fmo GS
to_rm_from_fmo <- c("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/Q9: <G575-A>Ð, <R660-A>+", 
          "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/Q10: <G575-A>+, <R660-A>+", 
          "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/Q11: <G575-A>+, <R660-A>Ð", 
          "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/Q12: <G575-A>Ð, <R660-A>Ð", 
          "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/Q13: <G575-A>Ð, <V655-A>+", 
          "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/Q14: <G575-A>+, <V655-A>+", 
          "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/Q15: <G575-A>+, <V655-A>Ð", 
          "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/Q16: <G575-A>Ð, <V655-A>Ð", 
          "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/Q17: <V655-A>Ð, <R660-A>+", 
          "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/Q18: <V655-A>+, <R660-A>+",
          "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/Q19: <V655-A>+, <R660-A>Ð", 
          "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/Q20: <V655-A>Ð, <R660-A>Ð")
for (g in to_rm_from_fmo) { gs_pop_remove(fmo_b1, g)}

# FMO: Set marker names where needed
pData(parameters(getData(fmo_b1[[1]])))[,c(1,2)]
fmo_b1_markernames <- c("Time", "FSC-A", "FSC-H", "SSC-A", "SSC-H", "CD8b", "L/D", "CD45RA", "CD3", 
                        "Unloaded CD1b", "CD4", "CD56", "GMM_R660", "CD14/CD19", "CCR7", "GMM_V655", "HLA-DR", "TRAV1-2", 
                        "CD38")
names(fmo_b1_markernames) <- c("Time", "FSC-A", "FSC-H", "SSC-A", "SSC-H", "<B710-A>", "<B515-A>", 
                               "<G780-A>", "<G610-A>", "<G575-A>", "<R780-A>", "<R710-A>", "<R660-A>", 
                               "<V780-A>", "<V710-A>", "<V655-A>", "<V610-A>", "<V510-A>", "<V450-A>")
markernames(fmo_b1) <- fmo_b1_markernames
markernames(fmo_b1)
pData(parameters(getData(fmo_b1[[1]])))[,c(1,2)]

# Save FMO as a GatingSet for now, with Krystle's gating tree
save_gs(fmo_b1, gs_out_dir_b1_fmo_krystle)

# Lymphocytes gate stats: over 10,000? (Yes)
getPopStats(gs_b1, subpopulation = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes")
# ranges between 226,181 and 842,897 Lymphocyte events

# Make a couple example plots of Krystle's gating
png(filename = file.path(gating_plots_dir_b1, "20180820_Krystle_GatingTree.png"), width = 1000, height = 1000)
plot(gs_b1, fontsize = 20)
dev.off()

png(filename = file.path(gating_plots_dir_b1, "20180820_Krystle_plotGate_Example.png"), width = 1000, height = 1000)
plotGate(gs_b1[[1]])
dev.off()

##################################

# Start autogating

# Remove Krystle's gates

nodes2rm <- c("/CD3-/Q1: CD38Ð, HLA-DR+", "/CD3-/Q2: CD38+, HLA-DR+", 
              "/CD3-/Q3: CD38+, HLA-DRÐ", "/CD3-/Q4: CD38Ð, HLA-DRÐ", "/Live-CD3+/CD14-CD19-/56+",
              "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/Q1: <V655-A>Ð, TRAV1-2+", 
              "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/Q2: <V655-A>+, TRAV1-2+", 
              "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/Q3: <V655-A>+, TRAV1-2Ð", 
              "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/Q4: <V655-A>Ð, TRAV1-2Ð", 
              "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/Q9: <G575-A>Ð, <R660-A>+", 
              "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/Q10: <G575-A>+, <R660-A>+", 
              "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/Q11: <G575-A>+, <R660-A>Ð", 
              "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/Q12: <G575-A>Ð, <R660-A>Ð", 
              "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/Q13: <G575-A>Ð, <V655-A>+", 
              "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/Q14: <G575-A>+, <V655-A>+", 
              "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/Q15: <G575-A>+, <V655-A>Ð", 
              "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/Q16: <G575-A>Ð, <V655-A>Ð", 
              "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/Q17: <V655-A>Ð, <R660-A>+", 
              "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/Q18: <V655-A>+, <R660-A>+",
              "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/Q19: <V655-A>+, <R660-A>Ð", 
              "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/Q20: <V655-A>Ð, <R660-A>Ð", 
              "/Q1: CD56Ð, CD3+", "/Q2: CD56+, CD3+", "/Q3: CD56+, CD3Ð", 
              "/Q4: CD56Ð, CD3Ð", "/S")
for(node in nodes2rm) {
  Rm(node, gs_b1)
}
plot(gs_b1, fontsize = 25, bool = T)

# Now the basic nodes are:
# c("root", "/CD3-", "/Live-CD3+", "/Live-CD3+/CD14-CD19-", "/Live-CD3+/CD14-CD19-/S", 
#   "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper", "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper", 
#   "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes")

# Gating Strategy:
# CD56 gated on all ungated cells or maybe singlet, live cells. 
# Gate HLA-DR and CD38 on CD3- as a quadrant gate.
# CCR7 and CD45RA on CD3- under the Lymphocytes gate.
# TRAV1-2 as a mindensity gate against FSC-A.

# Time
png(filename = file.path(gating_plots_dir_b1,
                         sprintf("%s_%s.png", "Time_vs_FSC-A",
                                 format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))),
    width = 1900, height = 1030)
ggcyto(gs_b1, mapping = aes(x = "Time", y = "FSC-A"),
       subset = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes") +
  geom_hex(bins=120) +
  ggcyto_par_set(limits = "instrument") +
  facet_grid(Individual ~ TimePoint)
dev.off()
png(filename = file.path(gating_plots_dir_b1,
                         sprintf("%s_%s.png", "Time_vs_CD56",
                                 format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))),
    width = 1900, height = 1030)
ggcyto(gs_b1, mapping = aes(x = "Time", y = "CD56"),
       subset = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes") +
  geom_hex(bins=120) +
  ggcyto_par_set(limits = "instrument") +
  facet_grid(Individual ~ TimePoint)
dev.off()
png(filename = file.path(gating_plots_dir_b1,
                         sprintf("%s_%s.png", "Time_vs_CD4",
                                 format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))),
    width = 1900, height = 1030)
ggcyto(gs_b1, mapping = aes(x = "Time", y = "CD4"),
       subset = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes") +
  geom_hex(bins=120) +
  ggcyto_par_set(limits = "instrument") +
  facet_grid(Individual ~ TimePoint)
dev.off()
png(filename = file.path(gating_plots_dir_b1,
                         sprintf("%s_%s.png", "Time_vs_CD8b",
                                 format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))),
    width = 1900, height = 1030)
ggcyto(gs_b1, mapping = aes(x = "Time", y = "CD8b"),
       subset = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes") +
  geom_hex(bins=120) +
  ggcyto_par_set(limits = "instrument") +
  facet_grid(Individual ~ TimePoint)
dev.off()
png(filename = file.path(gating_plots_dir_b1,
                         sprintf("%s_%s.png", "Time_vs_GMM_V655",
                                 format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))),
    width = 1900, height = 1030)
ggcyto(gs_b1, mapping = aes(x = "Time", y = "GMM_V655"),
       subset = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes") +
  geom_hex(bins=120) +
  ggcyto_par_set(limits = "instrument") +
  facet_grid(Individual ~ TimePoint)
dev.off()
png(filename = file.path(gating_plots_dir_b1,
                         sprintf("%s_%s.png", "Time_vs_GMM_R660",
                                 format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))),
    width = 1900, height = 1030)
ggcyto(gs_b1, mapping = aes(x = "Time", y = "GMM_R660"),
       subset = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes") +
  geom_hex(bins=120) +
  ggcyto_par_set(limits = "instrument") +
  facet_grid(Individual ~ TimePoint)
dev.off()
png(filename = file.path(gating_plots_dir_b1,
                         sprintf("%s_%s.png", "Time_vs_Unloaded_CD1b",
                                 format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))),
    width = 1900, height = 1030)
ggcyto(gs_b1, mapping = aes(x = "Time", y = "Unloaded CD1b"),
       subset = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes") +
  geom_hex(bins=120) +
  ggcyto_par_set(limits = "instrument") +
  facet_grid(Individual ~ TimePoint)
dev.off()
png(filename = file.path(gating_plots_dir_b1,
                         sprintf("%s_%s.png", "Time_vs_CD3",
                                 format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))),
    width = 1900, height = 1030)
ggcyto(gs_b1, mapping = aes(x = "Time", y = "CD3"),
       subset = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes") +
  geom_hex(bins=120) +
  ggcyto_par_set(limits = "instrument") +
  facet_grid(Individual ~ TimePoint)
dev.off()

#######################

# Gating
getNodes(gs_b1)

# CD56
ggcyto(gs_b1[[1]], mapping = aes(x = "CD3", y = "CD56"), subset = "root") + geom_hex(bins=120) +
  ggcyto_par_set(limits = "instrument")
ggcyto(gs_b1[[1]], mapping = aes(x = "CD3", y = "CD56"), subset = "CD3-") + geom_hex(bins=120) +
  ggcyto_par_set(limits = "instrument")

# using root instead of CD3- because the CD3+ CD56+ population is slightly higher than the CD3- CD56+ along the CD56 axis, it looks like
add_pop(gs_b1, alias = "CD56+", pop = "+", parent = "root",
        dims = "CD56", gating_method = "mindensity", gating_args = "gate_range=c(975,1500)",
        mc.cores = 14, parallel_type = "multicore", collapseDataForGating = T, groupBy = "Individual")
gc()

png(filename = file.path(gating_plots_dir_b1,
                         sprintf("%s_gate_vs_CD3_%s.png", "CD56+",
                                 format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))),
    width = 1900, height = 1030)
ggcyto(gs_b1, aes(x = "CD3", y = "CD56"), subset = "root") + geom_hex(bins=120) +
  ggcyto_par_set(limits = "instrument") + geom_gate("/CD56+") +
  facet_grid(Individual ~ TimePoint) + labs(title = "/CD56+   collapsed by Individual", caption = batch)
dev.off()

gc()

# HLA-DR and CD38
ggcyto(gs_b1[[1]], mapping = aes(x = "CD38", y = "HLA-DR"), subset = "/CD3-") + geom_hex(bins=120) +
  ggcyto_par_set(limits = "instrument")

add_pop(gs_b1, alias = "CD38+", pop = "+", parent = "/CD3-",
        dims = "CD38", gating_method = "mindensity2", gating_args = "min=1400,max=1850",
        mc.cores = 14, parallel_type = "multicore", collapseDataForGating = T, groupBy = "Individual")

add_pop(gs_b1, alias = "HLA-DR+", pop = "+", parent = "/CD3-",
        dims = "HLA-DR", gating_method = "mindensity2", gating_args = "min=900,max=1700",
        mc.cores = 14, parallel_type = "multicore", collapseDataForGating = T, groupBy = "Individual")

png(filename = file.path(gating_plots_dir_b1,
                         sprintf("%s_gate_%s.png", "CD38_and_HLA-DR",
                                 format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))),
    width = 1900, height = 1030)
ggcyto(gs_b1, mapping = aes(x = "CD38", y = "HLA-DR"), subset = "/CD3-") + geom_hex(bins=120) +
  ggcyto_par_set(limits = "instrument") + geom_gate(c("/CD3-/CD38+", "/CD3-/HLA-DR+")) +
  labs(title = "CD38 and HLA-DR   collapsed by Individual", caption = batch) + facet_grid(Individual ~ TimePoint)
dev.off()

gc()

# Copy the HLA-DR and CD338 gates over to the root node
add_pop(gs_b1, alias = "CD38+", pop = "+", parent = "root",
        dims = "CD38", gating_method = "refGate", gating_args = "/CD3-/CD38+",
        mc.cores = 14, parallel_type = "multicore")
add_pop(gs_b1, alias = "HLA-DR+", pop = "+", parent = "root",
        dims = "HLA-DR", gating_method = "refGate", gating_args = "/CD3-/HLA-DR+",
        mc.cores = 14, parallel_type = "multicore")

png(filename = file.path(gating_plots_dir_b1,
                         sprintf("%s_gate_on_Lymphocytes_%s.png", "CD38_and_HLA-DR",
                                 format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))),
    width = 1900, height = 1030)
ggcyto(gs_b1, mapping = aes(x = "CD38", y = "HLA-DR"), subset = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes") + geom_hex(bins=120) +
  ggcyto_par_set(limits = "instrument") + geom_gate(c("/CD38+", "/HLA-DR+")) +
  labs(title = "CD38 and HLA-DR    collapsed by Individual", caption = batch) + facet_grid(Individual ~ TimePoint)
dev.off()

# TRAV1-2
ggcyto(gs_b1[[1]], mapping = aes(x = "FSC-A", y = "TRAV1-2"),
       subset = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes") +
  geom_hex(bins=120) +
  ggcyto_par_set(limits = "instrument")

add_pop(gs_b1, alias = "TRAV1-2+", pop = "+", parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes",
        dims = "TRAV1-2", gating_method = "mindensity2", gating_args = "min=1000,max=2000",
        mc.cores = 14, parallel_type = "multicore", collapseDataForGating = T, groupBy = "Individual")

png(filename = file.path(gating_plots_dir_b1,
                         sprintf("%s_gate_%s.png", "TRAV1-2",
                                 format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))),
    width = 1900, height = 1030)
ggcyto(gs_b1, mapping = aes(x = "FSC-A", y = "TRAV1-2"),
       subset = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes") +
  geom_hex(bins=120) +
  ggcyto_par_set(limits = "instrument") +
  facet_grid(Individual ~ TimePoint) +
  geom_gate("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/TRAV1-2+") +
  labs(title = "TRAV1-2   collapsed by Individual", caption = batch) + facet_grid(Individual ~ TimePoint)
dev.off()

# CD4 vs CD8
# Note also that we are using CD8b here, which is more rigorous than CD8a (which is what we previously used for all the panels). CD8a on its own is sometimes expressed on cells as an activation marker, and doesn't necessarily appear on t cells. So the fact that we are using CD8b is good.
ggcyto(gs_b1[[1]], mapping = aes(x = "CD8b", y = "CD4"),
       subset = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes") + geom_hex(bins=120)

# write a function that will mimic quadGate.seq i.e. 3 mindensity gates

# alias dim1_alias and dim2_alias shouldn't contain + character
quadGate_seq_mimicer <- function(my_gs, parent, dim1, dim2, dim1_gating_args, dim2_gating_args, dim1_alias = dim1, dim2_alias = dim2) {
  add_pop(my_gs, alias = dim1_alias, pop = "+/-", parent = parent,
          dims = dim1, gating_method = "mindensity2", gating_args = dim1_gating_args,
          mc.cores = 14, parallel_type = "multicore", collapseDataForGating = T, groupBy = "Individual")
  
  dim1_pos <- paste0(parent, "/", dim1_alias, "+")
  dim1_neg <- paste0(parent, "/", dim1_alias, "-")
  
  add_pop(my_gs, alias = dim2_alias, pop = "+/-", parent = dim1_pos,
          dims = dim2, gating_method = "mindensity2", gating_args = dim2_gating_args,
          mc.cores = 14, parallel_type = "multicore", collapseDataForGating = T, groupBy = "Individual")
  
  add_pop(my_gs, alias = dim2_alias, pop = "+/-", parent = dim1_neg,
          dims = dim2, gating_method = "mindensity2", gating_args = dim2_gating_args,
          mc.cores = 14, parallel_type = "multicore", collapseDataForGating = T, groupBy = "Individual")
}

quadGate_seq_mimicer(my_gs = gs_b1, parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes",
                     dim1 = "CD4", dim2 = "CD8b",
                     dim1_gating_args = "min=1000,max=2800",
                     dim2_gating_args = "min=500,max=2000")

plot(gs_b1, fontsize = 25, bool = T)

png(filename = file.path(gating_plots_dir_b1,
                         sprintf("%s_gate_%s.png", "CD4_and_CD8",
                                 format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))),
    width = 1900, height = 1030)
ggcyto(gs_b1, mapping = aes(x = "CD4", y = "CD8b"),
       subset = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes") +
  geom_hex(bins=120) +
  geom_gate(c("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+",
              "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+/CD8b+",
              "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4-/CD8b+")) +
  facet_grid(Individual ~ TimePoint) +
  geom_stats(type = c("gate_name", "percent"), adjust = 0.5, fill = "white", alpha = 0.7) +
  xlim(c(-2000,3500)) + ylim(-1250,4000) +
  labs(title = "CD4 and CD8   collapsed by Individual", caption = batch) + facet_grid(Individual ~ TimePoint)
dev.off()

# Make boolean gates to identify the four quadrants (using the CD8b- gates seems to fail)
flowWorkspace::add(gs_b1, eval(substitute(flowWorkspace::booleanFilter(v),
                                          list(v = as.symbol("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+&/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+/CD8b+")))),
                   parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes", name = "CD4+ CD8+")
flowWorkspace::add(gs_b1, eval(substitute(flowWorkspace::booleanFilter(v),
                                          list(v = as.symbol("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+&!/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+/CD8b+")))),
                   parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes", name = "CD4+ CD8-")
flowWorkspace::add(gs_b1, eval(substitute(flowWorkspace::booleanFilter(v),
                                          list(v = as.symbol("!/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+&!/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4-/CD8b+")))),
                   parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes", name = "CD4- CD8-")
flowWorkspace::add(gs_b1, eval(substitute(flowWorkspace::booleanFilter(v),
                                          list(v = as.symbol("!/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+&/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4-/CD8b+")))),
                   parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes", name = "CD4- CD8+")
flowWorkspace::recompute(gs_b1, "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes")

getNodes(gs_b1)

plot(gs_b1, fontsize = 25, bool = T)

# GMM
# using FMOs
# Strategy: Gate the FMO (draw a single gate that works for all the FMO samples) and then copy it over to the clinical samples and visualize the gates to make sure they look legit.
# In order to gate the FMOs, what I want to do is draw a gate on the "Full stain sample" using the other FMO samples as the negative control depending on which gate I am drawing.
pData(parameters(getData(fmo_b1[[1]])))[,c(1,2)]
# For the GMM CD1b APC (or GMM_R660 / R660-A) gate, use the "GMM CD1b APC FMO" sample as the negative control for that channel
# For the GMM CD1b BV650 (or GMM_V655 / V655-A) gate, use the "GMM CD1b BV650 FMO" sample as the negative control for that channel
# For the Unloaded CD1b (or Unloaded CD1b / G575-A) gate, use the "Unloaded CD1b FMO" sample as the negative control for that channel. Should I be generous with this gate (i.e. draw a low gate in order to identify as many unloaded-cd1b-reactive cells as possible for removal?)

# Use tailgate on the "GMM CD1b APC FMO" sample and copy it over to the full stain sample
add_pop(subset(fmo_b1, `SAMPLE ID` %in% c("GMM CD1b APC FMO")),
        pop = "+", alias = "GMM_R660+", parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes",
        dims = "GMM_R660", gating_method = "tailgate",
        gating_args = "tol=0.1,bias=650", collapseDataForGating = F) # 
plotGate(subset(fmo_b1, `SAMPLE ID` %in% c("GMM CD1b APC FMO")),
         "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/GMM_R660+", xbins = 120)
# Copy the gate to the full stain sample
gmm_r660_gate <- getGate(subset(fmo_b1, `SAMPLE ID` %in% c("GMM CD1b APC FMO"))[[1]],
                         "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/GMM_R660+")
add(subset(fmo_b1, `SAMPLE ID` %in% c("Full Stain")), gmm_r660_gate,
    name = "GMM_R660+", parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes")
recompute(subset(fmo_b1, `SAMPLE ID` %in% c("GMM CD1b APC FMO", "Full Stain")),
          "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes")
plotGate(subset(fmo_b1, `SAMPLE ID` %in% c("GMM CD1b APC FMO", "Full Stain")),
         "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/GMM_R660+", xbins = 120)

# Rm("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/GMM_R660+",
#    subset(fmo_b1, `SAMPLE ID` %in% c("GMM CD1b APC FMO", "Full Stain")))

png(filename = file.path(gating_plots_dir_b1,
                         sprintf("%s_gate_%s.png", "GMM_R660_FMO",
                                 format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))),
    width = 600, height = 400)
ggcyto(subset(fmo_b1, `SAMPLE ID` %in% c("GMM CD1b APC FMO", "Full Stain")),
       mapping = aes(x = "Unloaded CD1b", y = "GMM_R660"),
       subset = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes") +
  geom_hex(bins=240) +
  geom_gate(c("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/GMM_R660+")) +
  facet_grid(. ~ `SAMPLE ID`) +
  geom_stats(type = c("gate_name", "percent"), adjust = 0.7, fill = "white", alpha = 0.7) +
  xlim(c(-300,1800)) + ylim(-1000,3250) +
  labs(caption = batch)
dev.off()

# Now copy the gate to the clinical samples
add(gs_b1, gmm_r660_gate,
    name = "GMM_R660+", parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes")
recompute(gs_b1,
          "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes")
# And save gating plot
png(filename = file.path(gating_plots_dir_b1,
                         sprintf("%s_gate_%s.png", "GMM_R660_Samples",
                                 format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))),
    width = 1900, height = 1030)
ggcyto(gs_b1,
       mapping = aes(x = "Unloaded CD1b", y = "GMM_R660"),
       subset = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes") +
  geom_hex(bins=240) +
  geom_gate(c("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/GMM_R660+")) +
  facet_grid(Individual ~ TimePoint) +
  geom_stats(type = c("gate_name", "percent"), adjust = 0.7, fill = "white", alpha = 0.7) +
  xlim(c(0,1500)) + ylim(0,3250) +
  labs(caption = batch)
dev.off()

# Use tailgate on the "GMM CD1b APC FMO" sample and copy it over to the full stain sample
add_pop(subset(fmo_b1, `SAMPLE ID` %in% c("GMM CD1b BV650 FMO")),
        pop = "+", alias = "GMM_V655+", parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes",
        dims = "GMM_V655", gating_method = "tailgate",
        gating_args = "tol=0.1,bias=700", collapseDataForGating = F) # 
plotGate(subset(fmo_b1, `SAMPLE ID` %in% c("GMM CD1b BV650 FMO")),
         "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/GMM_V655+", xbins = 120)
# Copy the gate to the full stain sample
gmm_v655_gate <- getGate(subset(fmo_b1, `SAMPLE ID` %in% c("GMM CD1b BV650 FMO"))[[1]],
                         "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/GMM_V655+")
add(subset(fmo_b1, `SAMPLE ID` %in% c("Full Stain")), gmm_v655_gate,
    name = "GMM_V655+", parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes")
recompute(subset(fmo_b1, `SAMPLE ID` %in% c("GMM CD1b BV650 FMO", "Full Stain")),
          "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes")
plotGate(subset(fmo_b1, `SAMPLE ID` %in% c("GMM CD1b BV650 FMO", "Full Stain")),
         "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/GMM_V655+", xbins = 120)

# Rm("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/GMM_V655+",
#    subset(fmo_b1, `SAMPLE ID` %in% c("GMM CD1b BV650 FMO", "Full Stain")))

png(filename = file.path(gating_plots_dir_b1,
                         sprintf("%s_gate_%s.png", "GMM_V655_FMO",
                                 format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))),
    width = 600, height = 400)
ggcyto(subset(fmo_b1, `SAMPLE ID` %in% c("GMM CD1b BV650 FMO", "Full Stain")),
       mapping = aes(x = "Unloaded CD1b", y = "GMM_V655"),
       subset = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes") +
  geom_hex(bins=240) +
  geom_gate(c("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/GMM_V655+")) +
  facet_grid(. ~ `SAMPLE ID`) +
  geom_stats(type = c("gate_name", "percent"), adjust = 0.7, fill = "white", alpha = 0.7) +
  xlim(c(-300,1800)) + ylim(-1000,3250) +
  labs(caption = batch)
dev.off()

# Now copy the gate to the clinical samples
add(gs_b1, gmm_v655_gate,
    name = "GMM_V655+", parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes")
recompute(gs_b1,
          "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes")
# And save gating plot
png(filename = file.path(gating_plots_dir_b1,
                         sprintf("%s_gate_%s.png", "GMM_V655_Samples",
                                 format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))),
    width = 1900, height = 1030)
ggcyto(gs_b1,
       mapping = aes(x = "Unloaded CD1b", y = "GMM_V655"),
       subset = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes") +
  geom_hex(bins=240) +
  geom_gate(c("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/GMM_V655+")) +
  facet_grid(Individual ~ TimePoint) +
  geom_stats(type = c("gate_name", "percent"), adjust = 0.7, fill = "white", alpha = 0.7) +
  xlim(c(0,1500)) + ylim(0,3250) +
  labs(caption = batch)
dev.off()



# Use tailgate on the "Unloaded CD1b FMO" sample and copy it over to the full stain sample
# I chose to use a generous gate (i.e. conservative low gate to remove false events)
add_pop(subset(fmo_b1, `SAMPLE ID` %in% c("Unloaded CD1b FMO")),
        pop = "+", alias = "Unloaded_CD1b+", parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes",
        dims = "Unloaded CD1b", gating_method = "tailgate",
        gating_args = "tol=0.1,bias=400", collapseDataForGating = F) # 
plotGate(subset(fmo_b1, `SAMPLE ID` %in% c("Unloaded CD1b FMO")),
         "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/Unloaded_CD1b+", xbins = 120)
# Copy the gate to the full stain sample
gmm_unloaded_gate <- getGate(subset(fmo_b1, `SAMPLE ID` %in% c("Unloaded CD1b FMO"))[[1]],
                             "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/Unloaded_CD1b+")
add(subset(fmo_b1, `SAMPLE ID` %in% c("Full Stain")), gmm_unloaded_gate,
    name = "Unloaded_CD1b+", parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes")
recompute(subset(fmo_b1, `SAMPLE ID` %in% c("Unloaded CD1b FMO", "Full Stain")),
          "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes")
plotGate(subset(fmo_b1, `SAMPLE ID` %in% c("Unloaded CD1b FMO", "Full Stain")),
         "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/Unloaded_CD1b+", xbins = 120)

# Rm("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/Unloaded_CD1b+",
#    subset(fmo_b1, `SAMPLE ID` %in% c("Unloaded CD1b FMO", "Full Stain")))

png(filename = file.path(gating_plots_dir_b1,
                         sprintf("%s_gate_%s.png", "Unloaded_CD1b_FMO",
                                 format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))),
    width = 600, height = 400)
ggcyto(subset(fmo_b1, `SAMPLE ID` %in% c("Unloaded CD1b FMO", "Full Stain")),
       mapping = aes(x = "Unloaded CD1b", y = "GMM_V655"),
       subset = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes") +
  geom_hex(bins=240) +
  geom_gate(c("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/Unloaded_CD1b+")) +
  facet_grid(. ~ `SAMPLE ID`) +
  geom_stats(type = c("gate_name", "percent"), adjust = 0.7, fill = "white", alpha = 0.7) +
  xlim(c(-300,1800)) + ylim(-1000,3250) +
  labs(caption = batch)
dev.off()

# Now copy the gate to the clinical samples
add(gs_b1, gmm_unloaded_gate,
    name = "Unloaded_CD1b+", parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes")
recompute(gs_b1,
          "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes")
# And save gating plot
png(filename = file.path(gating_plots_dir_b1,
                         sprintf("%s_gate_%s.png", "Unloaded_CD1b_Samples",
                                 format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))),
    width = 1900, height = 1030)
ggcyto(gs_b1,
       mapping = aes(x = "Unloaded CD1b", y = "GMM_V655"),
       subset = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes") +
  geom_hex(bins=240) +
  geom_gate(c("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/Unloaded_CD1b+")) +
  facet_grid(Individual ~ TimePoint) +
  geom_stats(type = c("gate_name", "percent"), adjust = 0.6, fill = "white", alpha = 0.7) +
  xlim(c(0,1500)) + ylim(0,3250) +
  labs(caption = batch)
dev.off()

flowWorkspace::add(gs_b1, eval(substitute(flowWorkspace::booleanFilter(v),
                                          list(v = as.symbol("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/GMM_R660+&/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/GMM_V655+&!/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/Unloaded_CD1b+")))),
                   parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes", name = "GMM_clean")
flowWorkspace::recompute(gs_b1, "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes")

ggcyto(gs_b1,
       mapping = aes(x = "GMM_R660", y = "GMM_V655"),
       subset = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes") +
  geom_overlay("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/GMM_clean",
               size = 0.3, alpha = 0.7, color = "red") +
  geom_hex(bins=240) +
  facet_grid(Individual ~ TimePoint) +
  xlim(c(-100,2500)) + ylim(-100,3250)

png(filename = file.path(gating_plots_dir_b1,
                         sprintf("%s_%s.png", "GMM_clean_on_GMM_R660_vs_GMM_V655",
                                 format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))),
    width = 1900, height = 1030)
ggcyto(gs_b1,
       mapping = aes(x = "GMM_R660", y = "GMM_V655"),
       subset = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes") +
  geom_overlay("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/GMM_clean",
               size = 0.3, alpha = 0.7, color = "red") +
  geom_hex(bins=240) +
  facet_grid(Individual ~ TimePoint) +
  xlim(c(-100,2500)) + ylim(-100,3250)
dev.off()

# Overlay the GMM_clean population on top of Lymphocytes along the Unloaded vs loaded channels
png(filename = file.path(gating_plots_dir_b1,
                         sprintf("%s_%s.png", "GMM_clean_on_Unloaded_vs_GMM_V655",
                                 format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))),
    width = 1900, height = 1030)
ggcyto(gs_b1,
       mapping = aes(x = "Unloaded CD1b", y = "GMM_V655"),
       subset = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes") +
  geom_overlay("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/GMM_clean",
               size = 0.3, alpha = 0.7, color = "red") +
  geom_hex(bins=240) +
  facet_grid(Individual ~ TimePoint) +
  xlim(c(-100,2500)) + ylim(-100,3250) +
  labs(title = "GMM_clean on Lymphocytes", caption = batch) +
  geom_gate(c("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/Unloaded_CD1b+"))
dev.off()

png(filename = file.path(gating_plots_dir_b1,
                         sprintf("%s_%s.png", "GMM_clean_on_Unloaded_vs_GMM_R660",
                                 format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))),
    width = 1900, height = 1030)
ggcyto(gs_b1,
       mapping = aes(x = "Unloaded CD1b", y = "GMM_R660"),
       subset = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes") +
  geom_overlay("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/GMM_clean",
               size = 0.3, alpha = 0.7, color = "red") +
  geom_hex(bins=240) +
  facet_grid(Individual ~ TimePoint) +
  xlim(c(-100,2500)) + ylim(-100,3250) +
  labs(title = "GMM_clean on Lymphocytes", caption = batch) +
  geom_gate(c("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/Unloaded_CD1b+"))
dev.off()

# CD45RA and CCR7

ggcyto(gs_b1, mapping = aes(x = "CCR7", y = "CD45RA"),
       subset = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+ CD8-") +
  geom_hex(bins=120) +
  ggcyto_par_set(limits = "instrument")
# For CD4+ CD8 I think I should use two mindensity gates (no quadGate.seq-type gate necessary)

ggcyto(gs_b1, mapping = aes(x = "CCR7", y = "CD45RA"),
       subset = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8+") +
  geom_hex(bins=120) +
  ggcyto_par_set(limits = "instrument")
# Use a mindensity gate along CD45RA and then two more along CCR7 (don't use the quadGate mimicer cause this is finicky)

ggcyto(gs_b1, mapping = aes(x = "CCR7", y = "CD45RA"),
       subset = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8-") +
  geom_hex(bins=120) +
  ggcyto_par_set(limits = "instrument")
# Use a mindensity gate along CD45RA and then two more along CCR7 (don't use the quadGate mimicer cause this is finicky)

ggcyto(gs_b1, mapping = aes(x = "CCR7", y = "CD45RA"),
       subset = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+ CD8+") +
  geom_hex(bins=180) +
  ggcyto_par_set(limits = "instrument")
# Copy over from the CD4- CD8- population

# CD4+ CD8-
add_pop(gs_b1, alias = "CD45RA+", pop = "+", parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+ CD8-",
        dims = "CD45RA", gating_method = "mindensity2", gating_args = "min=1100,max=2200",
        mc.cores = 14, parallel_type = "multicore", collapseDataForGating = T, groupBy = "Individual")
add_pop(gs_b1, alias = "CCR7+", pop = "+", parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+ CD8-",
        dims = "CCR7", gating_method = "mindensity2", gating_args = "min=900,max=1100",
        mc.cores = 14, parallel_type = "multicore", collapseDataForGating = T, groupBy = "Individual")
png(filename = file.path(gating_plots_dir_b1,
                         sprintf("%s_gate_%s.png", "CCR7_and_CD45RA_on_CD4+CD8-",
                                 format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))),
    width = 1900, height = 1030)
ggcyto(gs_b1, mapping = aes(x = "CCR7", y = "CD45RA"),
       subset = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+ CD8-") +
  geom_hex(bins=120) +
  geom_gate(c("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+ CD8-/CD45RA+",
              "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+ CD8-/CCR7+")) +
  facet_grid(Individual ~ TimePoint) +
  geom_stats(type = c("gate_name", "percent"), adjust = 0.5, fill = "white", alpha = 0.7) +
  ggcyto_par_set(limits = "instrument") +
  labs(title = "CD45RA and CCR7 on CD4+ CD8-   collapsed by Individual", caption = batch) + facet_grid(Individual ~ TimePoint)
dev.off()

# CD4- CD8+
add_pop(gs_b1, alias = "CD45RA", pop = "+/-", parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8+",
        dims = "CD45RA", gating_method = "mindensity2", gating_args = "gate_range=c(1100,1500)",
        mc.cores = 14, parallel_type = "multicore", collapseDataForGating = T, groupBy = "Individual")
add_pop(gs_b1, alias = "CCR7+", pop = "+", parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8+/CD45RA+",
        dims = "CCR7", gating_method = "mindensity2", gating_args = "gate_range=c(1350,1550)",
        mc.cores = 14, parallel_type = "multicore", collapseDataForGating = T, groupBy = "Individual")
add_pop(gs_b1, alias = "CCR7+", pop = "+", parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8+/CD45RA-",
        dims = "CCR7", gating_method = "mindensity2", gating_args = "gate_range=c(1800,2100)",
        mc.cores = 14, parallel_type = "multicore", collapseDataForGating = T, groupBy = "Individual")
png(filename = file.path(gating_plots_dir_b1,
                         sprintf("%s_gate_%s.png", "CCR7_and_CD45RA_on_CD4-CD8+",
                                 format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))),
    width = 1900, height = 1030)
ggcyto(gs_b1, mapping = aes(x = "CCR7", y = "CD45RA"),
       subset = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8+") +
  geom_hex(bins=120) +
  geom_gate(c("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8+/CD45RA+",
              "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8+/CD45RA+/CCR7+",
              "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8+/CD45RA-/CCR7+")) +
  facet_grid(Individual ~ TimePoint) +
  geom_stats(type = c("gate_name", "percent"), adjust = 0.5, fill = "white", alpha = 0.7) +
  ggcyto_par_set(limits = "instrument") +
  labs(title = "CD45RA and CCR7 on CD4- CD8+   collapsed by Individual", caption = batch) + facet_grid(Individual ~ TimePoint)
dev.off()

# CD4- CD8-
add_pop(gs_b1, alias = "CD45RA", pop = "+/-", parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8-",
        dims = "CD45RA", gating_method = "mindensity2", gating_args = "gate_range=c(1400,2150)",
        mc.cores = 14, parallel_type = "multicore", collapseDataForGating = T, groupBy = "Individual")
add_pop(gs_b1, alias = "CCR7+", pop = "+", parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8-/CD45RA+",
        dims = "CCR7", gating_method = "mindensity2", gating_args = "gate_range=c(1000,1400)",
        mc.cores = 14, parallel_type = "multicore", collapseDataForGating = T, groupBy = "Individual")
add_pop(gs_b1, alias = "CCR7+", pop = "+", parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8-/CD45RA-",
        dims = "CCR7", gating_method = "mindensity2", gating_args = "gate_range=c(1000,1450)",
        mc.cores = 14, parallel_type = "multicore", collapseDataForGating = T, groupBy = "Individual")
png(filename = file.path(gating_plots_dir_b1,
                         sprintf("%s_gate_%s.png", "CD45RA_on_CD4-CD8-",
                                 format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))),
    width = 1900, height = 1030)
ggcyto(gs_b1, mapping = aes(x = "CCR7", y = "CD45RA"),
       subset = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8-") +
  geom_hex(bins=120) +
  geom_gate(c("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8-/CD45RA+")) +
  facet_grid(Individual ~ TimePoint) +
  geom_stats(type = c("gate_name", "percent"), adjust = 0.5, fill = "white", alpha = 0.7) +
  ggcyto_par_set(limits = "instrument") +
  labs(title = "CD45RA on CD4- CD8-   collapsed by Individual", caption = batch) + facet_grid(Individual ~ TimePoint)
dev.off()
png(filename = file.path(gating_plots_dir_b1,
                         sprintf("%s_gate_%s.png", "CCR7_and_CD45RA_on_CD4-CD8-_CD45RA+",
                                 format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))),
    width = 1900, height = 1030)
ggcyto(gs_b1, mapping = aes(x = "CCR7", y = "CD45RA"),
       subset = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8-/CD45RA+") +
  geom_hex(bins=120) +
  geom_gate(c("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8-/CD45RA+",
              "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8-/CD45RA+/CCR7+")) +
  facet_grid(Individual ~ TimePoint) +
  geom_stats(type = c("gate_name", "percent"), adjust = 0.5, fill = "white", alpha = 0.7) +
  ggcyto_par_set(limits = "instrument") +
  labs(title = "CD45RA and CCR7 on CD4- CD8-/CD45RA+   collapsed by Individual", caption = batch) + facet_grid(Individual ~ TimePoint)
dev.off()
png(filename = file.path(gating_plots_dir_b1,
                         sprintf("%s_gate_%s.png", "CCR7_and_CD45RA_on_CD4-CD8-_CD45RA-",
                                 format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))),
    width = 1900, height = 1030)
ggcyto(gs_b1, mapping = aes(x = "CCR7", y = "CD45RA"),
       subset = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8-/CD45RA-") +
  geom_hex(bins=120) +
  geom_gate(c("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8-/CD45RA+",
              "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8-/CD45RA-/CCR7+")) +
  facet_grid(Individual ~ TimePoint) +
  geom_stats(type = c("gate_name", "percent"), adjust = 0.5, fill = "white", alpha = 0.7) +
  ggcyto_par_set(limits = "instrument") +
  labs(title = "CD45RA and CCR7 on CD4- CD8-/CD45RA-   collapsed by Individual", caption = batch) + facet_grid(Individual ~ TimePoint)
dev.off()

# CD4+ CD8+
add_pop(gs_b1, alias = "CD45RA", pop = "+/-", parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+ CD8+",
        dims = "CD45RA", gating_method = "refGate", gating_args = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8-/CD45RA+",
        mc.cores = 14, parallel_type = "multicore")
add_pop(gs_b1, alias = "CCR7", pop = "+/-", parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+ CD8+/CD45RA+",
        dims = "CCR7", gating_method = "refGate", gating_args = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8-/CD45RA+/CCR7+",
        mc.cores = 14, parallel_type = "multicore")
add_pop(gs_b1, alias = "CCR7", pop = "+/-", parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+ CD8+/CD45RA-",
        dims = "CCR7", gating_method = "refGate", gating_args = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8-/CD45RA-/CCR7+",
        mc.cores = 14, parallel_type = "multicore")
png(filename = file.path(gating_plots_dir_b1,
                         sprintf("%s_gate_%s.png", "CCR7_and_CD45RA_on_CD4+CD8+",
                                 format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))),
    width = 1900, height = 1030)
ggcyto(gs_b1, mapping = aes(x = "CCR7", y = "CD45RA"),
       subset = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+ CD8+") +
  geom_hex(bins=120) +
  geom_gate(c("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+ CD8+/CD45RA+",
              "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+ CD8+/CD45RA+/CCR7+",
              "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+ CD8+/CD45RA-/CCR7+")) +
  facet_grid(Individual ~ TimePoint) +
  geom_stats(type = c("gate_name", "percent"), adjust = 0.5, fill = "white", alpha = 0.7) +
  ggcyto_par_set(limits = "instrument") +
  labs(title = "CD45RA and CCR7 on CD4+ CD8+   collapsed by Individual", caption = batch) + facet_grid(Individual ~ TimePoint)
dev.off()

# Within each co-receptor group, create the memory population boolean gates
flowWorkspace::add(gs_b1, eval(substitute(flowWorkspace::booleanFilter(v),
                                          list(v = as.symbol("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+ CD8+/CD45RA+&/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+ CD8+/CD45RA+/CCR7+")))),
                   parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+ CD8+", name = "Naive")
flowWorkspace::add(gs_b1, eval(substitute(flowWorkspace::booleanFilter(v),
                                          list(v = as.symbol("!/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+ CD8+/CD45RA+&!/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+ CD8+/CD45RA-/CCR7+")))),
                   parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+ CD8+", name = "TEM")
flowWorkspace::add(gs_b1, eval(substitute(flowWorkspace::booleanFilter(v),
                                          list(v = as.symbol("!/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+ CD8+/CD45RA+&/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+ CD8+/CD45RA-/CCR7+")))),
                   parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+ CD8+", name = "TCM")
flowWorkspace::add(gs_b1, eval(substitute(flowWorkspace::booleanFilter(v),
                                          list(v = as.symbol("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+ CD8+/CD45RA+&!/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+ CD8+/CD45RA+/CCR7+")))),
                   parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+ CD8+", name = "TEMRA")

flowWorkspace::add(gs_b1, eval(substitute(flowWorkspace::booleanFilter(v),
                                          list(v = as.symbol("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8-/CD45RA+&/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8-/CD45RA+/CCR7+")))),
                   parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8-", name = "Naive")
flowWorkspace::add(gs_b1, eval(substitute(flowWorkspace::booleanFilter(v),
                                          list(v = as.symbol("!/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8-/CD45RA+&!/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8-/CD45RA-/CCR7+")))),
                   parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8-", name = "TEM")
flowWorkspace::add(gs_b1, eval(substitute(flowWorkspace::booleanFilter(v),
                                          list(v = as.symbol("!/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8-/CD45RA+&/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8-/CD45RA-/CCR7+")))),
                   parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8-", name = "TCM")
flowWorkspace::add(gs_b1, eval(substitute(flowWorkspace::booleanFilter(v),
                                          list(v = as.symbol("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8-/CD45RA+&!/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8-/CD45RA+/CCR7+")))),
                   parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8-", name = "TEMRA")

flowWorkspace::add(gs_b1, eval(substitute(flowWorkspace::booleanFilter(v),
                                          list(v = as.symbol("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8+/CD45RA+&/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8+/CD45RA+/CCR7+")))),
                   parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8+", name = "Naive")
flowWorkspace::add(gs_b1, eval(substitute(flowWorkspace::booleanFilter(v),
                                          list(v = as.symbol("!/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8+/CD45RA+&!/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8+/CD45RA-/CCR7+")))),
                   parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8+", name = "TEM")
flowWorkspace::add(gs_b1, eval(substitute(flowWorkspace::booleanFilter(v),
                                          list(v = as.symbol("!/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8+/CD45RA+&/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8+/CD45RA-/CCR7+")))),
                   parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8+", name = "TCM")
flowWorkspace::add(gs_b1, eval(substitute(flowWorkspace::booleanFilter(v),
                                          list(v = as.symbol("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8+/CD45RA+&!/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8+/CD45RA+/CCR7+")))),
                   parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8+", name = "TEMRA")

flowWorkspace::add(gs_b1, eval(substitute(flowWorkspace::booleanFilter(v),
                                          list(v = as.symbol("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+ CD8-/CD45RA+&/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+ CD8-/CCR7+")))),
                   parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+ CD8-", name = "Naive")
flowWorkspace::add(gs_b1, eval(substitute(flowWorkspace::booleanFilter(v),
                                          list(v = as.symbol("!/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+ CD8-/CD45RA+&!/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+ CD8-/CCR7+")))),
                   parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+ CD8-", name = "TEM")
flowWorkspace::add(gs_b1, eval(substitute(flowWorkspace::booleanFilter(v),
                                          list(v = as.symbol("!/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+ CD8-/CD45RA+&/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+ CD8-/CCR7+")))),
                   parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+ CD8-", name = "TCM")
flowWorkspace::add(gs_b1, eval(substitute(flowWorkspace::booleanFilter(v),
                                          list(v = as.symbol("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+ CD8-/CD45RA+&!/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+ CD8-/CCR7+")))),
                   parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+ CD8-", name = "TEMRA")

# Make boolean gates to identify the four quadrants
flowWorkspace::add(gs_b1, eval(substitute(flowWorkspace::booleanFilter(v),
                                          list(v = as.symbol("CD4+ CD8+/Naive|CD4+ CD8-/Naive|CD4- CD8-/Naive|CD4- CD8+/Naive")))),
                   parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes", name = "Naive")
flowWorkspace::add(gs_b1, eval(substitute(flowWorkspace::booleanFilter(v),
                                          list(v = as.symbol("CD4+ CD8+/TEM|CD4+ CD8-/TEM|CD4- CD8-/TEM|CD4- CD8+/TEM")))),
                   parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes", name = "TEM")
flowWorkspace::add(gs_b1, eval(substitute(flowWorkspace::booleanFilter(v),
                                          list(v = as.symbol("CD4+ CD8+/TCM|CD4+ CD8-/TCM|CD4- CD8-/TCM|CD4- CD8+/TCM")))),
                   parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes", name = "TCM")
flowWorkspace::add(gs_b1, eval(substitute(flowWorkspace::booleanFilter(v),
                                          list(v = as.symbol("CD4+ CD8+/TEMRA|CD4+ CD8-/TEMRA|CD4- CD8-/TEMRA|CD4- CD8+/TEMRA")))),
                   parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes", name = "TEMRA")

flowWorkspace::recompute(gs_b1, "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes")

getPopStats(gs_b1, subpopulations = c("Lymphocytes/Naive", "Lymphocytes/TEM", "Lymphocytes/TCM", "Lymphocytes/TEMRA"))

getNodes(gs_b1)

plot(gs_b1, fontsize = 25, bool = T)

# CD3_clean population

flowWorkspace::add(gs_b1, eval(substitute(flowWorkspace::booleanFilter(v),
                                          list(v = as.symbol("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/GMM_R660+|/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/GMM_V655+|/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/Unloaded_CD1b+&!/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/GMM_clean")))),
                   parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes", name = "GMM_dirty")
flowWorkspace::add(gs_b1, eval(substitute(flowWorkspace::booleanFilter(v),
                                          list(v = as.symbol("/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes&!/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/GMM_dirty")))),
                   parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes", name = "CD3_clean")
flowWorkspace::add(gs_b1, getGate(gs_b1, "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/GMM_clean"), parent = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD3_clean", name = "GMM_clean")
flowWorkspace::recompute(gs_b1, "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes")

###############################

# Make a plot of the automated gating tree
png(filename = file.path(gating_plots_dir_b1, "20180906_AutoGate_GatingTree.png"), width = 1000, height = 1000)
plot(gs_b1, fontsize = 25, bool = T)
dev.off()

# Save the GatingSets
save_gs(gs_b1, gs_out_dir_b1_autogate)
save_gs(fmo_b1, gs_out_dir_b1_fmo_autogate)

# GMM_clean count
gmm_clean_count <- dplyr::left_join(pData(gs_b1),
                                    getPopStats(gs_b1, subpopulation = "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD3_clean/GMM_clean"),
                                    by = c("rowname" = "name"))[, c("SAMPLE ID", "Individual", "TimePoint",
                                                                    "SampleDate", "TimePoint_v2", "Day",
                                                                    "Count", "ParentCount")]
gmm_clean_count
write.table(x = gmm_clean_count, file = file.path(gating_plots_dir_b1,
                                                  sprintf("%s_%s.tsv", "GMM_clean_stats",
                                                          format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))),
            sep = "\t", row.names = F)
