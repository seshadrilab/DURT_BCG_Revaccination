library(here)
library(openCyto)
library(ggcyto)

gs_dir_b1 <- here::here("out/GatingSets/20180908_GatingSet_Batch1")
gs_dir_b2 <- here::here("out/GatingSets/20180908_GatingSet_Batch2")
gs_dir_b3 <- here::here("out/GatingSets/20180908_GatingSet_Batch3")
gs_dir_b4 <- here::here("out/GatingSets/20180908_GatingSet_Batch4")

gs_combined_dir_tmp <- here::here("out/GatingSets/20180908_GatingSet_AllBatches_tmp")
gs_combined_dir <- here::here("out/GatingSets/20180908_GatingSet_AllBatches")

gs_b1 <- load_gs(gs_dir_b1)
gs_b2 <- load_gs(gs_dir_b2)
gs_b3 <- load_gs(gs_dir_b3)
gs_b4 <- load_gs(gs_dir_b4)

pData(gs_b1)$Batch <- "Batch 1"
pData(gs_b2)$Batch <- "Batch 2"
pData(gs_b3)$Batch <- "Batch 3"
pData(gs_b4)$Batch <- "Batch 4"

gs_list <- GatingSetList(list(gs_b1, gs_b2, gs_b3, gs_b4))
gs <- flowWorkspace::rbind2(gs_list)

save_gs(gs, gs_combined_dir_tmp)

# note that the ECD+ gate is for MA_G610
plotGate(gs_b1, "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/ECD+")
# APC+ gate is for MA_R660
plotGate(gs_b1, "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/APC+")
# Unloaded CD1b
plotGate(gs_b1, "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/Endo+")

# note that the ECD+ gate is for MA_G610
plotGate(gs_b2, "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/ECD+")
# APC+ gate is for MA_R660
plotGate(gs_b2, "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/APC+")
# Unloaded CD1b
plotGate(gs_b2, "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/Endo+")

# note that the ECD+ gate is for MA_G610
plotGate(gs_b3, "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/ECD+")
# APC+ gate is for MA_R660
plotGate(gs_b3, "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/APC+")
# Unloaded CD1b
plotGate(gs_b3, "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/Endo+")

# note that the ECD+ gate is for MA_G610
plotGate(gs_b4, "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/ECD+")
# APC+ gate is for MA_R660
plotGate(gs_b4, "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/APC+")
# Unloaded CD1b
plotGate(gs_b4, "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/Endo+")

# So, 3 tetramer positive populations to run tSNE on:
# Mycolic Acid: tSNE run population defined by ECD+ APC+, Endo-. Small population, so we want to be extra sure what we see is real. Also filter out aGC+ cells, but this is accomplished via CD3_clean
# MR1: include TRAV1-2- events for now
# aGalCer: Perform 2 tSNE runs: 1 including GammaDelta and one not including GammaDelta. This is because, although we are interested in the GammaDelta population, some of it may be artifact.

############

# Add boolean gates

# MA_clean defined by MA (ECD)+, MA (APC)+, Endo -, aGalCer-, MR1-
flowWorkspace::add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                                       list(v = as.symbol("/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/ECD+&/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/APC+&!/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/Endo+&!/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/aGalCer+&!/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/MR1+")))),
                                                          parent = "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes", name = "MA_clean")

# aGalCer_clean defined by aGalCer+, MA (ECD)-, MA (APC)-, Endo -, MR1-
flowWorkspace::add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                                       list(v = as.symbol("/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/aGalCer+&!/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/ECD+&!/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/APC+&!/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/Endo+&!/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/MR1+")))),
                   parent = "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes", name = "aGalCer_clean")

# MR1_clean defined by MR1+, aGalCer-, MA (ECD)-, MA (APC)-, Endo -
flowWorkspace::add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                                       list(v = as.symbol("/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/MR1+&!/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/aGalCer+&!/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/ECD+&!/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/APC+&!/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/Endo+")))),
                   parent = "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes", name = "MR1_clean")

# Define "dirty" gates. These events will be excluded from CD3_clean, so they should contain everything we need to throw away.
# Since there are 2 MA tetramers, MA_dirty defined by MA (ECD)+  OR  MA (APC)+  AND NOT  MA_clean
# Remember to exclude Endo+ from CD3_clean gate separately (there could be some Endo+ that is MA ECD/APC-, which wouldn't get removed by MA_dirty)
flowWorkspace::add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                                       list(v = as.symbol("/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/ECD+|/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/APC+&!/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/MA_clean")))),
                   parent = "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes", name = "MA_dirty")

flowWorkspace::add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                                       list(v = as.symbol("/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/aGalCer+&!/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/aGalCer_clean")))),
                   parent = "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes", name = "aGalCer_dirty")

flowWorkspace::add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                                       list(v = as.symbol("/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/MR1+&!/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/MR1_clean")))),
                   parent = "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes", name = "MR1_dirty")

# CD3_clean
flowWorkspace::add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                                       list(v = as.symbol("/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes&!/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/MA_dirty&!/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/aGalCer_dirty&!/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/MR1_dirty&!/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/Endo+")))),
                   parent = "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes", name = "CD3_clean")

plot(gs, bool = T, fontsize = 25)

# Define the MA "GEM" cells:
# remembering to include CD4 :)
flowWorkspace::add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                                       list(v = as.symbol("/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/MA_clean&/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/TRAV1-2+&/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD4+")))),
                   parent = "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes", name = "GEM_MA")

# Technically at this point all the *_clean gates should be included in CD3_clean, but let's still copy over all the nodes to go under CD3_clean in case we want to calculate arbitrary population proportions out of CD3_clean

plot(gs, bool = T, fontsize = 25)

# Copy original gates under CD3_clean
nodes2CopyCD3 <- c("/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD4+", 
                   "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD8+", 
                   "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/DN", 
                   "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/DP",
                   "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/GammaDelta", 
                   "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/MR1+",
                   "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/TRAV1-2+", 
                   "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/MA_clean", 
                   "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/aGalCer_clean", 
                   "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/MR1_clean",
                   "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/aGalCer+",
                   "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/GEM_MA")
nodes2CopyCD4 <- c("/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD4+/Q1: CCR7Ð, CD45RA+", 
                   "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD4+/Q2: CCR7+, CD45RA+", 
                   "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD4+/Q3: CCR7+, CD45RAÐ", 
                   "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD4+/Q4: CCR7Ð, CD45RAÐ")
nodes2CopyCD8 <- c("/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD8+/Q1: CCR7Ð, CD45RA+", 
                   "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD8+/Q2: CCR7+, CD45RA+", 
                   "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD8+/Q3: CCR7+, CD45RAÐ", 
                   "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD8+/Q4: CCR7Ð, CD45RAÐ")
nodes2CopyDN <- c("/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/DN/Q1: CCR7Ð, CD45RA+", 
                  "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/DN/Q2: CCR7+, CD45RA+", 
                  "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/DN/Q3: CCR7+, CD45RAÐ", 
                  "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/DN/Q4: CCR7Ð, CD45RAÐ")
nodes2CopyDP <- c("/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/DP/Q1: CCR7Ð, CD45RA+", 
                  "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/DP/Q2: CCR7+, CD45RA+", 
                  "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/DP/Q3: CCR7+, CD45RAÐ", 
                  "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/DP/Q4: CCR7Ð, CD45RAÐ")
for(node in nodes2CopyCD3) {
  gate <- getGate(gs, node)
  flowWorkspace::add(gs, gate, parent="/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD3_clean") 
}
for(node in nodes2CopyCD4) {
  gate <- getGate(gs, node)
  flowWorkspace::add(gs, gate, parent="/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD3_clean/CD4+") 
}
for(node in nodes2CopyCD8) {
  gate <- getGate(gs, node)
  flowWorkspace::add(gs, gate, parent="/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD3_clean/CD8+") 
}
for(node in nodes2CopyDN) {
  gate <- getGate(gs, node)
  flowWorkspace::add(gs, gate, parent="/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD3_clean/DN") 
}
for(node in nodes2CopyDP) {
  gate <- getGate(gs, node)
  flowWorkspace::add(gs, gate, parent="/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD3_clean/DP") 
}

# MAIT
flowWorkspace::add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                                       list(v = as.symbol("/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD3_clean/MR1_clean&/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD3_clean/TRAV1-2+")))),
                   parent = "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD3_clean", name = "MAIT")

# aGalCer_clean_NoGD.
flowWorkspace::add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                                       list(v = as.symbol("/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD3_clean/aGalCer_clean&!/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD3_clean/GammaDelta")))),
                   parent = "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD3_clean", name = "aGalCer_clean_NoGD")

gc()

# Memory populations
flowWorkspace::add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                                       list(v = as.symbol("/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD3_clean/CD4+/Q2: CCR7+, CD45RA+|/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD3_clean/CD8+/Q2: CCR7+, CD45RA+|/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD3_clean/DN/Q2: CCR7+, CD45RA+|/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD3_clean/DP/Q2: CCR7+, CD45RA+")))),
                   parent = "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD3_clean", name = "Naive")

flowWorkspace::add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                                       list(v = as.symbol("/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD3_clean/CD4+/Q4: CCR7Ð, CD45RAÐ|/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD3_clean/CD8+/Q4: CCR7Ð, CD45RAÐ|/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD3_clean/DN/Q4: CCR7Ð, CD45RAÐ|/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD3_clean/DP/Q4: CCR7Ð, CD45RAÐ")))),
                   parent = "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD3_clean", name = "TEM")

flowWorkspace::add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                                       list(v = as.symbol("/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD3_clean/CD4+/Q1: CCR7Ð, CD45RA+|/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD3_clean/CD8+/Q1: CCR7Ð, CD45RA+|/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD3_clean/DN/Q1: CCR7Ð, CD45RA+|/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD3_clean/DP/Q1: CCR7Ð, CD45RA+")))),
                   parent = "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD3_clean", name = "TEMRA")

flowWorkspace::add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                                       list(v = as.symbol("/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD3_clean/CD4+/Q3: CCR7+, CD45RAÐ|/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD3_clean/CD8+/Q3: CCR7+, CD45RAÐ|/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD3_clean/DN/Q3: CCR7+, CD45RAÐ|/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD3_clean/DP/Q3: CCR7+, CD45RAÐ")))),
                   parent = "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD3_clean", name = "TCM")

flowWorkspace::recompute(gs, "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes")
gc()

cd3stats <- getPopStats(gs, subpopulation = "CD3_clean")
write.table(x = cd3stats,
            file = here::here("out/GatingPlots/CD3_stats.tsv"), row.names = F, col.names = T, sep = "\t")

tetstats <- getPopStats(gs, subpopulation = c("CD3_clean/MA_clean", "CD3_clean/aGalCer_clean",
                                              "CD3_clean/MR1_clean", "CD3_clean/MAIT", "CD3_clean/aGalCer_clean_NoGD"))[order(Population, name),]
write.table(x = tetstats,
            file = here::here("out/GatingPlots/Tet_stats.tsv"), row.names = F, col.names = T, sep = "\t")


plot(gs, bool = T, fontsize = 25)

# get pop stats here

save_gs(gs, gs_combined_dir)

ggcyto(subset(gs, Batch == "Batch 1"), aes(x = "MA_G610", y = "MA_R660"), subset = "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/") +
  geom_hex(bins=60) +
  ggcyto_par_set(limits = "instrument") +
  facet_grid(Individual ~ TimePoint) +
  geom_overlay("/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD3_clean/MA_clean", size = 0.2, alpha = 0.1, col = "red")
ggcyto(subset(gs, Batch == "Batch 2"), aes(x = "MA_G610", y = "MA_R660"), subset = "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/") +
  geom_hex(bins=60) +
  ggcyto_par_set(limits = "instrument") +
  facet_grid(Individual ~ TimePoint) +
  geom_overlay("/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD3_clean/MA_clean", size = 0.2, alpha = 0.1, col = "red")
ggcyto(subset(gs, Batch == "Batch 3"), aes(x = "MA_G610", y = "MA_R660"), subset = "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/") +
  geom_hex(bins=60) +
  ggcyto_par_set(limits = "instrument") +
  facet_grid(Individual ~ TimePoint) +
  geom_overlay("/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD3_clean/MA_clean", size = 0.2, alpha = 0.1, col = "red")
ggcyto(subset(gs, Batch == "Batch 4"), aes(x = "MA_G610", y = "MA_R660"), subset = "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/") +
  geom_hex(bins=60) +
  ggcyto_par_set(limits = "instrument") +
  facet_grid(Individual ~ TimePoint) +
  geom_overlay("/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD3_clean/MA_clean", size = 0.2, alpha = 0.1, col = "red")


ggcyto(subset(gs, Batch %in% c("Batch 1", "Batch 2")), aes(x = Time, y = MA_R660), subset = "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/") +
  geom_hex(bins=60) +
  ggcyto_par_set(limits = "instrument") +
  facet_grid(Individual ~ TimePoint)
ggcyto(subset(gs, Batch %in% c("Batch 3", "Batch 4")), aes(x = Time, y = MA_R660), subset = "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/") +
  geom_hex(bins=60) +
  ggcyto_par_set(limits = "instrument") +
  facet_grid(Individual ~ TimePoint)

ggcyto(subset(gs, Batch %in% c("Batch 1", "Batch 2")), aes(x = Time, y = MA_G610), subset = "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/") +
  geom_hex(bins=60) +
  ggcyto_par_set(limits = "instrument") +
  facet_grid(Individual ~ TimePoint)
ggcyto(subset(gs, Batch %in% c("Batch 3", "Batch 4")), aes(x = Time, y = MA_G610), subset = "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/") +
  geom_hex(bins=60) +
  ggcyto_par_set(limits = "instrument") +
  facet_grid(Individual ~ TimePoint)

ggcyto(subset(gs, Batch %in% c("Batch 1", "Batch 2")), aes(x = Time, y = Unloaded_CD1b), subset = "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/") +
  geom_hex(bins=60) +
  ggcyto_par_set(limits = "instrument") +
  facet_grid(Individual ~ TimePoint)
ggcyto(subset(gs, Batch %in% c("Batch 3", "Batch 4")), aes(x = Time, y = Unloaded_CD1b), subset = "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/") +
  geom_hex(bins=60) +
  ggcyto_par_set(limits = "instrument") +
  facet_grid(Individual ~ TimePoint)

##############################

gs <- load_gs(gs_combined_dir)

# This is a really heavy plot... RStudio crashes
# ggcyto(gs, aes(x = "MA_R660", y = "MA_G610"), subset = "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD3_clean") +
#   geom_hex(bins=60) +
#   ggcyto_par_set(limits = "instrument") +
#   facet_wrap(. ~ Individual) +
#   geom_overlay("/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD3_clean/MA_clean", size = 0.2, alpha = 0.1, col = "red")

# ggcyto(gs[["T003009 W4_5284656"]], aes(x = "MA_R660", y = "MA_G610"), subset = "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD3_clean") +
#   geom_hex(bins=240) +
#   ggcyto_par_set(limits = "instrument") +
#   facet_wrap(. ~ Individual) +
#   geom_overlay("/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD3_clean", size = 0.2, alpha = 0.1, col = "red")
#
# Above plot shows that when geom_hex() is used, it performs a sampling of the underlying events

recompute(gs, "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/")

attributes(getGate(gs, "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/APC+")[[1]])$boundaries

ggcyto(gs, aes(x = "FSC-A", y = "MA_R660"), subset = "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/") +
  geom_hex(bins=240) +
  ggcyto_par_set(limits = "instrument") +
  facet_grid(Individual ~ TimePoint) +
  geom_gate("/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/APC+")

ggcyto(gs, aes(x = "FSC-A", y = "MA_G610"), subset = "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/") +
  geom_hex(bins=240) +
  ggcyto_par_set(limits = "instrument") +
  facet_grid(Individual ~ TimePoint) +
  geom_gate("/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/ECD+")

ggcyto(gs, aes(x = "<U395-A>", y = "Unloaded_CD1b"), subset = "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/") +
  geom_hex(bins=240) +
  ggcyto_par_set(limits = "instrument") +
  facet_grid(Individual ~ TimePoint) +
  geom_gate("/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/Endo+")

ggcyto(gs, aes(x = "MA_G610", y = "MA_R660"), subset = "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/") +
  geom_hex(bins=240) +
  ggcyto_par_set(limits = "instrument") +
  facet_grid(Individual ~ TimePoint) +
  geom_overlay("/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD3_clean/MA_clean", size = 0.2, alpha = 0.1, col = "red")

# plotGate looks really different from geom_hex. It seems to remove some of the data when there are more bins


ggcyto(gs, aes(x = "MA_G610", y = "MA_R660"), subset = "/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/") +
  geom_hex(bins=60) +
  ggcyto_par_set(limits = "instrument") +
  facet_grid(Individual ~ TimePoint) +
  geom_overlay("/Live CD3+/CD14-19-/Singlet/CCR7 Keeper/MR1 Keeper/Lymphocytes/CD3_clean/MA_clean", size = 0.2, alpha = 0.1, col = "red")
