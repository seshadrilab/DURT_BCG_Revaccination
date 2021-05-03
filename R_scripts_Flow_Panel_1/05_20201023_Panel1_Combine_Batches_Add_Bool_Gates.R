library(openCyto) # 1.24.0
library(CytoML) # 1.12.0
library(flowCore) # required for description()
library(flowWorkspace) # required for gh_pop_get_data()
library(here)
library(tidyverse)

# 20201023
# re-generate the GatingSet in order to complete analyses for the manuscript.

# Also, add the following boolean gates:
# GMM-CD1b+ TRAV1-2-
# GMM-CD1b+ CD8+
# GMM-CD1b+ Naive
# GMM+TRAV1-2+CD4+

gs_dir_b1 <- here::here("out/GatingSets/20180906_GatingSet_Batch1_AutoGate")
gs_dir_b2 <- here::here("out/GatingSets/20180906_GatingSet_Batch2_AutoGate")
gs_dir_b3 <- here::here("out/GatingSets/20180907_GatingSet_Batch3_AutoGate")
gs_dir_b4 <- here::here("out/GatingSets/20180907_GatingSet_Batch4_AutoGate")

gs_combined_dir_tmp <- here::here("out/GatingSets/20180907_GatingSets_Combined_AutoGate_tmp") # An intermediate combined GatingSet
gs_combined_dir <- here::here("out/GatingSets/20180907_GatingSets_Combined_AutoGate_mod20201023")

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
# gs <- load_gs(gs_combined_dir_tmp)

# Add gates of interest

plot(gs, fontsize = 25, bool = T)

lymphocytes_path <- "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes"
cd3_clean_path <- "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD3_clean"
gmm_clean_path_parent_lymph <- "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/GMM_clean"
trav1_2_path <- "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/TRAV1-2+"
cd8_path <- "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4- CD8+"
cd4_path <- "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD4+ CD8-"
naive_path <- "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/Naive"

# GMM+TRAV1-2+CD4+
# originally defined under lymphocytes_path, so do this again to be consistent with original GatingSet
gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("", gmm_clean_path_parent_lymph,
                                                         "&", trav1_2_path,
                                                         "&", cd4_path))))),
           parent = lymphocytes_path, name = "GEM_GMM")
gem_gmm_path_parent_lymph <- "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/GEM_GMM"

# GMM-CD1b+ TRAV1-2-
gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("", gmm_clean_path_parent_lymph,
                                                         "&!", trav1_2_path))))),
           parent = cd3_clean_path, name = "GMM_and_not_TRAV1_2")
GMM_and_not_TRAV1_2_path <- "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD3_clean/GMM_and_not_TRAV1_2"

# GMM-CD1b+ CD8+
gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("", gmm_clean_path_parent_lymph,
                                                         "&", cd8_path))))),
           parent = cd3_clean_path, name = "GMM_CD8")
GMM_CD8_path <- "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD3_clean/GMM_CD8"

# GMM-CD1b+ Naive
gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("", gmm_clean_path_parent_lymph,
                                                         "&", naive_path))))),
           parent = cd3_clean_path, name = "GMM_Naive")
GMM_Naive_path <- "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD3_clean/GMM_Naive"

# To avoid confusion, also copy the populations of interest to be children of "denominator" populations, where necessary
# # Copy GMM_clean under CD3_clean --> nvmd, this already exists
# gmm_clean_gate_parent_lymph <- getGate(gs, gmm_clean_path_parent_lymph)
# flowWorkspace::add(gs, gmm_clean_gate_parent_lymph, parent=cd3_clean_path) 
gmm_clean_path_parent_CD3_clean <- "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD3_clean/GMM_clean"
# Copy GEM_GMM under CD3_clean
gem_gmm_gate_parent_lymph <- getGate(gs, gem_gmm_path_parent_lymph)
flowWorkspace::add(gs, gem_gmm_gate_parent_lymph, parent=cd3_clean_path) 
gem_gmm_path_parent_CD3_clean <- "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD3_clean/GEM_GMM"
# Copy GEM_GMM under CD3_clean/GMM_clean
gem_gmm_gate_parent_CD3_clean <- getGate(gs, gem_gmm_path_parent_CD3_clean)
flowWorkspace::add(gs, gem_gmm_gate_parent_CD3_clean, parent=gmm_clean_path_parent_CD3_clean) 
gem_gmm_path_parent_GMM_clean <- "/Live-CD3+/CD14-CD19-/S/CCR7 Keeper/CD38 Keeper/Lymphocytes/CD3_clean/GMM_clean/GEM_GMM"

# Make a plot of the gating tree
plot(gs, bool = T, fontsize = 25)
# and save it
png(filename = here::here("out/GatingPlots/20201023_GatingTree_AutoGated_Allbatches_AddnlBoolGates.png"), width = 1500, height = 1000)
plot(gs, fontsize = 25, bool = T)
dev.off()

flowWorkspace::recompute(gs, lymphocytes_path)
gc()

save_gs(gs, gs_combined_dir)

##########################################

# While we're at it, extract cell counts for populations of interest, and save to disk
pop_paths_of_interest <- c(cd3_clean_path,
                           gmm_clean_path_parent_CD3_clean,
                           gem_gmm_path_parent_CD3_clean,
                           GMM_and_not_TRAV1_2_path,
                           GMM_CD8_path,
                           GMM_Naive_path)
new_names_pop_paths_of_interest <- c("CD3_clean",
                                     "GMM",
                                     "GEM_GMM",
                                     "GMM_Not_TRAV1_2",
                                     "GMM_CD8",
                                     "GMM_Naive")
pop_dat <- gs_pop_get_count_fast(gs, subpopulations = pop_paths_of_interest) %>% 
  dplyr::select(name, Population, Count) %>% 
  pivot_wider(names_from = Population, values_from = Count) %>% 
  rename_at(vars(all_of(pop_paths_of_interest)),
            ~ new_names_pop_paths_of_interest) %>% 
  left_join(pData(gs) %>% 
              dplyr::select("SAMPLE ID", "Individual", "TimePoint_v2", "Batch", "SampleDate", "Day", "rowname"),
            by = c("name" = "rowname")) %>% 
  rename("SAMPLE_ID" = `SAMPLE ID`) %>% 
  dplyr::select(SAMPLE_ID, "Individual", "TimePoint_v2", "Batch", "SampleDate", "Day",
                all_of(new_names_pop_paths_of_interest))
head(pop_dat)

if(save_output) {
  write.csv(pop_dat, here::here("out/GatingPlots/20201023_BCG_SATVI_TBRU_GMM_Cell_Counts.csv"), row.names = F)
  # pop_dat <- read.csv(here::here("out/20201023_BCG_SATVI_TBRU_GMM_Cell_Counts.csv"), check.names = F, stringsAsFactors = F)
}