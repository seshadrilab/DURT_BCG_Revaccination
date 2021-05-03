# library(here)
library(openCyto)
library(ggcyto)
# library(lubridate)

# pre_inh_metadata_file <- here::here("data/151127_PreINH_Manifest.xlsx")
# bcg_vacc_metadata_file <- here::here("data/150518_BCGVacc_Manifest.xlsx")
# pre_inh_metadata <- readxl::read_excel(pre_inh_metadata_file, sheet = 1, col_names = T, n_max = 20)
# pre_inh_metadata$Timepoint <- "Pre"
# bcg_vacc_metadata <- readxl::read_excel(bcg_vacc_metadata_file, sheet = 1, col_names = T)
# all_ext_metadata <- rbind(pre_inh_metadata[, c("Sample ID", "Sample Date", "Timepoint")],
#                           bcg_vacc_metadata[, c("Sample ID", "Sample Date", "Timepoint")])
# all_ext_metadata$`Sample ID` <- as.character(all_ext_metadata$`Sample ID`)
# colnames(all_ext_metadata)[which(colnames(all_ext_metadata) == "Sample Date")] <- "SampleDate"

#################
# Batch 1
#################

workspace_path_b1 <- here::here("data/20180712 DURT BCG B1/Concatenated Files/20180712 DURT BCG B1 XML .xml")
fcs_path_b1_controls <- here::here("data/20180712 DURT BCG B1")
gs_out_dir_b1_ctrl <- here::here("out/GatingSets/20190314_GatingSet_Batch1_Controls")

ws_b1 <- openWorkspace(workspace_path_b1)
plyr::join(getSampleGroups(ws_b1), getSamples(ws_b1), by = "sampleID")
getKeywords(ws_b1, "27637.fcs")[grep("\\$|P[0-9]", names(getKeywords(ws_b1, "27637.fcs")), invert = T, value = T)]
keywords2import=c("SAMPLE ID", "$DATE", "EXPORT TIME", "TUBE NAME", "WELL ID", "EXPERIMENT NAME", "PLATE NAME")

gs_b1_controls <- parseWorkspace(ws_b1, name = "Controls",
                        keywords = keywords2import,
                        path = fcs_path_b1_controls)
pData(gs_b1_controls)$filename <- sapply(rownames(pData(gs_b1_controls)), function(x) {
  gsub("/.+/", "", description(getData(gs_b1_controls[[x]]))$FILENAME)
}, USE.NAMES = F)
pData(gs_b1_controls)$rowname <- rownames(pData(gs_b1_controls))

pData(gs_b1_controls)

plot(gs_b1_controls, bool = T, fontsize = 25)

# tSNE requires all marker names be set.
markernames(gs_b1_controls)
pData(parameters(getData(gs_b1_controls[[1]])))[,c(1,2)]
# gs_markernames <- c("Time", "FSC-A", "FSC-H", "SSC-A", "SSC-H", "MR1", "Unloaded_CD1b", "TRAV1-2", 
#                     "CCR7", "CD14/CD19", "L/D", "CD8b", "CD3", "CD45RA", 
#                     "MA_R660", "CD4", "aGalCer", "MA_G610", "GD")
# names(gs_markernames) <- pData(parameters(getData(gs_b1_controls[[1]])))[,1]
# markernames(gs_b1_controls) <- gs_markernames
# markernames(gs_b1_controls)
# pData(parameters(getData(gs_b1_controls[[1]])))[,c(1,2)]

save_gs(gs_b1_controls, gs_out_dir_b1_ctrl)
