library(here)
library(openCyto)
library(ggcyto)
library(lubridate)

pre_inh_metadata_file <- here::here("data/151127_PreINH_Manifest.xlsx")
bcg_vacc_metadata_file <- here::here("data/150518_BCGVacc_Manifest.xlsx")
pre_inh_metadata <- readxl::read_excel(pre_inh_metadata_file, sheet = 1, col_names = T, n_max = 20)
pre_inh_metadata$Timepoint <- "Pre"
bcg_vacc_metadata <- readxl::read_excel(bcg_vacc_metadata_file, sheet = 1, col_names = T)
all_ext_metadata <- rbind(pre_inh_metadata[, c("Sample ID", "Sample Date", "Timepoint")],
                          bcg_vacc_metadata[, c("Sample ID", "Sample Date", "Timepoint")])
all_ext_metadata$`Sample ID` <- as.character(all_ext_metadata$`Sample ID`)
colnames(all_ext_metadata)[which(colnames(all_ext_metadata) == "Sample Date")] <- "SampleDate"

#################
# Batch 1
#################

workspace_path_b1 <- here::here("data/20180712 DURT BCG B1/Concatenated Files/20180712 DURT BCG B1 XML .xml")
fcs_path_b1 <- here::here("data/20180712 DURT BCG B1/Concatenated Files")
gs_out_dir_b1 <- here::here("out/GatingSets/20180908_GatingSet_Batch1")

ws_b1 <- openWorkspace(workspace_path_b1)
tail(plyr::join(getSampleGroups(ws_b1), getSamples(ws_b1), by = "sampleID"))
unlist(getKeywords(ws_b1, "86 V08")[grep("\\$|P[0-9]", names(getKeywords(ws_b1, "86 V08")), invert = T, value = T)])
keywords2import=c("SAMPLE ID", "$DATE", "EXPORT TIME", "TUBE NAME", "WELL ID", "EXPERIMENT NAME", "PLATE NAME")

gs_b1 <- parseWorkspace(ws_b1, name = "Samples",
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

plot(gs_b1, bool = T, fontsize = 25)

# tSNE requires all marker names be set.
markernames(gs_b1)
pData(parameters(getData(gs_b1[[1]])))[,c(1,2)]
gs_markernames <- c("Time", "FSC-A", "FSC-H", "SSC-A", "SSC-H", "MR1", "Unloaded_CD1b", "TRAV1-2", 
                    "CCR7", "CD14/CD19", "L/D", "CD8b", "CD3", "CD45RA", 
                    "MA_R660", "CD4", "aGalCer", "MA_G610", "GD")
names(gs_markernames) <- pData(parameters(getData(gs_b1[[1]])))[,1]
markernames(gs_b1) <- gs_markernames
markernames(gs_b1)
pData(parameters(getData(gs_b1[[1]])))[,c(1,2)]

save_gs(gs_b1, gs_out_dir_b1)

#################
# Batch 2
#################

workspace_path_b2 <- here::here("data/20180727 DURT BCG B2/Concatenated Files/20180727 DURT BCG B2 XML .xml")
fcs_path_b2 <- here::here("data/20180727 DURT BCG B2/Concatenated Files")
gs_out_dir_b2 <- here::here("out/GatingSets/20180908_GatingSet_Batch2")

ws_b2 <- openWorkspace(workspace_path_b2)
tail(plyr::join(getSampleGroups(ws_b2), getSamples(ws_b2), by = "sampleID"))
keywords2import=c("SAMPLE ID", "$DATE", "EXPORT TIME", "TUBE NAME", "WELL ID", "EXPERIMENT NAME", "PLATE NAME")

gs_b2 <- parseWorkspace(ws_b2, name = "Samples",
                        keywords = keywords2import,
                        path = fcs_path_b2)
pData(gs_b2)$filename <- sapply(rownames(pData(gs_b2)), function(x) {
  gsub("/.+/", "", description(getData(gs_b2[[x]]))$FILENAME)
}, USE.NAMES = F)
pData(gs_b2)$rowname <- rownames(pData(gs_b2))
pData(gs_b2)$Individual <- sapply(strsplit(pData(gs_b2)$`SAMPLE ID`, " "), function(x) { x[[1]] })
pData(gs_b2)$TimePoint <- sapply(strsplit(pData(gs_b2)$`SAMPLE ID`, " "), function(x) { x[[2]] })
gs_b2_new <- dplyr::left_join(pData(gs_b2), all_ext_metadata, by = c("Individual" = "Sample ID", "TimePoint" = "Timepoint"))
rownames(gs_b2_new) <- gs_b2_new$rowname
pData(gs_b2) <- gs_b2_new
pData(gs_b2)$TimePoint_v2 <- factor(plyr::revalue(pData(gs_b2)$TimePoint, c("Pre" = "Pre-INH",
                                                                            "V08" = "Pre-BCG",
                                                                            "V13" = "3 weeks",
                                                                            "V28" = "1 year")),
                                    levels = c("Pre-INH", "Pre-BCG", "3 weeks", "1 year"))
# Add a Day column which sets the Pre-BCG SampleDate as Day 0 and calculates all the other days based on that. Then I can possibly use a linear x-axis.
pre_bcg_indices <- which(pData(gs_b2)$TimePoint_v2 == "Pre-BCG")
pre_bcg_dates <- pData(gs_b2)[pre_bcg_indices, "SampleDate"]
names(pre_bcg_dates) <- pData(gs_b2)[pre_bcg_indices, "Individual"]
pData(gs_b2)$Day <- purrr::map2_dbl(pData(gs_b2)$Individual, pData(gs_b2)$SampleDate,
                                    function(ind, current_date) {
                                      ymd(current_date) - ymd(pre_bcg_dates[[ind]])
                                    })
pData(gs_b2)

plot(gs_b2, bool = T, fontsize = 25)

# tSNE requires all marker names be set.
markernames(gs_b2)
pData(parameters(getData(gs_b2[[1]])))[,c(1,2)]
gs_markernames <- c("Time", "FSC-A", "FSC-H", "SSC-A", "SSC-H", "MR1", "Unloaded_CD1b", "TRAV1-2", 
                    "CCR7", "CD14/CD19", "L/D", "CD8b", "CD3", "CD45RA", 
                    "MA_R660", "CD4", "aGalCer", "MA_G610", "GD")
names(gs_markernames) <- pData(parameters(getData(gs_b2[[1]])))[,1]
markernames(gs_b2) <- gs_markernames
markernames(gs_b2)
pData(parameters(getData(gs_b2[[1]])))[,c(1,2)]

save_gs(gs_b2, gs_out_dir_b2)

#################
# Batch 3
#################

workspace_path_b3 <- here::here("data/20180731 DURT BCG B3/Concatenated Files/20180731 DURT BCG B3 XML.xml")
fcs_path_b3 <- here::here("data/20180731 DURT BCG B3/Concatenated Files")
gs_out_dir_b3 <- here::here("out/GatingSets/20180908_GatingSet_Batch3")

ws_b3 <- openWorkspace(workspace_path_b3)
tail(plyr::join(getSampleGroups(ws_b3), getSamples(ws_b3), by = "sampleID"))
keywords2import=c("SAMPLE ID", "$DATE", "EXPORT TIME", "TUBE NAME", "WELL ID", "EXPERIMENT NAME", "PLATE NAME")

gs_b3 <- parseWorkspace(ws_b3, name = "Samples",
                        keywords = keywords2import,
                        path = fcs_path_b3)
pData(gs_b3)$filename <- sapply(rownames(pData(gs_b3)), function(x) {
  gsub("/.+/", "", description(getData(gs_b3[[x]]))$FILENAME)
}, USE.NAMES = F)
pData(gs_b3)$rowname <- rownames(pData(gs_b3))
pData(gs_b3)$Individual <- sapply(strsplit(pData(gs_b3)$`SAMPLE ID`, " "), function(x) { x[[1]] })
pData(gs_b3)$TimePoint <- sapply(strsplit(pData(gs_b3)$`SAMPLE ID`, " "), function(x) { x[[2]] })
gs_b3_new <- dplyr::left_join(pData(gs_b3), all_ext_metadata, by = c("Individual" = "Sample ID", "TimePoint" = "Timepoint"))
rownames(gs_b3_new) <- gs_b3_new$rowname
pData(gs_b3) <- gs_b3_new
pData(gs_b3)$TimePoint_v2 <- factor(plyr::revalue(pData(gs_b3)$TimePoint, c("Pre" = "Pre-INH",
                                                                            "V08" = "Pre-BCG",
                                                                            "V13" = "3 weeks",
                                                                            "V28" = "1 year")),
                                    levels = c("Pre-INH", "Pre-BCG", "3 weeks", "1 year"))
# Add a Day column which sets the Pre-BCG SampleDate as Day 0 and calculates all the other days based on that. Then I can possibly use a linear x-axis.
pre_bcg_indices <- which(pData(gs_b3)$TimePoint_v2 == "Pre-BCG")
pre_bcg_dates <- pData(gs_b3)[pre_bcg_indices, "SampleDate"]
names(pre_bcg_dates) <- pData(gs_b3)[pre_bcg_indices, "Individual"]
pData(gs_b3)$Day <- purrr::map2_dbl(pData(gs_b3)$Individual, pData(gs_b3)$SampleDate,
                                    function(ind, current_date) {
                                      ymd(current_date) - ymd(pre_bcg_dates[[ind]])
                                    })
pData(gs_b3)

plot(gs_b3, bool = T, fontsize = 25)

# tSNE requires all marker names be set.
markernames(gs_b3)
pData(parameters(getData(gs_b3[[1]])))[,c(1,2)]
gs_markernames <- c("Time", "FSC-A", "FSC-H", "SSC-A", "SSC-H", "MR1", "Unloaded_CD1b", "TRAV1-2", 
                    "CCR7", "CD14/CD19", "L/D", "CD8b", "CD3", "CD45RA", 
                    "MA_R660", "CD4", "aGalCer", "MA_G610", "GD")
names(gs_markernames) <- pData(parameters(getData(gs_b3[[1]])))[,1]
markernames(gs_b3) <- gs_markernames
markernames(gs_b3)
pData(parameters(getData(gs_b3[[1]])))[,c(1,2)]

save_gs(gs_b3, gs_out_dir_b3)

#################
# Batch 4
#################

workspace_path_b4 <- here::here("data/20180803 DURT BCG B4/Concatenated Files/20180803 DURT BCG B4 XML.xml")
fcs_path_b4 <- here::here("data/20180803 DURT BCG B4/Concatenated Files")
gs_out_dir_b4 <- here::here("out/GatingSets/20180908_GatingSet_Batch4")

ws_b4 <- openWorkspace(workspace_path_b4)
tail(plyr::join(getSampleGroups(ws_b4), getSamples(ws_b4), by = "sampleID"))
keywords2import=c("SAMPLE ID", "$DATE", "EXPORT TIME", "TUBE NAME", "WELL ID", "EXPERIMENT NAME", "PLATE NAME")

gs_b4 <- parseWorkspace(ws_b4, name = "Samples",
                        keywords = keywords2import,
                        path = fcs_path_b4)
pData(gs_b4)$filename <- sapply(rownames(pData(gs_b4)), function(x) {
  gsub("/.+/", "", description(getData(gs_b4[[x]]))$FILENAME)
}, USE.NAMES = F)
pData(gs_b4)$rowname <- rownames(pData(gs_b4))

# Fix Batch_4 mislabelling:
# "Sample 224 V28 is mislabeled, the file name is correct on all files, but the sample ID reads as 228 V28."
pData(gs_b4)[, c("name", "SAMPLE ID", "filename")]
pData(gs_b4)[pData(gs_b4)[, c("SAMPLE ID")] != pData(gs_b4)[, c("filename")], c("name", "SAMPLE ID", "filename")]
pData(gs_b4)[which(pData(gs_b4)$filename == "224 V28"), "SAMPLE ID"] <- "224 V28"
pData(gs_b4)[pData(gs_b4)[, c("SAMPLE ID")] != pData(gs_b4)[, c("filename")], c("name", "SAMPLE ID", "filename")]

pData(gs_b4)$Individual <- sapply(strsplit(pData(gs_b4)$`SAMPLE ID`, " "), function(x) { x[[1]] })
pData(gs_b4)$TimePoint <- sapply(strsplit(pData(gs_b4)$`SAMPLE ID`, " "), function(x) { x[[2]] })
gs_b4_new <- dplyr::left_join(pData(gs_b4), all_ext_metadata, by = c("Individual" = "Sample ID", "TimePoint" = "Timepoint"))
rownames(gs_b4_new) <- gs_b4_new$rowname
pData(gs_b4) <- gs_b4_new
pData(gs_b4)$TimePoint_v2 <- factor(plyr::revalue(pData(gs_b4)$TimePoint, c("Pre" = "Pre-INH",
                                                                            "V08" = "Pre-BCG",
                                                                            "V13" = "3 weeks",
                                                                            "V28" = "1 year")),
                                    levels = c("Pre-INH", "Pre-BCG", "3 weeks", "1 year"))
# Add a Day column which sets the Pre-BCG SampleDate as Day 0 and calculates all the other days based on that. Then I can possibly use a linear x-axis.
pre_bcg_indices <- which(pData(gs_b4)$TimePoint_v2 == "Pre-BCG")
pre_bcg_dates <- pData(gs_b4)[pre_bcg_indices, "SampleDate"]
names(pre_bcg_dates) <- pData(gs_b4)[pre_bcg_indices, "Individual"]
pData(gs_b4)$Day <- purrr::map2_dbl(pData(gs_b4)$Individual, pData(gs_b4)$SampleDate,
                                    function(ind, current_date) {
                                      ymd(current_date) - ymd(pre_bcg_dates[[ind]])
                                    })
pData(gs_b4)

plot(gs_b4, bool = T, fontsize = 25)

# tSNE requires all marker names be set.
markernames(gs_b4)
pData(parameters(getData(gs_b4[[1]])))[,c(1,2)]
gs_markernames <- c("Time", "FSC-A", "FSC-H", "SSC-A", "SSC-H", "MR1", "Unloaded_CD1b", "TRAV1-2", 
                    "CCR7", "CD14/CD19", "L/D", "CD8b", "CD3", "CD45RA", 
                    "MA_R660", "CD4", "aGalCer", "MA_G610", "GD")
names(gs_markernames) <- pData(parameters(getData(gs_b4[[1]])))[,1]
markernames(gs_b4) <- gs_markernames
markernames(gs_b4)
pData(parameters(getData(gs_b4[[1]])))[,c(1,2)]

save_gs(gs_b4, gs_out_dir_b4)
