# Scripts for Durable expansion of TCR-δ meta-clonotypes after BCG revaccination in humans
# James, C. et al. 
# Script Authored: January 28, 2021 (kmayerbl)
# Reviewed: April 30, 2021 (kmayerb/)

### THIS SCRIPT ESTIMATE BETA-BINOMIAL MODELS. MAKES VOLCANO PLOTS, AND V01,V28 PLOTS. 
require(readr)
require(dplyr)
require(tidyr)
require(ggplot2)
require(corncob)
set.seed(210128)

# PRELIMINARY FUNCTIONS ####
#' <f> lambda-like function for adding feature column to a dataframe using - applied with purrr::map2(dataframes, names)
f  <- function(df, name){
  df = as.data.frame(df)
  df$variable = rownames(df)
  df['feature'] = name
  return(df)
}

########################################################################################################################
### Beta Binom Maximum Likelihood 
########################################################################################################################
#' do_corncob 
#' 
#' Define the beta-binomial we are attempting to fit
#' 
#' @param mydata data.frame
do_corncob <- function(mydata, frm = as.formula('cbind(W, M - W) ~ visit')){
  cb1 = bbdml(formula = frm,
              phi.formula = ~ 1,
              data = mydata)
  return(cb1)
}
# Run corncob, returning NA for failures, such as all samples in one class
possibly_do_corncob = purrr::possibly(do_corncob, otherwise = NA)
#' get bbdml coefficients into a table
#' 
#' @param cb is object result of corncob::bbdml
parse_corncob <- function(cb,i =1){
  y = summary(cb)$coefficients
  rdf = as.data.frame(y) 
  rdf$param = rownames(rdf)
  rdf = rdf %>% mutate(estimate = Estimate,  se = `Std. Error`, tvalue =  `t value`, pvalue = `Pr(>|t|)`, param) %>% 
    mutate(type = ifelse(grepl(param, pattern = "phi"), "phi", "mu")) %>% 
    mutate(type2 = ifelse(grepl(param, pattern = "Intercept"), "intercept", "covariate")) 
  rdf$feature = i
  return(rdf)
}

# Tabulation results, counts the number of clones in each 
# bulk samples that are conformant in meta-clonotypes. 
d = readr::read_tsv("2021_01_13_delta_BCG_metaclones.tsv.tab.tsv")
# All meta_clones, 101 identified
mc = readr::read_tsv('2021_01_13_delta_BCG_metaclones.tsv')
# Write out the metaclnootypes file
mc %>%  mutate(feature = paste0(v_d_gene, "+", cdr3_d_aa,"+", radius, "+", regex)) %>%
  select(nsubject,v_d_gene,j_d_gene,cdr3_d_aa,radius,regex,feature) %>% 
  readr::write_csv("James_et_al_2021/2021_01_13_delta_BCG_metaclones.tsv")

# Before moving on the regressions, 
# We audit the input files. 
# That is, open  all files ending TRD.tcrdist.tsv, 
# and get the calculated sum of templates and gen usages. 
input_dir = '/Volumes/T7/Seshadri/sampleExport_2020-07-22_23-33-32/tcrdist3ready/'
fs = list.files(input_dir, pattern = 'TRD.tcrdist.tsv')

# Function opens a file from input_dir, 
# and compute gene usage used in each bulk file
audit<-function(input_dir, f){
  df = readr::read_tsv(paste0(input_dir, f))
  # How many templates, per sample
  calc_sum_temp = df %>% arrange(desc(productive_frequency)) %>% summarise(sum(templates))
  gene_usage = df %>% arrange(desc(productive_frequency)) %>% 
    group_by(v_d_gene,j_d_gene) %>% summarise(n=n(), s_templates = sum(templates)) %>% arrange(desc(n))
  return(list('gene_usage'=gene_usage, 'calc_sum_temp'=calc_sum_temp))
}

# apply the audit function to the file list <fs> 
audit_df = purrr::map(fs, ~audit(input_dir= input_dir, f = .x))
# <audit_df> is a list for each file ($gene_usage), ($calc_sum_temp)
# Contains a M, the total counts per file.
# We do this to know the number of total productive reads per sample
audit_df2 = tibble(M = purrr::map_dbl(audit_df, ~.x$calc_sum_temp[[1]]), file = fs) %>%
  mutate(sample_key = stringr::str_replace(file, pattern = '.tsv.TRD.tcrdist.tsv', replacement = ''))
# M     file                             sample_key  
# 14343 109-V01_TCRa.tsv.TRD.tcrdist.tsv 109-V01_TCRa
#  5756 109-V08_TCRa.tsv.TRD.tcrdist.tsv 109-V08_TCRa
# 11595 109-V13_TCRa.tsv.TRD.tcrdist.tsv 109-V13_TCRa
# 15845 109-V28_TCRa.tsv.TRD.tcrdist.tsv 109-V28_TCRa
#  9580 122-V01_TCRa.tsv.TRD.tcrdist.tsv 122-V01_TCRa

# <mc2> adds a feature string by combining TRDV+CDR3+RADIUS 
mc2 = mc %>%  mutate(feature = paste0(v_d_gene, "+", cdr3_d_aa,"+", radius, "+", regex)) %>% 
  group_by(v_d_gene, j_d_gene, feature) %>% summarise(n=n())

# setup A dataframe for regressions
# M, total counts
# W0, counts exactly matching meta-clonotype centroids
# W, counts within centroid radius
# WR, counts withint centroid radius and regex
d2 = d %>% mutate(feature = paste0(v_d_gene, "+", cdr3_d_aa,"+", radius, "+", regex)) %>% 
  mutate(sample_key = stringr::str_replace(sample, pattern = '.tsv.TRD.tab1TRD.bulk_tabulation.tsv', replacement ='')) %>% 
  left_join(audit_df2) %>% 
  dplyr::select(sample_key, feature, visit, subject, M, W0 =  bulk_sum_counts_tcrdist0, WR = bulk_sum_counts_regex_adj, W= bulk_sum_counts, v_d_gene, j_d_gene) %>%
  mutate(f0 = W0/M, fR= WR/M, f=W/M)

# COMPARE VISIT0 (PRE) TO VISIT 28 (1-YEAR), 
# FOR EACH MEATCLONOTPE, 
# BREAK INTO SEPERATE DATAFRAMES, 
# WHERE THE ROWS ARE COUNTS PER SAMPLE
dfs = d2 %>% group_by(feature) %>%
  filter(visit %in% c('V01','V28')) %>% 
  group_split()

# GET FEATURE NAME, MATCHING THE DFS list of DATAFRAMES
dfs_names = purrr::map_chr(dfs, ~.x$feature[[1]])
# NOW NAME THE LIST
names(dfs) = dfs_names

##################
### FIT MODELS ###
##################
# FIT MODEL FOR EXACT, RADIUS, AND RADIUS+MOTIF Meta-clonotypes
bbrt_W = purrr::map(dfs, ~possibly_do_corncob(mydata = .x, frm = as.formula('cbind(W, M - W) ~ visit')))
bbrt_W0 = purrr::map(dfs, ~possibly_do_corncob(mydata = .x, frm = as.formula('cbind(W0, M - W0) ~ visit')))
bbrt_WR = purrr::map(dfs, ~possibly_do_corncob(mydata = .x, frm = as.formula('cbind(WR, M - WR) ~ visit')))

# PARSE (RADIUS) RESULTS
rW = do.call(rbind, purrr::map2(bbrt_W,names(dfs), ~parse_corncob(.x, .y))) %>% 
  filter(type == "mu") %>% 
  filter(type2 == "covariate") %>%
  arrange(pvalue) %>%
  mutate(method= 'RADIUS')

# PARSE (RADIUS + REGEX RESULTS)
rWR = do.call(rbind, purrr::map2(bbrt_WR,names(dfs), ~parse_corncob(.x, .y))) %>% 
  filter(type == "mu") %>% 
  filter(type2 == "covariate") %>%
  arrange(pvalue) %>% 
  mutate(method= 'RADIUS+MOTIF')

# PARSE (EXACT RESULTS)
rW0 = do.call(rbind, purrr::map2(bbrt_W0,names(dfs), ~parse_corncob(.x, .y))) %>% 
  filter(type == "mu") %>% 
  filter(type2 == "covariate") %>%
  arrange(pvalue) %>% 
  mutate(method= 'EXACT')

# COMBINE ALL THE RESULTS
# ALSO NEGATIVE, LOG10 TRANFORM PVALUE FOR VOLCANO PLOTS
R = rbind(rW,rWR,rW0) %>% 
  mutate(negative_log10_pvalue = -1 * log10(pvalue)) %>% 
  filter(abs(Estimate) < 5)
# SORT DESCENDING BY -LOG10(P)
R %>% arrange(desc(negative_log10_pvalue))

# WRITE REGRSSION RESULTS
R %>% arrange(desc(negative_log10_pvalue)) %>% 
  readr::write_tsv("James_et_al_2021/2021_01_13_delta_BCG_beta_binomial_results.tsv")
  #readr::write_csv("bioRxiv/SI/2021_01_13_Beta_Binomial.results.tsv")

# DEFINE VISUAL SCHEMA
pub_theme = theme_classic()+
  theme(axis.title.x = element_text(size = 12, face = 'bold')) +
  theme(axis.title.y = element_text(size = 14, face = 'bold')) + 
  theme(axis.text.x = element_text(size = 12, face= "bold", angle = 90)) + 
  theme(axis.text.y = element_text(size = 12, face= "bold")) + 
  theme(legend.position = "bottom") + 
  theme(strip.text = element_text(size = 12, face= "bold"),
        strip.background = element_rect(fill = "#00000010", size=1, color = "white")) + 
  theme(panel.grid.minor = element_blank())

#############
### PLOTS ###
#############
# VOLCANO PLOT
#pdf("bioRxiv/Beta_Binomial2.pdf", width = 7, height = 3.5)
ggdf = R %>% filter(method != "RADIUS") %>% 
  mutate(method = ifelse(method == "RADIUS+MOTIF", "Meta-Clonotypes", "Clonotypes")) %>% 
  mutate(param = ifelse(param == "mu.visitV28", "Post (V28)", param))
stopifnot(dim(ggdf)[1] < dim(R)[1])
ggplot(ggdf, aes(x = Estimate, y = negative_log10_pvalue)) + 
  geom_point(size = 1, shape= 1) + facet_wrap(method~param) + 
  theme_bw() +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  pub_theme + 
  ylab(expression(-log[10](P[value])))
#dev.off()

# NOW WE WANT TO LOOK AT RAW DATA FOR THE MOST SIGNIFICANT HIT(S)
R %>% arrange(desc(negative_log10_pvalue)) %>% head(6)
# Estimate Std. Error   t value   Pr(>|t|)       param   estimate        se    tvalue     pvalue type     type2                                                   feature       method negative_log10_pvalue
# 1  1.2385507  0.5829102  2.124771 0.04195756 mu.visitV28  1.2385507 0.5829102  2.124771 0.04195756   mu covariate     TRDV2*01+CACDTLLGDTRTDKLIF+18+(D[ST][LV]LGDT[RG]TDKL) RADIUS+MOTIF             1.3771897
# 2  1.8244380  1.0940744  1.667563 0.10580968 mu.visitV28  1.8244380 1.0940744  1.667563 0.10580968   mu covariate     TRDV2*01+CACDTLLGDTRTDKLIF+18+(D[ST][LV]LGDT[RG]TDKL)        EXACT             0.9754746
# 3 -0.9106175  0.5526086 -1.647853 0.10981872 mu.visitV28 -0.9106175 0.5526086 -1.647853 0.10981872   mu covariate        TRDV2*01+CACDYLTGGSFTDKLIF+18+(D[AY]LTGGS[FY]TDKL)       RADIUS             0.9593236
# 4 -0.9729900  0.5948649 -1.635649 0.11236366 mu.visitV28 -0.9729900 0.5948649 -1.635649 0.11236366   mu covariate   TRDV2*01+CACDPVLGDTRYTDKLIF+18+(D[PS]V[LM]GDTR[Y]?TDKL) RADIUS+MOTIF             0.9493741
# 5  1.1195023  0.6862182  1.631409 0.11325919 mu.visitV28  1.1195023 0.6862182  1.631409 0.11325919   mu covariate        TRDV2*01+CACDTLLGDTGTDKLIF+18+(D[TV]LLGDT[RG]TDKL) RADIUS+MOTIF             0.9459265
# 6  0.9698537  0.6043657  1.604746 0.11902680 mu.visitV28  0.9698537 0.6043657  1.604746 0.11902680   mu covariate TRDV2*01+CACDSLLGDTRTDKLIF+18+(D[ST][LMV]LG[DE][KT]RTDKL) RADIUS+MOTIF             0.9243552
feature_to_examine= 'TRDV2*01+CACDTLLGDTRTDKLIF+18+(D[ST][LV]LGDT[RG]TDKL)'
visits_to_compare = c('V01','V28')

# WHAT IS HALF THE LOWEST DETECTION
k = filter(d2, feature =='TRDV2*01+CACDTLLGDTRTDKLIF+18+(D[ST][LV]LGDT[RG]TDKL)', visit %in% c('V01','V28')) %>% pull(fR)
half_min_detect = min(k[k!=0]) / 2

# ON A PER SAMPLE BASIS WHAT IS 1/10th the detection limit (1/10)M
ggdf = filter(d2, feature =='TRDV2*01+CACDTLLGDTRTDKLIF+18+(D[ST][LV]LGDT[RG]TDKL)', visit %in% c('V01','V28')) %>% 
  mutate(fR_low = ifelse(fR == 0, .1/M, fR)) %>% 
  mutate(fR_nd = ifelse(fR == 0, "ND", "D")) 
# <both_down> is a dataframe that captures that both pre and post were zero V01+V28 == 0,
both_down = ggdf %>% 
  select(subject, feature, fR, visit) %>% ungroup() %>%
  tidyr::spread(key = visit, value=fR) %>% 
  mutate(V28 = ifelse(is.na(V28), 0,V28)) %>% 
  mutate(V01 = ifelse(is.na(V01), 0,V01)) %>% 
  mutate(both = ifelse(V01+V28 == 0, "ND", "D"))%>%
  mutate(either = ifelse(V01+V28 > 0, 1, 0))
# subject feature                                                     V01       V28 both 
# <dbl> <chr>                                                     <dbl>     <dbl> <chr>
# 1      52 TRDV2*01+CACDTLLGDTRTDKLIF+18+(D[ST][LV]LGDT[RG]TDKL) 0         0.000260  D    
# 2      71 TRDV2*01+CACDTLLGDTRTDKLIF+18+(D[ST][LV]LGDT[RG]TDKL) 0         0.000139  D    
# 3      86 TRDV2*01+CACDTLLGDTRTDKLIF+18+(D[ST][LV]LGDT[RG]TDKL) 0         0         ND   
# 4      93 TRDV2*01+CACDTLLGDTRTDKLIF+18+(D[ST][LV]LGDT[RG]TDKL) 0.0000214 0         D    
# 5      96 TRDV2*01+CACDTLLGDTRTDKLIF+18+(D[ST][LV]LGDT[RG]TDKL) 0         0.000191  D    
# 6      109 TRDV2*01+CACDTLLGDTRTDKLIF+18+(D[ST][LV]LGDT[RG]TDKL) 0         0         ND   
# 7      122 TRDV2*01+CACDTLLGDTRTDKLIF+18+(D[ST][LV]LGDT[RG]TDKL) 0         0.000374  D    
# 8      131 TRDV2*01+CACDTLLGDTRTDKLIF+18+(D[ST][LV]LGDT[RG]TDKL) 0         0         ND   
# 9      155 TRDV2*01+CACDTLLGDTRTDKLIF+18+(D[ST][LV]LGDT[RG]TDKL) 0.000632  0.000272  D    
# 10     159 TRDV2*01+CACDTLLGDTRTDKLIF+18+(D[ST][LV]LGDT[RG]TDKL) 0         0         ND   
# 11     164 TRDV2*01+CACDTLLGDTRTDKLIF+18+(D[ST][LV]LGDT[RG]TDKL) 0         0.0000774 D    
# 12     173 TRDV2*01+CACDTLLGDTRTDKLIF+18+(D[ST][LV]LGDT[RG]TDKL) 0         0.000416  D    
# 13     189 TRDV2*01+CACDTLLGDTRTDKLIF+18+(D[ST][LV]LGDT[RG]TDKL) 0         0         ND   
# 14     224 TRDV2*01+CACDTLLGDTRTDKLIF+18+(D[ST][LV]LGDT[RG]TDKL) 0.00148   0.00203   D    
# 15     228 TRDV2*01+CACDTLLGDTRTDKLIF+18+(D[ST][LV]LGDT[RG]TDKL) 0.000614  0.000864  D    
# 16     238 TRDV2*01+CACDTLLGDTRTDKLIF+18+(D[ST][LV]LGDT[RG]TDKL) 0         0.000434  D    
# 17     241 TRDV2*01+CACDTLLGDTRTDKLIF+18+(D[ST][LV]LGDT[RG]TDKL) 0         0         ND  

# HOW MANY ARE DETECTED IN V01
sum(both_down$V01 > 0)
# 4
sum(both_down$V28 > 0)
# 10

# HOW MANY EXPANDED BETWEEN V01 to V28
sum(both_down$V01 < both_down$V28)
ggdf2 = ggdf %>% left_join(both_down)
stopifnot(dim(ggdf)[1] == dim(ggdf2)[1])

# SET A PUBLISHING THEME
pub_theme = theme_classic()+
  theme(axis.title.x = element_text(size = 10, face = 'bold')) +
  theme(axis.title.y = element_text(size = 10, face = 'bold')) + 
  theme(axis.text.x = element_text(size = 10, face= "bold")) + 
  theme(axis.text.y = element_text(size = 10, face= "bold")) + 
  theme(legend.position = "bottom") + 
  theme(strip.text = element_text(size = 10, face= "bold"),
        strip.background = element_rect(fill = "#0000FF10", size=1, color = "white")) + 
  theme(panel.grid.minor = element_blank())

# MAKE A PLOT OF CHANGE BETWEEN V01,V28
#pdf(paste("James_et_al_2021/fR_raw_changes_",feature_to_examine,".pdf"), width = 4, height = 3.5)
ggplot(ggdf2, 
       aes(x = visit,
           y = fR_low,
           group = subject,
           shape = fR_nd ,
           fill = fR_nd)) + 
  geom_point(size =2) + 
  geom_line(alpha= .9, aes(linetype = both)) +
  scale_shape_manual("",values = c(21,22)) +
  scale_fill_manual("", values = c("black","white")) +
  scale_x_discrete(expand = c(.05, .05)) + 
  scale_color_discrete('') + 
  scale_y_log10(expand = c(.01,.01)) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 20)) + 
  theme(axis.text.y = element_text(angle = 0, face = "bold", size = 20)) + ylab("") +
  annotation_logticks(side = 'lr') + 
  pub_theme + 
  xlab("") + ylab("") +
  theme(legend.position = "none")
#dev.off()

## END ##
