# =======================================================
# step 1 of revision re-analysis
# prepare updated annotations for GO enrichment analysis
# revision version
# =======================================================

library (GO.db)
library (data.table)
library (org.Mm.eg.db)
library (openxlsx)

set.seed(2401)

# =======================
# set up I/O folders
# =======================
main_base <- "./bulk_RNASeq/final/"

# ===============================================================
# read the manually curated GO term inclusion / exclusion terms
# ===============================================================
GO_list <- read.csv(file = paste(main_base, "/other_data/GO_Term_categories_v1.csv", sep=""))

# ============================================
# GO term extraction / filtering / reporting
# ============================================

filter_terms <-GO_list$Include[which( nchar(GO_list$Include) > 0) ]

# prepare the GO lists
includeIDX <- sapply( filter_terms, function (y) {
  which(grepl(pattern =  y, x = Term(GOTERM), ignore.case = F, fixed = T ))
})

includeIDX <- unique( includeIDX )
                 
excludeIDX <- sapply( GO_list$Exclude[which( nchar(GO_list$Exclude) > 0) ], function (y) {
  which(grepl(pattern =  y, x = Term(GOTERM), ignore.case = F, fixed = T ))
})

includeIDX <- sapply( includeIDX, function (x) setdiff( x, unique( unlist(excludeIDX) )))

GO_term_list <- sapply( includeIDX, function (x)  Term(GOTERM)[ x ] )
names(GO_term_list) <- filter_terms

# combine all TOR related terms
GO_term_list[["TOR"]] <- c( GO_term_list[["TOR"]],
                            GO_term_list[["phosphatidylinositol 3-kinase/protein kinase B signal"]] )
GO_term_list[["phosphatidylinositol 3-kinase/protein kinase B signal"]] <- NA

filter_terms[ which(filter_terms == "TOR")] <- "TOR/PI3K"
filter_terms[ which(filter_terms == "epidermal growth factor")] <- "Egfr"
filter_terms[ which(filter_terms == "reactive oxygen")] <- "ROS"

names(GO_term_list) <- filter_terms

# number of GO terms used for analysis
length(unique(unlist(GO_term_list)))
# 867

# write out to xlsx file
# write out GO classification determined above
wb <- createWorkbook()

for (name in names(GO_term_list)) {
  sheetName <- name
  if (!any(is.na(GO_term_list[[ name ]]))) {
    tmp_df <- as.data.frame(GO_term_list[[ name ]])
    tmp_df$GO <- rownames(tmp_df)
    colnames( tmp_df ) <- c("Description", "ID")
    
    addWorksheet(wb, sheetName)
    
    writeData( wb, sheetName, tmp_df )
    addFilter( wb, sheetName, row = 1, cols = 1:ncol( tmp_df ) )
    setColWidths( wb, sheetName, cols = 1:ncol( tmp_df ), widths="auto" )
  }
}

saveWorkbook(wb, file = paste(main_base, "Supplement/GO_terms_used.xlsx", sep=""), overwrite = T) 

# =============================================================================
# prepare a dataframe linking genes to GO term groups used for later analyses
# =============================================================================

sapply( includeIDX, function (y) {
  sapply(y, function (x){
    go_id <- GOID( GOTERM[ x ] )
    allegs <- tryCatch( get(go_id, org.Mm.egGO2ALLEGS), 
                        error = function (error) {
                          # message(error, "\n")
                          return (NA)
                        }  )
    
    # get symbols
    if (!any(is.na(allegs))) {
      genes <- unique(unlist(mget(allegs,org.Mm.egSYMBOL)))
      return(genes)
    }
  })

}) -> genes_per_group

genes_per_group <- sapply( 1:length(genes_per_group), function (x) {
  if ("list" %in% class(genes_per_group[[x]]) ) {
    unique(unlist(genes_per_group[[x]]))
  } else {
    return ( as.character(genes_per_group[[x]]) )
  }  
})

names( genes_per_group ) <- GO_list$Include[which( nchar(GO_list$Include) > 0) ]

df <- data.frame(
  name = rep(names(genes_per_group), lengths(genes_per_group)),
  value = unlist(genes_per_group),
  row.names = NULL,
  stringsAsFactors = FALSE
)

df$name[ which(df$name %in% c("TOR", "phosphatidylinositol 3-kinase/protein kinase B signal")) ] <- "TOR/PI3K"

# save the data frame
saveRDS(object = df, file = paste(main_base, "/RDS_files/GO_Term_df_collapsed.RDS", sep="") )


sessionInfo()
# R version 4.3.2 (2023-10-31)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 22.04.5 LTS
# 
# Matrix products: default
# BLAS:   /opt/R/4.3.2/lib/R/lib/libRblas.so 
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0
# 
# locale:
#   [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=de_AT.UTF-8       
# [4] LC_COLLATE=en_GB.UTF-8     LC_MONETARY=de_AT.UTF-8    LC_MESSAGES=en_GB.UTF-8   
# [7] LC_PAPER=de_AT.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
# [10] LC_TELEPHONE=C             LC_MEASUREMENT=de_AT.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Europe/Vienna
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] openxlsx_4.2.5.2     org.Mm.eg.db_3.18.0  data.table_1.15.0    GO.db_3.18.0        
# [5] AnnotationDbi_1.64.1 IRanges_2.36.0       S4Vectors_0.40.2     Biobase_2.62.0      
# [9] BiocGenerics_0.48.1 
# 
# loaded via a namespace (and not attached):
#   [1] bit_4.0.5               compiler_4.3.2          crayon_1.5.2            Rcpp_1.0.12            
# [5] zip_2.3.1               blob_1.2.4              bitops_1.0-7            Biostrings_2.70.3      
# [9] png_0.1-8               fastmap_1.1.1           R6_2.5.1                XVector_0.42.0         
# [13] igraph_2.1.4            GenomeInfoDb_1.38.8     GenomeInfoDbData_1.2.11 DBI_1.2.2              
# [17] rlang_1.1.3             KEGGREST_1.42.0         cachem_1.0.8            stringi_1.8.3          
# [21] bit64_4.0.5             RSQLite_2.3.6           memoise_2.0.1           cli_3.6.2              
# [25] magrittr_2.0.3          zlibbioc_1.48.2         rstudioapi_0.16.0       packrat_0.9.2          
# [29] lifecycle_1.0.4         vctrs_0.6.5             RCurl_1.98-1.14         httr_1.4.7             
# [33] tools_4.3.2             pkgconfig_2.0.3  
