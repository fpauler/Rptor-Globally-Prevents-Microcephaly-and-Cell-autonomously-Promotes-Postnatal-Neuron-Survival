# ====================================
# analyses a list of DEGs with iDEA
# arguments from command line:
# ====================================

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
	error ( "expected are 3 command line arguments ( in that exact order) : [annotation RDS] [DEG list RDS] [prefix char]" )
}

annoRDS <- args[1]
degRDS <- args[2]
prefix = args[3]

message( annoRDS, " ", degRDS, " ", prefix )

library (iDEA)

# a data.frame with gene names as rows and term names as columns
# cells show a 1 if a gene is associated with a trem or 0 if not
go_anno <- readRDS(file = annoRDS )

# this is a list of lists
# first level is genotype, second level is t-type
de_res <- readRDS(file = degRDS )

idea_GO_results <- list()

# the DEG dataframe has several levels - analyse here 
for ( level1 in names(de_res) ) {

  idea_GO_results[[ level1 ]] <- list()

  for ( level2 in  names(de_res[[ level1 ]]) ) {

    idea_GO_results[[ level1 ]][[ level2 ]] <- list()    
    
    #for (dir in c("neg", "pos")) {
    for (dir in c("all")) {

       message( level1, " ", level2, " ", dir )

       ## Assume you have obtained the DE results from i.e. zingeR, edgeR or MAST with the data frame res_DE (column: pvalue and LogFC)
       if (dir == "pos") {
          deg_df <- de_res[[ level1 ]][[ level2 ]][ which(de_res[[ level1 ]][[ level2 ]]$avg_log2FC > 0), ]
       } else if (dir == "neg") {
          deg_df <- de_res[[ level1 ]][[ level2 ]][ which(de_res[[ level1 ]][[ level2 ]]$avg_log2FC < 0), ]
       } else {
         deg_df <- de_res[[ level1 ]][[ level2 ]]
       }

       pvalue <- deg_df$p_val #### the pvalue column
       zscore <- qnorm(pvalue/2.0, lower.tail=FALSE) #### convert the pvalue to z-score
       beta <- deg_df$avg_log2FC ## effect size
       se_beta <- abs(beta/zscore) ## to approximate the standard error of beta
       beta_var = se_beta^2  ### square 
       summary = data.frame(beta = beta,beta_var = beta_var)
       ## add the gene names as the rownames of summary
       rownames(summary) = rownames(deg_df) ### or the gene id column in the res_DE results
    
       idea <- CreateiDEAObject(summary, go_anno, max_var_beta = 100, min_precent_annot = 0.0025, num_core=10)
    
       
       idea <- iDEA.fit(idea,
                      fit_noGS=FALSE,
                      init_beta=NULL, 
                      init_tau=c(-2,0.5),
                      min_degene=5,
                      em_iter=15,
                      mcmc_iter=1000, 
                      fit.tol=1e-5,
                      modelVariant = F,
                      verbose=TRUE)

       # deal with potential empty list entries that cause errors for BMA calculations
       empty_list_ids <- which( sapply( idea@de, function (x) length(x)) == 0)
       filled_list_ids <- which( sapply( idea@de, function (x) length(x)) > 0)
       idea@de[ empty_list_ids ] = NULL
       idea@annot_id <- idea@annot_id[ filled_list_ids ]
       idea@annotation <- idea@annotation[ filled_list_ids ]
        
        idea <- iDEA.louis(idea)
        idea <- iDEA.BMA(idea)
        
       idea_GO_results[[ level1 ]][[ level2 ]][[ dir ]] <- idea
    }
  }
}

# save the analysis
saveRDS( object = idea_GO_results, file = paste(prefix,"idea_results.RDS", sep="." ) )

# report the environment
sessionInfo()
