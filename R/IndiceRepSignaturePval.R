#' @title IndiceRepSignaturePval
#'
#' @description this function returns a named list : each element of this list is a numeric vector that gives values for each droplet of a specific metric except the last element, a character vector that #' gives the status of each cell. The fist element of list called ratio_interest_sign_vs_random provides for each droplet the ratio between the percentage of detected genes in the interest signature and #' the mean percentage of detected genes in random signatures of the same size than the interest signature.
#' The second element of the list called detection_rate_signature gives the percentage of detected genes in the interest signature.      
#' The third element of the list called mean_random_detection_rate_per_cell gives the mean percentage of detected genes in all the random signatures.
#' The fourth element called p-value gives the numeric value (p-value) that allow to assess if the detection rate of the cell is significantly higher or lower than expected randomly. To distinguish the #'  two cases (lower or higher) when the percentage of detection rate in the interest signature is higher than the median percentage of detection rate, the sign of the value is positive. In the 
#'  alternative case (percentage of detected genes in the interest signature is lower than the median percentage of detected genes in the random signatures)  
#' The fifth element called abs_pvalue is the same thing as the fourth element (p-value) but all the values are positive (original values).
#' The sixth element called statut_pvalue is a character vector that gives the status of each cell according to the value of the fourth element : POS_SIGNI (detection rate higher than expected randomly), #' NEG_SIGNI (detection rate lower than expected randomly) 
#' or NOT_SIGNI (detection rate not significantly different than expected randomly).   
#'
#'
#' @return (named list). This function returns a named list : each element of this list is a numeric vector that gives values for each droplet of a specific metric except the last element, a character #'   vector that gives the status of each cell.
#' The fist element of list called ratio_interest_sign_vs_random provides for each droplet the ratio between the percentage of detected genes in the interest signature and the mean percentage of detected #' genes in random signatures of the same size than the interest signature.
#' The second element of the list called detection_rate_signature gives the percentage of detected genes in the interest signature.      
#' The third element of the list called mean_random_detection_rate_per_cell gives the mean percentage of detected genes in all the random signatures.
#' The fourth element called p-value gives the numeric value (p-value) that allow to assess if the detection rate of the cell is significantly higher or lower than expected randomly. To distinguish the #'  two cases (lower or higher) when the percentage of detection rate in the interest signature is higher than the median percentage of detection rate, the sign of the value is positive. In the #'  
#' alternative case (percentage of detected genes in the interest signature is lower than the median percentage of detected genes in the random signatures)  
#' The fifth element called abs_pvalue is the same thing as the fourth element (p-value) but all the values are positive (original values).
#' The sixth element called statut_pvalue is a character vector that gives the status of each cell according to the value of the fourth element : POS_SIGNI (detection rate higher than expected randomly), #' NEG_SIGNI (detection rate lower than expected randomly) or NOT_SIGNI (detection rate not significantly different than expected randomly). 
#'
#' @param exprMat object of class matrix that corresponds to the single-cell expression matrix.
#'
#' @param biological_signature  a character vector that corresponds to a list of genes of interest. 
#'
#' @param min_prop_cells_detected Numerical value that corresponds to the minimal proportion of genes that have to be detected to take into account a gene in the matrix. 
#'
#' @param number_random_signatures integer value. The number of random signatures that will be used to generate the null distribution.
#'
#' @param thresholdpval numeric value which corresponds to the p-value threshold applied to determine whether or not a cell is significantly enriched for a given biological signature.
#' 
#' @param random_seed integer value random seed to ensure the reproducibility of the results
#'
#' @importFrom purrr pmap
#'
#' @importFrom dplyr case_when
#'
#' @importFrom Seurat Read10X CreateSeuratObject PercentageFeatureSet SCTransform RunPCA RunUMAP FindNeighbors FindClusters DimPlot FeaturePlot
#'
#' @export
#'
#' @examples 
#'
#' require(Seurat)
#'
#' 
#'
#' platelet <- read.table(GenicSignature_example("platelet_signature"),header=FALSE)[[1]]
#'
#' results_platelets <- IndiceRepSignaturePval(exprMat = processPBMC@assays$SCT@counts,
#'biological_signature = platelet,
#'min_prop_cells_detected = 0.02,
#'number_random_signatures = 1000,
#'thresholdpval = 0.01,
#'random_seed = 1)
#'
#' processPBMC@meta.data$status_platelets <- results_platelets$statut_fdr
#'
#' processPBMC$platelets_ratio_query_sign_vs_random <- results_platelets$ratio_interest_sign_vs_random
#'
#'
#' Seurat::DimPlot(processPBMC,group.by= "status_platelets")
#'
#'
#' Seurat::FeaturePlot(processPBMC,"platelets_ratio_query_sign_vs_random")
#'


# This function can be used several times with different value of the random seed in order to compute the mean of the pvalues.
IndiceRepSignaturePval <- function(exprMat,biological_signature,min_prop_cells_detected,number_random_signatures = 100,thresholdpval,random_seed){
#opening bracket of the function IndiceRepSignaturePval
  
#Setting of a random seed to ensure the reproducibility of the analysis.  
print("Setting of a random seed to ensure the reproducibility of the analysis.")
set.seed(random_seed) 

#Selection of the genes that will be used to test signatures. The resulting matrix that contains only these genes will be called used_exprMat 
print("Selection of the genes that will be used to test signatures. The resulting matrix that contains only these genes will be called used_exprMat")


#Selection of the genes that will be used to generate the random signatures. The resulting matrix only containing these genes will be called used_exprMat. 

used_exprMat <- exprMat[which(unlist(apply(exprMat,1,function(x){length(which(x>0))})) >= min_prop_cells_detected*ncol(exprMat)),]

 
########################### determination of the detection rate for the random signatures #########################

#Generation of 100 random signatures (we only keep the positions of genes in the matrix)

#We only keep from  the signature provided as input of the function (the set of genes provided as value of the argument biological_signature) genes that are in the matrix (that is those that are a detection rate higher than the minimal detection rate) in the expression matrix.
echantillonnage <- replicate(number_random_signatures,sample(rownames(used_exprMat),length(intersect(rownames(used_exprMat),biological_signature))),simplify = FALSE) 

echantillonnage <- lapply(echantillonnage,function(x){
#opening brackets of the function allowing to retrieve the positions of genes in the matrix.
  which(rownames(used_exprMat) %in% x)
  
  
  #closing brackets of the function that allows to retrieve positione of genes in the matrix.
})


#for each of the cell, we can test the detection rate associated to each of the random signatures generated above.
print("for each of the droplets, we apply for each of the signature.")

#for each of the droplets, we apply for each of the signature.
test <- apply(used_exprMat,2,function(z){
  

#for the z ieme droplet, we call the vector associated to the column corresponding to the droplet droplets.  
droplets <- z


#for each droplet the detection rate associated to each random signature generated above is computed.  
percentages <- lapply(echantillonnage,function(y){#computation of the percentage of the number of genes detected for each of the random signatures.
  
   (length(which(droplets[y]>0))/length(droplets[y]))*100
   
   
   
 } ## computation of the detection rate (percentage) for each droplet.
 ) 
  
  
  return(percentages)
  
  
})

print("for each droplet, we store in a vector the percentages of detected genes for each of the random signature.")

#for each droplet, we store in a vector the percentages of detected genes for each of the random signature.
test <- lapply(test,unlist)

#names of the elements of the storage list will correspond to names of the matrix columns.

print("names of the elements of the storage list will correspond to names of the matrix columns.")
names(test) <- colnames(used_exprMat)
#pour chacune des cellules, on détermine la moyenne du taux de détection des signatures aléatoires (en pourcentage)
print("for each cell, we determine the mean of the detection rate for each random signatures (percentage)")


mean_random_detection_rate_per_cell <- lapply(test,mean) 

mean_random_detection_rate_per_cell <- unlist(mean_random_detection_rate_per_cell)


############### Determination for the true biological signature of the detection rate associated to each droplet. #########################

print("Determination for the true biological signature of the detection rate associated to each droplet")
#computation for the real biological signature of the detection rate (percentage) for each droplet.
detection_rate_signature <- apply(used_exprMat,2,function(z){
  
    
  return(
  (length(which(z[which(rownames(used_exprMat) %in% intersect(rownames(used_exprMat),biological_signature))]>0))/length(z[which(rownames(used_exprMat) %in% intersect(rownames(used_exprMat),biological_signature))]))*100
  )
})

 
#computation of the ration detection_rate_signature on mean_random_detection_rate_per_cell
print("Calcul du ratio detection_rate_signature sur mean_random_detection_rate_per_cell") 

## ratio between the percentage of detected genes in the interest signature and the mean percentage of detected genes in all the random signatures generated. 
ratio_interest_sign_vs_random <- (detection_rate_signature/mean_random_detection_rate_per_cell)


### computation of the p-value : probability to obtain a given detection rate under the hypothesis that the signature has been generated  randomly (no biological meaning)
print("computation of the p-value : probability to obtain a given detection rate under the hypothesis that the signature has been generated  randomly (no biological meaning)")

pvaltoret <- purrr::pmap(.l = list(random = test , signature = detection_rate_signature),.f = function(random,signature){
  
return(pvaldistri(random,signature))
  
})


pvaltoret <- unlist(pvaltoret)

fdrtoret <- stats::p.adjust(abs(pvaltoret),method = "BH", n = length(pvaltoret))

####################################### Determination for each signature for each cell #######################

#Determination of the enrichment status for each signature for each cell.
print("Determination of the enrichment status for each signature for each cell.")
statutpval <- dplyr::case_when(pvaltoret < 0 & abs(pvaltoret) < thresholdpval ~ "NEG_SIGNI",
                                                                pvaltoret < 0 & abs(pvaltoret) >= thresholdpval ~ "NOT_SIGNI",
                                                                pvaltoret > 0 & abs(pvaltoret) < thresholdpval ~ "POS_SIGNI",
                                                                pvaltoret > 0 & abs(pvaltoret) >= thresholdpval ~ "NOT_SIGNI")


statut_fdr <- dplyr::case_when(pvaltoret < 0 & fdrtoret < thresholdpval ~ "NEG_SIGNI",
                                                                pvaltoret < 0 & fdrtoret >= thresholdpval ~ "NOT_SIGNI",
                                                                pvaltoret > 0 & fdrtoret < thresholdpval ~ "POS_SIGNI",
                                                                pvaltoret > 0 & fdrtoret >= thresholdpval ~ "NOT_SIGNI")



#We name the elements of the vector that give the status associated to each droplet in function of the p-value
names(statutpval) <- names(pvaltoret)


#We give a name to the elements of the vector that give the status associated to each droplet according to the fdr value.
names(statut_fdr) <- names(fdrtoret)
######################################### We return a list that give for each droplet values associated to parameters that allow to assess the enrichment of a signature in a given cell. ##  
    return(list("ratio_interest_sign_vs_random" = ratio_interest_sign_vs_random,
               "detection_rate_signature" = detection_rate_signature,
               "mean_random_detection_rate_per_cell" = mean_random_detection_rate_per_cell,
               "p-value" = pvaltoret,
               "abs_pvalue" = abs(pvaltoret), 
               "statut_pvalue" = statutpval,
               "FDR" = fdrtoret,
               "statut_fdr" = statut_fdr
              ))
  
## closing bracket of the function IndiceRepSignaturePval
}



# load(processPBMC_example(data = "processPBMC"))
# pbmc_data <- Seurat::Read10X(data.dir = pbmc_example("pbmc"))
# pbmc <- Seurat::CreateSeuratObject(counts = pbmc_data)
# pbmc <- Seurat::PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")
# pbmc <- Seurat::SCTransform(pbmc, verbose = FALSE)
# pbmc <- Seurat::RunPCA(pbmc, verbose = FALSE)
# pbmc <- Seurat::RunUMAP(pbmc, dims = 1:30, verbose = FALSE)
# pbmc <- Seurat::FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)
# pbmc <- Seurat::FindClusters(pbmc, verbose = FALSE)
# Seurat::DimPlot(pbmc, label = TRUE)


