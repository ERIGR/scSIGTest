
#' @title pvaldistri
#' @description compute the p-value associated to a detection rate for a biological signature under the null hypothesis that this signature is equivalent to a random signature.  
#' 
#' @param distridetecrate numeric vector of distribution rate of the set of random signature that have the same size that the signature of interest.
#'
#' @param detectionratesignature (numeric) detection rate of the signature of interest
#'
#' @return return a numeric value which represent the probability that this detection rate can be observed with a random signature.
#'  If the value is negative, it only means that the detection rate is lower than the median detection rate for random signatures of the same size than the tested signature. 
#' In the opposite case, the value is positive.
#' 
#' @import Seurat 
#'
#' @export 
#'
#' @examples
#' require(Seurat)
#' 
#'
#' platelet <- read.table(GenicSignature_example("platelet_signature"),header=FALSE)[[1]]
#'
#' used_expr_mat <- processPBMC@assays$SCT@counts
#' echantillonnage <- replicate(100,sample(rownames(used_expr_mat),
#' length(intersect(rownames(used_expr_mat),platelet))),simplify = FALSE)
#'
#' echantillonnage <- lapply(echantillonnage,function(x){which(rownames(used_expr_mat) %in% x)})
#'
#' 
#'
#' test <- apply(used_expr_mat,2,function(z){ droplets <- z ; 
#' percentages <- lapply(echantillonnage,function(y){
#' (length(which(droplets[y]>0))/length(droplets[y]))*100}) 
#' return(percentages)
#' })
#'
#' test <- lapply(test,unlist)
#'
#'
#' names(test) <- colnames(used_expr_mat)
#'
#'
#' mean_random_detection_rate_per_cell <- lapply(test,mean) 
#'
#' mean_random_detection_rate_per_cell <- unlist(mean_random_detection_rate_per_cell)
#' detection_rate_signature <- apply(used_expr_mat,2,
#' function(z){return((length(
#' which(z[which(rownames(used_expr_mat) %in% intersect(rownames(used_expr_mat),platelet))]>0)
#' )/length(z[which(rownames(used_expr_mat) %in% intersect(rownames(used_expr_mat),platelet))]))*100)})
#'
#'
#'pvaltoret <- purrr::pmap(.l = list(random = test , signature = detection_rate_signature),
#' .f = function(random,signature){return(pvaldistri(random,signature))})
#'
#'
#'
#'
#'
#' pvaltoret <- unlist(pvaltoret)
#'
#' fdrtoret <- stats::p.adjust(abs(pvaltoret),method = "BH", n = length(pvaltoret))
#'
#'



pvaldistri <- function(distridetecrate,detectionratesignature){
  
#setting of boundaries for the density.  
  verysmall <- 0
  verylarge <- 200000
  
#estimation fonction de densite de distribution des valeurs des taux de detection des signatures aleatoires de la meme taille que la signature d interet
# estimation of the density distribution of the values of the detection rate of random signatures of the same size than the interest signature.
#Ici on genere un vecteur X contenant les valeurs de depart dans l ordre croissant (vecteur distridetecrate dans l ordre croissant) mais bornees par les valeurs 0 et 200000. 
# on estime le nombre d'occurences 
X = c(verysmall, stats::density(distridetecrate)$x, verylarge)

# estimation de l'estimation  de la fonction de repartition des valeurs des taux de detection des signatures aleatoires 
#de la meme taille que la signature d'interet 

SUM = cumsum(stats::density(distridetecrate)$y)
# Vector normalized on 1 (normalisation des valeurs sur l'intervalle [0;1] pour avoir une distribution de probabilites cad que AUC = 1                                                                                                                                                      )
Y = c( 0, SUM/max(SUM) , 1)

# Valeur x de la distribution pour laquelle y = 0.5 (moitie  de l'AUC de la fonction de distribution 0.5). 
# Value x of the distribution for which y = 0.5 (half of the area under curve) 
# on supprime le premier element et le dernier element du vecteur Y d ou le [-c(1,length(Y))]. De meme pour X d'où le [-c(1,length(X))].
# On calcule la valeur de X pour laquelle on a la moitié de l aire sous la courbe.
#tester fonction avec le code ci-dessous car probablement le code actuel ne fonctionne que si la distribution des signatures aléatoires suit une loi normale
#if detectionratesignature est  superieure à la mediane
#revdensSUM = cumsum(rev(stats::density(distridetecrate)$y))
#revdensY = c(0,revdensSUM/max(revdensSUM),1)
#pval = stats::approxfun(rev(X),revdensY)(detectionrates)
#si c est le cas inverse pval <- stats::approxfun(X,Y)(detectionratesignature); return(-pval)

middistribution <- stats::approxfun(Y[-c(1,length(Y))],X[-c(1,length(X))])(0.5)
if(detectionratesignature >= middistribution){
  
  Yrev <- rev(Y)
  
  pval <- stats::approxfun(X,Yrev)(detectionratesignature)
  if(pval==0){#instructions si pvalue estimee à 0 par la fonction de densite 
    
    pval <- 10^-10

    
  }
  return(pval)
  
}
else{
    
  pval <- stats::approxfun(X,Y)(detectionratesignature)
    if(pval==0){
    
    pval <- 10^-10
    
  }
  return(-pval)
  
  
}
  
  
  
}



# pbmc_data <- Seurat::Read10X(data.dir = pbmc_example())
#
# 
# pbmc <- Seurat::CreateSeuratObject(counts = pbmc_data)
# pbmc <- Seurat::PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")
#' pbmc <- Seurat::SCTransform(pbmc, verbose = FALSE)
#' pbmc <- Seurat::RunPCA(pbmc, verbose = FALSE)
#' pbmc <- Seurat::RunUMAP(pbmc, dims = 1:30, verbose = FALSE)
#' pbmc <- Seurat::FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)
#' pbmc <- Seurat::FindClusters(pbmc, verbose = FALSE)
#' Seurat::DimPlot(pbmc, label = TRUE)
#'load(processPBMC_example(data= "processPBMC"))
