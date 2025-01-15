#' @title df_rep_of_each_level_a_factor_in_each_level_another_factor
#'
#' @description This function allows to make a barplot representation of the proportions of the different level of a factor into another factor.
#'
#' @param df a dataframe (class data.frame) with at least two columns of class character or class factor. These columns should give the level of a specific factor for each sample (row) of the dataframe. 
#'
#' @param statcol an element of class character that should correspond to the name of a column of class character or factor of the dataframe df.
#'
#' @param otherfaccol an element of class character that should correspond to the name of a column of class character or factor of the datarame df. The value of the parameter otherfaccol 
#' should not be the same as the value of the parameter statcol. 
#'
#' @param optionrep an element of class character that should correspond either to percentage or brut. In the first case the distribution of the levels of the factor otherfaccol inside the levels of
#' the factor statcol is represented in percentage otherwise when the value of optionrep is brut representing  the real number of samples.
#'
#' @param personnal_colors an element of class logical. if TRUE  the colors to represent the levels of the factor otherfaccol are provided by the user otherwise different colors choosen according to the number of levels of the factor otherfaccol will be used (determined with the package chameleon).   
#'
#' @param colchoose  a vector of class character. It provides the color that will be associated to each level of the factor otherfaccol if the parameter personnal_colors is set to TRUE.
#'
#' @return generate a ggplot2 barplot representation that gives the distribution of the levels of the factor otherfaccol into the level of the factor statcol. 
#'
#' @importFrom chameleon distinct_colors
#'
#' @importFrom ggplot2 ggplot
#' 
#' @importFrom ggplot2 geom_bar
#'
#' @importFrom ggplot2 theme_bw
#' 
#' @importFrom ggplot2 coord_flip
#'
#' @importFrom ggplot2 scale_fill_discrete 
#'
#' @importFrom ggplot2 scale_fill_manual
#'
#' @importFrom ggplot2 theme
#'
#' @export
#' 
#' @examples 
#' require(Seurat)
#' 
#'
#' 
#' df_rep_of_each_level_a_factor_in_each_level_another_factor(processPBMC@meta.data,
#'statcol = "seurat_clusters",
#'otherfaccol = "Phase",
#'personnal_colors = FALSE,
#'optionrep = "percentage")
#'
#' df_rep_of_each_level_a_factor_in_each_level_another_factor(processPBMC@meta.data,
#' statcol = "seurat_clusters",
#' otherfaccol = "Phase",
#' personnal_colors = TRUE,
#' colchoose = c("S" = "blue","G1" = "green","G2M" = "red"),
#' optionrep = "percentage")
#'
#'

df_rep_of_each_level_a_factor_in_each_level_another_factor <- function(df,statcol,otherfaccol,optionrep,personnal_colors=FALSE,colchoose = NULL){
##### opening bracket of the function df_rep_of_each_level_a_factor_in_each_level_another_factor

if(inherits(optionrep,"character") == FALSE){
  
  stop("The value given for the argument optionrep must be of class character (must be either brut or percentage)")
  
  
}  

if(optionrep %in% c("brut","percentage") == FALSE){
  
  stop("The value given for the argument optionrep must be either the character chain brut or percentage (must be of class character)")
  
  
}   
  
if(inherits(personnal_colors,"logical") == FALSE){
    
    stop("The value given for the argument personnal_colors must be of class logical")
    
  }
  
  

  
  # Generation of all possible combinations of levels of the two interest factors.
  comb <- base::expand.grid( unique(df[[statcol]]),unique(df[[otherfaccol]]))
  
  colnames(comb) <- c(statcol,otherfaccol)



#storage of the number of lines of the dataframe that correspond to each of the possible combination of the levels of the two interest factors. The result is stored in a vector that correspond to a column named effectif of the dataframe effectifcombi that corresponds to the dataframe comb with the column effectif added. 
  
  effectifcombi <- data.frame(
comb,
effectif = unlist(lapply(1:nrow(comb),function(z){
    
   # print(df[which(df[[otherfaccol]]==comb[[otherfaccol]][z] & comb[[statcol]] == df[[statcol]][z]),])  
  
    return(nrow(df[which(df[[otherfaccol]]==comb[[otherfaccol]][z] & df[[statcol]] == comb[[statcol]][z]),]))
  
}))
)
  
 
############ computation of the percentage of each level of the factor otherfaccol in each level of statcol 
effectifcombi <- data.frame(effectifcombi,
percentage = unlist(lapply(1:nrow(effectifcombi),function(z){
  
  return(
    ((effectifcombi$effectif[z])/(sum(effectifcombi$effectif[which(effectifcombi[[statcol]]==effectifcombi[[statcol]][z])])))*100
    )
  
  
  
}))

)
  
  
if(personnal_colors==FALSE){  

col_otherfaccol <- chameleon::distinct_colors(length(unique(df[[otherfaccol]])))[["name"]]

}else{
  
  if(length(colchoose) != length(unique(df[[otherfaccol]]))){
    
    stop(paste0("The number of colors in the vector colchoose must be equal to the number of levels of the factor ",otherfaccol, " in the dataframe given as value of the argument df "))
    }
  else{
    
  col_otherfaccol <- colchoose
  
  }
  
  
  
}
  
  
  
################################################## représentation ggplot ###################################################


effectifcombi[[otherfaccol]] <- as.character(effectifcombi[[otherfaccol]])

if(optionrep=="percentage"){


print(
     ggplot2::ggplot(effectifcombi,ggplot2::aes_string(x=statcol, y="percentage", fill=otherfaccol))+ ggplot2::geom_bar(stat="identity")+
     ggplot2::theme_bw()+
     ggplot2::coord_flip()+
     ggplot2::scale_fill_discrete(breaks = sort(unique(df[[otherfaccol]])))+
     ggplot2::scale_fill_manual(values= col_otherfaccol)+
     ggplot2::theme(axis.text= ggplot2::element_text(size=20))
 
)
}else{
  
print(
     ggplot2::ggplot(effectifcombi,ggplot2::aes_string(x=statcol, y="effectif", fill= otherfaccol))+ ggplot2::geom_bar(stat="identity")+
     ggplot2::theme_bw()+
     ggplot2::coord_flip()+
     ggplot2::scale_fill_discrete(breaks = sort(unique(df[[otherfaccol]])))+
     ggplot2::scale_fill_manual(values= col_otherfaccol)+
     ggplot2::theme(axis.text= ggplot2::element_text(size=20))
     )
  
  
}  
 
  
##### closing bracket of the function df_rep_each_level_a_factor_in_each_level_another_factor  
}
 


 
#pbmc_data <- Seurat::Read10X(data.dir = pbmc_example())
#pbmc <- Seurat::CreateSeuratObject(counts = pbmc_data)
#pbmc <- Seurat::PercentageFeatureSet(pbmc, pattern = "^MT-", 
#col.name = "percent.mt")
#pbmc <- Seurat::SCTransform(pbmc, verbose = FALSE)
#pbmc <- Seurat::RunPCA(pbmc, verbose = FALSE)
#pbmc <- Seurat::RunUMAP(pbmc, dims = 1:30, verbose = FALSE)
#pbmc <- Seurat::FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)
#pbmc <- Seurat::FindClusters(pbmc, verbose = FALSE)
#Seurat::DimPlot(pbmc, label = TRUE)
#s.genes <- Seurat::cc.genes$s.genes
#g2m.genes <- Seurat::cc.genes$g2m.genes
#pbmc <- Seurat::CellCycleScoring(pbmc, 
#s.features = s.genes, 
#g2m.features = g2m.genes, 
#set.ident = TRUE)
#load(processPBMC_example(data = "processPBMC")) 


