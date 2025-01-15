#' @title ProcessPBMC
#'
#' @description Perform a basic preprocessing of the PBMC dataset which is used in this package to provide some examples of use of some of the defined function.
#'
#' @importFrom Seurat Read10X CreateSeuratObject PercentageFeatureSet SCTransform RunPCA RunUMAP FindNeighbors FindClusters DimPlot CellCycleScoring
#' 
#' 
#'
#' @export
#'
#' @examples 
#'
#' require(Seurat)
#' 
#' 
#'
#' Seurat::DimPlot(processPBMC)
#' 
#'
#' @return A Seurat object containing the PBMC already basically processed and ready to use as an example to illustrate how some functions of the package work. 
#'
#' 
#'


ProcessPBMC <- function(){#opening bracket of the function ProcessPBMC

pbmc_data <- Seurat::Read10X(data.dir = pbmc_example())

pbmc <- Seurat::CreateSeuratObject(counts = pbmc_data)

pbmc <- Seurat::PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")

pbmc <- Seurat::SCTransform(pbmc, verbose = FALSE)

pbmc <- Seurat::RunPCA(pbmc, verbose = FALSE)

pbmc <- Seurat::RunUMAP(pbmc, dims = 1:30, verbose = FALSE)

pbmc <- Seurat::FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)

pbmc <- Seurat::FindClusters(pbmc, verbose = FALSE)

Seurat::DimPlot(pbmc, label = TRUE)

s.genes <- Seurat::cc.genes$s.genes

g2m.genes <- Seurat::cc.genes$g2m.genes

pbmc <- Seurat::CellCycleScoring(pbmc, 
s.features = s.genes, 
g2m.features = g2m.genes, 
set.ident = TRUE)

return(pbmc)

# closing bracket of the function ProcessPBMC
}


