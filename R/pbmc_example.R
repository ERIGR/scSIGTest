#' @title pbmc_example
#'
#' @description Get path to files to run the scSIGTest method on the pbmc dataset.
#'
#' scSIGTest comes bundled with some example files in its `inst/extdata`
#' directory. This function allows to load files to run the scSIGTest method for example on the pbmc dataset
#'
#' @param data  Name of the example dataset.
#' @export
#' @examples
#' 
#' pbmc_example("pbmc")
#'


pbmc_example <- function(data = "pbmc") {
  if (data == "pbmc"){
    print(system.file("extdata","filtered_gene_bc_matrices", package = "scSIGTest", mustWork = TRUE))
    system.file("extdata","filtered_gene_bc_matrices", package = "scSIGTest", mustWork = TRUE)
   }
}
