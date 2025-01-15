#' @title GenicSignature_example
#'
#' @description Get path to molecular signature files used as example to illustrate how scSIGTest works.
#'
#' scSIGTest comes bundled with some example files in its `inst/extdata`
#' directory. This function allows to get the path of files containing molecular signatures that can be used as input of the scSIGTest method.
#'
#' @param data  Name of the available signature that the user wants to test.
#' @export
#' @examples
#' 
#' GenicSignature_example("platelet_signature")
#'


GenicSignature_example <- function(data = "platelet_signature") {
  if (data == "platelet_signature"){
#opening bracket if data is equal to platelet_signature.
    
    print(system.file("extdata","HAY_BONE_MARROW_PLATELET.v2024.1.Hs.grp", package = "scSIGTest", mustWork = TRUE))
    
    system.file("extdata","HAY_BONE_MARROW_PLATELET.v2024.1.Hs.grp", package = "scSIGTest", mustWork = TRUE)


#closing bracket if data is equal to platelet signature.
   }
}
