# scSIGTest

## Author : Elie Robert

Basic package to test the enrichment or the depletion of a given signature in each cell of a scRNAseq dataset. It is based on the detection rate of genes of the signature. Is the number of detected genes in the signature higher or lower to what is expected randomly in randomly selected sets of genes of the same size than the  signature ? Random sets of genes of the same size than the interest set of genes are created to generate a null distribution of detection rate for each cell of the dataset. Then the detection rate of the interest signature is compared to this distribution and a kind of p-value is computed based on this distribution allow to cells in three categories : cells enriched for the signature, cells depleted for the signature or cells without significant enrichment or depetion of the signature.

If you use this package in a publication please cite : 
Arkoun B, Robert E, Boudia F, Mazzi S, Dufour V, Siret A, Mammasse Y, Aid Z, Vieira M, Imanci A, Aglave M, Cambot M, Petermann R, Souquere S, Rameau P, Catelain C, Diot R, Tachdjian G, Hermine O, Droin N, Debili N, Plo I, Malinge S, Soler E, Raslova H, Mercher T, Vainchenker W. Stepwise GATA1 and SMC3 mutations alter megakaryocyte differentiation in a Down syndrome leukemia model. J Clin Invest. 2022 Jul 15;132(14):e156290. doi: 10.1172/JCI156290. PMID: 35587378; PMCID: PMC9282925.
