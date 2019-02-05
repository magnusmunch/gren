#' 772 sequenced microRNAs in cervical cancer.
#'
#' A dataset containing 772 sequenced microRNA for 32 healthy women and 24 
#' women with pre-cursor lesions for cervical cancer. Additionally, 
#' the conservation levels of the microRNAs are given, coded in three levels:
#' "NotCons", "Mammals", "Broadly", that correspond to microRNAs only found in 
#' humans, found in most mammals, and found in most vertebrates, respectively.
#'
#' @format Contains the following:
#' \describe{
#'   \item{mirCerv}{A matrix with 56 rows and 772 columns}
#'   \item{respCerv}{Factor w/ 2 levels "CIN3","Normal"}
#'   \item{mirCons}{Factor w/ 3 levels "NotCons","Mammals","Broadly"}
#' }
#' @source Novianti, P.W., Snoek, B.C., Wilting, S.M., and van de Wiel, M.A. 
#'         (2017). Better diagnostic signatures from RNAseq data through use of 
#'         auxiliary co-data. Bioinformatics, 33, 1572–1574.
"dataCervical"

#' 2114 sequenced microRNAs, colon cancer treatment response, 
#' 4 clinical covariates, expression level of microRNAs.
#'
#' A dataset containing 2114 sequenced microRNA for 70 treatment responsive 
#' patients and 70 non-responders. Additionally, 4 clinical covariates for the
#' patients and the expression levels of the microRNAs, coded in three levels:
#' "nonExpr", "medExpr", "highExpr".
#'
#' @format Contains the following:
#' \describe{
#'   \item{mirCol}{A matrix with 88 rows and 2114 columns}
#'   \item{unpenCol}{data.frame with	88 observations of  4 variables}
#'   \item{respColl}{Factor w/ 2 levels "Progr","TherBenefit"}
#'   \item{mirExpr}{Factor w/ 3 levels "nonExpr","medExpr","highExpr"}
#' }
#' @source Neerincx, M., Poel, D., Sie, D.L.S., van Grieken, N.C.T., 
#'         Shankaraiah, R.C., van der Wolf - de Lijster, F.S.W., 
#'         van Waesberghe, J.H.T.M., Burggraaf, J.D., Eijk, P.P., Verhoef, C., 
#'         Ylstra, B., Meijer, G.A., van de Wiel, M.A., Buffart, T.E., 
#'         and others. (2018). Combination of a six microRNA expression profile 
#'         with four clinicopathological factors improves response prediction 
#'         to systemic treatment in patients with advanced colorectal cancer. 
#'         Submitted.
#'         
#'         Neerincx, M., Sie, D.L.S., van de Wiel, M.A., van Grieken, N.C.T., 
#'         Burggraaf, J.D., Dekker, H., Eijk, P.P., Ylstra, B., Verhoef, C., 
#'         Meijer, G.A., Buffart, T.E., and others. (2015). MiR expression 
#'         profiles of paired primary colorectal cancer and metastases by 
#'         next-generation sequencing. Oncogenesis, 4, e170.
"dataColon"