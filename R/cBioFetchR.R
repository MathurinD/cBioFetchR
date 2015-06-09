#' Visualize profiling data on a NaviCell map
#'
#' Provides functions to visualize any kind of biological profiling data on a NaviCell map in a biologically relevant manner.
#' Also contains functions to import datasets from cBioportal which format them for NaviCell export.
#'
#' @docType package
#' @name cBioFetchR-package
#' @author Mathurin Dorel \email{mathurin.dorel@@curie.fr}
#' @author Maintainer : Mathurin Dorel \email{mathurin.dorel@@curie.fr}
#' @keywords manip
#' @references TODO
#' @examples
#' conn = CGDS("http://www.cbioportal.org/")
#'
#' #Selection of the study id
#' studies = getCancerStudies(conn)
#' st_id = studies$cancer_study_id[1]
#' 
#' # Fetch data from cBioPortal for ACSN genes
#' fetcher= cBioNCviz(genes_list="http://acsn.curie.fr/files/acsn_v1.1.gmt", cell_type="Acute Myeloid Leukemia", cbio_data=genes_data)
#' # Save Data
#' saveDat(fetcher)
NULL
