#' @import cgdsr
#' @importFrom cgdsr CGDS getCancerStudies getProfileData getCaseList
#' @import RCurl

# Gene lists of the ACSN maps (.owl for the map, .gmt for the gene list)
ACSN_genes = list("global"="http://acsn.curie.fr/files/acsn_v1.0", "survival"="http://acsn.curie.fr/files/survival_v1.0", "apoptosis"="http://acsn.curie.fr/files/apoptosis_v1.0", "cell_motility"="http://acsn.curie.fr/files/emtcellmotility_v1.0", "cell_cycle"="http://acsn.curie.fr/files/cellcycle_v1.0", "DNA_repair"="http://acsn.curie.fr/files/dnarepair_v1.0")

# Retrieve the list of the genes available on the NaviCell map
# TODO replace by the getHugoList of NaviCell R API
# @param url URL pointing to a list of genes
getGenesList <- function(url="http://acsn.curie.fr/files/acsn_v1.0.gmt"){
    acsn_page = getURL(url)
    acsn_list = c()
    for ( ll in unlist(strsplit(acsn_page, "\n")) ) {
        acsn_list = c(acsn_list, unlist(strsplit(ll, "\t"))[-c(1, 2)])
    }
    acsn_genes = unique(acsn_list)
}

#' Retrieve a dataset from c-Bioportal
#'
#' Retrieve a dataset from c-Bioportal and select the genes that are present on the NaviCell map
#' @param conn A CGDS connexion to c-Bioportal
#' @param profile_id List of ids of the profiles we want to retrieve
#' @param case_id ID of the list of cases we want to retrieve
#' @param genes_url URL pointing to the list of genes of interest (in .gmt format)
#' @return A list of data.frames containing the data for each gene, for all combination experiment x sample
#' @export
#' @seealso \code{\link{cBioStudy}}, \code{\link{importDataSet}}, \code{\link{saveData}}
#' @author Mathurin Dorel \email{mathurin.dorel@@curie.fr}
# TODO replace by the url of the map
cBioDataSet <- function (conn, profile_id, case_id, genes_url="http://acsn.curie.fr/files/acsn_v1.0.gmt") {
    genes_data = list()
    for (gene in getGenesList(genes_url)) {
        dd = getProfileData(conn, gene, profile_id, case_id)
        if (nrow(dd) != 0) {
            colnames(dd) = gsub( gsub("[^_]$", "", case_id), "", colnames(dd) )
            genes_data[[gene]] = dd
            print(paste(gene, "included"))
        } else {
            print(paste(gene, "not included"))
        }
    }
    return(genes_data)
}

#' Retrieve the data from a c-Bioportal study
#'
#' Retrieve from a c-Bioportal study and select genes that are present on the NaviCell map
#' @param study_id ID of the study to retrieve.
#' @param genes_url URL pointing to the list of genes of interest (in .gmt format)
#' @return A list of data.frames containing the data for each gene, for all combination experiment x sample
#' @export
#' @seealso \code{\link{cBioDataSet}}, \code{\link{importDataSet}}, \code{\link{saveData}}
#' @author Mathurin Dorel \email{mathurin.dorel@@curie.fr}
cBioStudy <- function(study_id, genes_url="http://acsn.curie.fr/files/acsn_v1.0.gmt") {
    conn = CGDS("http://www.cbioportal.org/")

    # Retrieve genetic profiles ids of all profiles
    ca_id = paste0(study_id, "_all")
    profiles = getGeneticProfiles(conn, study_id)
    pr_id = profiles$genetic_profile_id

    genes_data = cBioDataSet(conn, pr_id, ca_id, genes_url)
    return(genes_data)
}

# Save data set in a file, avoid having to download everything every time
# TODO
saveDataSet <- function(fname, dataset) {
    print("Not implemented yet")
    for (gene in names(data_set)) {
    }
}

# Import a data set from a file
# TODO
importDataSet  <- function(fname) {
    fcontent = readLines(fname)
}

