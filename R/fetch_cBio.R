#' @import cgdsr
# @importFrom cgdsr CGDS getCancerStudies getProfileData getCaseLists
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

#' Create a connection to c-Bioportal API
#' @param url URL of the c-Bioportal website
#' @return A CGDS connexion object
#' @export
#' @seealso \code{\link{cgdsr::CGDS}}
cBioConnect <- function(url="http://www.cbioportal.org/") {
    return(CGDS(url))
}

#' List of c-Bioportal studies
#'
#' Get the list of c-Bioportal cancer studies
#' @param conn A CGDS connexion object
#' @return A data.frame with c-Bioportal studies ids (cancer_study_id), name (name) and description (description)
#' @export
#' @seealso \code{\link{cBioConnect}}, \code{\link{cgdsr::CGDS}}, \code{\link{cgdsr::getCancerStudies}}
listStudies <- function(conn="http://www.cbioportal.org/") {
    if (is.character(conn)) {
        return(getCancerStudies(CGDS(conn)))
    }
    return(getCancerStudies(conn))
}

#' Retrieve a dataset from c-Bioportal
#'
#' Retrieve a dataset from c-Bioportal for all genes that are present on the NaviCell map. Displays whether the NaviCell map genes have been included or not.
#' @param conn A CGDS connexion object
#' @param profile_id List of ids of the profiles we want to retrieve
#' @param case_id ID of the list of cases we want to retrieve
#' @param genes_url URL pointing to the list of genes of interest (in .gmt format)
#' @param method String, either "genes" or "profiles", specifying whether the data must be fetched by genes or by profiles. The result is the same but the "genes" version is more detailed.
#' @return A list of data.frames containing the data for each gene, for all combination experiment x sample
#' @export
#' @seealso \code{\link{cBioStudy}}, \code{\link{importDataSet}}, \code{\link{saveData}}
#' @author Mathurin Dorel \email{mathurin.dorel@@curie.fr}
# TODO replace by the url of the map
cBioDataSet <- function (conn, profile_ids, case_id, genes_url="http://acsn.curie.fr/files/acsn_v1.0.gmt", method=method) {
    if (method == "genes") {
        genes_data = list()
        for (gene in getGenesList(genes_url)) {
            dd = getProfileData(conn, gene, profile_ids, case_id)
            if (nrow(dd) != 0) {
                colnames(dd) = gsub( paste0(st_id, "_"), "", colnames(dd) )
                genes_data[[gene]] = dd
                print(paste(gene, "included"))
            } else {
                print(paste(gene, "not included"))
            }
        }
        print("------------------ Import finished -------------------------")
        return(genes_data)
    } else if (method == "profiles") {
        profiles_data = list()
        genes = getGenesList(genes_url)
        for (prof in profile_ids) {
            pr_code = gsub(paste0(st_id, "_"), "", prof)
            print(paste("Importing", pr_code, "data"))
            profiles_data[[pr_code]] = t(getProfileData(conn, genes, prof, case_id))
        }
        print("------------------ Import finished -------------------------")
        return(profiles_data)
    } else {
        stop("Invalid method, valids are 'profiles' and 'genes'")
    }
}

#' Create NCViz object from a c-Bioportal study
#'
#' Retrieve data and annotations from a c-Bioportal study and select genes that are present on the NaviCell map
#' @param study_id ID of the study to retrieve.
#' @param genes_url URL pointing to the list of genes of interest (in .gmt format)
#' @param nc_url URL of the NaviCell map
#' @param name Name of the dataset. If not provided, the name of the study provided by cBioPortal will be used
#' @param method String, either "genes" or "profiles", specifying whether the data must be fetched by genes or by profiles. The result is the same, however the "genes" version is more detailed but slower and uses more memory.
#' @return An NCviz object containing the data of the study 
#' @export
#' @seealso \code{\link{listStudies}}, \code{\link{cBioDataSet}}, \code{\link{cBioStudy}}
#' @author Mathurin Dorel \email{mathurin.dorel@@curie.fr}
cBioNCviz <- function(study_id, genes_url="http://acsn.curie.fr/files/acsn_v1.0.gmt", nc_url="http://acsn.curie.fr/files/acsn_v1.0.owl", name="", method="genes") {
    all_data = cBioStudy(study_id, genes_url, method=method)
    clinical_data = all_data$annotations
    genes_data = all_data$data

    studies = listStudies(conn)
    if (!study_id %in% studies$cancer_study_id) { stop(paste("The study", study_id, "does not exits on cBioPortal")) }
    if (name == "") { name = studies$name[which(studies$cancer_study_id==study_id)] }
    if (method == "genes") {
        ncviz = NCviz(nc_url=nc_url, cell_type=name, cbio_gene_data=genes_data, annotations=clinical_data)
    } else {
        ncviz = NCviz(nc_url=nc_url, cell_type=name, nc_data=genes_data, annotations=clinical_data)
    }

    return(ncviz)
}

#' Retrieve data and annotations from a c-Bioportal study
#'
#' Retrieve data and annotations from a c-Bioportal study and select genes that are present on the NaviCell map
#' @param study_id ID of the study to retrieve.
#' @param genes_url URL pointing to the list of genes of interest (in .gmt format)
#' @param method String, either "genes" or "profiles", specifying whether the data must be fetched by genes or by profiles. The result is the same but the "genes" version is more detailed.
#' @return A list containing the genes data as a list of dataframe (sample_id * profiling method, list indexed by gene names), and the annotations in a dataframe (sample_id * annotation_type)
#' @export
#' @seealso \code{\link{listStudies}}, \code{\link{cBioDataSet}}, \code{\link{cBioNCviz}}
#' @author Mathurin Dorel \email{mathurin.dorel@@curie.fr}
cBioStudy <- function(study_id, genes_url="http://acsn.curie.fr/files/acsn_v1.0.gmt", method="genes") {
    conn = CGDS("http://www.cbioportal.org/")

    # Retrieve genetic profiles ids of all profiles
    ca_id = paste0(study_id, "_all")
    profiles = getGeneticProfiles(conn, study_id)
    pr_id = profiles$genetic_profile_id

    clinical_data = getClinicalData(conn, ca_id)
    genes_data = cBioDataSet(conn, pr_id, ca_id, genes_url, method=method)

    return(list(data=genes_data, annotations=clinical_data))
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

