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
    for (ii in 1:length(acsn_genes)) {
        acsn_genes[ii] = gsub(" ", "_", acsn_genes[ii])
    }
    return(acsn_genes)
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
#' Get the list of c-Bioportal cancer studies. Allow the specification of cancer types to limit the search.
#' @param conn A CGDS connexion object
#' @param query A string specifying a subtype of cancer to look for (leave to "all" to list all available studies)
#' @param case_sensitive A boolean specifying whether the query should be case sensitive (does not apply if query == "all")
#' @return A data.frame with c-Bioportal studies ids (cancer_study_id), name (name) and description (description)
#' @export
#' @seealso \code{\link{cBioConnect}}, \code{\link{cgdsr::CGDS}}, \code{\link{cgdsr::getCancerStudies}}
listStudies <- function(conn="http://www.cbioportal.org/", query="all", case_sensitive=FALSE) {
    if (is.character(conn)) {
        studies = getCancerStudies(CGDS(conn))
    } else {
        studies = getCancerStudies(conn)
    }

    if (query == "all") {
        return(studies)
    } else {
        result = data.frame()

        for ( rr in as.numeric(rownames(studies)) ) {
            if ( (case_sensitive && grepl(query, studies$name[rr])) || grepl(tolower(query), tolower(studies$name[rr])) ) {
                result = rbind(result, studies[rr,])
            }
        }
        return(result)
    }
}

#' Retrieve a dataset from c-Bioportal
#'
#' Retrieve a dataset from c-Bioportal for all genes that are present on the NaviCell map. Displays whether the NaviCell map genes have been included or not.
#' @param conn A CGDS connexion object
#' @param study_id ID of the study to retrieve.
#' @param profile_id List of ids of the profiles we want to retrieve
#' @param case_id ID of the list of cases we want to retrieve
#' @param genes_list URL pointing to the list of genes of interest (in .gmt format), or a list of genes HUGO identifiers
#' @param method String, either "genes" or "profiles", specifying whether the data must be fetched by genes or by profiles. The result is the same, however the "genes" version is more detailed but much slower and uses more memory.
#' @details The method "profiles" does not work it there are too many samples
#' @return The format depends on the method :\cr
#' "genes" : a list (indexed by gene names) of dataframes (sample_id * profiling method)\cr
#' "profiles" : a list (indexed by profiling methods) of dataframes (genes * samples), and the annotations in a dataframe (sample_id * annotation_type)\cr
#' @export
#' @seealso \code{\link{cBioStudy}}, \code{\link{importDataSet}}, \code{\link{saveData}}
#' @author Mathurin Dorel \email{mathurin.dorel@@curie.fr}
# TODO replace by the url of the map
cBioDataSet <- function (conn, study_id, profile_ids, case_id, genes_list="http://acsn.curie.fr/files/acsn_v1.0.gmt", method="profiles") {
    if (length(genes_list) == 1) {
        genes_list = getGenesList(genes_list)
    } else {
        genes_list = genes_list
    }
    if (method == "genes") {
        imported = 0
        genes_data = list()
        for (gene in genes_list) {
            dd = getProfileData(conn, gene, profile_ids, case_id)
            if (nrow(dd) != 0) {
                colnames(dd) = gsub( paste0(study_id, "_"), "", colnames(dd) )
                genes_data[[gene]] = dd
                print(paste(gene, "included"))
                imported = imported+1
            } else {
                print(paste(gene, "not included"))
            }
        }
        print("------------------ Import finished -------------------------")
        print(paste0(imported, "/", length(genes_list), " genes successfully imported"))
        return(genes_data)
    } else if (method == "profiles") {
        profiles_data = list()
        for (prof in profile_ids) {
            pr_code = gsub(paste0(study_id, "_"), "", prof)
            print(paste("Importing", pr_code, "data"))
            caseList = getCaseLists(conn, study_id)
            nSamples = length(unlist(strsplit(caseList$case_ids[caseList$case_list_id == case_id], " ")))
            MAX_GENES = floor(90000 / nSamples) # Maximum number of genes per request, otherwise the server raises an error
            if (MAX_GENES < 1) {
                stop("Too many samples to use \"profiles\" method to get the data, use \"genes\" instead")
            }
            if (length(genes_list) < MAX_GENES) {
                dd = getProfileData(conn, genes_list, prof, case_id)
            } else {
                # If the genes list is too long, the URL gets too long and a 414 error is raised
                sub_lists = list()
                sub_lists[[ceiling(length(genes_list)/MAX_GENES)]] = 0
                lsl = 0
                while (length(genes_list) > (lsl+1) * MAX_GENES) {
                    lsl = lsl+1
                    sub_lists[[lsl]] = genes_list[(1+MAX_GENES*(lsl-1)):(MAX_GENES*lsl)]
                }
                sub_lists[[lsl+1]] = genes_list[(1+MAX_GENES*lsl):length(genes_list)]
                dd = matrix()
                for (subl in sub_lists) {
                    pdd = getProfileData(conn, subl, prof, case_id)
                    dd = cbind(dd, pdd)
                }
                dd = dd[,-1]
            }
            profiles_data[[pr_code]] = t(dd)
        }
        print("------------------ Import finished -------------------------")
        print(paste0(nrow(profiles_data[[pr_code]]), "/", length(genes_list), " genes successfully imported"))
        return(profiles_data)
    } else {
        stop("Invalid method, valids are 'profiles' and 'genes'")
    }
}

#' Create NCviz object from a c-Bioportal study
#'
#' Retrieve data and annotations from a c-Bioportal study and select genes that are present on the NaviCell map
#' @param study_id ID of the study to retrieve.
#' @param genes_list URL pointing to the list of genes of interest (in .gmt format), or a list of genes HUGO identifiers
#' @param nc_url URL of the NaviCell map
#' @param name Name of the dataset. If not provided, the name of the study provided by cBioPortal will be used
#' @param method String, either "genes" or "profiles", specifying whether the data must be fetched by genes or by profiles. The result is the same, however the "genes" version is more detailed but much slower and uses more memory.
#' @param url URL to the CGDS API
#' @return An NCviz object containing the data of the study 
#' @export
#' @seealso \code{\link{listStudies}}, \code{\link{cBioDataSet}}, \code{\link{cBioStudy}}
#' @author Mathurin Dorel \email{mathurin.dorel@@curie.fr}
cBioNCviz <- function(study_id, genes_list="http://acsn.curie.fr/files/acsn_v1.0.gmt", nc_url="http://acsn.curie.fr/files/acsn_v1.0.owl", name="", method="profiles", url="http://www.cbioportal.org/") {
    all_data = cBioStudy(study_id, genes_list, method=method, url=url)
    clinical_data = all_data$annotations
    genes_data = all_data$data
    studies = getCancerStudies(cBioConnect(url))

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
#' @param genes_list URL pointing to the list of genes of interest (in .gmt format), or a list of genes HUGO identifiers
#' @param method String, either "genes" or "profiles", specifying whether the data must be fetched by genes or by profiles. The result is the same but the "genes" version is more detailed.
#' @param url URL to the CGDS API
#' @return A list containing the annotations in a dataframe (sample_id * annotation_type) and the data. The format of the data depends on the method :\cr
#' "genes" : a list (indexed by gene names) of dataframes (sample_id * profiling method)\cr
#' "profiles" : a list (indexed by profiling methods) of dataframes (genes * samples), and the annotations in a dataframe (sample_id * annotation_type)\cr
#' @export
#' @seealso \code{\link{listStudies}}, \code{\link{cBioDataSet}}, \code{\link{cBioNCviz}}
#' @author Mathurin Dorel \email{mathurin.dorel@@curie.fr}
cBioStudy <- function(study_id, genes_list="http://acsn.curie.fr/files/acsn_v1.0.gmt", method="profiles", url="http://www.cbioportal.org/") {
    conn = cBioConnect(url)
    studies = listStudies(conn)
    if (!study_id %in% studies$cancer_study_id) { stop(paste("The study", study_id, "does not exist on", url)) }

    # Retrieve genetic profiles ids of all profiles
    ca_id = paste0(study_id, "_all")
    profiles = getGeneticProfiles(conn, study_id)
    pr_id = profiles$genetic_profile_id

    clinical_data = list()
    tryCatch({
        clinical_data = getClinicalData(conn, ca_id)
    }, error = function(e) {
        warning("Error when trying to retrieve clinical data, result is empty clinical data")
    })
    genes_data = cBioDataSet(conn, study_id, pr_id, ca_id, genes_list, method=method)

    return(list(data=genes_data, annotations=clinical_data))
}

getLine <- function(ll, split_char="\t") {
    return(gsub('\"', '', unlist(strsplit(ll, split_char))))
}

#' Import a study from a file
#'
#' Import a study from a file with its annotations. A study consists of several experiments on the same set of genes and samples.
#' For a study on m patients and n genes, the file must be formated as "M experiment_name", followed by a line with the list of samples ids, followed by n lines "genes_name data".
#' The annotations start by "ANNOTATIONS study_name", followed by a line with the list of annotations names, followed by m lines "sample_name annotation".
#' @param fname Name of the file
#' @return A list (indexed by profiling methods) of dataframes (genes * samples), and the annotations in a dataframe (sample_id * annotation_type).
#' @export
#' @seealso \code{\link{cBioStudy}}, \code{\link{saveData}}
#' @author Mathurin Dorel \email{mathurin.dorel@@curie.fr}
importStudy <- function(fname) {
    ff = readLines(fname)
    ll = 1

    all_data = list()
    while (ll <= length(ff)) {
        if (grepl("^M ", ff[ll])) {
            # Import data
            method = gsub("^M ", "", ff[ll])
            samples = getLine(ff[ll+1])

            ll = ll+2
            profile_data = data.frame()
            genes = c()
            while(!grepl("^M ", ff[ll]) && !grepl("^ANNOTATIONS", ff[ll]) && ll <= length(ff)) {
                dl = getLine(ff[ll])
                genes = c(genes, dl[1])
                profile_data = rbind( profile_data, suppressWarnings(as.numeric(dl[-1])) )

                ll = ll+1
            }
            colnames(profile_data) = samples
            rownames(profile_data) = genes
            all_data[[method]] = profile_data
        } else if (grepl("^ANNOTATIONS", ff[ll])) {
            # Import annotations
            annotations_names = getLine(ff[ll+1])
            ## Each annotations has to be imported as an independent list, because rbind on data.frame creates factors and puts values different from row 1 to NA
            annot_lists = list()
            for (name in annotations_names) {
                annot_lists[[name]] = c(0) # c() would not create the list slot
            }

            ll = ll+2
            samples = c()
            while(!grepl("^M ", ff[ll]) && !grepl("^ANNOTATIONS", ff[ll]) && ll <= length(ff)) {
                al = getLine(ff[ll])
                samples = c(samples, al[1])
                # Gather each annotation for this sample
                spl = al[1]
                for (id in 2:length(al)) {
                    annot = al[id]
                    if ( !is.na(suppressWarnings(as.numeric(annot))) ) { annot = as.numeric(annot) }
                    annot_lists[[id-1]] = c(annot_lists[[id-1]], annot)
                }

                ll = ll+1
            }
            for (nn in names(annot_lists)) {
                annot_lists[[nn]] = annot_lists[[nn]][-1]
            }
            annotations = data.frame(annot_lists)
            rownames(annotations) = samples
            colnames(annotations) = annotations_names
        }
    }
    return(list(data=all_data, annotations=annotations))
}

#' Import a study from a file to an NCviz object
#'
#' Import a study from a file in an NCviz object.
#' @param nc_url URL of the NaviCell server
#' @param cell_type Name of the data. If cell_type == "guess", it will be deduced from the file name by default.
#' @return An NCviz object.
#' @export
#' @rdname importStudy
importNCviz <- function(fname, nc_url="http://acsn.curie.fr/files/acsn_v1.0.owl", cell_type="guess") {
    dd = importStudy(fname)

    # Try to guess cell_type from the name of the file
    if (cell_type == "guess") {
        cell_type = gsub("_", " ", sub(".txt$", "", gsub("([^/]/)*", "", fname)))
    }

    ncviz = NCviz(nc_url=nc_url, cell_type=cell_type, annotations=dd$annotations, nc_data=dd$data)
    return(ncviz)
}

