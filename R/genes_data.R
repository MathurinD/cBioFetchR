#' @import methods
# @import RNaviCell

check_dataset <- function(object) {
    errors = c() 

    cbio_gene = (length(object@cbio_gene_data) != 0)
    cbio_profile = (length(object@cbio_profile_data) != 0)
    nc_data = (length(object@nc_data) != 0)
    if ((cbio_gene && cbio_profile) || (cbio_gene && nc_data) || (cbio_profile && nc_data)) {
        errors = c(errors, "Only one dataset must be provided")
    }
    if (!cbio_gene && !cbio_profile && !nc_data) {
        errors = c(errors, "One dataset must be provided")
    }

    if (length(object@cell_type) == 0) {
        errors = c(errors, "A cell_type argument must be provided to name the dataset")
    }
    if (length(object@nc_url) == 0) {
        errors = c(errors, "A NaviCell map URL must be provided")
    }

    if (length(errors) == 0) TRUE else errors
}

#' Biological profilin data visualisation
#'
#' S4 class for biological profiling data visualisation. Contains statistics and functions to communicate with a NaviCell API.
setClass("NCviz",
     representation(
            nc_url="character",
            cell_type="character", # Name of the dataset
            annotations="list",
            cbio_gene_data="list", # Data indexed by genes
            cbio_profile_data="list", # Data indexed by profiles
            nc_data="list", # Data indexed by experiment type
            verbose="logical"
            ),
     #prototype(nc_url="", cell_type="", cbio_data=list(), nc_data=list()),
     validity=check_dataset
     )

#' Constructor for the class NCviz. NCviz objects provide functions to visualize profiling data on a NaviCell map
#' @param nc_url URL of the NaviCell server
#' @param cell_type Name of the dataset to pass to the NaviCell server
#' @param annotations Cases annotations (clinical data)
#' @param cbio_data Data loaded with cBioStudy (incompatible with nc_data)
#' @param nc_data Data under the NaviCell format, a list indexed by experiment type containing data arranged in a dataframe with genes HUGO indentifiers as rownames and samples id as colnames (incompatible with cbio_data)
#' @export
#' @seealso cBioStudy
#' @author Mathurin Dorel \email{mathurin.dorel@@curie.fr}
#' @rdname NCviz-class
NCviz <- function(nc_url="", cell_type="", annotations=list(), cbio_gene_data=list(), cbio_profile_data=list(), nc_data=list(), verbose=TRUE) {
    return(new( "NCviz", nc_url=nc_url, cell_type=cell_type, annotations=annotations, cbio_gene_data=cbio_gene_data, cbio_profile_data=cbio_profile_data, nc_data=nc_data, verbose=verbose ))
}

# TODO Determine if some extra data are recquired
setMethod("initialize",
          "NCviz",
            function(.Object, nc_url, cell_type, annotations, cbio_gene_data, cbio_profile_data, nc_data, verbose=TRUE) {
                print(paste("Creation of an", "NCviz", "object"))
                .Object@nc_url = nc_url
                .Object@cell_type = cell_type
                .Object@cbio_gene_data = cbio_gene_data
                .Object@cbio_profile_data = cbio_profile_data
                .Object@nc_data = nc_data
                .Object@verbose = verbose
                .Object@annotations = annotations
                validObject(.Object) # Validate object

                # Copy data in the other format
                if (length(.Object@nc_data) == 0) {
                    if (length(.Object@cbio_gene_data) != 0) {
                        print("Converting form c-Bioportal gene list format to NaviCell format")
                        for (experiment in colnames(.Object@cbio_gene_data[[1]])) {
                            print(paste("Inserting experiment", experiment))
                            .Object@nc_data[[experiment]] = list()
                            for (gene in names(.Object@cbio_gene_data)) {
                                .Object@nc_data[[experiment]][[gene]] = .Object@cbio_gene_data[[gene]][[experiment]]
                            }
                            .Object@nc_data[[experiment]] = as.data.frame(t( as.data.frame(.Object@nc_data[[experiment]]) ))
                            colnames(.Object@nc_data[[experiment]]) = rownames(.Object@cbio_gene_data[[1]])
                        }
                    } else {
                        print("Converting form c-Bioportal profiles list format to NaviCell format")
                        for (experiment in names(.Object@cbio_profile_data)) {
                            .Object@nc_data[[experiment]] = t(.Object@cbio_profile_data[[experiment]])
                        }
                    }
                } else {
                    .Object@cbio_gene_data[[length(rownames(.Object@nc_data[[1]]))]] = 0
                    print("Data already in NaviCell format, no conversion needed") # as it is useless
                }

                # Add a group with all patients for NaviCell group visualisation, if not already present
                # Add NAs when annotations is not present for a patient
                if (! "all" %in% colnames(.Object@annotations) ) {
                    patients = colnames(.Object@nc_data[[1]])
                    group_all = data.frame( "all"=rep(1, length(patients)) )
                    if (length(annotations) == 0) {
                        rownames(group_all) = patients
                        .Object@annotations = group_all
                    } else {
                        if (all( patients %in% rownames(.Object@annotations) )) {
                            .Object@annotations = data.frame(group_all, annotations)
                        } else { # Add missing patients annotations as NA
                            to_add = patients[!patients %in% rownames(.Object@annotations)]
                            missing_annots = matrix(rep(NA, length(to_add) * ncol(.Object@annotations)), nrow=length(to_add))
                            rownames(missing_annots) = to_add
                            colnames(missing_annots) = colnames(.Object@annotations)
                            rownames(group_all) = c(rownames(.Object@annotations), rownames(missing_annots))
                            .Object@annotations = data.frame(group_all, rbind(.Object@annotations, missing_annots))
                        }
                    }
                }

                return(.Object)
            }
          )


#setGeneric(addData) # Add data tables
#setGeneric()

#TODO
NCdisplay <- function(obj) {
}
setGeneric("NCdisplay", NCdisplay)
#' Display data on a NaviCell map
#'
#' Display the data on a NaviCell map using an optimized visualisation setup
#' @param obj NCviz object
#' @rdname NCdisplay
##' @export
setMethod("NCdisplay", "NCviz", NCdisplay)

toFileName <- function(string) {
    string = gsub(",? ", "_", string)
    string = gsub("\\(|\\)", "", string)
    string = gsub(",", "_", string)
    return(string)
}

saveInFilesF <- function(obj, path="./", suffix="") {
    # Save data
    for (method in names(obj@nc_data)) {
        print(paste("Saving", method))
        ff = file(paste0(path, toFileName(obj@cell_type), "_", method, ifelse(suffix=="", "", "_"), suffix, ".tsv"), "w")
        writeLines(paste0(c("GENE", colnames(obj@nc_data[[method]])), collapse="\t"), ff)
        #write.table(obj@nc_data[[method]], ff, sep="\t", col.names=FALSE)
        for (gene in rownames(obj@nc_data[[method]])) {
            writeLines(paste0(c(gene, as.character(obj@nc_data[[method]][gene,])), collapse="\t"), ff)
         }
        close(ff)
    }
    # Save annotations
    ff = file(paste0(path, toFileName(obj@cell_type), "_Annotations", ifelse(suffix=="", "", "_"), suffix, ".tsv"), "w")
    writeLines(paste0(c("NAME", colnames(obj@annotations)), collapse="\t"), ff)
    for (spl in rownames(obj@annotations)) {
        writeLines(paste0(c(spl, as.character(obj@annotations[spl,])), collapse="\t"), ff)
    }
    close(ff)
}
setGeneric("saveInFiles", saveInFilesF)
#' Save the data in several files.
#'
#' Save the data from a NCviz object in several files, one .tsv file per profiling method. Each file can then be imported into NaviCell.
#' Also produces a file with annotations of all samples, also importable in NaviCell.
#'
#' @param obj NCviz object
#' @param path Folder where the files must be save, can be used to append a prefix to the filename
#' @param suffix Suffix appended to the filename
#' @return Produces .tsv files, no R output.
#' @author Mathurin Dorel \email{mathurin.dorel@@curie.fr}
#' @export
#' @rdname saveInFiles
setMethod("saveInFiles", "NCviz", saveInFilesF)


saveDataF <- function(obj, path="./", suffix="") {
    ff = file(paste0(path, toFileName(obj@cell_type), ifelse(suffix=="", "", "_"), suffix, ".txt"), "w")
    for (method in names(obj@nc_data)) {
        # Save data
        print(paste("Saving", method))
        writeLines(paste0("M ", method), ff)
        write.table(obj@nc_data[[method]], ff, sep="\t")
    }
    # Save annotations
    writeLines(paste0("ANNOTATIONS"), ff)
    write.table(obj@annotations, ff, sep="\t")
    close(ff)
}
setGeneric("saveData", saveDataF)
#' Save the data in one files.
#'
#' Save the data in a .txt files. The file cannot be directly exported to NaviCell but can be imported with the RncMapping package. For easily commucating data.
#'
#' @name saveData
#' @param obj NCviz object
#' @param path Folder where the files must be save, can be used to append a prefix to the filename
#' @param suffix Suffix appended to the filename
#' @return Produces a .txt file, no R output.
#' @author Mathurin Dorel \email{mathurin.dorel@@curie.fr}
#' @export
#' @rdname saveData
setMethod("saveData", "NCviz", saveDataF)

