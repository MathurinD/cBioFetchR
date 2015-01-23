#' @import methods
# @import RNaviCell

check_dataset <- function(object) {
    errors = c() 

    no_cbio = (length(object@cbio_data) == 0)
    no_nc_data = (length(object@nc_data) == 0)
    if (!no_cbio && !no_nc_data) {
        errors = c(errors, "Only one dataset must be provided")
    }
    if (no_cbio && no_nc_data) {
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
            cbio_data="list", # Data indexed by genes
            nc_data="list", # Data indexed by experiment type
            verbose="logical"
            ),
     #prototype(nc_url="", cell_type="", cbio_data=list(), nc_data=list()),
     validity=check_dataset
     )

#' Constructor for the class NCviz. NCviz objects provide functions to visualize profiling data on a NaviCell map
#' @param nc_url URL of the NaviCell server
#' @param cell_type Name of the dataset to pass to the NaviCell server
#' @param cbio_data Data loaded with cBioStudy (incompatible with nc_data)
#' @param nc_data Data under the NaviCell format, a list indexed by experiment type containing data arranged in a dataframe with genes HUGO indentifiers as rownames and samples id as colnames (incompatible with cbio_data)
#' @export
#' @seealso cBioStudy
#' @author Mathurin Dorel \email{mathurin.dorel@@curie.fr}
#' @rdname NCviz-class
NCviz <- function(nc_url="", cell_type="", cbio_data=list(), nc_data=list(), verbose=TRUE) {
    return(new( "NCviz", nc_url=nc_url, cell_type=cell_type, cbio_data=cbio_data, nc_data=nc_data, verbose=verbose ))
}

# TODO Determine if some extra data are recquired
setMethod("initialize",
          "NCviz",
            function(.Object, nc_url, cell_type, cbio_data, nc_data, verbose=TRUE) {
                print(paste("Creation of an", "NCviz", "object"))
                .Object@nc_url = nc_url
                .Object@cell_type = cell_type
                .Object@cbio_data = cbio_data
                .Object@nc_data = nc_data
                .Object@verbose = verbose
                validObject(.Object) # Validate object

                # Copy data in the other format
                if (length(.Object@nc_data) == 0) {
                    print("Converting form c-Bioportal format to NaviCell format")
                    for (experiment in colnames(.Object@cbio_data[[1]])) {
                        print(paste("Inserting experiment", experiment))
                        .Object@nc_data[[experiment]] = list()
                        for (gene in names(.Object@cbio_data)) {
                            .Object@nc_data[[experiment]][[gene]] = .Object@cbio_data[[gene]][[experiment]]
                        }
                        .Object@nc_data[[experiment]] = as.data.frame(t( as.data.frame(.Object@nc_data[[experiment]]) ))
                        colnames(.Object@nc_data[[experiment]]) = rownames(.Object@cbio_data[[1]])
                    }
                } else {
                    .Object@cbio_data[[length(rownames(.Object@nc_data[[1]]))]] = 0
                    print("No conversion from navicell format to cbioportal format") # as it is useless
                    for (gene in rownames(.Object@nc_data[[1]])) {
                        # TODO ? copy in not navicell format
                    }
                }

                return(.Object)
            }
          )


#setGeneric(addData) # Add data tables
#setGeneric()

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

saveInFilesF <- function(obj, path="./", suffix="") {
    for (method in names(obj@nc_data)) {
        print(paste("Saving", method))
        ff = file(paste0(path, gsub(" ", "_", obj@cell_type), "_", method, suffix, ".tsv"), "w")
        writeLines(c("GENE", colnames(obj@nc_data[[method]])), sep="\t")
        for (gene in rownames(obj@nc_data[[method]])) {
            writeLines(c(gene, obj@nc_data[[method]][gene,]), sep="\t")
        }
        close(ff)
    }
}
setGeneric("saveInFiles", saveInFilesF)
#' Save the data in several files.
#'
#' Save the data in several files, one .tsv file per profiling method. Each file can then be imported into NaviCell.
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
    ff = file(paste0(path, gsub(" ", "_", obj@cell_type), "_", suffix, ".tsv"), "w")
    for (method in names(obj@nc_data)) {
        print(paste("Saving", method))
        writeLines(paste0("M ", method))
        writeLines(c("GENE", colnames(obj@nc_data[[method]])), sep="\t")
        for (gene in rownames(obj@nc_data[[method]])) {
            writeLines(c(gene, obj@nc_data[[method]][gene,]), sep="\t")
        }
    }
    close(ff)
}
setGeneric("saveData", saveDataF)
#' Save the data in one files.
#'
#' Save the data in a .txt files. The file cannot be directly exported to NaviCell but can be imported with the RncMapping package.
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

