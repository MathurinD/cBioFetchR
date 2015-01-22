# S4 class for genes data
# Contains statistics and functions to communicate with a NaviCell API

library(method)

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

NCclass = "NCviz"
setClass(NCclass,
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

NCviz <- function(nc_url="", cell_type="", cbio_data=list(), nc_data=list(), verbose=TRUE) {
    return(new( NCclass, nc_url=nc_url, cell_type=cell_type, cbio_data=cbio_data, nc_data=nc_data, verbose=verbose ))
}

# TODO Determine if some data are recquired
setMethod("initialize",
          NCclass,
            function(.Object, nc_url, cell_type, cbio_data, nc_data, verbose=TRUE) {
                print(paste("Creation of an", NCclass, "object"))
                .Object@nc_url = nc_url
                .Object@cell_type = cell_type
                .Object@cbio_data = cbio_data
                .Object@nc_data = nc_data
                .Object@verbose = verbose
                validObject(.Object) # Validate object

                # Copy data in the other format
                if (length(.Object@nc_data) == 0) {
                    print("Converting form c-Bioportal format to NaviCell format")
                    .Object@nc_data[[length(colnames(.Object@cbio_data[[1]]))]] = 0
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
                    print("No conversion from navicell format to cbioportal format, as it is useless")
                    for (gene in rownames(.Object@nc_data[[1]])) {
                        # TODO ? copy in not navicell format
                    }
                }

                return(.Object)
            }
          )


setGeneric(addData) # Add data tables
setGeneric()

display_function="NCdisplay" # Display the data in a good manner in NaviCell
setGeneric(display_function, function(obj){return(standardGeneric(display_function))})
setMethod(display_function, NCclass,
            function(obj) {
            # TODO display in NaviCell
            }
         )

