#######################################################################################################################
#This is the class used to create an object which will contains every dataframe and list created when creating submaps#
#######################################################################################################################

setClassUnion("listOrNULL",members=c("list", "NULL"))
setClassUnion("dataframeOrNULL",members=c("data.frame", "NULL"))
setClassUnion("matrixOrNULL",members=c("matrix", "NULL"))
setClassUnion("characterOrNULL",members = c("character", "NULL"))
setClassUnion("doubleOrNULL",members = c("numeric", "NULL"))
setClassUnion("logicalOrNULL",members = c("numeric", "NULL"))

#' Class atlas
#'
#' Class \code{atlas} This is the class used to create an object which will contains every dataframe and list created when creating submaps.  
#'
#' @rdname atlas-class
#' @exportClass atlas
#' @slot segments_list the list of segments 
#' @slot submaps_list a list of submaps
#' @slot likelihood_summary a dataframe with both likelihood0 and likelihood1 over the submaps
#' @slot estimation_summary a dataframe with both a and f estimation over the submaps
#' @slot submap_summary a dataframe with summary statitistics about the submaps
#' @slot HBD_recap a dataframe with for one individual and for one marker a mean computation of all the HBD probabilities computed, on every individuals.
#' @slot FLOD_recap a dataframe with for one individual and for one marker a mean computation of all the FLOD scores computed, on every individuals.
#' @slot HBD_segments a list of dataframe with the HBD segments, on all individuals.
#' @slot bedmatrix  a bed.matrix object (refer to gaston package)
#' @slot bySegments a boolean indicating wheter the creation of summary statistics was made by segments (see documentation of Fantasio function)
#' @slot unit   the unit of the markers (cM or Bp).
#' @slot gap   the value of the gap used to pick marker when doing submaps by snps. (see function Fantasio for more infos)

setClass("atlas", representation(
        bedmatrix            = 'bed.matrix', 
        seeds                = 'matrixOrNULL',
        epsilon              = 'numeric',
        segments_list        = 'listOrNULL',
#       submaps_list         = 'listOrNULL', 
#       likelihood_summary   = 'listOrNULL',
        estimations          = 'listOrNULL', 
        submap_summary       = 'dataframeOrNULL',
        recap                = "characterOrNULL", 
        q                    = "doubleOrNULL",
        HBD_recap            = 'matrixOrNULL',
        FLOD_recap           = 'matrixOrNULL',  
        HBD_segments          = 'listOrNULL',
        unit                 = "characterOrNULL", 
        gap                  = "doubleOrNULL"
))

#' Show method of class atlas.
#'
#' @param object an atlas object
setMethod('show', signature("atlas"), 
  function(object){
       cat('An atlas of', ncol(object@seeds), 'submaps\n ')
  })


#' Constructor of class atlas.
#'
#' @param .Object the object 
#' @param bedmatrix a bed.matrix object
#' @param seeds
#' @param epsilon
#' @param segments.list a list of segments object
#' @param estimations
#' @param submap.summary

setMethod('initialize', signature='atlas', definition=function(.Object, bedmatrix, seeds, epsilon, segments.list, estimations, submap.summary)
{
  .Object@bedmatrix       <- bedmatrix
  .Object@seeds           <- seeds
  .Object@epsilon         <- epsilon
  .Object@segments_list   <- segments.list
  .Object@estimations     <- estimations
  .Object@submap_summary  <- submap.summary
  .Object
})


