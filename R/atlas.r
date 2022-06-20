#######################################################################################################################
#This is the class used to create an object which will contains every dataframe and list created when creating submaps#
#######################################################################################################################

setClassUnion("listOrNULL",members=c("list", "NULL"))
setClassUnion("dataframeOrNULL",members=c("data.frame", "NULL"))
setClassUnion("matrixOrNULL",members=c("matrix", "NULL"))
setClassUnion("characterOrNULL",members = c("character", "NULL"))
setClassUnion("doubleOrNULL",members = c("numeric", "NULL"))

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
#' @slot HBDsegments a list of dataframe with the HBDsegments, on all individuals.
#' @slot HFLOD a dataframe with the value of HFLOD scores for every markers through all submaps.
#' @slot bedmatrix  a bed.matrix object (refer to gaston package)
#' @slot bySegments a boolean indicating wheter the creation of summary statistics was made by segments (see documentation of Fantasio function)
#' @slot unit   the unit of the markers (cM or Bp).
#' @slot gap   the value of the gap used to pick marker when doing submaps by snps. (see function Fantasio for more infos)
#' @slot logisticRegression the results of logistic regression

setClass("atlas", representation(
        bedmatrix            = 'bed.matrix', 
        seeds                = 'matrixOrNULL',
        submaps              = 'submapsOrNULL',
#       segments_list        = 'listOrNULL',
#       submaps_list         = 'listOrNULL', 
#       likelihood_summary   = 'listOrNULL',
#       estimation_summary   = 'listOrNULL', 
        submap_summary       = 'dataframeOrNULL',
        HBD_recap            = 'matrixOrNULL',
        FLOD_recap           = 'matrixOrNULL',  
        HBDsegments          = 'listOrNULL',
        HFLOD                = 'dataframeOrNULL',
        bySegments           = "logical", 
        unit                 = "characterOrNULL", 
        gap                  = "doubleOrNULL",
        logisticRegression   = "listOrNULL"
))

#' Constructor method of atlas.
#'
#' @param .Object the object type
#' @param submaps a list of submaps
#' @param bedmatrix a bed.matrix object
#' @param segments_list a list of segments object
#' @param bySegments a boolean
#' @param unit the unit in which are markers (cM or Bp)
#' @param gap the gap used to create segments in the By Distance method
#' @param logReg list of logistic regression results on adjusted and unadjusted data

setMethod('initialize', signature='atlas', definition=function(.Object, submaps, bedmatrix, segments_list, bySegments, unit=NULL, gap=NULL)
{
  .Object@submaps_list        <- submaps
  .Object@bedmatrix           <- bedmatrix
  .Object@segments_list       <- segments_list
  .Object@bySegments          <- bySegments
  .Object@unit                <- unit
  .Object@gap                 <- gap
  .Object
})


#' Show method of atlas.
#'
#' @param object an atlas object
setMethod('show', signature("atlas"), 
  function(object){
       cat('An atlas of', length(object@submaps_list), 'submaps\n ')
  })



