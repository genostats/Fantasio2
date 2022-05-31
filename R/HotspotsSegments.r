##################################################################################
#This is the class that encapsulate the list of segments created by using the    #
#hotspots in the genome                                                          #
#                                                                                #
##################################################################################
#' Class HostspotsSegments
#'
#' Class \code{HostspotsSegments} This is the class that encapsulate the list of segments created by using the hotspots in the genome.
#' @exportClass HostspotsSegments
#' 
setClass("HostspotsSegments", representation = (
  HostspotsSegments = 'list'
))

#' Show method of HostspotsSegments.
#'
#' @param object the HostspotsSegments
setMethod('show', signature("HostspotsSegments"), 
          function(object){
            cat('A  list of ', length(object), 'chromosoms segments created using hotspots in the genome\n')
          })