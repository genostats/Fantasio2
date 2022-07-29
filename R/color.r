##################################################################################
#This function allows the algorithms for plotting HBD to chose a color           #
#                                                                                #
#!!! val : the phenotype of the individual                                       #
#                                                                                #
#*** return the value of the color                                               #
##################################################################################

color  <- function (val) {
  if(is.na(val))
    "black"
  else if (val==1) 
    "skyblue3"
  else if (val==2)
    "tomato"
  else
    "grey"
}