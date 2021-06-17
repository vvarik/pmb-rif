#' Add some variables
#'
#' @export
addIrisVars = function (dat) {
  #Set the right plate number, replicate number and colony
  dat[,media:=my.strsplit(filename, "-", 1)]
  dat[,antibiotic:=my.strsplit(filename, "-", 2)]
  dat[,condition:=paste(media, antibiotic, sep="-")]
  
  dat[,plate.id:=
                   filename %>% 
                   my.strsplit(".JPG.iris", 1) %>% 
                   my.strsplit("-", 3) %>% 
                   as.numeric()]
  
  #get library plate number and replicate numbers
  #modulo divisions work on zero-based numbers, so removing one to start with 
  #and adding it back later since both library plate and replicate are 1-based numbers
  dat[,plate.number:=(plate.id-1)%/%5 + 1]
  dat[,biorep:=(plate.id-1)%%5 + 1]
  
  dat[,colony:=10000*plate.number + 100*row + column]
  dat[,plate.id:=NULL]
  
}
