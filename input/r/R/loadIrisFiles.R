#' Load iris files 
#'
#' Load all iris files
#' @param folder.path Path to folder
#' @param keep.folder.name Boolean, defaults to FALSE
#' @export
loadIrisFiles = function(folder.path, keep.folder.name = F){

  loadOneIrisFile = function(iris.filename){
    iris.file.data = data.table::fread(iris.filename)
    
    iris.file.data[,filename:=iris.filename]
    
    return(iris.file.data)
  }
  
  all.iris.filenames = paste(folder.path, 
    list.files(path = folder.path, pattern = "\\.iris$", all.files = T),
    sep="/")
  
  all.iris.files.data =
    data.table::rbindlist(
      lapply(all.iris.filenames, loadOneIrisFile)
    )
  
  
  if(!keep.folder.name){
    all.iris.files.data[,filename:=basename(filename)]
  }
  
  return(all.iris.files.data)
}
