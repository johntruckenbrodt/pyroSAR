# this script is intended for cropping multiple images to their common extent
# John Truckenbrodt 2015

suppressMessages(library(raster))

raster::rasterOptions(format="ENVI",overwrite=T,datatype="FLT4S",setfileext=F)

rasterlist=lapply(commandArgs(T),function(x)raster::raster(x))

ext=raster::extent(rasterlist[[1]])
for(i in 2:length(rasterlist)){
  ext=raster::intersect(ext,raster::extent(rasterlist[[i]]))
}

for(item in rasterlist){
  outname=paste(tools::file_path_sans_ext(item@file@name),"sub1",sep="_")
  if(!file.exists(outname)){
    raster::writeRaster(raster::crop(item,ext),filename=outname,bandorder="BSQ",NAflag=0)
  }  
}




#/pvdata2/john/Landsat/patrick/LT51130262010246.envi /pvdata2/john/Landsat/patrick/LT51130272010246.envi


test1=raster::stack("D:/Landsat/test1")
test2=raster::stack("D:/Landsat/test2")

inter=raster::intersect(raster::extent(test1),raster::extent(test2))



mat1=raster::as.matrix(test1@layers[[1]])
mat2=raster::as.matrix(test2@layers[[1]])

if(length(mat1[is.na(mat1)&!is.na(mat2)])==0){
  stop("no intersect (yet)")
}




test=raster::stack("D:/Landsat/test")

test=raster::stack(test1)@layers
test3=raster::as.matrix(test[[1]])

