#!/bin/r
# function: the whole process we do analysis with gf gene offset

require(raster)
require(geosphere)
require(gdm)
require(foreach)
require(parallel)
require(doParallel)
require(gradientForest)
require(fields)

# set the working directory
setwd("/home/chichi/data/china/china3/future/gene_offset")

##############
#Read in population locations & climate data
##############

pops <- read.csv("/home/chichi/data/china/china3/future/gene_offset/pop.csv")
predNames <- c("Band_1","Band_2","Band_3","Band_4","Band_5","Band_6","Band_7","Band_8","Band_9","Band_10","Band_11","Band_12","Band_13","Band_14","Band_15","Band_16","Band_17","Band_18","Band_19")
presClim <- stack("/home/chichi/data/china/china3/future/gene_offset/bios_2.tif")
presClim <- presClim[[predNames]]
pred <- data.frame(pop=pops$code,long=pops$long, lat=pops$lat, extract(presClim, y=pops[,c("long","lat")]),
                   stringsAsFactors=FALSE)

######################
#GF model
######################
#read in maf data
snps <- read.csv("/home/chichi/data/china/china3/future/gene_offset/maf.csv", stringsAsFactors=FALSE)
snpsAll <- colnames(snps)[-c(1:3)]
snps <- merge(pred, snps, by.x="pop", by.y="code", all.x=TRUE)

gfMod <- gradientForest(data=snps, predictor.vars=predNames,response.vars=snpsAll,
                        ntree=500, 
                        maxLevel=log2(0.368*nrow(snps)/2), trace=T, 
                        corr.threshold=0.70)

shp <- shapefile("/home/chichi/data/china/china3/future/gene_offset/TlatifoliaCurrent.shp")

#load future climate data
futClims <- stack("/home/chichi/data/china/china3/future/gene_offset/wc2.1_2.5m_bioc_MIROC6_ssp370_2081-2100_2.tif") #stack future climate layers
#predNames <- c("wc2.1_2.5m_bioc_MIROC6_ssp126_2081.2100_1","wc2.1_2.5m_bioc_MIROC6_ssp126_2081.2100_2","wc2.1_2.5m_bioc_MIROC6_ssp126_2081.2100_3","wc2.1_2.5m_bioc_MIROC6_ssp126_2081.2100_4")
futClims <- futClims[[predNames]]
futClimDat <- as.data.frame(futClims, xy=TRUE, na.rm=TRUE)

futClimDatGF <- data.frame(futClimDat[,c("x","y")], predict(gfMod,futClimDat[,predNames])) 

popDatGF <- na.omit(as.data.frame(mask(presClim, shp), xy=TRUE))
#after this step the x, y are different from the presClim, and presClim is the same as the futClimDat in x, y

popDatGF <- data.frame(popDatGF[,c("x","y")], predict(gfMod, popDatGF[,predNames]))
popDatGF <- split(popDatGF, seq(nrow(popDatGF)))

###############
#Forward offset calculation
##############
cl <- makeCluster(8)
registerDoParallel(cl)
print("start parallel analysis")
forwardOffsetGF <- foreach(i = 1:length(popDatGF), .packages=c("fields","gdm","geosphere")) %dopar%{
  
  #get the focal population
  onePopGF <- popDatGF[[i]]
  
  #get destination populations and add gf distance
  combinedDatGF <- futClimDatGF[,c("x","y")]
  combinedDatGF["gfOffset"] <- c(rdist(onePopGF[,predNames], futClimDatGF[,predNames]))
  ###########
  # set the distance to 3 degree add by chichi
  # 3 degree is about 333 km
  combinedDatGF <- combinedDatGF[which(combinedDatGF$x > onePopGF$x-3 & combinedDatGF$x < onePopGF$x+3 & combinedDatGF$y > onePopGF$y-3 & combinedDatGF$y < onePopGF$y+3),]
  #############
  
  ##Get metrics for the focal population
  #coordinate of focal population
  coordGF <- onePopGF[,c("x","y")]
  # as the onePopGF X,Y in not the same as the combinedDatGF X,Y.
  
  #choose the pixels with the minimum fst
  minCoordsGF <- combinedDatGF[which(combinedDatGF$gfOffset == min(combinedDatGF$gfOffset)),]
  
  #calculate the distance to the sites with minimum fst, and selct the one with the shortest distance
  minCoordsGF["dists"] <- distGeo(p1=coordGF, p2=minCoordsGF[,1:2])
  minCoordsGF <- minCoordsGF[which(minCoordsGF$dist == min(minCoordsGF$dists)),]
  
  #if multiple sites have the same fst, and same distance, one is randomly chosen
  minCoordsGF <- minCoordsGF[sample(1:nrow(minCoordsGF),1),]
  
  #get local offset
  # due to floating point precision issues, we can't find the exact match
  #the following is the solution
  tolerance <- 0.00001
  offsetGF <- combinedDatGF[which(abs(combinedDatGF$x - coordGF$x) < tolerance & abs(combinedDatGF$y - coordGF$y) < tolerance),"gfOffset"]
  #offsetGF <- combinedDatGF[which(combinedDatGF$x == coordGF$x & combinedDatGF$y ==coordGF$y),"gfOffset"]
  
  #get the minimum predicted fst - forward offset in this case
  minValGF <- minCoordsGF$gfOffset
  
  #get distance and coordinates of site that minimizes fst
  toGoGF <- minCoordsGF$dists
  minPtGF <- minCoordsGF[,c("x","y")]
  
  #get bearing to the site that minimizes fst
  bearGF <- bearing(coordGF, minPtGF)
  
  #write out
  outGF <- c(x1=coordGF[[1]], y1=coordGF[[2]], local=offsetGF, forwardOffset=minValGF, predDist=toGoGF, bearing=bearGF,x2=minPtGF[[1]],y2=minPtGF[[2]])
  
}
stopCluster(cl)


#in this resultant dataframe the columns are:
#x1/y1: focal coordinates
#local: local offset
#forwardOffset: forward offset
#predDist: distance to site of forward offset
#bearing: bearing to site of forward offset
#x2/y2: coordinate of site of forward offset
forwardOffsetGF <- do.call(rbind, forwardOffsetGF)

write.csv(forwardOffsetGF,paste0("./forwardOffsetGF.csv"), row.names=FALSE)


###############
#Reverse offset calculation
##############
#Getting all coordinates in the range in current climate and transform using gf
popDatGF <- na.omit(as.data.frame(mask(presClim, shp), xy=TRUE))
popDatGF <- data.frame(popDatGF[,c("x","y")], predict(gfMod, popDatGF[,predNames]))

#Gets climate data from the range in future climate and transform using gf
futClimMask <- mask(futClims, mask=shp)
futClimDat <- as.data.frame(futClimMask, xy=TRUE, na.rm=TRUE)
futClimDatGF <- data.frame(futClimDat[,c("x","y")],predict(gfMod,futClimDat[,predNames]))


###############
#Reverse offset calculation
##############
cl <- makeCluster(32)
registerDoParallel(cl)
reverseOffsetGF <- foreach(i = 1:nrow(futClimDatGF), .packages=c("fields","gdm","geosphere")) %dopar%{
  
  #get the focal population in future climate
  onePopGF <- futClimDatGF[i,]
  
  #make prediction between focal population and current climate
  combinedDatGF <- popDatGF[,c("x","y")]
  combinedDatGF["gfOffset"] <- c(rdist(onePopGF[,predNames], popDatGF[,predNames]))

  ###########
  # set the distance to 3 degree add by chichi
  # 3 degree is about 333 km
  combinedDatGF <- combinedDatGF[which(combinedDatGF$x > onePopGF$x-3 & combinedDatGF$x < onePopGF$x+3 & combinedDatGF$y > onePopGF$y-3 & combinedDatGF$y < onePopGF$y+3),]
  #############
  
  ##Get metrics for the focal population
  #coordinate of focal population
  coordGF <- onePopGF[,c("x","y")]
  
  #choose the pixels with the minimum offset
  minCoordsGF <- combinedDatGF[which(combinedDatGF$gfOffset == min(combinedDatGF$gfOffset)),]
  
  #calculate the distance to the sites with minimum fst, and selct the one with the shortest distance
  minCoordsGF["dists"] <- distGeo(p1=coordGF, p2=minCoordsGF[,1:2])
  minCoordsGF <- minCoordsGF[which(minCoordsGF$dists == min(minCoordsGF$dists)),]
  
  #if multiple sites have the same fst, and same distance, one is randomly chosen
  minCoordsGF <- minCoordsGF[sample(1:nrow(minCoordsGF),1),]
  
  #get local offset
  # due to floating point precision issues, we can't find the exact match
  #the following is the solution
  tolerance <- 0.00001
  offsetGF <- combinedDatGF[which(abs(combinedDatGF$x - coordGF$x) < tolerance & abs(combinedDatGF$y - coordGF$y) < tolerance),"gfOffset"]
  # offsetGF <- combinedDatGF[which(combinedDatGF$x == coordGF$x & combinedDatGF$y == coordGF$y),"gfOffset"]
  
  #get the minimum predicted offset - reverse offset in this case
  minValGF <- minCoordsGF$gfOffset
  
  #get distance and coordinates of site that minimizes fst
  toGoGF <- minCoordsGF$dists
  minPtGF <- minCoordsGF[,c("x","y")]
  
  #get bearing to the site that minimizes fst
  bearGF <- bearing(coordGF, minPtGF)
  
  #write out
  outGF <- c(x1=coordGF[[1]], y1=coordGF[[2]],local=offsetGF, reverseOffset=minValGF, predDist=toGoGF, bearing=bearGF, x2=minPtGF[[1]],y2=minPtGF[[2]])
  
}

stopCluster(cl)


#in this resultant dataframe the columns are:
#x1/y1: focal coordinates
#local: local offset
#reverseOffset: reverse offset
#predDist: distance to site of reverse offset
#bearing: bearing to site of reverse offset
#x2/y2: coordinate of site of reverse offset
reverseOffsetGF <- do.call(rbind, reverseOffsetGF)

write.csv(reverseOffsetGF,paste0("./reverseOffsetGF.csv"), row.names=FALSE)