#my librarier
library(ENMeval)
library(raster)
library(MASS)
library(rgdal)

#load environmental variables in.asc extension 
bio1 <- raster("bio1_big.asc")
bio2 <- raster("bio2_big.asc")
bio3 <- raster("bio3_big.asc")
bio4 <- raster("bio4_big.asc")
bio5 <- raster("bio5_big.asc")
bio6 <- raster("bio6_big.asc")
bio7 <- raster("bio7_big.asc")
bio8 <- raster("bio8_big.asc")
bio9 <- raster("bio9_big.asc")
bio10 <- raster("bio10_big.asc")
bio11 <- raster("bio11_big.asc")
bio12 <- raster("bio12_big.asc")
bio13 <- raster("bio13_big.asc")
bio14 <- raster("bio14_big.asc")
bio15 <- raster("bio15_big.asc")
bio16 <- raster("bio16_big.asc")
bio17 <- raster("bio17_big.asc")
bio18 <- raster("bio18_big.asc")
bio19 <- raster("bio19_big.asc")
dem <- raster("el_big.asc")
bio_cat1 <- raster("land_big.asc")

#tell my code that this is a categorical variable
bio_cat2 <- as.factor(bio_cat1)

#stack everything
env <- stack(bio1, bio2, bio3, bio4, bio5, bio6, bio7, bio8,bio9,bio10,bio11,bio12,bio13,bio14,bio15,bio16,bio17,bio18,bio19,dem, bio_cat1)



#load observations
occs <- read.csv("san_tx_fin.csv")[,-1]

#load file with testiing observations
occs.testing <- read.csv("test_san.csv")[,-1]

#plot them
occur.ras <- rasterize(occs, env, 1)
plot(occur.ras)

presences <- which(values(occur.ras) == 1)
pres.locs <- coordinates(occur.ras)[presences, ]

#create bias file
dens <- kde2d(pres.locs[,1], pres.locs[,2], n = c(nrow(occur.ras), ncol(occur.ras)), lims = c(extent(env)[1], extent(env)[2], extent(env)[3], extent(env)[4]))
dens.ras <- raster(dens, env)
dens.ras2 <- resample(dens.ras, env)
plot(dens.ras2)
dens.ras

writeRaster(dens.ras2, "biastexas_san.asc", overwrite = TRUE)
length(which(!is.na(values(subset(env, 1)))))

#create background
bg <- xyFromCell(dens.ras2, sample(which(!is.na(values(subset(env, 1)))), 5000, prob=values(dens.ras2)[!is.na(values(subset(env, 1)))]))
colnames(bg) <- colnames(occs)

#line with categorical variables
#enmeval_results <- ENMevaluate(occ, env, bg, tune.args = list(fc = c("L","LQ","H", "LQH", "LQHP", "LQHPT"), rm = 1:5), partitions = "randomkfold", partition.settings = list(kfolds = 10), algorithm = "maxnet",categoricals="bio_cat2")

#line without categorical variables
enmeval_results <- ENMevaluate(occs, env, bg, tune.args = list(fc = c("L","LQ","H", "LQH", "LQHP", "LQHPT"), rm = 1:5), partitions = "randomkfold", partition.settings = list(kfolds = 10), algorithm = "maxnet")

#line with test file (% data used for evaluating fitting)
#enmeval_results <- ENMevaluate(occs, occs.testing=testfile,env, bg, tune.args = list(fc = c("L","LQ","H", "LQH", "LQHP", "LQHPT"), rm = 1:5), partitions = "testing", partition.settings = list(kfolds = 10), algorithm = "maxnet")

#save results
enmeval_results@results
write.csv(enmeval_results@results, "enmeval_results_san.csv")

