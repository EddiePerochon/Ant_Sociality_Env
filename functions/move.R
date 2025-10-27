#### Packages and internal function
pkg <- c("vegan")
lapply(pkg, function(x) require(x, quietly = TRUE, character.only = TRUE))

#--> vegan: for procrustes analysis

### Associate a colour of the colour cloud to each of the 657 colours of the colours() vector
quiet <- TRUE
if(!quiet){
	library(rgl)	
	library(chroma)	#--> chroma: calculate CIE2000 distance between colours (instead of spaceXYZ::DeltaE used previously)
	nam <- colours()
	ind <- NULL
	ddd <- NULL
	for(n in nam){
		d <- chroma::deltaE(n, hex_LAB)
		i <- which(d == min(d))
		ind <- c(ind, i)
		ddd <- c(ddd, d[i])
	print(match(n, nam))
	}
	range(ddd)
	col_LAB <- cbind.data.frame(nam, ind, ddd)
	save(col_LAB, file = "col_LAB.RData")
	rgl::plot3d(coords_LAB, type = "n")
	rgl::spheres3d(coords_LAB[col_LAB$ind,], col = hex_LAB[col_LAB$ind], r = 35)
	rgl::text3d(coords_LAB[col_LAB$ind,], text = col_LAB$nam, col = hex_LAB[col_LAB$ind])
}


####	####	####	####
#--> Function to move a cloud according two (or more points) points

"move" <- function(coords, coords_LAB, hex_LAB, col_LAB, from, to){
#--> look at select3d to select point interactively

	#--> coords: cloud of data points (after initiation)
	#--> coords_LAB: cloud of color points in euclidean space (after initiation) 
	#--> hex_LAB: colors (hex code) associated to coords_LAB points 
	#--> col_LAB: correspondance between 657 names of colours and colours of coords_LAB
	#--> from: vector of two (or more) points defining the vector of origin
	#--> to: vector of two (or more) points defining the vector of arrival

	### Asymetric procruste transformation to move (translation + rotation + scaling) Y into X
	ind <- col_LAB[match(to, col_LAB$nam),"ind"]
	X <- coords_LAB[ind,]
	Y <- coords[from,]
	pro <- vegan::procrustes(X,Y)
	mov <- stats::predict(pro, coords)

return(mov)
}


####	####	####	####
#--> description file
quiet <- TRUE
if(!quiet){
	library(terra)
	library(FNN)
	library(colorspace)
	source("C:/Users/Seb/Documents/Ma recherche/Papers/2024_Methods in Ecology and Evolution_project/functions/initiate.R")
	source("C:/Users/Seb/Documents/Ma recherche/Papers/2024_Methods in Ecology and Evolution_project/functions/move.R")

	dir <- "C:/Users/Seb/Documents/Ma recherche/Papers/2024_Methods in Ecology and Evolution_project/RData/"
	path <- paste0(dir, "LAB.RData")#
	load(path)
	coords_LAB <- as.matrix(coords_LAB)

	dir <- "C:/Users/Seb/Documents/Ma recherche/Papers/2024_Methods in Ecology and Evolution_project/RData/"
	path <- paste0(dir, "col_LAB.RData")#
	load(path)

	#### Test the function on beck categories
	### Load data to test the function
	path <-  paste0("C:/Users/Seb/Documents/Ma recherche/Data/Climatic/modified/beck_10km.tif")
	beck_10km <- rast(path)
	path <- paste0("C:/Users/Seb/Documents/Ma recherche/Data/Climatic/modified/beck_df.txt")
	beck_df <- read.table(path, h = T)
	cat <- beck_df$c # the climate categories
	ncat <- length(cat) # number of categories
	pal_beck <- beck_df$col_cc # the beck colours 
	names(pal_beck) <- cat

	path <-  paste0("C:/Users/Seb/Documents/Ma recherche/Data/Climatic/modified/worldclim_pca_10km.tif")
	worldclim_pca_10km <- rast(path)

	pca_li <- terra::values(worldclim_pca_10km)
	fac <- terra::values(beck_10km)[,1]
	index <- which(is.nan(fac)) 
	li <- as.data.frame(pca_li[-index,1:3]) # all the worldclim pca first 3 axis 
	fac <- fac[-index] 
	fac <- as.factor(fac) 
	levels(fac) <- cat # categories of koppen-geiger
	local <- terra::split(li, fac)  
	local <- lapply(local, function(x) apply(x, 2, mean))  
	coords <- as.data.frame(t(as.data.frame(local))) 
	
	### Run the initiate function
	ini <- initiate(coords, coords_LAB, hex_LAB, ds = 0.1)
	plot(beck_10km, col = ini$pal)

	### Move the data cloud between green & cyan	
		# Af --> green
		# EF --> blue
	from <- c("Af", "EF")
	to <- c("green", "blue")
	mov <- move(ini$coords, ini$coords_LAB, hex_LAB, col_LAB, from = from, to = to)
	cat <- row.names(ini$coords)
	rgl::text3d(ini$coords, texts = cat)
	ind <- match(from, cat)
	rgl::spheres3d(ini$coords[ind,], col = to, r = 0.25)
	nam <- col_LAB$nam
	ind <- match(to, nam)
	rgl::spheres3d(ini$coords_LAB[col_LAB[ind,"ind"],], col = to, r = 0.15)

	rgl::spheres3d(mov, r = 0.1)
	#--> rotation +/- pi with two points

	### Move the data cloud with deformation (3 points)
		# Af --> green
		# EF --> cyan
		# BWh --> red
	from <- c("Af", "BWh", "EF")
	to <- c("green", "red", "blue")
	mov <- move(ini$coords, ini$coords_LAB, hex_LAB, col_LAB, from = from, to = to)
	cat <- row.names(ini$coords)
	rgl::text3d(ini$coords, texts = cat)
	ind <- match(from, cat)
	rgl::spheres3d(ini$coords[ind,], col = to, r = 0.25)
	nam <- col_LAB$nam
	ind <- match(to, nam)
	rgl::spheres3d(ini$coords_LAB[col_LAB[ind,"ind"],], col = to, r = 0.15)

	rgl::spheres3d(mov, r = 0.1)

	knn <- FNN::get.knnx(ini$coords_LAB, mov, k = 1)
	w <- knn$nn.index[,1]
	pal <- hex_LAB[w]
	names(pal) <- cat	
	plot(beck_10km, col = pal)
	x11()
	plot(beck_10km, col = pal_beck)
}

