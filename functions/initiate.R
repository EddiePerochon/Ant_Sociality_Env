#### Packages and internal function
pkg <- c("rgl", "geometry", "FNN", "colorspace", "ade4")
lapply(pkg, function(x) require(x, quietly = TRUE, character.only = TRUE))

#--> rgl: could be replaced by plotly
#--> geometry: for convexhull
#--> FNN: for k_nearest neighbor
#--> colorspace: for plotting palette
#--> ade4: for pca in order_pal() internal function

"order_pal" <- function(pal){
#--> order the color of the palette
	hsv <- as(hex2RGB(pal), "HSV")
	pca <- dudi.pca(hsv@coords[,1:2], scannf = F)
	hcl <- hclust(dist(pca$li), method = "ward.D2")
	fac <- cutree(hcl, h = 2.25)
	# li1 <- coords_LAB[match(pal, hex_LAB),1]
	lab <- decode_colour(pal, to = "lab")
	knn <- FNN::get.knnx(coords_LAB, lab, k = 1)
	w <- knn$nn.index[,1]
	li1 <- coords_LAB[w,1]
	loc <- unlist(lapply(split(li1, fac), mean))
	v <- hsv@coords[,3]
	names(v) <- pal
	ord <- split(v, fac)[order(loc)]
	ord <- unlist(lapply(lapply(ord, sort), names))
	pal_ord <- pal[ord]
return(pal_ord)
}

"order_pal" <- function(pal){
      hsv <- as(hex2RGB(pal), "HSV")
      pca <- dudi.pca(hsv@coords[,1:2], scannf = F)
      hcl <- hclust(dist(pca$li), method = "ward.D2")
      fac <- cutree(hcl, h = 2.25)
      li1 <- coords_LAB[match(pal, hex_LAB),1]
      loc <- unlist(lapply(split(li1, fac), mean))
      ord <- split(hsv@coords[,3], fac)[order(loc)]
      ord <- unlist(lapply(lapply(ord, sort), base::names))
      pal_ord <- pal[ord]
return(pal_ord)
}

####	####	####	####
#--> Function to initiate the shiny app: simplified version of coords2pal
	#--> adjust the cloud of colours to cloud of data

"initiate" <- function(coords, coords_LAB, hex_LAB, inside = TRUE, s = 1, ds = 0.05, plot = TRUE, reorder = TRUE){

	#--> coords: cloud of data points (1D, 2D, or 3D - if more, only the first 3 dimensions are retained)
	#--> coords_LAB: cloud of color points in euclidean space (the first 3D are used)
	#--> hex_LAB: colors (hex code) associated to coords_LAB points 
	#--> inside: if TRUE, test if the data cloud is inside the colour cloud and repeat resizing of colour cloud until it falls inside
	#--> s: initial value to resize the colour cloud if inside is TRUE
	#--> ds: increement to resize colour cloud after scaling and translation
	#--> plot: if TRUE, plot the two clouds after transformation in 3D space
	#--> reorder: if TRUE, try to reorder color before plotting the palette

	### Class of coords: work with matrix
	if(!is.matrix(coords)){
		local <- as.matrix(coords)
	}	else{
			local <- coords
		}

	### Axes
	nc <- ncol(local)
	cn <- colnames(local)
	if(is.null(cn))	colnames(local) <- paste0("Ax", seq(1,nc))
	if(nc < 3){
		repeat{
			local <- cbind(local, 0)
			nc <- ncol(local)
			if(nc == 3){
				break
			}
		}
	}

	## Reorder axes according to ranges
	index <- apply(local, 2, sd)
	index <- order(index, decreasing = TRUE)
	test <- identical(index, seq(1,nc))
	if(!test){
		print("Axes have been reordered according to standard deviation")
		local <- local[,index]
	}
	
	## Keep only the first 3 axes
	if(nc > 3){
		print("Only the first 3 axes have been conserved")
		local <- local[,1:3]
	}

	## Center cloud of data around the mean and scale the axes according to sd of axis 1
	local <- apply(local, 2, function(x) x-mean(x))
	s1 <- sd(local[,1])
	local <- apply(local, 2, function(x) x/s1)

	### Categories
	nr <- nrow(local)
	rn <- rownames(local)
	if(is.null(rn) | identical(rn, as.character(1:nr))){
		rn <- paste0("C", 1:nr)
		rownames(local) <- rn
	}

	### Move the scaled colour cloud [according to sd of the first axis] in data cloud space 
	if(inside){
		repeat{
			local_LAB <- apply(coords_LAB, 2, function(x) x-mean(x))
			s1_LAB <- sd(local_LAB[,1]) / s
			local_LAB <- apply(local_LAB, 2, function(x) x/s1_LAB)
			chull_LAB <- geometry::convhulln(local_LAB)
			test <- geometry::inhulln(chull_LAB, local)
			if(sum(!test) == 0){
				break
			}
			s <- s + ds
		}
		print(s)
	}	else{
			local_LAB <- apply(coords_LAB, 2, function(x) x-mean(x))
			s1_LAB <- sd(local_LAB[,1]) / s
			local_LAB <- apply(local_LAB, 2, function(x) x/s1_LAB)
			chull_LAB <- geometry::convhulln(local_LAB)
		}

	### Plot the result
	## Build palette
	knn <- FNN::get.knnx(local_LAB[,1:nc], local, k = 1)
	w <- knn$nn.index[,1]
	pal <- hex_LAB[w]
	names(pal) <- rn	

	## Plot the results
	if(plot){
		# Plot the cloud
		rgl::view3d(theta = 25, phi = 25, zoom = 1)
		r <- 1.5 * range(local[,1])
		rgl::plot3d(local, type = "n", xlim = r, ylim = r, zlim = r)
		if(nr <= 300){
			rgl::spheres3d(local, r = 0.1, col = pal)
		}	else{
				rgl::points3d(local, col = pal)
			}
		m <- apply(local, 2, mean)
		rgl::spheres3d(m, r = 0.05)
		chull_LAB.t <- t(chull_LAB)
		rgl::triangles3d(local_LAB[chull_LAB.t,1],local_LAB[chull_LAB.t,2],local_LAB[chull_LAB.t,3],
		 col = "grey", alpha = 0.2)
		rgl::points3d(local_LAB, r = 0.01, col = hex_LAB, alpha = 0.1)
		rgl::cur3d()

		# Plot the palette
		if(nr <= 300){
			if(reorder){
				p <- order_pal(pal)
			}	else{
					p <- pal
				}
			colorspace::swatchplot(p, cvd = TRUE)
			text(seq(0+1/nr/2, 1-1/nr/2, le = nr), 0.9575, names(p), cex = 0.75)
		}
	}
res <- list(pal = pal, coords = local, coords_LAB = local_LAB)
return(res)
}


####	####	####	####
#--> description file
quiet <- TRUE
if(!quiet){
	library(terra)
	source("C:/Users/Seb/Documents/Ma recherche/Papers/2024_Methods in Ecology and Evolution_project/functions/initiate.R")

	dir <- "C:/Users/Seb/Documents/Ma recherche/Papers/2024_Methods in Ecology and Evolution_project/RData/"
	path <- paste0(dir, "LAB.RData")#
	load(path)
	coords_LAB <- as.matrix(coords_LAB)

	#### Test the function on random colors
	### Create data to test the function
	nr <- 10
	coords <- cbind.data.frame(rnorm(nr, 1, 1), rnorm(nr, 0, 10), rnorm(nr, 8, 5))
	names(coords) <- NULL

	### Test the function
	res <- initiate(coords, coords_LAB, hex_LAB, ds = 0.2)

	
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
	
	## Test the function
	res <- initiate(coords, coords_LAB, hex_LAB, ds = 0.1)
	pol <- as.polygons(beck_10km)
	plot(beck_10km, col = res$pal)
	plot(pol, col = res$pal, border = "white")

	### Load data to test the function
	path <-  paste0("C:/Users/Seb/Documents/Ma recherche/Data/Climatic/modified/worldclim_bio_10km.tif")
	worldclim_10km <- rast(path)
	tem_pre <- terra::values(worldclim_10km)[,c(1,12)]
	# tem_pre <- terra::values(worldclim_pca_10km)[,1:3]
	coords <- as.data.frame(tem_pre[-index,])
	coords[,2] <- -coords[,2]**(1/2.25)
	# coords <- apply(coords, 2, scale)
	
	## Test function
	res <- initiate(coords, coords_LAB, hex_LAB, inside = F, s = 1.25)
	plot(coords, pch = 20, cex = 0.1, col = res$pal)
	pal <- rep(NA, nrow(tem_pre))
	pal[-index] <- res$pal
	r <- rast(ext(worldclim_10km), resolution = res(worldclim_10km))
	terra::values(r) <- pal
	plot(r)
#--> try with ecozones
}

