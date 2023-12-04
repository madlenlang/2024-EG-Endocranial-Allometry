###### RUNNING ANALYSES WITH SLIDING LANDMARKS ######

library(geomorph)
library(phytools) 
library(Arothron)

setwd("D:/Lang - GMM Endocast Project/Extant Paper/Extant GM Data")

## LOAD FILEs ## 
	Classifier <- read.csv("EG-GMM-Classifier.csv") 
	mytree <- read.nexus("mytree.nex")	
	coordsraw <- read.amira.dir("EG_Landmarks_CombinedSets",42)	# read raw coordinate data from file
	omit <- c(5,14,23,25,26,37)	# removing duplicate landmarks
	coords.all <-coordsraw[-omit,,]
	
## Procurstes Transformation/Data set up ##	
	sliders <- rbind(define.sliders(26:31), define.sliders(31:36))	
	Y.gpa <- gpagen(coords.all, curves = sliders)	
	GDF <- geomorph.data.frame(coords = coords.all, landmarks = Y.gpa$coords, Csize = Y.gpa$Csize, species = Classifier$Species, order = Classifier$Order, clade = Classifier$Clade, lago = Classifier$Lago)

## DERMOPTERA and TARSIERS removed for ANOVAs ##
	class2 <- Classifier[-c(33,51,124,125),]	# removing dermoptera and tarsiers from classifier
	omit <- c(33, 51, 124, 125)
	coords.omit <- Y.gpa$coords[,,-omit]		# removing dermoptera and tarsiers  from Y.gpa
	Csize.omit <- Y.gpa$Csize[-omit]
	GDF2 <- geomorph.data.frame(landmarks2 = coords.omit, Csize2 = Csize.omit, species2 = class2$Species, order2 = class2$Order, clade2 = class2$Clade, lago2 = class2$Lago)

## RUNNING PANOVA and PAIRWISE TESTS - CLADE ##
	PANOVA <- procD.lm(landmarks2~clade2, data = GDF2, iter = 9999, RRPP = TRUE)
		summary(PANOVA)
	PW <- pairwise(PANOVA, groups = GDF2$clade2, covariate = NULL)
		summary(PW)
		
## Supplementary Materials - RUNNING ANOVA and PAIRWISE TESTS - CLADE [leporidae and ochotonidae seperated] ##
	PANOVA.lag <- procD.lm(landmarks2~lago2, data = GDF2, iter = 9999, RRPP = TRUE)
		summary(PANOVA.lag)
	PW.lag <- pairwise(PANOVA.lag, groups = GDF2$lago2, covariate = NULL)
		summary(PW.lag)
		
## Supplementary Materials - RUNNING ANOVA and PAIRWISE TESTS - ORDER ##
	PANOVA.ord <- procD.lm(landmarks2~order2, data = GDF2, iter = 9999, RRPP = TRUE)
		summary(PANOVA.ord)
	PW.ord <- pairwise(PANOVA.ord, groups = GDF2$order2, covariate = NULL)
		summary(PW.ord)

## Homogeneity of Slopes Test - Clade ##
	fit.unique <- procD.lm(landmarks2~log(Csize.omit)*clade2,data = GDF2,iter = 9999, RRPP = FALSE) 
		summary(fit.unique)
	fit.common <- procD.lm(landmarks2~log(Csize.omit)+clade2,data = GDF2,iter = 9999, RRPP = FALSE)
		summary(fit.common)
	anova(fit.common, fit.unique, print.progress = FALSE)
	
	PW<- pairwise(fit = fit.unique, fit.null = fit.common, groups = GDF2$clade2, covariate = Csize.omit)
		summary(PW)
		summary(PW, test.type = "DL", confidence = 0.95, stat.table = FALSE)	# diffs in slope vector lengths (magnitude)
		summary(PW, test.type = "VC", angle.type = "deg", stat.table = FALSE)	# diffs in slipe vector orientation (direction)

## Homogeneity of Slopes Test with Phylo context - Clade ##
	Species.tar <- Classifier$Species[Classifier$Clade == "Tarsioidea"]
	Species.der <- Classifier$Species[Classifier$Clade == "Dermoptera"]
	mytree2 <- drop.tip(mytree, c(Species.der, Species.tar))		# dropping dermoptera and tariers from tree
		
	names(coords.omit) = GDF2$species2								# assign species names to Y.gpa
	dimnames(coords.omit)[[3]] <- gsub (".txt", "", GDF2$species2)	# assigning species names to the Y.GPA 
	sort(dimnames(coords.omit)[[3]])==sort(mytree2$tip.label)		# sort and identify which names are not aligned

	GDFphy2 <- geomorph.data.frame(landmarks2 = coords.omit, phy=mytree2, Csize2 = Csize.omit, clade2 = class2$Clade)

	fitpgls.unique <- procD.pgls(GDFphy2$landmarks2~log(GDFphy2$Csize2)*clade2, phy = GDFphy2$phy, data = GDFphy2, iter = 9999)	
		anova(fitpgls.unique )	
	fitpgls.common <- procD.pgls(GDFphy2$landmarks2~log(GDFphy2$Csize2)+clade2, phy = GDFphy2$phy, data = GDFphy2, iter = 9999)	
		anova(fitpgls.common)	
	
	anova(fitpgls.common, fitpgls.unique, print.progress = FALSE)
	
	PWphy<- pairwise(fit = fitpgls.unique, fit.null = fitpgls.common, groups = GDF2$clade2, covariate = Csize.omit)
		summary(PWphy)
		summary(PWphy, test.type = "DL", confidence = 0.95, stat.table = FALSE)		# diffs in slope vector lengths (magnitude)
		summary(PWphy, test.type = "VC", angle.type = "deg", stat.table = FALSE)	# diffs in slipe vector orientation (direction)
		
## Overall Phylogenetic Signal in Shape ##
	names(Y.gpa$Csize) = GDF$species								
	dimnames(Y.gpa$coords)[[3]] <- gsub (".txt", "", GDF$species)	
	sort(dimnames(Y.gpa$coords)[[3]])==sort(mytree$tip.label)		
	
	phypl <- physignal(A = Y.gpa$coords,phy = mytree,iter = 9999)
		summary(phypl)

## Overall Phylogenetic Signal in Size ##
	phypl <- physignal(A = Y.gpa$Csize,phy = mytree,iter = 9999)
		summary(phypl)
		
## Setting up data for phylogenetic and other analyses ##
	lm.lag <- GDF$landmarks[ , , GDF$clade == "Lagomorpha"]		# seperate coordinates from common Y.gpa for each clade of interest 
	lm.pla <- GDF$landmarks[ , , GDF$clade == "Platyrrhini"]
	lm.rod <- GDF$landmarks[ , , GDF$clade == "Rodentia"]	
	lm.sca <- GDF$landmarks[ , , GDF$clade == "Scandentia"]
	lm.str <- GDF$landmarks[ , , GDF$clade == "Strepsirrhini"]
	
	csize.pla <- GDF$Csize[GDF$clade == "Platyrrhini"]	# seperate Csize from common Y.gpa coordinates for each clade of interest
	csize.lag <- GDF$Csize[GDF$clade == "Lagomorpha"]	
	csize.rod <- GDF$Csize[GDF$clade == "Rodentia"]	
	csize.sca <- GDF$Csize[GDF$clade == "Scandentia"]
	csize.str <- GDF$Csize[GDF$clade == "Strepsirrhini"]	
		
	Species.der <- Classifier$Species[Classifier$Clade == "Dermoptera"]		# assign species label used to match tip names
	Species.lag <- Classifier$Species[Classifier$Clade == "Lagomorpha"]		
	Species.pla <- Classifier$Species[Classifier$Clade == "Platyrrhini"]
	Species.rod <- Classifier$Species[Classifier$Clade == "Rodentia"]
	Species.sca <- Classifier$Species[Classifier$Clade == "Scandentia"]
	Species.str <- Classifier$Species[Classifier$Clade == "Strepsirrhini"]
	Species.tar <- Classifier$Species[Classifier$Clade == "Tarsioidea"]
		
	lag.tree <- drop.tip(mytree, c(Species.pla, Species.tar, Species.rod, Species.sca, Species.str, Species.der))	# drop tips from tree for groups not of interest
	pla.tree <- drop.tip(mytree, c(Species.lag, Species.rod, Species.sca, Species.tar, Species.der, Species.str))
	rod.tree <- drop.tip(mytree, c(Species.lag, Species.tar, Species.pla, Species.sca, Species.str, Species.der))
	sca.tree <- drop.tip(mytree, c(Species.lag, Species.tar, Species.rod, Species.pla, Species.str, Species.der))
	str.tree <- drop.tip(mytree, c(Species.lag, Species.tar, Species.rod, Species.sca, Species.pla, Species.der))
	
	dimnames(lm.lag)[[3]] <- gsub (".txt", "", Species.lag)	# assigning species names to the Y.GPA 
	dimnames(lm.pla)[[3]] <- gsub (".txt", "", Species.pla)	
	dimnames(lm.rod)[[3]] <- gsub (".txt", "", Species.rod)	
	dimnames(lm.sca)[[3]] <- gsub (".txt", "", Species.sca)	
	dimnames(lm.str)[[3]] <- gsub (".txt", "", Species.str)	
	
	names(csize.lag) = Species.lag	# assigning species names to Csize
	names(csize.pla) = Species.pla
	names(csize.rod) = Species.rod
	names(csize.sca) = Species.sca
	names(csize.rod) = Species.rod
	names(csize.str) = Species.str
	
## Clade Specific Physignals in Shape from common y.gpa ##
	phy.lm.lag <- physignal(A = lm.lag,phy = lag.tree,iter = 9999)
		summary(phy.lm.lag)
	phy.lm.pla <- physignal(A = lm.pla,phy = pla.tree,iter = 9999)
		summary(phy.lm.pla)
	phy.lm.rod <- physignal(A = lm.rod,phy = rod.tree,iter = 9999)
		summary(phy.lm.rod)
	phy.lm.sca <- physignal(A = lm.sca,phy = sca.tree,iter = 9999)
		summary(phy.lm.sca)
	phy.lm.str <- physignal(A = lm.str,phy = str.tree,iter = 9999)
		summary(phy.lm.str)

## Clade Specific Physignals in Size from common y.gpa ##	
	phy.cs.lag <- physignal(A = csize.lag,phy = lag.tree,iter = 9999)
		summary(phy.cs.lag)	
	phy.cs.pla <- physignal(A = csize.pla,phy = pla.tree,iter = 9999)
		summary(phy.cs.pla )
	phy.cs.rod <- physignal(A = csize.rod,phy = rod.tree,iter = 9999)
		summary(phy.cs.rod)
	phy.cs.sca <- physignal(A = csize.sca,phy = sca.tree,iter = 9999)
		summary(phy.cs.sca)
	phy.cs.str <- physignal(A = csize.str,phy = str.tree,iter = 9999)
		summary(phy.cs.str)

## Common ANOVAs (linear model and PGLS) of shape and size ##
	fit <- procD.lm(landmarks~log(Y.gpa$Csize),data = GDF,iter = 9999, RRPP = TRUE) 
		summary(fit)
	
	names(Y.gpa$Csize) = GDF$species								# assign species names to Y.gpa$Csize
	dimnames(Y.gpa$coords)[[3]] <- gsub (".txt", "", GDF$species)	# assigning species names to the Y.GPA 
	sort(dimnames(Y.gpa$coords)[[3]])==sort(mytree$tip.label)		# identify which names are not aligned
	GDFphy <- geomorph.data.frame(landmarks = Y.gpa$coords, phy=mytree, Csize = Y.gpa$Csize, clade = Classifier$Clade)

	fitpgls <- procD.pgls(GDFphy$landmarks~log(GDFphy$Csize), phy = GDFphy$phy, data = GDFphy, iter = 9999)	
		anova(fitpgls)

## Clade ANOVAs (linear model and PGLS) of shape and size using clade specific procrustes transformations ##
## Lagomorpha ##
	coords.lag <- GDF$coords[ , , GDF$clade == "Lagomorpha"]
	Y.gpa.lag <- gpagen(coords.lag, curves = sliders)	
	
	names(Y.gpa.lag$Csize) = Species.lag								
	dimnames(Y.gpa.lag$coords)[[3]] <- gsub (".txt", "", Species.lag)	
	sort(dimnames(Y.gpa.lag$coords)[[3]])==sort(lag.tree$tip.label)		
	gdf.lag <- geomorph.data.frame(phy=lag.tree, lmlag = Y.gpa.lag$coords, Csize.ord = Y.gpa.lag$Csize, species.ord = Species.lag)
	
	fit.unique.lag <- procD.lm(Y.gpa.lag$coords~log(Y.gpa.lag$Csize),data = gdf.lag,iter = 9999, RRPP = TRUE) 
		summary(fit.unique.lag)
			
	fitpgls.lag <- procD.pgls(Y.gpa.lag$coords~log(Y.gpa.lag$Csize), phy = lag.tree, data = gdf.lag, iter = 9999)	
		anova(fitpgls.lag)	
			
## Platyrrhini ##
	coords.pla <- GDF$coords[ , , GDF$clade == "Platyrrhini"]
	Y.gpa.pla <- gpagen(coords.pla, curves = sliders)	
	
	names(Y.gpa.pla$Csize) = Species.pla								
	dimnames(Y.gpa.pla$coords)[[3]] <- gsub (".txt", "", Species.pla)	
	sort(dimnames(Y.gpa.pla$coords)[[3]])==sort(pla.tree$tip.label)		
	gdf.pla <- geomorph.data.frame(phy=pla.tree, lmpla = Y.gpa.pla$coords, Csize.ord = Y.gpa.pla$Csize, species.ord = Species.pla)
	
	fit.unique.pla <- procD.lm(Y.gpa.pla$coords~log(Y.gpa.pla$Csize), data = gdf.pla, iter = 9999, RRPP = TRUE) 
		summary(fit.unique.pla)
			
	fitpgls.pla <- procD.pgls(Y.gpa.pla$coords~log(Y.gpa.pla$Csize), phy = pla.tree, data = gdf.pla, iter = 9999)	
		anova(fitpgls.pla)	

## Rodentia ##
	coords.rod <- GDF$coords[ , , GDF$clade == "Rodentia"]
	Y.gpa.rod <- gpagen(coords.rod, curves = sliders)	

	names(Y.gpa.rod$Csize) = Species.rod								
	dimnames(Y.gpa.rod$coords)[[3]] <- gsub (".txt", "", Species.rod)	
	sort(dimnames(Y.gpa.rod$coords)[[3]])==sort(rod.tree$tip.label)		
	gdf.rod <- geomorph.data.frame(phy=rod.tree, lmrod = Y.gpa.rod$coords, Csize.ord = Y.gpa.rod$Csize, species.ord = Species.rod)
	
	fit.unique.rod <- procD.lm(Y.gpa.rod$coords~log(Y.gpa.rod$Csize), data = gdf.rod, iter = 9999, RRPP = TRUE) 
		summary(fit.unique.rod)
			
	fitpgls.rod <- procD.pgls(Y.gpa.rod$coords~log(Y.gpa.rod$Csize), phy = rod.tree, data = gdf.rod, iter = 9999)	
		anova(fitpgls.rod)	

## Scandentia ##
	coords.sca <- GDF$coords[ , , GDF$clade == "Scandentia"]
	Y.gpa.sca <- gpagen(coords.sca, curves = sliders)	
	
	names(Y.gpa.sca$Csize) = Species.sca								
	dimnames(Y.gpa.sca$coords)[[3]] <- gsub (".txt", "", Species.sca)	
	sort(dimnames(Y.gpa.sca$coords)[[3]])==sort(sca.tree$tip.label)		
	gdf.sca <- geomorph.data.frame(phy=sca.tree, lmsca = Y.gpa.sca$coords, Csize.ord = Y.gpa.sca$Csize, species.ord = Species.sca)
	
	fit.unique.sca <- procD.lm(Y.gpa.sca$coords~log(Y.gpa.sca$Csize), data = gdf.sca, iter = 9999, RRPP = TRUE) 
		summary(fit.unique.sca)
			
	fitpgls.sca <- procD.pgls(Y.gpa.sca$coords~log(Y.gpa.sca$Csize), phy = sca.tree, data = gdf.sca, iter = 9999)	
		anova(fitpgls.sca)	

## Strepsirrhini ##
	coords.str <- GDF$coords[ , , GDF$clade == "Strepsirrhini"]
	Y.gpa.str <- gpagen(coords.str, curves = sliders)	
	
	names(Y.gpa.str$Csize) = Species.str								
	dimnames(Y.gpa.str$coords)[[3]] <- gsub (".txt", "", Species.str)	
	sort(dimnames(Y.gpa.str$coords)[[3]])==sort(str.tree$tip.label)		
	gdf.str <- geomorph.data.frame(phy=str.tree, lmstr = Y.gpa.str$coords, Csize.ord = Y.gpa.str$Csize, species.ord = Species.str)
	
	fit.unique.str <- procD.lm(Y.gpa.str$coords~log(Y.gpa.str$Csize), data = gdf.str, iter = 9999, RRPP = TRUE) 
		summary(fit.unique.str)
			
	fitpgls.str <- procD.pgls(Y.gpa.str$coords~log(Y.gpa.str$Csize), phy = str.tree, data = gdf.str, iter = 9999)	
		anova(fitpgls.str)

## Common ANCOVAs (linear model and PGLS) of shape and size with clade as a covariate - Dermoptera and Tarsiers removed ##
	fit <- procD.lm(landmarks2~log(Csize2)*clade2,data = GDF2,iter = 9999, RRPP = TRUE) 
		summary(fit)
	
	fitpgls <- procD.pgls(GDFphy2$landmarks2~log(GDFphy2$Csize2)*clade2, phy = GDFphy2$phy, data = GDFphy2, iter = 9999)	
		anova(fitpgls)

## RUNNING SIMPLE PCA ##
	PCA <- gm.prcomp(Y.gpa$coords)

## Broken Stick Model (Wilson et al., 2023; https://doi.org/10.1038/s41467-023-38365-0) 
	ev<-PCA$sdev^2
	
	bsm <- function(ev) { 
		# Broken stick model (MacArthur 1957)
		n = length(ev)
			bsm = data.frame(j=seq(1:n), p=0)
		bsm$p[1] = 1/n
		for (i in 2:n) bsm$p[i] = bsm$p[i-1] + (1/(n + 1 - i))
		bsm$p = 100*bsm$p/n
		# Plot eigenvalues and % of variation for each axis
		op = par(mfrow=c(2,1),omi=c(0.1,0.3,0.1,0.1), mar=c(1, 1, 1, 1))
		barplot(ev, main="Eigenvalues", col="lightskyblue", las=2)
		abline(h=mean(ev), col="red")
		legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")
		barplot(t(cbind(100*ev/sum(ev), bsm$p[n:1])), beside=TRUE, 
			main="% variation", col=c("lightskyblue",2), las=2)
		legend("topright", c("% eigenvalue", "Broken stick model"), 
			pch=15, col=c("lightskyblue",2), bty="n")
		par(op)
	}
	bsm(ev)

## PCA Plot - CLADE ##
	#Sequence: Dermoptera, Lagomorpha, Platyrrhini, Rodentia, Scandentia, Strepsirrhini, Tarsioidea
	cols <- c("goldenrod2", "brown1", "steelblue1", "mediumpurple3", "yellowgreen", "midnightblue", "violetred")	
	shapes <- c(25, 22, 23, 24, 21, 23, 23)
				
	plot(PCA, pch = shapes[as.factor(Classifier$Clade)], cex = 1.25, col = cols[as.factor(Classifier$Clade)], bg = cols[as.factor(Classifier$Clade)])
	legend("topleft", pch = shapes[as.factor(Classifier$Clade)], cex = 1.25,  legend = unique(Classifier$Clade), col = cols[as.factor(Classifier$Clade)], bty = "n")

## Supplementary Materials - PCAs SHOWING OTHER AXES - CLADE ##
	plot(PCA, axis1 = 3, axis2 = 4, pch = shapes[as.factor(Classifier$Clade)], cex = 1.25, col = cols[as.factor(Classifier$Clade)], bg = cols[as.factor(Classifier$Clade)])
	legend("topleft", pch = shapes[as.factor(Classifier$Clade)], cex = 1.25,  legend = unique(Classifier$Clade), col = cols[as.factor(Classifier$Clade)], bty = "n")

## VISUALIZATION of allometry with with PredLine - nonphylo ##
	fit <- procD.lm(landmarks~log(Y.gpa$Csize)*clade,data = GDF,iter = 9999) 
	plotAllometry(fit, size = Y.gpa$Csize, logsz = TRUE, method = "PredLine", pch = shapes[as.factor(Classifier$Clade)], col = cols[as.factor(Classifier$Clade)], bg = cols[as.factor(Classifier$Clade)], cex = 1.5)

## VISUALIZATION of allometry with RegScore - nonphylo ##
	fit <- procD.lm(landmarks~log(Y.gpa$Csize),data = GDF,iter = 9999) 
	plotAllometry(fit, size = Y.gpa$Csize, logsz = TRUE, method = "RegScore", pch = shapes[as.factor(Classifier$Clade)], col = cols[as.factor(Classifier$Clade)], bg = cols[as.factor(Classifier$Clade)], cex = 1.5)

## Wireframe ##
	wf <- read.morphologika("wireframe.txt")	# define a wireframe
	links <- wf$wireframe

## Wireframes warps for allometry based on RegScore per clade ##
## Lagomorpga ##
	fit.lag <- procD.lm(lm.lag~log(csize.lag),iter = 9999, RRPP = TRUE) 
	pAllo.lag <- plotAllometry(fit.lag, size = csize.lag, logsz = TRUE, method = "RegScore", pch = 19, col = "black", cex = 1.5)											# plot allometry using RegScore	
	preds.lag <- shape.predictor(pAllo.lag $GM$fitted, x= pAllo.lag$RegScore, Intercept = FALSE, predmin = min(pAllo.lag$RegScore),  predmax = max(pAllo.lag$RegScore))		# preditct shape from RegScore
	plotRefToTarget(preds.lag$predmin, preds.lag$predmax, method="points", links = links)																					# show wireframe warp; M1 is grey, M2 is black

## Platyrrhini ##
	fit.pla <- procD.lm(lm.pla~log(csize.pla),iter = 9999, RRPP = TRUE) 
	pAllo.pla <- plotAllometry(fit.pla, size = csize.pla, logsz = TRUE, method = "RegScore", pch = 19, col = "black", cex = 1.5)						
	preds.pla <- shape.predictor(pAllo.pla$GM$fitted, x= pAllo.pla$RegScore, Intercept = FALSE, predmin = min(pAllo.pla$RegScore),  predmax = max(pAllo.pla$RegScore))		
	plotRefToTarget(preds.pla$predmin, preds.pla$predmax, method="points", links = links)	

## Rodentia ##
	fit.rod <- procD.lm(lm.rod~log(csize.rod),iter = 9999, RRPP = TRUE) 
	pAllo.rod <- plotAllometry(fit.rod, size = csize.rod, logsz = TRUE, method = "RegScore", pch = 19, col = "black", cex = 1.5)						
	preds.rod <- shape.predictor(pAllo.rod$GM$fitted, x= pAllo.rod$RegScore, Intercept = FALSE, predmin = min(pAllo.rod$RegScore),  predmax = max(pAllo.rod$PredLine))		
	plotRefToTarget(preds.rod$predmin, preds.rod$predmax, method="points", links = links)	
		
## Scandentia ##
	fit.sca <- procD.lm(lm.sca~log(csize.sca),iter = 9999, RRPP = TRUE) 
	pAllo.sca <- plotAllometry(fit.sca, size = csize.sca, logsz = TRUE, method = "RegScore", pch = 19, col = "black", cex = 1.5)						
	preds.sca <- shape.predictor(pAllo.sca$GM$fitted, x= pAllo.sca$RegScore, Intercept = FALSE, predmin = min(pAllo.sca$RegScore),  predmax = max(pAllo.sca$RegScore))		
	plotRefToTarget(preds.sca$predmin, preds.sca$predmax, method="points", links = links)	
	
## Strepsirrhini ##
	fit.str <- procD.lm(lm.str~log(csize.str),iter = 9999, RRPP = TRUE) 
	pAllo.str <- plotAllometry(fit.str, size = csize.str, logsz = TRUE, method = "RegScore", pch = 19, col = "black", cex = 1.5)						
	preds.str <- shape.predictor(pAllo.str$GM$fitted, x= pAllo.str$RegScore, Intercept = FALSE, predmin = min(pAllo.str$RegScore),  predmax = max(pAllo.str$RegScore))		
	plotRefToTarget(preds.str$predmin, preds.str$predmax, method="points", links = links)	
	
	#rgl.snapshot("Allo-CLADE-lat.png", fmt = "png")
	#rgl.snapshot("Allo-CLADE-dor.png", fmt = "png")

## Surface warps for PCA ##
    ref <-mshape(Y.gpa$coords)	# get coords for the mean shape from the transformed landmark matrix
    var <- Y.gpa[["coords"]]
    avg <-mshape(var)			# calculate the mean shape
    findMeanSpec(var)			# identify species closest to mean
   
    mean.ply <- read.ply("Paraxerus_cepapi.ply", ShowSpecimen=FALSE)    # load the .ply for the mean species
    mean.coords <- coords.all[,,94]										# get the untransformed coordinates for the mean specimen/species (38 is the sequence number for the mean species (i.e. Tamiasciurus hudsonicus))
    ygpa.mc <- Y.gpa$coords[,,94]										# get transformed coords (from Y.gpa) for the mean specimen/species
       
    averagemesh <- warpRefMesh(mean.ply, mean.coords, avg, color=NULL, centered=FALSE)	# create a mesh of the average shape
   
    PC1min <- PCA$shapes$shapes.comp1$min	# identify the coords for the min of PC1 (comp1)
	PC1max <- PCA$shapes$shapes.comp1$max
	PC2min <- PCA$shapes$shapes.comp2$min
	PC2max <- PCA$shapes$shapes.comp2$max
   
	PC1minwarp <- plotRefToTarget(M1 = ygpa.mc, M2 = PC1min, mesh = averagemesh, method = "surface", color = "lightseagreen")	
		rgl.snapshot( "PC1minwarp", fmt = "png", top = TRUE )
	PC1maxwarp <- plotRefToTarget(M1 = ygpa.mc, M2 = PC1max, mesh = averagemesh, method = "surface", color = "lightseagreen")	
		rgl.snapshot( "PC1maxwarp", fmt = "png", top = TRUE )
	PC2minwarp <- plotRefToTarget(M1 = ygpa.mc, M2 = PC2min, mesh = averagemesh, method = "surface", color = "lightseagreen")	
		rgl.snapshot( "PC2minwarp", fmt = "png", top = TRUE )
	PC2maxwarp <- plotRefToTarget(M1 = ygpa.mc, M2 = PC2max, mesh = averagemesh, method = "surface", color = "lightseagreen")	
		rgl.snapshot( "PC2maxwarp", fmt = "png", top = TRUE )
		
