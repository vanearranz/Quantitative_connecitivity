#### Step 0. Set up the working directory and load all the functions and packages #### 

## Load functions from functions.R
- LoadSpecies3
- pelc.significant.func
- sig.plot 
- sig.plot2

require(dplyr)

#Create an output folder
dir.create("output")

#### Step 1. Simulation of the real data #### 

## X number of similation of n random breaks for each Species (according to the real number of breaks that it was found by the authors)

Asqu2 <- LoadSpecies("Asqu", "Species/Asqu.csv", 10000, 4)
Astu2 <- LoadSpecies("Astu", "Species/Astu.csv", 10000, 7)
Aten2 <- LoadSpecies("Aten", "Species/Aten.csv", 10000, 6)
CglaspA2 <- LoadSpecies("CglapA", "Species/CglaspA.csv", 10000, 2)
CglaspC2 <- LoadSpecies("CglaspC", "Species/CglaspC.csv", 10000, 4)
Corn2 <- LoadSpecies("Corn", "Species/Corn.csv", 10000, 5)
Crad2 <- LoadSpecies("Crad", "Species/Crad.csv", 10000, 6)
Echlo2 <- LoadSpecies("Echlo", "Species/Echlo.csv", 10000, 4)
Hiris2 <- LoadSpecies("Hiris", "Species/Hiris.csv", 10000, 3)
Hsex2 <- LoadSpecies("Hsex", "Species/Hsex.csv", 10000, 3)
Lsma2 <- LoadSpecies("Lsma", "Species/Lsma.csv", 10000, 9)
Paust2 <- LoadSpecies("Paust", "Species/Paust.csv", 10000, 4)
Pcan2 <- LoadSpecies("Pcan", "Species/Pcan.csv", 10000, 2)
Pelo2 <- LoadSpecies("Pelo", "Species/Pelo.csv", 10000, 7)
Pnov2 <- LoadSpecies("Pnov", "Species/Pnov.csv", 10000, 4)
Preg2 <- LoadSpecies("Preg", "Species/Preg.csv", 10000, 2)
Psub2 <- LoadSpecies("Psub", "Species/Psub.csv", 10000, 6)
Sbre2 <- LoadSpecies("Sbre", "Species/Sbre.csv", 10000, 2)
Spel2 <- LoadSpecies("Spel", "Species/Spel.csv", 10000, 8)
Zlut2 <- LoadSpecies("Zlut", "Species/Zlut.csv", 10000, 3)
Zsub2 <- LoadSpecies("Zsub", "Species/Zsub.csv", 10000, 5)

#### Steo 2. Create the matrix (Bmat) for each subset of analysis ####
# Example with all species with Fst data #

Segments_Asqu<- read.csv("output/Segments_Asqu.csv")
Asqu_F <- as.matrix(Segments_Asqu[-(1:4)])
rm("Segments_Asqu")

Segments_Astu <- read.csv("output/Segments_Astu.csv")
Astu_F <- as.matrix(Segments_Astu[-(1:4)])
rm("Segments_Astu")

Segments_Aten <- read.csv("output/Segments_Aten.csv")
Aten_F <- as.matrix(Segments_Aten[-(1:4)])
rm("Segments_Aten")

Segments_CglaspA <- read.csv("output/Segments_CglapA.csv")
CglaspA_F <- as.matrix(Segments_CglaspA[-(1:4)])
rm("Segments_CglaspA")

Segments_CglaspC <- read.csv("output/Segments_CglaspC.csv")
CglaspC_F <- as.matrix(Segments_CglaspC[-(1:4)])
rm("Segments_CglaspC")

Segments_Corn <- read.csv("output/Segments_Corn.csv")
Corn_F <- as.matrix(Segments_Corn[-(1:4)])
rm("Segments_Corn")

Segments_Crad <- read.csv("output/Segments_Crad.csv")
Crad_F <- as.matrix(Segments_Crad[-(1:4)])
rm("Segments_Crad")

Segments_Echlo <- read.csv("output/Segments_Echlo.csv")
Echlo_F <- as.matrix(Segments_Echlo[-(1:4)])
rm("Segments_Echlo")

Segments_Hiri <- read.csv("output/Segments_Hiris.csv")
Hiris_F <- as.matrix(Segments_Hiri[-(1:4)])
rm("Segments_Hiri")

Segments_Hsex <- read.csv("output/Segments_Hsex.csv")
Hsex_F <- as.matrix(Segments_Hsex[-(1:4)])
rm("Segments_Hsex")

Segments_Lsma <- read.csv("output/Segments_Lsma.csv")
Lsma_F <- as.matrix(Segments_Lsma[-(1:4)])
rm("Segments_Lsma")

Segments_Paust <- read.csv("output/Segments_Paust.csv")
Paust_F <- as.matrix(Segments_Paust[-(1:4)])
rm("Segments_Paust")

Segments_Pcan <- read.csv("output/Segments_Pcan.csv")
Pcan_F <- as.matrix(Segments_Pcan[-(1:4)])
rm("Segments_Pcan")

Segments_Pelo <- read.csv("output/Segments_Pelo.csv")
Pelo_F <- as.matrix(Segments_Pelo[-(1:4)])
rm("Segments_Pelo")

Segments_Pnov <- read.csv("output/Segments_Pnov.csv")
Pnov_F <- as.matrix(Segments_Pnov[-(1:4)])
rm("Segments_Pnov")

Segments_Preg <- read.csv("output/Segments_Preg.csv")
Preg_F <- as.matrix(Segments_Preg[-(1:4)])
rm("Segments_Preg")

Segments_Psub <- read.csv("output/Segments_Psub.csv")
Psub_F <- as.matrix(Segments_Psub[-(1:4)])
rm("Segments_Psub")

Segments_Sbre <- read.csv("output/Segments_Sbre.csv")
Sbre_F <- as.matrix(Segments_Sbre[-(1:4)])
rm("Segments_Sbre")

Segments_Spel <- read.csv("output/Segments_Spel.csv")
Spel_F <- as.matrix(Segments_Spel[-(1:4)])
rm("Segments_Spel")

Segments_Zlut <- read.csv("output/Segments_Zlut.csv")
Zlut_F <- as.matrix(Segments_Zlut[-(1:4)])
rm("Segments_Zlut")

Segments_Zsub <- read.csv("output/Segments_Zsub.csv")
Zsub_F <- as.matrix(Segments_Zsub[-(1:4)])
rm("Segments_Zsub")

# Add all the species matrices relevant to the desire analysis 
## In this example all the 21 species using Fst data

Bmat_F_all <- Asqu_F+Astu_F+Aten_F+CglaspA_F+CglaspC_F+Corn_F+Crad_F+Echlo_F+Hiris_F+Hsex_F+Lsma_F+Paust_F+Pcan_F+Pelo_F+Pnov_F+Preg_F+Psub_F+Sbre_F+Spel_F+Zlut_F+Zsub_F


#### Step 3. Test for statistical significance of the real break for each segment ####

## Add a vector with the real data 
realvec <- read.csv("realvec/real_Fst_all.csv", header = FALSE)

res.out.allF2 <- pelc.significant.func(mat= Bmat_F_all, rvec = realvec, siglevel = c(0.8, 0.9, 0.95))


#### Step 4. Plot the significance levels in each segment #### 

## Add the geographic information of each segment to plot in a map
segdat <- read.csv("Species/segdat.csv", as.is=T)

# Plot with sig.plot or sig.plot2 with different levels of significance
sig.plot(res.out.allF2$sig, segdat = segdat)
sig.plot2(res.out.allF2$sig, segdat = segdat)


### EXTRA. More Bmat for all the Fst analysis ####

Bmat_F_low <- Asqu_F+Aten_F+Corn_F+Crad_F+Hiris_F+Lsma_F+Zsub_F+Spel_F+CglaspA_F+CglaspC_F+Sbre_F+Zlut_F 
Bmat_F_high <- Hsex_F+Astu_F+Echlo_F+Paust_F+Pcan_F+Pelo_F+Pnov_F+Preg_F+Psub_F

Bmat_Hab_low <- Aten_F+Asqu_F+Pelo_F+Hsex_F+Preg_F+Crad_F+CglaspA_F+CglaspC_F+Lsma_F+Pcan_F+Sbre_F+Zlut_F+Zsub_F
Bmat_Hab_sub <- Echlo_F+Hiris_F+Paust_F+Psub_F+Pnov_F+Astu_F
Bmat_Hab_high <- Spel_F+Corn_F

Bmat_microsat <- Aten_F+Paust_F+Psub_F+Pnov_F+Echlo_F
Bmat_mitocon <- Pelo_F+Asqu_F+Preg_F+Astu_F+Corn_F+Crad_F+CglaspA_F+CglaspC_F+Hiris_F+Lsma_F+Pcan_F+Sbre_F+Spel_F+Zlut_F+Zsub_F



