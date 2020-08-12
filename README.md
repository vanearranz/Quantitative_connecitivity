# Quantitative_connecitivity
This repository contains the data and scripts we used to process the comparative phylogeographic analysis in the manuscript "Geographic concordance of genetic barriers in New Zealand coastal marine species". 

Those scripts are based on Pelc et al. 2009 paper with our own dataset from previously published and new generated population genetic studies of 21 marine species from New Zealand coast. We generated our custom scripts and statistical test provided in this repository. 

# Input files 

### Species folder

- Segments.csv with a total of 92 rarified sampling locations therefore divided the New Zealand (NZ) coastline into 92 between-site coastal segments. The New Zealand coastline was then linearized using simplified coastal distances between all rarified sampling points. The starting point (0 km) was Cape Reinga at the northern tip of the North Island. 
- segdat.csv with the geographic coordinates (latitude and longitude) of each of the 92 segments to plot the results in a NZ map. 
- 21 "Species_name".csv files with the starting and end point of each of the sampling segments of each population genetic study. Those starting and end points are based on the previously linearized NZ coast. 
  
### realvec folder

The number of genetic breaks in each of the 92 segments across all species was recorded based on different fixation indexes (Fst and Phist). Different .csv files correspond to the different groups for analysing the data. 

- real_Fst_all.csv : genetic breaks found including all 21 species based on Fst analysis
- real_Fst_high : genetic breaks found including 9 high-dispersal species based on Fst analysis
- real_Fst_low : genetic breaks found including 12 low-dispersal species based on Fst analysis
- real_Hab_high : genetic breaks found including 2 species inhabiting the high-intertidal zone based on Fst analysis
- real_Hab_low : genetic breaks found including 13 species inhabiting the low-intertidal zone based on Fst analysis
- real_hab_sub : genetic breaks found including 6 species inhabiting the subtidal zone based on Fst analysis
- real_microsats : genetic breaks found including 5 species studied using microsatellites DNA markers based on Fst analysis
- real_mitocon : genetic breaks found including 16 species studied using mitochondrial DNA markers based on Fst analysis
- real_Phist_all : genetic breaks found including 13 species based on Phist analysis
- real_Phist_high : genetic breaks found including 4 high-dispersal species based on Phist analysis
- real_Phist_low : genetic breaks found including 9 low-dispersal species based on Phist analysis

# Scripts 
- functions.R : functions to load in the R environment to generate the random simulations of the data set, perform the randomisation test and generate the plots for the Figures found in the manuscript. 

- NZ_quantitative_connectivity.R : R Script example to perform the analysis and reproduce the results.  

# Additional analysis 

### Complete_analysis folder
Use the scripts in this folder to generate all the simulation matrices performed and analysed in the manuscript. 


## If you use these data files or scripts in future work, please cite:
XXXXXX

# References 
Pelc, R. A., R. R. Warner, and S. D. Gaines. "Geographical patterns of genetic structure in marine species with contrasting life histories." Journal of Biogeography 36, no. 10 (2009): 1881-1890.
