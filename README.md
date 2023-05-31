#### Authors:  

Andrii Zaiats, Megan E. Cattau, Rongsong Liu, David S. Pilliod, Patricia Kaye T. Dumandan, Ahmad Hojatimalekshah, Donna Delparte, T. Trevor Caughlin

---  

#### Overview:

This repository stores the scripts that were used for the article "Multi-scale structural heterogeneity reflects disturbance and recovery in arid shrublands".  

___  

The repository contains data processing, tabular data, figures, modeling results presented in the manuscript. 

Folders:\
    - **data**: locations of randomized plot centroids, burnt/unburnt polygons, original and scaled resolutions of canopy height models, site-level data, rugosity and heterogeneity meteric calculated with DWT.\
    - **figures**: includes .png/ figures from _figures.R_
    
R scripts:\
    - *get_cell_rugosity.R*: extracts canopy rugosity from CHMs.  
    - *get_chm_res.R*: extracts resolutions form CHM and DWT-transformed CHMs.  
    - *scale_decompose.R*: tiles and applies DWT to CHMs.  
    - *scale_depend_disturbance.R*: wildfire effect on structural heterogeneity.  
    - *scale_demography.R*: the effect of structural heterogeneity on recruitment.  
    - *figures.R*: generates figures for the manuscript.
