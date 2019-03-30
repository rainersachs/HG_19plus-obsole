# all 6 core files (2 .csv and 4 .R double checked and functional 1/14/2019)
# Copyright:    (C) 2017-2019 Sachs Undergraduate Research Apprentice Program (URAP) class  at UCB
#               This program and its accompanying materials are distributed 
#               under the terms of the GNU General Public License v3. As detailed in that
#               license no warranty, explicit or implied, comes with this suite of R scripts
# Filename:     dataAndInfo.R 
# Purpose:      Concerns radiogenic mouse Harderian gland (HG) tumorigenesis. Loads 
#               ion and tumor prevalence data from CSV files. It is part of the 
#               customized source code for the NASAmouseHG project.
# Contact:      Rainer K. Sachs 
# Website:      https://github.com/rainersachs/mouseHG_Chang_2019plus
# Attribution:  This R script was developed at UC Berkeley. Written by Dae Woong 
#               Ham Summer 2017. Additions, corrections, changes, quality 
#               control, reorganization by Edward Huang, Yimin Lin, Mark Ebert,
#               Yunzhi Zhang and Ray Sachs 2017-2019.

# Relevant references and abbreviations:
#
#   ".93Alp" = Alpen et al. "Tumorigenic potential of high-Z, high-LET charged-
#                           particle radiations." Rad Res 136:382-391 (1993).
#
#   ".94Alp" = Alpen et al. "Fluence-based relative biological effectiveness for
#                           charged particle carcinogenesis in mouse Harderian 
#                           gland." Adv Space Res 14(10): 573-581. (1994).  
#
#   "16Chang" = Chang et al. "Harderian Gland Tumorigenesis: Low-Dose and LET 
#                            Response." Radiat Res 185(5): 449-460. (2016). 
#
#   "16Srn" = Siranart et al. "Mixed Beam Murine Harderian Gland Tumorigenesis: 
#                             Predicted Dose-Effect Relationships if neither 
#                             Synergism nor Antagonism Occurs." 
#                             Radiat Res 186(6): 577-591 (2016).  
#
#   "17Cuc" = Cucinotta & Cacao. "Non-Targeted Effects Models Predict 
#                                Significantly Higher Mars Mission Cancer Risk 
#                                than Targeted Effects Models." 
#                                Sci Rep 7(1): 1832. (2017). PMC5431989
#
#   "DER"     = Dose-effect relation(ship)"
#   "HZE"     = High atomic number Z and high energy
#   "HG"      = Harderian Gland
#   "IEA"     = Incremental Effect Additivity
#   "LET"     = Linear energy transfer; stopping power
#   "NTE"     = Non-targeted effects
#   "SEA"     = Simple Effect Additivity
#   "TE"      = Targeted effects
#   "cGy"     = Centigray

rm(list=ls()) 

# Data used is that in 16Chang plus new data becoming available in the summer of 2018 and later.
# rks_raw_data_ordered.csv includes data analyzed in .93Alp and .94Alp. Does not 
# include gamma-ray data. Includes LET=100 keV/micron for Ti, an ad-hoc compromise
# between lower value at beam entry and higher value at mouse cage.

ion_data <- data.frame(read.csv("one-ion.csv")) 
mix_data <- data.frame(read.csv("mix_ion.csv")) 
# The two .csv files contain all input HG data. Changes in the data set such as additions, corrections, or deletions should be made only in the .csv files. Such changes will then be implemented automatically by the scripts everywhere else.

# The following, which shows how to compute ion speed and the Katz amorphous track structure parameter,
# may be used for adding Cucinotta's models in 16 Chang to our scripts and comparing them to our more parsimonious models.
# GeVu is kinetic energy per atomic mass unit. An example for 670Ne20 is GeVu =10^-3*670.
# The calculations here can and will approximate Z_eff by Z, e.g. Z_eff = 10 for Ne.
#Katz = 1/round(Z^2 * (2.57 * GeVu ^2 + 4.781 * GeVu + 2.233) / (2.57 * GeVu ^2 + 4.781 * GeVu), 3) 
#special relativistic calculation of Z^2/beta^2. The numerics include conversion from GeV to joules and from u to kg.
#beta_star =Z*round(sqrt(1 / Katz), 3) #  i.e. beta = Z*sqrt(beta^2/Z^2). 
