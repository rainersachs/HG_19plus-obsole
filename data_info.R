# Copyright:    (C) 2017-2019 Sachs Undergraduate Research Apprentice Program 
#               (URAP) class at University of California, Berkeley.
#               This program and its accompanying materials are distributed 
#               under the terms of the GNU General Public License v3. As 
#               detailed in that license no warranty, explicit or implied, 
#               comes with this suite of R scripts.
# Filename:     data_info.R 
# Purpose:      Concerns radiogenic mouse Harderian gland (HG) tumorigenesis. Loads 
#               ion and tumor prevalence data from CSV files. It is part of the 
#               customized source code for the Chang 2019 HG project.
# Contact:      Rainer K. Sachs 
# Website:      https://github.com/rainersachs/mouseHG_Chang_2019plus
# Mod history:  04 Apr 2019
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
#   "NSRL"    = NASA Space Radiation Laboratory in Brookhaven NY
#   "NTE"     = Non-targeted effects
#   "SEA"     = Simple Effect Additivity
#   "TE"      = Targeted effects
#   "cGy"     = Centigray

rm(list=ls()) 

# Data used is that in 16Chang plus later NSRL data scored before 3/31/2019.
ion_data <- data.frame(read.csv("one_ion.csv")) 
mix_data <- data.frame(read.csv("mix_ion.csv"))

# Filter NA entries, determined by the existence of a dose value.
ion_data <- ion_data[!(is.na(ion_data$dose)), ]
mix_data <- mix_data[!(is.na(mix_data$dose)), ]
# The two .csv files contain all input HG data except Y_0
