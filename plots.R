# Copyright:    (C) 2017-2019 Sachs Undergraduate Research Apprentice Program
#               This program and its accompanying materials are distributed 
#               under the terms of the GNU General Public License v3.
# Filename:     plots.R 
# Purpose:      Concerns radiogenic mouse Harderian gland tumorigenesis. 
#               Contains code to generate many of the LSSR paper figures. It's
#               part of the source code for the Chang 2019 HG project.
# Contact:      Rainer K. Sachs 
# Website:      https://github.com/rainersachs/mouseHG_Chang_2019plus
# Mod history:  04 Apr 2019
# Details:      See data_info.R for further licensing, attribution, 
#               references, and abbreviation information.

source("monteCarlo.R") # Load Monte Carlo.
library(Hmisc) # Error bars.

#= shift needed due to the synergy-theory-oriented way DERs are calibrated =#
dat_down <- ion_data 
# dat_down stands for ion data shifted downward (by subtracting background) 
vvector <- ion_data[, "Prev"] # vvector is used only in this chunk
vvector <- vvector - Y_0
dat_down[, "Prev"] <- vvector

# EGH: Probably can be rewritten as follows, will check later
# dat_down[, "Prev"] <- ion_data[, "Prev"] - Y_0

# Edward: above chunk added 6/2/2019 

# ==========================================================#
#====== Fig. 1_final. Convex, Concave, Standard ============#
#=== Figs. 1&2_final are schematic. need not be redrawn ====#
# ==========================================================#
# d2 <- 0.01 * 0:200
# a <- 2; b <- .6; c <- 4
# E1 <- a * d2 + b * d2 ^ 2  # Convex
# E2 <- a * d2  # Linear no-threshold (LNT) same initial slope
# E3 <- c * (1 - (exp(- a * d2 / c))) # Concave, same initial slope
# plot(d2, E1, type = 'l', lwd = 3, bty = 'l', ann = FALSE)
# lines(d2, E2, lwd = 3)
# lines(d2, E3, lwd = 3)
# 
# a <- 0.45; b <- 1 / 8; c <- 0.8
# E1 <- b * d2 + 0.35 * b * d2 ^ 2  # Convex
# E2 <- 0.7 * a * d2  # Linear no-threshold (LNT) 
# E3 <- c * (1 - (exp(- 2 * a * d2 / c))) # Concave
# plot(d2, E3, type = 'l', lwd = 3, bty = 'l', ann = FALSE)
# lines(d2, E1, lwd = 3)
# lines(d2, E2, lwd = 3)
#============================================================#
#= Fig.3_final. Low LET Data SHIFTED, Error Bars, and DER. ==#
#============================================================#
ddose <- 0:701 
plot(c(0, 701), c(-.02, .75), pch = 19, col = 'white', ann = FALSE, bty = 'u') #RKS to RKS: force 7 ticks on y axis, see par arguments to plot()
lines(ddose, 1 - exp(- coef(summary(low_LET_model, correlation = TRUE))[1] * ddose), lwd = 2)# next is alpha particle data

errbar(dat_down[dat_down$Beam == "He", ][, "dose"],
       dat_down[dat_down$Beam == "He", ][, "Prev"], 
       yplus  = dat_down[dat_down$Beam == "He", ][, "Prev"] + dat_down[
         dat_down$Beam == "He", ][, "SD"],
       yminus = dat_down[dat_down$Beam == "He", ][, "Prev"] - dat_down[
         dat_down$Beam == "He", ][, "SD"], pch = 19, cap = 0.02, add = TRUE, 
       col = 'red', errbar.col = 'red', lwd = 1) # Alpha particle data

errbar(dat_down[dat_down$Beam == "p", ][, "dose"], 
       dat_down[dat_down$Beam == "p", ][, "Prev"],
       yplus  = dat_down[dat_down$Beam == "p", ][, "Prev"] + dat_down[
         dat_down$Beam == "p", ][, "SD"], 
       yminus = dat_down[dat_down$Beam == "p", ][, "Prev"] - dat_down[
         dat_down$Beam == "p", ][, "SD"], pch = 19, cap = 0.02, add = TRUE,
       col = 'black', errbar.col = 'black', lwd = 1) # Proton data
legend(x = "bottomright", legend = "SD", cex=0.6)
print(Y_0) # amount data points shifted down compared to DER curves etc.

##==================================================#
#== Fig.4A_final Fe DERs. points, error bars, SHIFTED==#
##==================================================#
d1_Fe <- c(0.01 * 0:9, 0.1 * 1:9, 1:80)
fe_six = calibrated_HZE_te_der(dose = d1_Fe, L = 193) #TE-only Fe600 1-ion DER
plot(c(0, 80.1), c(0, .65), col = "white", bty = 'L', ann = FALSE) # Set plot area
lines(d1_Fe, fe_six, col = 'black') # Fe TE-only 1-ion DER, no background
fe_six_nte = calibrated_HZE_nte_der(dose = d1_Fe, L = 193) # same for NTE-also
lines(d1_Fe, fe_six_nte, col = 'red', lwd =2)

errbar(dat_down[dat_down$LET==193,][,"dose"], dat_down[dat_down$LET==193,][,"Prev"],
       yplus  = dat_down[dat_down$LET == 193, ][, "Prev"] + dat_down[dat_down$LET == 193, ][, "SD"], 
       yminus = dat_down[dat_down$LET == 193, ][, "Prev"] - dat_down[dat_down$LET == 193, ][, "SD"],
       pch = 19, cap = 0.04, add = TRUE, col = 'black', errbar.col = 'black', lwd = 1)
legend(x = "topleft", legend = "Fe 600, not 95%CI, SD, Y_0=", cex=0.49)
print(Y_0)

## =============================================================#
##= Fig4_B final: Fe 1-ion NTE-also der "mixture" + both ribbons =# 
##==============================================================#
# Declare ratios and LET values for plot; Fe with LET 193 
dd0 <- c(0.01 * 0:9, 0.1 * 1:9, 1:81)

# We use the plot that neglects adjustable parameter correlations
uncorr_fig_0 <- simulate_monte_carlo(n = 500, dd0, 193, 1, vcov = FALSE)

ci_data <- data.frame(dose = dd0,
                      # Monte Carlo values
                      uncorrBottom = uncorr_fig_0$monte_carlo[1, ],
                      uncorrTop = uncorr_fig_0$monte_carlo[2, ], 
                      
                      # one-ion DER for comparison
                      fe_six = calibrated_HZE_nte_der(dose = dd0, L = 193),
                      # IEA baseline mixture DER I(d), denoted by id below
                      i = calculate_id(dd0, 193, 1, model = "NTE")[, 2])

plot(c(0, 81), c(0, .62), col = "white", bty = 'L', ann = FALSE) # Set plot area

polygon(x = c(dd0, rev(dd0)), 
        y = c(ci_data[, "uncorrTop"], rev(ci_data[, "uncorrBottom"])), 
        col = "aquamarine2", lwd = .4, border = "aquamarine2") # CI ribbon

# We use the plot that takes adjustable parameter correlations into account
corr_fig_0 <- simulate_monte_carlo(n = 500, dd0, 193, 1)
ci_data <- data.frame(dose = dd0,
                      # Monte Carlo values
                      corrBottom = corr_fig_0$monte_carlo[1, ],
                      corrTop = corr_fig_0$monte_carlo[2, ], #
                      
                      # one-ion DERs for comparison
                      fe_six = calibrated_HZE_nte_der(dose = dd0, L = 193),
                      # p = calibrated_low_LET_der(dose = dd0, L = .4),
                      
                      # IEA baseline mixture DER I(d), denoted by id below
                      i = calculate_id(dd0, 193, 1, model = "NTE")[, 2])

polygon(x = c(dd0, rev(dd0)), 
        y = c(ci_data[, "corrTop"], rev(ci_data[, "corrBottom"])),
        col = "yellow", lwd = .4, border = "orange") # CI ribbon
lines(dd0, calibrated_HZE_nte_der(dose = dd0, L = 193), col = 'red', lwd = 2)

# ========================================================== #
#= Fig. 5A_final p-Fe mix, NTE, points+error, narrow ribbon =# 
# ========================================================== #
# Declare ratios and LET values for plot
ratios <- c(3/7, 4/7) # for Fe-p 
LET_vals <- c(193, .4)
d5_new <- c(0.01 * 0:9, 0.1 * 1:9, 1:70)
d5_Fe = c(0.01 * 0:9, 0.1 * 1:9, 1:30)
d5_p <- c(0.01 * 0:9, 0.1 * 1:9, 1:40)
# We use the plot that takes adjustable parameter correlations into account
corr_fig_5A <- simulate_monte_carlo(n = 500, d5_new, LET_vals, ratios, model = "NTE")
# The first argument, n, is the number of Monte Carlo repeats. Increase for
# greater accuracy. Decrease to speed up the program.
ci_data <- data.frame(dose = d5_new,
                      # Monte Carlo values
                      corrBottom = corr_fig_5A$monte_carlo[1, ],
                      corrTop = corr_fig_5A$monte_carlo[2, ], #
                      
                      # one-ion DERs 
                      fe_six = calibrated_HZE_nte_der(dose = d5_new, L = 193),
                      p = calibrated_low_LET_der(dose = d5_new, L = .4),
                      
                      # IEA baseline mixture DER I(d), denoted by id below
                      i = calculate_id(d5_new, LET_vals, ratios, model = "NTE")[, 2])

plot(c(0, 71), c(0,.62), col = "white", bty = 'L', ann = FALSE) # Set plot area
polygon(x = c(d5_new, rev(d5_new)), y = c(ci_data[, "corrTop"], rev(ci_data[, "corrBottom"])),
        xpd = -1, col = "yellow", lwd = .4, border = "orange") # Narrow CI ribbon
fe_six_one = calibrated_HZE_nte_der(dose = d5_Fe, L = 193)
p_one = calibrated_low_LET_der(dose = d5_p, L = .4)
lines(d5_p, p_one, col = 'brown', lwd = 2) # p DER
lines(d5_Fe, fe_six_one, col = 'blue') # Fe NTE-also 1-ion DER
lines(ci_data[, "dose"], ci_data[, "i"], col = 'red', lwd = 1) # I(d)

# The following chunk is just asking for trouble with unamed vectors
# Edward please fix
Prev=as.numeric(mix_data$Prev)-Y_0; SD=as.numeric(mix_data$SD)
errbar(70,Prev[3], yplus  = Prev[3] +  SD[3],
       yminus = Prev[3] -SD[3],
       pch = 19, cap = 0.05, add = TRUE, col = 'black', errbar.col = 'black', lwd = 1)
#legend(x = "topleft", legend = "not 95%CI, SD", cex=0.3)

# ========================================================== #
#= Fig. 5B_final p-Fe mix, TE, point & error, narrow ribbon =# 
# ========================================================== #
# We use the plot that takes adjustable parameter correlations into account
corr_fig_5B <- simulate_monte_carlo(n = 500, d5_new, LET_vals, ratios, model = "TE")
# The first argument, n, is the number of Monte Carlo repeats. Increase for
# greater accuracy. Decrease to speed up the program.
ci_data <- data.frame(dose = d5_new,
                      # Monte Carlo values
                      corrBottom = corr_fig_5B$monte_carlo[1, ],
                      corrTop = corr_fig_5B$monte_carlo[2, ], #
                      
                      # one-ion DERs 
                      fe_six = calibrated_HZE_te_der(dose = d5_new, L = 193),
                      p = calibrated_low_LET_der(dose = d5_new, L = .4),
                      
                      # IEA baseline mixture DER I(d), denoted by id below
                      i = calculate_id(d5_new, LET_vals, ratios, model = "TE")[, 2])

plot(c(0, 71), c(0,.62), col = "white", bty = 'L', ann = FALSE) # Set plot area
polygon(x = c(d5_new, rev(d5_new)), y = c(ci_data[, "corrTop"], rev(ci_data[, "corrBottom"])),
        xpd = -1, col = "yellow", lwd = .4, border = "orange") # Narrow CI ribbon
fe_six_te_one = calibrated_HZE_te_der(dose = d5_Fe, L = 193)
p_one = calibrated_low_LET_der(dose = d5_p, L = .4)
lines(d5_p, p_one, col = 'brown', lwd = 2) # p DER
lines(d5_Fe, fe_six_te_one, col = 'blue') # Fe TE-only 1-ion DER
lines(ci_data[, "dose"], ci_data[, "i"], col = 'red', lwd = 1) # I(d)


Prev=as.numeric(mix_data$Prev)-Y_0; SD=as.numeric(mix_data$SD)
errbar(70,Prev[3], yplus  = Prev[3] +  SD[3],
       yminus = Prev[3] -SD[3],
       pch = 19, cap = 0.05, add = TRUE, col = 'black', errbar.col = 'black', lwd = 1)


# ================================================================ #
#=== Figs. 6&7_A&B_final. 2 other mixtures. NTE-also, TE-only. ====#
#== Export as images width 400, height 350, font arial. then jpg ==#
#= Fig.6A_final. Fe-Si 50-50 total 40 cGy, NTE-also. like Fig. 5A =#
#==================================================================#
# Declare ratios and LET values for plot
ratios <- c(1/2, 1/2)
LET_vals <- c(193, 70)
d6 <- c(0.01 * 0:9, 0.1 * 1:9, 1:40)
d6_ion <- c(0.01 * 0:9, 0.1 * 1:9, 1:20)
# We use the plot that takes adjustable parameter correlations into account
corr_fig_6A <- simulate_monte_carlo(n = 500, d6, LET_vals, ratios, model = "NTE")
# The first argument, n, is the number of Monte Carlo repeats. Increase for
# greater accuracy. Decrease to speed up the program.
ci_data <- data.frame(dose = d6,
                      # Monte Carlo values
                      corrBottom = corr_fig_6A$monte_carlo[1, ],
                      corrTop = corr_fig_6A$monte_carlo[2, ], #
                      
                      # one-ion DERs for comparison
                      fe_six = calibrated_HZE_nte_der(dose = d6, L = 193),
                      si = calibrated_HZE_nte_der(dose = d6, L = 70),
                      # does d6_ion above rather than d_6 work?
                      # IEA baseline mixture DER I(d), denoted by id below
                      i = calculate_id(d6, LET_vals, ratios, model = "NTE")[, 2])

#We make the ribbon plot for correlated parameters
plot(c(0, 41), c(0, .36), col = "white", bty = 'L', ann = FALSE) # Set plot area
polygon(x = c(d6, rev(d6)), y = c(ci_data[, "corrTop"], rev(ci_data[, "corrBottom"])),
        xpd = -1, col = "yellow", lwd = .4, border = "orange") # Narrow CI ribbon

fe_six_one = calibrated_HZE_nte_der(dose = d6_ion, L = 193)
si_one = calibrated_HZE_nte_der(dose = d6_ion, L = 70)
lines(d6_ion, si_one, col = 'brown', lwd=2)
lines(d6_ion, fe_six_one, col = 'blue', lwd=2) # Fe TE-only 1-ion DER
lines(ci_data[, "dose"], ci_data[, "i"], col = 'red', lwd = 2) # IEA NSNA I(d)

Prev=as.numeric(mix_data$Prev)-Y_0; SD=as.numeric(mix_data$SD)
errbar(40,Prev[5], yplus  = Prev[5] +  SD[5],
       yminus = Prev[5] -SD[5],
       pch = 19, cap = 0.05, add = TRUE, col = 'black', errbar.col = 'black', lwd = 1)

#=========== 6B_final. repeat Fig. 6A_final but for TE-only DERs ============#
# We use the plot that takes adjustable parameter correlations into account
corr_fig_6B <- simulate_monte_carlo(n = 500, d6, LET_vals, ratios, model = "TE")
# The first argument, n, is the number of Monte Carlo repeats. Increase for
# greater accuracy. Decrease to speed up the program.
ci_data <- data.frame(dose = d6,
                      # Monte Carlo values
                      corrBottom = corr_fig_6B$monte_carlo[1, ],
                      corrTop = corr_fig_6B$monte_carlo[2, ], #
                      
                      # one-ion DERs for comparison
                      fe_six = calibrated_HZE_te_der(dose = d6, L = 193),
                      si = calibrated_HZE_te_der(dose = d6, L = 70),
                      # does d6_ion above rather than d_6 work?
                      # IEA baseline mixture DER I(d), denoted by id below
                      i = calculate_id(d6, LET_vals, ratios, model = "TE")[, 2])

#We make the ribbon plot for correlated parameters
plot(c(0, 41), c(0, .36), col = "white", bty = 'L', ann = FALSE) # Set plot area
polygon(x = c(d6, rev(d6)), y = c(ci_data[, "corrTop"], rev(ci_data[, "corrBottom"])),
        xpd = -1, col = "yellow", lwd = .4, border = "orange") # Narrow CI ribbon

fe_six_one = calibrated_HZE_te_der(dose = d6_ion, L = 193)
si_one = calibrated_HZE_te_der(dose = d6_ion, L = 70)
lines(d6_ion, si_one, col = 'brown', lwd=2)
lines(d6_ion, fe_six_one, col = 'blue', lwd=2) # Fe TE-only 1-ion DER
lines(ci_data[, "dose"], ci_data[, "i"], col = 'red', lwd = 2) # IEA NSNA I(d)

Prev=as.numeric(mix_data$Prev)-Y_0; SD=as.numeric(mix_data$SD)
errbar(40,Prev[6], yplus  = Prev[6] +  SD[6],
       yminus = Prev[6] -SD[6],
       pch = 19, cap = 0.05, add = TRUE, col = 'black', errbar.col = 'black', lwd = 1)

# ========================================================== #
#= Fig. 7A_final p-Si mix, NTE, points+error, narrow ribbon =# 
# ========================================================== #
# Declare ratios and LET values for plot
ratios <- c(.6, .4) # for Fe-p 
LET_vals <- c(.4, 70)
d7_new <- c(0.01 * 0:9, 0.1 * 1:9, 1:100)
d7_Si = c(0.01 * 0:9, 0.1 * 1:9, 1:40)
d7_p <- c(0.01 * 0:9, 0.1 * 1:9, 1:60)
# We use the plot that takes adjustable parameter correlations into account
corr_fig_7A <- simulate_monte_carlo(n = 500, d7_new, LET_vals, ratios, model = "NTE")
# The first argument, n, is the number of Monte Carlo repeats. Increase for
# greater accuracy. Decrease to speed up the program.
ci_data <- data.frame(dose = d7_new,
                      # Monte Carlo values
                      corrBottom = corr_fig_7A$monte_carlo[1, ],
                      corrTop = corr_fig_7A$monte_carlo[2, ], #
                      
                      # one-ion DERs 
                      Si = calibrated_HZE_nte_der(dose = d7_new, L = 70),
                      p = calibrated_low_LET_der(dose = d7_new, L = .4),
                      
                      # IEA baseline mixture DER I(d), denoted by id below
                      i = calculate_id(d7_new, LET_vals, ratios, model = "NTE")[, 2])

plot(c(0, 101), c(0,.36), col = "white", bty = 'L', ann = FALSE) # Set plot area
polygon(x = c(d7_new, rev(d7_new)), y = c(ci_data[, "corrTop"], rev(ci_data[, "corrBottom"])),
        xpd = -1, col = "yellow", lwd = .4, border = "orange") # Narrow CI ribbon
si_one = calibrated_HZE_nte_der(dose = d7_Si, L = 70)
p_one = calibrated_low_LET_der(dose = d7_p, L = .4)
lines(d7_p, p_one, col = 'brown', lwd = 2) # p DER
lines(d7_Si, si_one, col = 'blue') # Si NTE-also 1-ion DER
lines(ci_data[, "dose"], ci_data[, "i"], col = 'red', lwd = 1) # I(d)

# The following chunk is just asking for trouble with unamed vectors
# Edward please fix
Prev=as.numeric(mix_data$Prev)-Y_0; SD=as.numeric(mix_data$SD)
errbar(100,Prev[1], yplus  = Prev[1] +  SD[1],
       yminus = Prev[1] -SD[1],
       pch = 19, cap = 0.05, add = TRUE, col = 'black', errbar.col = 'black', lwd = 1)
#legend(x = "topleft", legend = "not 95%CI, SD", cex=0.3)

# ========================================================== #
#= Fig. 7B_final p-Si mix, TE, point & error, narrow ribbon =# 
# ========================================================== #
# We use the plot that takes adjustable parameter correlations into account
corr_fig_7B <- simulate_monte_carlo(n = 500, d7_new, LET_vals, ratios, model = "TE")
# The first argument, n, is the number of Monte Carlo repeats. Increase for
# greater accuracy. Decrease to speed up the program.
ci_data <- data.frame(dose = d7_new,
                      # Monte Carlo values
                      corrBottom = corr_fig_7B$monte_carlo[1, ],
                      corrTop = corr_fig_7B$monte_carlo[2, ], #
                      
                      # one-ion DERs 
                      si = calibrated_HZE_te_der(dose = d7_new, L = 70),
                      p = calibrated_low_LET_der(dose = d7_new, L = .4),
                      
                      # IEA baseline mixture DER I(d), denoted by id below
                      i = calculate_id(d7_new, LET_vals, ratios, model = "TE")[, 2])

plot(c(0, 100), c(0,.36), col = "white", bty = 'L', ann = FALSE) # Set plot area
polygon(x = c(d7_new, rev(d7_new)), y = c(ci_data[, "corrTop"], rev(ci_data[, "corrBottom"])),
        xpd = -1, col = "yellow", lwd = .4, border = "orange") # Narrow CI ribbon
si_te_one = calibrated_HZE_te_der(dose = d7_Si, L = 70)
p_one = calibrated_low_LET_der(dose = d7_p, L = .4)
lines(d7_p, p_one, col = 'brown', lwd = 2) # p DER
lines(d7_Si, si_te_one, col = 'blue') # Si TE-only 1-ion DER
lines(ci_data[, "dose"], ci_data[, "i"], col = 'red', lwd = 1) # I(d)


Prev=as.numeric(mix_data$Prev)-Y_0; SD=as.numeric(mix_data$SD)
errbar(100,Prev[1], yplus  = Prev[1] +  SD[1],
       yminus = Prev[1] -SD[1],
       pch = 19, cap = 0.05, add = TRUE, col = 'black', errbar.col = 'black', lwd = 1)

#==================================================================#
#== Fig8_final. 3-ion mix. Correlated vs Uncorr CI Overlay Plot ===#
#==================================================================#
# Consists of 3 HZE ions in our 5/20/2019 data set.
# Declare ratios and LET values for plot
ratios <- rep(1/3,3)
LET_vals <- c(70, 100, 193)
d8 <- c(0.1 * 0:9, 1:60)
d8t <- c(0.1 * 0:9, 1:20)
# We begin with the correlated plot
corr_fig_8 <- simulate_monte_carlo(n = 500, d8, LET_vals, ratios, model = "NTE")
# Comments for Fig. 10 apply with minor changes here and in some other lines
# We now calculate the uncorrelated Monte Carlo
uncorr_fig_8 <- simulate_monte_carlo(n = 500, d8, LET_vals, ratios, model = "NTE", vcov = FALSE)

ci_data <- data.frame(dose = d8,
                      # Monte Carlo values
                      corrBottom = corr_fig_8$monte_carlo[1, ],
                      corrTop = corr_fig_8$monte_carlo[2, ],
                      uncorrBottom = uncorr_fig_8$monte_carlo[1, ],
                      uncorrTop = uncorr_fig_8$monte_carlo[2, ],
                      
                      # DER values
                      si = calibrated_HZE_nte_der(dose = d8, L = 70),
                      ti = calibrated_HZE_nte_der(dose = d8, L = 100),
                      fe_six = calibrated_HZE_nte_der(dose = d8, L = 193),
                      i = calculate_id(d8, LET_vals, ratios, model = "NTE")[, 2])

# Plotting call. 
plot(c(0, 60), c(0, .450), col = "white", bty = 'l', ann = FALSE) # Next is broad CI ribbon
polygon(x = c(d8, rev(d8)), y = c(ci_data[, "uncorrTop"], rev(ci_data[, "uncorrBottom"])),
        xpd = -1, col = "aquamarine2", lwd = .5, border = "aquamarine2") # Wide CI

polygon(x = c(d8, rev(d8)), y = c(ci_data[, "corrTop"], rev(ci_data[, "corrBottom"])),
        xpd = -1, col = "yellow", border = "orange", lwd = .2) # Narrow CI

lines(d8t, calibrated_HZE_nte_der(dose = d8t, L = 70), col = 'black', lwd = 2)
lines(d8t, calibrated_HZE_nte_der(dose = d8t, L = 100), col = 'black', lwd = 2)
lines(d8t, calibrated_HZE_nte_der(dose = d8t, L = 193), col = 'black', lwd = 2)
lines(ci_data[, "dose"], ci_data[, "i"], col = 'red', lwd = 3,lty=3) # I(d)
#errbar(49,.35, yplus  = .35 +  .05,
#      yminus = .35 -.05,
#     pch = 19, cap = 0.05, add = TRUE, col = 'black', errbar.col = 'black', lwd = 1) 

#==================================================================#
#== Fig9_final. 8-ion mix. Correlated vs Uncorr CI Overlay Plot ===#
#==================================================================#
# Consists of all 8 HZE ions in our 5/20/2019 data set.
# Declare ratios and LET values for plot
ratios <- rep(1/8,8)
LET_vals <- c(20, 25, 70, 100, 193, 250, 464, 953)
d9 <- c(.001*0:9,0.1 * 1:9, 1:50)
d9t <- c(.001*0:9,0.1 * 1:9, 1:6,6.1,6.2,6.25)
d9B = c(.001*0:9,0.1 * 1:9, 1:10)
# We begin with the correlated plot
corr_fig_9 <- simulate_monte_carlo(n = 500, d9, LET_vals, ratios, model = "NTE")
# We now calculate the uncorrelated Monte Carlo
uncorr_fig_9 <- simulate_monte_carlo(n = 500, d9, LET_vals, ratios, model = "NTE", vcov = FALSE)

ci_data <- data.frame(dose = d9,
                      # Monte Carlo values
                      corrBottom = corr_fig_9$monte_carlo[1, ],
                      corrTop = corr_fig_9$monte_carlo[2, ],
                      uncorrBottom = uncorr_fig_9$monte_carlo[1, ],
                      uncorrTop = uncorr_fig_9$monte_carlo[2, ],
                      
                      # DER values
                      o = calibrated_HZE_nte_der(dose = d9, L = 20),
                      ne = calibrated_HZE_nte_der(dose = d9, L = 25),
                      si = calibrated_HZE_nte_der(dose = d9, L = 70),
                      ti = calibrated_HZE_nte_der(dose = d9, L = 100),
                      fe_six = calibrated_HZE_nte_der(dose = d9, L = 193),
                      fe_three = calibrated_HZE_nte_der(dose = d9, L = 250),
                      nb = calibrated_HZE_nte_der(dose = d9, L = 464),
                      la = calibrated_HZE_nte_der(dose = d9, L = 953),
                      i = calculate_id(d9, LET_vals, ratios, model = "NTE")[, 2])

# Plotting call. 
plot(c(0, 50), c(0, .40), col = "white", bty = 'l', ann = FALSE) # Next is broad CI ribbon
polygon(x = c(d9, rev(d9)), y = c(ci_data[, "uncorrTop"], rev(ci_data[, "uncorrBottom"])),
        xpd = -1, col = "aquamarine2", lwd = .5, border = "aquamarine2") # Wide CI

polygon(x = c(d9, rev(d9)), y = c(ci_data[, "corrTop"], rev(ci_data[, "corrBottom"])),
        xpd = -1, col = "yellow", border = "orange", lwd = .2) # Narrow CI

lines(d9t, calibrated_HZE_nte_der(dose = d9t, L = 20), col = 'blue', lwd = 1)
lines(d9t, calibrated_HZE_nte_der(dose = d9t, L = 25), col = 'blue', lwd = 1)
lines(d9t, calibrated_HZE_nte_der(dose = d9t, L = 70), col = 'blue', lwd = 1)
lines(d9t, calibrated_HZE_nte_der(dose = d9t, L = 100), col = 'blue', lwd = 1)
lines(d9t, calibrated_HZE_nte_der(dose = d9t, L = 193), col = 'blue', lwd = 1)
lines(d9t, calibrated_HZE_nte_der(dose = d9t, L = 250), col = 'blue', lwd = 1)
lines(d9t, calibrated_HZE_nte_der(dose = d9t, L = 464), col = 'blue', lwd = 1)
lines(d9t, calibrated_HZE_nte_der(dose = d9t, L = 953), col = 'blue', lwd = 1)
lines(ci_data[, "dose"], ci_data[, "i"], col = 'red', lwd = 2) # I(d)
errbar(49,.35, yplus  = .35 +  .05,
       yminus = .35 -.05,
       pch = 19, cap = 0.05, add = TRUE, col = 'black', errbar.col = 'black', lwd = 1) #OK up to here
#=====================================#
# panel B ============================#
#=====================================#
ratios <- rep(1/8,8)
LET_vals <- c(20, 25, 70, 100, 193, 250, 464, 953)
d9 <- c(.001*0:9,0.1 * 1:9, 1:50)
d9t <- c(.001*0:9,0.1 * 1:9, 1:6,6.1,6.2,6.25)
d9B = c(.001*0:9,0.1 * 1:9, 1:10)
# We begin with the correlated plot
corr_fig_9B <- simulate_monte_carlo(n = 500, d9B, LET_vals, ratios, model = "NTE")
# We now calculate the uncorrelated Monte Carlo
uncorr_fig_9B <- simulate_monte_carlo(n = 500, d9B, LET_vals, ratios, model = "NTE", vcov = FALSE)

ci_data <- data.frame(dose = d9B,
                      # Monte Carlo values
                      corrBottom = corr_fig_9B$monte_carlo[1, ],
                      corrTop = corr_fig_9B$monte_carlo[2, ],
                      uncorrBottom = uncorr_fig_9B$monte_carlo[1, ],
                      uncorrTop = uncorr_fig_9B$monte_carlo[2, ],
                      
                      # DER values
                      o = calibrated_HZE_nte_der(dose = d9B, L = 20),
                      ne = calibrated_HZE_nte_der(dose = d9B, L = 25),
                      si = calibrated_HZE_nte_der(dose = d9B, L = 70),
                      ti = calibrated_HZE_nte_der(dose = d9B, L = 100),
                      fe_six = calibrated_HZE_nte_der(dose = d9B, L = 193),
                      fe_three = calibrated_HZE_nte_der(dose = d9B, L = 250),
                      nb = calibrated_HZE_nte_der(dose = d9B, L = 464),
                      la = calibrated_HZE_nte_der(dose = d9B, L = 953),
                      i = calculate_id(d9B, LET_vals, ratios, model = "NTE")[, 2])

# Plotting call. 
plot(c(0, 10), c(0, .12), col = "white", bty = 'l', ann = FALSE) # Next is broad CI ribbon
polygon(x = c(d9B, rev(d9B)), y = c(ci_data[, "uncorrTop"], rev(ci_data[, "uncorrBottom"])),
        xpd = -1, col = "aquamarine2", lwd = .5, border = "aquamarine2") # Wide CI

polygon(x = c(d9B, rev(d9B)), y = c(ci_data[, "corrTop"], rev(ci_data[, "corrBottom"])),
        xpd = -1, col = "yellow", border = "orange", lwd = .2) # Narrow CI

lines(d9t, calibrated_HZE_nte_der(dose = d9t, L = 20), col = 'blue', lwd = 1)
lines(d9t, calibrated_HZE_nte_der(dose = d9t, L = 25), col = 'blue', lwd = 1)
lines(d9t, calibrated_HZE_nte_der(dose = d9t, L = 70), col = 'blue', lwd = 1)
lines(d9t, calibrated_HZE_nte_der(dose = d9t, L = 100), col = 'blue', lwd = 1)
lines(d9t, calibrated_HZE_nte_der(dose = d9t, L = 193), col = 'blue', lwd = 1)
lines(d9t, calibrated_HZE_nte_der(dose = d9t, L = 250), col = 'blue', lwd = 1)
lines(d9t, calibrated_HZE_nte_der(dose = d9t, L = 464), col = 'blue', lwd = 1)
lines(d9t, calibrated_HZE_nte_der(dose = d9t, L = 953), col = 'blue', lwd = 1)
lines(ci_data[, "dose"], ci_data[, "i"], col = 'red', lwd = 2, lty = 2) # I(d)

#errbar(49,.35, yplus  = .35 +  .05,
#yminus = .35 -.05,
#pch = 19, cap = 0.05, add = TRUE, col = 'black', errbar.col = 'black', lwd = 1)

# #==========================================================#
# #===== ALL the rest is probably redundant or obsolete =====#
# #========== But some might be useful for WebSup ===========#
# #======= So RKS merely commented it out  5/16/2019 ========#
# #==========================================================#
# 
# # =========================================================#
# #= Fe 1-ion NTE-also der "mixture" 
# #==========================================================#
# # Declare ratios and LET values for plot; Fe with LET 193 only
# dd0 <- c(0.01 * 0:9, 0.1 * 1:9, 1:81)
# 
# # We use the plot that neglects adjustable parameter correlations
# uncorr_fig_0 <- simulate_monte_carlo(n = 500, dd0, 193, 1, vcov = FALSE)
# 
# ci_data <- data.frame(dose = dd0,
#                       # Monte Carlo values
#                       uncorrBottom = uncorr_fig_0$monte_carlo[1, ],
#                       uncorrTop = uncorr_fig_0$monte_carlo[2, ], 
#                       
#                       # one-ion DERs for comparison
#                       fe_six = calibrated_HZE_nte_der(dose = dd0, L = 193),
#                       # p = calibrated_low_LET_der(dose = dd0, L = .4),
#                       
#                       # IEA baseline mixture DER I(d), denoted by id below
#                       i = calculate_id(dd0, 193, 1, model = "NTE")[, 2])
# 
# plot(c(0, 81), c(0, .62), col = "white", bty = 'L', ann = FALSE) # Set plot area
# 
# polygon(x = c(dd0, rev(dd0)), 
#         y = c(ci_data[, "uncorrTop"], rev(ci_data[, "uncorrBottom"])), 
#         col = "aquamarine2", lwd = .4, border = "aquamarine2") # CI ribbon
# 
# # We use the plot that takes adjustable parameter correlations into account
# corr_fig_0 <- simulate_monte_carlo(n = 500, dd0, 193, 1)
# ci_data <- data.frame(dose = dd0,
#                       # Monte Carlo values
#                       corrBottom = corr_fig_0$monte_carlo[1, ],
#                       corrTop = corr_fig_0$monte_carlo[2, ], #
#                       
#                       # one-ion DERs for comparison
#                       fe_six = calibrated_HZE_nte_der(dose = dd0, L = 193),
#                       # p = calibrated_low_LET_der(dose = dd0, L = .4),
#                       
#                       # IEA baseline mixture DER I(d), denoted by id below
#                       i = calculate_id(dd0, 193, 1, model = "NTE")[, 2])
# 
# polygon(x = c(dd0, rev(dd0)), 
#         y = c(ci_data[, "corrTop"], rev(ci_data[, "corrBottom"])),
#         col = "yellow", lwd = .4, border = "orange") # CI ribbon
# lines(dd0, calibrated_HZE_nte_der(dose = dd0, L = 193), col = 'red', lwd = 2)
# 
# observed <- select(filter(Fe_600, dose < 90), c(dose, Prev, SD))
# 
# errbar(observed[["dose"]], observed[["Prev"]], 
#        yplus = observed[["Prev"]] + observed[["SD"]],
#        yminus = observed[["Prev"]] - observed[["SD"]],
#        pch = 20, cap = 0.03, add = TRUE, col = 'black', errbar.col = 'black', lwd = 1)
# 
# legend(x = "topleft", legend = "Fe600", cex=0.51)
# 
# 
# #===============================================================================#
# #== Fig. 1REBP. SEA Synergy Theory. Not needed for LSSR paper but maybe WebSup==#
# #===============================================================================#
# # d1 <- 0.01 * 0:100
# # d <- 0.5 * d1
# # E1 <- d1 ^ 2
# # E2 <- 2 * d1 ^ 2
# # SEA <- d ^ 2 + 2 * d ^ 2
# # plot(d1, E2, type = 'l', lwd = 3, bty = 'l', ann = FALSE)
# # lines(d1, E1, lwd = 3)
# # lines(d1, SEA, lwd = 3, lty = 2)
# 
# #==================================================================#
# #====== Fig. 3. Shape of DER for Fe 600 MeV/u. probably redundant==#
# #==================================================================#
# d3A <- c(0.01 * 0:9, 0.1 * 1:9, 1:150) 
# prevalence <- calibrated_HZE_nte_der(d3A, 193) 
# plot(d3A, prevalence, type = 'l', bty = 'u', lwd = 3, ann = FALSE)
# # legend(x = "bottomright", 
# #        legend = "dose in centiGy; Fe 193 zoom in twice",
# #        cex = 0.4, inset = 0.025)
# 
# d3B <- 2 * 10 ^ -5 * 0:1600 # Zoom in by a factor of 10^4
# prevalence <- calibrated_HZE_nte_der(d3B, 193)
# plot(d3B, prevalence, type = 'l', bty = 'u', lwd = 3, ann = FALSE)
# 
# d3C <- 10 ^ -6 * 0:1600 # Zoom in by another factor of 20
# prevalence <- calibrated_HZE_nte_der(d3C, 193)
# plot(d3C, prevalence, type = 'l', bty = 'u', lwd = 3)
# 
# #=============== Sample Fig. for visual check of models ===================#
# dvchk <- 0.2*0:200 
# plot(c(-.01, 40), c(-.02, .4), pch = 19, col = 'white', ann = FALSE, bty = 'u')
# LET = 193
# lines(dvchk, Y_0 + calibrated_HZE_nte_der(dvchk, LET))
# lines(dvchk, Y_0 + calibrated_HZE_te_der(dvchk, LET), col= 'purple')
# # mdf %>% filter(Z < 3 & d>0) # shows .R grammar for subsetting data base
# visual_data <- ion_data %>% filter (LET == 193) # needs number 193 not name??
# errbar(visual_data[, "dose"], visual_data[, "Prev"], 
#        yplus  = visual_data[, "Prev"] +  visual_data[, "SD"],
#        yminus = visual_data[, "Prev"] -  visual_data[, "SD"], pch = 19,
#        cap = 0.02, add = TRUE, col = 'blue', errbar.col = 'blue', lwd = 2)
# legend(x = "topleft", legend = "LET = 193",cex=0.6) 
# #============= This plot runs. Looking at the low doses suggests that perhaps if we
# # confined attention to data for doses in the interval (0, 40] cGy nte4 would do much better. 
# # Also, because the line defining visual_data uses filter which only accepts a numerical
# # entry for LET I don't see how to get corresponding plots for all subdata sets of interest
# # without a huge hassle. For both these reasons this plot chunk needs work. ============#
# 
# # =============================================================== #
# #=========== OMIT? Fig. 3B_final SLI der ribbon,  is CI? =========# 
# #=================================================================#
# # Declare ratios and LET values for plot; SLI
# dd0 <- 0:701
# # Ribbon plot for this 1-parameter model
# uncorr_fig_0 <- simulate_monte_carlo(n = 500, dd0, .4, 1, vcov = FALSE)
# 
# ci_data <- data.frame(dose = dd0,
#                       # Monte Carlo values
#                       uncorrBottom = uncorr_fig_0$monte_carlo[1, ],
#                       uncorrTop = uncorr_fig_0$monte_carlo[2, ], 
#                       
#                       # one-ion DERs for comparison
#                       # fe_six = calibrated_HZE_nte_der(dose = dd0, L = .4),
#                       #p = calibrated_low_LET_der(dose = dd0, L = .4),
#                       
#                       # IEA baseline mixture DER I(d), denoted by id below
#                       i = calculate_id(dd0, .4, 1, model = "NTE")[, 2])
# # Edward: why does this chunk work even though the above line says "NTE"?
# 
# plot(c(0, 701), c(0, .7), col = "white", bty = 'L', ann = FALSE) # Set plot area
# 
# polygon(x = c(dd0, rev(dd0)), 
#         y = c(ci_data[, "uncorrTop"], rev(ci_data[, "uncorrBottom"])), 
#         col = "yellow", lwd = .4, border = "orange") # CI ribbon
# lines(ddose, 1 - exp(- coef(summary(low_LET_model, correlation = TRUE))[1] * ddose), lwd = 2, col='red')
# errbar(ion_data[ion_data$Beam == "He", ][, "dose"], ion_data[ion_data$Beam == "He", ][, "Prev"], 
#        yplus  = ion_data[ion_data$Beam == "He", ][, "Prev"] + ion_data[ion_data$Beam == "He", ][, "SD"],
#        yminus = ion_data[ion_data$Beam == "He", ][, "Prev"] - ion_data[ion_data$Beam == "He", ][, "SD"],
#        pch = 20, cap = 0.02, add = TRUE, col = 'black', errbar.col = 'black', lwd = 1) # Alpha particle data
# 
# errbar(ion_data[ion_data$Beam == "p", ][, "dose"], ion_data[ion_data$Beam == "p", ][, "Prev"],
#        yplus  = ion_data[ion_data$Beam == "p", ][, "Prev"] + ion_data[ion_data$Beam == "p", ][, "SD"], 
#        yminus = ion_data[ion_data$Beam == "p", ][, "Prev"] - ion_data[ion_data$Beam == "p", ][, "SD"], 
#        pch = 20, cap = 0.02, add = TRUE, col = 'black', errbar.col = 'black', lwd = 1) # Proton data
# 
# #========================================================#
# #=============Fig. 6. One HZE Ion, One Low-LET ===========#
# #========================================================#
# # We will always use NTE & TE for HZE in the Chang 2019 paper. 
# d <- 1 * 0:302
# r <- c(.2, .8) # Dose proportions.
# plot(x = d, y = calibrated_HZE_nte_der(dose = d, L = 193), type = "l", #  Now plot DERs 
#      xlab = "dose", ylab = "HG", bty = 'u', col = 'green', lwd = 2) # Fe DER
# lines(x = d, y = calibrated_low_LET_der(d, 0), col = 'orange', lwd = 2) # Low LET DER 
# lines(x = d, y = calculate_id(d, c(193, 0), r)[, 2], col = "red", lwd = 3) # I(d)
# lines(x = d, y = calibrated_HZE_nte_der(dose = .2 * d, L = 193) + 
#         calibrated_low_LET_der(.8 * d, 0), lty = 2, lwd = 2) # SEA mixture baseline S(d)
# 
# r <- c(.8, .2) # Panel B, proportions reversed with low LET small  
# plot(x = d, y = calibrated_HZE_nte_der(dose = d, L = 193), type = "l", 
#      xlab = "dose", ylab = "HG", bty = 'u', col = 'green', lwd = 2) 
# lines(x = d, y = calibrated_low_LET_der(d, 0), col = 'orange', lwd = 2) 
# lines(x = d, y = calculate_id(d, c(193, 0), r)[, 2], col = "red", lwd = 3) # I(d)
# lines(x = d, y = calibrated_HZE_nte_der(dose = .8 * d, L = 193) +
#         calibrated_low_LET_der(.2 *d, 0), lty= 2, lwd=2) 
# 
# #===================================================#
# #======== Fig. 7. 80% LowLET and Four HZE Ions =====#
# #===============================================#
# d7 <- c(0.01*0:9, 0.1*1:9, 1:41)
# plot(c(0,41),c(0,0.35), col="white", xlab = "Dose (cGy)", ylab = "HG", bty = 'l')
# lines(x = d7, y = calculate_SEA(d7, c(0.4, 40, 110, 180, 250), c(.8, rep(.05,4)),
#                                 lowLET = TRUE), col = "black", lwd = 2, lty = 2)
# lines(x = d7, y = calibrated_low_LET_der(dose = d7, L = 0.4),
#       col = "orange", lwd = 2) # low LET DER
# lines(x = d7, y = calibrated_HZE_nte_der(dose = d7, L = 40),
#       col = "green", lwd = 2) # HZE DER
# lines(x = d7, y = calibrated_HZE_nte_der(dose = d7, L = 110),
#       col = "purple", lwd = 2)
# lines(x = d7, y = calibrated_HZE_nte_der(dose = d7, L = 180),
#       col = "blue", lwd = 2)
# lines(x = d7, y = calibrated_HZE_nte_der(dose = d7, L = 250),
#       col = "aquamarine2", lwd = 2)
# lines(x = d7, y = calculate_id(d7, c(0.4, 40, 110, 180, 250),
#                                c(.8, rep(.05, 4)), model = "NTE")[, 2], col = "red", lwd = 3) # I(d)
# legend(x = "topleft", legend = c("Low-LET","L=40", "L=110", "L=180", 
#                                  "L=250", "I(d)", "S(d)"),
#        col = c("orange", "green","purple","blue", "aquamarine2", "red", "black"), 
#        lwd = c(2, 2, 2, 2, 2, 3, 2), 
#        lty = c(1, 1, 1, 1,  1, 1, 2), cex = 0.3, inset = 0.025)
# 
# #==================================================================================#
# #=== Fe-Si mix, IEA, SEA 06/2018 mix data point probably obsolete even for WebSup==#
# #==================================================================================#
# # d8 <- c(0.01 * 0:9, 0.1 * 1:9, 1:41)
# # plot(c(0, 40), c(0, 0.4), col = "white", bty = 'l', xlab = "Dose (cGy)", ylab = "HG", ann=F)
# # #lines(x = d8, y = calibrated_HZE_nte_der(dose = d8, L = 70), col = "cyan", lwd = 2)
# # #lines(x = d8, y = calibrated_HZE_nte_der(dose = d8, L = 193), col = "orange", lwd = 2)
# # lines(x = d8/2, y = calibrated_HZE_nte_der(dose = d8/2, L = 70), col = "cyan", lwd = 2) # for HRS 2019
# # lines(x = d8/2, y = calibrated_HZE_nte_der(dose = d8/2, L = 193), col = "orange", lwd = 2)
# # lines(x = d8, y = calculate_id(d8, c(70, 193), c(0.5 , 0.5), model = "NTE")[, 2],
# #       col = "red", lwd = 3) # I(d)
# # lines(x = d8, y = calculate_SEA(d8, c(70, 193), c(1/2, 1/2), n = 2), col = "black",
# #       lwd = 2, lty = 2)
# # # RKS: above is baselines; below is new unpublished mix data June 2018
# # # new works but numbers like "40" and "1" will need to be replaced by names eventually
# # errbar(40, mix_data[6,"Prev"],
# #        yplus  = mix_data[6,"Prev"] + mix_data[6,"SD"] , 
# #        yminus = mix_data[6,"Prev"] -  mix_data[6,"SD"], 
# #        pch = 19, cap = 0.02, add = TRUE, col = 'black', errbar.col = 'black', lwd = 2)
# # 
# # legend(x = "topleft", legend = c("Fe56 (600 MeV/u)", "Si28", "IEA", "SEA"),
# #        col = c("orange", "cyan", "red", "black"), lwd = c(2, 2, 2, 2),
# #        lty = c(1, 1, 1, 2), cex = 0.6, inset = 0.05)
# #legend(x="bottomright",legend="SD not 95%CI",cex=0.6, inset=0.05)
# #=======================================================#
# #=== Fig. 9 7 HZE Ion mix. Each Ion Contributes 7 cGy ==#
# #=======================================================#
# d9 <- c(0.01 * 0:9, 0.1 * 1:9, 0.5 * 2:100)
# plot(x = d9, y = calculate_SEA(d9, c(25, 70, 100, 193, 250, 464, 953), rep(1/7, 7)), 
#      type = "l", xlab = "Dose (cGy)", ylab = "HG", bty = 'u', col = "black", lwd = 2, lty = 2)
# lines(x = d9, y = calibrated_HZE_nte_der(dose = d9, L = 25), col = "pink", lwd = 2)
# lines(x = d9, y = calibrated_HZE_nte_der(dose = d9, L = 70), col = "orange", lwd = 2)
# lines(x = d9, y = calibrated_HZE_nte_der(dose = d9, L = 100), col = "aquamarine2", lwd = 4)
# lines(x = d9, y = calibrated_HZE_nte_der(dose = d9, L = 193), col = "green", lwd = 2)
# lines(x = d9, y = calibrated_HZE_nte_der(dose = d9, L = 250), col = "blue", lwd = 2)
# lines(x = d9, y = calibrated_HZE_nte_der(dose = d9, L = 464), col = "purple", lwd = 2)
# lines(x = d9, y = calibrated_HZE_nte_der(dose = d9, L = 953), col = "violet", lwd = 2)
# lines(x = d9, y = calculate_id(d9, c(25, 70, 100, 193, 250, 464, 953),
#                                rep(1/7, 7), model = "NTE")[, 2], col = "red", lwd = 3) # I(d)
# legend(x = "topleft", legend = c("Ne20 NTE-TE IDER", "Si28 NTE-TE IDER", 
#                                  "Ti48 NTE-TE IDER", "Fe56 (600 MeV/u) NTE-TE IDER", 
#                                  "Fe56 (300 MeV/u) NTE-TE IDER", "Nb93 NTE-TE IDER",
#                                  "La139 NTE-TE IDER", "IEA MIXDER (Equally Distributed)",
#                                  "SEA MIXDER (Equally Distributed)"),
#        col = c("pink", "orange", "aquamarine2", "green", "blue",
#                "purple", "violet", "red", "black"), 
#        lwd = c(2, 2, 2, 2, 2, 2, 2, 2, 2), 
#        lty = c(1, 1, 1, 1, 1, 1, 1, 1, 2), cex = 0.3, inset = 0.0125)
# 
# 
# #+++++++++++++++++++++ Confidence Interval Ribbon Plots +++++++++++++++++++++++#
# 
# #===========================================================#
# #=== Fig. 10. Fe 600 MeV/u, Si, Equal Doses, CI, mix point====#
# #===========================================================#
# # Fe56 (600 MeV/u) and Si28 in equal proportions for a total of 40 cGy
# # Declare ratios and LET values for plot
# ratios <- c(1/2, 1/2)
# LET_vals <- c(193, 70)
# d10 <- c(0.01 * 0:9, 0.1 * 1:9, 1:40)
# d10_one <- c(0.01 * 0:9, 0.1 * 1:9, 1:20)
# # We use the plot that takes adjustable parameter correlations into account
# corr_fig_10 <- simulate_monte_carlo(n = 500, d10, LET_vals, ratios, model = "NTE")
# # The first argument, n, is the number of Monte Carlo repeats. Increase for
# # greater accuracy. Decrease to speed up the program.
# ci_data <- data.frame(dose = d10,
#                        # Monte Carlo values
#                        corrBottom = corr_fig_10$monte_carlo[1, ],
#                        corrTop = corr_fig_10$monte_carlo[2, ], #
# 
#                        # one-ion DERs for comparison
#                        #fe_six = calibrated_HZE_nte_der(dose = d10, L = 193),
#                        #si = calibrated_HZE_nte_der(dose = d10, L = 70),
# 
#                        # IEA baseline mixture DER I(d), denoted by id below
#                        i = calculate_id(d10, LET_vals, ratios, model = "NTE")[, 2])
# 
# # We make the ribbon plot for correlated parameters
# plot(c(0, 41), c(0, .40), col = "white", bty = 'L', ann = FALSE) # Set plot area
# polygon(x = c(d10, rev(d10)), y = c(ci_data[, "corrTop"], rev(ci_data[, "corrBottom"])),
#          xpd = -1, col = "yellow", lwd = .4, border = "orange") # Narrow CI ribbon
# fe_six_one = calibrated_HZE_nte_der(dose = d10_one, L = 193)
# si_one = calibrated_HZE_nte_der(dose = d10_one, L = 70)
# lines(d10_one, si_one, col = 'brown', lwd = 2) # Si DER
# lines(d10_one, fe_six_one, col = 'blue') # Fe
# lines(ci_data[, "dose"], ci_data[, "i"], col = 'red', lwd = 3) # I(d)
# 
# observed=select(filter(mix_data, Beam == "Si.Fe600"), c(dose,Prev, SD))
# #dd=as.numeric(observed[1]) or as.numeric(observed["dose"]) should also work
# errbar(40,.325, yplus  = .325 +  observed[["SD"]],
#          yminus = .325 -observed[["SD"]],
#          pch = 19, cap = 0.02, add = TRUE, col = 'black', errbar.col = 'black', lwd = 2)
# # fix above so called by name; see Fig. 1 new
# 
# legend(x = "bottomright", legend = "95%CI not SD", cex=0.6)
# # pch = c(19,19), cex = 1, inset = 0.025)
# #=============================================================#
# #== UNDER CONSTRUCTION Fig. 7A. 80% LowLET and 4 HZE Ions ribbon=====#
# #=============================================================#
# ratios <- .01 * c(10, 2.5, 2.5, 5, 80) 
# LET_vals <- c(17, 70, 100, 193, 0.4) 
# d7A <- c(0.01 * 0:9, 0.1 * 1:9, 1:60) # wrong order for low LET?
# corr_fig_7A <- simulate_monte_carlo(n = 500, d7A, LET_vals, ratios, model="NTE")
# ci_data <- data.frame(dose = d7A, 
#                       corrBottom = corr_fig_7A$monte_carlo[1, ],
#                       corrTop = corr_fig_7A$monte_carlo[2, ],
#                       H = calibrated_low_LET_der(dose = d7A, L = .4),
#                       O = calibrated_HZE_nte_der(dose = d7A, L = 17),
#                       Si = calibrated_HZE_nte_der(dose = d7A, L = 70),
#                       Ti = calibrated_HZE_nte_der(dose = d7A, L = 100),
#                       Fe = calibrated_HZE_nte_der(dose = d7A, L = 193),
#                       i = calculate_id(d7A,LET_vals, ratios, model = "NTE")[,2])
#                      
# plot(c(0,60),c(0,1), col = "white", ann = FALSE)
# polygon(x = c(d7A, rev(d7A)), y = c(ci_data[, "corrTop"], rev(ci_data[, "corrBottom"])),
#         xpd = -1, col = "aquamarine2", lwd = .5, border = "aquamarine2") # CI
# lines(ci_data[,"dose"], ci_data[, "Fe"],
#       lwd = 2) # Fe DER. resume here after 7.26.2018
# lines(ci_data[,"dose"], ci_data[, "i"],col = "red",
#       lwd = 2)
# lines(ci_data[,"dose"], ci_data[, "O"],
#       lwd = 2, col = "blue")
# #       col = "green", lwd = 2) # HZE DER
# # lines(x = d7, y = calibrated_HZE_nte_der(dose = d7, L = 110),
# #       col = "purple", lwd = 2)
# # lines(x = d7, y = calibrated_HZE_nte_der(dose = d7, L = 180),
# #       col = "blue", lwd = 2)
# # lines(x = d7, y = calibrated_HZE_nte_der(dose = d7, L = 250),
# #       col = "aquamarine2", lwd = 2)
# # lines(x = d7, y = calculate_id(d7, c(0.4, 40, 110, 180, 250),
# #                                c(.8, rep(.05, 4)), model = "NTE")[, 2], col = "red", lwd = 3) # I(d)
# # legend(x = "topleft", legend = c("Low-LET","L=40", "L=110", "L=180", 
# #                                  "L=250", "I(d)", "S(d)"),
# #        col = c("orange", "green","purple","blue", "aquamarine2", "red", "black"), 
# #        lwd = c(2, 2, 2, 2, 2, 3, 2), 
# #        lty = c(1, 1, 1, 1,  1, 1, 2), cex = 0.3, inset = 0.025)
# 
# #================================================#
# #===== Fig. 7A. 80% LowLET and 4 HZE Ions CI included END =====#
# #================================================#
# 
# 
# 
# 
# #==============================================================================#
# #======================= Fig. 3. Fe 600 MeV/u dual CI plot ====================#
# #==============================================================================#
# # EGH: Need to change double x-axis label to single.
# 
# dose_40 <- 0:40
# dose_120 <- 0:120
# 
# fig_3_corr_1 <- simulate_monte_carlo(n = 500, dose_40, LET = 193, ratios = 1, 
#                                        model = "NTE")
# fig_3_uncorr_1 <- simulate_monte_carlo(n = 500, dose_40, LET = 193, ratios = 1, 
#                                        model = "NTE", vcov = FALSE)
# fig_3_corr_2 <- simulate_monte_carlo(n = 500, dose_120, LET = 193, ratios = 1, 
#                                        model = "NTE")
# fig_3_uncorr_2 <- simulate_monte_carlo(n = 500, dose_120, LET = 193, ratios = 1, 
#                                        model = "NTE", vcov = FALSE)
# 
# par(mfrow = c(1, 2)) # Split plot window.
# plot(dose_40, fig_3_uncorr_1$monte_carlo[2, ], type = "n", # Panel 1.
#      ylab = "Prevalence (%)",
#      xlab = "Dose (cGy)")
# 
# ci_data <- data.frame(dose = dose_40,
#                       # Monte Carlo values
#                       uncorrBottom = fig_3_uncorr_1$monte_carlo[1, ],
#                       uncorrTop = fig_3_uncorr_1$monte_carlo[2, ],
#                       corrBottom = fig_3_corr_1$monte_carlo[1, ],
#                       corrTop = fig_3_corr_1$monte_carlo[2, ])
# polygon(x = c(dose_40, rev(dose_40)), y = c(ci_data[, "uncorrTop"], 
#                                     rev(ci_data[, "uncorrBottom"])),
#         xpd = -1, col = "lightblue", lwd = .5, border = "lightblue")
# polygon(x = c(dose_40, rev(dose_40)), y = c(ci_data[, "corrTop"], 
#                                             rev(ci_data[, "corrBottom"])),
#         xpd = -1, col = "aquamarine2", lwd = .5, border = "aquamarine2")
# lines(dose_40, calibrated_HZE_nte_der(dose_40, 193), col = 'blue', lwd = 1.5)
# 
# 
# plot(dose_120, fig_3_uncorr_2$monte_carlo[2, ], type = "n", # Panel 2.
#      xlab = "Dose (cGy)",
#      ylab = "")
# 
# ci_data <- data.frame(dose = dose_120,
#                       # Monte Carlo values
#                       uncorrBottom = fig_3_uncorr_2$monte_carlo[1, ],
#                       uncorrTop = fig_3_uncorr_2$monte_carlo[2, ],
#                       corrBottom = fig_3_corr_2$monte_carlo[1, ],
#                       corrTop = fig_3_corr_2$monte_carlo[2, ])
# 
# polygon(x = c(dose_120, rev(dose_120)), y = c(ci_data[, "uncorrTop"], 
#                                             rev(ci_data[, "uncorrBottom"])),
#         xpd = -1, col = "lightblue", lwd = .5, border = "lightblue")
# polygon(x = c(dose_120, rev(dose_120)), y = c(ci_data[, "corrTop"], 
#                                               rev(ci_data[, "corrBottom"])),
#         xpd = -1, col = "aquamarine2", lwd = .5, border = "aquamarine2")
# lines(dose_120, calibrated_HZE_nte_der(dose_120, 193), col = 'blue', lwd = 1.5)



