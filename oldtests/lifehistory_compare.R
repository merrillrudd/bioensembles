### from life history studies, Adams and Rios 2015
study_vec <- c(NA, "Puerto Rico", "Cuba", "GOM", "GOM", "GOM", "GOM", NA, "Dry Tortugas", "FL Keys", "S. Florida", "GOM", "Marquesas Key FL", "S. Florida", "GOM", "GOM", "Florida", "W. Florida", "SE US")
loo_vec <- c(566, 913, 850, 966, NA, 381, 896, 566, 651, 336, 428, 917, 397, 448, 939, NA, NA, 849, NA)
vbk_vec <- c(0.19, 0.08, 0.10, 0.08, NA, 0.56, 0.09, 0.19, 0.11, 0.56, 0.26, 0.08, 0.17, 0.23, 0.07, NA, NA, 0.11, NA)
t0_vec <- c(-0.78, -1.78, -1.38, -1.77, NA, -0.16, -1.98, -2.33, -3.34, -0.19, -0.93, -1.84, -3.74, -1.02, -2.01, NA, NA, -1.33, NA)
M_vec <- c(0.25, 0.13, rep(NA, 16), mean(c(0.16,0.29)))
Lm_vec <- c(NA, 249, NA, 152, 152, NA, NA, 198, NA, NA, NA, NA, NA, 163, 166, 177, NA, NA, mean(c(152, 193)))

lhstudies <- data.frame("Location"=study_vec, "Loo"=loo_vec, "K"=vbk_vec, "t0"=t0_vec, "M"=M_vec, "Lm"=Lm_vec)

a_vec <- c(2.55e-5, 1.52e-2, 2.37e-2, 9.5e-5)
b_vec <- c(2.97, 3.11, 2.95, 2.75)