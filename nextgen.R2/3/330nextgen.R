# if (!caracas::has_sympy()) {
#   caracas::install_sympy() 
# }
library(caracas)
library(ggplot2)
#with 1 control/intervention
#construct F (transmission matrix)
Fmat_1 <- matrix_(c("Bc*v", "Bc*v", "Bn*(1-v)", "Bn*(1-v)"), 2, 2)

#construct V (transition matrix)
Vmat_1 <- matrix_(c("-y", "0", "0", "-y"), 2, 2)

#find -V^-1
Vmat_1_inv <- -1*inv(Vmat_1)

#find K = -FV-1
K_1 <- Fmat_1 %*% Vmat_1_inv
#Latex
#tex(K)

#find eigenvalues
K_eigen_1 <- eigenval(K_1)
#store in vector
eigen_vec <- numeric(length(K_eigen_1))
for(i in 1:length(K_eigen_1)){
  eigen_vec[i] <- K_eigen_1[[i]]$eigval
}

K_eigen1 <- K_eigen_1[[1]]$eigval
#inspect eigen_vec manually to find spectral radius (max absolute val)

###combining interventions
#oli's suggestion -- instead of multiplying coverage for both, use a separate parameter

#with 2 controls/interventions
#construct F (transmission matrix)
Fmat_2 <- matrix_(c("Bn*(1-v1)*(1-v2)", "Bn*(1-v1)*(1-v2)", "Bn*(1-v1)*(1-v2)", "Bn*(1-v1)*(1-v2)",
                  "Bc1*v1", "Bc1*v1", "Bc1*v1", "Bc1*v1",
                  "Bc2*v2", "Bc2*v2", "Bc2*v2", "Bc2*v2",
                  "Bc12*v12", "Bc12*v12", "Bc12*v12", "Bc12*v12"), 4, 4)

#construct V (transition matrix)
Vmat_2 <- matrix_(c("-y", "0", "0", "0",
                  "0", "-y", "0", "0",
                  "0", "0", "-y", "0",
                  "0", "0", "0", "-y"), 4, 4)

#find -V^-1
Vmat_2_inv <- -1*inv(Vmat_2)

#find K = -FV-1
K_2 <- Fmat_2 %*% Vmat_2_inv
#Latex
#tex(K)

#find eigenvalues
K_eigen_2 <- eigenval(K_2)
#store in vector
eigen_vec <- numeric(length(K_eigen_2))
for(i in 1:length(K_eigen_2)){
  eigen_vec[i] <- K_eigen_2[[i]]$eigval
}
K_eigen2 <- K_eigen_2[[1]]$eigval

##replace parameters with values
#1 = SRs
#2 = TIRS


epsilon <- 0.04
#mosquito emergence rate
a <- 0.5
#biting rate
b <- 0.38
#Pr(human infected | bitten by infectious mosquito)
c <- 0.38
#Pr(mosquito infected | biting infected human)
g <- 0.07
#daily mortality due to fogging
g_12 <- (g_1^2 + g_2^2)/(g_1 + g_2)
#daily mosquito mortality rate due to SRs + TIRS -- competing hazards
n <- 10
gamma <- 1/7

#beta for intervention 1
Bc1_val <- (epsilon*(a^2)*b*c*exp(-g*n))/(g^2)
#beta for intervention 2
Bc2_val <- (epsilon*(a^2)*b*c*exp(-g*n))/(g^2)
#beta for interventions 1 and 2
Bc12_val <- (epsilon*(a^2)*b*c*exp(-g*n))/(g^2)
#beta for no intervention
Bn_val <- (epsilon*(a^2)*b*c*exp(-g*n))/(g^2)
n <- 17.5
gamma <- 1/7


#coverage for intervention 1
v1_val <- 0.5
#coverage for intervention 2
v2_val <- 0.5
v12_val <- v1_val*v2_val
#less arbitrary way to indicate non-independence?

#R0 when combining interventions 1 and 2
R0_12 <- subs(K_eigen2, list("Bn" = Bn_val,
                                 "Bc1" = Bc1_val,
                                 "Bc2" = Bc2_val,
                                 "Bc12" = Bc12_val,
                                 "v1" = v1_val,
                                 "v2" = v2_val,
                                 "v12" = v12_val,
                                 "y" = gamma))

##show with 1 each and compare
#same level of coverage (50%)

#intervention 1 only
R0_1 <- subs(K_eigen1, list("Bn" = Bn_val,
                                "Bc" = Bc1_val,
                                "v" = v1_val,
                                "y" = gamma))

#intervention 2 only
R0_2 <- subs(K_eigen1, list("Bn" = Bn_val,
                                "Bc" = Bc2_val,
                                "v" = v2_val,
                                "y" = gamma))

##plot combos of coverages -- SRs = v1, TIRS = v2
v_vals <- seq(0, 1, 0.05)

df <- expand.grid(v_vals, v_vals)

names(df)[names(df) == "Var1"] <- "v1"
names(df)[names(df) == "Var2"] <- "v2"

df$R0 <- numeric(length(df$v1))

for(i in 1:length(df$v1)){
  df$R0[i] <- as.numeric(as_expr(subs(K_eigen2, list("Bn" = Bn_val,
                                                         "Bc1" = Bc1_val,
                                                         "Bc2" = Bc2_val,
                                                         "Bc12" = Bc12_val,
                                                         "v1" = df$v1[i],
                                                         "v2" = df$v2[i],
                                                         "v12" = df$v1[i]*df$v2[i],
                                                         "y" = gamma))))
}

#ggplot2 heatmap
ggplot(df, aes(v1, v2, fill= R0)) + 
  geom_tile() +  
  xlab("Larviciding") + ylab("Fogging") +
  scale_fill_gradient(low="khaki1", high="red") +
  theme_bw()

df[which.max(df$R0),]
df[which.min(df$R0),]
