# ------------------------ Setup ------------------------
library(ggplot2)
library(plyr)
library(mvtnorm)

# Common parameters
beta11 <- 0.5
beta12 <- -0.5
beta21 <- 0.5
beta22 <- -0.5
LD2 <- 0.5
vare <- 100
alpha <- 0.05
m <- 2
nrep <- 10000

n1 <- 1000
n2 <- 1000

r1_vec <- c(0, 0.2, 0.4, 0.6, 0.8, 0.99)
font_family <- "Arial"

# ------------------------ Theoretical power ------------------------
res.theory <- data.frame()

for (LD1 in r1_vec) {
  R1 <- matrix(c(1, LD1, LD1, 1), nrow = 2)
  R2 <- matrix(c(1, LD2, LD2, 1), nrow = 2)
  
  eig1 <- eigen(R1, symmetric = TRUE)
  U1 <- eig1$vectors
  D1 <- eig1$values
  
  varg1 <- beta11^2 + beta12^2 + 2 * LD1 * beta11 * beta12
  vary1 <- varg1 + vare
  b1_1 <- beta11 + LD1 * beta12
  b2_1 <- beta12 + LD1 * beta11
  q2_1 <- b1_1^2 / vary1
  q2_2 <- b2_1^2 / vary1
  mu1 <- sqrt(n1 / vary1 / (1 - q2_1)) * b1_1
  mu2 <- sqrt(n1 / vary1 / (1 - q2_2)) * b2_1
  mu_vec1 <- c(mu1, mu2)
  RinvSqrt1 <- U1 %*% diag(sqrt(1 / D1)) %*% t(U1)
  mu_transf1 <- t(U1) %*% RinvSqrt1 %*% mu_vec1
  ncp1 <- sum(mu_transf1^2)
  
  eig2 <- eigen(R2, symmetric = TRUE)
  U2 <- eig2$vectors
  D2 <- eig2$values
  varg2 <- beta21^2 + beta22^2 + 2 * LD2 * beta21 * beta22
  vary2 <- varg2 + vare
  b1_2 <- beta21 + LD2 * beta22
  b2_2 <- beta22 + LD2 * beta21
  q2_3 <- b1_2^2 / vary2
  q2_4 <- b2_2^2 / vary2
  mu3 <- sqrt(n2 / vary2 / (1 - q2_3)) * b1_2
  mu4 <- sqrt(n2 / vary2 / (1 - q2_4)) * b2_2
  mu_vec2 <- c(mu3, mu4)
  RinvSqrt2 <- U2 %*% diag(sqrt(1 / D2)) %*% t(U2)
  mu_transf2 <- t(U2) %*% RinvSqrt2 %*% mu_vec2
  ncp2 <- sum(mu_transf2^2)
  
  mu_comb <- c(mu_vec1, mu_vec2)
  Sigma_comb <- matrix(0, 4, 4)
  Sigma_comb[1:2, 1:2] <- R1
  Sigma_comb[3:4, 3:4] <- R2
  eig_comb <- eigen(Sigma_comb, symmetric = TRUE)
  U_comb <- eig_comb$vectors
  D_comb <- eig_comb$values
  RinvSqrt_comb <- U_comb %*% diag(sqrt(1 / D_comb)) %*% t(U_comb)
  mu_transf_comb <- t(U_comb) %*% RinvSqrt_comb %*% mu_comb
  ncp_comb <- sum(mu_transf_comb^2)
  
  crit2 <- qchisq(1 - alpha, df = 2)
  crit4 <- qchisq(1 - alpha, df = 4)
  power1 <- 1 - pchisq(crit2, df = 2, ncp = ncp1)
  power2 <- 1 - pchisq(crit2, df = 2, ncp = ncp2)
  power_comb <- 1 - pchisq(crit4, df = 4, ncp = ncp_comb)
  
  res.theory <- rbind(
    res.theory,
    data.frame(r1 = LD1, Method = "mbat1", Power = power1, Type = "Theory"),
    data.frame(r1 = LD1, Method = "mbat2", Power = power2, Type = "Theory"),
    data.frame(r1 = LD1, Method = "mbat_combined", Power = power_comb, Type = "Theory")
  )
}

# ------------------------ Empirical power ------------------------
res.emp <- data.frame()

for (LD1 in r1_vec) {
  R1 <- matrix(c(1, LD1, LD1, 1), nrow = 2)
  eig1 <- eigen(R1, symmetric = TRUE)
  z1_list <- vector("list", nrep)
  se1_list <- vector("list", nrep)
  
  for (k in 1:nrep) {
    X1 <- rmvnorm(n1, mean = c(0, 0), sigma = R1)
    y1 <- X1 %*% c(beta11, beta12) + rnorm(n1, 0, sqrt(vare))
    fit1 <- apply(X1, 2, function(x) {
      summary(lm(y1 ~ x))$coefficients[2, 1:2]
    })
    z1_list[[k]] <- fit1[1, ] / fit1[2, ]
    se1_list[[k]] <- fit1[2, ]
  }
  
  R2 <- matrix(c(1, LD2, LD2, 1), nrow = 2)
  eig2 <- eigen(R2, symmetric = TRUE)
  
  for (k in 1:nrep) {
    z1 <- z1_list[[k]]
    
    X2 <- rmvnorm(n2, mean = c(0, 0), sigma = R2)
    y2 <- X2 %*% c(beta21, beta22) + rnorm(n2, 0, sqrt(vare))
    fit2 <- apply(X2, 2, function(x) {
      summary(lm(y2 ~ x))$coefficients[2, 1:2]
    })
    z2 <- fit2[1, ] / fit2[2, ]
    
    z_combined <- c(z1, z2)
    Sigma <- matrix(0, 4, 4)
    Sigma[1:2, 1:2] <- R1
    Sigma[3:4, 3:4] <- R2
    eig_comb <- eigen(Sigma, symmetric = TRUE)
    
    UPz_comb <- crossprod(eig_comb$vectors, z_combined)
    chisq_comb <- crossprod(UPz_comb, 1 / eig_comb$values * UPz_comb)
    P_mbat_combined <- pchisq(chisq_comb, df = 4, lower.tail = FALSE)
    
    UPz1 <- crossprod(eig1$vectors, z1)
    chisq1 <- crossprod(UPz1, 1 / eig1$values * UPz1)
    P_mbat1 <- pchisq(chisq1, df = m, lower.tail = FALSE)
    
    UPz2 <- crossprod(eig2$vectors, z2)
    chisq2 <- crossprod(UPz2, 1 / eig2$values * UPz2)
    P_mbat2 <- pchisq(chisq2, df = m, lower.tail = FALSE)
    
    res.emp <- rbind(
      res.emp,
      data.frame(r1 = LD1, Method = "mbat1", P = P_mbat1),
      data.frame(r1 = LD1, Method = "mbat2", P = P_mbat2),
      data.frame(r1 = LD1, Method = "mbat_combined", P = P_mbat_combined)
    )
  }
}

dd.emp <- ddply(res.emp, .(r1, Method), summarise, Power = mean(P < alpha))
dd.emp$Type <- "Empirical"

# ------------------------ Merge and plot ------------------------
combined <- rbind(
  res.theory[, c("r1", "Method", "Power", "Type")],
  dd.emp[, c("r1", "Method", "Power", "Type")]
)

plot_power <- function(data, 
                       linewidth_empirical = 1.2,
                       linewidth_theory = 0.8, 
                       point_size = 2.5) {
  
  ggplot(data, aes(x = r1, y = Power, color = Method, 
                   linetype = Type, shape = Type, 
                   group = interaction(Method, Type))) + 
    
    geom_point(size = point_size) + 
    geom_line(aes(linewidth = Type)) + 
    
    scale_linetype_manual(values = c("Empirical" = "22", "Theory" = "solid")) + 
    
    scale_color_manual(values = c("mbat1" = "#66ccff", 
                                  "mbat2" = "#DBB71D", 
                                  "mbat_combined" = "#873186"), 
                       breaks = c("mbat1", "mbat2", "mbat_combined"), 
                       labels = c("mBAT Population 1", 
                                  "mBAT Population 2", 
                                  "multiBAT")) +  
    
    scale_linewidth_manual(values = c("Empirical" = linewidth_empirical, 
                                      "Theory" = linewidth_theory)) + 
    
    scale_shape_manual(values = c("Empirical" = 17, "Theory" = 16)) + 
    
    scale_y_continuous(limits = c(0, 1), 
                       breaks = c(0, 0.25, 0.50, 0.75, 1.00),
                       expand = c(0, 0.02)) +
    
    scale_x_continuous(breaks = r1_vec) +
    
    labs(title = " ", x = expression(r[1]), y = "Power") + 
    
    theme_classic(base_size = 18, base_family = font_family) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", size = 0.8),
      axis.ticks = element_line(color = "black", size = 0.8),
      axis.text = element_text(color = "black", size = 18),
      axis.title = element_text(color = "black", size = 18),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 22),
      legend.title = element_blank(),
      legend.position = "right",
      legend.text = element_text(size = 18),
      plot.margin = margin(t = 10, r = 15, b = 10, l = 10)
    )
}

p <- plot_power(combined, linewidth_empirical = 1.2, linewidth_theory = 0.8)

print(p)

ggsave("sim_power_r1_style_matched.png", p, width = 16, height = 7, dpi = 600)
ggsave("sim_power_r1_style_matched.pdf", p, width = 16, height = 7, device = cairo_pdf)
