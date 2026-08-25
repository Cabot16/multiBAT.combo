# ------------------------ Setup ------------------------

library(ggplot2)
library(plyr)
library(mvtnorm)

# Common parameters
beta11 <- 1
beta12 <- -1
beta21 <- 1
beta22 <- -1
LD2 <- 0.6
vare <- 100
alpha <- 0.05
m <- 2
nrep <- 5000

n1 <- 500
n2 <- 500

r1_vec <- c(0, 0.2, 0.4, 0.6, 0.8, 0.99)

font_family <- "Arial"

BASE_SIZE <- 18
TITLE_SIZE <- 22
TAG_SIZE <- 28
TICK_SIZE <- BASE_SIZE + 2
AXIS_TITLE_SIZE <- BASE_SIZE + 2

DATE_TAG <- format(Sys.Date(), "%m%d")
OUTPUT_DIR <- "."
OUTPUT_FILE_STUB <- "FigureS1"

CURVE_SPAN <- 0.9

method_levels <- c(
  "Population 1",
  "Population 2",
  "Single-ancestry meta",
  "multiBAT",
  "multiBAT-combo"
)

method_colors <- c(
  "Population 1"         = "#66ccff",
  "Population 2"         = "#DBB71D",
  "Single-ancestry meta" = "#A1C28E",
  "multiBAT"             = "#873186",
  "multiBAT-combo"       = "#39C5BB"
)

# ------------------------ Helper functions ------------------------

ACAT_weighted <- function(p, w) {
  valid <- !is.na(p)
  p <- p[valid]
  w <- w[valid]
  
  if (length(p) == 0) {
    return(NA_real_)
  }
  
  w <- w / sum(w)
  p[p == 1] <- 1 - 1e-16
  is.small <- p < 1e-16
  
  if (sum(is.small) == 0) {
    cct.stat <- sum(w * tan((0.5 - p) * pi))
  } else {
    cct1 <- sum(w[is.small] * ((1 / p[is.small]) / pi))
    cct2 <- sum(w[!is.small] * tan((0.5 - p[!is.small]) * pi))
    cct.stat <- cct1 + cct2
  }
  
  if (cct.stat > 1e15) {
    pval <- (1 / cct.stat) / pi
  } else {
    pval <- 0.5 - atan(cct.stat) / pi
  }
  
  return(pval)
}

save_plot_pair <- function(plot_obj, file_stub, width, height) {
  if (!dir.exists(OUTPUT_DIR)) {
    dir.create(OUTPUT_DIR, recursive = TRUE)
  }
  
  ggsave(
    filename = file.path(OUTPUT_DIR, paste0(file_stub, "_", DATE_TAG, ".png")),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 600,
    bg = "white"
  )
  
  ggsave(
    filename = file.path(OUTPUT_DIR, paste0(file_stub, "_", DATE_TAG, ".pdf")),
    plot = plot_obj,
    width = width,
    height = height,
    device = cairo_pdf,
    bg = "white"
  )
}

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
    data.frame(r1 = LD1, Method = "Population 1", Power = power1, Type = "Theory"),
    data.frame(r1 = LD1, Method = "Population 2", Power = power2, Type = "Theory"),
    data.frame(r1 = LD1, Method = "multiBAT", Power = power_comb, Type = "Theory")
  )
}

# ------------------------ Empirical power ------------------------

res.emp <- data.frame()

for (LD1 in r1_vec) {
  R1 <- matrix(c(1, LD1, LD1, 1), nrow = 2)
  eig1 <- eigen(R1, symmetric = TRUE)
  z1_list <- vector("list", nrep)
  
  for (k in 1:nrep) {
    X1 <- rmvnorm(n1, mean = c(0, 0), sigma = R1)
    y1 <- X1 %*% c(beta11, beta12) + rnorm(n1, 0, sqrt(vare))
    fit1 <- apply(X1, 2, function(x) {
      summary(lm(y1 ~ x))$coefficients[2, 1:2]
    })
    z1_list[[k]] <- fit1[1, ] / fit1[2, ]
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
    
    w_pop1_meta <- n1 / (n1 + n2)
    w_pop2_meta <- n2 / (n1 + n2)
    
    P_single_meta <- ACAT_weighted(
      p = c(P_mbat1, P_mbat2),
      w = c(w_pop1_meta, w_pop2_meta)
    )
    
    w_joint <- 1 - abs(n1 - n2) / (n1 + n2)
    w_joint <- max(0, min(1, w_joint))
    
    w_pop1 <- (1 - w_joint) * n1 / (n1 + n2)
    w_pop2 <- (1 - w_joint) * n2 / (n1 + n2)
    
    P_multibat_combo <- ACAT_weighted(
      p = c(P_mbat1, P_mbat2, P_mbat_combined),
      w = c(w_pop1, w_pop2, w_joint)
    )
    
    res.emp <- rbind(
      res.emp,
      data.frame(r1 = LD1, Method = "Population 1", P = P_mbat1),
      data.frame(r1 = LD1, Method = "Population 2", P = P_mbat2),
      data.frame(r1 = LD1, Method = "Single-ancestry meta", P = P_single_meta),
      data.frame(r1 = LD1, Method = "multiBAT", P = P_mbat_combined),
      data.frame(r1 = LD1, Method = "multiBAT-combo", P = P_multibat_combo)
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

combined$Method <- factor(combined$Method, levels = method_levels)
combined$Type <- factor(combined$Type, levels = c("Theory", "Empirical"))

plot_power <- function(data,
                       linewidth_empirical = 1.2,
                       linewidth_theory = 0.8,
                       point_size = 2.8,
                       curve_span = CURVE_SPAN) {
  
  ggplot(
    data,
    aes(
      x = r1,
      y = Power,
      colour = Method,
      group = interaction(Method, Type)
    )
  ) +
    geom_point(
      aes(shape = Type),
      size = point_size
    ) +
    geom_smooth(
      aes(linetype = Type, linewidth = Type),
      method = "loess",
      formula = y ~ x,
      se = FALSE,
      span = curve_span
    ) +
    
    scale_color_manual(
      values = method_colors,
      breaks = method_levels,
      labels = method_levels,
      drop = FALSE
    ) +
    
    scale_linetype_manual(
      values = c("Theory" = "solid", "Empirical" = "22"),
      breaks = c("Theory", "Empirical"),
      labels = c("Theory", "Empirical")
    ) +
    
    scale_shape_manual(
      values = c("Theory" = 16, "Empirical" = 17),
      breaks = c("Theory", "Empirical"),
      labels = c("Theory", "Empirical")
    ) +
    
    scale_linewidth_manual(
      values = c("Theory" = linewidth_theory, "Empirical" = linewidth_empirical),
      guide = "none"
    ) +
    
    scale_y_continuous(
      breaks = c(0, 0.25, 0.50, 0.75, 1.00),
      labels = sprintf("%.2f", c(0, 0.25, 0.50, 0.75, 1.00)),
      expand = c(0, 0.02)
    ) +
    
    coord_cartesian(
      ylim = c(0, 1)
    ) +
    
    scale_x_continuous(
      breaks = r1_vec,
      labels = sprintf("%.2f", r1_vec)
    ) +
    
    labs(
      title = " ",
      x = expression(r[1]),
      y = "Power",
      colour = NULL,
      linetype = NULL,
      shape = NULL
    ) +
    
    guides(
      colour = guide_legend(order = 1),
      linetype = guide_legend(order = 2),
      shape = guide_legend(order = 2)
    ) +
    
    theme_classic(base_size = BASE_SIZE, base_family = font_family) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_blank(),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 1.0
      ),
      axis.ticks = element_line(colour = "black", linewidth = 1.0),
      axis.text = element_text(
        colour = "black",
        size = TICK_SIZE,
        face = "plain"
      ),
      axis.title = element_text(
        colour = "black",
        size = AXIS_TITLE_SIZE,
        face = "bold"
      ),
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = TITLE_SIZE
      ),
      legend.title = element_blank(),
      legend.position = "right",
      legend.text = element_text(
        size = BASE_SIZE,
        face = "plain"
      ),
      plot.margin = margin(t = 10, r = 15, b = 10, l = 10)
    )
}

p <- plot_power(
  combined,
  linewidth_empirical = 1.2,
  linewidth_theory = 0.8,
  point_size = 2.8,
  curve_span = CURVE_SPAN
)

print(p)

save_plot_pair(
  plot_obj = p,
  file_stub = OUTPUT_FILE_STUB,
  width = 16,
  height = 7
)
