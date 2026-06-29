# =======================================================
# 1. Setup & Parameters
# =======================================================
library(ggplot2)
library(plyr)
library(mvtnorm)
library(dplyr)
library(tidyr)
library(patchwork)
library(rstatix)

# Common parameters
beta11 <- 0.5; beta12 <- -0.5
beta21 <- 0.5; beta22 <- -0.5
LD1 <- 0.5; LD2 <- 0.5
vare <- 100
alpha <- 0.05
m <- 2
nrep <- 10000

nvec1 <- seq(100, 1000, by = 100)
n2 <- 1000

# Global Color Mapping
method_colors <- c(
  "pop1"           = "#66ccff",
  "pop2"           = "#DBB71D",
  "multiBAT"       = "#873186",
  "multiBAT-combo" = "#39C5BB"
)
# Global Font
font_family <- "Arial"
# =======================================================
# 2. Functions (ACAT & Simulation)
# =======================================================
ACAT_weighted <- function(p, w) {
  valid <- !is.na(p)
  p <- p[valid]
  w <- w[valid]
  if (length(p) == 0) return(NA_real_)
  
  w <- w / sum(w)
  p[p == 1] <- 1 - 1e-16
  is.small <- (p < 1e-16)
  
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

simulate_power <- function(n1, n2 = 1000, nrep_inner = 100, nrep_outer = 100,
                           beta1_val = 0.8, beta2_val = -0.8) {
  res.box <- data.frame()
  
  R1 <- matrix(c(1, LD1, LD1, 1), nrow = 2)
  R2 <- matrix(c(1, LD2, LD2, 1), nrow = 2)
  eig1 <- eigen(R1, symmetric = TRUE)
  eig2 <- eigen(R2, symmetric = TRUE)
  
  get_z_se <- function(X, y) {
    res <- apply(X, 2, function(x) {
      fit <- summary(lm(y ~ x))$coefficients[2, 1:2]
      c(z = fit[1] / fit[2], se = fit[2])
    })
    list(z = res[1, ], se = res[2, ])
  }
  
  for (outer in 1:nrep_outer) {
    P_mbat1 <- numeric(nrep_inner)
    P_mbat2 <- numeric(nrep_inner)
    P_mbat_combined <- numeric(nrep_inner)
    P_ACATO_weighted <- numeric(nrep_inner)
    
    for (k in 1:nrep_inner) {
      X1 <- rmvnorm(n1, mean = c(0, 0), sigma = R1)
      y1 <- X1 %*% c(beta1_val, beta2_val) + rnorm(n1, 0, sqrt(vare))
      zse1 <- get_z_se(X1, y1)
      z1 <- zse1$z
      
      X2 <- rmvnorm(n2, mean = c(0, 0), sigma = R2)
      y2 <- X2 %*% c(beta1_val, beta2_val) + rnorm(n2, 0, sqrt(vare))
      zse2 <- get_z_se(X2, y2)
      z2 <- zse2$z
      
      z_combined <- c(z1, z2)
      Sigma <- matrix(0, 4, 4)
      Sigma[1:2, 1:2] <- R1
      Sigma[3:4, 3:4] <- R2
      eig_comb <- eigen(Sigma, symmetric = TRUE)
      
      UPz_comb <- crossprod(eig_comb$vectors, z_combined)
      chisq_comb <- crossprod(UPz_comb, (1 / eig_comb$values) * UPz_comb)
      P_mbat_combined[k] <- pchisq(chisq_comb, df = 4, lower.tail = FALSE)
      
      UPz1 <- crossprod(eig1$vectors, z1)
      chisq1 <- crossprod(UPz1, (1 / eig1$values) * UPz1)
      P_mbat1[k] <- pchisq(chisq1, df = m, lower.tail = FALSE)
      
      UPz2 <- crossprod(eig2$vectors, z2)
      chisq2 <- crossprod(UPz2, (1 / eig2$values) * UPz2)
      P_mbat2[k] <- pchisq(chisq2, df = m, lower.tail = FALSE)
      
      w_joint <- 1 - abs(n1 - n2) / (n1 + n2)
      w_joint <- max(0, min(1, w_joint))
      w_pop1 <- (1 - w_joint) * n1 / (n1 + n2)
      w_pop2 <- (1 - w_joint) * n2 / (n1 + n2)
      
      P_ACATO_weighted[k] <- ACAT_weighted(
        p = c(P_mbat1[k], P_mbat2[k], P_mbat_combined[k]),
        w = c(w_pop1, w_pop2, w_joint)
      )
    }
    
    res.box <- rbind(
      res.box,
      data.frame(N1 = n1, Replicate = outer, Method = "pop1",           Power = mean(P_mbat1 < alpha)),
      data.frame(N1 = n1, Replicate = outer, Method = "pop2",           Power = mean(P_mbat2 < alpha)),
      data.frame(N1 = n1, Replicate = outer, Method = "multiBAT",       Power = mean(P_mbat_combined < alpha)),
      data.frame(N1 = n1, Replicate = outer, Method = "multiBAT-combo", Power = mean(P_ACATO_weighted < alpha))
    )
  }
  
  return(res.box)
}

# =======================================================
# 3. Data Generation
# =======================================================

# ------------------------
# 3A. Theoretical power
# ------------------------
res.theory <- data.frame()

for (n1 in nvec1) {
  R1 <- matrix(c(1, LD1, LD1, 1), nrow = 2)
  R2 <- matrix(c(1, LD2, LD2, 1), nrow = 2)
  
  eig1 <- eigen(R1, symmetric = TRUE)
  varg1 <- beta11^2 + beta12^2 + 2 * LD1 * beta11 * beta12
  vary1 <- varg1 + vare
  b1_1  <- beta11 + LD1 * beta12
  b2_1  <- beta12 + LD1 * beta11
  mu_vec1 <- c(
    sqrt(n1 / vary1 / (1 - (b1_1^2 / vary1))) * b1_1,
    sqrt(n1 / vary1 / (1 - (b2_1^2 / vary1))) * b2_1
  )
  ncp1 <- as.numeric(t(mu_vec1) %*% solve(R1) %*% mu_vec1)
  
  eig2 <- eigen(R2, symmetric = TRUE)
  varg2 <- beta21^2 + beta22^2 + 2 * LD2 * beta21 * beta22
  vary2 <- varg2 + vare
  b1_2  <- beta21 + LD2 * beta22
  b2_2  <- beta22 + LD2 * beta21
  mu_vec2 <- c(
    sqrt(n2 / vary2 / (1 - (b1_2^2 / vary2))) * b1_2,
    sqrt(n2 / vary2 / (1 - (b2_2^2 / vary2))) * b2_2
  )
  ncp2 <- as.numeric(t(mu_vec2) %*% solve(R2) %*% mu_vec2)
  
  ncp_comb <- ncp1 + ncp2
  
  power1     <- pchisq(qchisq(1 - alpha, df = 2), df = 2, ncp = ncp1, lower.tail = FALSE)
  power2     <- pchisq(qchisq(1 - alpha, df = 2), df = 2, ncp = ncp2, lower.tail = FALSE)
  power_comb <- pchisq(qchisq(1 - alpha, df = 4), df = 4, ncp = ncp_comb, lower.tail = FALSE)
  
  res.theory <- rbind(
    res.theory,
    data.frame(N1 = n1, Method = "pop1",     Power = power1,     Type = "Theory"),
    data.frame(N1 = n1, Method = "pop2",     Power = power2,     Type = "Theory"),
    data.frame(N1 = n1, Method = "multiBAT", Power = power_comb, Type = "Theory")
  )
}

# ------------------------
# 3B. Boxplot data
# ------------------------
res_null <- simulate_power(420, 1000, nrep_inner = 100, nrep_outer = 50, beta1_val = 0, beta2_val = 0)
res_200  <- simulate_power(200, 1000, nrep_inner = 100, nrep_outer = 50)
res_1000 <- simulate_power(1000, 1000, nrep_inner = 100, nrep_outer = 50)

# ------------------------
# 3C. Simulation line data for new D panel
# ------------------------
res.emp.line <- data.frame()

get_z_only <- function(X, y) {
  apply(X, 2, function(x) {
    fit <- summary(lm(y ~ x))$coefficients[2, 1:2]
    fit[1] / fit[2]
  })
}

for (n1 in nvec1) {
  R1 <- matrix(c(1, LD1, LD1, 1), nrow = 2)
  R2 <- matrix(c(1, LD2, LD2, 1), nrow = 2)
  
  eig1 <- eigen(R1, symmetric = TRUE)
  eig2 <- eigen(R2, symmetric = TRUE)
  
  P1 <- numeric(nrep)
  P2 <- numeric(nrep)
  Pc <- numeric(nrep)
  Pcombo <- numeric(nrep)
  
  for (k in 1:nrep) {
    X1 <- rmvnorm(n1, mean = c(0, 0), sigma = R1)
    y1 <- X1 %*% c(beta11, beta12) + rnorm(n1, 0, sqrt(vare))
    z1 <- get_z_only(X1, y1)
    
    X2 <- rmvnorm(n2, mean = c(0, 0), sigma = R2)
    y2 <- X2 %*% c(beta21, beta22) + rnorm(n2, 0, sqrt(vare))
    z2 <- get_z_only(X2, y2)
    
    UPz1 <- crossprod(eig1$vectors, z1)
    P1[k] <- pchisq(sum(UPz1^2 / eig1$values), df = 2, lower.tail = FALSE)
    
    UPz2 <- crossprod(eig2$vectors, z2)
    P2[k] <- pchisq(sum(UPz2^2 / eig2$values), df = 2, lower.tail = FALSE)
    
    z_combined <- c(z1, z2)
    Sigma <- matrix(0, 4, 4)
    Sigma[1:2, 1:2] <- R1
    Sigma[3:4, 3:4] <- R2
    eig_comb <- eigen(Sigma, symmetric = TRUE)
    
    UPzc <- crossprod(eig_comb$vectors, z_combined)
    Pc[k] <- pchisq(sum(UPzc^2 / eig_comb$values), df = 4, lower.tail = FALSE)
    
    w_joint <- 1 - abs(n1 - n2) / (n1 + n2)
    w_joint <- max(0, min(1, w_joint))
    w_pop1 <- (1 - w_joint) * n1 / (n1 + n2)
    w_pop2 <- (1 - w_joint) * n2 / (n1 + n2)
    
    Pcombo[k] <- ACAT_weighted(
      p = c(P1[k], P2[k], Pc[k]),
      w = c(w_pop1, w_pop2, w_joint)
    )
  }
  
  res.emp.line <- rbind(
    res.emp.line,
    data.frame(N1 = n1, Method = "pop1",           Power = mean(P1 < alpha),     Type = "Simulation"),
    data.frame(N1 = n1, Method = "pop2",           Power = mean(P2 < alpha),     Type = "Simulation"),
    data.frame(N1 = n1, Method = "multiBAT",       Power = mean(Pc < alpha),     Type = "Simulation"),
    data.frame(N1 = n1, Method = "multiBAT-combo", Power = mean(Pcombo < alpha), Type = "Simulation")
  )
}

# =======================================================
# 4. Significance Calculation
# =======================================================
prepare_sig_data <- function(dat) {
  dat_filtered <- dat %>%
    filter(Method %in% c("pop2", "multiBAT", "multiBAT-combo")) %>%
    mutate(Method = factor(Method, levels = c("pop2", "multiBAT", "multiBAT-combo")))
  
  stat.test <- dat_filtered %>%
    pairwise_wilcox_test(Power ~ Method, paired = TRUE) %>%
    add_significance() %>%
    mutate(stars = p.adj.signif)
  
  tag_data <- data.frame(
    start   = c(2, 3, 2),
    end     = c(3, 4, 4),
    y       = c(0.95, 1.00, 1.05),
    p.stars = as.character(stat.test$stars[c(1, 3, 2)]),
    h       = c(0.015, 0.015, 0.015)
  )
  
  return(list(tag_data = tag_data, max_y = 1.10))
}

sig_200  <- prepare_sig_data(res_200)
sig_1000 <- prepare_sig_data(res_1000)

global_ymax <- 1.10

# =======================================================
# 5. Final Plotting
# =======================================================
make_boxplot <- function(dat, title_text = NULL, y_label = "Power",
                         drop_y = FALSE, add_hline = FALSE, sig_data = NULL) {
  dat$Method <- factor(dat$Method, levels = names(method_colors))
  
  p <- ggplot(dat, aes(x = Method, y = Power, fill = Method)) +
    geom_boxplot(size = 1.1, alpha = 0.9) +
    scale_fill_manual(
      values = method_colors,
      breaks = c("pop1", "pop2", "multiBAT", "multiBAT-combo"),
      labels = c("mBAT Population 1", "mBAT Population 2", "multiBAT", "multiBAT-combo")
    ) +
    theme_bw(base_size = 18, base_family = font_family) +
    labs(title = title_text, y = y_label, x = NULL) +
    scale_y_continuous(
      breaks = c(0, 0.25, 0.50, 0.75, 1.00)
    ) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 22)
    ) +
    coord_cartesian(ylim = c(0, global_ymax))
  
  if (add_hline) {
    p <- p + geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", size = 1.1)
  }
  
  if (drop_y) {
    p <- p + theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  }
  
  if (!is.null(sig_data)) {
    p <- p +
      geom_segment(
        data = sig_data$tag_data,
        aes(x = start, xend = end, y = y, yend = y),
        inherit.aes = FALSE,
        linewidth = 0.8
      ) +
      geom_segment(
        data = sig_data$tag_data,
        aes(x = start, xend = start, y = y - h, yend = y),
        inherit.aes = FALSE,
        linewidth = 0.8
      ) +
      geom_segment(
        data = sig_data$tag_data,
        aes(x = end, xend = end, y = y - h, yend = y),
        inherit.aes = FALSE,
        linewidth = 0.8
      ) +
      geom_text(
        data = sig_data$tag_data,
        aes(
          x = (start + end) / 2,
          y = y + 0.02,
          label = p.stars,
          size = ifelse(p.stars == "ns", 5, 8)   # ⭐关键
        ),
        inherit.aes = FALSE
      ) +
      scale_size_identity()
  }
  
  return(p)
}

# ------------------------
# 5A. Panel A
# ------------------------
p_null <- make_boxplot(res_null, title_text = " ", y_label = "FPR", add_hline = TRUE)
p_1000 <- make_boxplot(res_1000, title_text = "Balanced Samples", sig_data = sig_1000)
p_200  <- make_boxplot(res_200,  title_text = "Imbalanced Samples", drop_y = TRUE, sig_data = sig_200)

panel_A <- (p_null | p_1000 | p_200) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

# ------------------------
# 5B. New Panel D / panel_B
# Theory + Simulation + multiBAT-combo
# ------------------------
combined_line <- bind_rows(res.theory, res.emp.line) %>%
  mutate(
    Method = factor(Method, levels = c("pop1", "pop2", "multiBAT", "multiBAT-combo")),
    Type   = factor(Type, levels = c("Theory", "Simulation"))
  )

panel_B <- ggplot(
  combined_line,
  aes(
    x = N1,
    y = Power,
    color = Method,
    linetype = Type,
    shape = Type,
    group = interaction(Method, Type)
  )
) +
  geom_point(size = 3) +
  geom_line(aes(linewidth = Type)) +
  scale_color_manual(
    values = method_colors,
    breaks = c("pop1", "pop2", "multiBAT", "multiBAT-combo"),
    labels = c("mBAT Population 1", "mBAT Population 2", "multiBAT", "multiBAT-combo")
  ) +
  scale_linetype_manual(
    values = c("Theory" = "solid", "Simulation" = "22")
  ) +
  scale_shape_manual(
    values = c("Theory" = 16, "Simulation" = 17)
  ) +
  scale_linewidth_manual(
    values = c("Theory" = 0.8, "Simulation" = 1.2)
  ) +
  scale_y_continuous(
    limits = c(0, 1.05),
    breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    expand = c(0, 0.02)
  ) +
  labs(
    x = expression(bold(N[1])),
    y = "Power"
  ) +
  theme_classic(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    axis.text = element_text(color = "black"),
    legend.title = element_blank(),
    legend.position = "right"
  )

# ------------------------
# 5C. Final Combined Plot
# ------------------------
final_plot <- (panel_A / panel_B) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 28))

# =======================================================
# 6. Output
# =======================================================

ggsave("individual_panel_A.png", panel_A, width = 16, height = 7, dpi = 600)
ggsave("individual_panel_A.pdf", panel_A, width = 16, height = 7, device = cairo_pdf)

ggsave("individual_panel_B.png", panel_B, width = 16, height = 7, dpi = 600)
ggsave("individual_panel_B.pdf", panel_B, width = 16, height = 7, device = cairo_pdf)

ggsave("final_combined_analysis.png", final_plot, width = 16, height = 14, dpi = 600)
ggsave("final_combined_analysis.pdf", final_plot, width = 16, height = 14, device = cairo_pdf)
