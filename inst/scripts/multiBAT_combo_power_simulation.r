# =======================================================
# 1. Setup and parameters
# =======================================================

library(ggplot2)
library(plyr)
library(mvtnorm)
library(dplyr)
library(tidyr)
library(patchwork)
library(rstatix)
library(grid)

beta11 <- 1
beta12 <- -1
beta21 <- 1
beta22 <- -1

LD1 <- 0.6
LD2 <- 0.6
vare <- 100
alpha <- 0.05
m <- 2

QUICK_TEST <- FALSE

if (QUICK_TEST) {
  NREP_LINE <- 200
  NREP_INNER <- 10
  NREP_OUTER <- 5
} else {
  NREP_LINE <- 10000
  NREP_INNER <- 100
  NREP_OUTER <- 50
}

nrep <- NREP_LINE

nvec1 <- seq(100, 1000, by = 100)
n2 <- 1000

BASE_SIZE <- 18
TITLE_SIZE <- 22
TAG_SIZE <- 28
TICK_SIZE <- BASE_SIZE + 2
AXIS_TITLE_SIZE <- BASE_SIZE + 2

font_family <- "Arial"

DATE_TAG <- format(Sys.Date(), "%m%d")
OUTPUT_DIR <- "."
OUTPUT_FILE_STUB <- "Figure2"

OUTPUT_WIDTH <- 16
OUTPUT_HEIGHT <- 14

FPR_YMIN <- 0.00
FPR_YMAX <- 0.15
FPR_YBREAKS <- c(0.00, 0.05, 0.10, 0.15)

POWER_YMAX <- 1.20
POWER_BALANCED_YMIN <- 0.25
POWER_IMBALANCED_YMIN <- 0.00

POWER_BALANCED_YBREAKS <- c(0.25, 0.50, 0.75, 1.00)
POWER_IMBALANCED_YBREAKS <- c(0.00, 0.25, 0.50, 0.75, 1.00)

SCATTER_YMIN <- 0.00
SCATTER_YMAX <- 1.00
SCATTER_YBREAKS <- seq(0, 1.00, by = 0.20)

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

# =======================================================
# 2. Functions
# =======================================================

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

get_z_se <- function(X, y) {
  res <- apply(X, 2, function(x) {
    fit <- summary(lm(y ~ x))$coefficients[2, 1:2]
    c(z = fit[1] / fit[2], se = fit[2])
  })
  
  list(z = res[1, ], se = res[2, ])
}

get_z_only <- function(X, y) {
  apply(X, 2, function(x) {
    fit <- summary(lm(y ~ x))$coefficients[2, 1:2]
    fit[1] / fit[2]
  })
}

p_to_stars <- function(p) {
  ifelse(
    is.na(p), "ns",
    ifelse(
      p <= 0.0001, "****",
      ifelse(
        p <= 0.001, "***",
        ifelse(
          p <= 0.01, "**",
          ifelse(p <= 0.05, "*", "ns")
        )
      )
    )
  )
}

simulate_power <- function(n1, n2 = 1000, nrep_inner = 100, nrep_outer = 50,
                           beta1_val = beta11, beta2_val = beta12) {
  res.box <- data.frame()
  
  R1 <- matrix(c(1, LD1, LD1, 1), nrow = 2)
  R2 <- matrix(c(1, LD2, LD2, 1), nrow = 2)
  
  eig1 <- eigen(R1, symmetric = TRUE)
  eig2 <- eigen(R2, symmetric = TRUE)
  
  for (outer in 1:nrep_outer) {
    P_mbat1 <- numeric(nrep_inner)
    P_mbat2 <- numeric(nrep_inner)
    P_single_meta <- numeric(nrep_inner)
    P_mbat_combined <- numeric(nrep_inner)
    P_multibat_combo <- numeric(nrep_inner)
    
    for (k in 1:nrep_inner) {
      X1 <- rmvnorm(n1, mean = c(0, 0), sigma = R1)
      y1 <- X1 %*% c(beta1_val, beta2_val) + rnorm(n1, 0, sqrt(vare))
      zse1 <- get_z_se(X1, y1)
      z1 <- zse1$z
      
      X2 <- rmvnorm(n2, mean = c(0, 0), sigma = R2)
      y2 <- X2 %*% c(beta1_val, beta2_val) + rnorm(n2, 0, sqrt(vare))
      zse2 <- get_z_se(X2, y2)
      z2 <- zse2$z
      
      UPz1 <- crossprod(eig1$vectors, z1)
      chisq1 <- crossprod(UPz1, (1 / eig1$values) * UPz1)
      P_mbat1[k] <- pchisq(chisq1, df = m, lower.tail = FALSE)
      
      UPz2 <- crossprod(eig2$vectors, z2)
      chisq2 <- crossprod(UPz2, (1 / eig2$values) * UPz2)
      P_mbat2[k] <- pchisq(chisq2, df = m, lower.tail = FALSE)
      
      w_pop1_meta <- n1 / (n1 + n2)
      w_pop2_meta <- n2 / (n1 + n2)
      
      P_single_meta[k] <- ACAT_weighted(
        p = c(P_mbat1[k], P_mbat2[k]),
        w = c(w_pop1_meta, w_pop2_meta)
      )
      
      z_combined <- c(z1, z2)
      Sigma <- matrix(0, 4, 4)
      Sigma[1:2, 1:2] <- R1
      Sigma[3:4, 3:4] <- R2
      
      eig_comb <- eigen(Sigma, symmetric = TRUE)
      UPz_comb <- crossprod(eig_comb$vectors, z_combined)
      chisq_comb <- crossprod(UPz_comb, (1 / eig_comb$values) * UPz_comb)
      P_mbat_combined[k] <- pchisq(chisq_comb, df = 4, lower.tail = FALSE)
      
      w_joint <- 1 - abs(n1 - n2) / (n1 + n2)
      w_joint <- max(0, min(1, w_joint))
      w_pop1 <- (1 - w_joint) * n1 / (n1 + n2)
      w_pop2 <- (1 - w_joint) * n2 / (n1 + n2)
      
      P_multibat_combo[k] <- ACAT_weighted(
        p = c(P_mbat1[k], P_mbat2[k], P_mbat_combined[k]),
        w = c(w_pop1, w_pop2, w_joint)
      )
    }
    
    res.box <- rbind(
      res.box,
      data.frame(N1 = n1, Replicate = outer, Method = "Population 1", Power = mean(P_mbat1 < alpha)),
      data.frame(N1 = n1, Replicate = outer, Method = "Population 2", Power = mean(P_mbat2 < alpha)),
      data.frame(N1 = n1, Replicate = outer, Method = "Single-ancestry meta", Power = mean(P_single_meta < alpha)),
      data.frame(N1 = n1, Replicate = outer, Method = "multiBAT", Power = mean(P_mbat_combined < alpha)),
      data.frame(N1 = n1, Replicate = outer, Method = "multiBAT-combo", Power = mean(P_multibat_combo < alpha))
    )
  }
  
  return(res.box)
}

prepare_sig_data <- function(dat, compare_multibat = TRUE) {
  dat_filtered <- dat %>%
    filter(Method %in% c(
      "Population 2",
      "Single-ancestry meta",
      "multiBAT",
      "multiBAT-combo"
    )) %>%
    mutate(Method = factor(Method, levels = method_levels))
  
  if (compare_multibat) {
    comparisons <- data.frame(
      group1 = c(
        "Population 2",
        "Single-ancestry meta",
        "multiBAT"
      ),
      group2 = c(
        "multiBAT-combo",
        "multiBAT-combo",
        "multiBAT-combo"
      ),
      y = c(0.98, 1.06, 1.14),
      stringsAsFactors = FALSE
    )
  } else {
    comparisons <- data.frame(
      group1 = c(
        "Population 2",
        "Single-ancestry meta"
      ),
      group2 = c(
        "multiBAT-combo",
        "multiBAT-combo"
      ),
      y = c(1.02, 1.10),
      stringsAsFactors = FALSE
    )
  }
  
  calc_p <- function(g1, g2) {
    x <- dat_filtered %>%
      filter(Method == g1) %>%
      arrange(Replicate) %>%
      pull(Power)
    
    y <- dat_filtered %>%
      filter(Method == g2) %>%
      arrange(Replicate) %>%
      pull(Power)
    
    wilcox.test(x, y, paired = TRUE, exact = FALSE)$p.value
  }
  
  comparisons$p.raw <- mapply(calc_p, comparisons$group1, comparisons$group2)
  comparisons$p.adj <- p.adjust(comparisons$p.raw, method = "BH")
  comparisons$p.stars <- p_to_stars(comparisons$p.adj)
  
  tag_data <- comparisons %>%
    mutate(
      start = match(group1, method_levels),
      end = match(group2, method_levels),
      h = 0.015
    )
  
  list(tag_data = tag_data, max_y = POWER_YMAX)
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

# =======================================================
# 3. Data generation
# =======================================================

res.theory <- data.frame()

for (n1 in nvec1) {
  R1 <- matrix(c(1, LD1, LD1, 1), nrow = 2)
  R2 <- matrix(c(1, LD2, LD2, 1), nrow = 2)
  
  varg1 <- beta11^2 + beta12^2 + 2 * LD1 * beta11 * beta12
  vary1 <- varg1 + vare
  b1_1 <- beta11 + LD1 * beta12
  b2_1 <- beta12 + LD1 * beta11
  
  mu_vec1 <- c(
    sqrt(n1 / vary1 / (1 - (b1_1^2 / vary1))) * b1_1,
    sqrt(n1 / vary1 / (1 - (b2_1^2 / vary1))) * b2_1
  )
  
  ncp1 <- as.numeric(t(mu_vec1) %*% solve(R1) %*% mu_vec1)
  
  varg2 <- beta21^2 + beta22^2 + 2 * LD2 * beta21 * beta22
  vary2 <- varg2 + vare
  b1_2 <- beta21 + LD2 * beta22
  b2_2 <- beta22 + LD2 * beta21
  
  mu_vec2 <- c(
    sqrt(n2 / vary2 / (1 - (b1_2^2 / vary2))) * b1_2,
    sqrt(n2 / vary2 / (1 - (b2_2^2 / vary2))) * b2_2
  )
  
  ncp2 <- as.numeric(t(mu_vec2) %*% solve(R2) %*% mu_vec2)
  ncp_comb <- ncp1 + ncp2
  
  power1 <- pchisq(
    qchisq(1 - alpha, df = 2),
    df = 2,
    ncp = ncp1,
    lower.tail = FALSE
  )
  
  power2 <- pchisq(
    qchisq(1 - alpha, df = 2),
    df = 2,
    ncp = ncp2,
    lower.tail = FALSE
  )
  
  power_comb <- pchisq(
    qchisq(1 - alpha, df = 4),
    df = 4,
    ncp = ncp_comb,
    lower.tail = FALSE
  )
  
  res.theory <- rbind(
    res.theory,
    data.frame(N1 = n1, Method = "Population 1", Power = power1, Type = "Theory"),
    data.frame(N1 = n1, Method = "Population 2", Power = power2, Type = "Theory"),
    data.frame(N1 = n1, Method = "multiBAT", Power = power_comb, Type = "Theory")
  )
}

res_null <- simulate_power(
  1000,
  1000,
  nrep_inner = NREP_INNER,
  nrep_outer = NREP_OUTER,
  beta1_val = 0,
  beta2_val = 0
)

res_200 <- simulate_power(
  200,
  1000,
  nrep_inner = NREP_INNER,
  nrep_outer = NREP_OUTER
)

res_1000 <- simulate_power(
  1000,
  1000,
  nrep_inner = NREP_INNER,
  nrep_outer = NREP_OUTER
)

res.emp.line <- data.frame()

R2 <- matrix(c(1, LD2, LD2, 1), nrow = 2)
eig2 <- eigen(R2, symmetric = TRUE)

P2_fixed <- numeric(NREP_INNER * NREP_OUTER)

for (k in seq_along(P2_fixed)) {
  X2 <- rmvnorm(n2, mean = c(0, 0), sigma = R2)
  y2 <- X2 %*% c(beta21, beta22) + rnorm(n2, 0, sqrt(vare))
  z2 <- get_z_only(X2, y2)
  
  UPz2 <- crossprod(eig2$vectors, z2)
  P2_fixed[k] <- pchisq(sum(UPz2^2 / eig2$values), df = 2, lower.tail = FALSE)
}

power2_fixed <- mean(P2_fixed < alpha)

for (n1 in nvec1) {
  R1 <- matrix(c(1, LD1, LD1, 1), nrow = 2)
  R2 <- matrix(c(1, LD2, LD2, 1), nrow = 2)
  
  eig1 <- eigen(R1, symmetric = TRUE)
  eig2 <- eigen(R2, symmetric = TRUE)
  
  P1 <- numeric(nrep)
  Pmeta <- numeric(nrep)
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
    P2_current <- pchisq(sum(UPz2^2 / eig2$values), df = 2, lower.tail = FALSE)
    
    w_pop1_meta <- n1 / (n1 + n2)
    w_pop2_meta <- n2 / (n1 + n2)
    
    Pmeta[k] <- ACAT_weighted(
      p = c(P1[k], P2_current),
      w = c(w_pop1_meta, w_pop2_meta)
    )
    
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
      p = c(P1[k], P2_current, Pc[k]),
      w = c(w_pop1, w_pop2, w_joint)
    )
  }
  
  res.emp.line <- rbind(
    res.emp.line,
    data.frame(N1 = n1, Method = "Population 1", Power = mean(P1 < alpha), Type = "Simulation"),
    data.frame(N1 = n1, Method = "Population 2", Power = power2_fixed, Type = "Simulation"),
    data.frame(N1 = n1, Method = "Single-ancestry meta", Power = mean(Pmeta < alpha), Type = "Simulation"),
    data.frame(N1 = n1, Method = "multiBAT", Power = mean(Pc < alpha), Type = "Simulation"),
    data.frame(N1 = n1, Method = "multiBAT-combo", Power = mean(Pcombo < alpha), Type = "Simulation")
  )
}

# =======================================================
# 4. Significance calculation
# =======================================================

sig_200 <- prepare_sig_data(
  res_200,
  compare_multibat = TRUE
)

sig_1000 <- prepare_sig_data(
  res_1000,
  compare_multibat = FALSE
)

# =======================================================
# 5. Plotting
# =======================================================

tag_theme <- theme(
  plot.tag = element_text(
    face = "bold",
    size = TAG_SIZE,
    family = font_family
  )
)

make_boxplot <- function(dat, title_text = NULL, y_label = "Power",
                         drop_y = FALSE, add_hline = FALSE,
                         sig_data = NULL, y_min = 0.00,
                         y_max = POWER_YMAX,
                         y_breaks = POWER_IMBALANCED_YBREAKS) {
  dat$Method <- factor(dat$Method, levels = method_levels)
  
  p <- ggplot(dat, aes(x = Method, y = Power, fill = Method)) +
    geom_boxplot(linewidth = 1.1, alpha = 0.9, show.legend = FALSE) +
    scale_fill_manual(
      values = method_colors,
      breaks = method_levels,
      labels = method_levels,
      drop = FALSE,
      guide = "none"
    ) +
    theme_bw(base_size = BASE_SIZE, base_family = font_family) +
    labs(title = title_text, y = y_label, x = NULL) +
    scale_y_continuous(
      breaks = y_breaks,
      labels = sprintf("%.2f", y_breaks)
    ) +
    coord_cartesian(ylim = c(y_min, y_max)) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(face = "bold", size = AXIS_TITLE_SIZE),
      axis.text = element_text(face = "plain", size = TICK_SIZE, colour = "black"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_line(colour = "black", linewidth = 1.0),
      plot.title = element_text(hjust = 0.5, face = "bold", size = TITLE_SIZE),
      plot.subtitle = element_text(face = "bold", size = BASE_SIZE),
      strip.text = element_text(face = "bold", size = BASE_SIZE),
      legend.position = "none"
    )
  
  if (add_hline) {
    p <- p + geom_hline(
      yintercept = 0.05,
      linetype = "dashed",
      colour = "red",
      linewidth = 1.1
    )
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
          label = p.stars
        ),
        inherit.aes = FALSE,
        size = 6,
        fontface = "plain",
        family = font_family
      )
  }
  
  p
}

p_null <- make_boxplot(
  res_null,
  title_text = " ",
  y_label = "False positive rate",
  add_hline = TRUE,
  y_min = FPR_YMIN,
  y_max = FPR_YMAX,
  y_breaks = FPR_YBREAKS
) +
  labs(tag = "A") +
  tag_theme

p_200 <- make_boxplot(
  res_200,
  title_text = "Imbalanced samples",
  y_label = "Power",
  drop_y = FALSE,
  sig_data = sig_200,
  y_min = POWER_IMBALANCED_YMIN,
  y_max = POWER_YMAX,
  y_breaks = POWER_IMBALANCED_YBREAKS
) +
  labs(tag = "C") +
  tag_theme


p_1000 <- make_boxplot(
  res_1000,
  title_text = "Balanced samples",
  y_label = "Power",
  drop_y = TRUE,
  sig_data = sig_1000,
  y_min = POWER_IMBALANCED_YMIN,
  y_max = POWER_YMAX,
  y_breaks = POWER_IMBALANCED_YBREAKS
)

combined_line <- bind_rows(res.theory, res.emp.line) %>%
  mutate(
    SampleRatio = N1 / n2,
    Method = factor(Method, levels = method_levels),
    Type = factor(Type, levels = c("Theory", "Simulation"))
  )

panel_scatter <- ggplot(
  combined_line,
  aes(
    x = SampleRatio,
    y = Power,
    colour = Method,
    group = interaction(Method, Type)
  )
) +
  geom_point(
    aes(shape = Type),
    size = 3,
    show.legend = FALSE
  ) +
  geom_smooth(
    aes(linetype = Type),
    method = "loess",
    formula = y ~ x,
    se = FALSE,
    span = 0.9,
    linewidth = 1.1,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = method_colors,
    guide = "none"
  ) +
  scale_linetype_manual(
    values = c("Theory" = "solid", "Simulation" = "44"),
    guide = "none"
  ) +
  scale_shape_manual(
    values = c("Theory" = 16, "Simulation" = 17),
    guide = "none"
  ) +
  scale_x_continuous(
    limits = c(0.1, 1.0),
    breaks = c(0.25, 0.50, 0.75, 1.00),
    labels = sprintf("%.2f", c(0.25, 0.50, 0.75, 1.00)),
    expand = c(0.01, 0.01)
  ) +
  scale_y_continuous(
    breaks = SCATTER_YBREAKS,
    labels = sprintf("%.2f", SCATTER_YBREAKS),
    expand = c(0, 0.02)
  ) +
  coord_cartesian(
    ylim = c(SCATTER_YMIN, SCATTER_YMAX)
  ) +
  labs(
    x = "Sample ratio",
    y = "Power",
    tag = "B"
  ) +
  theme_classic(base_size = BASE_SIZE, base_family = font_family) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.8
    ),
    axis.ticks = element_line(colour = "black", linewidth = 1.0),
    axis.title = element_text(face = "bold", size = AXIS_TITLE_SIZE),
    axis.text.x = element_text(face = "plain", size = TICK_SIZE, colour = "black"),
    axis.text.y = element_text(face = "plain", size = TICK_SIZE, colour = "black"),
    plot.title = element_text(face = "bold", size = TITLE_SIZE),
    plot.subtitle = element_text(face = "bold", size = BASE_SIZE),
    strip.text = element_text(face = "bold", size = BASE_SIZE),
    legend.position = "none",
    plot.tag = element_text(
      face = "bold",
      size = TAG_SIZE,
      family = font_family
    )
    
  )

# =======================================================
# 6. Manual legends
# =======================================================

source_legend_df <- data.frame(
  Source = factor(c("Theory", "Simulation"), levels = c("Theory", "Simulation")),
  y = c(0.53, 0.45)
)

source_legend_plot <- ggplot(source_legend_df) +
  geom_segment(
    aes(x = 0.00, xend = 0.95, y = y, yend = y, linetype = Source),
    colour = "black",
    linewidth = 1.3
  ) +
  geom_point(
    aes(x = 0.475, y = y, shape = Source),
    colour = "black",
    size = 3.8
  ) +
  geom_text(
    aes(x = 1.10, y = y, label = Source),
    hjust = 0,
    vjust = 0.5,
    family = font_family,
    fontface = "plain",
    size = BASE_SIZE / ggplot2::.pt
  ) +
  scale_linetype_manual(
    values = c("Theory" = "solid", "Simulation" = "44"),
    guide = "none"
  ) +
  scale_shape_manual(
    values = c("Theory" = 16, "Simulation" = 17),
    guide = "none"
  ) +
  coord_cartesian(
    xlim = c(0.00, 2.10),
    ylim = c(0.20, 0.90),
    clip = "off"
  ) +
  theme_void(base_family = font_family) +
  theme(
    plot.margin = margin(5, 5, 5, 5)
  )

method_legend_df <- data.frame(
  Method = factor(method_levels, levels = method_levels),
  y = rev(seq(1, by = 0.35, length.out = length(method_levels)))
)

method_legend_plot <- ggplot(method_legend_df, aes(x = 0, y = y)) +
  annotate(
    "text",
    x = 0.00,
    y = max(method_legend_df$y) + 0.45,
    label = "Method",
    hjust = 0,
    vjust = 0.5,
    family = font_family,
    fontface = "bold",
    size = BASE_SIZE / ggplot2::.pt
  ) +
  geom_point(
    aes(fill = Method),
    shape = 22,
    size = 6.5,
    colour = "black",
    stroke = 0.25
  ) +
  geom_text(
    aes(x = 0.28, label = Method),
    hjust = 0,
    vjust = 0.5,
    family = font_family,
    fontface = "plain",
    size = BASE_SIZE / ggplot2::.pt
  ) +
  scale_fill_manual(
    values = method_colors,
    guide = "none"
  ) +
  coord_cartesian(
    xlim = c(-0.05, 2.55),
    ylim = c(0.20, 6.60),
    clip = "off"
  ) +
  theme_void(base_family = font_family) +
  theme(
    plot.margin = margin(5, 5, 5, 5)
  )

right_legend_column <- source_legend_plot / method_legend_plot +
  plot_layout(
    heights = c(0.65, 1.35)
  )

# =======================================================
# 7. Final layout and output
# =======================================================

panel_power_pair <- (p_200 | p_1000) +
  plot_layout(widths = c(1, 1))

panel_C <-  panel_power_pair

row_gap <- plot_spacer()

layout_design <- "
ABG
DDG
CCG
"

final_plot <- wrap_plots(
  A = p_null,
  B = panel_scatter,
  D = row_gap,
  C = panel_C,
  G = right_legend_column,
  design = layout_design,
  widths = c(1.00, 2.60, 1.15),
  heights = c(1.00, 0.02, 1.00)
)

save_plot_pair(final_plot, OUTPUT_FILE_STUB, width = OUTPUT_WIDTH, height = OUTPUT_HEIGHT)
