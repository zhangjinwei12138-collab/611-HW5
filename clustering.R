suppressPackageStartupMessages({
  library(cluster)
  library(ggplot2)
  library(dplyr)
})

# =========================================
# 0) Output directory
# =========================================
OUT_DIR <- "C:/Users/22103/Desktop/611/hw5/"  
dir.create(file.path(OUT_DIR, "figs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "artifacts"), recursive = TRUE, showWarnings = FALSE)
cat("Output directories created at:", OUT_DIR, "\n")

# -------------------------------
# 1) Data generation function
# -------------------------------
generate_hypercube_clusters <- function(n, k, side_length, noise_sd = 1.0) {
  centers <- diag(side_length, n, n)
  data_list <- lapply(1:n, function(i) {
    matrix(rnorm(k * n, mean = 0, sd = noise_sd), nrow = k, ncol = n) +
      matrix(rep(centers[i, ], each = k), nrow = k)
  })
  data <- do.call(rbind, data_list)
  labels <- rep(1:n, each = k)
  return(list(x = data, y = labels))
}

# -------------------------------
# 2) Simulation parameters
# -------------------------------
dimensions   <- c(6, 5, 4, 3, 2)
side_lengths <- 10:1
k <- 100
noise_sd <- 1.0

# Reproducibility settings
GLOBAL_SEED <- 123
Kmax_offset <- 2
B_boot      <- 20
KM_nstart   <- 20
KM_itermax  <- 50

# -------------------------------
# 3) Main loop
# -------------------------------
results <- data.frame()

for (n in dimensions) {
  for (L in side_lengths) {
    set.seed(GLOBAL_SEED)
    dat <- generate_hypercube_clusters(n, k, L, noise_sd)
    
    gap_result <- clusGap(
      dat$x,
      FUN = kmeans,
      K.max   = n + Kmax_offset,  # still n + 2
      B       = B_boot,           # still 20 bootstrap samples
      nstart  = KM_nstart,        # still 20
      iter.max= KM_itermax        # still 50 iterations
    )
    
    best_k <- maxSE(gap_result$Tab[, "gap"],
                    gap_result$Tab[, "SE.sim"],
                    method = "Tibs2001SEmax")
    
    results <- rbind(results, data.frame(
      Dimension = n, SideLength = L, EstimatedK = best_k, TrueK = n
    ))
    cat("Finished n =", n, "L =", L, "=> K^ =", best_k, "\n")
  }
}

# Save raw results (for reporting)
write.csv(results, file.path(OUT_DIR, "artifacts", "task1_gap_results.csv"), row.names = FALSE)

# -------------------------------
# 4) Basic summary statistics 
# -------------------------------
# For each (Dimension, L), mark whether the estimate equals the true K
by_L <- results %>%
  mutate(Success = as.integer(EstimatedK == TrueK)) %>%
  group_by(Dimension, SideLength) %>%
  summarise(
    TrueK = first(TrueK),
    SuccessRate = mean(Success),         # Single seed → only 0 or 1
    EstimatedK = first(EstimatedK),
    .groups = "drop"
  ) %>%
  arrange(Dimension, desc(SideLength))   # Start with well-separated clusters → move closer

write.csv(by_L, file.path(OUT_DIR, "artifacts", "task1_gap_success_rates.csv"), row.names = FALSE)

# For each dimension, find the first L where EstimatedK < TrueK
summary_by_dim <- by_L %>%
  group_by(Dimension) %>%
  arrange(desc(SideLength), .by_group = TRUE) %>%
  summarise(
    TrueK = first(TrueK),
    # Range of L where EstimatedK == TrueK
    plateau_L_max = {
      tmp <- SideLength[EstimatedK == TrueK]
      if (length(tmp)) max(tmp) else NA_integer_
    },
    plateau_L_min = {
      tmp <- SideLength[EstimatedK == TrueK]
      if (length(tmp)) min(tmp) else NA_integer_
    },
    # Failure point: largest L below the correct range where EstimatedK < TrueK
    fail_start_L = {
      if (!length(SideLength)) NA_integer_ else {
        drops <- by_L %>%
          filter(Dimension == first(Dimension),
                 SideLength < ifelse(length(SideLength[EstimatedK == TrueK])>0,
                                     max(SideLength[EstimatedK == TrueK]), Inf),
                 EstimatedK < TrueK) %>%
          pull(SideLength)
        if (length(drops)) max(drops) else NA_integer_
      }
    },
    plateau_range = ifelse(
      is.na(plateau_L_min) | is.na(plateau_L_max),
      NA_character_,
      sprintf("[%d … %d]", plateau_L_min, plateau_L_max)
    ),
    .groups = "drop"
  ) %>% arrange(desc(Dimension))

write.csv(summary_by_dim, file.path(OUT_DIR, "artifacts", "task1_gap_summary_by_dimension.csv"), row.names = FALSE)

# -------------------------------
# 5) Visualization 
# -------------------------------
p <- ggplot(results, aes(x = SideLength, y = EstimatedK, color = as.factor(Dimension))) +
  geom_line(linewidth = 1.2) +
  geom_point() +
  geom_hline(aes(yintercept = Dimension, color = as.factor(Dimension)),
             linetype = "dashed", linewidth = 0.7) +
  scale_x_reverse() +
  labs(
    title = "Estimated Number of Clusters vs. Side Length",
    subtitle = sprintf("B=%d; nstart=%d; iter.max=%d; K.max=n+%d; seed=%d",
                       B_boot, KM_nstart, KM_itermax, Kmax_offset, GLOBAL_SEED),
    x = "Side Length (cluster separation)",
    y = "Estimated number of clusters (Gap Statistic)",
    color = "Dimension (True # Clusters)"
  ) +
  theme_minimal(base_size = 14)

print(p)
ggsave(file.path(OUT_DIR, "figs", "task1_gap_vs_side_length.png"),
       p, width = 9, height = 6, dpi = 300)

# -------------------------------
# 6) Save run settings and session information
# -------------------------------
settings_txt <- sprintf(
  paste(
    "Task 1 run settings:",
    "OUT_DIR: %s",
    "Dimensions: %s",
    "Side lengths: %s",
    "k per cluster: %d",
    "noise_sd: %.2f",
    "seed: %d",
    "Gap: B=%d, K.max=n+%d",
    "kmeans: nstart=%d, iter.max=%d",
    sep = "\n"
  ),
  OUT_DIR,
  paste(dimensions, collapse = ", "),
  paste(rev(side_lengths), collapse = " .. "),
  k, noise_sd, GLOBAL_SEED,
  B_boot, Kmax_offset, KM_nstart, KM_itermax
)
writeLines(settings_txt, file.path(OUT_DIR, "artifacts", "task1_run_settings.txt"))
capture.output(sessionInfo(), file = file.path(OUT_DIR, "artifacts", "task1_sessionInfo.txt"))

cat("\nArtifacts written to:\n",
    file.path(OUT_DIR, "artifacts", "task1_gap_results.csv"), "\n",
    file.path(OUT_DIR, "artifacts", "task1_gap_success_rates.csv"), "\n",
    file.path(OUT_DIR, "artifacts", "task1_gap_summary_by_dimension.csv"), "\n",
    file.path(OUT_DIR, "artifacts", "task1_run_settings.txt"), "\n",
    file.path(OUT_DIR, "artifacts", "task1_sessionInfo.txt"), "\n",
    "Figure:\n",
    file.path(OUT_DIR, "figs", "task1_gap_vs_side_length.png"), "\n", sep = "")


# ================================
# BIOS611: Clustering – Task 2
# ================================
library(plotly)
library(ggplot2)
library(cluster)    # clusGap
library(stats)

# -------------------------------
# 1) Data generator
# -------------------------------
generate_shell_clusters <- function(n_shells, k_per_shell, max_radius, noise_sd = 0.1,
                                    inner_radius = 0.5) {
  if (max_radius <= inner_radius) {
    radii <- rep(inner_radius, n_shells)
  } else {
    radii <- seq(inner_radius, max_radius, length.out = n_shells)
  }
  X <- matrix(NA_real_, nrow = n_shells * k_per_shell, ncol = 3)
  labs <- integer(n_shells * k_per_shell)
  
  row_idx <- 1
  for (s in seq_len(n_shells)) {
    r0 <- radii[s]
    for (j in 1:k_per_shell) {
      v <- rnorm(3); v <- v / sqrt(sum(v^2))
      r <- max(r0 + rnorm(1, 0, noise_sd), 0)
      X[row_idx, ] <- r * v
      labs[row_idx] <- s
      row_idx <- row_idx + 1
    }
  }
  list(x = X, y = labs, radii = radii)
}

# -------------------------------
# 2) Quick interactive check (one example dataset)
# -------------------------------
set.seed(1)
demo <- generate_shell_clusters(n_shells = 4, k_per_shell = 100, max_radius = 6, noise_sd = 0.1)
plot_ly(
  x = demo$x[,1], y = demo$x[,2], z = demo$x[,3],
  type = "scatter3d", mode = "markers",
  marker = list(size = 3),
  color = as.factor(demo$y)
)

# -------------------------------
# 3) Spectral clustering wrapper for clusGap
# -------------------------------
spectral_wrapper_factory <- function(d_threshold = 1) {
  function(x, k) {
    n <- nrow(x)
    dmat <- as.matrix(dist(x, method = "euclidean"))
    A <- (dmat < d_threshold) * 1
    diag(A) <- 0
    
    deg <- rowSums(A)
    deg[deg == 0] <- 1e-8
    Dm12 <- diag(1 / sqrt(deg), n, n)
    L <- diag(deg, n, n) - A
    Lsym <- Dm12 %*% L %*% Dm12
    
    eig <- eigen((Lsym + t(Lsym)) / 2, symmetric = TRUE)
    U <- eig$vectors[, (n - k + 1):n, drop = FALSE]
    
    row_norms <- sqrt(rowSums(U^2))
    row_norms[row_norms == 0] <- 1
    U_norm <- U / row_norms
    
    km <- kmeans(U_norm, centers = k, nstart = 20, iter.max = 50)
    list(cluster = km$cluster)
  }
}
spectral_FUN <- spectral_wrapper_factory(d_threshold = 1)

# -------------------------------
# 4) Simulation loop with clusGap
# -------------------------------
n_shells    <- 4
k_per_shell <- 100
noise_sd2   <- 0.1
max_radii   <- 10:0  # 10 down to 0

sim_results <- data.frame()

set.seed(123)
for (R in max_radii) {
  dat <- generate_shell_clusters(n_shells, k_per_shell, R, noise_sd2)
  gap <- clusGap(
    dat$x,
    FUN = spectral_FUN,
    K.max = 8,
    B = 20
  )
  best_k <- maxSE(gap$Tab[, "gap"], gap$Tab[, "SE.sim"], method = "Tibs2001SEmax")
  sim_results <- rbind(sim_results, data.frame(MaxRadius = R, EstimatedK = best_k))
  cat("Finished max_radius =", R, " => Estimated K =", best_k, "\n")
}


write.csv(sim_results, file.path(OUT_DIR, "artifacts", "task2_spectral_gap_results.csv"), row.names = FALSE)

# -------------------------------
# 5) Visualization
# -------------------------------
p2 <- ggplot(sim_results, aes(x = MaxRadius, y = EstimatedK)) +
  geom_line(linewidth = 1.2) +
  geom_point() +
  geom_hline(yintercept = 4, linetype = "dashed") +
  scale_x_continuous(breaks = max_radii) +
  labs(
    title = "Spectral Clustering (Gap Statistic) on Concentric Shells",
    subtitle = "d_threshold = 1; n_shells = 4; k_per_shell = 100; noise_sd = 0.1",
    x = "max_radius",
    y = "Estimated number of clusters"
  ) +
  theme_minimal(base_size = 14)

print(p2)
ggsave(file.path(OUT_DIR, "figs", "task2_gap_vs_maxradius.png"),
       p2, width = 9, height = 6, dpi = 300)
