

## ============================================
## Disjoint, non-nested circle manifolds in R
## ============================================

set.seed(123)

## ---- Parameters ----
n_per_circle <- 500

radii <- c(
  small  = 0.5,
  medium = 1.5,
  large  = 3.0
)

noise_frac <- 0   # relative noise (0 for exact manifolds)
margin     <- 0.5    # extra separation buffer

## ---- Define non-overlapping centers ----
centers <- list(
  small  = c(0, 0),
  medium = c(6, 0),
  large  = c(0, 7)
)

## ---- Sampling function ----
sample_noisy_circle <- function(n, R, center, noise_frac = 0) {
  theta <- runif(n, 0, 2 * pi)
  r_noise <- if (noise_frac > 0) {
    rnorm(n, mean = R, sd = noise_frac * R)
  } else {
    rep(R, n)
  }
  
  x <- center[1] + r_noise * cos(theta)
  y <- center[2] + r_noise * sin(theta)
  cbind(x, y)
}

## ---- Generate data ----
circle_data <- do.call(rbind, lapply(names(radii), function(name) {
  R <- radii[name]
  ctr <- centers[[name]]
  
  pts <- sample_noisy_circle(
    n = n_per_circle,
    R = R,
    center = ctr,
    noise_frac = noise_frac
  )
  
  data.frame(
    x = pts[, 1],
    y = pts[, 2],
    manifold = name,
    radius = R
  )
}))

## ---- Visualization ----
library(ggplot2)

ggplot(circle_data, aes(x = x, y = y)) +
  geom_point(size = 1.2, alpha = 0.7) +
  coord_equal() +
  theme_bw() +
  xlab("") +
  ylab("") +
  theme(
    text = element_text(size = 30),
    legend.title = element_blank(),
    legend.position = c(.85, .85),
    legend.background = element_blank(),
    legend.box.background = element_rect(color = "black")
  )



# Save the data to csv 
write.csv(circle_data,
          file = "disjoint_circle_manifolds.csv",
          row.names = FALSE)