# Create bound for suspicious features.

boundFind <- function(
    ind_suspicious,
    Diag_split,
    X_split,
    m,
    gap1,
    gap2,
    eps,
    maxdimension,
    maxscale
) {
  if (!is.numeric(m) || length(m) != 1L || is.na(m) ||
      m < 2L || m != as.integer(m)) {
    stop("m must be one integer-like value greater than or equal to 2.")
  }

  if (!is.numeric(eps) || length(eps) != 1L ||
      !is.finite(eps) || eps <= 0) {
    stop("eps must be one finite positive number.")
  }

  if (!is.numeric(maxdimension) || length(maxdimension) != 1L ||
      is.na(maxdimension) || maxdimension < 1L ||
      maxdimension != as.integer(maxdimension)) {
    stop("maxdimension must be at least 1 because boundFind recovers H1 cycles.")
  }

  if (!is.numeric(maxscale) || length(maxscale) != 1L ||
      !is.finite(maxscale) || maxscale <= 0) {
    stop("maxscale must be one finite positive number.")
  }

  n_blocks <- as.integer(m - 1L)

  if (length(gap1) < m || length(gap2) < m) {
    stop("gap1 and gap2 must each contain at least m entries.")
  }

  expected_dim <- c(n_blocks, n_blocks)
  if (!identical(dim(ind_suspicious), expected_dim) ||
      !identical(dim(Diag_split), expected_dim) ||
      !identical(dim(X_split), expected_dim)) {
    stop("ind_suspicious, Diag_split, and X_split must be (m - 1) x (m - 1) objects.")
  }

  bound <- matrix(
    vector("list", n_blocks * n_blocks),
    nrow = n_blocks,
    ncol = n_blocks
  )

  Suspicious <- matrix(
    vector("list", n_blocks * n_blocks),
    nrow = n_blocks,
    ncol = n_blocks
  )

  empty_bound <- function() {
    matrix(
      Inf,
      nrow = 2L,
      ncol = 2L,
      dimnames = list(NULL, c("x", "y"))
    )
  }

  for (i in seq_len(n_blocks)) {
    for (j in seq_len(n_blocks)) {
      ind_suspicious_ij <- ind_suspicious[[i, j]]
      n_local_features <- length(ind_suspicious_ij)

      if (n_local_features == 0L) {
        next
      }

      X_split_ij <- X_split[[i, j]]
      Diag_split_ij <- Diag_split[[i, j]]

      if (is.null(X_split_ij) || NCOL(X_split_ij) != 2L) {
        stop(sprintf(
          "X_split[[%d, %d]] must be a non-NULL object with exactly two columns.",
          i, j
        ))
      }

      if (is.null(Diag_split_ij) ||
          is.null(Diag_split_ij[["cycleLocation"]])) {
        stop(sprintf(
          "Diag_split[[%d, %d]] does not contain cycleLocation.",
          i, j
        ))
      }

      gap3 <- c(gap1[i], gap1[i + 1L])
      gap4 <- c(gap2[j], gap2[j + 1L])
      lengthX <- NROW(X_split_ij)

      # Preserve the original nested structure expected by Projected_Merge().
      bound[[i, j]] <- matrix(
        vector("list", n_local_features),
        nrow = 1L,
        ncol = n_local_features
      )

      Suspicious[[i, j]] <- matrix(
        vector("list", n_local_features),
        nrow = 1L,
        ncol = n_local_features
      )

      for (ij in seq_along(ind_suspicious_ij)) {
        diag_id <- ind_suspicious_ij[[ij]]

        if (is.na(diag_id) || diag_id < 1L ||
            diag_id > length(Diag_split_ij[["cycleLocation"]])) {
          stop(sprintf(
            "Invalid diagram index for suspicious feature %d in block (%d, %d).",
            ij, i, j
          ))
        }

        data1 <- Diag_split_ij[["cycleLocation"]][[diag_id]]
        IN <- as.integer(as.vector(data1))
        IN <- IN[!is.na(IN)]

        # Supplemental boundary vertices follow the real observations.
        boundary_codes <- unique(IN[IN > lengthX] - lengthX)

        # Retain unique real-data indices and preserve matrix dimensions.
        data_indices <- unique(IN[IN >= 1L & IN <= lengthX])
        data2 <- X_split_ij[data_indices, , drop = FALSE]
        colnames(data2) <- c("x", "y")

        local_width <- abs(gap3[2L] - gap3[1L])
        lengthBound <- max(
          2L,
          floor(lengthX / 5) + 2L,
          floor(local_width / eps)
        )
        lengthBound <- as.integer(lengthBound)

        bound1 <- matrix(
          numeric(0),
          nrow = 0L,
          ncol = 2L,
          dimnames = list(NULL, c("x", "y"))
        )

        # Boundary code 1.
        if (1L %in% boundary_codes) {
          if (i != 1L) {
            bound1 <- rbind(
              bound1,
              cbind(
                rep(gap3[1L], lengthBound),
                seq(gap4[1L], gap4[2L], length.out = lengthBound)
              )
            )
          } else {
            bound1 <- rbind(
              bound1,
              cbind(
                rep(gap3[2L], lengthBound),
                seq(gap4[1L], gap4[2L], length.out = lengthBound)
              )
            )
          }
        }

        # Boundary code 2.
        if (2L %in% boundary_codes) {
          if (i > 1L && i < n_blocks) {
            bound1 <- rbind(
              bound1,
              cbind(
                rep(gap3[2L], lengthBound),
                seq(gap4[1L], gap4[2L], length.out = lengthBound)
              )
            )
          } else if (j > 1L) {
            bound1 <- rbind(
              bound1,
              cbind(
                seq(gap3[1L], gap3[2L], length.out = lengthBound),
                rep(gap4[1L], lengthBound)
              )
            )
          } else {
            bound1 <- rbind(
              bound1,
              cbind(
                seq(gap3[1L], gap3[2L], length.out = lengthBound),
                rep(gap4[2L], lengthBound)
              )
            )
          }
        }

        # Boundary code 3.
        if (3L %in% boundary_codes) {
          if (i == 1L || i == n_blocks || j == 1L) {
            bound1 <- rbind(
              bound1,
              cbind(
                seq(gap3[1L], gap3[2L], length.out = lengthBound),
                rep(gap4[2L], lengthBound)
              )
            )
          } else {
            bound1 <- rbind(
              bound1,
              cbind(
                seq(gap3[1L], gap3[2L], length.out = lengthBound),
                rep(gap4[1L], lengthBound)
              )
            )
          }
        }

        # Boundary code 4.
        if (4L %in% boundary_codes) {
          bound1 <- rbind(
            bound1,
            cbind(
              seq(gap3[1L], gap3[2L], length.out = lengthBound),
              rep(gap4[2L], lengthBound)
            )
          )
        }

        bound1 <- unique(bound1)
        colnames(bound1) <- c("x", "y")

        rips_input <- rbind(bound1, data2)

        # Fewer than three points cannot support a one-dimensional cycle.
        if (NROW(rips_input) < 3L) {
          bound[[i, j]][[ij]] <- empty_bound()
          next
        }

        DiagProj <- ripsDiag(
          X = rips_input,
          maxdimension = maxdimension,
          maxscale = maxscale,
          library = "Dionysus",
          location = TRUE,
          printProgress = FALSE
        )

        h1_indices <- which(DiagProj$diagram[, 1] == 1)

        if (length(h1_indices) == 0L) {
          bound[[i, j]][[ij]] <- empty_bound()
          next
        }

        lifetimes <-
          DiagProj$diagram[h1_indices, 3] -
          DiagProj$diagram[h1_indices, 2]
        selected_id <- h1_indices[[which.max(lifetimes)]]
        selected_cycle <- DiagProj$cycleLocation[[selected_id]]

        if (is.null(selected_cycle) || length(selected_cycle) == 0L) {
          bound[[i, j]][[ij]] <- empty_bound()
          next
        }

        Suspicious[[i, j]][[ij]] <- selected_cycle

        cycle_points <- unique(matrix(
          as.vector(selected_cycle),
          ncol = 2L,
          byrow = FALSE
        ))
        colnames(cycle_points) <- c("x", "y")

        on_boundary <-
          cycle_points[, 1] == gap3[1L] |
          cycle_points[, 1] == gap3[2L] |
          cycle_points[, 2] == gap4[1L] |
          cycle_points[, 2] == gap4[2L]

        recovered_bound <- cycle_points[
          on_boundary,
          ,
          drop = FALSE
        ]

        if (NROW(recovered_bound) == 0L) {
          recovered_bound <- empty_bound()
        }

        bound[[i, j]][[ij]] <- recovered_bound
      }
    }
  }

  return(list(
    bound = bound,
    Suspicious = Suspicious
  ))
}