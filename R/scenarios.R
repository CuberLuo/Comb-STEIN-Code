## Modified simulation scenarios based on manuscript Table 1.
##
## For scenarios whose largest and second-largest true utilities were separated
## by less than 0.10, the efficacy probability at each utility-maximizing dose
## was increased by the minimum amount needed to create a 0.10 gap. These
## changes make selection probabilities interpretable at practical simulation
## sizes while leaving all other efficacy cells unchanged. In S7, the toxicity
## probability at its OBDC (2,3) was also reduced from 0.38 to 0.30: the original
## value was almost equal to phi and caused frequent false safety elimination.
## The modified S7 toxicity surface remains monotone in both agents.

mk <- function(...) {
  matrix(c(...), ncol = 2, byrow = TRUE)
}

scenarios <- list(
  S1 = mk(
    0.10,0.15, 0.15,0.25, 0.30,0.30, 0.35,0.35, 0.45,0.45,
    0.20,0.30, 0.30,0.35, 0.40,0.40, 0.45,0.45, 0.50,0.40,
    0.30,0.75, 0.35,0.60, 0.50,0.45, 0.60,0.40, 0.65,0.35
  ),
  S2 = mk(
    0.05,0.05, 0.10,0.20, 0.15,0.40, 0.20,0.767, 0.30,0.70,
    0.20,0.15, 0.15,0.35, 0.20,0.45, 0.40,0.60, 0.55,0.55,
    0.15,0.30, 0.20,0.40, 0.25,0.50, 0.45,0.35, 0.60,0.25
  ),
  S3 = mk(
    0.05,0.10, 0.10,0.15, 0.15,0.25, 0.20,0.30, 0.30,0.40,
    0.15,0.20, 0.20,0.45, 0.25,0.7835, 0.30,0.60, 0.35,0.45,
    0.20,0.35, 0.30,0.70, 0.35,0.45, 0.45,0.35, 0.50,0.25
  ),
  S4 = mk(
    0.05,0.05, 0.10,0.20, 0.15,0.40, 0.20,0.7835, 0.25,0.70,
    0.20,0.15, 0.15,0.35, 0.20,0.45, 0.25,0.50, 0.30,0.70,
    0.25,0.30, 0.30,0.40, 0.35,0.50, 0.45,0.55, 0.60,0.65
  ),
  S5 = mk(
    0.10,0.10, 0.14,0.15, 0.16,0.25, 0.20,0.30, 0.23,0.45,
    0.13,0.30, 0.15,0.45, 0.25,0.8335, 0.30,0.75, 0.35,0.70,
    0.20,0.35, 0.25,0.40, 0.30,0.75, 0.40,0.50, 0.45,0.50
  ),
  S6 = mk(
    0.10,0.30, 0.15,0.767, 0.20,0.43, 0.25,0.40,
    0.15,0.35, 0.18,0.55, 0.25,0.70, 0.35,0.35,
    0.18,0.40, 0.20,0.45, 0.35,0.48, 0.40,0.37
  ),
  S7 = mk(
    0.05,0.20, 0.10,0.35, 0.25,0.60, 0.30,0.50,
    0.15,0.30, 0.24,0.60, 0.30,0.7462, 0.45,0.70,
    0.25,0.40, 0.30,0.50, 0.40,0.55, 0.50,0.45
  ),
  S8 = mk(
    0.05,0.20, 0.08,0.25, 0.12,0.35, 0.18,0.40,
    0.10,0.25, 0.14,0.40, 0.16,0.7873, 0.25,0.70,
    0.13,0.35, 0.15,0.45, 0.20,0.70, 0.35,0.75
  ),
  S9 = mk(
    0.01,0.20, 0.05,0.30, 0.10,0.40,
    0.05,0.30, 0.08,0.35, 0.20,0.65,
    0.10,0.35, 0.18,0.65, 0.25,0.7731
  ),
  S10 = mk(
    0.10,0.20, 0.15,0.30, 0.18,0.6769,
    0.15,0.30, 0.18,0.6769, 0.25,0.60,
    0.20,0.45, 0.26,0.50, 0.35,0.55
  ),
  S11 = mk(
    0.10,0.20, 0.20,0.35, 0.25,0.6835,
    0.15,0.30, 0.30,0.60, 0.40,0.80
  ),
  S12 = mk(
    0.05,0.20, 0.10,0.35, 0.15,0.667,
    0.15,0.30, 0.20,0.45, 0.25,0.60
  )
)

ncol_map <- c(
  S1 = 5, S2 = 5, S3 = 5, S4 = 5, S5 = 5, S6 = 4,
  S7 = 4, S8 = 4, S9 = 3, S10 = 3, S11 = 3, S12 = 3
)
nrow_map <- c(
  S1 = 3, S2 = 3, S3 = 3, S4 = 3, S5 = 3, S6 = 3,
  S7 = 3, S8 = 3, S9 = 3, S10 = 3, S11 = 2, S12 = 2
)

obdc <- function(...) {
  matrix(
    c(...),
    ncol = 2,
    byrow = TRUE,
    dimnames = list(NULL, c("i", "j"))
  )
}

table1.obdcs <- list(
  S1  = obdc(3,1),
  S2  = obdc(1,4, 1,5),
  S3  = obdc(2,3, 3,2),
  S4  = obdc(1,4, 1,5, 2,5),
  S5  = obdc(2,3, 2,4, 3,3),
  S6  = obdc(1,2, 2,3),
  S7  = obdc(1,3, 2,2),
  S8  = obdc(2,3, 2,4, 3,3),
  S9  = obdc(2,3, 3,2),
  S10 = obdc(1,3, 2,2),
  S11 = obdc(1,3, 2,2),
  S12 = obdc(1,3, 2,3)
)

build.pTpE <- function(name) {
  mat <- scenarios[[name]]
  I <- nrow_map[[name]]
  J <- ncol_map[[name]]
  pT <- matrix(mat[, 1], nrow = I, ncol = J, byrow = TRUE)
  pE <- matrix(mat[, 2], nrow = I, ncol = J, byrow = TRUE)
  list(pT = pT, pE = pE)
}
