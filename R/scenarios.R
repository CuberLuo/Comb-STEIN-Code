## Diversified simulation scenarios derived from the manuscript setting.
## All toxicity and efficacy probabilities are stored to at most two decimal
## places. Toxicity is nondecreasing in both agents. Under the utility rule,
## S1-S4 have one OBDC, S5-S8 have two OBDCs, and S9-S12 have three OBDCs;
## the complete OBDC location set is different in every scenario.

## Convert row-wise toxicity/efficacy value pairs into a two-column matrix.
mk <- function(...) {
  matrix(c(...), ncol = 2, byrow = TRUE)
}

scenarios <- list(
  S1 = mk(
    0.05,0.35, 0.10,0.80, 0.15,0.65, 0.20,0.50, 0.30,0.40,
    0.10,0.40, 0.15,0.60, 0.20,0.55, 0.30,0.45, 0.40,0.35,
    0.15,0.45, 0.20,0.50, 0.30,0.40, 0.40,0.35, 0.50,0.30
  ),
  S2 = mk(
    0.05,0.25, 0.10,0.35, 0.15,0.50, 0.20,0.65, 0.30,0.70,
    0.10,0.30, 0.15,0.45, 0.20,0.60, 0.25,0.85, 0.35,0.65,
    0.15,0.35, 0.20,0.50, 0.25,0.55, 0.35,0.60, 0.45,0.50
  ),
  S3 = mk(
    0.05,0.20, 0.10,0.30, 0.15,0.40, 0.25,0.50, 0.35,0.45,
    0.10,0.30, 0.15,0.45, 0.20,0.60, 0.30,0.65, 0.40,0.55,
    0.15,0.40, 0.25,0.85, 0.30,0.70, 0.40,0.55, 0.50,0.40
  ),
  S4 = mk(
    0.05,0.20, 0.10,0.35, 0.15,0.50, 0.20,0.70, 0.30,0.90,
    0.10,0.25, 0.15,0.40, 0.20,0.55, 0.30,0.65, 0.40,0.75,
    0.15,0.30, 0.20,0.45, 0.30,0.60, 0.40,0.70, 0.50,0.80
  ),
  S5 = mk(
    0.05,0.25, 0.10,0.45, 0.20,0.80, 0.30,0.65, 0.40,0.50,
    0.10,0.40, 0.20,0.80, 0.30,0.70, 0.40,0.55, 0.50,0.40,
    0.20,0.55, 0.30,0.65, 0.40,0.60, 0.50,0.45, 0.60,0.30
  ),
  S6 = mk(
    0.05,0.25, 0.10,0.45, 0.20,0.65, 0.30,0.85,
    0.10,0.40, 0.20,0.60, 0.30,0.85, 0.40,0.70,
    0.20,0.50, 0.30,0.65, 0.40,0.70, 0.50,0.55
  ),
  S7 = mk(
    0.05,0.25, 0.10,0.45, 0.20,0.60, 0.30,0.55,
    0.10,0.40, 0.20,0.80, 0.30,0.70, 0.40,0.60,
    0.20,0.80, 0.30,0.65, 0.40,0.55, 0.50,0.45
  ),
  S8 = mk(
    0.05,0.20, 0.10,0.35, 0.15,0.50, 0.25,0.60,
    0.10,0.30, 0.15,0.50, 0.25,0.65, 0.35,0.85,
    0.15,0.40, 0.25,0.60, 0.35,0.85, 0.45,0.70
  ),
  S9 = mk(
    0.05,0.30, 0.20,0.80, 0.30,0.60,
    0.20,0.80, 0.20,0.80, 0.35,0.65,
    0.30,0.50, 0.35,0.60, 0.45,0.55
  ),
  S10 = mk(
    0.05,0.30, 0.10,0.55, 0.25,0.85,
    0.10,0.55, 0.25,0.85, 0.35,0.70,
    0.25,0.85, 0.35,0.70, 0.45,0.60
  ),
  S11 = mk(
    0.10,0.75, 0.10,0.75, 0.25,0.60,
    0.10,0.75, 0.25,0.60, 0.40,0.50
  ),
  S12 = mk(
    0.05,0.25, 0.15,0.50, 0.30,0.85,
    0.15,0.45, 0.30,0.85, 0.30,0.85
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

## Convert row/column index pairs into a named two-column OBDC matrix.
obdc <- function(...) {
  matrix(
    c(...),
    ncol = 2,
    byrow = TRUE,
    dimnames = list(NULL, c("i", "j"))
  )
}

table1.obdcs <- list(
  S1  = obdc(1,2),
  S2  = obdc(2,4),
  S3  = obdc(3,2),
  S4  = obdc(1,5),
  S5  = obdc(1,3, 2,2),
  S6  = obdc(1,4, 2,3),
  S7  = obdc(2,2, 3,1),
  S8  = obdc(2,4, 3,3),
  S9  = obdc(1,2, 2,1, 2,2),
  S10 = obdc(1,3, 2,2, 3,1),
  S11 = obdc(1,1, 1,2, 2,1),
  S12 = obdc(1,3, 2,2, 2,3)
)

## Reshape a named scenario into toxicity and efficacy dose-combination grids.
## Returns a list containing matrices `pT` and `pE`.
build.pTpE <- function(name) {
  mat <- scenarios[[name]]
  I <- nrow_map[[name]]
  J <- ncol_map[[name]]
  pT <- matrix(mat[, 1], nrow = I, ncol = J, byrow = TRUE)
  pE <- matrix(mat[, 2], nrow = I, ncol = J, byrow = TRUE)
  list(pT = pT, pE = pE)
}
