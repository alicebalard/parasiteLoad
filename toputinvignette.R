toyData <- read.csv("https://raw.githubusercontent.com/alicebalard/Article_IntensityEimeriaHMHZ/master/data/cleanedData.csv")
toyData <-toyData[toyData$Aspiculuris_Syphacia != 0,]

m <- mlHyb(Aspiculuris_Syphacia~., data = toyData, model = "negbin", hybridEffect = TRUE)
m

m0 <- mlHyb(Aspiculuris_Syphacia~., data = toyData, model = "negbin", hybridEffect = FALSE)
m0

lrt.mlHyb(m)
##### to do

# groups
mlHyb(Aspiculuris_Syphacia~Sex, data = toyData, model = "negbin", hybridEffect = TRUE)
mlHyb(Aspiculuris_Syphacia~., data = toyData, model = "negbin", hybridEffect = TRUE)


all.vars(Aspiculuris_Syphacia~.)[2]
mytest <- function(formula, data, model, hybridIndex = "HI", myparamBounds = "default",
         hybridEffect = TRUE,
         config = list(optimizer = "optimx", method = c("L-BFGS-B", "bobyqa"), control = list(follow.on = TRUE))){
  # extract response from formula
  response <- all.vars(formula)[1]
  # so far, implemented for 1 categorical group (e.g. sex)
  group <- all.vars(formula)[2]
  group}

mlHyb(Aspiculuris_Syphacia~Sex, data = toyData, model = "negbin", hybridEffect = TRUE)

# extract factor

length(all.vars(Aspiculuris_Syphacia~Sex))
length(all.vars(Aspiculuris_Syphacia~.))
length(all.vars(Aspiculuris_Syphacia~Sex*group))

# fix sides?

# print table 4 hypothesis

# plots
# plot(lm(Aspiculuris_Syphacia~Sex, data = toyData)) Ca ira


# > lm.fit
# function (x, y, offset = NULL, method = "qr", tol = 1e-07, singular.ok = TRUE,
#           ...)
# {
#   if (is.null(n <- nrow(x)))
#     stop("'x' must be a matrix")
#   if (n == 0L)
#     stop("0 (non-NA) cases")
#   p <- ncol(x)
#   if (p == 0L) {
#     return(list(coefficients = numeric(), residuals = y,
#                 fitted.values = 0 * y, rank = 0, df.residual = length(y)))
#   }
#   ny <- NCOL(y)
#   if (is.matrix(y) && ny == 1)
#     y <- drop(y)
#   if (!is.null(offset))
#     y <- y - offset
#   if (NROW(y) != n)
#     stop("incompatible dimensions")
#   if (method != "qr")
#     warning(gettextf("method = '%s' is not supported. Using 'qr'",
#                      method), domain = NA)
#   chkDots(...)
#   z <- .Call(C_Cdqrls, x, y, tol, FALSE)
#   if (!singular.ok && z$rank < p)
#     stop("singular fit encountered")
#   coef <- z$coefficients
#   pivot <- z$pivot
#   r1 <- seq_len(z$rank)
#   dn <- colnames(x)
#   if (is.null(dn))
#     dn <- paste0("x", 1L:p)
#   nmeffects <- c(dn[pivot[r1]], rep.int("", n - z$rank))
#   r2 <- if (z$rank < p)
#     (z$rank + 1L):p
#   else integer()
#   if (is.matrix(y)) {
#     coef[r2, ] <- NA
#     if (z$pivoted)
#       coef[pivot, ] <- coef
#     dimnames(coef) <- list(dn, colnames(y))
#     dimnames(z$effects) <- list(nmeffects, colnames(y))
#   }
#   else {
#     coef[r2] <- NA
#     if (z$pivoted)
#       coef[pivot] <- coef
#     names(coef) <- dn
#     names(z$effects) <- nmeffects
#   }
#   z$coefficients <- coef
#   r1 <- y - z$residuals
#   if (!is.null(offset))
#     r1 <- r1 + offset
#   if (z$pivoted)
#     colnames(z$qr) <- colnames(x)[z$pivot]
#   qr <- z[c("qr", "qraux", "pivot", "tol", "rank")]
#   c(z[c("coefficients", "residuals", "effects", "rank")], list(fitted.values = r1,
#                                                                assign = attr(x, "assign"), qr = structure(qr, class = "qr"),
#                                                                df.residual = n - z$rank))
# }

