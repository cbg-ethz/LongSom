beta_binom_params <- function (x, y){
			suppressPackageStartupMessages(library(VGAM))
			fit = invisible(vglm(cbind(x,y) ~ 1, betabinomialff, trace = FALSE))
			alpha <- Coef(fit)[[1]]
			beta <- Coef(fit)[[2]]
			return(c(alpha,beta))
		}