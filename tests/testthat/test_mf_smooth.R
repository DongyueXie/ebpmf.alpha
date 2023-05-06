set.seed(12345)
n <- 50
p <- 256
K <- 3
snr <- 1

## Step 1: sample U, an orthogonal matrix
rand_semdef_sym_mat <- crossprod(matrix(runif(n * n), n, n))
rand_ortho_mat <- eigen(rand_semdef_sym_mat)$vector[, 1:K]
u_1 <- rand_ortho_mat[, 1]
u_2 <- rand_ortho_mat[, 2]
u_3 <- rand_ortho_mat[, 3]

f1 = c(rep(0,p/8), rep(1, p/4), rep(0, p/4), rep(0, p/4),rep(0,p/8))
f2 = c(rep(0,p/8), rep(0, p/4), rep(1, p/4), rep(0, p/4),rep(0,p/8))
f3 = c(rep(0,p/8), rep(0, p/4), rep(0, p/4), rep(1, p/4),rep(0,p/8))

L = cbind(u_1,u_2,u_3)
FF=cbind(f1,f2,f3)

M = n / 3 * u_1 %*% t(f1) +
  n / 5 * u_2 %*% t(f2) +
  n / 6 * u_3 %*% t(f3)
v = var(c(M))/snr
X = M + matrix(rnorm(n*p,0,sqrt(v)),nrow=n,ncol=p)

library(flashier)
fl <- flash(X, greedy.Kmax = 10L, backfit = TRUE,
            ebnm.fn = c(ebnm::ebnm_point_normal,
                        ebnm_dwt_haar))

for (k in 1:fl$n.factors) {
  plot(fl$F.pm[,k],type='l')
}



###### ebpmf + dwt prior #########

Y = matrix(rpois(n*p,exp(X)),nrow=n,ncol = p)
fit = ebpmf_log(Y,l0=0,f0=0,
                flash_control = list(fix_f0=T,ebnm.fn=c(ebnm_point_normal,ebnm_dwt_haar)),
                var_type = 'constant')
plot(fit$elbo_trace)
plot(fit$fit_flash$F.pm[,5],type='l')

###### ebpmf + ndwt prior #########

fit = ebpmf_log(Y,l0=0,f0=0,
                flash_control = list(fix_f0=T,ebnm.fn=c(ebnm_point_normal,ebnm_ndwt),add_greedy_warmstart=F,backfit_warmstart=F),
                var_type = 'constant')
plot(fit$elbo_trace)
plot(fit$fit_flash$F.pm[,5],type='l')

###########MoMA example#########
get.X <- function(n=50,p=200,snr=1) {
  #n <- 50
  #p <- 200
  K <- 3
  #snr <- 1

  ## Step 1: sample U, an orthogonal matrix
  rand_semdef_sym_mat <- crossprod(matrix(runif(n * n), n, n))
  rand_ortho_mat <- eigen(rand_semdef_sym_mat)$vector[, 1:K]
  u_1 <- rand_ortho_mat[, 1]
  u_2 <- rand_ortho_mat[, 2]
  u_3 <- rand_ortho_mat[, 3]

  ## Step 2: generate V, the signal
  set_zero_n_scale <- function(x, index_set) {
    x[index_set] <- 0
    x <- x / sqrt(sum(x^2))
    x
  }

  b_1 <- 7 / 20 * p
  b_2 <- 13 / 20 * p

  x <- as.vector(seq(p))

  # Sinusoidal signal
  v_1 <- sin((x + 15) * pi / 17)
  v_1 <- set_zero_n_scale(v_1, b_1:p)

  # Gaussian-modulated sinusoidal signal
  v_2 <- exp(-(x - 100)^2 / 650) * sin((x - 100) * 2 * pi / 21)
  v_2 <- set_zero_n_scale(v_2, c(1:b_1, b_2:p))

  # Sinusoidal signal
  v_3 <- sin((x - 40) * pi / 30)
  v_3 <- set_zero_n_scale(v_3, 1:b_2)

  ## Step 3, the noise
  eps <- matrix(rnorm(n * p), n, p)

  ## Step 4, put the pieces together
  X <- n / 3 * u_1 %*% t(v_1) +
    n / 5 * u_2 %*% t(v_2) +
    n / 6 * u_3 %*% t(v_3) +
    eps

  # Print the noise-to-signal ratio
  cat(paste("norm(X) / norm(noise) = ", norm(X) / norm(eps)))

  # Plot the signals
  yrange <- max(c(v_1, v_2, v_3))
  plot(v_1,
       type = "l",
       ylim = c(-yrange, yrange),
       ylab = "v", xlab = "i",
       main = "Plot of Signals"
  )
  lines(v_2, col = "blue")
  lines(v_3, col = "red")
  legend(0, 0.25,
         legend = expression(v[1], v[2], v[3]),
         lty = 1,
         col = c("black", "blue", "red"),
         cex = 0.6
  )

  return(list(X=X,L = cbind(u_1,u_2,u_3),FF=cbind(v_1,v_2,v_3)))
}

set.seed(12345)
Y = get.X(n=100,p=256)
fl <- flash(Y$X, greedy.Kmax = 10L, backfit = TRUE,
            ebnm.fn = c(ebnm::ebnm_point_normal,
                        ebnm_dwt))

for (k in 1:fl$n.factors) {
  plot(fl$F.pm[,k],type='l')
}

