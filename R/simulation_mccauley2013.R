#' Check McCauley Parameter Vector
#'
#' Checks the pvec contains required parameters.
#'
#' @param pvec The pvec to check contains all required three process model parameters
#'
#' @return logical
mccauley2013_check_pvec <- function(pvec) {
  accepted_names <- c("alpha_w","alpha_s","beta_w","beta_s",
                      "eta_w","Wc","Tp",                # allostatic rates (eta_s derived)
                      "mu_w","mu_s","phi",         # circadian
                      "lambda_w","lambda_s","xi","p0","u0","k0","pf")
  diffvals <- setdiff(names(pvec), accepted_names)
  if (length(diffvals) > 0) {
    msg <- sprintf("mccauley2013_check_pvec model halted!:\n [%s] \n is/are unsupported parameters\nPlease remove these before continuing.",
                   diffvals)
    stop(call. = FALSE, msg)
  }
  if (!all(accepted_names %in% names(pvec))) {
    stop(call. = FALSE,
         "mccauley2013_check_pvec model halted: missing parameters in supplied pvec.\nSee help(ODEM_make_pvec) for required names.")
  }
  TRUE
}

#' Make McCauley Model Parameter Vector
#'
#' Creates a parameter vector with defaults drawn from the published model. The
#' values provided here are placeholders that approximate those reported in the
#' paper and can be altered as required.
#'
#' @param alpha_w Homeostatic dissipation rate for performance during wakefulness
#' @param alpha_s Homeostatic dissipation rate for performance during sleep
#' @param beta_w Scaling factor for \code{u} during wakefulness
#' @param beta_s Scaling factor for \code{u} during sleep
#' @param kappa Build-up rate for \code{u} during wakefulness
#' @param mu_w Offset of circadian process during wakefulness
#' @param mu_s Offset of circadian process during sleep
#' @param phi Phase position of circadian process
#' @param lambda Rate constant for modulation of circadian amplitude
#' @param A Asymptotic amplitude of circadian modulation
#' @param p0 Initial value for the performance process
#' @param u0 Initial value for the homeostatic process
#' @param k0 Initial value for the circadian amplitude
#'
#' @export
mccauley2013_make_pvec <- function(alpha_w = 0.028,   
                               alpha_s = 0.260,   
                               beta_w  = 0.260,   
                               beta_s  = 0.260,   
                               eta_w   = 0.0126,  
                               Wc      =  20.2, 
                               Tp      = 24,
                               mu_w    = 0.466,
                               mu_s    = -1.50,
                               phi     = 21.2,     
                               lambda_w=  0.49,    
                               lambda_s=  -0.49, # note kept this flipped compared to McCauley (2021) and later papers    
                               xi      =  1.09,
                               p0=5.22,u0=38.5,k0=.0212,pf=0) {
  
  pvec <- c(alpha_w=alpha_w, alpha_s=alpha_s, beta_w=beta_w, beta_s=beta_s,
            eta_w=eta_w, Wc=Wc, Tp=Tp,mu_w=mu_w, mu_s=mu_s, phi=phi,
            lambda_w=lambda_w, lambda_s=lambda_s, xi=xi,p0=p0,u0=u0,k0=k0,pf=pf)
  mccauley2013_check_pvec(pvec)
  pvec
}

#' Default parameter vector for the Differential Performance Model
#'
#' @export
mccauley2013_pvec <- mccauley2013_make_pvec()

# Single-harmonic circadian, time t in hours (absolute), phi in hours
# @return Circadian modulation value
mccauley2013_cfun <- function(t, phi=21.2, tau=24) {
  sin( 2*pi * (t - phi)/tau )
}

derive_eta_s <- function(eta_w, Wc, Tp) {
  Ac=Wc/Tp
  if (Ac <= 0 || Ac >= 1) stop("Ac must be in (0,1)")
  eta_w * (Ac / (Ac-1))
}

# Internal helper to append model columns to a FIPS_df
mccauley2013_cols  <- c("p", "u", "k", "c")
mccauley2013_append_model_cols  <- function(.FIPS_df) {
  .FIPS_df[, mccauley2013_cols] <- NA
  .FIPS_df
}

# One RK4 step of the 2013 ODE over step h (hours)
.mccauley2013_derivs <- function(p,u,k, t_abs, wake, pvec) {
  with(as.list(pvec), {
    c_t  <- mccauley2013_cfun(t_abs, phi)
    #c_t  <- unified_Cfun(t_abs%%Tp,phi) # optionally allow use of 5-harmonic Circadian function
    eta_s <- derive_eta_s(eta_w, Wc, Tp)
    g_w <- k * (c_t + mu_w)
    g_s <- k * (c_t + mu_s) + ((alpha_s*beta_s)/eta_w)*((Wc-Tp)/Wc)
    if (wake) {

      dp <- -alpha_w * (p - pf - beta_w * u) + g_w
      du <- eta_w *u
      dk <- lambda_w * k * (1 - (k/xi))
      
    } else {
      
      dp <- -alpha_s * (p - pf - beta_s *u) + g_s
      du <- eta_s * u + 1
      dk <- lambda_s * k
      
    }
    c(dp = unname(dp), du = unname(du), dk = unname(dk), c_t = unname(c_t))
  })
}

.mccauley2013_rk4  <- function(p,u,k,t0, hour_step, wake, pvec) {
  k1 <- .mccauley2013_derivs(p,u,k,t0,wake,pvec)
  k1dp=k1[["dp"]];k1du=k1[["du"]];k1dk=k1[["dk"]]
  k2 <- .mccauley2013_derivs(p + 0.5*hour_step*k1dp, u + 0.5*hour_step*k1du, k + 0.5*hour_step*k1dk, t0 + 0.5*hour_step, wake, pvec)
  k2dp=k2[["dp"]];k2du=k2[["du"]];k2dk=k2[["dk"]]
  k3 <- .mccauley2013_derivs(p + 0.5*hour_step*k2dp, u + 0.5*hour_step*k2du, k + 0.5*hour_step*k2dk, t0 + 0.5*hour_step, wake, pvec)
  k3dp=k3[["dp"]];k3du=k3[["du"]];k3dk=k3[["dk"]]
  k4 <- .mccauley2013_derivs(p +     hour_step*k3dp, u +     hour_step*k3du, k +     hour_step*k3dk, t0 +     hour_step, wake, pvec)
  k4dp=k4[["dp"]];k4du=k4[["du"]];k4dk=k4[["dk"]];k4c_t=k4[["c_t"]]
  p_next <- p + (hour_step/6)*(k1dp + 2*k2dp + 2*k3dp + k4dp)
  u_next <- u + (hour_step/6)*(k1du + 2*k2du + 2*k3du + k4du)
  k_next <- k + (hour_step/6)*(k1dk + 2*k2dk + 2*k3dk + k4dk)
  c_next <- as.numeric(k4c_t)  # last circadian sample (for logging); not critical
  # keep k within [0, xi]
  xi <- pvec[["xi"]]
  if (!is.na(xi)) k_next <- max(0, min(xi, k_next))
  list(p=p_next, u=u_next, k=k_next,c=c_next)
}

#' Simulate: McCauley Model
#'
#' Runs the McCauley et al. (2013) PVT model over the supplied FIPS_df.
#'
#' @param pvec Parameter vector, see [mccauley2013_pvec]
#' @param dat Input dataframe (ensure this is a FIPS_df)
#'
#' @return Dataframe with simulation results
#' @export
#' 
mccauley2013_simulate  <- function(pvec, dat, substep_minutes = 2) {
  # initial state
  p <- pvec[["p0"]]; u <- pvec[["u0"]]; k <- pvec[["k0"]]
  
  # Choose substep
  hour_step <- substep_minutes / 60  # hours
  
  for (i in seq_len(nrow(dat))) {
    if (i>1) {
      dt_row <- time_length(dat$datetime[i] - dat$datetime[i-1], "hour")
      if (is.na(dt_row) || dt_row <= 0) {
        # Fill circadian at row start for completeness
        dat$c[i] <- mccauley2013_cfun(dat$time[i], pvec[["phi"]])
        dat$p[i] <- p; dat$u[i] <- u; dat$k[i] <- k
        next
      }
      wake   <- isTRUE(dat$wake_status[i])
      
      t_abs  <- as.decimaltime(dat$datetime[1])+(dat$sim_hours[i]-dat$sim_hours[1])  # absolute clock in hours
      # Integrate in small steps to avoid large-step error
      remaining <- dt_row
      while (remaining > 0) {
        step <- min(hour_step, remaining)
        res  <- .mccauley2013_rk4(p, u, k, t_abs, step, wake, pvec)
        p <- res$p; u <- res$u; k <- res$k
        t_abs <- t_abs + step
        remaining <- remaining - step
      }
      # Log end-of-row state and circadian (c at row end)
      dat$p[i] <- p
      dat$u[i] <- u
      dat$k[i] <- k
      dat$c[i] <- mccauley2013_cfun(t_abs, pvec[["phi"]])
    } else {
      dat$p[i] <- pvec["p0"]
      dat$u[i] <- pvec["u0"]
      dat$k[i] <- pvec["k0"]
      dat$c[i] <- mccauley2013_cfun(dat$time[i], pvec[["phi"]])
    }
    
  }
  dat
}

#' McCauley Simulation Dispatcher
#'
#' @param dat A FIPS_df
#' @param pvec Parameter vector
#' @param model_formula Optional formula for custom predictions
#' 
mccauley2013_simulation_dispatch  <- function(dat, pvec, model_formula = NULL, substep_minutes = 2) {
  mccauley2013_check_pvec(pvec)
  dat <- mccauley2013_append_model_cols(dat)
  dat <- mccauley2013_simulate(pvec, dat, substep_minutes = substep_minutes)
  dat <- FIPS_simulation(dat, modeltype = "mccauley", pvec = pvec, pred_stat = "fatigue", pred_cols = mccauley2013_cols)
  if (is.null(model_formula)) {
    dat <- process_bmm_formula(dat, fatigue~p, pvec)
  } else {
    dat <- process_bmm_formula(dat, model_formula, pvec)
  }
  dat
}

## Code below computes steady-state parameters from McCauley (2013) Appendix
# State-transition matrices for upper-triangular linear systems
# For A = [[a, b],[0, d]], exp(A*t) = [[e^{a t}, b*(e^{a t}-e^{d t})/(a-d)],[0, e^{d t}]]
.stm_upper <- function(a, b, d, t) {
  ea <- exp(a*t); ed <- exp(d*t)
  off <- if (abs(a - d) < 1e-12) b * t * ea else b * (ea - ed) / (a - d)
  rbind(c(ea, off), c(0, ed))
}
Phi_w <- function(t, alpha_w, beta_w, eta_w) .stm_upper(alpha_w, alpha_w*beta_w, eta_w, t)
Psi_s <- function(t, alpha_s, beta_s, eta_s) .stm_upper(alpha_s, alpha_s*beta_s, eta_s, t)

# Steady-state kappa at wake onset (t0), from piecewise logistic (wake) + exponential (sleep)
# Wake duration W, sleep duration S = T - W
# logistic wake: dκ/dt = λ_w κ (1 - κ/ξ); sleep: dκ/dt = λ_s κ
# steady κ0 solves: κ0 = e^{λ_s S} * ξ / (1 + (ξ/κ0 - 1) e^{-λ_w W})
# => κ0 = ξ * (exp(λ_s*S) - exp(-λ_w*W)) / (1 - exp(-λ_w*W))
kappa0_steady <- function(W, Time, lambda_w, lambda_s, xi) {
  S <- Time - W
  r <- exp(-lambda_w * W)     # wake logistic factor
  s <- exp(lambda_s * S)     # sleep exp factor (since dk/dt = -lambda_s * k)
  xi * (s - r) / (1 - r)
}
# κ(t) over the day in steady state, given κ0 at wake onset t0
kappa_ss_fun <- function(t_abs, t0, W, Time, lambda_w, lambda_s, xi, kappa0) {
  rel <- ((t_abs - t0) %% Time)
  if (rel < W) {
    r <- exp(-lambda_w * rel)
    xi / (1 + (xi / kappa0 - 1) * r)
  } else {
    S <- rel - W
    rW <- exp(-lambda_w * W)
    kw_end <- xi / (1 + (xi / kappa0 - 1) * rW)
    kw_end * exp(lambda_s * S)  # <-- sign
  }
}

# Inhomogeneity vectors b_w(t) and b_s(t):
# Wake: dx/dt = A_w x + [g_w(t); 0],  g_w = κ(t) * (c(t) + mu_w)
# Sleep: dx/dt = A_s x + [g_s(t); 1], g_s = κ(t) * (c(t) + mu_s) + (alpha_s*beta_s)/eta_s
b_w <- function(t, kappa_t, c_t, mu_w) c(kappa_t * (c_t + mu_w), 0)
b_s <- function(t, kappa_t, c_t, mu_s, alpha_s, beta_s, eta_s) {
  c(kappa_t * (c_t + mu_s) + (alpha_s*beta_s)/eta_s, 1)
}

.convolve_forcing <- function(L, t0, dt_target,
                              stm_fun, b_fun,
                              circ_fun, kappa_fun) {
  n  <- max(1L, ceiling(L / dt_target))
  h  <- L / n
  tau <- seq(0, L, by = h)
  m  <- length(tau)
  
  vals <- matrix(0, nrow = 2, ncol = m)
  for (j in seq_len(m)) {
    t   <- t0 + tau[j]
    STM <- stm_fun(L - tau[j])
    b   <- b_fun(t, kappa_fun(t), circ_fun(t))
    vals[, j] <- STM %*% b
  }
  
  middle <- if (m > 2) rowSums(vals[, 2:(m-1), drop = FALSE]) else c(0, 0)
  area   <- h * (0.5 * vals[, 1] + middle + 0.5 * vals[, m])
  as.numeric(area)  # length 2
}

# Returns p0,u0,k0 at wake onset, plus state at sleep onset (end of wake) for convenience
#' @export
#' 
mccauley_initial_values <- function(pvec, W, Time = 24, t0 = 0, dt_minutes = 1) {
  # unpack
  alpha_w <- -pvec[["alpha_w"]]; beta_w <- -pvec[["beta_w"]]
  alpha_s <- -pvec[["alpha_s"]]; beta_s <- -pvec[["beta_s"]]
  eta_w   <- pvec[["eta_w"]];   Wc     <- pvec[["Wc"]]
  Tp      <- pvec[["Tp"]];
  mu_w    <- pvec[["mu_w"]];    mu_s   <- pvec[["mu_s"]]
  phi     <- pvec[["phi"]]
  lambda_w<- pvec[["lambda_w"]]; lambda_s <- pvec[["lambda_s"]]
  xi      <- pvec[["xi"]]

  eta_s <- derive_eta_s(eta_w, Wc, Tp)
  S <- Time - W
  dt <- dt_minutes / 60  # hours
  
  k0 <- kappa0_steady(W, Time, lambda_w, lambda_s, xi)
  rW <- exp(-lambda_w * W)
  k_sleep_start <- xi / (1 + (xi / k0 - 1) * rW)
  k_fun <- function(t) kappa_ss_fun(t, t0, W, Time, lambda_w, lambda_s, xi, k0)
  c_fun <- function(t) mccauley2013_cfun(t, phi)  # Eq. (6)
  
  Phi <- function(tau) Phi_w(tau, alpha_w, beta_w, eta_w)
  Psi <- function(tau) Psi_s(tau, alpha_s, beta_s, eta_s)
  
  wake_forcing <- function(t, kappa_t, c_t) b_w(t, kappa_t, c_t, mu_w)
  sleep_forcing<- function(t, kappa_t, c_t) b_s(t, kappa_t, c_t, mu_s, alpha_s, beta_s, eta_s)
  
  I_w <- .convolve_forcing(W, t0, dt, Phi, wake_forcing, c_fun, k_fun)
  I_s <- .convolve_forcing(S, t0 + W, dt, Psi, sleep_forcing, c_fun, k_fun)
  
  M <- Psi(S) %*% Phi(W)
  I2 <- diag(2)
  rhs <- I_w + I_s
  x0 <- solve(I2 - M, rhs)
  p0 <- unname(x0[1]); u0 <- unname(x0[2])
  
  # Also report state at sleep onset (end of wake): xW = Phi(W) x0 + I_w
  xW <- Phi(W) %*% x0 + I_w
  list(p0 = p0, u0 = u0, k0 = k0,h0=0,
       p_sleep_start = unname(xW[1]), u_sleep_start = unname(xW[2]),k_sleep_start=k_sleep_start,
       eta_s = eta_s)
}
