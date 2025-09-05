#' Check McCauley Parameter Vector
#'
#' Checks the pvec contains required parameters.
#'
#' @param pvec The pvec to check contains all required model parameters
#'
#' @return logical
mccauley2024_check_pvec <- function(pvec) {
  accepted_names <- c("alpha_w","alpha_s","beta_w","beta_s","eta_s",
                      "eta_w","Wc","Tp",               
                      "mu_w","mu_s","phi",         
                      "lambda_w","lambda_s","xi_u","xi_k","xi_h","p0","u0","k0","zeta_w","zeta_s","pf","v_w","v_s","gamma_h","h0")
  diffvals <- setdiff(names(pvec), accepted_names)
  if (length(diffvals) > 0) {
    msg <- sprintf("mccauley2024_check_pvec model halted!:\n [%s] \n is/are unsupported parameters\nPlease remove these before continuing.",
                   diffvals)
    stop(call. = FALSE, msg)
  }
  if (!all(accepted_names %in% names(pvec))) {
    stop(call. = FALSE,
         "mccauley2024_check_pvec model halted: missing parameters in supplied pvec.\nSee help(ODEM_make_pvec) for required names.")
  }
  TRUE
}

# Single-harmonic circadian, time t in hours (absolute), phi in hours
# @return Circadian modulation value
mccauley2024_cfun <- function(t, phi=21.2, tau=24) {
  sin( 2*pi * (t - phi)/tau )
}

derive_eta_s_2024 <- function(eta_w, Wc, Tp) {
  eta_w * (Wc / (Tp-Wc))
}


#' Make McCauley (2024) Model Parameter Vector
#'
#' Creates a parameter vector with defaults drawn from the published model.
#'
#' @param alpha_w Homeostatic build-up rate for performance during wakefulness
#' @param alpha_s Homeostatic dissipation rate for performance during sleep
#' @param beta_w Scaling factor for impact of u(t) onto p(t) during wakefulness
#' @param beta_s Scaling factor for impact of u(t) onto p(t) during sleep
#' @param eta_w Build-up rate for u(t) during wakefulness
#' @param eta_s Dissipation rate for u(t) during sleep
#' @param Wc Critical threshold (for bifurcation)
#' @param Tp Time period, usually 24hr
#' @param mu_w Offset of circadian process during wakefulness
#' @param mu_s Offset of circadian process during sleep
#' @param phi Phase position of circadian process
#' @param lambda Rate constant for modulation of circadian amplitude
#' @param xi Metric-dependent scaling factor for the effects of u(t) k(t) and h(t)
#' @param zeta Scaling factor for impact of p(t) onto u(t)
#' @param gamma_h Scaling factor for dynamic asymptote of h(t) during sleep
#' @param v Decay rate for h(t) during wakefulness/sleep
#' @param p0 Initial value for the performance process
#' @param u0 Initial value for the homeostatic process
#' @param k0 Initial value for the circadian amplitude
#'
#' @export
mccauley2024_make_pvec <- function(alpha_w = 0.028,   
                               alpha_s = 0.260,   
                               beta_w  = 0.260,   
                               beta_s  = 0.260,   
                               eta_w   = 0.0126,  
                               Wc      =  20.2, 
                               Tp      = 24,
                               mu_w    = 0.46,
                               mu_s    = -1.50,
                               phi     = 21.2,     
                               lambda_w=  0.49,    
                               lambda_s=  0.49,
                               xi_u      =  1.09,
                               xi_k      =  1.09,
                               xi_h      =  1.09,
                               p0=5.22,u0=38.5,k0=.0212,pf=0,h0=0.02, # placeholders, use inits function
                               # NEW PARAMETERS from 2024 papers
                               zeta_w  = 1.31,    
                               zeta_s  = 1.31,
                               v_w=1.37,
                               v_s=1.37,
                               gamma_h=0.71) {
  eta_s <- derive_eta_s_2024(eta_w, Wc, Tp)
  pvec <- c(alpha_w=alpha_w, alpha_s=alpha_s, beta_w=beta_w, beta_s=beta_s,
            eta_w=eta_w, eta_s=eta_s,Wc=Wc, Tp=Tp,mu_w=mu_w, mu_s=mu_s, phi=phi,
            lambda_w=lambda_w, lambda_s=lambda_s, xi_u=xi_u,xi_k=xi_k,xi_h=xi_h,p0=p0,u0=u0,k0=k0,h0=h0,pf=pf,
            zeta_w=zeta_w,zeta_s=zeta_s,gamma_h=gamma_h,v_w=v_w,v_s=v_s)
  mccauley2024_check_pvec(pvec)
  pvec
}

#' Default parameter vector for the McCauley 2024 model
#'
#' @export
mccauley2024_pvec <- mccauley2024_make_pvec()


# Internal helper to append model columns to a FIPS_df
mccauley2024_cols  <- c("p", "u", "k", "h","c")
mccauley2024_append_model_cols  <- function(.FIPS_df) {
  .FIPS_df[, mccauley2024_cols] <- NA
  .FIPS_df
}

# One RK4 step of the 2024 ODE over step h (hours)
.mccauley2024_derivs <- function(p,u,k,h, t_abs, wake, pvec) {
  with(as.list(pvec), {
    if (wake) {
      #c_t  <- unified_Cfun(t_abs%%Tp,Tp-phi) # optionally allow use of 5-harmonic Circadian function
      c_t  <- mccauley2024_cfun(t_abs, phi)
      g_w <- k * (c_t + mu_w)
      dp <- -alpha_w * (p - pf - beta_w*u) + xi_h*(alpha_w -v_w)*h + xi_k*g_w
      du <- -eta_w * (zeta_w * xi_u * (p - pf) - u) + eta_w*zeta_w*xi_u*xi_h*h
      dk <- lambda_w * k * (1 - k)
      dh <- -v_w*h
    } else {
      #c_t  <- unified_Cfun(t_abs%%Tp,(Tp-phi)+12)
      c_t  <- mccauley2024_cfun(t_abs, phi)
      g_s <- k * (c_t + mu_s)
      dp <- -alpha_s * ( (1-  ((xi_h*gamma_h*v_s)/alpha_s))*(p - pf) - beta_s*(u - (1/eta_s))) + xi_h*(alpha_s - v_s*(1+gamma_h*xi_h))*h + xi_k*g_s
      du <- eta_s*(zeta_s*xi_u *(p-pf) - u) - eta_s*zeta_s*xi_u*xi_h*h + 1
      dk <- -lambda_s * k
      dh <- -v_s*( (1+gamma_h*xi_h)*h - gamma_h*(p-pf))
    }
    c(dp = unname(dp), du = unname(du), dk = unname(dk), dh = unname(dh), c_t = unname(c_t))
  })
}

.mccauley2024_rk4  <- function(p,u,k,h,t0, hour_step, wake, pvec) {
  k1 <- .mccauley2024_derivs(p,u,k,h,t0,wake,pvec)
  k1dp=k1[["dp"]];k1du=k1[["du"]];k1dk=k1[["dk"]];k1dh=k1[["dh"]]
  k2 <- .mccauley2024_derivs(p + 0.5*hour_step*k1dp, u + 0.5*hour_step*k1du, k + 0.5*hour_step*k1dk, h + 0.5*hour_step*k1dh, t0 + 0.5*hour_step, wake, pvec)
  k2dp=k2[["dp"]];k2du=k2[["du"]];k2dk=k2[["dk"]];k2dh=k2[["dh"]]
  k3 <- .mccauley2024_derivs(p + 0.5*hour_step*k2dp, u + 0.5*hour_step*k2du, k + 0.5*hour_step*k2dk, h + 0.5*hour_step*k2dh, t0 + 0.5*hour_step, wake, pvec)
  k3dp=k3[["dp"]];k3du=k3[["du"]];k3dk=k3[["dk"]];k3dh=k3[["dh"]]
  k4 <- .mccauley2024_derivs(p +     hour_step*k3dp, u +     hour_step*k3du, k +     hour_step*k3dk, h +     hour_step*k3dh, t0 +     hour_step, wake, pvec)
  k4dp=k4[["dp"]];k4du=k4[["du"]];k4dk=k4[["dk"]];k4dh=k4[["dh"]];k4c_t=k4[["c_t"]]
  p_next <- p + (hour_step/6)*(k1dp + 2*k2dp + 2*k3dp + k4dp)
  u_next <- u + (hour_step/6)*(k1du + 2*k2du + 2*k3du + k4du)
  k_next <- k + (hour_step/6)*(k1dk + 2*k2dk + 2*k3dk + k4dk)
  h_next <- h + (hour_step/6)*(k1dh + 2*k2dh + 2*k3dh + k4dh)
  c_next <- as.numeric(k4[["c_t"]])  # last circadian sample (for logging); not critical
  # keep k within [0, xi]
  xi <- pvec[["xi_k"]]
  if (!is.na(xi)) k_next <- max(0, min(xi, k_next))
  list(p=p_next, u=u_next, k=k_next,c=c_next, h=h_next)
}

#' Simulate: McCauley Model
#' Runs the McCauley et al. (2024) model over the supplied FIPS_df.
#'
#' @param pvec Parameter vector, see [mccauley2024_pvec]
#' @param dat Input dataframe (ensure this is a FIPS_df)
#'
#' @return Dataframe with simulation results
#' @export
#' 
mccauley2024_simulate  <- function(pvec, dat, substep_minutes = 2) {
  # initial state
  p <- pvec[["p0"]]; u <- pvec[["u0"]]; k <- pvec[["k0"]]; h <- pvec[["h0"]]
  
  # Choose substep
  hour_step <- substep_minutes / 60  # hours
  
  for (i in seq_len(nrow(dat))) {
    if (i>1) {
      dt_row <- time_length(dat$datetime[i] - dat$datetime[i-1], "hour")
      if (is.na(dt_row) || dt_row <= 0) {
        # Fill circadian at row start for completeness
        dat$c[i] <- mccauley2024_cfun(dat$time[i], pvec[["phi"]])
        dat$p[i] <- p; dat$u[i] <- u; dat$k[i] <- k; dat$h[i] <- h
        next
      }
      wake   <- dat$wake_status[i]
      
      t_abs  <- dat$time[i]
      # Integrate in small steps to avoid large-step error
      remaining <- dt_row
      while (remaining > 0) {
        step <- min(hour_step, remaining)
        res  <- .mccauley2024_rk4(p, u, k, h, t_abs, step, wake, pvec)
        p <- res$p; u <- res$u; k <- res$k; h <- res$h; c <-res$c
        t_abs <- t_abs + step
        remaining <- remaining - step
      }
      # Log end-of-row state and circadian (c at row end)
      dat$p[i] <- p
      dat$u[i] <- u
      dat$k[i] <- k
      dat$h[i] <- h
      dat$c[i] <- c
     
    } else {
      dat$p[i] <- pvec["p0"]
      dat$u[i] <- pvec["u0"]
      dat$k[i] <- pvec["k0"]
      dat$h[i] <- pvec["h0"]
      dat$c[i] <- mccauley2024_cfun(dat$time[i], pvec[["phi"]])
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
mccauley2024_simulation_dispatch  <- function(dat, pvec, model_formula = NULL, substep_minutes = 2) {
  mccauley2024_check_pvec(pvec)
  dat <- mccauley2024_append_model_cols(dat)
  dat <- mccauley2024_simulate(pvec, dat, substep_minutes = substep_minutes)
  dat <- FIPS_simulation(dat, modeltype = "mccauley", pvec = pvec, pred_stat = "fatigue", pred_cols = mccauley2024_cols)
  if (is.null(model_formula)) {
    dat <- process_bmm_formula(dat, fatigue~p, pvec)
  } else {
    dat <- process_bmm_formula(dat, model_formula, pvec)
  }
  dat
}

## Code to create an equilibrium state for initialisation
#' Compute steady-state initial conditions for the McCauley (2024) model
#'
#' Simulates a stable schedule for `ndays` to approximate equilibrium, then
#' returns p0/u0/k0/h0 at the requested simulation start time-of-day.
#'
#' @param sim_start POSIXct. Timestamp of the *real* simulation start; only the
#'   time-of-day is used to pick initial states on the final equilibrium day.
#' @param sleep_hrs Numeric hours of nightly sleep (default 8).
#' @param wake_time Character "HH:MM:SS" (local clock) for daily wake time (default "07:30:00").
#' @param ndays Integer number of adaptation days (default 30).
#' @param tz Timezone (default taken from `sim_start`, else "UTC").
#' @param round_minutes Rounding used by `parse_sleeptimes` (default 30). Must divide 60.
#' @param pvec Parameter vector used for the adaptation run. `p0/u0/k0/h0` in
#'   `pvec` are ignored; the returned inits will replace them. Defaults to `mccauley2024_pvec`.
#' @param substep_minutes Internal ODE substep (default 2) passed to the simulator.
#' @param return_equilibrium If TRUE, also return the final-day dataframe.
#' @param save_equilibrium_path Optional path to save the final-day equilibrium (RDS).
#'
#' @return Named vector c(p0=..., u0=..., k0=..., h0=...) or a list with
#'   $inits and $equilibrium if `return_equilibrium=TRUE`.
#' @export
mccauley2024_make_inits <- function(
    sim_start,
    sleep_hrs = 8,
    wake_time = "07:30:00",
    ndays = 30,
    tz = NULL,
    round_minutes = 30,
    pvec = mccauley2024_pvec,
    substep_minutes = 2,
    return_equilibrium = FALSE,
    return_adaptation = FALSE,
    save_equilibrium_path = NULL
) {
  stopifnot(inherits(sim_start, "POSIXct"))
  if (is.null(tz)) tz <- attr(sim_start, "tzone") %||% "UTC"
  if (!mccauley2024_check_pvec(pvec)) stop("Invalid `pvec`.")
  
  # Anchor the adaptation period `ndays` before sim_start, at the daily wake time.
  sim_start_date <- lubridate::as_date(sim_start, tz = tz)
  anchor_wake <- lubridate::ymd_hms(
    paste0(sim_start_date - lubridate::days(ndays), " ", wake_time),
    tz = tz
  )
  
  cycles <- tibble::tibble(
    n_days    = ndays,
    sleep_hrs = sleep_hrs,
    wake_time = wake_time,
    start_idx = 0L
  )
  
  # Build a regular schedule and parse into a FIPS_df
  sleeptimes <- build_sleeptimes(anchor_wake = anchor_wake, cycles = cycles, tz = tz)
  
  series_start <- lubridate::ymd_hms(sleeptimes$sleep.start[1], tz = tz)
  series_end   <- lubridate::ymd_hms(sleeptimes$sleep.end[nrow(sleeptimes)], tz = tz)
  series_end   <- lubridate::floor_date(series_end, "day") + lubridate::hours(23) + lubridate::minutes(59)
  
  adaptation_df <- parse_sleeptimes(
    sleeptimes      = sleeptimes,
    series.start    = series_start,
    series.end      = series_end,
    sleep.start.col = "sleep.start",
    sleep.end.col   = "sleep.end",
    sleep.id.col    = "sleep.id",
    roundvalue      = round_minutes
  )
  
  # Run the model across the adaptation period
  adaptation_run <- FIPS_simulate(
    FIPS_df   = adaptation_df,
    modeltype = "mccauley2024",
    pvec      = pvec
  )
  
  # Final day of adaptation
  last_day_start <- lubridate::floor_date(adaptation_run$datetime[nrow(adaptation_run)], "day")
  eq <- dplyr::filter(adaptation_run, datetime >= last_day_start)
  
  # Target time-of-day = time-of-day of sim_start
  target_tod <- as.decimaltime(sim_start)
  # Find closest row (tolerant to rounding)
  idx <- which.min(abs(eq$time - target_tod))
  if (length(idx) == 0L || is.na(idx)) stop("Failed to locate target time-of-day in equilibrium day.")
  
  inits <- c(
    p0 = unname(eq$p[idx]),
    u0 = unname(eq$u[idx]),
    k0 = unname(eq$k[idx]),
    h0 = unname(eq$h[idx])
  )
  
  if (!is.null(save_equilibrium_path)) {
    saveRDS(eq, file = save_equilibrium_path)
  }
  
  if (return_equilibrium & !return_adaptation) {
    return(list(inits = inits, equilibrium = eq))
  }
  if (!return_equilibrium & return_adaptation) {
    return(list(inits = inits, adaptation = adaptation_run))
  }
  if (return_equilibrium & return_adaptation) {
    return(list(inits = inits, equilibrium = eq, adaptation = adaptation_run))
  }
  inits
}

#' Convenience: inject steady-state inits into a parameter vector
#'
#' @inheritParams mccauley2024_make_inits
#' @return Updated `pvec` with p0/u0/k0/h0 set to steady-state values.
#' @export
mccauley2024_init_pvec <- function(
    sim_start,
    sleep_hrs = 8,
    wake_time = "07:30:00",
    ndays = 30,
    tz = NULL,
    round_minutes = 30,
    pvec = mccauley2024_pvec,
    substep_minutes = 2
) {
  inits <- mccauley2024_make_inits(
    sim_start = sim_start,
    sleep_hrs = sleep_hrs,
    wake_time = wake_time,
    ndays = ndays,
    tz = tz,
    round_minutes = round_minutes,
    pvec = pvec,
    substep_minutes = substep_minutes,
    return_equilibrium = FALSE
  )
  pvec[c("p0", "u0", "k0", "h0")] <- inits
  pvec
}


