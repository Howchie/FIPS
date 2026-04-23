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
mccauley_cfun <- function(t, phi=21.2, tau=24) {
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
                               p0=5.22,u0=28.5,k0=.0212,pf=0,h0=0.02, # placeholders, use inits function
                               # NEW PARAMETERS from 2024 papers
                               zeta_w  = 1.31,    
                               zeta_s  = 1.31,
                               v_w=1.37,
                               v_s=1.37,
                               gamma_h=0.71
                               ) {
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
mccauley2024_cols <- c("p","u","k","h")
# Root function: return 0 at switching moments
# check transition to/from bed
.mccauley_rootfun <- function(t, y, parms, pvec_C, df, change_points, dynamic_vars,...) {
  if (t > 1) {
    g_change <- as.numeric(all((t - change_points) != 0))
    g_CBTmin <- 1
    if ("x" %in% names(y)) { # if we're fitting the dynamic state variable, use that
      if (dynamic_vars$model == "K99") {
        cur_phase = atan2(y["xc"], y["x"])
        g_CBTmin = -2.98 - cur_phase # cross the CBT-min threshold (Ref Woelders et al 2017, citing May et al 2002)
      } else if (dynamic_vars$model == "FJK") {
        g_CBTmin = unlist(.FJK_derivs(t, y, pvec_C, df, dynamic_vars = dynamic_vars))[2]
      }
    }
  } else {g_change=1;g_CBTmin=1}
  return(c(g_change,g_CBTmin))
}

.mccauley_eventfun <- function(t, y, pvec, pvec_C, df, change_points, dynamic_vars,...) {
  g = .mccauley_rootfun(t, y, pvec, pvec_C, df, change_points, dynamic_vars)
  if (g[1] == 0) {
    j = which(df$t_abs == t)
    if (df$switch_direction[j] == "Wake") {
      y["S"] = 0
    } else {
      y["S"] = 1
    }
  } else if (abs(g[2])<1e-4) {
    if (dynamic_vars$model == "K99") {
      cur_phase = atan2(y["xc"], y["x"])
      if (cur_phase < 0) { # needed to ensure we don't pick up the flip to the positive angle
        y["CBT_min"] = (t %% 24) + pvec_C["phi_adj"]
      }
    } else if (dynamic_vars$model == "FJK") {
      if (y["x"] < 0) {
        y["CBT_min"] = (t %% 24)
      }
    }
  }
  return(y)
}

# One RK4 step of the 2024 ODE over step h (hours)
.mccauley2024_derivs <- function(t, state, parms, df, ...) {
  with(as.list(c(state, parms)), {
    j <- findInterval(t, df$t_abs, rightmost.closed = TRUE) # use the row corresponding to t<=t_abs
    if (!("x" %in% names(state))) {
      C <- mccauley_cfun(t, phi)
    } else {
      C=-x # McCauley's model uses a peaking at night function
    }
    if (S==0) {
      g_w <- k * (C + mu_w)
      dp <- -alpha_w * (p - pf - beta_w * u) + xi_h * (alpha_w - v_w) * h + xi_k * g_w
      du <- -eta_w * (zeta_w * xi_u * (p - pf) - u) + eta_w * zeta_w * xi_u * xi_h * h
      dk <- lambda_w * k * (1 - k)
      dh <- -v_w * h
    } else {
      g_s <- k * (-C + mu_s)
      dp <- -alpha_s * ((1 - ((xi_h * gamma_h * v_s) / alpha_s)) * (p - pf) - beta_s * (u - (1 / eta_s))) + xi_h * (alpha_s - v_s * (1 + gamma_h * xi_h)) * h + xi_k * g_s
      du <- eta_s * (zeta_s * xi_u * (p - pf) - u) - eta_s * zeta_s * xi_u * xi_h * h + 1
      dk <- -lambda_s * k
      dh <- -v_s * ((1 + gamma_h * xi_h) * h - gamma_h * (p - pf))
    }
    c(unname(dp), unname(du), unname(dk), unname(dh),0)
  })
}

## Coupling for McCauley (2024) model with dynamic circadian (FJK - Forger et al 1999)
## This keeps each subsystem modular while integrating as one ODE when needed.
# Coupled function that shares S and C across subsystems if desired
.mccauley2024_derivs_coupled <- function(t, y, parms, df, pvec_C,use_dynamic,dynamic_vars,...) {
   # Phillips derivatives depend on dynamic C
   d_mccauley <- .mccauley2024_derivs(
      t = t,
      state = y,
      parms = parms,
      df = df
  )

  if (use_dynamic) {
    if (dynamic_vars$model == "FJK") {
      d_c <- unlist(.FJK_derivs(
      t = t,
      state = y,
      parms = pvec_C,
      df = df,
      dynamic_vars = dynamic_vars
      ))
    } else if (dynamic_vars$model=="K99") {
        d_c <- unlist(.K99_derivs(
          t = t,
          state = y,
          parms = pvec_C,
          df = df,
          dynamic_vars = dynamic_vars
        ))
    }
  } else {
     return(list(d_mccauley))
  }
  return(list(c(d_mccauley, d_c)))
}

mccauley2024_simulation_dispatch  <- function(
  dat, 
  pvec=mccauley2024_pvec,  pvec_C=FJK_pvec,
  method = "lsoda",
  use_dynamic = TRUE,
  dynamic_vars = fips_default_dynamic_vars(model = "FJK"),
  model_formula=NULL) {
  dynamic_vars <- fips_normalize_dynamic_vars(dynamic_vars)
  # Ensure required inputs exist
  mccauley2024_check_pvec(pvec)
  if (dynamic_vars$model == "FJK") {
    FJK_check_pvec(pvec_C)
  } else if (dynamic_vars$model=="K99") {
    K99_check_pvec(pvec_C)
  }
  pvec_C <- .fips_adjust_pvec_for_light_metric(pvec_C, dynamic_vars, model_name = dynamic_vars$model)


  # Prepare time and schedule columns
  dat$t_abs <- dat$sim_hours + dat$time[1]
  change_points <- dat$t_abs[dat$change_point == 1]

  dat <- .fips_build_light_drive(dat, dynamic_vars)
  S0 = as.numeric(dat$wake_status[1] == 0)
  y0 <- c(
    p = pvec[["p0"]],
    u = pvec[["u0"]],
    k = pvec[["k0"]],
    h = pvec[["h0"]],
    S = S0
  )
  if (use_dynamic) {
    y0 = c(y0,
      xc = pvec_C[["xc0"]],
      x = pvec_C[["x0"]]
    ) 
    if (!dynamic_vars$approx_N) {
      y0 = c(y0, L = pvec_C[["L0"]])
    }
    y0 = c(y0, CBT_min = pvec_C[["CBT_min0"]])
  }

  # Integrate the coupled system
  ODE_results <- as.data.frame(deSolve::ode(
    y = y0,
    times = dat$t_abs,
    func = .mccauley2024_derivs_coupled,
    rootfun = .mccauley_rootfun,
    parms = pvec,
    method = method,
    events = list(func = .mccauley_eventfun, root = TRUE, maxroot = 1e6),
    df = dat,
    change_points = change_points,
    pvec_C = pvec_C,
    use_dynamic = use_dynamic,
    dynamic_vars=dynamic_vars
  )) %>% dplyr::rename(t_abs = time)

  # Merge back to data; compute Phillips outputs using dynamic C from the path
  dat <- dat %>%
    dplyr::left_join(ODE_results, by = "t_abs")
  
  # Keep downstream processing consistent with Phillips interface
  if (use_dynamic & dynamic_vars$model=="FJK") {
    pred_cols = unique(c(mccauley2024_cols, FJK_cols))
  } else if (use_dynamic & dynamic_vars$model=="K99") {
    pred_cols = unique(c(mccauley2024_cols, K99_cols))
  } else {
    pred_cols = mccauley2024_cols
    dat$C = mccauley_cfun(dat$t_abs,pvec["phi"])
  }
  dat <- FIPS_simulation(dat,
    modeltype = "mccauley2024", pvec = pvec, pred_stat = "fatigue",
    pred_cols = pred_cols
  )
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
    series_start_time = 23.5,
    ndays = 30,
    tz = NULL,
    round_minutes = 5,
    pvec = mccauley2024_pvec,
    pvec_C = FJK_pvec,
    dynamic_C = FALSE,
    dynamic_vars=fips_default_dynamic_vars(model = "FJK"),
    return_equilibrium = FALSE,
    return_adaptation = FALSE,
    save_equilibrium_path = NULL
) {
  dynamic_vars <- fips_normalize_dynamic_vars(dynamic_vars)
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
    pvec      = pvec,
    pvec_C = pvec_C,
    dynamic_C = dynamic_C,
    dynamic_vars=dynamic_vars
  )
  
  # Final day of adaptation
  last_day_start <- lubridate::floor_date(adaptation_run$datetime[nrow(adaptation_run)], "day")
  eq <- dplyr::filter(adaptation_run, datetime >= last_day_start)
  
  # Target time-of-day = time-of-day of sim_start
  target_tod <- series_start_time
  # Find closest row (tolerant to rounding)
  idx <- which.min(abs(eq$time - target_tod))
  if (length(idx) == 0L || is.na(idx)) stop("Failed to locate target time-of-day in equilibrium day.")
  
  if (dynamic_C & !dynamic_vars$approx_N) {
    inits <- c(
      p0 = unname(eq$p[idx]),
      u0 = unname(eq$u[idx]),
      k0 = unname(eq$k[idx]),
      h0 = unname(eq$h[idx]),
      xc0 = unname(eq$xc[idx]),
      x0 = unname(eq$x[idx]),
      L0 = unname(eq$L[idx]),
      CBT_min0 = unname(eq$CBT_min[idx])
    )
  } else if (dynamic_C & dynamic_vars$approx_N) {
    inits <- c(
      p0 = unname(eq$p[idx]),
      u0 = unname(eq$u[idx]),
      k0 = unname(eq$k[idx]),
      h0 = unname(eq$h[idx]),
      xc0 = unname(eq$xc[idx]),
      x0 = unname(eq$x[idx]),
      CBT_min0 = unname(eq$CBT_min[idx])
    )
  } else {
    inits <- c(
      p0 = unname(eq$p[idx]),
      u0 = unname(eq$u[idx]),
      k0 = unname(eq$k[idx]),
      h0 = unname(eq$h[idx])
    )
  }
  
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
    series_start_time = 23.5,
    ndays = 30,
    tz = NULL,
    round_minutes = 30,
    pvec = mccauley2024_pvec,
    pvec_C = FJK_pvec,
    dynamic_C = FALSE,
    dynamic_vars=fips_default_dynamic_vars(model = "FJK")
) {
  dynamic_vars <- fips_normalize_dynamic_vars(dynamic_vars)
  inits <- mccauley2024_make_inits(
    sim_start = sim_start,
    sleep_hrs = sleep_hrs,
    wake_time = wake_time,
    series_start_time = series_start_time,
    ndays = ndays,
    tz = tz,
    round_minutes = round_minutes,
    pvec = pvec,
    pvec_C = pvec_C,
    dynamic_C = dynamic_C,
    dynamic_vars=dynamic_vars,
    return_equilibrium = FALSE
  )
  pvec[c("p0", "u0", "k0", "h0")] <- inits[c("p0", "u0", "k0", "h0")]
  if (dynamic_C & !dynamic_vars$approx_N) {
    pvec_C[c("xc0","x0","L0","CBT_min0")] = inits[c("xc0","x0","L0","CBT_min0")]
  } else if (dynamic_C & dynamic_vars$approx_N) {
    pvec_C[c("xc0","x0","CBT_min0")] = inits[c("xc0","x0","CBT_min0")]
  }
  return(list(pvec = pvec, pvec_C = pvec_C))
}

