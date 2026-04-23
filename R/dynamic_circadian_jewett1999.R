## Dynamic Oscillator for Circadian Pacemaker
# Primary equations reflect the Jewett et al (1999) formulation of the model, with St Hilaire (2007) adjustments 
# to Process L, and non-photic component, and parameters (alpha0, Beta, and G). CBT_min adjustment value of 0.97 also
# taken from St Hilaire (2007). value of pow taken from Kronauer (2000) (see Rea et al 2022)

K99_check_pvec <- function(pvec) {
  accepted_names <- c("mu","q","rho","tau_c","kappa","beta",
                      "G","K","alpha0", "I0","pow","x0","xc0","L0","rho_n","CBT_min0","phi_adj")
  diffvals <- setdiff(names(pvec), accepted_names)
  if (length(diffvals) > 0) {
    msg <- sprintf("K99_check_pvec model halted!:\n [%s] \n is/are unsupported parameters\nPlease remove these before continuing.",
                   diffvals)
    stop(call. = FALSE, msg)
  }
  if (!all(accepted_names %in% names(pvec))) {
    stop(call. = FALSE,
         "K99_check_pvec model halted: missing parameters in supplied pvec.\nSee help(ODEM_make_pvec) for required names.")
  }
  TRUE
}

#' @export
#'
K99_make_pvec <- function(
			mu = 0.13, q = 1/3, rho = 0.99729, tau_c = 24.2, kappa = 0.55,
      G = 37, K = 0.4, alpha0 = 0.1, I0 = 9500, pow = 0.6, beta = 0.007, x0=-0.25,xc0=0.5,L0=.3, rho_n = 0.032,CBT_min0=4,phi_adj=0.97
) {
  pvec <- c(mu=mu, q=q, rho=rho, tau_c=tau_c, kappa=kappa, beta=beta,
            G=G, K=K, alpha0=alpha0, I0=I0, pow=pow,
            x0=x0, xc0=xc0,L0=L0,rho_n=rho_n,CBT_min0=CBT_min0,phi_adj=phi_adj)
  K99_check_pvec(pvec)
  pvec
}

#' @export
#'
K99_pvec <- K99_make_pvec()

# Internal helper to append model columns to a FIPS_df
K99_cols  <- c("x","xc","L")

.K99_rootfun <- function(t, state, parms, df, ...) {
    cur_phase = atan2(state["xc"], state["x"])
    deriv = -2.98 - cur_phase # cross the CBT-min threshold (Ref Woelders et al 2017, citing May et al 2002)
    return(deriv)
}

.K99_eventfun <- function(t, state, parms, df, ...) {
  cur_phase = atan2(state["xc"], state["x"])
  if (cur_phase < 0) { # needed to ensure we don't pick up the flip to the positive angle
    state["CBT_min"] = (t %% 24) + parms["phi_adj"]
  }
  return(state)
}


.K99_derivs <- function(t, state, parms, df, dynamic_vars = NULL) {
  with(as.list(c(state, parms)), {
    
    j <- findInterval(t, df$t_abs, rightmost.closed = TRUE) # use the row corresponding to t<=t_abs
    if ("S" %in% names(state)) { # if we're jointly fitting the Phillips model, use the actual S state
      sleeping=S
    } else { # fall-back to sleep=in_bed
      sleeping = !df$out_of_bed[j]
    }
    light_col <- if ("light_drive" %in% names(df)) "light_drive" else "light_lux"
    I_t <- if (sleeping) 0 else df[[light_col]][j]
    f_t <- if (sleeping) 0 else 1
    
    aI <- (alpha0 * (I_t / I0)^pow) * (I_t/(I_t+100))
    # 3D (A, xc, L) version; if L not present, treat as 2D
    if ("L" %in% names(state)) {
      B <- G * aI * f_t * (1 - L) * (1 - K * xc) * (1 - K * x)
      w2 <- (24 / (rho * tau_c))^2
      # St Hilaire (2007) non-photic component
      if (as.numeric(sleeping) == 1 & (((t%%24 - CBT_min) %% 24) < 16.5 | ((t - CBT_min) %% 24) > 21)) {
          N_s_hat <- rho_n * (1 / 3 - 1)
      } else {
         N_s_hat <- rho_n * (1 / 3)
      }
      N_s    <- N_s_hat * (1 - tanh(10 * x)) 
      dx <- (pi / 12) * (xc + mu * ((x / 3) + (4*x^3)/3 - (256/105) * x^7 ) + B + N_s)
      dxc <- (pi/12) * (q*B*xc - ( w2 + kappa*B)*x)
      dL <- 60 * (aI * f_t * (1 - L) - beta * L)

      return(list(c(dxc, dx, dL, 0)))
    } else {
      n_inf <- aI / (beta + aI)
      if (as.numeric(sleeping) == 1 & (((t - CBT_min) %% 24) < 16.5 | ((t - CBT_min) %% 24) > 21)) {
          N_s_hat <- rho_n * (1 / 3 - 1)
      } else {
         N_s_hat <- rho_n * (1 / 3)
      }
      N_s    <- N_s_hat * (1 - tanh(10 * x)) 
      B <- G * aI * f_t * (1 - n_inf) * (1 - K * xc) * (1 - K * x)
      w2 <- (24 / (rho * tau_c))^2
      dx <- (pi / 12) * (xc + mu * ((x / 3) + (4*x^3)/3 - (256/105) * x^7 ) + B + N_s)
      dxc <- (pi/12) * (q*B*xc - ( w2 + kappa*B)*x)
      return(list(c(dxc, dx, 0)))
    }
  })
}

K99_simulation_dispatch  <- function(dat, pvec, model_formula = NULL, 
                                     dynamic_vars=fips_default_dynamic_vars(model = "K99"),
                                     method="lsoda") {
  dynamic_vars <- fips_normalize_dynamic_vars(dynamic_vars, model = "K99")
  pvec <- .fips_adjust_pvec_for_light_metric(pvec, dynamic_vars, model_name = "K99")
  K99_check_pvec(pvec)

  # Prepare time and schedule columns
  dat <- dat %>% dplyr::rename(out_of_bed = wake_status)
  dat$t_abs <- dat$sim_hours + dat$time[1]
  dat <- .fips_build_light_drive(dat, dynamic_vars)

  # Toggle between Creaser (2021) 2D/3D representation for L state
  if (dynamic_vars$approx_N) {
    y0 = c(xc = pvec[["xc0"]], x = pvec[["x0"]], CBT_min=pvec[["CBT_min0"]])
  } else {
    y0 = c(xc = pvec[["xc0"]], x = pvec[["x0"]], L = pvec[["L0"]], CBT_min=pvec[["CBT_min0"]])
  }
  # Integrate the ODE
  ODE_results = as.data.frame(deSolve::odeode(
    y = y0,
    times = dat$t_abs,
    func = .K99_derivs,
    parms = pvec, method = method, df = dat, dynamic_vars = dynamic_vars, rootfun = .K99_rootfun,
    events = list(func = .K99_eventfun, root = TRUE, maxroot = 1e6)
  )) %>% dplyr::rename(t_abs = time)
  
   # Merge back to data
   dat <- dat %>%
       dplyr::left_join(ODE_results, by = "t_abs") %>%
     dplyr::rename(wake_status=out_of_bed)
  dat <- FIPS_simulation(dat, modeltype = "K99", pvec = pvec, pred_stat = "x", pred_cols = K99_cols)
  dat
}

#' @export
#'
K99_make_inits <- function(
    sim_start,
    sleep_hrs = 8,
    wake_time = "07:30:00",
    series_start_time = 23.5,
    ndays = 30,
    tz = NULL,
    round_minutes = 5,
    pvec = K99_pvec,
    dynamic_vars=fips_default_dynamic_vars(model = "K99"),
    return_equilibrium = FALSE,
    return_adaptation = FALSE,
    save_equilibrium_path = NULL
) {
  dynamic_vars <- fips_normalize_dynamic_vars(dynamic_vars, model = "K99")
  stopifnot(inherits(sim_start, "POSIXct"))
  if (is.null(tz)) tz <- attr(sim_start, "tz") %||% "UTC"
  if (!K99_check_pvec(pvec)) stop("Invalid `pvec`.")
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
    modeltype = "K99",
    pvec      = pvec,
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
  if (dynamic_vars$approx_N) {
    inits <- c(
      x0 = unname(eq$x[idx]),
      xc0 = unname(eq$xc[idx]),
      CBT_min0 = unname(eq$CBT_min[idx])
    )
  } else {
    inits <- c(
      x0 = unname(eq$x[idx]),
      xc0 = unname(eq$xc[idx]),
      L0 = unname(eq$L[idx]),
      CBT_min0 = unname(eq$CBT_min[idx])
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

#' @export
K99_init_pvec <- function(
    sim_start,
    sleep_hrs = 8,
    wake_time = "07:30:00",
    series_start_time = 23.5,
    ndays = 30,
    tz = NULL,
    round_minutes = 30,
    pvec = K99_pvec, 
    dynamic_vars=fips_default_dynamic_vars(model = "K99")
) {
  dynamic_vars <- fips_normalize_dynamic_vars(dynamic_vars, model = "K99")
  inits <- K99_make_inits(
    sim_start = sim_start,
    sleep_hrs = sleep_hrs,
    wake_time = wake_time,
    series_start_time = series_start_time,
    ndays = ndays,
    tz = tz,
    round_minutes = round_minutes,
    pvec = pvec,
    return_equilibrium = FALSE,
    dynamic_vars=dynamic_vars
  )
  if (dynamic_vars$approx_N) {
    pvec[c("x0", "xc0","CBT_min0")] <- inits[c("x0", "xc0","CBT_min0")]
  } else {
    pvec[c("x0", "xc0", "L0", "CBT_min0")] <- inits[c("x0", "xc0", "L0", "CBT_min0")]
  }
  pvec
}


