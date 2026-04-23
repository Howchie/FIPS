## Dynamic Circadian Oscillator function
# Primary equations from Forger et al (1999)
# Parameter adjustments, and 2D reduction, from Creaser et al (2021) (see also O'Diekman, 2022)

FJK_check_pvec <- function(pvec) {
  accepted_names <- c("mu","q","rho","tau_c","kappa","beta",
                      "G","K","alpha0", "I0","pow","xc0","x0","L0","CBT_min0","phi_adj")
  diffvals <- setdiff(names(pvec), accepted_names)
  if (length(diffvals) > 0) {
    msg <- sprintf("FJK_check_pvec model halted!:\n [%s] \n is/are unsupported parameters\nPlease remove these before continuing.",
                   diffvals)
    stop(call. = FALSE, msg)
  }
  if (!all(accepted_names %in% names(pvec))) {
    stop(call. = FALSE,
         "FJK_check_pvec model halted: missing parameters in supplied pvec.\nSee help(ODEM_make_pvec) for required names.")
  }
  TRUE
}

#' @export
#'
FJK_make_pvec <- function(
			mu = 0.23, q = 4/3, rho = 0.99669, tau_c = 24.2, kappa = 0.55,
      G = 33.75, K = 0.4, alpha0 = 0.05, I0 = 9500, pow = 0.5, beta = 0.0075, x0=-0.25,xc0=0.5,L0=.3,CBT_min0=4,phi_adj=0
) {
  pvec <- c(mu=mu, q=q, rho=rho, tau_c=tau_c, kappa=kappa, beta=beta,
            G=G, K=K, alpha0=alpha0, I0=I0, pow=pow,
            xc0=xc0, x0=x0,L0=L0,CBT_min0=CBT_min0,phi_adj=phi_adj)
  FJK_check_pvec(pvec)
  pvec
}

#' @export
#'
FJK_pvec <- FJK_make_pvec()

# Internal helper to append model columns to a FIPS_df
FJK_cols  <- c("xc","x","L")

.FJK_rootfun <- function(t, state, parms, df) {
    deriv = unlist(.FJK_derivs(t, state, parms, df))
    return(deriv[2]) 
}

.FJK_eventfun <- function(t, state, parms, df) {
  if (state["x"] < 0) { # only adjust for the derivative for state x<0, where the minimum is CBT
    state["CBT_min"] = (t %% 24)
  }
  return(state)
}


.FJK_derivs <- function(t, state, parms,df) {
  with(as.list(c(state, parms)), {
    
    j <- findInterval(t, df$t_abs, rightmost.closed = TRUE) # use the row corresponding to t<=t_abs
    if ("S" %in% names(state)) { # if we're jointly fitting the Phillips model, use the actual S state
      sleeping=S
    } else { # fall-back to sleep=in_bed
      sleeping = !df$out_of_bed[j]
    }
    I_t <- if (sleeping) 0 else df$light_lux[j]
    f_t <- if (sleeping) 0 else 1
    aI <- alpha0 * (I_t / I0)^pow * (I_t/(I_t+100))
    # 3D (xc, x, L) version; if L not present, treat as 2D
    if ("L" %in% names(state)) {
      B <- G * aI * f_t * (1 - L) * (1 - K * x) * (1 - K * xc)
      w2 <- (24 / (rho * tau_c))^2
      dxc <- (pi/12) * (mu * (xc - q * xc^3) - x * (w2 + kappa * B))
      dx <- (pi/12) * (xc + B)
      dL <- 60 * (aI * f_t * (1 - L) - beta * L)
      return(list(c(dxc, dx, dL, 0)))
    } else {
      n_inf <- aI / (beta + aI)
      B <- G * aI * f_t * (1 - n_inf) * (1 - K * x) * (1 - K * xc)
      w2 <- (24 / (rho * tau_c))^2
      dxc <- (pi/12) * (mu * (xc - q * xc^3) - x * (w2 + kappa * B))
      dx <- (pi/12) * (xc + B)
      return(list(c(dxc, dx, 0)))
    }
  })
}

FJK_simulation_dispatch  <- function(dat, pvec, model_formula = NULL, 
                                     dynamic_vars=list(approx_N=FALSE,model="FJK",light_levels=lux_levels()),
                                     method="lsoda") {
  FJK_check_pvec(pvec)

  # Prepare time and schedule columns
  dat <- dat %>% dplyr::rename(out_of_bed = wake_status)
  dat$t_abs <- dat$sim_hours + dat$time[1]
  if (!"light_lux" %in% names(dat)) {
    light_levels = dynamic_vars$light_levels
    if(!"func"%in%names(light_levels)) {light_levels$func=NULL}
    dat = add_light_lux(dat,level_windows=light_levels$windows, fallback=light_levels$fallback,func=light_levels$func)
  }
  # Toggle between Creaser (2021) 2D/3D representation for L state
  if (dynamic_vars$approx_N) {
    y0 = c(xc = pvec[["xc0"]], x = pvec[["x0"]], CBT_min = pvec[["CBT_min0"]])
  } else {
    y0 = c(xc = pvec[["xc0"]], x = pvec[["x0"]], L = pvec[["L0"]], CBT_min = pvec[["CBT_min0"]])
  }
  
  # Integrate the ODE
  ODE_results = as.data.frame(deSolve::ode(
    y = y0,
    times = dat$t_abs,
    func = .FJK_derivs,
    parms = pvec, method = method, df = dat, rootfun = .FJK_rootfun,
    events = list(func = .FJK_eventfun, root = TRUE, maxroot = 1e6)
  )) %>% dplyr::rename(t_abs = time)
  
   # Merge back to data
  dat <- dat %>%
    dplyr::left_join(ODE_results,by="t_abs") %>%
       dplyr::rename(wake_status=out_of_bed)
  dat <- FIPS_simulation(dat, modeltype = "FJK", pvec = pvec, pred_stat = "x", pred_cols = FJK_cols)
  dat
}

#' @export
#'
FJK_make_inits <- function(
    sim_start,
    sleep_hrs = 8,
    wake_time = "07:30:00",
    series_start_time = 23.5,
    ndays = 30,
    tz = "Australia/Perth",
    round_minutes = 5,
    pvec = FJK_pvec,
    dynamic_vars=list(approx_N=FALSE,model="FJK",light_levels=lux_levels()),
    return_equilibrium = FALSE,
    return_adaptation = FALSE,
    save_equilibrium_path = NULL
) {
  stopifnot(inherits(sim_start, "POSIXct"))
  if (is.null(tz)) tz <- attr(sim_start, "tz") %||% "UTC"
  if (!FJK_check_pvec(pvec)) stop("Invalid `pvec`.")

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
    modeltype = "FJK",
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
FJK_init_pvec <- function(
    sim_start,
    sleep_hrs = 8,
    wake_time = "07:30:00",
    series_start_time = 23.5,
    ndays = 30,
    tz = NULL,
    round_minutes = 30,
    pvec = FJK_pvec,
    dynamic_vars=list(approx_N=FALSE,model="FJK",light_levels=lux_levels())
) {
  inits <- FJK_make_inits(
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
    pvec[c("xc0", "x0","CBT_min0")] <- inits[c("xc0", "x0","CBT_min0")]
  } else {
    pvec[c("xc0", "x0", "L0", "CBT_min0")] <- inits[c("xc0", "x0", "L0", "CBT_min0")]
  }
  pvec
}


