phillips2017_check_pvec <- function(pvec) {
  accepted_names <- c("mu_wake","mu_sleep","chi_wake","chi_sleep","lambda","gammaR",
                      "Kd1","Kd2","R1_tot_0", "R2_tot","D_wake","D_sleep",
                      "a","omega","phi","A_tot_0",
                      "D_mid","D_s","p_max")
  diffvals <- setdiff(names(pvec), accepted_names)
  if (length(diffvals) > 0) {
    msg <- sprintf("phillips2017_check_pvec model halted!:\n [%s] \n is/are unsupported parameters\nPlease remove these before continuing.",
                   diffvals)
    stop(call. = FALSE, msg)
  }
  if (!all(accepted_names %in% names(pvec))) {
    stop(call. = FALSE,
         "phillips2017_check_pvec model halted: missing parameters in supplied pvec.\nSee help(ODEM_make_pvec) for required names.")
  }
  TRUE
}

# Single-harmonic circadian, time t in hours (absolute), phi in hours
phillips2017_Cfun <- function(t, phi = 7.95, a = 3.25) {
  a * cos(((2 * pi) / 24) * (t - phi))
}

#' @export
#'
phillips2024_lux <- function(t,s1=7.5,s2=16.5,c=0.6,sunrise=6,sunset=19) {
  out = numeric(length(t))
  mask = t>sunrise & t< sunset
  out[!mask]=50
  out[mask] = 150 * tanh(0.6 * (t[mask] - s1)) - 150 * tanh(c * (t[mask] - s2))
  out
}


#' @export
#'
phillips2017_make_pvec <- function(
						  mu_wake = 869.5, mu_sleep = 596.4,
						  chi_wake = 18.18, chi_sleep = 4.20,
						  lambda = 291, gammaR = 0.9677,
						  Kd1 = 1, Kd2 = 100,
						  R1_tot_0 = 582, R2_tot = 300,A_tot_0 = 711.5,
						  a = 3.25, omega = 2*pi/24, phi = 7.95,
						  D_mid = 583.2, D_s = 5.872, p_max = 60,
						  D_wake=555.4,D_sleep=572.7
) {
  pvec <- c(mu_wake=mu_wake, mu_sleep=mu_sleep, chi_wake=chi_wake, chi_sleep=chi_sleep,lambda=lambda,gammaR=gammaR,
            Kd1=Kd1, Kd2=Kd2,R1_tot_0=R1_tot_0, R2_tot=R2_tot,a=a, omega=omega, phi=phi,
            D_mid=D_mid, D_s=D_s, p_max=p_max,A_tot_0=A_tot_0,D_wake=D_wake,D_sleep=D_sleep)
  phillips2017_check_pvec(pvec)
  pvec
}

#' @export
#'
phillips2017_pvec <- phillips2017_make_pvec()


# Internal helper to append model columns to a FIPS_df
phillips2017_cols  <- c("p", "A_tot", "A_u","R1_tot","R1_b", "D", "C","S")

# Binding helper identical to Phillips, without the C term (for reuse)
.phillips2017_compute_binding_noC <- function(A_tot, R1_tot, pvec) {
  with(as.list(pvec), {
    # Approximate B = R2_u/(R2_u + Kd2) assuming R2_u ≪ R2_tot
    B <- R2_tot / (R2_tot + Kd2)
    term <- A_tot + R1_tot + (Kd1 / (1 - B))
    discriminant <- term^2 - 4 * A_tot * R1_tot
    if(any(discriminant<0)) stop("Discriminant cannot be <0")
    A1_b <- 0.5 * (term - sqrt(discriminant))
    A2_b <- B * (A_tot - A1_b)
    A_u  <- A_tot - A1_b - A2_b
    list(A1_b = A1_b, A2_b = A2_b, A_u = A_u)
  })
}

# Phillips binding but allowing a direct C_value override
.phillips2017_compute_binding_with_C <- function(A_tot, R1_tot, pvec, C_value) {
  b <- .phillips2017_compute_binding_noC(A_tot = A_tot, R1_tot = R1_tot, pvec = pvec)
  D <- b$A1_b + C_value
  list(A1_b = b$A1_b, A2_b = b$A2_b, A_u = b$A_u, D = D, c_t = C_value)
}

.phillips2017_compute_C <- function(y,dynamic_vars) {
  if (dynamic_vars$model=="K99") {
    C = -0.608 * y[["x"]] + 0.685 * y[["xc"]] + 0.075
  } else if (dynamic_vars$model=="FJK") {
    C = -0.569 * y[["x"]] + 0.708 * y[["xc"]] + 0.075
  }
  return(C)
}

# Root functions: return 0 at switching moments
#  g1: fall asleep (only allowed in bed) when D rises through D_on
#  g2: wake when D falls through D_off
#  g3: check transition to/from bed (either wake up, or check if D>threshold upon entry to bed)
.phillips2017_rootfun <- function(t, y, parms, pvec_C, df, change_points,dynamic_vars,...) {
  g_sleep <- 1
  g_wake <- 1
  g_change <- as.numeric(all((t - change_points) != 0))
  g_CBTmin <- 1
  j <- findInterval(t, df$t_abs, rightmost.closed = TRUE) # use the row corresponding to t<=t_abs
  out_of_bed = df$out_of_bed[j]
  S <- y["S"]
  if ("x" %in% names(y)) { # if we're fitting the dynamic state variable, use that
    c_t =  parms["a"]*.phillips2017_compute_C(y, dynamic_vars)
    if (dynamic_vars$model == "K99") {
      cur_phase = atan2(y["xc"], y["x"])
      g_CBTmin = -2.98 - cur_phase # cross the CBT-min threshold (Ref Woelders et al 2017, citing May et al 2002)
    } else if (dynamic_vars$model=="FJK") {
      g_CBTmin = unlist(.FJK_derivs(t, y, pvec_C, df, dynamic_vars = dynamic_vars))[2]
    }
  } else { # fallback to static function from Phillips (2017)
    c_t = phillips2017_Cfun(t,parms[["phi"]],parms[["a"]])
  }

  if (!out_of_bed & g_change != 0) {
    D <- unname(.phillips2017_compute_binding_with_C(y["A_tot"], y["R1_tot"], parms, c_t)$D)
    if (S < 0.5) {
      g_sleep = D - parms[["D_sleep"]] # root 1: sleep onset
    } else {
      g_wake = D - parms[["D_wake"]] # root 2: wake
    }
  }
  
  return(c(g_sleep,g_wake,g_change,g_CBTmin))
}

## This is the event handler for when the root function is triggered. It checks sleep/wake traansitions:
# 1. When transitioning to bed (check if D>sleep threshold)
# 2. When tranisitioning out of bed (wake up)
# 3. While in bed, check whether sleep/wake thresholds are crossed and update state accordingly
.phillips2017_eventfun <- function(t, y, pvec, pvec_C,df, change_points,dynamic_vars,...) {
  j <- findInterval(t, df$t_abs, rightmost.closed = TRUE) # use the row corresponding to t<=t_abs
  out_of_bed = df$out_of_bed[j]
  g = .phillips2017_rootfun(t, y, pvec, pvec_C, df, change_points, dynamic_vars)
  if ("x" %in% names(y)) { # if we're fitting the dynamic state variable, use that
    c_t =  pvec["a"]*.phillips2017_compute_C(y,dynamic_vars)
  } else { # fallback to static function from Phillips (2017)
    c_t = phillips2017_Cfun(t,pvec[["phi"]],pvec[["a"]])
  }
  D <- unname(.phillips2017_compute_binding_with_C(y["A_tot"], y["R1_tot"], pvec, c_t)$D)
  S <- y["S"]
  if (g[3] == 0) {
    if (df$switch_direction[j] == "Wake" & S == 1) {
      y["S"] = 0
    } else if (df$switch_direction[j] == "Sleep" & D >= pvec[["D_sleep"]]) {
      y["S"] = 1
    }
  } else if (abs(g[1]) < 1e-4 | abs(g[2]) < 1e-4) {
    if (S == 1 & abs(D - pvec[["D_wake"]]) < 1) {
      y["S"] = 0
    } else if (S == 0 & (abs(D - pvec[["D_sleep"]]) < 1)) {
      y["S"] = 1
    }
  } else if (g[4]< 1e-4) {
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
.phillips2017_derivs <- function(t,state,parms,...) {
  with(as.list(c(state, parms)), {
    # Binding uses only A_tot, R1_tot; C contributes to D = A1_b + C_value
    bind <- .phillips2017_compute_binding_noC(A_tot = A_tot, R1_tot = R1_tot, pvec = parms)
    mu  <- ifelse(S == 1, mu_sleep, mu_wake)
    chi <- ifelse(S == 1, chi_sleep, chi_wake)

    dA_tot <- (mu - A_tot) / chi
    dR1_tot <- (bind$A1_b - gammaR * R1_tot) / lambda
    # dS is handled via events; derivative is 0 between events
    c(unname(dA_tot), unname(dR1_tot), 0)
  })
}

.phillips2017_add_states <- function(dat, pvec,dynamic_vars,pvec_C) {
  with(as.list(c(dat, pvec)), {
    # Check sleep state, only allow transition to sleep if in bed
    if (!("x" %in% names(dat))) { # if we're fitting the dynamic state variable, use that
      dat$C = phillips2017_Cfun(dat$t_abs,phi,a)
      c_t = dat$C
    } else {
      c_t =  pvec["a"]*.phillips2017_compute_C(list(xc = dat$xc, x = dat$x), dynamic_vars)
      dat$C=c_t
    }
    binding  <- .phillips2017_compute_binding_with_C(A_tot, R1_tot, pvec, c_t)
    p <- p_max / (1 + exp((D_mid - binding$D) / D_s))
    dat$p = p
    dat$D = binding$D
    dat$R1_b = binding$A1_b
    dat$A_u = binding$A_u
    dat
  })
}

## Coupling for Phillips (2017) model with dynamic circadian (FJK)
## This keeps each subsystem modular while integrating as one ODE when needed.
# Coupled function that shares S and C across subsystems if desired
.phillips2017_derivs_coupled <- function(t, y, parms, df, pvec_C,use_dynamic,dynamic_vars,...) {
   # Phillips derivatives depend on dynamic C
   d_ph <- .phillips2017_derivs(
     t = t,
     state = y,
     parms = parms
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
     return(list(d_ph))
  }
  return(list(c(d_ph, d_c)))
}

# - If use_dynamic = FALSE, falls back to Phillips-only with static C(t)
# - If use_dynamic = TRUE, integrates Phillips + FJK together, where
#   Phillips depends on C(t) from FJK, and FJK can depend on S (sleep state)
#   from Phillips rather than the schedule.
#   NB I use the transformation from Phillips et al (2023)
phillips2017_simulation_dispatch <- function(
  dat,
  pvec=phillips2017_pvec,
  pvec_C=FJK_pvec,
  method = "lsoda",
  use_dynamic = TRUE,
  dynamic_vars = fips_default_dynamic_vars(model = "FJK"),
  model_formula=NULL
) {
  dynamic_vars <- fips_normalize_dynamic_vars(dynamic_vars)

  # Ensure required inputs exist
  phillips2017_check_pvec(pvec)
  if (dynamic_vars$model == "FJK") {
    FJK_check_pvec(pvec_C)
  } else if (dynamic_vars$model=="K99") {
    K99_check_pvec(pvec_C)
  }
  pvec_C <- .fips_adjust_pvec_for_light_metric(pvec_C, dynamic_vars, model_name = dynamic_vars$model)

  # Prepare time and schedule columns
  dat <- dat %>% dplyr::rename(out_of_bed = wake_status)
  dat$t_abs <- dat$sim_hours + dat$time[1]
  change_points <- dat$t_abs[dat$change_point == 1]

  dat <- .fips_build_light_drive(dat, dynamic_vars)

  # Initial S: use Phillips thresholds with initial C from circadian pvec if dynamic_C, otherwise static c_t from Phillips (2017)
  c_t = ifelse(use_dynamic, .phillips2017_compute_C(c(x = pvec_C[["x0"]], xc = pvec_C[["xc0"]]), dynamic_vars),
    phillips2017_Cfun(dat$time[1], pvec[["phi"]], pvec[["a"]])
  )
  binding0 <- .phillips2017_compute_binding_with_C(
    A_tot = pvec[["A_tot_0"]],
    R1_tot = pvec[["R1_tot_0"]],
    pvec = pvec,
    C_value = c_t
  )
  S0 <- unname(ifelse(binding0$D >= pvec[["D_sleep"]], 1, 0))

  y0 <- c(
    A_tot = pvec[["A_tot_0"]],
    R1_tot = pvec[["R1_tot_0"]],
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
    func = .phillips2017_derivs_coupled,
    rootfun = .phillips2017_rootfun,
    parms = pvec,
    method = method,
    events = list(func = .phillips2017_eventfun, root = TRUE, maxroot = 1e6),
    df = dat,
    change_points = change_points,
    pvec_C = pvec_C,
    use_dynamic = use_dynamic,
    dynamic_vars=dynamic_vars
  )) %>% dplyr::rename(t_abs = time)
  # Merge back to data; compute Phillips outputs using dynamic C from the path
  dat <- dat %>%
    dplyr::left_join(ODE_results, by = "t_abs") %>%
    dplyr::rename(wake_status=out_of_bed)
  
  dat <- .phillips2017_add_states(dat, pvec, dynamic_vars, pvec_C)

  # Keep downstream processing consistent with Phillips interface
  if (use_dynamic & dynamic_vars$model == "FJK") {
    pred_cols = unique(c(phillips2017_cols, FJK_cols))
  } else if (use_dynamic & dynamic_vars$model == "K99") {
    pred_cols = unique(c(phillips2017_cols, K99_cols))
  } else {
    pred_cols = phillips2017_cols
  }
  dat <- FIPS_simulation(dat,
    modeltype = "phillips2017", pvec = pvec, pred_stat = "fatigue",
    pred_cols = pred_cols
  )
  if (is.null(model_formula)) {
    dat <- process_bmm_formula(dat, fatigue~p, pvec)
  } else {
    dat <- process_bmm_formula(dat, model_formula, pvec)
  }
  dat
}

#' @export
#'
phillips2017_make_inits <- function(
    sim_start,
    sleep_hrs = 8,
    wake_time = "07:30:00",
    series_start_time = 23.5,
    ndays = 30,
    tz = NULL,
    round_minutes = 5,
    fast_adapt = FALSE,
    pvec = phillips2017_pvec,
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
  if (!phillips2017_check_pvec(pvec)) stop("Invalid `pvec`.")

  # Anchor the adaptation period `ndays` before sim_start, at the daily wake time.
  sim_start_date <- lubridate::as_date(sim_start, tz = tz)
  anchor_wake <- lubridate::ymd_hms(
    paste0(sim_start_date - lubridate::days(ndays), " ", wake_time),
    tz = tz
  )
  attach(as.list(pvec))
    T_awake = sleep_hrs-24
    term_wake=mu_wake*(1-exp(T_awake/chi_wake)) + mu_sleep*(1-exp(-sleep_hrs/chi_sleep))*exp(T_awake/chi_wake)
    term_sleep=mu_sleep*(1-exp(-sleep_hrs/chi_sleep)) + mu_wake*(1-exp(T_awake/chi_wake))*exp(-sleep_hrs/chi_sleep)
    denom=(1-exp((T_awake/chi_wake)*(sleep_hrs/chi_sleep)))
  detach()
  pvec["A_tot_0"] = term_wake / denom
  if (fast_adapt) { # set lambda super short for initialisation
    pvec["lambda"] = 48
  } 
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
    modeltype = "phillips2017",
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
      A_tot_0 = unname(eq$A_tot[idx]),
      R1_tot_0 = unname(eq$R1_tot[idx]),
      xc0 = unname(eq$xc[idx]),
      x0 = unname(eq$x[idx]),
      L0 = unname(eq$L[idx]),
      CBT_min0 = unname(eq$CBT_min[idx])
    )
  } else if (dynamic_C & dynamic_vars$approx_N) {
    inits <- c(
      A_tot_0 = unname(eq$A_tot[idx]),
      R1_tot_0 = unname(eq$R1_tot[idx]),
      xc0 = unname(eq$xc[idx]),
      x0 = unname(eq$x[idx]),
      CBT_min0 = unname(eq$CBT_min[idx])
    )
  } else {
    inits <- c(
      A_tot_0 = unname(eq$A_tot[idx]),
      R1_tot_0 = unname(eq$R1_tot[idx])
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
#' 
phillips2017_init_pvec <- function(
    sim_start,
    sleep_hrs = 8,
    wake_time = "07:30:00",
    series_start_time = 23.5,
    ndays = 30,
    tz = NULL,
    round_minutes = 30,
    pvec = phillips2017_pvec,
    pvec_C = FJK_pvec,
    dynamic_C = FALSE,
    fast_adapt = TRUE,
    dynamic_vars=fips_default_dynamic_vars(model = "FJK")
) {
  dynamic_vars <- fips_normalize_dynamic_vars(dynamic_vars)
  inits <- phillips2017_make_inits(
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
    fast_adapt=fast_adapt,
    return_equilibrium = FALSE
  )
  pvec[c("A_tot_0", "R1_tot_0")] <- inits[c("A_tot_0", "R1_tot_0")]
  if (dynamic_C & !dynamic_vars$approx_N) {
    pvec_C[c("xc0","x0","L0","CBT_min0")] = inits[c("xc0","x0","L0","CBT_min0")]
  } else if (dynamic_C & dynamic_vars$approx_N) {
    pvec_C[c("xc0","x0","CBT_min0")] = inits[c("xc0","x0","CBT_min0")]
  }
  return(list(pvec = pvec, pvec_C = pvec_C))
}


