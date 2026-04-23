
#' Check McCauley (2013) Parameter Vector
#'
#' Validates that the supplied parameter vector contains the required entries
#' for the McCauley (2013) model.
#'
#' @param pvec Named numeric vector of model parameters.
#'
#' @return logical
mccauley2013_check_pvec <- function(pvec) {
  required <- c(
    "alpha_w", "alpha_s", "beta_w", "beta_s", "eta_w",
    "Tp", "mu_w", "mu_s", "phi", "lambda_w", "lambda_s",
    "xi", "p0", "u0", "k0"
  )
  optional <- c("eta_s", "Ac", "Wc")
  accepted <- c(required, optional)

  diffvals <- setdiff(names(pvec), accepted)
  if (length(diffvals) > 0) {
    msg <- sprintf(
      "mccauley2013_check_pvec model halted!:
 [%s] 
 is/are unsupported parameters
Please remove these before continuing.",
      diffvals
    )
    stop(call. = FALSE, msg)
  }

  missing <- setdiff(required, names(pvec))
  if (length(missing) > 0) {
    stop(
      call. = FALSE,
      sprintf(
        "mccauley2013_check_pvec model halted: missing parameters in supplied pvec.
Required names include: %s",
        paste(missing, collapse = ", ")
      )
    )
  }

  if (!any(c("eta_s", "Ac", "Wc") %in% names(pvec))) {
    stop(call. = FALSE,
         "mccauley2013_check_pvec model halted: supply either `eta_s`, `Ac`, or `Wc` to determine sleep dissipation.")
  }
  TRUE
}

derive_eta_s_2013 <- function(eta_w, Ac) {
  eta_w * (Ac / (Ac - 1))
}

#' Make McCauley (2013) Model Parameter Vector
#'
#' Creates a parameter vector with defaults matching the McCauley et al. (2013)
#' model. Note that the published formulation uses negative coefficients for
#' alpha and beta. Here we retain those published signs so the differential
#' equations can be implemented exactly as in the original study.
#'
#' @param alpha_w Homeostatic build-up rate for performance during wakefulness.
#' @param alpha_s Homeostatic dissipation rate for performance during sleep.
#' @param beta_w Scaling factor for the effect of u(t) during wakefulness.
#' @param beta_s Scaling factor for the effect of u(t) during sleep.
#' @param eta_w Build-up rate for the homeostatic process during wake.
#' @param eta_s Dissipation rate for the homeostatic process during sleep. If
#'   `NA`, it will be derived from `eta_w` and `Ac` via Eq. (4) in McCauley et al. (2013).
#' @param Ac Critical wake ratio (Wc / T). Defaults to 0.84 as reported in the paper.
#' @param Wc Critical duration of wakefulness (hours). Only used when `Ac` is not supplied.
#' @param Tp Period of the wake/sleep cycle (usually 24 h).
#' @param mu_w Circadian offset during wakefulness.
#' @param mu_s Circadian offset during sleep.
#' @param phi Circadian phase (hours).
#' @param lambda_w Rate constant for circadian amplitude during wake.
#' @param lambda_s Rate constant for circadian amplitude during sleep (typically negative).
#' @param xi Asymptotic maximum of the circadian modulation.
#' @param p0,u0,k0 Initial conditions for the state variables.
#'
#' @export
mccauley2013_make_pvec <- function(
    alpha_w = -0.028,
    alpha_s = -0.260,
    beta_w = -0.260,
    beta_s = -0.260,
    eta_w = 0.0074,
    eta_s = NA_real_,
    Ac = NA_real_,
    Wc = 20.2,
    Tp = 24,
    mu_w = 0.33,
    mu_s = -1.50,
    phi = 21.2,
    lambda_w = 0.49,
    lambda_s = -0.49,
    xi = 1.09,
    p0 = 5,
    u0 = 25,
    k0 = xi
) {
  if (is.na(Ac)) {
    if (is.na(Wc)) stop("Supply either `Ac` or `Wc` so that eta_s can be derived.")
    Ac <- Wc / Tp
  }
  if (is.na(eta_s)) {
    eta_s <- derive_eta_s_2013(eta_w, Ac)
  }

  pvec <- c(
    alpha_w = alpha_w,
    alpha_s = alpha_s,
    beta_w = beta_w,
    beta_s = beta_s,
    eta_w = eta_w,
    eta_s = eta_s,
    Ac = Ac,
    Wc = Wc,
    Tp = Tp,
    mu_w = mu_w,
    mu_s = mu_s,
    phi = phi,
    lambda_w = lambda_w,
    lambda_s = lambda_s,
    xi = xi,
    p0 = p0,
    u0 = u0,
    k0 = k0
  )

  mccauley2013_check_pvec(pvec)
  pvec
}

#' Default parameter vector for the McCauley (2013) model
#'
#' @export
mccauley2013_pvec <- mccauley2013_make_pvec()

mccauley2013_cols <- c("p", "u", "k")

.mccauley2013_derivs <- function(t, state, parms, df, ...) {
  with(as.list(c(state, parms)), {
    j <- findInterval(t, df$t_abs, rightmost.closed = TRUE) # use the row corresponding to t<=t_abs
    if (!("x" %in% names(state))) {
      C <- mccauley_cfun(t, phi, Tp)
    } else {
      C=-x # McCauley's model uses a peaking at night function
    }

    if (S == 0) {
      g_w <- k * (C + mu_w)
      dp <- alpha_w * (p + beta_w * u) + g_w
      du <- eta_w * u
      dk <- lambda_w * k * (1 - k / xi)
    } else {
      g_s <- k * (-C + mu_s) + (alpha_s * beta_s / eta_s)
      dp <- alpha_s * (p + beta_s * u) + g_s
      du <- eta_s * u + 1
      dk <- lambda_s * k
    }

    c(dp, du, dk, 0)
  })
}

.mccauley2013_derivs_coupled <- function(t, y, parms, df, pvec_C,use_dynamic,dynamic_vars,...) {
   # Phillips derivatives depend on dynamic C
   d_mccauley <- .mccauley2013_derivs(
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
      df = df
      ))
    } else if (dynamic_vars$model=="K99") {
        d_c <- unlist(.K99_derivs(
          t = t,
          state = y,
          parms = pvec_C,
          df = df
        ))
    }
  } else {
     return(list(d_mccauley))
  }
  return(list(c(d_mccauley, d_c)))
}

mccauley2013_simulation_dispatch <- function(
    dat,
    pvec = mccauley2013_pvec,pvec_C=FJK_pvec,
    method = "lsoda",
    use_dynamic = TRUE,
    dynamic_vars = list(approx_N = FALSE,model="FJK"),
    model_formula = NULL
) {
  # Ensure required inputs exist
  if (!"eta_s" %in% names(pvec)) {
    if ("Ac" %in% names(pvec)) {
      pvec["eta_s"] <- derive_eta_s_2013(pvec[["eta_w"]], pvec[["Ac"]])
    } else if ("Wc" %in% names(pvec)) {
      Ac <- pvec[["Wc"]] / pvec[["Tp"]]
      pvec["Ac"] <- Ac
      pvec["eta_s"] <- derive_eta_s_2013(pvec[["eta_w"]], Ac)
    } else {
      stop("`pvec` must include `eta_s`, `Ac`, or `Wc`.")
    }
  }
  mccauley2013_check_pvec(pvec)
  if (dynamic_vars$model == "FJK") {
    FJK_check_pvec(pvec_C)
  } else if (dynamic_vars$model=="K99") {
    K99_check_pvec(pvec_C)
  }
  dat$t_abs <- dat$sim_hours + dat$time[1]
  change_points <- dat$t_abs[dat$change_point == 1]

  # Ensure light is present for FJK
  if (!"light_lux" %in% names(dat)) {
    light_levels = dynamic_vars$light_levels
    if(!"func"%in%names(light_levels)) {light_levels$func=NULL}
    dat = add_light_lux(dat,level_windows=light_levels$windows, fallback=light_levels$fallback,func=light_levels$func)
  }

  S0 <- as.numeric(dat$wake_status[1] == 0)
  y0 <- c(
    p = pvec[["p0"]],
    u = pvec[["u0"]],
    k = pvec[["k0"]],
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

  ODE_results <- as.data.frame(
    deSolve::ode(
      y = y0,
      times = dat$t_abs,
      func = .mccauley2013_derivs_coupled,
      rootfun = .mccauley_rootfun,
      parms = pvec,
      method = method,
      events = list(func = .mccauley_eventfun, root = TRUE, maxroot = 1e6),
      df = dat,
      change_points = change_points,
      pvec_C = pvec_C,
      use_dynamic = use_dynamic,
      dynamic_vars=list(approx_N = dynamic_vars$approx_N, model=dynamic_vars$model)
    )
  ) %>% dplyr::rename(t_abs = time)

  dat <- dat %>% dplyr::left_join(ODE_results, by = "t_abs")
    # Keep downstream processing consistent with Phillips interface
  if (use_dynamic & dynamic_vars$model=="FJK") {
    pred_cols = unique(c(mccauley2013_cols, FJK_cols))
  } else if (use_dynamic & dynamic_vars$model=="K99") {
    pred_cols = unique(c(mccauley2013_cols, K99_cols))
  } else {
    pred_cols = mccauley2013_cols
    dat$C = mccauley_cfun(dat$t_abs,pvec["phi"])
  }

  dat <- FIPS_simulation(
    dat,
    modeltype = "mccauley2013",
    pvec = pvec,
    pred_stat = "fatigue",
    pred_cols = pred_cols
  )

  if (is.null(model_formula)) {
    dat <- process_bmm_formula(dat, fatigue ~ p, pvec)
  } else {
    dat <- process_bmm_formula(dat, model_formula, pvec)
  }
  dat
}

#' Compute steady-state initial conditions for the McCauley (2013) model
#'
#' Uses a regular sleep schedule to find an approximate equilibrium for the
#' McCauley (2013) system.
#'
#' @inheritParams mccauley2024_make_inits
#' @param pvec Parameter vector used for the adaptation run. `p0/u0/k0` in
#'   `pvec` are ignored and replaced by the returned values.
#'
#' @return Named vector `c(p0=..., u0=..., k0=...)` or a list containing
#'   `inits` and `equilibrium`/`adaptation` when requested.
#' @export
mccauley2013_make_inits <- function(
    sim_start,
    sleep_hrs = 8,
    wake_time = "07:30:00",
    series_start_time = 23.5,
    ndays = 30,
    tz = NULL,
    round_minutes = 5,
    pvec = mccauley2013_pvec,
    pvec_C = FJK_pvec,
    dynamic_C = FALSE,
    dynamic_vars=list(approx_N=FALSE,model="FJK"),
    return_equilibrium = FALSE,
    return_adaptation = FALSE,
    save_equilibrium_path = NULL
) {
  stopifnot(inherits(sim_start, "POSIXct"))
  if (is.null(tz)) tz <- attr(sim_start, "tzone") %||% "UTC"
  if (!"eta_s" %in% names(pvec)) {
    if ("Ac" %in% names(pvec)) {
      pvec["eta_s"] <- derive_eta_s_2013(pvec[["eta_w"]], pvec[["Ac"]])
    } else if ("Wc" %in% names(pvec)) {
      Ac <- pvec[["Wc"]] / pvec[["Tp"]]
      pvec["Ac"] <- Ac
      pvec["eta_s"] <- derive_eta_s_2013(pvec[["eta_w"]], Ac)
    } else {
      stop("`pvec` must include `eta_s`, `Ac`, or `Wc`.")
    }
  }
  if (!mccauley2013_check_pvec(pvec)) stop("Invalid `pvec`.")

  sim_start_date <- lubridate::as_date(sim_start, tz = tz)
  anchor_wake <- lubridate::ymd_hms(
    paste0(sim_start_date - lubridate::days(ndays), " ", wake_time),
    tz = tz
  )

  cycles <- tibble::tibble(
    n_days = ndays,
    sleep_hrs = sleep_hrs,
    wake_time = wake_time,
    start_idx = 0L
  )

  sleeptimes <- build_sleeptimes(anchor_wake = anchor_wake, cycles = cycles, tz = tz)

  series_start <- lubridate::ymd_hms(sleeptimes$sleep.start[1], tz = tz)
  series_end <- lubridate::ymd_hms(sleeptimes$sleep.end[nrow(sleeptimes)], tz = tz)
  series_end <- lubridate::floor_date(series_end, "day") + lubridate::hours(23) + lubridate::minutes(59)

  adaptation_df <- parse_sleeptimes(
    sleeptimes = sleeptimes,
    series.start = series_start,
    series.end = series_end,
    sleep.start.col = "sleep.start",
    sleep.end.col = "sleep.end",
    sleep.id.col = "sleep.id",
    roundvalue = round_minutes
  )

  adaptation_run <- FIPS_simulate(
    FIPS_df = adaptation_df,
    modeltype = "mccauley2013",
    pvec = pvec,
    pvec_C = pvec_C,
    dynamic_C = dynamic_C,
    dynamic_vars=dynamic_vars
  )

  last_day_start <- lubridate::floor_date(adaptation_run$datetime[nrow(adaptation_run)], "day")
  eq <- dplyr::filter(adaptation_run, datetime >= last_day_start)

  target_tod <- series_start_time
  idx <- which.min(abs(eq$time - target_tod))
  if (length(idx) == 0L || is.na(idx)) stop("Failed to locate target time-of-day in equilibrium day.")
  if (dynamic_C & !dynamic_vars$approx_N) {
    inits <- c(
      p0 = unname(eq$p[idx]),
      u0 = unname(eq$u[idx]),
      k0 = unname(eq$k[idx]),
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
      xc0 = unname(eq$xc[idx]),
      x0 = unname(eq$x[idx]),
      CBT_min0 = unname(eq$CBT_min[idx])
    )
  } else {
    inits <- c(
      p0 = unname(eq$p[idx]),
      u0 = unname(eq$u[idx]),
      k0 = unname(eq$k[idx])
    )
  }

  if (!is.null(save_equilibrium_path)) {
    saveRDS(eq, file = save_equilibrium_path)
  }

  if (return_equilibrium && !return_adaptation) {
    return(list(inits = inits, equilibrium = eq))
  }
  if (!return_equilibrium && return_adaptation) {
    return(list(inits = inits, adaptation = adaptation_run))
  }
  if (return_equilibrium && return_adaptation) {
    return(list(inits = inits, equilibrium = eq, adaptation = adaptation_run))
  }
  inits
}

#' Convenience helper to inject steady-state inits into a parameter vector
#'
#' @inheritParams mccauley2013_make_inits
#' @return Updated `pvec` with fresh initial conditions.
#' @export
mccauley2013_init_pvec <- function(
    sim_start,
    sleep_hrs = 8,
    wake_time = "07:30:00",
    series_start_time = 23.5,
    ndays = 30,
    tz = NULL,
    round_minutes = 5,
    pvec = mccauley2013_pvec,
    pvec_C = FJK_pvec,
    dynamic_C = FALSE,
    dynamic_vars=list(approx_N=FALSE,model="FJK")
) {
  inits <- mccauley2013_make_inits(
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
  pvec[c("p0", "u0", "k0")] <- inits[c("p0", "u0", "k0")]
  if (dynamic_C & !dynamic_vars$approx_N) {
    pvec_C[c("xc0","x0","L0","CBT_min0")] = inits[c("xc0","x0","L0","CBT_min0")]
  } else if (dynamic_C & dynamic_vars$approx_N) {
    pvec_C[c("xc0","x0","CBT_min0")] = inits[c("xc0","x0","CBT_min0")]
  }
  return(list(pvec = pvec, pvec_C = pvec_C))
}
