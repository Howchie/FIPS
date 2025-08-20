## Generate 7 days of sleep deprived (4hr per day) + 7 days recovery (8hrs per day) then a period of several days without sleep
library(FIPS)
library(dplyr)
library(lubridate)
library(tidyr)
library(ggplot2)
# Anchor (first wake) and series bounds
wake_datetime <- lubridate::ymd_hms('2025-05-03 05:00:00', tz = "Australia/Perth")

# Define any number of cycles here (add rows as needed).
cycles <- tibble::tribble(
  ~n_days, ~sleep_hrs, ~wake_time,
  7,      4,          "03:00:00",
  7,      8,          "07:00:00",
  3,       5,          "18:00:00",
  4,       10,         "08:30:00"
)

sleeptimes <- build_sleeptimes(anchor_wake = wake_datetime, cycles = cycles, tz = "Australia/Perth")

# Feed into FIPS parser
simulation_df <- parse_sleeptimes(
  sleeptimes = sleeptimes,
  series.start = lubridate::ymd_hms('2025-05-02 23:00:00', tz = "Australia/Perth"),
  series.end   = lubridate::ymd_hms('2025-05-23 23:00:00', tz = "Australia/Perth"),
  sleep.start.col = "sleep.start",
  sleep.end.col   = "sleep.end",
  sleep.id.col    = "sleep.id",
  roundvalue = 5
) %>%
  dplyr::filter(datetime >= sleeptimes$sleep.start[1])


## Simulations
# Simulate from Unified Model
simulation_unified <- FIPS_simulate(
  FIPS_df = simulation_df ,
  modeltype = "unified",
  pvec = unified_make_pvec(
    U0 = 24.12,
    L0 = 0,
    S0 = 0,
    phi = 2.02,
    kappa = 4.13,
    tau_s = 1,
    tau_w = 40,
    tau_la = 4.06*24,
    sigma = 1,
    wc = 1.14,
    wd = -0.46
  )
)

# Simulate from Three Process Model
simulation_tpm <- FIPS_simulate(
  FIPS_df = simulation_df ,
  modeltype = "TPM",
  pvec = TPM_make_pvec(
    la = 2.4,
    ha = 14.3,
    d = -0.0353,
    g = log((14.3 - 14) / (14.3 - 7.96)) / 8,
    bl = 12.2,
    Cm = 0,
    Ca = 2.5,
    p = 16.8,
    Um = -0.5,
    Ua = 0.5,
    Wc = -5.72,
    Wd = -1.51,
    S0 = 7.96,
    KSS_intercept = 10.6,
    KSS_beta = -0.6
  )
)

# Compute initial states for ODE (can be either at start of sleep or start of rest, assumes a steady state -- might not be perfect)
ODE_init = mccauley_initial_values(mccauley2013_make_pvec(), W=16, T=24, t0=8.5)

## Simulate from the McCauley (2013) model (uses the later 2021 mappings with positive parameters, also essentially nested within the 2024 model)
simulation_mccauley2013 <- FIPS_simulate(
  FIPS_df = simulation_df ,
  modeltype = "mccauley2013",
  pvec = mccauley2013_make_pvec(
    alpha_w = 0.028,   
    alpha_s = 0.260, 
    beta_w =  0.260,  
    beta_s =  0.260,   
    eta_w = 0.0074, 
    Wc =  20.2,
    Tp = 24,
    mu_w =  .33,
    mu_s = -1.50,
    phi = 21.02,     
    lambda_w = 0.49,
    lambda_s= -0.49,
    xi = 1.09,
    p0=ODE_init$p_sleep_start,
    u0=ODE_init$u_sleep_start,
    k0=ODE_init$k_sleep_start,
    pf=0
  )
)

# Simulate from the McCauley (2024) model that includes sleep inertia and bidirectional feedback
simulation_mccauley2024 <- FIPS_simulate(
  FIPS_df = simulation_df ,
  modeltype = "mccauley2024",
  pvec = mccauley2024_make_pvec(
    alpha_w = 0.028,   
    alpha_s = 0.260,   
    beta_w  = 0.26,   
    beta_s  = 0.26,   
    eta_w   = 0.0126,  
    Wc      =  20.2, 
    Tp      = 24,
    mu_w    = .466,
    mu_s    = -1.50,
    phi     = 21.2,     
    lambda_w =  0.49,    
    lambda_s =  0.49,
    xi_u      =  1.09,
    xi_k      =  1.09,
    xi_h      =  1.09,
    p0=ODE_init$p_sleep_start,
    u0=ODE_init$u_sleep_start,
    k0=ODE_init$k_sleep_start,
    pf=0,h0=0,
    zeta_w  = 1.31,    
    zeta_s  = 1.31,
    v_w=1.37,
    v_s=1.37,
    gamma=.71
  )
)

# Simulate from the McCauley (2024) model with KSS parameterization
simulation_mccauley2024_KSS <- FIPS_simulate(
  FIPS_df = simulation_df ,
  modeltype = "mccauley2024",
  pvec = mccauley2024_make_pvec(
    alpha_w = 0.022,   
    alpha_s = 0.037,   
    beta_w  = 0.26,   
    beta_s  = 0.26,   
    eta_w   = 0.0126,  
    Wc      =  22.02, 
    Tp      = 24,
    mu_w    = .82,
    mu_s    = -1.50,
    phi     = 21.2,     
    lambda_w =  0.49,    
    lambda_s =  0.49,
    xi_u      =  0.51,
    xi_k      =  0.51,
    xi_h      =  0.51,
    p0=ODE_init$p_sleep_start,
    u0=ODE_init$u_sleep_start,
    k0=ODE_init$k_sleep_start,
    pf=1,h0=0,
    zeta_w  = 1.31,    
    zeta_s  = 1.31,
    v_w=1.37,
    v_s=1.37,
    gamma=.71
  )
)

# Showcase new overlayed plot
plots <- FIPS_plot_overlay(
  list(
    unified=simulation_unified,
    mccauley2024 = simulation_mccauley2024
  ),
  plot_stat = "fatigue"
)
plots

# Haven't figured out how (or if) to show different stats on same figure
plots <- FIPS_plot_overlay(
  list(
    TPM=simulation_tpm,
    mccauley2024 = simulation_mccauley2024_KSS
  ),
  plot_stat = "KSS"
)
plots

## Demo shift-tagging feature to e.g. extract mean fatigue per shift across a time period
# anchor for the first day we want to generate shifts for
anchor <- as_date(min(simulation_mccauley2024$datetime))  # or a fixed Date
tz <- "Australia/Perth"

shift_schedule <- tibble::tribble(
  ~n_days, ~start_time, ~end_time, ~label, ~start_offset,
  13L,      "09:00:00",  "17:00:00", "Day", 0,
  3,        "20:00:00",  "04:00:00", "Night", 13,
  4L,      "09:00:00",  "17:00:00", "Day", 17
)

shifts <- build_shifts(anchor, shift_schedule, tz)

df_tagged <- tag_shifts(simulation_mccauley2024, shifts)

# Now you can summarise by shift
summary_by_shift <- df_tagged %>%
  group_by(shift_id, label) %>%
  summarise(
    n = n(),
    mean_fatigue = mean(fatigue, na.rm = TRUE),
    .groups = "drop"
  ) %>% filter(!is.na(label))
plot(summary_by_shift$mean_fatigue)

library(rlang)

# df:       time series with a POSIXct 'datetime' column + your metric
# metric:   bare column name in df (e.g., fatigue)
# shifts:   data frame with shift_start, shift_end, label (from build_shifts)
# series:   optional bare column in df to color lines/points (e.g., Label)
# shade_df: optional intervals to shade (cols: start, end, label). If NULL, uses shifts.
# x_mode:   "hours" for hours since start, or "datetime" for real time axis
plot_shift_metric <- function(df,
                              metric,
                              shifts,
                              series = NULL,
                              shade_df = NULL,
                              x_mode = c("hours", "datetime"),
                              shade_alpha = 0.12,
                              thresholds=NULL) {
  x_mode <- match.arg(x_mode)
  metric_quo <- rlang::enquo(metric)
  series_quo <- rlang::enquo(series)
  
  stopifnot(inherits(df$datetime, "POSIXct"))
  stopifnot(all(c("shift_start","shift_end") %in% names(shifts)))
  
  # Tag rows with their shift (dplyr non-equi join, right-open [start,end))
  df_tag <- df %>% 
    left_join(
      shifts %>% select(shift_id, label, shift_start, shift_end),
      by = join_by(datetime >= shift_start, datetime < shift_end)
    )
  
  # Summaries per shift (and per series if provided)
  group_vars <- c("shift_id", "label", "shift_start", "shift_end")
  if (!rlang::quo_is_null(series_quo)) {
    df_tag <- df_tag %>% mutate(.series = !!series_quo)
    group_vars <- c(group_vars, ".series")
  }
  
  shift_summ <- df_tag %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      y_mean = mean(!!metric_quo, na.rm = TRUE),
      y_min  = min(!!metric_quo, na.rm = TRUE),
      y_max  = max(!!metric_quo, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(shift_mid = shift_start + (shift_end - shift_start)/2)
  
  # X mapping (hours since start or datetime)
  t0 <- min(df$datetime, shifts$shift_start, na.rm = TRUE)
  
  map_x <- function(t) {
    if (x_mode == "hours") as.numeric(difftime(t, t0, units = "hours")) else t
  }
  
  df_plot <- df %>%
    mutate(.x = map_x(datetime),
           .y = !!metric_quo)
  
  summ_plot <- shift_summ %>%
    mutate(.x_mid = map_x(shift_mid))
  
  # Shading intervals: use provided shade_df if given, otherwise use shifts
  if (is.null(shade_df)) {
    shade_df <- shifts %>%
      transmute(start = shift_start, end = shift_end, label = label %||% "Shift")
  } else {
    stopifnot(all(c("start","end") %in% names(shade_df)))
  }
  shade_plot <- shade_df %>%
    transmute(
      xmin = map_x(start),
      xmax = map_x(end),
      label = label
    )
  
  # Day break lines
  x_min_dt <- min(df$datetime, na.rm = TRUE)
  x_max_dt <- max(df$datetime, na.rm = TRUE)
  day_breaks <- seq(floor_date(x_min_dt, "day"),
                    ceiling_date(x_max_dt, "day"),
                    by = "1 day")
  vline_x <- map_x(day_breaks)
  # Build plot ---------------------------------------------------------------
  p <- ggplot()
  
  # background shading
  p <- p + geom_rect(
    data = shade_plot,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = "grey50", alpha = shade_alpha
  )
  
  # day grid
  p <- p + geom_vline(xintercept = vline_x, linetype = 3, linewidth = 0.3) +
    geom_hline(yintercept = thresholds, linetype = 2, linewidth = 0.3)
  
  
  # shift summaries: error bars (min–max) + mean point
  if (!rlang::quo_is_null(series_quo)) {
    p <- p +
      geom_errorbar(data = summ_plot, aes(x = .x_mid, ymin = y_min, ymax = y_max, color = .series), width = 0) +
      geom_point    (data = summ_plot, aes(x = .x_mid, y = y_mean, color = .series), size = 2)
  } else {
    p <- p +
      geom_errorbar(data = summ_plot, aes(x = .x_mid, ymin = y_min, ymax = y_max), width = 0) +
      geom_point    (data = summ_plot, aes(x = .x_mid, y = y_mean), size = 2)
  }
  
  # labels & theme
  xlab_txt <- if (x_mode == "hours") "Simulation Hours" else "Time"
  p +
    labs(x = xlab_txt, y = rlang::as_name(metric_quo)) +
    theme_classic(base_size = 12)
}

# Assume you already have:
#   simulation_df: with columns datetime (POSIXct), fatigue, Label (model name), etc.
#   shifts:        from build_shifts(...)

# Shade sleep instead of shifts? Build a generic interval df:
sleep_shade <- sleeptimes %>%
  transmute(start = sleep.start, end = sleep.end, label = "Sleep")

plot_shift_metric(
  df        = simulation_mccauley2024,
  metric    = fatigue,
  shifts    = shifts,          # used for per-shift summaries
  shade_df  = sleep_shade,     # background shading (optional)
  x_mode    = "datetime",
  thresholds=c(11,20)# or "datetime"
)

