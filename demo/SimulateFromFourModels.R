## Generate 7 days of sleep deprived (4hr per day) + 7 days recovery (8hrs per day) then a period of several days without sleep
library(FIPS)
library(dplyr)
library(lubridate)
sleeptimes_common <- tibble::tibble(
  sleep.start = seq(
    from = lubridate::ymd_hms('2018-05-03 01:00:00', tz = "Australia/Perth"),
    to = lubridate::ymd_hms('2018-05-11 05:00:00', tz = "Australia/Perth"),
    by = '24 hours'),
  sleep.end = sleep.start + lubridate::dhours(8),
  sleep.id = rank(sleep.start)+7)

wake_datetime = lubridate::ymd_hms('2018-05-03 05:00:00', tz = "Australia/Perth")
ndays_deprived =7; n_days_recover=7
sleep_hrs_deprived = 4; sleep_hrs_recover=8

sleeptimes <- rbind(
  tibble(
    sleep.start = seq(
      from = wake_datetime - hours(sleep_hrs_deprived),
      to = wake_datetime - hours(sleep_hrs_deprived) + days(ndays_deprived - 1),
      by = '24 hours'),
    sleep.end = sleep.start + lubridate::dhours(sleep_hrs_deprived)
  ),
  tibble(
    sleep.start = seq(
      from =wake_datetime - hours(sleep_hrs_deprived) + days(ndays_deprived),
      to = wake_datetime - hours(sleep_hrs_recover)  + days(ndays_deprived) + days(n_days_recover),
      by = '24 hours'),
    sleep.end = sleep.start + lubridate::dhours(sleep_hrs_recover)
  )
) %>%
  mutate(sleep.id = rank(sleep.start))

simulation_df = parse_sleeptimes(
  sleeptimes = sleeptimes ,
  series.start = lubridate::ymd_hms('2018-05-02 21:00:00', tz = "Australia/Perth"),
  series.end = lubridate::ymd_hms('2018-05-19 23:00:00', tz = "Australia/Perth"),
  sleep.start.col = "sleep.start",
  sleep.end.col = "sleep.end",
  sleep.id.col = "sleep.id",
  roundvalue = 5
)

simulation_df <- simulation_df %>%
  dplyr::filter(datetime>=sleeptimes$sleep.start[1])

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

# Showcase new overlayed plot
plots <- FIPS_plot_overlay(
  list(
    unified=simulation_unified,
    mccauley2013 = simulation_mccauley2013,
    mccauley2024 = simulation_mccauley2024
  ),
  plot_stat = "fatigue"
)

plots

# Haven't figured out how (or if) to show different stats on same figure
plot(simulation_tpm,plot_stat="KSS")

