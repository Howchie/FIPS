.fips_default_light_config <- function() {
  list(
    metric = "lux",
    input = "lux_levels",
    cs = list(enabled = FALSE, force_I0 = TRUE, cla_params = list()),
    spd = list(
      source = "D65",
      cct = NULL,
      basis = NULL,
      sens = NULL,
      ybar = NULL,
      wavelength = NULL,
      relative_spd = NULL,
      sources = NULL,
      sources_path = NULL
    ),
    fallback_lux = 150
  )
}

fips_default_dynamic_vars <- function(model = "FJK") {
  list(
    approx_N = FALSE,
    model = model,
    light_levels = lux_levels(),
    light = .fips_default_light_config()
  )
}

.fips_merge_list <- function(base, extra) {
  if (is.null(extra)) return(base)
  nms <- names(extra)
  if (is.null(nms)) return(base)
  for (nm in nms) {
    val <- extra[[nm]]
    if (is.list(base[[nm]]) && is.list(val)) {
      base[[nm]] <- .fips_merge_list(base[[nm]], val)
    } else {
      base[[nm]] <- val
    }
  }
  base
}

fips_normalize_dynamic_vars <- function(dynamic_vars = NULL, model = NULL) {
  if (is.null(dynamic_vars)) {
    return(fips_default_dynamic_vars(model = if (is.null(model)) "FJK" else model))
  }
  if (!is.list(dynamic_vars)) stop("dynamic_vars must be a list")

  use_model <- if (!is.null(dynamic_vars$model)) dynamic_vars$model else if (!is.null(model)) model else "FJK"
  out <- fips_default_dynamic_vars(model = use_model)
  out <- .fips_merge_list(out, dynamic_vars)

  if (!is.list(out$light_levels) || is.null(out$light_levels$windows)) {
    out$light_levels <- lux_levels()
  }
  if (is.null(out$light)) out$light <- .fips_default_light_config()
  out$light <- .fips_merge_list(.fips_default_light_config(), out$light)

  if (isTRUE(out$light$cs$enabled)) out$light$metric <- "cs"
  out$light$metric <- tolower(as.character(out$light$metric)[1])
  if (!(out$light$metric %in% c("lux", "cs"))) {
    stop("dynamic_vars$light$metric must be one of 'lux' or 'cs'")
  }

  out$light$input <- tolower(as.character(out$light$input)[1])
  out
}

.fips_light_levels_apply <- function(dat, dynamic_vars) {
  if ("light_lux" %in% names(dat)) return(dat)
  light_levels <- dynamic_vars$light_levels
  if (!"func" %in% names(light_levels)) light_levels$func <- NULL
  add_light_lux(
    dat,
    level_windows = light_levels$windows,
    fallback = light_levels$fallback,
    func = light_levels$func
  )
}

.fips_adjust_pvec_for_light_metric <- function(pvec, dynamic_vars, model_name = "dynamic circadian") {
  metric <- dynamic_vars$light$metric
  if (!identical(metric, "cs")) return(pvec)
  if (!"I0" %in% names(pvec)) return(pvec)

  if (isTRUE(dynamic_vars$light$cs$force_I0)) {
    pvec[["I0"]] <- 0.7
  } else if (!isTRUE(all.equal(as.numeric(pvec[["I0"]]), 0.7, tolerance = 1e-8))) {
    warning(sprintf("%s using CS typically expects I0 = 0.7; received I0 = %s", model_name, as.character(pvec[["I0"]])))
  }
  pvec
}

.fips_resolve_sources_path <- function(dynamic_vars) {
  explicit <- dynamic_vars$light$spd$sources_path
  if (!is.null(explicit) && nzchar(explicit) && file.exists(explicit)) {
    return(explicit)
  }

  ext <- tryCatch(system.file("extdata", "Sources.RData", package = "FIPS"), error = function(e) "")
  if (nzchar(ext) && file.exists(ext)) return(ext)

  NULL
}

.fips_load_spd_library <- function(dynamic_vars) {
  provided <- dynamic_vars$light$spd$sources
  if (is.list(provided) && length(provided) > 0) return(provided)

  path <- .fips_resolve_sources_path(dynamic_vars)
  if (!is.null(path)) return(readRDS(path))

  stop(
    "CS mode requires packaged spectral source data at 'inst/extdata/Sources.RData', ",
    "or user-supplied dynamic_vars$light$spd$sources / sources_path."
  )
}

.fips_trapz <- function(x, y) {
  if (length(x) != length(y) || length(x) < 2) return(0)
  o <- order(x)
  x <- as.numeric(x[o])
  y <- as.numeric(y[o])
  dx <- diff(x)
  sum(dx * (head(y, -1) + tail(y, -1)) / 2)
}

.fips_scale_max <- function(x) {
  m <- suppressWarnings(max(x, na.rm = TRUE))
  if (!is.finite(m) || m <= 0) return(rep(0, length(x)))
  as.numeric(x) / m
}

.fips_cie_daylight_xy <- function(cct) {
  if (!is.finite(cct) || cct <= 0) stop("CCT must be finite and > 0")
  if (cct < 4000 || cct > 25000) {
    warning("CCT outside CIE daylight recommended range (4000K-25000K)")
  }
  x <- if (cct <= 7000) {
    -4.6070e9 / cct^3 + 2.9678e6 / cct^2 + 0.09911e3 / cct + 0.244063
  } else {
    -2.0064e9 / cct^3 + 1.9018e6 / cct^2 + 0.24748e3 / cct + 0.237040
  }
  y <- -3.0 * x^2 + 2.87 * x - 0.275
  c(x = x, y = y)
}

.fips_cie_daylight_M1M2 <- function(x, y) {
  denom <- 0.0241 + 0.2562 * x - 0.7341 * y
  M1 <- (-1.3515 - 1.7703 * x + 5.9114 * y) / denom
  M2 <- (0.0300 - 31.4424 * x + 30.0717 * y) / denom
  c(M1 = M1, M2 = M2)
}

.fips_cie_daylight_spd <- function(cct, basis, wl_min = 380, wl_max = 730) {
  need <- c("wl", "S0", "S1", "S2")
  if (!all(need %in% names(basis))) {
    stop("CIE daylight basis must contain columns: wl, S0, S1, S2")
  }
  xy <- .fips_cie_daylight_xy(cct)
  M <- .fips_cie_daylight_M1M2(xy[["x"]], xy[["y"]])
  b <- basis[basis$wl >= wl_min & basis$wl <= wl_max, need, drop = FALSE]
  spd <- b$S0 + M[["M1"]] * b$S1 + M[["M2"]] * b$S2
  data.frame(Wavelength = b$wl, S = .fips_scale_max(spd))
}

.fips_rescale_spd_to_lux <- function(target_lux, relative_spd, ybar, K = 683) {
  if (!is.finite(target_lux) || target_lux < 0) return(rep(NA_real_, length(relative_spd)))
  denom <- K * sum(as.numeric(ybar) * as.numeric(relative_spd), na.rm = TRUE)
  if (!is.finite(denom) || denom <= 0) return(rep(0, length(relative_spd)))
  (target_lux / denom) * as.numeric(relative_spd)
}

.fips_align_to_wavelength <- function(values, wl_values, wl_target, name = "signal") {
  if (length(values) != length(wl_values)) {
    stop(sprintf("%s and wavelength vectors must have equal length", name))
  }
  o <- order(wl_values)
  wl_values <- as.numeric(wl_values[o])
  values <- as.numeric(values[o])
  keep <- !(duplicated(wl_values) | is.na(wl_values))
  wl_values <- wl_values[keep]
  values <- values[keep]
  if (length(wl_values) < 2) stop(sprintf("%s has insufficient wavelength support", name))
  stats::approx(x = wl_values, y = values, xout = wl_target, rule = 2)$y
}

.fips_cla2_calc <- function(wl, E_lambda, sens, cla_params = list()) {
  defaults <- list(
    cla_norm = 1548,
    by_yellow_weight = 0.2616,
    by_gain = 0.21,
    rod_to_melanopsin_gain = 2.30,
    rod_to_opponent_gain = 1.6,
    rod_saturation = 6.5215,
    cone_mix_melanopsin = 1,
    cone_mix_opponent = 0.16,
    cs_max = 0.7,
    cla_to_cs_scale = 355.7,
    cs_exp = 1.1026,
    temporal_factor = 1,
    field_factor = 1
  )
  cfg <- .fips_merge_list(defaults, cla_params)

  required <- c("Melanopsin", "Vlambda.mac", "Scone.mac", "Vprime")
  miss <- setdiff(required, names(sens))
  if (length(miss) > 0) stop("sens is missing: ", paste(miss, collapse = ", "))

  melanopsin_w <- .fips_trapz(wl, sens$Melanopsin * E_lambda)
  photopic_w <- .fips_trapz(wl, sens$Vlambda.mac * E_lambda)
  scone_w <- .fips_trapz(wl, sens$Scone.mac * E_lambda)
  rod_w <- .fips_trapz(wl, sens$Vprime * E_lambda)

  blue_minus_yellow <- scone_w - cfg$by_yellow_weight * photopic_w
  rod_adaptation <- 1 - exp(-rod_w / cfg$rod_saturation)

  rod_inhib_ipRGC <- cfg$rod_to_melanopsin_gain *
    (rod_w / (photopic_w + cfg$cone_mix_melanopsin * scone_w)) * rod_adaptation

  opponent_drive <- max(0, cfg$by_gain * blue_minus_yellow)

  rod_inhib_opponent <- cfg$rod_to_opponent_gain *
    (rod_w / (photopic_w + cfg$cone_mix_opponent * scone_w)) * rod_adaptation

  melanopsin_drive <- max(0, melanopsin_w)
  is_cool_state <- blue_minus_yellow >= 0

  cla_raw <- melanopsin_drive - rod_inhib_ipRGC +
    (is_cool_state * (opponent_drive - rod_inhib_opponent))

  cla <- max(0, cfg$cla_norm * cla_raw)
  cla_eff <- cfg$temporal_factor * cfg$field_factor * cla
  cs <- cfg$cs_max * (1 - (1 / (1 + (cla_eff / cfg$cla_to_cs_scale)^cfg$cs_exp)))

  list(CLA = cla, CS = cs)
}

.fips_get_spd_context <- function(dynamic_vars) {
  spd_cfg <- dynamic_vars$light$spd
  sources <- .fips_load_spd_library(dynamic_vars)

  wl <- spd_cfg$wavelength
  if (is.null(wl) && !is.null(sources[["Wavelength"]])) {
    wl <- as.numeric(sources[["Wavelength"]])
  }
  if (is.null(wl) || length(wl) < 2) stop("Unable to determine wavelength grid for CS computation")

  sens <- spd_cfg$sens
  if (is.null(sens)) sens <- sources[["Sens"]]
  if (is.null(sens)) stop("Missing spectral sensitivity table 'Sens' for CS computation")
  sens_wl_col <- intersect(c("wl", "Wavelength", "wavelength", "lambda", "Lambda"), names(sens))
  sens_wl <- if (length(sens_wl_col) > 0) sens[[sens_wl_col[1]]] else NULL
  if (is.null(sens_wl)) sens_wl <- wl

  ybar <- spd_cfg$ybar
  if (is.null(ybar) && is.list(sources[["ybar"]]) && !is.null(sources[["ybar"]]$values)) ybar <- sources[["ybar"]]$values
  if (is.null(ybar) && !is.list(sources[["ybar"]]) && !is.null(sources[["ybar"]])) ybar <- sources[["ybar"]]
  if (is.null(ybar) && !is.null(sources[["Sens"]][["Vlambda.mac"]])) ybar <- sources[["Sens"]][["Vlambda.mac"]]
  if (is.null(ybar)) stop("Missing ybar values for SPD scaling")

  input_mode <- dynamic_vars$light$input
  rel_spd <- spd_cfg$relative_spd
  src_name <- NA_character_

  if (is.null(rel_spd)) {
    if (identical(input_mode, "cct")) {
      cct <- spd_cfg$cct
      if (is.null(cct) || !is.finite(cct)) stop("light$input='cct' requires light$spd$cct")
      basis <- spd_cfg$basis
      if (is.null(basis)) basis <- sources[["CIE_Lambdas"]]
      if (is.null(basis)) stop("Missing CIE daylight basis (CIE_Lambdas)")
      day_spd <- .fips_cie_daylight_spd(cct, basis)
      rel_spd <- day_spd$S
      spd_wl <- day_spd$Wavelength
      src_name <- paste0("CCT_", format(round(cct, 2), scientific = FALSE))
    } else {
      src_name <- as.character(spd_cfg$source)[1]
      if (!nzchar(src_name)) src_name <- "D65"
      spd_entry <- sources[[src_name]]
      if (is.null(spd_entry)) stop(sprintf("SPD source '%s' not found in source library", src_name))
      rel_spd <- if (is.list(spd_entry) && !is.null(spd_entry$values)) spd_entry$values else spd_entry
      spd_wl <- if (!is.null(sources[["Wavelength"]])) as.numeric(sources[["Wavelength"]]) else wl
    }
  } else {
    spd_wl <- wl
    src_name <- "user_relative_spd"
  }

  rel_spd <- as.numeric(rel_spd)
  rel_spd[!is.finite(rel_spd)] <- 0
  if (any(rel_spd < 0, na.rm = TRUE)) {
    warning("Negative SPD values detected; clamping to 0")
    rel_spd[rel_spd < 0] <- 0
  }

  rel_spd <- .fips_align_to_wavelength(rel_spd, spd_wl, wl, name = "SPD")
  sens2 <- sens
  for (nm in c("Melanopsin", "Vlambda.mac", "Scone.mac", "Vprime")) {
    sens2[[nm]] <- .fips_align_to_wavelength(sens[[nm]], sens_wl, wl, name = nm)
  }
  ybar_wl <- if (!is.null(sources[["Wavelength"]])) sources[["Wavelength"]] else sens_wl
  ybar2 <- .fips_align_to_wavelength(ybar, ybar_wl, wl, name = "ybar")

  list(wl = wl, rel_spd = rel_spd, sens = sens2, ybar = ybar2, source_name = src_name)
}

.fips_compute_cs_from_lux <- function(lux_values, dynamic_vars) {
  ctx <- .fips_get_spd_context(dynamic_vars)
  out <- rep(NA_real_, length(lux_values))

  good <- is.finite(lux_values) & lux_values >= 0
  if (!any(good)) return(out)

  unique_lux <- sort(unique(as.numeric(lux_values[good])))
  cs_map <- setNames(numeric(length(unique_lux)), as.character(unique_lux))

  for (k in seq_along(unique_lux)) {
    this_lux <- unique_lux[k]
    E_lambda <- .fips_rescale_spd_to_lux(this_lux, ctx$rel_spd, ctx$ybar)
    if (all(E_lambda == 0)) {
      cs_map[k] <- 0
    } else {
      cs_map[k] <- .fips_cla2_calc(ctx$wl, E_lambda, ctx$sens, cla_params = dynamic_vars$light$cs$cla_params)$CS
    }
  }

  out[good] <- as.numeric(cs_map[as.character(lux_values[good])])
  attr(out, "spd_source") <- ctx$source_name
  out
}

.fips_build_light_drive <- function(dat, dynamic_vars) {
  use_model <- if (is.list(dynamic_vars) && !is.null(dynamic_vars$model)) dynamic_vars$model else "FJK"
  dynamic_vars <- fips_normalize_dynamic_vars(dynamic_vars, model = use_model)
  dat <- .fips_light_levels_apply(dat, dynamic_vars)

  metric <- dynamic_vars$light$metric
  if (identical(metric, "lux")) {
    dat$light_drive <- dat$light_lux
    attr(dat, "light_metric") <- "lux"
    return(dat)
  }

  dat$light_cs <- .fips_compute_cs_from_lux(dat$light_lux, dynamic_vars)
  dat$light_drive <- dat$light_cs
  attr(dat, "light_metric") <- "cs"
  dat
}
