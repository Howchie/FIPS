#' Parse Sleep Times to FIPS_df
#'
#' This function parses a standardised sleeptime dataframe into the full FIPS format, ready for simulation and modelling.
#' The sleeptime format requires a sleep.id column (vector), a series of sleep times, and a series of corresponding wake times.
#' This format is the simplest to work with for human-readable or human-generated dataframes. See [parse_sleepwake_sequence] for
#' binary input methods.
#'
#' It is crucial that that following conditions are met for all arguments:
#' * Ensure that all specified datetimes for all datetime arguments are in an identical timezone.
#' * Ensure that the minimum sleep start time is >= series.start
#' * Ensure that the maximum wake time (sleep end) is <= series.end
#' * Ensure that each sleep start is < the corresponding sleep.end
#'
#' @param sleeptimes A dataframe in the sleep time format (see help for more info)
#' @param series.start A POSIXct object indicating the start datetime of the simulation (i.e., pre-first sleep waking duration)
#' @param series.end  A POSIXct object indicating the end datetime of the simulation
#' @param roundvalue An whole numeric (integer) value to round the sleep times to in minutes (`default = 5 minutes`). Second precision not supported.
#' @param sleep.start.col [string] The column in the dataframe containing the sleep start times
#' @param sleep.end.col [string] The column name in the dataframe containing the sleep end times
#' @param sleep.id.col [string] A column name specifying the sleep id sequence (i.e., `1:n()`)
#'
#' @examples
#'
#'  my_sleeptimes = tibble::tribble(
#'    ~sleep.id,          ~sleep.start,            ~sleep.end,
#'    1L, "2018-05-21 01:00:00", "2018-05-21 07:00:00",
#'    2L, "2018-05-21 23:00:00", "2018-05-22 04:00:00",
#'    3L, "2018-05-23 01:00:00", "2018-05-23 09:00:00") %>%
#'    dplyr::mutate(
#'      sleep.start = lubridate::ymd_hms(sleep.start),
#'      sleep.end = lubridate::ymd_hms(sleep.end))
#'
#'  my_simstart = lubridate::ymd_hms('2018-05-20 22:00:00')
#'  my_simend   = lubridate::ymd_hms('2018-05-23 10:00:00')
#'
#'  my_FIPS_df = parse_sleeptimes(
#'    sleeptimes = my_sleeptimes,
#'    series.start = my_simstart,
#'    series.end = my_simend,
#'    sleep.start.col = "sleep.start",
#'    sleep.end.col = "sleep.end",
#'    sleep.id.col = "sleep.id",
#'    roundvalue = 5)
#'
#' @seealso
#' For binary input parsing see: [parse_sleepwake_sequence]
#'
#' @return FIPS_df
#'
#' @export
#' @md
#'
#' @importFrom rlang := .data
parse_sleeptimes <- function(sleeptimes, series.start, series.end,
                             roundvalue = 5, sleep.start.col, sleep.end.col, sleep.id.col) {
  
  # Assert all colnames specified are actually in sleeptimes
  valid_names = checkmate::check_names(
    names(sleeptimes), permutation.of = c(sleep.id.col, sleep.start.col, sleep.end.col), what = "colnames")
  
  # Let's give a very useful message for this one.
  if(valid_names != T) {
    stop("At least one of the column strings you have specified does not exist in the sleeptimes dataframe. ",
         valid_names)}
  
  # Assert that series.start <= min(sleep.start.col) & length 1 & is a datetime & same timezones
  checkmate::assert_posixct(series.start, upper = min(sleeptimes[[sleep.start.col]]),
                            len = 1, .var.name = "series start datetime")
  # Assert that simulation end time >= max(sleep.end.col) & length 1 & is a datetime & same timezones
  checkmate::assert_posixct(series.end, lower = max(sleeptimes[[sleep.end.col]]),
                            len = 1, .var.name = "series end datetime")
  # Assert all sleep ends are less than simulation end times
  checkmate::assert_posixct(sleeptimes[[sleep.end.col]], upper = series.end, .var.name = "sleep.end datetimes")
  # Assert all sleep start times are less than sleep end times
  checkmate::assert_true(all(sleeptimes[[sleep.start.col]] < sleeptimes[[sleep.end.col]]))
  # Assert that ALL sleep end times <= series.end & all are a datetime & same timezone
  checkmate::assert_posixct(sleeptimes[[sleep.start.col]], lower = series.start, .var.name = "sleep.start datetimes")
  
  # roundvalue checks for whole number, as second precision not supported
  if(roundvalue < 1) stop("roundvalue must be a whole number > 1. FIPS does not support second-level precision")
  if(!(roundvalue%%1==0)) {
    # Print warning
    warning(paste("roundvalue must be a whole number > 1. FIPS does not support second-level precision.",
                  "User set round value of", roundvalue, "now set to", round(roundvalue)))
    # round userinput
    roundvalue = round(roundvalue)}
  
  
  
  # Now rename the user supplied sleeptime columns to "sleep.id", "sleep.start", and "sleep.end".
  sleeptimes = sleeptimes %>%
    dplyr::rename("sleep.id" := !!sym(sleep.id.col),
                  "sleep.start" := !!sym(sleep.start.col),
                  "sleep.end" := !!sym(sleep.end.col))
  
  # Round sleep and wake times to the desired epoch value
  rounded.sleeptimes <- sleeptimes %>%
    round_times("sleep.start", round_by = roundvalue) %>%
    round_times("sleep.end", round_by = roundvalue)
  
  # This makes the end of the sleep period occur 5 mins prior so that wake period starts at correct epoch
  rounded.sleeptimes <- rounded.sleeptimes %>%
    dplyr::mutate(sleep.end = .data$sleep.end - lubridate::minutes(roundvalue))
  
  # Assign minimum sleep start
  minimum.sleepstart = min(rounded.sleeptimes[["sleep.start"]])
  maximum.sleepend = max(rounded.sleeptimes[["sleep.end"]])
  
  # Now expand out the series of sleep wake times
  processed.sleeptimes <- expand_sleep_series(rounded.sleeptimes, expand_by = roundvalue)
  
  presleep.times <- NULL
  postwake.times <- NULL
  
  if(series.start < minimum.sleepstart) {
    presleep.times <- generate_presleep_times(series.start, minimum.sleepstart, roundvalue)
  }
  if(series.end > maximum.sleepend) {
    postwake.times <- generate_postwake_times(series.end, maximum.sleepend, roundvalue)
  }
  
  joined.times <- dplyr::bind_rows(presleep.times, processed.sleeptimes) %>%
    dplyr::bind_rows(postwake.times) %>%
    dplyr::mutate(wake_status_int = as.integer(.data$wake_status)) %>%
    dplyr::mutate(change_point = change_points(.data$wake_status_int)) %>%
    dplyr::mutate(switch_direction = status_dir(.data$wake_status_int, .data$change_point)) %>%
    dplyr::mutate(status_duration = time_in_status(.data$wake_status, roundvalue)) %>%
    dplyr::mutate(time_awake = ifelse(.data$wake_status, .data$status_duration, NA)) %>%
    dplyr::mutate(total_prev = shifted_time_status(.data$wake_status, .data$change_point, roundvalue)) %>%
    generate_decimal_timeunits("datetime")
  
  return(as_FIPS_df(joined.times))
}
#' @export
#' @rdname parse_sleeptimes
sleeptimes_to_FIPSdf = parse_sleeptimes
#' Fill pre-observation wake times
#'
#' The first sleep is unlikely to also be the start of the mission simulation
#' Thus, this function fills the start of the tibble with the all times between
#' The mission start time and the first instance of sleep, intervaled by X minutes
#'
#' @param simulationstart start of simulation
#' @param firstsleep first sleep in the sleep dataframe
#' @param expand_by expand
#' @return returns expanded tibble containing sleep.id = NA (due to waking) and wake_status = T
#' @keywords internal
generate_presleep_times <- function(simulationstart, firstsleep, expand_by = 5) {
  if (simulationstart >= firstsleep)
    stop("[Developer] Simulation Start must before first sleep if using this function")
  emins = paste(expand_by, "mins")
  tibble::tibble(
    datetime = seq(simulationstart, firstsleep - lubridate::minutes(expand_by), by = emins),
    sleep.id = NA,
    wake_status = T
  )
}
#' Fill post-observation wake times
#'
#' The last wake moment is unlikely to also be the end of the series.
#' This function fills constructs a tibble with the all times between
#' the final wake episode and the end of the series, intervaled by `expand_by` minutes
#'
#' @param simulationend start of simulation
#' @param lastwake first sleep in the sleep dataframe
#' @param expand_by expand value
#'
#' @return returns expanded tibble containing sleep.id = NA (due to waking) and wake_status = T
#' @keywords internal
#' @importFrom tibble tibble
generate_postwake_times <- function(simulationend, lastwake, expand_by = 5) {
  if (simulationend <= lastwake)
    stop("[Developer] Simulation end must after last sleep if using this function")
  emins = paste(expand_by, "mins")
  tibble::tibble(
    datetime = seq(lastwake + lubridate::minutes(expand_by), simulationend, by = emins),
    sleep.id = NA,
    wake_status = T
  )
}
#' Default light levels and cutoffs (lux)
#'
#' Provides a simple set of default illuminance levels (lux) for day, evening,
#' and night, plus time-of-day cutoffs used by [add_light_lux()].
#'
#' @param day Daytime illuminance in lux (default 500).
#' @param evening Evening illuminance in lux (default 120).
#' @param night Nighttime illuminance in lux (default 5).
#' @param morning Morning illuminance in lux (default 60).
#' @param morning_start Hour of day that morning begins (0-24, default 4).
#' @param day_start Hour of day that daytime begins (0-24, default 6).
#' @param evening_start Hour of day that evening begins (0-24, default 19).
#' @param night_start Hour of day that night begins (0-24, default 21).
#' @param fallback Fallback illuminance applied when no window matches (default 150).
#'
#' @return A list with `levels`, `cutoffs`, `windows`, and `fallback`.
#'
#' @examples
#' lux_defaults <- lux_levels()
#' lux_defaults$levels
#' lux_defaults$cutoffs
#' lux_defaults$windows
#' lux_defaults$fallback
#'
#' @export
lux_levels <- function(labels = c("morning", "day", "evening", "night"),
                       levels = c(60,250,50,10),
                       start = c(4,7,18,21), stop = c(7,18,21,4),
                       fallback = 150) {
  windows <- data.frame(
    labels = labels,
    levels = levels,
    start = start,
    stop = stop,
    stringsAsFactors = FALSE
  )
  list(
    windows = windows,
    fallback = as.numeric(fallback)
  )
}
#' Add a `light_lux` column to a FIPS dataframe
#'
#' Adds illuminance (lux) to your FIPS data as a new `light_lux` column. Users can:
#'   1) supply a measured light series to be resampled onto the FIPS timeline,
#'   2) provide a custom function that returns lux values, or
#'   3) use default named levels based on time-of-day.
#'
#' Priority is `series` > `func` > `levels/cutoffs` > package defaults.
#'
#' The resulting `FIPS_df` gains attributes `light_source` and `light_meta`
#' describing provenance for reproducibility.
#'
#' @param FIPS_df A FIPS-formatted dataframe (tibble). Must contain either
#'   a `datetime` column (`POSIXct`) or a numeric `time_h` column (continuous hours).
#' @param series Optional dataframe of measured light with either a `datetime`
#'   (`POSIXct`) or `time_h` column and a numeric `lux` column. Values are mapped
#'   to the FIPS timeline using last-observation-carried-forward (step function).
#' @param func Optional function to generate lux. It will be called either as
#'   `fun(FIPS_df)` and must return a numeric vector of `nrow(FIPS_df)`, or as
#'   `fun(x)` where `x` is the chosen time key (`datetime` or `time_h`).
#' @param levels Optional named numeric vector like `c(day=250, evening=30, night=3)`
#'   to override defaults.
#' @param cutoffs Optional named numeric vector like
#'   `c(day_start=7, evening_start=18, night_start=22)` to override defaults.
#' @param level_windows Optional specification of start/stop hours (0-24) for each level.
#'   Supply either a named list (e.g., `list(day = c(6, 19))`) or a data frame with
#'   columns `level`, `start`, and `stop`. Overrides defaults and any `cutoffs`.
#' @param fallback Numeric lux value applied when no window matches (default `150`).
#' @param join Which time key to use when both exist. One of `c("datetime", "time_h", "time")`.
#'   When omitted, the first available column in `FIPS_df` is used.
#' @param overwrite If `FALSE` and `light_lux` already exists, returns `FIPS_df` unchanged.
#' @param tz Timezone to assume when working with `datetime` (default `"UTC"`).
#'
#' @return The input `FIPS_df` with a new numeric column `light_lux` and
#'   attributes `light_source` and `light_meta`.
#'
#' @examples
#' # Suppose `df` is your FIPS_df with `datetime`, `time_h`, and `sleep` columns
#' # 1) Measured series (irregular, multi-day) -> resampled onto df
#' set.seed(1)
#' meas <- tibble::tibble(
#'   datetime = seq(min(df$datetime), max(df$datetime), by = "30 min"),
#'   lux = pmax(0, 50 + 200 * sin(2*pi*(as.numeric(lubridate::hour(datetime)) +
#'                                  lubridate::minute(datetime)/60)/24) + rnorm(dplyr::n(), 0, 20))
#' )
#' df1 <- add_light_lux(df, series = meas)
#'
#' # 2) Custom function using entire df (e.g., sinusoidal daylight pattern)
#' my_fun <- function(df) 150 + 100 *
#'                   sin(2*pi*((df$time_h %% 24) - 10)/24)
#' df2 <- add_light_lux(df, fun = my_fun)
#'
#' # 3) Named levels + custom cutoffs (no series/function)
#' df3 <- add_light_lux(
#'   df,
#'   levels = c(day = 300, evening = 20, night = 2),
#'   cutoffs = c(day_start = 6, evening_start = 19, night_start = 22)
#' )
#'
#' # Using with deSolve inside your RHS:
#' # light_fun <- stats::approxfun(df1$time_h, df1$light_lux, method = "constant", rule = 2)
#' # derivs <- function(t, y, p) { L <- light_fun(t); ... }
#'
#' @export
add_light_lux <- function(FIPS_df,
                          series = NULL,
                          func = NULL,
                          level_windows = NULL,
                          fallback = 150,
                          join = c("datetime", "time"),
                          overwrite = TRUE,
                          tz = "UTC") {
  stopifnot(is.data.frame(FIPS_df))
  if (!overwrite && "light_lux" %in% names(FIPS_df)) return(FIPS_df)
  join_opts <- c("datetime", "time")
  join_supplied <- !missing(join)
  join <- match.arg(join, join_opts)
  if (join %in% names(FIPS_df)) {
    key <- join
  } else {
    available <- intersect(join_opts, names(FIPS_df))
    if (length(available) == 0) {
      stop("FIPS_df must contain 'datetime' or 'time'")
    }
    if (join_supplied) {
      stop(sprintf("Requested join column '%s' not found in FIPS_df", join))
    }
    key <- available[1]
  }
  to_num_time <- function(x) {
    as.numeric(lubridate::with_tz(x, tzone = tz))
  }
  step_resample <- function(x_src, y_src, x_tgt) {
    o <- order(x_src)
    x_src <- x_src[o]
    y_src <- y_src[o]
    idx <- findInterval(x_tgt, x_src, left.open = FALSE, rightmost.closed = TRUE)
    idx[idx < 1] <- 1L
    idx[idx > length(y_src)] <- length(y_src)
    y_src[idx]
  }

  normalize_windows <- function(win) {
    if (is.null(win)) return(NULL)
    if (is.data.frame(win)) {
      has_level <- "level" %in% names(win) || "levels" %in% names(win)
      has_label <- "label" %in% names(win) || "labels" %in% names(win)
      if (!all(c("start", "stop") %in% names(win)) || !has_level) {
        stop("level_windows data.frame must contain 'start', 'stop', and one of 'level'/'levels' columns")
      }
      if ("levels" %in% names(win) && !"level" %in% names(win)) win$level <- win$levels
      if ("labels" %in% names(win) && !"label" %in% names(win)) win$label <- win$labels
      if (!has_label) win$label <- NA_character_
    } else {
      stop("level_windows must be a data.frame")
    }

    win$start <- win$start %% 24
    win$stop <- win$stop %% 24
    win
  }
  assign_by_windows <- function(tod, windows_df, fallback_val) {
    n <- length(tod)
    lux <- rep(fallback_val, n)
    assigned <- rep(FALSE, n)
    for (i in seq_len(nrow(windows_df))) {
      lvl <- windows_df$label[i]
      lvl_val <- windows_df$level[i]
      if (is.null(lvl_val) || !is.finite(lvl_val)) next
      start <- as.numeric(windows_df$start[i]) %% 24
      stop <- as.numeric(windows_df$stop[i]) %% 24
      if (is.na(start) || is.na(stop)) next
      if (abs(start - stop) < .Machine$double.eps) {
        mask <- rep(TRUE, n)
      } else if (start < stop) {
        mask <- tod >= start & tod < stop
      } else {
        mask <- tod >= start | tod < stop
      }
      if (any(mask)) {
        lux[mask] <- lvl_val
        assigned[mask] <- TRUE
      }
    }
    list(lux = lux, unmatched = sum(!assigned))
  }
  provenance <- list()
  if (!is.null(series)) {
    has_dt <- "datetime" %in% names(series)
    has_time <- "time" %in% names(series)
    if (!(has_dt || has_time)) stop("series must contain 'datetime', 'time_h', or 'time'")
    if (!"lux" %in% names(series)) stop("series must contain a numeric 'lux' column")
    if (!is.numeric(series$lux)) stop("series$lux must be numeric")
    if (key == "datetime") {
      if (has_dt) {
        x_src <- to_num_time(series$datetime)
      } else {
        origin <- min(FIPS_df$datetime, na.rm = TRUE)
        time_col <- if (has_th) "time_h" else "time"
        x_src <- to_num_time(origin + as.numeric(series[[time_col]]) * 3600)
      }
      x_tgt <- to_num_time(FIPS_df$datetime)
    } else {
      time_col <- if (key %in% names(series)) {
        key
      } else if (has_time) {
        "time"
      } else {
        stop(sprintf("series must include a '%s' column when join = '%s'", key, key))
      }
      x_src <- as.numeric(series[[time_col]])
      x_tgt <- as.numeric(FIPS_df[[key]])
    }
    FIPS_df$light_lux <- as.numeric(step_resample(x_src, as.numeric(series$lux), x_tgt))
    provenance$light_source <- "user_series"
    provenance$light_meta <- list(key = key, series_cols = names(series))
  } else if (!is.null(func)) {
    vals <- try(func(FIPS_df), silent = TRUE)
    if (inherits(vals, "try-error") || length(vals) != nrow(FIPS_df)) {
      vals <- func(FIPS_df[[key]])
    }
    if (!is.numeric(vals) || length(vals) != nrow(FIPS_df)) stop("func must return numeric vector of length nrow(FIPS_df)")
    FIPS_df$light_lux <- as.numeric(vals)
    provenance$light_source <- "function"
    provenance$light_meta <- list(key = key)
  } else {
    defs <- lux_levels()
    windows <- normalize_windows(defs$windows)
    if (!is.null(level_windows)) {
      windows <- normalize_windows(level_windows)
    }
    if (is.null(windows) || nrow(windows) == 0) {
      stop("No level windows defined")
    }
    fallback_val <- fallback
    if (is.null(fallback_val) || length(fallback_val) == 0) {
      fallback_val <- defs$fallback
    }
    if (length(fallback_val) > 1) {
      fallback_val <- fallback_val[1]
    }
    fallback_val <- as.numeric(fallback_val)
    if (!is.finite(fallback_val)) {
      fallback_val <- if (!is.null(defs$fallback)) defs$fallback else 150
    }
    tod <- FIPS_df$time
    assign_res <- assign_by_windows(tod, windows, fallback_val)
    FIPS_df$light_lux <- as.numeric(assign_res$lux)
    provenance$light_source <- "defaults"
    provenance$light_meta <- list(
      windows = windows,
      fallback = fallback_val,
      key = key,
      unmatched = assign_res$unmatched
    )
  }
  if (any(!is.finite(FIPS_df$light_lux))) warning("Non-finite values produced in light_lux; check inputs")
  if (any(FIPS_df$light_lux < 0, na.rm = TRUE)) warning("Negative lux values detected; setting to 0")
  FIPS_df$light_lux[FIPS_df$light_lux < 0] <- 0
  attr(FIPS_df, "light_source") <- provenance$light_source
  attr(FIPS_df, "light_meta") <- provenance$light_meta
  FIPS_df
}
#' Round times by column
#'
#' @param .stdata The sleeptimes dataframe
#' @param colname the column required to be rounded
#' @param round_by Amount (in minutes) to round sleep times to
#'
#' @return The sleep dataframe with all sleep.start and sleep.end rounded to X minute interval
#' @importFrom dplyr mutate
#' @importFrom lubridate round_date
#'
#'
#' @export
round_times <- function(.stdata, colname, round_by = 5) {
  if(round_by < 5) warning("Rounding less than 5 will result in an excessively large dataframe for long series")
  
  colname = rlang::ensym(colname)
  .stdata %>%
    dplyr::mutate(!!colname := lubridate::round_date(!!colname, paste(round_by, "mins")))
}
#' Expand Sleep Times to full vector
#'
#' Turns the paired sleeptimes into a long single vectored datetime sequence
#'
#' @param .stdata A sleeptimes dataframe
#' @param expand_by Amount (in minutes) to expand sleep times by
#'
#' @importFrom rlang .data
#'
#' @return Sleeptimedataframe with single columns vector for datetime and wake status
#' @keywords internal
expand_sleep_series <- function(.stdata, expand_by = 5) {
  
  emins = paste(expand_by, "mins")
  
  .stdata %>%
    dplyr::group_by(.data$sleep.id) %>%
    tidyr::expand(datetime = seq(min(.data$sleep.start), max(.data$sleep.end), by = emins)) %>%
    dplyr::mutate(wake_status = F) %>%
    dplyr::ungroup() %>%
    tidyr::complete(datetime = seq(min(.data$datetime), max(.data$datetime), by = emins), fill = list(wake_status = T))
}
# Build N arbitrary cycles of sleep from an anchor wake datetime.
# `cycles` must have columns:
#   - n_days (integer): how many days in this cycle
#   - sleep_hrs (numeric): sleep duration (hours) for this cycle
#' @export
#' 
build_sleeptimes <- function(anchor_wake, cycles, tz = NULL) {
  stopifnot(inherits(anchor_wake, "POSIXct"))
  tz <- tz %||% lubridate::tz(anchor_wake)
  default_wake_time <- format(lubridate::with_tz(anchor_wake, tz), "%H:%M:%S")
  anchor_date <- as.Date(lubridate::with_tz(anchor_wake, tz),tz=tz)
  if (!"start_idx"%in%names(cycles)) {
    cycles <- cycles %>%
      mutate(start_idx = dplyr::lag(cumsum(n_days), default = 0L))
  }
  cycles = cycles %>%
    mutate(
      n_days    = as.integer(n_days),
      sleep_hrs = as.numeric(sleep_hrs),
      wake_time = if ("wake_time" %in% names(.)) wake_time else NA_character_,
      wake_time = tidyr::replace_na(wake_time, default_wake_time),
      cycle     = dplyr::row_number()
    ) 
  per_day <- purrr::map_dfr(seq_len(nrow(cycles)), function(i) {
    tibble(
      day_index = seq.int(from = cycles$start_idx[i],
                          length.out = cycles$n_days[i]),
      sleep_hrs = cycles$sleep_hrs[i],
      wake_time = cycles$wake_time[i]
    )
  })
  # Compose wake datetimes: (anchor_date + day_index) + cycle-specific time-of-day
  wake_date <- anchor_date + per_day$day_index
  tod <- lubridate::hms(per_day$wake_time)             # "HH:MM:SS" -> Period
  wake_dt   <- ymd_hms(paste(wake_date, per_day$wake_time), tz = tz)
  tibble(
    sleep.start = wake_dt - lubridate::dhours(per_day$sleep_hrs),
    sleep.end   = wake_dt
  ) %>%
    arrange(sleep.start) %>%
    mutate(sleep.id = row_number())
}
# schedule columns (per row):
#   n_days (int)            : number of calendar days this shift appears
#   start_time, end_time    : "HH:MM:SS"
#   label (chr, optional)   : name for the shift
#   start_offset (int, opt) : how many days after anchor to start (default 0)
#' @export
#' 
build_shifts <- function(anchor_date, schedule, tz) {
  # normalize anchor_date to Date in the target tz
  if (inherits(anchor_date, "POSIXt")) {
    anchor_date <- as_date(with_tz(anchor_date, tz))
  } else {
    anchor_date <- as_date(anchor_date)
  }
  
  if (!"label" %in% names(schedule))       schedule$label <- NA_character_
  if (!"start_offset" %in% names(schedule)) schedule$start_offset <- 0L
  
  schedule %>%
    mutate(
      n_days       = as.integer(n_days),
      start_offset = as.integer(start_offset),
      label        = coalesce(label, paste0("shift_", row_number())),
      start_time   = as.character(start_time),
      end_time     = as.character(end_time)
    ) %>%
    tidyr::uncount(n_days, .remove = FALSE, .id = "day_idx") %>%
    mutate(
      # independent day indices per row
      day_index   = start_offset + day_idx,
      date        = anchor_date + day_index,
      shift_start = ymd_hms(paste(date, start_time), tz = tz),
      raw_end     = ymd_hms(paste(date, end_time),   tz = tz),
      # roll over midnight if needed
      shift_end   = if_else(raw_end <= shift_start, raw_end + days(1), raw_end)
    ) %>%
    arrange(shift_start) %>%
    mutate(shift_id = row_number()) %>%
    select(shift_id, label, date, day_index, shift_start, shift_end)
}
# Tag a point-in-time dataframe (duplicates if shifts overlap by design)
#' @export
#' 
tag_shifts <- function(df, shifts) {
  stopifnot(inherits(df$datetime, "POSIXct"),
            inherits(shifts$shift_start, "POSIXct"),
            inherits(shifts$shift_end, "POSIXct"))
  
  # optional sanity check: timezones should match
  if (!identical(lubridate::tz(df$datetime[1]), lubridate::tz(shifts$shift_start[1]))) {
    warning("Time zones differ between df$datetime and shifts; align before joining.")
  }
  
  df %>%
    dplyr::left_join(
      shifts %>% dplyr::select(shift_id, label, shift_start, shift_end),
      by = dplyr::join_by(datetime >= shift_start, datetime < shift_end)
    )
}
