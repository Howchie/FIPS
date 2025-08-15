#' FIPS Time Series Plot
#'
#' @param dats A FIPS_simulation object (i.e., FIPS_df with simulation results)
#' @param from The starting datetime to be plotted
#' @param to The ending datetime to be plotted
#' @param plot_stat Which variables to plot
#' @param fatigue_CIs A logical indicating whether uncertainty intervals on fatigue should be plotted
#' @param observed A data frame with any observed sleepiness ratings or objective indicators to plot against predictions
#' @param observed_y The name of the observed sleepiness ratings in the observed data frame
#' @param ... additional arguments passed to ggplot2 theme_classic
#' @importFrom rlang .data
#' @return A ggplot2 object displaying fatigue and other requested processes over time
#' @md
#' @export
FIPS_plot <- function(dats,
                      from = NULL,
                      to = NULL,
                      plot_stat = NULL,
                      fatigue_CIs = FALSE,
                      observed = NULL,
                      observed_y = NULL,
                      ...) {


  if(FIPS_Simulation_lost_attributes(dats)) {
    stop("Your FIPS_Simulation object has lost attributes (have you wrangled the dataframe with dplyr?).
          No plot method availble. Please manually save attributes if plotting essential.")
  }


  if(!is_FIPS_simulation(dats)) {
    stop("Requires a FIPS_df which has had model simulation run on it")
  }

  modeltype = get_FIPS_modeltype(dats)
  if(! modeltype %in% c("TPM", "unified","mccauley")) {
    warning("You supplied a modeltype argument that doesn't match the model type specified in your
             FIPS_simulation FIPS_df. Defaulting to using one specified in the FIPS_df.")
  }

  # Figure out appropriate default plot_stat and observed_y based on modeltype
  if(is.null(plot_stat)) plot_stat <- get_FIPS_pred_stat(dats)
  # Make observation variable plot_stat unless otherwise specified
  if(is.null(observed_y)) observed_y <- plot_stat

  # Filter based on selected datetimes from and to
  if (!is.null(from)) {
    dats <- dats %>% dplyr::filter(.data$datetime > from)
    if (!is.null(observed))
      observed <- observed %>% dplyr::filter(.data$datetime > from)
    }

  if (!is.null(to)) {
    dats <- dats %>% dplyr::filter(.data$datetime < to)
    if (!is.null(observed))
      observed <- observed %>% dplyr::filter(.data$datetime < to)
  }

  if(!any((get_FIPS_pred_stat(dats) %in% plot_stat)) & fatigue_CIs == TRUE){
    warning("Will not plot fatigue CIs without a predicted model value (alertness/fatigue)")
    fatigue_CIs = FALSE
  }

  # Get start and end of sleeptimes for plotting sleep as rectangles
  sim_results <- dats  %>%
    dplyr::group_by(.data$sleep.id) %>% dplyr::mutate(
      sleepstart = ifelse(is.na(.data$sleep.id), NA, min(.data$sim_hours)),
      sleepend = ifelse(is.na(.data$sleep.id), NA, max(.data$sim_hours))) %>%
    #Also get end of each day for dashed lines indicating day end
    dplyr::group_by(day) %>%
    dplyr::mutate(eod = .data$sim_hours[which(time == max(time))] + 24 - max(time))

  # Filter out any end of days after specified date range
  sim_results$eod[sim_results$eod > max(sim_results$sim_hours)] <- NA

  plot_out <- ggplot2::ggplot(sim_results, aes(x = .data$sim_hours)) +
    geom_rect(aes(xmin = .data$sleepstart, xmax = .data$sleepend,
                  ymin = -Inf, ymax = Inf, fill = 'Sleep'), alpha = 0.1, na.rm = T) +
    scale_fill_manual('Sleep', name = "", values = 'grey80', guide = guide_legend(override.aes = list(alpha = 1))) +
    geom_vline(aes(xintercept = .data$eod), linetype = 2, na.rm = T) +
    theme_classic(...) +
    xlab ("Simulation Hours") +
    ylab("")


  long_data <- tidyr::pivot_longer(sim_results, !!plot_stat, names_to = "stat")

  # Change factor order to put fatigue first
  if("fatigue" %in% plot_stat){
    fac_levels <- c("fatigue", unique(long_data$stat)[!unique(long_data$stat) == "fatigue"])
    long_data$stat <- factor(long_data$stat, levels = fac_levels)
  }

  plot_out <- plot_out +
    geom_line(data = long_data, aes(y = .data$value, color = stat), size = 1) +
    labs(colour = "Predicted Value")

  if (fatigue_CIs) {
    # Figure out which fill is appropriate
    correct_fill <- which(levels(long_data$stat) == "fatigue") + 1
    plot_out <- plot_out +
      geom_ribbon(aes(ymin = .data$fatigue_lower, ymax = .data$fatigue_upper), alpha = 0.2, fill = correct_fill)
    }

  if (!is.null(observed)) {
    plot_out <-
      plot_out + geom_point(data = observed, aes(y = !!as.name(observed_y)))
  }

  return(plot_out)

}


#' plot.FIPS_simulation
#'
#' S3 plot method for FIPS_simulation
#'
#' @param x A valid .FIPS_simulation series that has been simulated
#' @param from The starting datetime to be plotted
#' @param to The ending datetime to be plotted
#' @param plot_stat Which variables to plot
#' @param fatigue_CIs A logical indicating whether uncertainty intervals on fatigue should be plotted
#' @param observed A data frame with any observed sleepiness ratings or objective indicators to plot against predictions
#' @param observed_y The name of the observed sleepiness ratings in the observed data frame
#' @param ... additional arguments passed to ggplot2 theme_classic
#'
#' @export
plot.FIPS_simulation <- function(
  x,
  from = NULL,
  to = NULL,
  plot_stat = NULL,
  fatigue_CIs = FALSE,
  observed = NULL,
  observed_y = NULL,
  ...) {
  FIPS_plot(
    dats = x,
    from = from,
    to = to,
    plot_stat = plot_stat,
    fatigue_CIs = fatigue_CIs,
    observed = observed,
    observed_y = observed_y,
    ...
  )
}

#' @export
#' 
FIPS_plot_overlay <- function(dats_list,
                              from = NULL,
                              to = NULL,
                              plot_stat = NULL,
                              fatigue_CIs = FALSE,
                              labels = NULL,
                              facet = NULL,
                              show_sleep = TRUE,
                              ...) {
  
  # Helper to check one FIPS object
  .check_fips <- function(d) {
    if (FIPS_Simulation_lost_attributes(d)) {
      stop("A FIPS_Simulation object lost attributes (was it wrangled?). Save attributes before plotting.")
    }
    if (!is_FIPS_simulation(d)) {
      stop("All elements must be FIPS_df objects with simulation results.")
    }
  }
  
  # If they accidentally pass a single object, hand off to original
  if (!is.list(dats_list)) {
    return(FIPS_plot(dats_list, from = from, to = to, plot_stat = plot_stat,
                     fatigue_CIs = fatigue_CIs, ...))
  }
  
  if (length(dats_list) == 0) stop("Empty list supplied.")
  lapply(dats_list, .check_fips)
  
  # Labels: prefer names(dats_list), else user-supplied, else sim1/sim2...
  if (is.null(labels)) labels <- names(dats_list)
  if (is.null(labels) || any(labels == "")) {
    labels <- paste0("sim", seq_along(dats_list))
  }
  
  # Default plot_stat: use the first object’s default if not provided
  # TODO make this per-list item
  if (is.null(plot_stat)) plot_stat <- get_FIPS_pred_stat(dats_list[[1]])
  plot_stat <- as.character(plot_stat)
  
  # Build a base copy (for sleep rects / eod)
  base_d <- dats_list[[1]]
  if (!is.null(from)) base_d <- dplyr::filter(base_d, .data$datetime > from)
  if (!is.null(to))   base_d <- dplyr::filter(base_d, .data$datetime < to)
  
  base_results <- base_d %>%
    dplyr::group_by(.data$sleep.id) %>%
    dplyr::mutate(
      sleepstart = ifelse(is.na(.data$sleep.id), NA, min(.data$sim_hours)),
      sleepend   = ifelse(is.na(.data$sleep.id), NA, max(.data$sim_hours))
    ) %>%
    dplyr::group_by(.data$day) %>%
    dplyr::mutate(eod = .data$sim_hours[which(time == max(time))] + 24 - max(time)) %>%
    dplyr::ungroup()
  
  base_results$eod[base_results$eod > max(base_results$sim_hours, na.rm = TRUE)] <- NA
  
  # Stack all datasets with a label
  make_long <- function(d, lab) {
    if (!is.null(from)) d <- dplyr::filter(d, .data$datetime > from)
    if (!is.null(to))   d <- dplyr::filter(d, .data$datetime < to)
    
    sim_results <- d %>%
      dplyr::group_by(.data$sleep.id) %>%
      dplyr::mutate(
        sleepstart = ifelse(is.na(.data$sleep.id), NA, min(.data$sim_hours)),
        sleepend   = ifelse(is.na(.data$sleep.id), NA, max(.data$sim_hours))
      ) %>%
      dplyr::group_by(.data$day) %>%
      dplyr::mutate(eod = .data$sim_hours[which(time == max(time))] + 24 - max(time)) %>%
      dplyr::ungroup()
    
    long <- tidyr::pivot_longer(
      sim_results,
      cols = dplyr::all_of(plot_stat),
      names_to = "stat",
      values_to = "value"
    )
    long$label <- lab
    long
  }
  
  long_all <- purrr::map2_dfr(dats_list, labels, make_long)
  
  # Faceting default: facet if >1 stat requested
  if (is.null(facet)) facet <- length(unique(long_all$stat)) > 1
  
  # Put fatigue first if present
  if ("fatigue" %in% long_all$stat) {
    fac_levels <- c("fatigue", setdiff(unique(long_all$stat), "fatigue"))
    long_all$stat <- factor(long_all$stat, levels = fac_levels)
  }
  
  # Build plot
  p <- ggplot2::ggplot()
  
  if (show_sleep) {
    p <- p +
      ggplot2::geom_rect(
        data = base_results,
        ggplot2::aes(xmin = .data$sleepstart, xmax = .data$sleepend,
                     ymin = -Inf, ymax = Inf, fill = "Sleep"),
        alpha = 0.1,
        inherit.aes = FALSE,
        na.rm = TRUE
      ) +
      ggplot2::scale_fill_manual(
        name = "", values = "grey80",
        guide = ggplot2::guide_legend(override.aes = list(alpha = 1))
      )
  }
  
  p <- p +
    ggplot2::geom_vline(
      data = base_results,
      ggplot2::aes(xintercept = .data$eod),
      linetype = 2,
      na.rm = TRUE
    )
  
  # Lines: color by label
  if (!facet && length(unique(long_all$stat)) == 1) {
    p <- p +
      ggplot2::geom_line(
        data = long_all,
        ggplot2::aes(x = .data$sim_hours, y = .data$value, color = .data$label),
        size = 1
      ) +
      ggplot2::labs(color = "Scenario")
  } else {
    p <- p +
      ggplot2::geom_line(
        data = long_all,
        ggplot2::aes(x = .data$sim_hours, y = .data$value,
                     color = .data$label),
        size = 0.9
      ) +
      ggplot2::facet_wrap(~stat, scales = "free_y", ncol = 1) +
      ggplot2::labs(color = "Label")
  }
  
  # Optional CIs for fatigue (neutral fill so color can stay on label)
  if (fatigue_CIs && "fatigue" %in% unique(long_all$stat)) {
    fat_dat <- long_all[long_all$stat == "fatigue", , drop = FALSE]
    has_bounds <- all(c("fatigue_lower", "fatigue_upper") %in% names(fat_dat))
    if (has_bounds) {
      p <- p +
        ggplot2::geom_ribbon(
          data = fat_dat,
          ggplot2::aes(x = .data$sim_hours,
                       ymin = .data$fatigue_lower, ymax = .data$fatigue_upper,
                       group = .data$label),
          alpha = 0.15,
          inherit.aes = FALSE
        )
    } else {
      warning("fatigue_CIs requested, but fatigue_lower/upper not found; skipping ribbons.")
    }
  }
  
  p +
    ggplot2::theme_classic(...) +
    ggplot2::xlab("Simulation Hours") +
    ggplot2::ylab("")
}

#' Build a collection of FIPS simulations for overlay plotting
#' @param ... either named FIPS_simulation objects, or a single named list of them
#' @return an object of class FIPS_collection
#' @export
#' 
FIPS_collect <- function(...) {
  x <- list(...)
  # allow FIPS_collect(list(name=sim, ...)) or FIPS_collect(name=sim, ...)
  if (length(x) == 1L && is.list(x[[1]]) && is.null(attr(x[[1]], "class"))) x <- x[[1]]
  
  if (!length(x)) stop("No simulations supplied.")
  if (is.null(names(x)) || any(names(x) == "")) {
    stop("Please provide names, e.g., FIPS_collect(baseline = sim1, nap = sim2).")
  }
  bad <- which(!vapply(x, is_FIPS_simulation, logical(1)))
  if (length(bad)) {
    stop("All elements must be FIPS_simulation objects. Problem at: ",
         paste0(names(x)[bad], collapse = ", "))
  }
  class(x) <- c("FIPS_collection", "list")
  x
}

#' @export
plot.FIPS_collection <- function(x,
                                 from = NULL,
                                 to = NULL,
                                 plot_stat = NULL,
                                 fatigue_CIs = FALSE,
                                 observed = NULL,
                                 observed_y = NULL,
                                 ...) {
  FIPS_plot_overlay(
    dats_list   = x,
    from        = from,
    to          = to,
    plot_stat   = plot_stat,
    fatigue_CIs = fatigue_CIs,
    # optional: pass-through observed once (shared points)
    labels      = names(x),
    ...
  )
}
