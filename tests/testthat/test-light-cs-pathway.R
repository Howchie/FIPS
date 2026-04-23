test_that("add_light_lux accepts both levels/labels and level/label window columns", {
  df <- data.frame(time = c(5, 8, 19, 22))

  defaults <- add_light_lux(df)
  expect_equal(defaults$light_lux, c(60, 250, 50, 10))

  wins <- data.frame(
    labels = c("day", "night"),
    levels = c(100, 10),
    start = c(0, 12),
    stop = c(12, 0)
  )
  custom <- add_light_lux(df, level_windows = wins)
  expect_equal(custom$light_lux, c(100, 100, 10, 10))
})

test_that("FJK default light pathway remains lux-driven", {
  df <- simulation_df[1:250, ]

  dv <- fips_default_dynamic_vars(model = "FJK")
  sim <- FIPS_simulate(
    FIPS_df = df,
    modeltype = "FJK",
    pvec = FJK_pvec,
    dynamic_vars = dv
  )

  expect_true("light_lux" %in% names(sim))
  expect_true("light_drive" %in% names(sim))
  expect_false("light_cs" %in% names(sim))
  expect_equal(sim$light_drive, sim$light_lux)
})

test_that("FJK CS pathway creates light_cs and maps to light_drive", {
  df <- simulation_df[1:250, ]

  dv <- fips_default_dynamic_vars(model = "FJK")
  dv$light$metric <- "cs"
  dv$light$spd$source <- "D65"

  sim <- FIPS_simulate(
    FIPS_df = df,
    modeltype = "FJK",
    pvec = FJK_make_pvec(I0 = 9500),
    dynamic_vars = dv
  )

  expect_true("light_cs" %in% names(sim))
  expect_true("light_drive" %in% names(sim))
  expect_equal(sim$light_drive, sim$light_cs)

  cs <- sim$light_cs[is.finite(sim$light_cs)]
  expect_true(length(cs) > 0)
  expect_true(all(cs >= 0))
  expect_true(all(cs <= 0.7 + 1e-8))
})

test_that("CS mode warns when force_I0 is FALSE and I0 differs from 0.7", {
  df <- simulation_df[1:250, ]

  dv <- fips_default_dynamic_vars(model = "FJK")
  dv$light$metric <- "cs"
  dv$light$cs$force_I0 <- FALSE
  dv$light$spd$source <- "D65"

  expect_warning(
    FIPS_simulate(
      FIPS_df = df,
      modeltype = "FJK",
      pvec = FJK_make_pvec(I0 = 9500),
      dynamic_vars = dv
    ),
    regexp = "I0 = 0.7"
  )
})
