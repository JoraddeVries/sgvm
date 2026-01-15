yin_beta_growth <- function(t, ym,  tm, te) {
    #calculate the growth rate at time t for parameters tm (midpoint), te (endpoint) and ym (maximum)
    cm = ym * ((2*te-tm) / (te * (te-tm))) * (tm/te)**(tm/(te-tm))
    r = t<te ? cm * ((te-t) / (te-tm)) * (t/tm)**(tm/(te-tm)) : 0
    return(r)
}

calc_woody_growth <- function(dt, par) {

  n <- nrow(dt)

  cell_cohorts <- data.table(
    cohort_id = integer(),
    birth_step = integer(),  # index in dt
    age = numeric(),         # age in growing degree days
    size = numeric(),        # cell size in micons
    wall = numeric()         # cell wall thickness in microns
  )

  next_cohort_id <- 1
  emergence_buffer <- 0

  for (i in seq_len(n)) {

    # Environmental drivers at this step
    env <- as.list(dt[i])

    # thermal time
    tbase = 5 #degrees C
    gdd = max(0, env.Temp-tbase)

    # ---- 1. Cohort emergence ----
    emergence_buffer <- emergence_buffer + gdd/100

    while (emergence_buffer >= 1) {
      cell_cohorts <- rbind(
        cell_cohorts,
        data.table(
          cohort_id = next_cohort_id,
          birth_step = i,
          age = 0,
          size = 0,
          wall = 0
        )
      )
      next_cohort_id <- next_cohort_id + 1
      emergence_buffer <- emergence_buffer - 1
    }

    # ---- 2. Update all existing cohorts in place ----
    if (nrow(cell_cohorts) > 0) {

      # Age update
      cell_cohorts[, age := age + gdd]

      # update states
      cell_cohorts[, size := size + yin_beta_growth(gdd, 1, 500, 1000)]
      cell_cohorts[, wall := wall + yin_beta_growth(gdd, 1, 800, 1500)]
    }
  }

  return(cell_cohorts)
}