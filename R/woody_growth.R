yin_beta_growth <- function(t, ym, tm, te) {
    #calculate the growth rate at time t for parameters tm (midpoint), te (endpoint) and ym (maximum)
    cm = ym * ((2*te-tm) / (te * (te-tm))) * (tm/te) ^ (tm/(te-tm))
    r = pmax(0, cm * ((te - t) / (te - tm)) * ((t / tm) ^ (tm/(te-tm))))
    return(r)
}

calc_woody_growth <- function(dt) {

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
    gdd = max(0, env$Temp-tbase) * env$time_step/86400

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
      # update cohort id and the emergence buffer
      next_cohort_id <- next_cohort_id + 1
      emergence_buffer <- emergence_buffer - 1
    }

    # ---- 2. Update all existing cohorts in place ----
    if (nrow(cell_cohorts) > 0) {

      # Age update
      cell_cohorts[, age := age + gdd]

      # Growth limitations of temperature and water
      f_W = env$f_Tr
      f_C = pmin(1,pmax(0,env$NPP))

      # update states
      cell_cohorts[, size := size + f_W * yin_beta_growth(t=age, ym=1, tm=100, te=200) * gdd]
      cell_cohorts[, wall := wall + f_C * yin_beta_growth(t=age, ym=1, tm=200, te=400)  * gdd]
    }
  }

  return(cell_cohorts)
}