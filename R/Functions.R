set_environment <- function(dt, input_data, par, cloud_cover = 0.2) {
 
  # Interpolate input data
  input_daily <- interpolate_data(input_data, par()$INTER)
  
  # Join monthly climate variables
  dt <- merge(dt, input_daily, by = "doy", all.x = TRUE, sort = FALSE)
  
  # Row-wise calling of clear_sky
  cs <- mapply(
    function(lat, doy, tod) {
      clear_sky(lat = (lat/180*pi), DOY = doy, f = tod)
    },
    par()$latitude, dt$doy, dt$tod,
    SIMPLIFY = FALSE
  )
  
  # Extract clear-sky components
  dt[, c("Ig", "Idir", "Idif", "theta", "phi", "dayLength") := 
       .(sapply(cs, `[[`, "Ig"),
         sapply(cs, `[[`, "Idir"),
         sapply(cs, `[[`, "Idif"),
         sapply(cs, `[[`, "theta"),
         sapply(cs, `[[`, "phi"),
         sapply(cs, `[[`, "dayLength"))
  ]

  # set timestep length in seconds
  dt[, time_step := dayLength/par()$n_steps * 3600]
  
  # calculate the fraction incoming radiation. Unit conversions so irradiance is expressed in kJ/m2/day
  # interpolated, or monthly
  if(par()$INTER) {
    # 1. Monthly aggregated values
    monthly <- dt[, .(
      fIrr_month = mean(srad) /
        (sum(Ig * time_step / 1000, na.rm = TRUE) /
        30.5 / (par()$n_cohorts + 1))
    ), by = month]

    # 2. Mid-month doy for interpolation
    monthly[, doy_mid := c(15, 45, 74, 105, 135, 166,
                          196, 227, 258, 288, 319, 349)]

    # 3. Interpolate to daily (and therefore to all tod rows)
    dt[, fIrr := approx(
            x = monthly$doy_mid,
            y = monthly$fIrr_month,
            xout = doy,
            rule = 2
          )$y]
  } else {
    dt[, fIrr := mean(srad) /
          (sum(Ig * time_step / 1000, na.rm = TRUE) /
           30.5 / (par()$n_cohorts + 1)),
   by = month]
  }

  # Row-wise call to cloudy_sky
  cls <- mapply(
    function(Ig, fIrr, lat, doy, tod) {
      cloudy_sky(Ig = Ig * fIrr, lat = (lat/180*pi), DOY = doy, f = tod)
    },
    dt$Ig, dt$fIrr, par()$latitude, dt$doy, dt$tod,
    SIMPLIFY = FALSE
  )
  
  # Extract cloudy-sky components
  dt[, c("Ig", "Idir", "Idif", "theta", "phi", "dayLength") := 
       .(sapply(cls, `[[`, "Ig"),
         sapply(cls, `[[`, "Idir"),
         sapply(cls, `[[`, "Idif"),
         sapply(cls, `[[`, "theta"),
         sapply(cls, `[[`, "phi"),
         sapply(cls, `[[`, "dayLength"))
  ]
  
  # Add PAR conversion factors to the data.table
  dt[, cdir := waveband_conversion(Itype = "direct",  waveband = "PAR", mode = "flux")]
  dt[, cdif := waveband_conversion(Itype = "diffuse", waveband = "PAR", mode = "flux")]
  
  # atmospheric COS concentration in ppm
  dt[, Ca := par()$Ca]
  
  # air temperature
  dt[, Temp := tmin + (tmax - tmin) * sin((pi/2) * (tod / 0.6))]
  
  # Soil Water content
  dt[, Water := par()$Wmax * par()$Winit]
  
  # Snow
  dt[, Snow := par()$Sinit]

  return(dt)
}

get_data <- function(lat, lon, data) {
  
  vars <- c("prec","srad","tmin","tmax","vapr","lai")
  
  # initialise output
  dt <- data.table(
    month = 1:12
  )
  
  # loop over the different vars in the column, these are provided in stacked rasters
  for (v in vars) {
    dt[, (v) := NA_real_]

    path <- if(v == "lai") {
      "data/LAI_MonthlyMean_2011_2020_0.5deg.tif"
    } else {
      paste0("data/worldclim/",data,"/10m/wc2.1_",data,"_10m_", v, ".tif")
    }
    
    # throw an error if the path is empty
    full_path <- system.file(path, package = "SGVM")

    if (full_path == "") {
      stop(
        "Raster not found in sgvm package:\n",
        path
      )
    }

    # load stacked raster (12 layers)
    r <- rast(full_path)
    
    # loop over the layers, representing months
    for (m in 1:12) {
      
      rm <- r[[m]]
      px <- extract(rm, cbind(lon, lat))[, 1]
      
      # fallback to nearest pixel if NA
      if (is.na(px)) {
        cell <- terra::cellFromXY(rm, cbind(lon, lat))
        cell <- terra::adjacent(
          rm,
          cell,
          directions = 8,
          pairs = FALSE
        )[1]
        
        px <- values(rm)[cell]
        
        warning(
          sprintf(
            "Data missing at (lat=%.3f, lon=%.3f) for %s month %d; using nearest pixel.",
            lat, lon, v, m
          ),
          call. = FALSE
        )
      }
      
      dt[month == m, (v) := px]
    }
  }
    
  #load Biomass raster
  r_bio <- rast(
    system.file(
      "data/ESACCI-BIOMASS-L4-AGB-MERGED-50000m-fv6.0.tif",
      package = "SGVM"
    )
  )
  bio <- extract(r_bio, cbind(lon, lat))[,1] * 100 # from Mg/ha to g/m2
  dt[, biomass := bio]
  
  return(dt[])
}

interpolate_data <- function(dt_month, INTER) {
  
  # Ensure df_month is data.table
  dt_month <- as.data.table(dt_month)
  
  # Middle-of-month days of year
  month_doy <- c(15, 45, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349)
  
  # Prepare output data.table
  dt_daily <- data.table(doy = 1:365)
  dt_daily[, month := ceiling(doy / 30.5)]
  
  # Climate variable names
  vars <- c("prec","srad","tmin","tmax","vapr","lai","biomass")
  
  # Loop through variables and interpolate
  for (col in vars) {
    if (INTER) {

    # --- Interpolated version ---
    dt_daily[, (col) :=
               approx(
                 x = month_doy,
                 y = dt_month[[col]],
                 xout = dt_daily$doy,
                 rule = 2
               )$y
    ]

  } else {

    # --- Monthly constant version ---
    dt_daily[, (col) :=
               dt_month[[col]][month]
    ]
  }
  }
  
  return(dt_daily)
}

calc_assimilation <- function(dt, par, kdif = 0.7) {
  
  stopifnot(is.data.table(dt))
  
  setDT(dt)
  
  # ----- Ensure sorted for cumulative LAI -----
  setorder(dt, day, tod, cohort)
  
  # ----- Cumulative LAI within day -----
  dt[, cumLAI := cumsum(lai_coh), by = .(day,tod)]
  dt[, cumLAI_above := shift(cumLAI, fill = 0), by = .(day,tod)]
  
  # Direct extinction coefficient
  dt[, kdir := sin(theta)]
  
  # ----- Light entering each cohort -----
  dt[, Idir_in := Idir * exp(-kdir * cumLAI_above)]
  dt[, Idif_in := Idif * exp(-kdif * cumLAI_above)]
  
  # ----- Light exiting each cohort -----
  dt[, Idir_out := ifelse(dt$cohort>par()$n_cohorts, 0, Idir * exp(-kdir * cumLAI))]
  dt[, Idif_out := ifelse(dt$cohort>par()$n_cohorts, 0, Idif * exp(-kdif * cumLAI))]
  
  # ----- Light intercepted -----
  dt[, intercepted_dir := Idir_in - Idir_out]
  dt[, intercepted_dif := Idif_in - Idif_out]
  
  # ----- Sunlit fraction -----
  dt[, f_sun := exp(-kdir * cumLAI_above)]
  
  # --- Assimilation: sun leaves ---
  dt[dt$cohort<=par()$n_cohorts, c("A_sun", "Tr_sun", "gs", "vpd") := {
    res <- mapply(
      calcA,
      PPFD  = intercepted_dir*cdir + intercepted_dif*cdif,
      Ca    = Ca,
      TleafC= Temp,
      VP    = vapr*1000,
      O2    = 210,
      LN    = 1,
      gs    = ifelse(tod<par()$st_closure, NA, par()$gs_md),
      SIMPLIFY = FALSE
    )
    
    .( sapply(res, `[[`, "An"),
       sapply(res, `[[`, "Tr"),
       sapply(res, `[[`, "gs"),
       sapply(res, `[[`, "vpd") )
  }]
  
  # --- Assimilation: shade leaves ---
  dt[dt$cohort<=par()$n_cohorts, c("A_shade", "Tr_shade", "gs", "vpd") := {
    res <- mapply(
      calcA,
      PPFD  = intercepted_dif*cdif,
      Ca    = Ca,
      TleafC= Temp,
      VP    = vapr*1000,
      O2    = 210,
      LN    = 1,
      gs    = ifelse(tod<par()$st_closure, NA, par()$gs_md),
      SIMPLIFY = FALSE
    )
    
    .( sapply(res, `[[`, "An"),
       sapply(res, `[[`, "Tr"),
       sapply(res, `[[`, "gs"),
       sapply(res, `[[`, "vpd") )
  }]

  # --- Potential assimilation per cohort (g CO2 m-2 ground day-1) ---
  dt[, Assim_pot :=
       (A_sun  * lai_coh * f_sun +
          A_shade * lai_coh * (1 - f_sun))
          * time_step * 1e-6 * 12.011 #convert from umol CO2/m2 ground/s to g C/m2 ground/timestep
  ]
  # Potential transpiration in L/m2
  dt[dt$cohort<=par()$n_cohorts, Tr := 
       (Tr_sun  * lai_coh * f_sun +
        Tr_shade * lai_coh * (1 - f_sun)) 
        * time_step
  ]

  # The last cohort represents the undergrowth and soil, from which water transpires and evaporates, 
  # calculated using the Priestly Taylor equation, a simplified version of Pennman Monteith
  # We assume that 1% of all soil water can be evaporated in a day, as drying out of the topsoil layer causes resistance, limiting evaporation
  dt[cohort > par()$n_cohorts,
    Tr := pmin(
      par()$Wmax * 0.01 / par()$n_steps,
      calcPET(
        Irr  = intercepted_dir + intercepted_dif,
        Temp = Temp,
        time = time_step
      )$PET
    )
  ]
  
  # calculate the water balance
  dt <- calc_water(dt, par)
  
  # calculate GPP
  dt[, Assim_wlim := Assim_pot * f_Tr] # gC / m2

  #get first tod value to add the night's respiration to the first day
  first_tod <- dt$tod[1]

  #calculate average yearly temperature, which will used as the baseline for the temperature response of maintenance respiration,
  #thereby assumining long term accliation to temperature (see Atkin 2015)
  Tavg <- pmax(par()$rmMin, (mean(dt$tmin) + mean(dt$tmax)) / 2)
  
  #calculate maintenance respiration
  dt[dt$cohort<=par()$n_cohorts, rm := biomass / par()$n_cohorts # divide the respiration costs over the cohorts
        * (1-par()$fHW) * par()$rmAvg*par()$rmQ10**((Temp-Tavg)/10) # maintenance respiration rate based on temperature
        * fifelse(tod == first_tod, (time_step + (24-dayLength) * 3600) / 86400, time_step / 86400)] # gC/gC per timestep, add the night to the first time step
  
  #calculate total ecosystem respiration = 0.69 accounts for the conversion of glucose to biomass (Poorter 1997)
  dt[, re := rm + pmax(0, (Assim_wlim - rm)) * (1-par()$ccBIO)]
  
  #calculate net primary productivity
  dt[, NPP := Assim_wlim - re]
  
  return(dt)
}

calc_water <- function(dt, par) {

  # paramters
  p <- par()
  
  # Ensure order by time
  setorder(dt, day, tod, cohort)
  
  # Construct continuous time vector
  steps <- unique(dt[, .(day, tod)])
  setorder(steps, day, tod)
  
  n_steps <- nrow(steps)
  n_rows  <- nrow(dt)
  
  # Preallocate
  f_Tr_vec   <- numeric(n_rows)
  Transp_vec <- numeric(n_rows)
  Evap_vec <- numeric(n_rows)
  Leaching_vec <- numeric(n_rows)
  Evap_buffer_vec <- numeric(n_rows)
  buffer_size <- p$Wmax * 0.02 #assuming a buffer zone of 2cm in a bucket of 1m
  Water_vec  <- dt$Water   # shared bucket
  Snow_vec   <- dt$Snow
  
  # Index list in matching order
  step_index <- lapply(seq_len(n_steps), function(k) {
    dt[day == steps$day[k] & tod == steps$tod[k], which = TRUE]
  })
  
  # Loop through time
  for (k in seq_len(n_steps)) {
    
    idx   <- step_index[[k]]
    W_now <- Water_vec[idx[1]]
    S_now  <- Snow_vec[idx[1]]
    B_now <- Evap_buffer_vec[idx[1]]
    temp = dt$Temp[idx[1]]
    snow_fall = 0
    rain_fall = 0

    #Let's assume that rain only falls once a week, at night, so at the lowest value of tod
    if(steps$day[k]%%p$rain_freq == 0 && steps$tod[k] == steps$tod[1]) {
      if(temp>0) {
        # precipitation falls in the form of water, corrected for timestep, note that precipitation data is /month
        rain_fall <- dt$prec[idx[1]] * p$rain_freq * 86400/2635200
      } else {
        # precipitation falls in the form of snow, corrected for timestep, note that precipitation data is /month
        snow_fall <- dt$prec[idx[1]] * p$rain_freq * 86400/2635200
      }
    }
    
    # snowmelt, assuming a melting rate of 2.75 mm snow per degree day >0C, and a snow density of 10mm snow = 1 mm water
    snow_melt <- min(S_now, max(0, temp) * 0.275 * dt$time_step[idx[1]] / 86400)
    
    # First, add precip and snowmelt to the buffer layer
    B_new <- B_now + rain_fall + snow_melt

    # Then, calculate infiltration when buffer water content is larger than the max buffer size, and reduce B_new
    infiltration <- pmax(0, B_new - buffer_size)
    B_new <- B_new - infiltration

    # potential evaporation only for ground cohort
    ground <- idx[dt$cohort[idx] >  p$n_cohorts] 
    ev_sum <- sum(dt$Tr[ground])

    # Evaporation is limited either by water in the buffer layer, or by potential evap
    evap <- pmin(B_new, ev_sum)
    Evap_vec[ground] <- evap

    # subtract evap from buffer layer
    B_new <- B_new - evap

    # calculate the new water content past the buffer layer
    W_new <- pmax(0, W_now + infiltration)

    # calculate Leaching
    leach <- pmax(0, W_new - p$Wmax)
    Leaching_vec[ground] <- leach

    # update the water content
    W_new <- W_new - leach

    # potential transpiration only for leafy cohorts
    leafy  <- idx[dt$cohort[idx] <= p$n_cohorts]
    tr_vec <- dt$Tr[leafy]          # per-cohort potential transpiration
    tr_sum <- sum(tr_vec)

    # Shared scaling factor: how much of total demand can be met
    transp_total <- if(temp > 0) min(tr_sum, W_new) else 0
    f_tr         <- if(tr_sum > 0) transp_total / tr_sum else 0

    # Apply shared f_tr to each cohort individually
    f_Tr_vec[leafy]   <- f_tr
    Transp_vec[leafy] <- tr_vec * f_tr   # per-cohort actual transpiration

    # Total water withdrawn from soil
    W_new <- max(0, W_new - transp_total)
    
    # compute Water and Snow in the next time step
    if (k < n_steps) {
      Water_vec[ step_index[[k+1]] ] <- pmin(p$Wmax, W_new)
      Snow_vec[ step_index[[k+1]] ] <- pmax(0, S_now + snow_fall - snow_melt)
      Evap_buffer_vec[ step_index[[k+1]] ] <- B_new
    }
  }
  
  # Write back
  dt[, f_Tr := f_Tr_vec]
  dt[, Transp := Transp_vec]
  dt[, Evap := Evap_vec]
  dt[, Buffer := Evap_buffer_vec]
  dt[, Water := Water_vec]
  dt[, Snow := Snow_vec]
  dt[, Leaching := Leaching_vec]
  
  dt
}
