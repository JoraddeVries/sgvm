
set_environment <- function(dt, input_data, par, cloud_cover = 0.2) {

  # Add month column to dt
  dt[, month := ceiling(doy / 30.5)]  # approximate month from doy
  
  # Interpolate input data
  input_daily <- interpolate_data(input_data)
  
  # Join monthly precipitation
  dt <- merge(dt, input_daily, by = "doy", all.x = TRUE, sort = FALSE)
  
  # Row-wise calling of clear_sky
  cs <- mapply(
    function(lat, doy, tod) {
      clear_sky(lat = (lat/180*pi), DOY = doy, f = tod)
    },
    par$latitude, dt$doy, dt$tod,
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
  
  # Row-wise call to cloudy_sky
  cls <- mapply(
    function(Ig, lat, doy, tod) {
      cloudy_sky(Ig = Ig * (1-cloud_cover), lat = (lat/180*pi), DOY = doy, f = tod)
    },
    sapply(cs, `[[`, "Ig"),  # use clear-sky Ig as input
    par$latitude, dt$doy, dt$tod,
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
  
  # set timestep length in seconds
  dt[, time_step := dayLength/par$n_steps * 3600]
  
  # Get conversion factors
  cdir <- waveband_conversion(Itype = "direct",  waveband = "PAR", mode = "flux")
  cdif <- waveband_conversion(Itype = "diffuse", waveband = "PAR", mode = "flux")
  
  # Add PAR columns to the data.table
  dt[, PAR_dir := Idir * cdir]
  dt[, PAR_dif := Idif * cdif]
  
  # Optionally, total PAR
  dt[, PAR_total := PAR_dir + PAR_dif]
  
  # atmospheric COS concentration in ppm
  dt[, Ca := par$Ca]
  
  # air temperature
  dt[, Temp := tmin + (tmax - tmin) * sin((pi/2) * (tod / 0.6))]
  
  # Soil Water content
  dt[, Water := par$Wmax * par$Winit]
  
  # Snow
  dt[, Snow := par$Sinit]
}

get_data <- function(lat, lon, data) {
  
  vars <- c("prec","srad","tmin","tmax","vapr","lai_tot")
  
  # initialise output
  dt <- data.table(
    month = 1:12
  )
  
  # loop over the different vars in the column, these are provided in stacked rasters
  for (v in vars) {
    dt[, (v) := NA_real_]

    path <- if(v == "lai_tot") {
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

interpolate_data <- function(dt_month) {
  
  # Ensure df_month is data.table
  dt_month <- as.data.table(dt_month)
  
  # Middle-of-month days of year
  month_doy <- c(15, 45, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349)
  
  # Prepare output data.table
  dt_daily <- data.table(
    doy = 1:365
  )
  
  # Climate variable names
  vars <- c("prec","srad","tmin","tmax","tavg","vapr","wind","lai_tot","biomass")
  
  # Loop through variables and interpolate
  for (col in vars) {
    dt_daily[, (col) :=
               approx(
                 x = month_doy,
                 y = dt_month[[col]],
                 xout = dt_daily$doy,
                 rule = 2
               )$y
    ]
  }
  
  return(dt_daily)
}

calc_assimilation <- function(dt, par, kdif = 0.7) {
  
  stopifnot(is.data.table(dt))
  
  setDT(dt)
  
  # ----- Ensure sorted for cumulative LAI -----
  setorder(dt, doy, tod, cohort)
  
  # ----- Cumulative LAI within day -----
  dt[, cumLAI := cumsum(LAI), by = .(doy,tod)]
  dt[, cumLAI_above := shift(cumLAI, fill = 0), by = .(doy,tod)]
  
  # Direct extinction coefficient
  dt[, kdir := sin(theta)]
  
  # ----- Light entering each cohort -----
  dt[, Idir_in := Idir * exp(-kdir * cumLAI_above)]
  dt[, Idif_in := Idif * exp(-kdif * cumLAI_above)]
  
  # ----- Light exiting each cohort -----
  dt[, Idir_out := Idir * exp(-kdir * cumLAI)]
  dt[, Idif_out := Idif * exp(-kdif * cumLAI)]
  
  # ----- Light intercepted -----
  dt[, intercepted_dir := Idir_in - Idir_out]
  dt[, intercepted_dif := Idif_in - Idif_out]
  
  # ----- Sunlit fraction -----
  dt[, f_sun := exp(-kdir * cumLAI_above)]
  
  # --- Assimilation: sun and shade for each ROW (each cohort) ---
  dt[, c("A_sun", "Tr_sun") := {
    res <- mapply(
      calcA,
      PPFD  = intercepted_dir + intercepted_dif,
      Ca    = Ca,
      TleafC= Temp,
      VP    = vapr*1000,
      O2    = 210,
      LN    = 1,
      SIMPLIFY = FALSE
    )
    
    .( sapply(res, `[[`, "An"),
       sapply(res, `[[`, "Tr") )
  }]
  
  
  dt[, c("A_shade", "Tr_shade") := {
    res <- mapply(
      calcA,
      PPFD  = intercepted_dif,
      Ca    = Ca,
      TleafC= Temp,
      VP    = vapr*1000,
      O2    = 210,
      LN    = 1,
      SIMPLIFY = FALSE
    )
    
    .( sapply(res, `[[`, "An"),
       sapply(res, `[[`, "Tr") )
  }]
  
  # Potential transpiration
  dt[,Tr := 
       (Tr_sun  * LAI * f_sun +
        Tr_shade * LAI * (1 - f_sun)) 
        * time_step
  ]
  
  dt <- calc_water(dt, par)
  
  # --- Total assimilation per cohort (g CO2 m-2 ground per day) ---
  dt[, Assim_pot :=
       (A_sun  * LAI * f_sun +
          A_shade * LAI * (1 - f_sun))
          * time_step * 1e-6 * 44.0095 #convert from umol/m2 ground/s to g/m2 ground/timestep
  ]
  
  dt[, Assim_wlim := Assim_pot * f_Tr] # gC / m2
  
  dt[, rm := 0.0005*1.7**((Temp-15)/10) * time_step/86400] # gC/gC per timestep
  
  dt[, NPP := {
    net <- Assim_wlim - rm * biomass
    ifelse(net > 0, net * 0.69, net) # * 0.69 accounts for the conversion of glucose to biomass (Poorter 1997)
  }]
  
  return(dt)
}

calc_water <- function(dt, par) {
  
  # Ensure order by time
  setorder(dt, doy, tod, cohort)
  
  # Construct continuous time vector
  steps <- unique(dt[, .(doy, tod)])
  setorder(steps, doy, tod)
  
  n_steps <- nrow(steps)
  n_rows  <- nrow(dt)
  
  # Preallocate
  f_Tr_vec   <- numeric(n_rows)
  Uptake_vec <- numeric(n_rows)
  Water_vec  <- dt$Water   # shared bucket
  Snow_vec   <- dt$Snow
  
  # Index list in matching order
  step_index <- lapply(seq_len(n_steps), function(k) {
    dt[doy == steps$doy[k] & tod == steps$tod[k], which = TRUE]
  })
  
  # Loop through time
  for (k in seq_len(n_steps)) {
    
    idx   <- step_index[[k]]
    W_now <- Water_vec[idx[1]]
    S_now  <- Snow_vec[idx[1]]
    tr_sum <- sum(dt$Tr[idx])
    temp = dt$Temp[idx[1]]
    snow_fall = 0
    rain_fall = 0
    
    #Let's assume that rain only falls at night, so at the lowest value of tod
    if(steps$tod[k] == steps$tod[1]) {
      if(temp>0) {
        # precipitation falls in the form of water, corrected for timestep, note that precipitation data is /month
        rain_fall <- dt$prec[idx[1]] * 86400/2635200#* dt$time_step[idx[1]] / 2635200
      } else {
        # precipitation falls in the form of snow, corrected for timestep, note that precipitation data is /month
        snow_fall <- dt$prec[idx[1]] * 86400/2635200#* dt$time_step[idx[1]] / 2635200
      }
    }
    
    # snowmelt, assuming a melting rate of 2.75 mm snow per degree day >0C, and a snow density of 10mm snow = 1 mm water
    snow_melt <- pmin(S_now, pmax(0, temp) * 0.275 * dt$time_step[idx[1]] / 86400)
    
    # uptake
    uptake <- if(temp>0) pmin(tr_sum, W_now + rain_fall + snow_melt) else 0
    f_tr   <- if (tr_sum>0) uptake / tr_sum else 0
    
    # store
    f_Tr_vec[idx]   <- f_tr
    Uptake_vec[idx] <- uptake
    
    # compute Water and Snow in the next time step
    if (k < n_steps) {
      Water_vec[ step_index[[k+1]] ] <- pmin(par$Wmax, pmax(0, W_now - tr_sum + rain_fall + snow_melt))
      Snow_vec[ step_index[[k+1]] ] <- pmax(0, S_now + snow_fall - snow_melt)
    }
  }
  
  # Write back
  dt[, f_Tr := f_Tr_vec]
  dt[, Uptake := Uptake_vec]
  dt[, Water := Water_vec]
  dt[, Snow := Snow_vec]
  
  dt
}



