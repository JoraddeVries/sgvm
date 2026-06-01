# Dev notes
# - penology (leaf-on leaf-off slider)
# - PFTs (evergreen-gymno, temperate-deciduous angiosperm, tropical-evergreen-angiosperm)
# - transpiratie output
# - GPP, NPP en hout als output
# - incoming irradiance worldclim -> cloudy_sky

library(data.table)
library(terra)
library(DT)
library(ggplot2)
library(cowplot)

# -------------------------------------------------------
# Set model parameters
# -------------------------------------------------------

# fraction of LAI divided over the different cohorts
cohort_default <- data.table(
  lai_coh = c(0.33, 0.33, 0.33,0)
)

# =======================================================
#                 SHINY APP
# =======================================================
ui <- fluidPage(

  tags$head(
    tags$style("
      /* Inverted leaf-on slider: blue on the tails, grey in the selected middle */
      .inverted-slider .irs-line {
        background: #428bca !important;
      }
      .inverted-slider .irs-bar {
        background: #dee2e6 !important;
        border-top-color:    #dee2e6 !important;
        border-bottom-color: #dee2e6 !important;
      }
    "),
    tags$script("
      $(document).on('change', '#lai_invert', function() {
        $('#lai_range').closest('.shiny-input-container').toggleClass('inverted-slider', this.checked);
      });
    ")
  ),

  titlePanel("Static Global Vegetation Model"),
  
  sidebarLayout(
    
    sidebarPanel(
      width = 4,
      
      tags$h4("Climate Data"),

      br(),
      radioButtons(
        inputId = "climate_scenario",
        label = NULL,
        choices = c(
          "Historic 1970-2000" = "historic",
          "SSP1-2.6 2081–2100" = "UKESM_ssp126_2081-2100",
          "SSP2-4.5 2081–2100" = "UKESM_ssp245_2081-2100",
          "SSP3-7.0 2081–2100" = "UKESM_ssp370_2081-2100"
        ),
        selected = "historic"
      ),
      sliderInput("Ca", "Atmospheric CO2 (ppm)", min = 200, max = 1200, value = 400, step = 10),

      br(),
      radioButtons(
        inputId = "location",
        label = NULL,
        choices = c(
          "Netherlands" = "ned",
          "Canada" = "can",
          "Brazil" = "bra"
        ),
        selected = "ned"
      ),
      sliderInput("latitude", "Latitude (deg)", min = -90, max = 90, value = 52, step = 1),
      sliderInput("longitude", "Longitude (deg)", min = -180, max = 180, value = 6, step = 1),
      sliderInput(
        inputId = "lai_range",
        label   = "Leaf-on period (day of year)",
        min     = 1,
        max     = 365,
        value   = c(1, 365),  # default start / end
        step    = 1,
        sep     = ""           # no thousands separator
      ),
      checkboxInput(
        inputId = "lai_invert",
        label   = "Invert (southern hemisphere)",
        value   = FALSE
      ),
      sliderInput(
        inputId = "lai_coh_split",
        label   = "LAI distribution over cohorts (%, top → bottom)",
        min     = 0,
        max     = 100,
        value   = c(33, 67),
        step    = 1,
        sep     = ""
      ),
      textOutput("lai_coh_display"),

      br(),
      actionButton("load_data", "Load Climate and Vegetation Data", class = "btn-warning"),
      
      br(), br(),
      sliderInput("runtime", "Runtime (years)", min = 1, max = 5, value = 1, step = 1),
      
      br(),
      actionButton("run_btn", "Run Model", class = "btn-primary"),
      
      br(), br(),
      downloadButton("download_output",      "Download output (CSV)"),
      downloadButton("download_output_xlsx", "Download output (XLSX)")
    ),
    
    # -------------------------------------------------------
    # MAIN PANEL WITH TABS
    # -------------------------------------------------------
    mainPanel(
      width = 8,
      tabsetPanel(
        
        tabPanel("Input table", 
          DTOutput("clim_table"),
          br(),
          uiOutput("cell_editor"),
          br(),
          selectInput("adjust_col", "Change all values of variable:", choices = c("prec", "tmin", "tmax", "vapr","lai","biomass")),
          numericInput("adjust_val", "by:", value = 0, step = 0.1),
          actionButton("adjust_btn", "Apply change", class = "btn-primary"),
          tags$h4("Notes"),
          textAreaInput(
            "input_notes",
            label = NULL,
            placeholder = "prec is precipitation in mm/m2/month
srad is incoming solar radiation in kJ/m2/day
tmin and tmax are daily minimum and maximum temperatures in dC
vapr is water vapour pressure in kPa
lai is leaf area index in m2 leaf area / m2 ground area
biomass is total vegetation biomass in g/m2",
            width = "100%",
            height = "150px"
          ),
        ),
        tabPanel("Map", plotOutput("clim_map", height = "500px")),
        tabPanel("PP output", plotOutput("assim_plot", height = "500px")),
        tabPanel("Water output", plotOutput("water_plot", height = "500px")),
        # tabPanel(
        #   "Woody growth output", plotOutput("wood_growth_plot", height = "500px")
        # ),
        tabPanel("Output Summary",
          br(),
          DTOutput("monthly_summary"),
          br(),
          tags$h4("Notes"),
          textAreaInput(
            "summary_notes",
            label = NULL,
            placeholder = "The output table contains:
- four input variables (prec, tmin, tmax, biomass)
- two interpolated variables (tavg, lai)
- seven output variables (Ev, Tr, Le, RE, TR, GPP, NPP)

Month 13 are the yearly sums/averages

Units:
  prec in mm/m2/month(year)
  tmin, tmax, and tavg in dC
  biomass in g
  lai in m2 leaf area / m2 ground area
  RE is ecosystem respiration (maintenance + growth) in gC/month(year)
  Ev, Tr, and Le are Evaporation, Transpiration and Leaching in mm/month(year)
  GPP is gross primary productivity in gC/m2/month(year)
  NPP is net primary productivity (GPP - RE) in gBiomass/m2/month(year)",
            width = "100%",
            height = "150px"
          ),
          br(),
          downloadButton("download_summary", "Download summary (CSV)")
        ),
        tabPanel("Physiology",
          br(),
          fluidRow(
            column(3,
              tags$h4("Vary driver"),
              radioButtons(
                "phys_driver",
                label    = NULL,
                choices  = c("Light" = "light", "Temperature" = "temp", "CO2" = "co2"),
                selected = "light"
              ),
              hr(),
              tags$h4("Hold constant"),
              uiOutput("phys_controls")
            ),
            column(9,
              plotOutput("phys_plot", height = "500px")
            )
          )
        )
      )
    )
  )
)

# =======================================================
#                    SERVER
# =======================================================
server <- function(input, output, session) {

  # initiate parameters
  par <- reactiveVal(init_parameters())
  
  # Holds climate table (editable) and the original loaded values for change highlighting
  clim_data          <- reactiveVal(NULL)
  original_clim_data <- reactiveVal(NULL)
  clim_raster <- reactiveVal(NULL)
  current_month_str <- reactive({as.integer(format(Sys.Date(), "%m"))})   # 1-12
  
  # When "load data" is pressed:
  observeEvent(input$load_data, {

    # update par
    p <- par()
    p$latitude  <- input$latitude
    p$longitude <- input$longitude
    p$Ca        <- input$Ca
    par(p)

    # load climate data
    w <- SGVM::get_data(input$latitude, input$longitude, input$climate_scenario)
    clim_data(w)
    original_clim_data(copy(w))
    
    # Load a representative WorldClim raster to plot (example: January tavg)
    r <- rast(
      system.file(
        paste0("data/LAI_AnnualMaxMean_2011_2020_0.5deg.tif"),
        package = "SGVM"
      )
    )[[1]]
    
    clim_raster(r)
  })

  # Update Ca slider when the climate scenario is changed
  observeEvent(input$climate_scenario, {

    new_Ca <- switch(input$climate_scenario,
      "historic" = 400,
      "UKESM_ssp126_2081-2100" = 450,
      "UKESM_ssp245_2081-2100" = 550,
      "UKESM_ssp370_2081-2100" = 850
    )

    updateSliderInput(
      session,
      "Ca",
      value = new_Ca
    )

  })

  # Update lat-lon and initial snowpack sliders when the location is changed
  observeEvent(input$location, {

    new_lat <- switch(input$location,
      "ned" = 52,
      "can" = 66,
      "bra" = -2
    )

    new_lon <- switch(input$location,
      "ned" = 6,
      "can" = -138,
      "bra" = -59
    )

    updateSliderInput(
      session,
      "latitude",
      value = new_lat
    )

    updateSliderInput(
      session,
      "longitude",
      value = new_lon
    )

  })

  # -------------------------------------------------------
  # Climate and vegetation table
  # -------------------------------------------------------
  
  output$clim_table <- renderDT({
    req(clim_data())

    editable_cols <- c("prec", "tmin", "tmax", "vapr", "lai", "biomass")
    df   <- round(copy(clim_data()), 2)
    orig <- original_clim_data()

    if (!is.null(orig)) {
      orig_r <- round(copy(orig), 2)
      for (col in editable_cols) {
        if (col %in% names(df)) {
          df[[paste0(col, "_color")]] <- ifelse(
            df[[col]] > orig_r[[col]], "up",
            ifelse(df[[col]] < orig_r[[col]], "down", "")
          )
        }
      }
    }

    color_col_idx <- which(grepl("_color$", names(df))) - 1  # 0-based for DT

    tbl <- datatable(
      df,
      rownames = FALSE,
      editable = FALSE,
      selection = list(mode = "single", target = "cell"),
      options = list(
        pageLength = 12,
        dom = "t",
        columnDefs = if (length(color_col_idx) > 0)
          list(list(visible = FALSE, targets = color_col_idx))
        else
          list()
      )
    )

    for (col in editable_cols) {
      color_col <- paste0(col, "_color")
      if (color_col %in% names(df)) {
        tbl <- formatStyle(tbl, col, valueColumns = color_col,
                           backgroundColor = styleEqual(
                             c("up", "down"),
                             c("#c3e6cb", "#f5c6cb")
                           ))
      }
    }

    tbl
  })

  # Cohort LAI fractions derived from the two partition handles
  cohort_fracs <- reactive({
    s <- input$lai_coh_split
    c(s[1], s[2] - s[1], 100 - s[2], 0) / 100
  })

  output$lai_coh_display <- renderText({
    f <- round(cohort_fracs() * 100)
    paste0("C1 (top): ", f[1], "%  |  C2: ", f[2], "%  |  C3 (bottom): ", f[3], "%")
  })

  #Table editors
  selected_cell <- reactiveVal(NULL)

  observeEvent(input$clim_table_cells_selected, {
    selected_cell(input$clim_table_cells_selected)
  })

  output$cell_editor <- renderUI({
    cell <- selected_cell()
    dt   <- clim_data()

    # Nothing clicked yet
    if (is.null(cell) || length(cell) != 2) {
      return(tags$em("Click a cell in the climate table to edit it."))
    }

    row <- cell[1]
    col <- cell[2] + 1   # 🔑 convert DT → R indexing

    colname <- names(dt)[col]

    editable_cols <- c("prec", "tmin", "tmax", "vapr", "lai", "biomass")

    # Non-editable column → show message instead of error
    if (!colname %in% editable_cols) {
      return(
        tags$div(
          tags$hr(),
          tags$h4("Selected variable cannot be edited."),
          tags$p(
            paste(
              "Valiable:", colname,
              "| Month:", dt$month[row]
            )
          )
        )
      )
    } 
    
    tagList(
      tags$hr(),
      tags$h4(
        paste(
          "Variable", colname,
          "| Month", dt$month[row]
        )
      ),
      numericInput(
        inputId = "cell_value",
        label   = NULL,
        value   = dt[[colname]][row],
        step    = 0.1
      ),
      actionButton(
        "apply_cell_edit",
        "Apply change",
        class = "btn-primary"
      )
    )
  })

  observeEvent(input$apply_cell_edit, {
    cell <- selected_cell()
    req(cell)

    df <- clim_data()

    row <- cell[1]
    col <- cell[2] + 1   # 🔑 DT → R

    colname <- names(df)[col]

    val <- suppressWarnings(as.numeric(input$cell_value))
    if (is.na(val)) return()

    # Enforce non-negativity
    if (colname %in% c("prec", "vapr", "lai", "biomass")) {
      val <- max(0, val)
    }

    df[row, colname] <- val
    clim_data(df)
  })

  # Column-wise edits
  observeEvent(input$adjust_btn, {
    req(clim_data())
    
    col <- input$adjust_col
    val <- input$adjust_val
    
    df <- clim_data()
    
    if (col %in% names(df)) {
      if (col %in% c("prec", "vapr", "lai", "biomass")) { # precipitation and vp can't become negative
        df[[col]] <- pmax(0, as.numeric(df[[col]]) + val)
      } else {                          #temperature can become negative
        df[[col]] <- df[[col]] + val 
      }
      clim_data(df)          # update reactive value
    }
  })
  
  # -------------------------------------------------------
  # Run Model (button)
  # -------------------------------------------------------
  model_results <- eventReactive(input$run_btn, {

    p <- par()
    p$n_days    <- input$runtime * 365
    par(p)
    
    # --- Build base dt ---
    dt <- data.table(
      day = rep(1:par()$n_days, each = (par()$n_cohorts+1) * par()$n_steps),
      tod = rep(1:par()$n_steps, each = par()$n_cohorts+1) / par()$n_steps - 0.5/par()$n_steps,
      cohort = rep(1:(par()$n_cohorts+1), times = par()$n_days * par()$n_steps)
    )

    # add doy
    dt[, doy := ((day - 1) %% 365) + 1]
    
    # make sure the model runs also if the data wasn't loaded
    clim_table <- clim_data()
    if (is.null(clim_table)) {
      clim_table <- SGVM::get_data(input$latitude, input$longitude, input$climate_scenario)
    }
    
    # interpolate the climate and vegetation data to days
    dt <- SGVM::set_environment(dt, clim_table, par)

    # modify lai based on phenology input
    dt[, phenology := fifelse(xor(doy >= input$lai_range[1] & doy <= input$lai_range[2], input$lai_invert), 1, 0)]
    dt[, lai_coh := cohort_fracs()[cohort] * lai * phenology]

    # calculate assimilation
    dt <- SGVM::calc_assimilation(dt, par)
    
    # Return full dt so all tabs can use it
    return(dt)
  })
  
  # -------------------------------------------------------
  # Output summary (monthly)
  # -------------------------------------------------------
  monthly_summary <- reactive({
    req(model_results())
    
    #filter the results to only include the second year of simulation
    dt <- model_results()
    clim <- clim_data()
    
    # --- Add month to dt ---
    dt[, month := as.integer((doy - 1) / 30.4375) + 1]
    dt[month > 12, month := 12]
    
    # --- Assimilation & NPP summaries ---
    monthly <- dt[, .(
      tavg        = mean(.SD$Temp,    na.rm = TRUE),
      lai        = mean(.SD$lai,    na.rm = TRUE),
      Ev         = sum(Evap,  na.rm = TRUE),
      Tr         = sum(Transp,  na.rm = TRUE),
      Le         = sum(Leaching,  na.rm = TRUE),
      RE         = sum(re,  na.rm = TRUE),
      GPP        = sum(Assim_wlim, na.rm = TRUE),
      NPP        = sum(NPP,        na.rm = TRUE)
    ), by = month]
    
    # --- Climate summaries ---
    clim_month <- clim[, .(
      prec = sum(prec, na.rm = TRUE),
      tmin = mean(tmin, na.rm = TRUE),
      tmax = mean(tmax, na.rm = TRUE),
      biomass = mean(biomass,   na.rm = TRUE)
    ), by = month]
    
    # --- Merge ---
    summary <- merge(clim_month, monthly, by = "month", all = TRUE)
    
    # --- Yearly row (month 13) ---
    yearly <- summary[, .(
      month      = 13,
      prec       = sum(prec,       na.rm = TRUE),
      tmin       = mean(tmin,      na.rm = TRUE),
      tmax       = mean(tmax,      na.rm = TRUE),
      biomass    = mean(biomass,   na.rm = TRUE),
      tavg       = mean(tavg,      na.rm = TRUE),
      lai        = mean(lai,       na.rm = TRUE),
      Ev         = sum(Ev,         na.rm = TRUE),
      Tr         = sum(Tr,         na.rm = TRUE),
      Le         = sum(Le,         na.rm = TRUE),
      RE         = sum(RE,         na.rm = TRUE),
      GPP        = sum(GPP,        na.rm = TRUE),
      NPP        = sum(NPP,        na.rm = TRUE)
    )]
    
    rbind(summary, yearly)
  })
  
  output$monthly_summary <- renderDT({
    req(monthly_summary())
    
    datatable(
      round(monthly_summary(), 2),
      rownames = FALSE,
      options = list(
        dom = "t",
        pageLength = 13
      )
    )
  })
  
  output$download_summary <- downloadHandler(
    filename = function() {
      paste0("monthly_summary_", Sys.Date(), ".csv")
    },
    content = function(file) {
      fwrite(monthly_summary(), file)
    }
  )
  
  # -------------------------------------------------------
  # Assimilation graph (daily)
  # -------------------------------------------------------
  output$assim_plot <- renderPlot({
    dt <- model_results()
    req(dt)
    
    dt_day <- dt[, .(
      RE = sum(re, na.rm = TRUE),
      GPP = sum(Assim_wlim, na.rm = TRUE),
      NPP = sum(NPP, na.rm = TRUE)
    ), by = day]
    
    # Reshape to long format for easy legend handling
    dt_long <- melt(
      dt_day,
      id.vars = "day",
      measure.vars = c("GPP", "RE", "NPP"),
      variable.name = "Type",
      value.name = "value"
    )
    
    ggplot(dt_long, aes(x = day, y = value, color = Type)) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept=0, linetype="dashed") +
      scale_color_manual(
        values = c("GPP" = "blue", "RE" = "red", "NPP" = "black"),
        labels = c("GPP", "RE", "NPP")
      ) +
      labs(
        x = "Day",
        y = "Carbon flux (g/m²)",
        title = "Daily Primary Productivity",
        color = ""   # legend title
      ) +
      theme_cowplot() +
      theme(
        legend.position = "bottom",
        panel.background = element_rect(fill = NA, color = "black")
      )
  })
  
  # -------------------------------------------------------
  # Water graph (daily mean water)
  # -------------------------------------------------------
  output$water_plot <- renderPlot({
    dt <- model_results()
    req(dt)
    
    dt_water <- dt[, .(
      Water_mean = mean(Water),
      Snow_mean = mean(Snow),
      Tr_mean = sum(Tr) * 10 #increase value of transpiration for visualisation
    ), by = day]
    
    # Reshape to long format for easy legend handling
    dt_water_long <- melt(
      dt_water,
      id.vars = "day",
      measure.vars = c("Water_mean", "Snow_mean", "Tr_mean"),
      variable.name = "Type",
      value.name = "value"
    )
    
    ggplot(dt_water_long, aes(x = day, y = value, color = Type)) +
      geom_line(linewidth = 1) +
      scale_color_manual(
        values = c("Snow_mean" = "black", "Water_mean" = "blue", "Tr_mean" = "green"),
        labels = c("Soil Water", "Snowpack", "Evapotranspiration (x10)")
      ) +
      labs(
        x = "Day",
        y = "Water content of Soil/Snowpack (L/m²)",
        title = "Water",
        color = ""   # legend title
      ) +
      theme_cowplot() +
      theme(
        legend.position = "bottom",
        panel.background = element_rect(fill = NA, color = "black")
      )
  })

  # -------------------------------------------------------
  # Wood growth tab
  # -------------------------------------------------------   
  # wood_growth <- reactive({
  #   req(model_results())
  #   calc_woody_growth(model_results())
  # })

  # output$wood_growth_plot <- renderPlot({

  #   df <- wood_growth()
  #   req(nrow(df) > 0)

  #   # reshape to long format
  #   plot_dt <- melt(
  #     df,
  #     id.vars = c("cohort_id", "age"),
  #     measure.vars = c("size", "wall"),
  #     variable.name = "variable",
  #     value.name = "value"
  #   )

  #   ggplot(plot_dt, aes(x = cohort_id, y = value)) +
  #     geom_line(alpha = 0.4) +
  #     geom_point() +
  #     facet_wrap(
  #       ~ variable,
  #       scales = "free_y",
  #       labeller = as_labeller(
  #         c(
  #           size = "Cell size (µm)",
  #           wall = "Cell wall thickness (µm)"
  #         )
  #       )
  #     ) +
  #     labs(
  #       x = "Cell cohort",
  #       y = NULL
  #     ) +
  #     theme_minimal() +
  #     theme(
  #       strip.text = element_text(face = "bold"),
  #       panel.grid.minor = element_blank()
  #     )
  # })

 
  # -------------------------------------------------------
  # Raster + lat lon 
  # -------------------------------------------------------   
  output$clim_map <- renderPlot({
    req(clim_raster())
    
    r <- clim_raster()
    
    plot(
      r,
      main = paste0("Annual max LAI"),
      col = terrain.colors(100)
    )
    
    points(
      x = input$longitude,
      y = input$latitude,
      pch = 19,
      col = "red",
      cex = 1.5
    )
  })
  
  # -------------------------------------------------------
  # Download the output table
  # -------------------------------------------------------
  output$download_output <- downloadHandler(
    filename = function() paste0("model_output_", Sys.Date(), ".csv"),
    content  = function(file) fwrite(model_results(), file)
  )

  output$download_output_xlsx <- downloadHandler(
    filename = function() paste0("model_output_", Sys.Date(), ".xlsx"),
    content  = function(file) writexl::write_xlsx(model_results(), file)
  )

  # -------------------------------------------------------
  # Physiology response curves
  # -------------------------------------------------------
  output$phys_controls <- renderUI({
    switch(input$phys_driver,
      light = tagList(
        sliderInput("phys_temp", "Temperature (deg C)", min = -5,  max = 45,   value = 20,  step = 1),
        sliderInput("phys_ca",   "CO2 (ppm)",           min = 200, max = 1000, value = 415, step = 5)
      ),
      temp = tagList(
        sliderInput("phys_ppfd", "Light (umol/m2/s)",   min = 0,   max = 2000, value = 800, step = 50),
        sliderInput("phys_ca",   "CO2 (ppm)",           min = 200, max = 1000, value = 415, step = 5)
      ),
      co2 = tagList(
        sliderInput("phys_ppfd", "Light (umol/m2/s)",   min = 0,   max = 2000, value = 800, step = 50),
        sliderInput("phys_temp", "Temperature (deg C)", min = -5,  max = 45,   value = 20,  step = 1)
      )
    )
  })

  phys_data <- reactive({
    driver <- req(input$phys_driver)
    ppfd   <- if (!is.null(input$phys_ppfd)) input$phys_ppfd else 800
    temp   <- if (!is.null(input$phys_temp)) input$phys_temp else 20
    ca     <- if (!is.null(input$phys_ca))   input$phys_ca   else 415

    run_calcA <- function(p, t, c) {
      tryCatch(
        SGVM:::calcA(PPFD = p, Ca = c, TleafC = t, VP = 1000, O2 = 210000, LN = 2),
        error = function(e) list(An = NA_real_, Tr = NA_real_)
      )
    }

    x <- switch(driver,
      light = seq(0,   2000, by = 20),
      temp  = seq(-5,  45,   by = 1),
      co2   = seq(200, 1000, by = 10)
    )

    results <- switch(driver,
      light = lapply(x, function(p) run_calcA(p, temp, ca)),
      temp  = lapply(x, function(t) run_calcA(ppfd, t,   ca)),
      co2   = lapply(x, function(c) run_calcA(ppfd, temp, c))
    )

    An <- sapply(results, function(r) r$An)
    Tr <- sapply(results, function(r) r$Tr) * 55508  # L/m2/s -> mmol/m2/s

    # Maintenance respiration, Tavg fixed at 15 deg C as temperate acclimation reference
    Rm_temps <- if (driver == "temp") x else rep(temp, length(x))
    Rm <- SGVM:::calc_rm_rate(Rm_temps, par(), 15) * 1000  # g/g/day -> mg/g/day

    data.frame(x = x, An = An, Tr = pmax(0, Tr), Rm = Rm)
  })

  output$phys_plot <- renderPlot({
    df     <- phys_data()
    driver <- input$phys_driver

    x_label <- switch(driver,
      light = "Light (umol m-2 s-1)",
      temp  = "Temperature (deg C)",
      co2   = "CO2 (ppm)"
    )

    th <- theme_classic(base_size = 13) +
      theme(plot.title = element_text(face = "bold"))

    p1 <- ggplot(df, aes(x = x, y = An)) +
      geom_line(colour = "#2c7bb6", linewidth = 1) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
      labs(title = "Photosynthesis", x = x_label,
           y = expression(mu*mol~CO[2]~m^{-2}~s^{-1})) +
      th

    p2 <- ggplot(df, aes(x = x, y = Tr)) +
      geom_line(colour = "#1a9641", linewidth = 1) +
      labs(title = "Transpiration", x = x_label,
           y = expression(mmol~H[2]*O~m^{-2}~s^{-1})) +
      th

    p3 <- ggplot(df, aes(x = x, y = Rm)) +
      geom_line(colour = "#d7191c", linewidth = 1) +
      labs(title = "Maintenance respiration", x = x_label,
           y = expression(mg~g^{-1}~day^{-1})) +
      th

    cowplot::plot_grid(p1, p2, p3, nrow = 1)
  })

 }

# Run app
shinyApp(ui, server)
