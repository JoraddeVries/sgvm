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
# Default model parameters
# -------------------------------------------------------
par_default <- list(
  n_steps = 6,
  n_cohorts = 3,
  latitude = 60,
  longitude = 100,
  Ca = 400,
  Wmax = 300,
  Winit = 0.1,
  Sinit = 80
)

# fraction of LAI divided over the different cohorts
cohort_default <- data.table(
  LAI = c(0.33, 0.33, 0.33)
)

# =======================================================
#                 SHINY APP
# =======================================================
ui <- fluidPage(
  
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
          "SSP5-8.5 2081–2100" = "UKESM_ssp585_2081-2100"
        ),
        selected = "historic"
      ),
      sliderInput("latitude", "Latitude", min = -90, max = 90, value = par_default$latitude, step = 1),
      sliderInput("longitude", "Longitude", min = -180, max = 180, value = par_default$longitude, step = 1),
      sliderInput("Ca", "Atmospheric CO2", min = 200, max = 1200, value = par_default$Ca, step = 10),
      sliderInput(
        inputId = "lai_range",
        label   = "Leaf-on period (day of year)",
        min     = 1,
        max     = 365,
        value   = c(1, 365),  # default start / end
        step    = 1,
        sep     = ""           # no thousands separator
      ),
      
      actionButton("load_data", "Load Climate and Vegetation Data", class = "btn-warning"),
      
      br(), br(),
      tags$h4("Initial Conditions"),
      
      sliderInput("Wmax", "Maximum Soil Water Content (L/m²)",min = 50, max = 800, value = par_default$Wmax, step = 10),
      sliderInput("Winit", "Initial Soil Water Content (fraction)", min = 0, max = 1, value = par_default$Winit, step = 0.1),
      sliderInput("Sinit", "Initial Snowpack (L water/m²)", min = 0, max = 100, value = par_default$Sinit, step = 5),
      
      br(),
      actionButton("run_btn", "Run Model", class = "btn-primary"),
      
      br(), br(),
      downloadButton("download_output", "Save Model Output")
    ),
    
    # -------------------------------------------------------
    # MAIN PANEL WITH TABS
    # -------------------------------------------------------
    mainPanel(
      width = 8,
      tabsetPanel(
        
        tabPanel("Input Table", 
          DTOutput("clim_table"),
          br(),
          uiOutput("cell_editor"),
          br(),
          selectInput("adjust_col", "Change all values of variable:", choices = c("prec", "tmin", "tmax", "vapr","lai_tot","biomass")),
          numericInput("adjust_val", "by:", value = 0, step = 0.1),
          actionButton("adjust_btn", "Apply change", class = "btn-primary"),
          tags$h4("Notes"),
          textAreaInput(
            "input_notes",
            label = NULL,
            placeholder = "prec is precipitation in mm/m2/month\ntmin and tmax are daily minimum and maximum temperatures in dC\nvapr is water vapour pressure in kPa\nlai_tot is total leaf area index in m2 leaf area / m2 ground area\nbiomass is total vegetation biomass in g/m2",
            width = "100%",
            height = "150px"
          ),
        ),
        tabPanel("Map", plotOutput("clim_map", height = "500px")),
        tabPanel("Assimilation", plotOutput("assim_plot", height = "500px")),
        tabPanel("Water", plotOutput("water_plot", height = "500px")),
        tabPanel("Output Summary",
          br(),
          DTOutput("monthly_summary"),
          br(),
          tags$h4("Notes"),
          textAreaInput(
            "summary_notes",
            label = NULL,
            placeholder = "Month 13 are the yearly sums/averages\nprec in mm/m2/month(year)\ntmin and tmax in dC\nAssim_pot in gC/m2/month(year)\nAssim_wlim in gC/m2/month(year)\nNPP in gBiomass/m2/month(year)",
            width = "100%",
            height = "150px"
          ),
          br(),
          downloadButton("download_summary", "Download summary (CSV)")
        ),
        tabPanel(
          "Wood Growth", plotOutput("wood_growth_plot", height = "500px")
        )
      )
    )
  )
)


# =======================================================
#                    SERVER
# =======================================================
server <- function(input, output, session) {
  
  # Holds climate table (editable)
  clim_data <- reactiveVal(NULL)
  clim_raster <- reactiveVal(NULL)
  current_month_str <- reactive({as.integer(format(Sys.Date(), "%m"))})   # 1-12
  
  # When "load data" is pressed:
  observeEvent(input$load_data, {
    w <- SGVM::get_data(input$latitude, input$longitude, input$climate_scenario)
    clim_data(w)
    
    # Load a representative WorldClim raster to plot (example: January tavg)
    r <- rast(
      system.file(
        paste0("data/worldclim/",input$climate_scenario,"/10m/wc2.1_",input$climate_scenario,"_10m_tmax.tif"),
        package = "SGVM"
      )
    )[[current_month_str()]]
    
    clim_raster(r)
  })

  # -------------------------------------------------------
  # Climate and vegetation table
  # -------------------------------------------------------
  
  output$clim_table <- renderDT({
    req(clim_data())
    
    df <- copy(clim_data())
    
    datatable(
      round(df, 2),
      rownames = FALSE,
      editable = FALSE,
      selection = list(
        mode = "single",
        target = "cell"
      ),
      options = list(
        pageLength = 12,
        dom = "t"
      )
    )
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

    editable_cols <- c("prec", "tmin", "tmax", "vapr", "lai_tot", "biomass")

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
    if (colname %in% c("prec", "vapr", "lai_tot", "biomass")) {
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
      if (col %in% c("prec", "vapr", "lai_tot", "biomass")) { # precipitation and vp can't become negative
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
    
    # --- Build parameter set ---
    par <- list(
      n_steps = par_default$n_steps,
      n_cohorts = par_default$n_cohorts,
      latitude = input$latitude,
      longitude = input$longitude,
      Ca = input$Ca,
      Wmax = input$Wmax,
      Winit = input$Winit,
      Sinit = input$Sinit
    )
    
    # --- Build base dt ---
    dt <- data.table(
      doy = rep(1:365, each = par$n_cohorts * par$n_steps),
      tod = rep(1:par$n_steps, each = par$n_cohorts) / par$n_steps - 0.5/par$n_steps,
      cohort = rep(1:par$n_cohorts, times = 365 * par$n_steps)
    )
    
    # make sure the model runs also if the data wasn't loaded
    clim_table <- clim_data()
    if (is.null(clim_table)) {
      clim_table <- SGVM::get_data(input$latitude, input$longitude, input$climate_scenario)
    }
    
    # interpolate the climate and vegetation data to days
    dt <- SGVM::set_environment(dt, clim_table, par)

    # modify lai based on phenology input
    dt[, phenology := fifelse(doy >= input$lai_range[1] & doy <= input$lai_range[2], 1,0)]
    dt[, LAI := cohort_default$LAI[cohort] * lai_tot * phenology]

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
    
    dt <- model_results()
    clim <- clim_data()
    
    # --- Add month to dt ---
    dt[, month := as.integer((doy - 1) / 30.4375) + 1]
    dt[month > 12, month := 12]
    
    # --- Assimilation & NPP summaries ---
    assim_month <- dt[, .(
      LAI        = mean(.SD$lai_tot,    na.rm = TRUE),
      biomass    = mean(biomass,   na.rm = TRUE),
      Assim_pot  = sum(Assim_pot,  na.rm = TRUE),
      Assim_wlim = sum(Assim_wlim, na.rm = TRUE),
      NPP        = sum(NPP,        na.rm = TRUE)
    ), by = month]
    
    # --- Climate summaries ---
    clim_month <- clim[, .(
      prec = sum(prec, na.rm = TRUE),
      tmin = mean(tmin, na.rm = TRUE),
      tmax = mean(tmax, na.rm = TRUE)
    ), by = month]
    
    # --- Merge ---
    summary <- merge(clim_month, assim_month, by = "month", all = TRUE)
    
    # --- Yearly row (month 13) ---
    yearly <- summary[, .(
      month      = 13,
      prec       = sum(prec,       na.rm = TRUE),
      tmin       = mean(tmin,      na.rm = TRUE),
      tmax       = mean(tmax,      na.rm = TRUE),
      LAI        = mean(LAI,       na.rm = TRUE),
      biomass    = mean(biomass,   na.rm = TRUE),
      Assim_pot  = sum(Assim_pot,  na.rm = TRUE),
      Assim_wlim = sum(Assim_wlim, na.rm = TRUE),
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
      Assim_pot_day = sum(Assim_pot),
      Assim_wlim_day = sum(Assim_wlim),
      NPP_day = sum(NPP)
    ), by = doy]
    
    # Reshape to long format for easy legend handling
    dt_long <- melt(
      dt_day,
      id.vars = "doy",
      measure.vars = c("Assim_pot_day", "Assim_wlim_day", "NPP_day"),
      variable.name = "Type",
      value.name = "value"
    )
    
    ggplot(dt_long, aes(x = doy, y = value, color = Type)) +
      geom_line(linewidth = 1) +
      scale_color_manual(
        values = c("Assim_pot_day" = "black", "Assim_wlim_day" = "blue", "NPP_day" = "red"),
        labels = c("Potential assimilation", "Water-limited assimilation", "NPP")
      ) +
      labs(
        x = "Day of Year",
        y = "Daily Assimilation (g/m²)",
        title = "Daily Canopy Assimilation",
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
      Tr_mean = sum(Tr)
    ), by = doy]
    
    # Reshape to long format for easy legend handling
    dt_water_long <- melt(
      dt_water,
      id.vars = "doy",
      measure.vars = c("Water_mean", "Snow_mean", "Tr_mean"),
      variable.name = "Type",
      value.name = "value"
    )
    
    ggplot(dt_water_long, aes(x = doy, y = value, color = Type)) +
      geom_line(linewidth = 1) +
      scale_color_manual(
        values = c("Snow_mean" = "black", "Water_mean" = "blue", "Tr_mean" = "green"),
        labels = c("Soil Water", "Snow", "Transpiration")
      ) +
      labs(
        x = "Day of Year",
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
  wood_growth <- reactive({
    req(model_results())
    calc_woody_growth(model_results())
  })

  output$wood_growth_plot <- renderPlot({

    df <- wood_growth()
    req(nrow(df) > 0)

    # reshape to long format
    plot_dt <- melt(
      df,
      id.vars = c("cohort_id", "age"),
      measure.vars = c("size", "wall"),
      variable.name = "variable",
      value.name = "value"
    )

    ggplot(plot_dt, aes(x = cohort_id, y = value)) +
      geom_line(alpha = 0.4) +
      geom_point() +
      facet_wrap(
        ~ variable,
        scales = "free_y",
        labeller = as_labeller(
          c(
            size = "Cell size (µm)",
            wall = "Cell wall thickness (µm)"
          )
        )
      ) +
      labs(
        x = "Cell cohort",
        y = NULL
      ) +
      theme_minimal() +
      theme(
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      )
  })

 
  # -------------------------------------------------------
  # Raster + lat lon 
  # -------------------------------------------------------   
  output$clim_map <- renderPlot({
    req(clim_raster())
    
    r <- clim_raster()
    
    plot(
      r,
      main = paste0("WorldClim (tmax, month ",current_month_str(),")"),
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
    
    filename = function() {
      paste0("model_output_", Sys.Date(), ".csv")
    },
    
    content = function(file) {
      # Get the latest model results
      dt <- model_results()
      
      # Save to CSV
      fwrite(dt, file)  # data.table::fwrite for speed
    }
  )
  
 }

# Run app
shinyApp(ui, server)
