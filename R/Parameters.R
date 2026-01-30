get_parameters <- function() {

par <- list(

    #simulation config
    n_steps = 6,
    n_cohorts = 3,
    latitude = 60,
    longitude = 100,

    #inititation
    Ca = 400,
    Wmax = 300,
    Winit = 0.1,
    Sinit = 80,

    #plant function
    rm15 = 0.001,           # base mantenance respiration rate at 15 degrees C (g/g)
    rmQ10 = 1.7,            # Q10 for maintenance respiration
    fHW = 0.38,             # fraction of total biomass that is heartwood (see BAAD dataset; Falster et.al. 2015)
    ccBIO = 0.69,           # construction cost of biomass; amount of biomass constructed from 1 g of glucose (g/g)
    gs_md = 0.1             # stomatal conductance when stomates close at midday (mol/m2/s)
    )

return(par)

}