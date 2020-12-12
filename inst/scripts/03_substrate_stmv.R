
## NOTE:: substrate size is really only relevant for SSE/snowcrab domain right now as no
##        other data source has been found/identified
##        but working at the size of canada.east.highres for compatibility with bathymetry

## TODO:: add data collected by snow crab survey and any others for that matter

## TODO --- once settled, these params should be moved into substrate_parameters as defaults

require(aegis.substrate)

p = substrate_parameters(
  project_class="stmv",
  stmv_nmin = 100, # stmv_nmin/stmv_nmax changes with resolution
  stmv_nmax = 1000 # numerical time/memory constraint -- anything larger takes too much time .. anything less .. errors
)


p$stmv_distance_prediction_limits = p$stmv_distance_statsgrid * c( 1/2, 5 ) # range of permissible predictions km (i.e 1/2 stats grid to upper limit based upon data density)
p$stmv_distance_scale = p$stmv_distance_statsgrid * c( 1, 2, 3, 4, 5, 10, 20, 40) # km ... approx guesses of 95% AC range
p$stmv_distance_interpolation = p$stmv_distance_statsgrid * c( 1/2, 1, 2, 3, 4, 5, 10, 20, 40)  # range of permissible predictions km (i.e 1/2 stats grid to upper limit) .. in this case 5, 10, 20
p$stmv_distance_interpolate_predictions = p$stmv_distance_statsgrid * c( 1/2, 1, 2, 3, 4, 8) # finalizing preds using linear interpolation

use_parallel_mode = FALSE
if (use_parallel_mode) {
    # default is serial mode .. to enable parallel processing, pick and choose:
    scale_ncpus = ram_local( "ncores", ram_main=10, ram_process=4 ) # in GB; about 24  hr
    interpolate_ncpus = ram_local( "ncores", ram_main=2, ram_process=2 ) # nn hrs

    if (!exists("stmv_runmode", p)) p$stmv_runmode = list()

    p$stmv_runmode$globalmodel = TRUE

    p$stmv_runmode$scale = rep("localhost", scale_ncpus)

    p$stmv_runmode$interpolate = list(
      cor_0.25 = rep("localhost", interpolate_ncpus),
      cor_0.1  = rep("localhost", interpolate_ncpus),
      cor_0.05 = rep("localhost", interpolate_ncpus),
      cor_0.01 = rep("localhost", interpolate_ncpus)
    )
    p$stmv_runmode$interpolate_predictions = list(
      c1 = rep("localhost", interpolate_ncpus),
      c2 = rep("localhost", interpolate_ncpus),
      c3 = rep("localhost", interpolate_ncpus),
      c4 = rep("localhost", interpolate_ncpus),
      c5 = rep("localhost", interpolate_ncpus),
      c6 = rep("localhost", interpolate_ncpus),
      c7 = rep("localhost", interpolate_ncpus)
    )

    p$stmv_runmode$save_intermediate_results = TRUE
    p$stmv_runmode$save_completed_data = TRUE

}


stmv( p=p )


# quick look of data
DATA = substrate_db( p=p, DS="stmv_inputs" )
dev.new(); surface( as.image( Z=DATA$input$substrate.grainsize, x=DATA$input[, c("plon", "plat")], nx=p$nplons, ny=p$nplats, na.rm=TRUE) )

predictions = stmv_db( p=p, DS="stmv.prediction", ret="mean" )
statistics  = stmv_db( p=p, DS="stmv.stats" )

# locations = DATA$output$LOCS # these are the prediction locations
locations = bathymetry_db(spatial_domain=p$spatial_domain, DS="baseline")

# comparisons
dev.new(); surface( as.image( Z=log(predictions), x=locations, nx=p$nplons, ny=p$nplats, na.rm=TRUE) )

# stats
# statsvars = c( "sdTotal", "rsquared", "ndata", "sdSpatial", "sdObs", "phi", "nu", "localrange" )
statsvars = dimnames(statistics)[[2]]

dev.new(); levelplot( (predictions) ~ locations[,1] + locations[,2], aspect="iso" )
dev.new(); levelplot( statistics[,match("nu", statsvars)]  ~ locations[,1] + locations[,2], aspect="iso" ) # nu
dev.new(); levelplot( statistics[,match("sdTot", statsvars)]  ~ locations[,1] + locations[,2], aspect="iso" ) #sd total
dev.new(); levelplot( statistics[,match("localrange", statsvars)]  ~ locations[,1] + locations[,2], aspect="iso" ) #localrange


# as the interpolation process is so expensive, regrid based off the above run
substrate_db( p=p, DS="complete.redo" )


# quick map
b = bathymetry_db(spatial_domain=p$spatial_domain, DS="baseline")
o = substrate_db( p=p, DS="complete" )
lattice::levelplot( log(o$substrate.grainsize) ~ plon +plat, data=b, aspect="iso")


# or a cleaner map:
# p = aegis_parameters()
substrate_figures( p=p, varnames=c( "s.ndata", "s.sdTotal", "s.sdSpatial", "s.sdObs" ), logyvar=FALSE, savetofile="png" )
substrate_figures( p=p, varnames=c( "substrate.grainsize", "s.localrange", "s.nu", "s.phi"), logyvar=TRUE, savetofile="png" )


# to summarize just the global model
o = stmv_global_model( p=p, DS="global_model" )
summary(o)
plot(o)
AIC(o)  # [1]  3263839.33


# Global model results:
Family: gaussian
Link function: log
Family: gaussian
Link function: log

Formula:
substrate.grainsize ~ s(b.sdSpatial, k = 3, bs = "ts") + s(b.localrange,
    k = 3, bs = "ts") + s(log(z), k = 3, bs = "ts") + s(log(dZ),
    k = 3, bs = "ts") + s(log(ddZ), k = 3, bs = "ts")

Parametric coefficients:
              Estimate Std. Error  t value   Pr(>|t|)
(Intercept) -0.9053906  0.0111006 -81.5626 < 2.22e-16

Approximate significance of smooth terms:
                    edf Ref.df           F    p-value
s(b.sdSpatial)  1.99628      2   304.92667 < 2.22e-16
s(b.localrange) 1.99992      2  4461.65359 < 2.22e-16
s(log(z))       1.99666      2 18144.34273 < 2.22e-16
s(log(dZ))      1.98167      2    26.51703 2.2327e-12
s(log(ddZ))     1.97150      2     7.92948 0.00031715

R-sq.(adj) =   0.14   Deviance explained = 13.8%
GCV = 5.6951  Scale est. = 5.695     n = 713021

