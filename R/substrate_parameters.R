

substrate_parameters = function( p=list(), project_name="substrate", project_class="core", ... ) {

  p = parameters_add(p, list(...) ) # add passed args to parameter list, priority to args

  # ---------------------

  # create/update library list
#   p$libs = unique( c( p$libs, RLibrary ( "colorspace",  "fields", "geosphere", "lubridate",  "lattice",
#    "maps", "mapdata", "maptools", "parallel",  "rgdal", "rgeos",  "sp", "splancs", "GADMTools", "INLA" ) ) )
  p$libs = unique( c( p$libs, RLibrary ( "colorspace",  "lubridate", "lattice",
    "parallel",  "sf", "sp", "GADMTools", "INLA" ) ) )

  p$libs = unique( c( p$libs, project.library ( "aegis", "aegis.bathymetry", "aegis.coastline", "aegis.polygons", "aegis.substrate" ) ) )

  p = parameters_add_without_overwriting( p, project_name = project_name )
  p = parameters_add_without_overwriting( p, data_root = project.datadirectory( "aegis", p$project_name ) )
  p = parameters_add_without_overwriting( p, datadir  = file.path( p$data_root, "data" ) )
  p = parameters_add_without_overwriting( p, modeldir = file.path( p$data_root, "modelled" ) )

  if ( !file.exists(p$datadir) ) dir.create( p$datadir, showWarnings=FALSE, recursive=TRUE )
  if ( !file.exists(p$modeldir) ) dir.create( p$modeldir, showWarnings=FALSE, recursive=TRUE )


  p = parameters_add_without_overwriting( p,
    variabletomodel = "substrate.grainsize",
    spatial_domain = "canada.east.highres",
    spatial_domain_subareas = c( "canada.east",  "SSE", "SSE.mpa" , "snowcrab"),  # this is for bathymetry_db, not stmv
    aegis_dimensionality="space"
  )

  p = spatial_parameters( p=p)  # default (= only supported resolution of 0.2 km discretization)  .. do NOT change

  p = parameters_add_without_overwriting( p, inputdata_spatial_discretization_planar_km = p$pres/2 )
   #  controls resolution of data prior to modelling (km .. ie 20 linear units smaller than the final discretization pres)


  # ---------------------

  if (project_class=="core") {
    p$project_class="core"
    return(p)
  }

  # ---------------------

  if (project_class %in% c("carstm") ) {
    # simple run of carstm. There are two types:
    #   one global, run directly from  polygons defined in aegis.bathymetry/inst/scripts/99.bathymetry.carstm.R
    #   and one that is called secondarily specific to a local project's polygons (eg. snow crab)
    p$libs = c( p$libs, project.library ( "carstm", "INLA"  ) )
    p$project_class="carstm"

    p = parameters_add_without_overwriting( p,
      areal_units_xydata = "substrate_db(p=p, DS='areal_units_input')",
      areal_units_type = "lattice", # "stmv_fields" to use ageis fields instead of carstm fields ... note variables are not the same
      areal_units_resolution_km = 25, # default in case not provided ... 25 km dim of lattice ~ 1 hr; 5km = 79hrs; 2km = ?? hrs
      areal_units_proj4string_planar_km = p$aegis_proj4string_planar_km,  # coord system to use for areal estimation and gridding for carstm
      # areal_units_proj4string_planar_km = projection_proj4string("omerc_nova_scotia")  # coord system to use for areal estimation and gridding for carstm
      areal_units_overlay = "none",
      areal_units_timeperiod = "none",
      tus="none",
      fraction_todrop = 1/5,
      fraction_cv = 1.0,
      fraction_good_bad = 0.9,
      nAU_min = 30,
      carstm_modelengine = "inla",  # {model engine}.{label to use to store}
      carstm_model_label = "default",
      carstm_inputs_aggregated = TRUE
    )

    if ( !exists("carstm_inputdata_model_source", p))  p$carstm_inputdata_model_source = list()
    p$carstm_inputdata_model_source = parameters_add_without_overwriting( p$carstm_inputdata_model_source,
      bathymetry = "stmv"  # "stmv", "hybrid", "carstm"
    )


    if ( grepl("inla", p$carstm_modelengine) ) {
      if ( !exists("carstm_model_formula", p)  ) {
        p$carstm_model_formula = as.formula( paste(
         p$variabletomodel, ' ~ 1',
             '+ f( inla.group(z, method="quantile", n=9),  model="rw2", scale.model=TRUE, hyper=H$rw2)',
             '+ f(auid, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE, constr=TRUE, hyper=H$bym2)'
         ) )
      }
      if ( !exists("carstm_model_family", p)  )  p$carstm_model_family = "lognormal"
    }

    p = carstm_parameters( p=p )  # fill in anything missing with defaults and do some checks

    if ( p$inputdata_spatial_discretization_planar_km >= p$areal_units_resolution_km ) {
      warning( "p$inputdata_spatial_discretization_planar_km >= p$areal_units_resolution_km " )
    }
    message ("p$areal_units_resolution_km: ", p$areal_units_resolution_km)

    return(p)
  }


  # ---------------------

  if (project_class %in% c("stmv" ,"default") ) {

    p$libs = c( p$libs, project.library ( "stmv" ) )
    p$project_class="stmv"

    p = parameters_add_without_overwriting( p,
      DATA = 'substrate_db( p=p, DS="stmv_inputs" )',  # _highres
      stmv_model_label="default",
      stmv_variables = list(Y="substrate.grainsize", LOCS=c("plon", "plat")),  # required as fft has no formulae
      stmv_global_modelengine = "gam",  # only marginally useful .. consider removing it and use "none",
      stmv_local_modelengine="fft",
      stmv_variogram_method = "fft",
      stmv_filter_depth_m = 0.5,  # > 0.5 m
      stmv_rsquared_threshold = 0.01, # lower threshold  .. ignore
      stmv_distance_statsgrid = 5, # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
      stmv_nmin = 90, # min number of data points req before attempting to model in a localized space
      stmv_nmax = 1000, # no real upper bound.. just speed /RAM
      stmv_force_complete_method = "linear_interp"
    )

    p = parameters_add_without_overwriting( p,
      stmv_distance_prediction_limits = p$stmv_distance_statsgrid * c( 1/2, 5 ), # range of permissible predictions km (i.e 1/2 stats grid to upper limit based upon data density)
      stmv_distance_scale = p$stmv_distance_statsgrid * c( 1, 2, 3, 4, 5, 10, 20, 40), # km ... approx guesses of 95% AC range
      stmv_distance_interpolation = p$stmv_distance_statsgrid * c( 1/2, 1, 2, 3, 4, 5, 10, 20, 40),  # range of permissible predictions km (i.e 1/2 stats grid to upper limit) .. in this case 5, 10, 20
      stmv_distance_interpolate_predictions = p$stmv_distance_statsgrid * c( 1/2, 1, 2, 3, 4, 8) # finalizing preds using linear interpolation
    )


    # some tweaked options for substrate
    if ( p$stmv_local_modelengine %in% c("krige" )) {
      # nothing to do ... faster than gstat
    }

    if (p$stmv_global_modelengine == "gam") {
      p$libs = unique( c( p$libs, RLibrary ("mgcv")))
      p = parameters_add_without_overwriting( p,
        stmv_global_modelformula = formula( paste(
          p$variabletomodel, ' ~  1 + s( log(z), k=3, bs="ts") + s( log(dZ), k=3, bs="ts" ) + s( log(ddZ), k=3, bs="ts" )+ s( b.sdSpatial, k=3, bs="ts") + s( b.localrange, k=3, bs="ts")
          ') ),
        stmv_global_family = gaussian(link="log"),
        stmv_gam_optimizer = c("outer", "bfgs")
      )
    }

    if (p$stmv_local_modelengine == "gam") {
      p$libs = unique( c( p$libs, RLibrary ("mgcv")))
      p = parameters_add_without_overwriting( p,
        stmv_local_modelformula = formula( paste(
          p$variabletomodel, ' ~ s(plon,k=3, bs="ts") + s(plat, k=3, bs="ts") + s(plon, plat, k=200, bs="ts")') ),
        stmv_local_model_distanceweighted = TRUE,
        stmv_gam_optimizer = c("outer", "bfgs")
      )
    }

    if ( p$stmv_local_modelengine == "bayesx" ) {
      p$libs = unique( c( p$libs, RLibrary ("bayesx")) )
      if (!exists("stmv_local_modelformula", p)) p$stmv_local_modelformula = formula( paste(
        p$variabletomodel, ' ~ s(plon, bs="ps") + s(plat, bs="ps") + s(plon, plat, bs="te")') )  # more detail than "gs" .. "te" is preferred
      if (!exists("stmv_local_model_bayesxmethod", p)) p$stmv_local_model_bayesxmethod="MCMC"  # REML actually seems to be the same speed ... i.e., most of the time is spent in thhe prediction step ..
      if (!exists("stmv_local_model_distanceweighted", p)) p$stmv_local_model_distanceweighted = TRUE

    }

    if ( p$stmv_local_modelengine == "fft" ) {
      nu = 0.5  # exponential smoothing
      ac_local = 0.1  # ac at which to designate "effective range"
      p = parameters_add_without_overwriting( p,
        stmv_fft_filter = "matern tapered lowpass modelled fast_predictions", #  act as a low pass filter first before matern with taper
        # stmv_fft_filter = "matern_tapered_modelled",
        stmv_autocorrelation_fft_taper = 0.75,  # benchmark from which to taper
        stmv_autocorrelation_localrange = ac_local,  # for output to stats
        stmv_autocorrelation_interpolation = c(0.25, 0.1, 0.05, 0.01),
        stmv_lowpass_nu = nu, # exp
        stmv_lowpass_phi = stmv::matern_distance2phi( distance=p$pres/2, nu=nu, cor=ac_local )
      )
    }

    # default to serial mode
    p = parameters_add_without_overwriting( p,
      stmv_runmode = list(
        globalmodel = TRUE,
        scale = rep("localhost", 1),
        interpolate_correlation_basis = list(
          cor_0.25 = rep("localhost", 1),
          cor_0.1  = rep("localhost", 1),
          cor_0.05 = rep("localhost", 1),
          cor_0.01 = rep("localhost", 1)
        ),
        interpolate_predictions = list(
          c1 = rep("localhost", 1),
          c2 = rep("localhost", 1),
          c3 = rep("localhost", 1),
          c4 = rep("localhost", 1),
          c5 = rep("localhost", 1),
          c6 = rep("localhost", 1),
          c7 = rep("localhost", 1)
        ),
        save_intermediate_results = TRUE,
        save_completed_data = TRUE # just a dummy variable with the correct name
      )
    )

    p = aegis_parameters( p=p, DS="stmv" )  # get defaults

    if ( p$inputdata_spatial_discretization_planar_km >= p$pres ) {
      warning( "p$inputdata_spatial_discretization_planar_km >= p$pres " )
    }
  #    message ("p$stmv_distance_statsgrid: ", p$stmv_distance_statsgrid)

    return(p)
  }


  # ---------------------

  if (project_class %in% c("hybrid") ) {

    p$project_class="hybrid"
    p$libs = unique( c( p$libs, RLibrary ("carstm", "stmv", "INLA" )) )

    p = parameters_add_without_overwriting( p,
      DATA = 'substrate_db( p=p, DS="stmv_inputs" )',  # _highres
      stmv_variables = list(Y="substrate.grainsize", LOCS=c("plon", "plat")),  # required as fft has no formulae
      stmv_model_label="default",
      stmv_global_modelengine = "gam",  # only marginally useful .. consider removing it and use "none",
      stmv_global_modelformula = formula( paste(
        p$variabletomodel, ' ~  1
          + s( log(z), k=3, bs="ts") + s( log(dZ), k=3, bs="ts" ) + s( log(ddZ), k=3, bs="ts" )
          + s( b.sdSpatial, k=3, bs="ts") + s( b.localrange, k=3, bs="ts")
        ') ),
      stmv_global_family = gaussian(link="log"),
      stmv_local_modelengine="carstm",
      stmv_local_covariates_carstm = "",  # only model covariates globally
      stmv_local_all_carstm = "",  # ignoring au
      stmv_local_modelcall = paste(
        'inla(
          formula = z ~ 1
            + f(auid, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE, constr=TRUE, hyper=H$bym2),
          family = "normal",
          data= dat,
          control.compute=list(dic=TRUE, waic=TRUE, cpo=FALSE, config=FALSE),  # config=TRUE if doing posterior simulations
          control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
          control.predictor=list(compute=FALSE, link=1 ),
          control.fixed=H$fixed,  # priors for fixed effects, generic is ok
          verbose=FALSE
        ) '
      ),   # NOTE:: this is a local model call
      stmv_distance_statsgrid = 1, # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
      stmv_nmin = 100, # min number of data points req before attempting to model in a localized space
      stmv_nmax = 1000 # no real upper bound.. just speed /RAM
    )


    p = parameters_add_without_overwriting( p,
      stmv_distance_prediction_limits = p$stmv_distance_statsgrid * c( 1, 5 ), # range of permissible predictions km (i.e  stats grid to upper limit based upon data density)
      stmv_distance_interpolation = p$stmv_distance_statsgrid * c( 1/2, 1, 2, 3, 4, 5),  # range of permissible predictions km (i.e 1/2 stats grid to upper limit) .. in this case 5, 10, 20
      stmv_distance_interpolate_predictions = p$stmv_distance_statsgrid * c( 1/2, 1, 2, 3, 4, 5) # finalizing preds using linear interpolation
    )

    p = parameters_add_without_overwriting( p,
      stmv_runmode = list(
        carstm = rep("localhost", 1),
        globalmodel = FALSE,
        save_intermediate_results = TRUE,
        save_completed_data = TRUE
      )
    )

    p = aegis_parameters( p=p, DS="stmv" )  # get defaults


    if ( p$inputdata_spatial_discretization_planar_km >= p$pres ) {
      warning( "p$inputdata_spatial_discretization_planar_km >= p$pres " )
    }
    message ("p$stmv_distance_statsgrid: ", p$stmv_distance_statsgrid)
    return(p)

  }



}
