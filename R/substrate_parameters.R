

substrate_parameters = function( p=list(), project_name="substrate", project_class="default", ... ) {


  p = parameters_add(p, list(...) ) # add passed args to parameter list, priority to args


  # ---------------------

  # create/update library list
  p$libs = unique( c( p$libs, RLibrary ( "colorspace",  "fields", "geosphere", "lubridate",  "lattice",
    "maps", "mapdata", "maptools", "parallel",  "rgdal", "rgeos",  "sp", "splancs", "GADMTools", "INLA" ) ) )
  p$libs = unique( c( p$libs, project.library ( "aegis", "aegis.bathymetry", "aegis.coastline", "aegis.polygons", "aegis.substrate" ) ) )


  p = parameters_add_without_overwriting( p, project_name = project_name )
  p = parameters_add_without_overwriting( p, data_root = project.datadirectory( "aegis", p$project_name ) )
  p = parameters_add_without_overwriting( p, datadir  = file.path( p$data_root, "data" ) )
  p = parameters_add_without_overwriting( p, modeldir = file.path( p$data_root, "modelled" ) )

  if ( !file.exists(p$datadir) ) dir.create( p$datadir, showWarnings=FALSE, recursive=TRUE )
  if ( !file.exists(p$modeldir) ) dir.create( p$modeldir, showWarnings=FALSE, recursive=TRUE )


  p = parameters_add_without_overwriting( p,
    variabletomodel = "substrate.grainsize",
    spatial_domain = "SSE",
    spatial_domain_subareas = c( "canada.east",  "SSE", "SSE.mpa" , "snowcrab"),  # this is for bathymetry_db, not stmv
    aegis_dimensionality="space"
  )

  p = spatial_parameters( p=p)  # default (= only supported resolution of 0.2 km discretization)  .. do NOT change

  p = parameters_add_without_overwriting( p,
      inputdata_spatial_discretization_planar_km = 1  #  controls resolution of data prior to modelling (km .. ie 20 linear units smaller than the final discretization pres)
  )


  # ---------------------

  if (project_class=="default") return(p)

  # ---------------------

  if (project_class=="carstm")  {
    # simple run of carstm. There are two types:
    #   one global, run directly from  polygons defined in aegis.bathymetry/inst/scripts/99.bathymetry.carstm.R
    #   and one that is called secondarily specific to a local project's polygons (eg. snow crab)
    p$libs = c( p$libs, project.library ( "carstm", "INLA"  ) )

    p = parameters_add_without_overwriting( p,
      areal_units_source = "lattice", # "stmv_fields" to use ageis fields instead of carstm fields ... note variables are not the same
      areal_units_resolution_km = 5, # default in case not provided ... 25 km dim of lattice ~ 1 hr; 5km = 79hrs; 2km = ?? hrs
      areal_units_proj4string_planar_km = aegis::projection_proj4string("utm20"),  # coord system to use for areal estimation and gridding for carstm
      # areal_units_proj4string_planar_km = projection_proj4string("omerc_nova_scotia")  # coord system to use for areal estimation and gridding for carstm
      areal_units_overlay = "none",
      carstm_modelengine = "inla",  # {model engine}.{label to use to store}
      carstm_model_label = "default",
      carstm_inputs_aggregated = FALSE
    )

    if ( !exists("carstm_model_call", p)  ) {
      if ( grepl("inla", p$carstm_modelengine) ) {
        p$carstm_model_call = paste(
         'inla(
            formula =', p$variabletomodel, ' ~ 1
              + f( inla.group(z, method="quantile", n=9),  model="rw2", scale.model=TRUE, hyper=H$rw2)
              + f(auid, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE, constr=TRUE, hyper=H$bym2),
            family = "lognormal",
            data= M,
            control.compute=list(dic=TRUE, waic=TRUE, cpo=TRUE, config=TRUE),
            control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
            control.predictor=list(compute=FALSE, link=1 ),
            control.fixed=H$fixed,  # priors for fixed effects, generic is ok
            control.inla = list(h=1e-4, tolerance=1e-9, cmin=0), # restart=3), # restart a few times in case posteriors are poorly defined
            verbose=TRUE
          ) ' )
      }
    }

    p = carstm_parameters( p=p )  # fill in anything missing with defaults and do some checks

    return(p)
  }


  # ---------------------

  if (project_class=="stmv") {
    p$libs = unique( c( p$libs, project.library ( "stmv" ) ) )

    if (!exists("stmv_variables", p)) p$stmv_variables = list()
    if (!exists("LOCS", p$stmv_variables)) p$stmv_variables$LOCS=c("plon", "plat")
    if (!exists("stmv_local_modelengine", p)) p$stmv_local_modelengine="fft"  # currently the perferred approach
    if (!exists("stmv_global_modelengine", p)) p$stmv_global_modelengine = "gam"
    if (!exists("stmv_global_modelformula", p)) p$stmv_global_modelformula = formula( paste(
      p$variabletomodel, ' ~  s( log(z), k=3, bs="ts") + s( log(dZ), k=3, bs="ts" ) + s( log(ddZ), k=3, bs="ts" ) ') )
    if (!exists("stmv_global_family", p)) p$stmv_global_family = gaussian(link="log")

    if (p$stmv_global_modelengine == "gam") {
      p$libs = unique( c( p$libs, RLibrary ("mgcv")) )
      if (!exists("stmv_gam_optimizer", p)) p$stmv_gam_optimizer = c("outer", "bfgs")
      if (!exists("stmv_local_modelformula", p)) p$stmv_local_modelformula = formula( paste(
        p$variabletomodel, ' ~ s(plon,k=3, bs="ts") + s(plat, k=3, bs="ts") + s(plon, plat, k=200, bs="ts")' ) )
    }

    # some tweaked options for substrate
    if ( p$stmv_local_modelengine %in% c("krige" )) {
      # nothing to do ... faster than gstat
    }

    if (p$stmv_local_modelengine == "gam") {
      p$libs = unique( c( p$libs, RLibrary ("mgcv")))
      if (!exists("stmv_local_model_distanceweighted", p)) p$stmv_local_model_distanceweighted = TRUE
      if (!exists("stmv_gam_optimizer", p)) p$stmv_gam_optimizer = c("outer", "bfgs")
    }

    if ( p$stmv_local_modelengine == "bayesx" ) {
      p$libs = unique( c( p$libs, RLibrary ("bayesx")) )
      if (!exists("stmv_local_modelformula", p)) p$stmv_local_modelformula = formula( paste(
        p$variabletomodel, ' ~ s(plon, bs="ps") + s(plat, bs="ps") + s(plon, plat, bs="te")') )  # more detail than "gs" .. "te" is preferred
      if (!exists("stmv_local_model_bayesxmethod", p)) p$stmv_local_model_bayesxmethod="MCMC"  # REML actually seems to be the same speed ... i.e., most of the time is spent in thhe prediction step ..
      if (!exists("stmv_local_model_distanceweighted", p)) p$stmv_local_model_distanceweighted = TRUE

    }

    if ( p$stmv_local_modelengine == "fft" ) {
      # definitely a cleaner (not overly smoothed) image than a GAM
      # NOTE that  p$stmv_lowpass_phi and  p$stmv_lowpass_nu are very critical choices
      # if (!exists("stmv_fft_filter", p)) p$stmv_fft_filter = "lowpass" # only act as a low pass filter .. depth has enough data for this. Otherwise, use:
      # if (!exists("stmv_fft_filter", p)) p$stmv_fft_filter = "matern" to ~ krige
      if (!exists("stmv_lowpass_phi", p)) p$stmv_lowpass_phi = p$pres*2 # FFT based method when operating gloablly
      if (!exists("stmv_lowpass_nu", p)) p$stmv_lowpass_nu = 0.5 # this is exponential covar
    }

    p = aegis_parameters( p=p, DS="stmv"  )
    return(p)
  }



}
