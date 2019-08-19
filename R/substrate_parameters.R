

substrate_parameters = function( p=NULL, project.name=NULL, project.mode="default", ... ) {

  # ---------------------
  # deal with additional passed parameters
  if ( is.null(p) ) p=list()
  p_add = list(...)
  if (length(p_add) > 0 ) p = c(p, p_add)
  i = which(duplicated(names(p), fromLast = TRUE ))
  if ( length(i) > 0 ) p = p[-i] # give any passed parameters a higher priority, overwriting pre-existing variable


  # ---------------------

  # create/update library list
  p$libs = c( p$libs, RLibrary ( "colorspace",  "fields", "geosphere", "lubridate",  "lattice",
    "maps", "mapdata", "maptools", "parallel",  "rgdal", "rgeos",  "sp", "splancs", "GADMTools" ) )
  p$libs = c( p$libs, project.library ( "aegis", "aegis.bathymetry", "aegis.coastline", "aegis.polygons", "aegis.substrate" ) )

  p$project.name = ifelse ( !is.null(project.name), project.name, "substrate" )

  if ( !exists("data_root", p) ) p$data_root = project.datadirectory( "aegis", p$project.name )
  if ( !exists("datadir", p) )   p$datadir  = file.path( p$data_root, "data" )
  if ( !exists("modeldir", p) )  p$modeldir = file.path( p$data_root, "modelled" )

  if ( !file.exists(p$datadir) ) dir.create( p$datadir, showWarnings=F, recursive=T )
  if ( !file.exists(p$modeldir) ) dir.create( p$modeldir, showWarnings=F, recursive=T )
  if (!exists("spatial.domain", p) ) p$spatial.domain = "canada.east.highres"
  if (!exists("spatial.domain.subareas", p)) p$spatial.domain.subareas = c( "canada.east", "SSE", "snowcrab", "SSE.mpa" )

  p = spatial_parameters( p=p)  # default (= only supported resolution of 0.2 km discretization)  .. do NOT change


  if (project.mode=="default") {
    return(p)
  }

  if (project.mode=="stmv") {
    p$libs = c( p$libs, project.library ( "stmv" ) )

    if (!exists("variables", p)) p$variables = list()
    if (!exists("LOCS", p$variables)) p$variables$LOCS=c("plon", "plat")
    if (!exists("stmv_local_modelengine", p)) p$stmv_local_modelengine="fft"  # currently the perferred approach
    if (!exists("stmv_global_modelengine", p)) p$stmv_global_modelengine = "gam"
    if (!exists("stmv_global_modelformula", p)) p$stmv_global_modelformula = formula( paste(
      'substrate.grainsize ~  s( log(z), k=3, bs="ts") + s( log(dZ), k=3, bs="ts" ) + s( log(ddZ), k=3, bs="ts" ) ') )
    if (!exists("stmv_global_family", p)) p$stmv_global_family = gaussian(link="log")

    if (p$stmv_global_modelengine == "gam") {
      p$libs = unique( c( p$libs, RLibrary ("mgcv")) )
      if (!exists("stmv_gam_optimizer", p)) p$stmv_gam_optimizer = c("outer", "bfgs")
      if (!exists("stmv_local_modelformula", p)) p$stmv_local_modelformula = formula( paste(
        'substrate.grainsize ~ s(plon,k=3, bs="ts") + s(plat, k=3, bs="ts") + s(plon, plat, k=200, bs="ts")' ) )
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
        'substrate.grainsize ~ s(plon, bs="ps") + s(plat, bs="ps") + s(plon, plat, bs="te")') )  # more detail than "gs" .. "te" is preferred
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

    if (!exists("stmv_dimensionality", p)) p$stmv_dimensionality="space"

    p = aegis_parameters( p=p, DS="stmv_spatial_model"  )
    return(p)
  }



  if (project.mode=="carstm") {
    p$libs = c( p$libs, project.library ( "carstm" ) )

    return(p)
  }

}
