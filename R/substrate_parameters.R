

substrate_parameters = function( p=NULL, project_name=NULL, project_class="default", ... ) {

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

  p$project_name = ifelse ( !is.null(project_name), project_name, "substrate" )

  if ( !exists("data_root", p) ) p$data_root = project.datadirectory( "aegis", p$project_name )
  if ( !exists("datadir", p) )   p$datadir  = file.path( p$data_root, "data" )
  if ( !exists("modeldir", p) )  p$modeldir = file.path( p$data_root, "modelled" )

  if ( !file.exists(p$datadir) ) dir.create( p$datadir, showWarnings=F, recursive=T )
  if ( !file.exists(p$modeldir) ) dir.create( p$modeldir, showWarnings=F, recursive=T )
  if (!exists("spatial_domain", p) ) p$spatial_domain = "canada.east.highres"
  if (!exists("spatial_domain_subareas", p)) p$spatial_domain_subareas = c( "canada.east", "SSE", "snowcrab", "SSE.mpa" )

  p = spatial_parameters( p=p)  # default (= only supported resolution of 0.2 km discretization)  .. do NOT change


  if (project_class=="default") {
    return(p)
  }

  if (project_class=="stmv") {
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

    if (!exists("aegis_dimensionality", p)) p$aegis_dimensionality="space"

    p = aegis_parameters( p=p, DS="stmv"  )
    return(p)
  }


  # ----------------------------------


  if (project_class=="carstm") {

    p$libs = c( p$libs, project.library ( "spatialreg", "INLA", "raster", "mgcv",  "carstm" ) )

    if ( !exists("project_name", p)) p$project_name = "substrate"

    p = aegis_parameters( p=p, DS="carstm" )

    # if ( !exists("spatial_domain", p)) p$spatial_domain = "snowcrab"  # defines spatial area, currenty: "snowcrab" or "SSE"
    if ( !exists("spatial_domain", p)) p$spatial_domain = "SSE"  # defines spatial area, currenty: "snowcrab" or "SSE"
    if ( !exists("areal_units_strata_type", p)) p$areal_units_strata_type = "lattice" # "stmv_lattice" to use ageis fields instead of carstm fields ... note variables are not the same

    if ( p$spatial_domain == "SSE" ) {
      if ( !exists("areal_units_overlay", p)) p$areal_units_overlay = "groundfish_strata" #.. additional polygon layers for subsequent analysis for now ..
      if ( !exists("areal_units_resolution_km", p)) p$areal_units_resolution_km = 25 # km dim of lattice ~ 1 hr
      if ( !exists("areal_units_proj4string_planar_km", p)) p$areal_units_proj4string_planar_km = projection_proj4string("utm20")  # coord system to use for areal estimation and gridding for carstm
      if ( !exists("inputdata_spatial_discretization_planar_km", p)) p$inputdata_spatial_discretization_planar_km = 1  # 1 km .. requires 32 GB RAM and limit of speed -- controls resolution of data prior to modelling to reduce data set and speed up modelling
    }

    if ( p$spatial_domain == "snowcrab" ) {
      if ( !exists("areal_units_overlay", p)) p$areal_units_overlay = "snowcrab_managementareas" # currently: "snowcrab_managementareas",  "groundfish_strata" .. additional polygon layers for subsequent analysis for now ..
      if ( !exists("areal_units_resolution_km", p)) p$areal_units_resolution_km = 25 # km dim of lattice ~ 1 hr
      if ( !exists("areal_units_proj4string_planar_km", p)) p$areal_units_proj4string_planar_km = projection_proj4string("utm20")  # coord system to use for areal estimation and gridding for carstm
      # if ( !exists("areal_units_proj4string_planar_km", p)) p$areal_units_proj4string_planar_km = projection_proj4string("omerc_nova_scotia")  # coord system to use for areal estimation and gridding for carstm
      if ( !exists("inputdata_spatial_discretization_planar_km", p)) p$inputdata_spatial_discretization_planar_km = 1  # 1 km .. requires 32 GB RAM and limit of speed -- controls resolution of data prior to modelling to reduce data set and speed up modelling
    }

    if ( !exists("carstm_modelengine", p)) p$carstm_modelengine = "inla.default"  # {model engine}.{label to use to store}

    if ( !exists("carstm_modelcall", p)) {
      if ( grepl("inla", p$carstm_modelengine) ) {
        p$libs = c( p$libs, RLibrary ( "INLA" ) )
        p$carstm_modelcall = '
          inla(
            formula = substrate.grainsize ~ 1
              + f(zi, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
              + f(strata, model="bym2", graph=sppoly@nb, scale.model=TRUE, constr=TRUE, hyper=H$bym2)
              + f(iid_error, model="iid", hyper=H$iid),
            family = "lognormal",
            data= M,
            control.compute=list(dic=TRUE, config=TRUE),
            control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
            control.predictor=list(compute=FALSE, link=1 ),
            control.fixed=H$fixed,  # priors for fixed effects, generic is ok
            # control.inla=list(int.strategy="eb") ,# to get empirical Bayes results much faster.
            # control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
            num.threads=4,
            blas.num.threads=4,
            verbose=TRUE
          ) '
      }
      if ( grepl("glm", p$carstm_modelengine) ) {
        p$carstm_modelcall = 'glm( formula = z ~ 1 + StrataID,  family = gaussian(link="log"), data= M[ which(M$tag=="observations"), ], family=gaussian(link="log")  ) '  # for modelengine='glm'
      }
      if ( grepl("gam", p$carstm_modelengine) ) {
        p$libs = c( p$libs, RLibrary ( "mgcv" ) )
        p$carstm_modelcall = 'gam( formula = z ~ 1 + StrataID,  family = gaussian(link="log"), data= M[ which(M$tag=="observations"), ], family=gaussian(link="log")  ) '  # for modelengine='gam'
      }
    }
    return(p)
  }

}
