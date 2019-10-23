

  substrate_carstm = function( p=NULL, DS=NULL, varnames=NULL, redo=FALSE ) {

    if ( is.null(p)) p = substrate_parameters(...)

    if ( !exists("project_name", p)) p$project_name = "substrate"
    if ( !exists("data_root", p) ) p$data_root = project.datadirectory( "aegis", p$project_name )
    if ( !exists("datadir", p) )   p$datadir  = file.path( p$data_root, "data" )
    if ( !exists("modeldir", p) )  p$modeldir = file.path( p$data_root, "modelled" )


  # ------------------


  if (DS =="carstm_auid") {
    # translate param values from one project to a unified representation
    # must be first to catch p
    pc = substrate_carstm(
      DS = "parameters",
      project_class = "carstm", # defines which parameter class / set to load
      project_name = "substrate",
      variabletomodel = "substrate.grainsize",
      spatial_domain = p$spatial_domain,  # defines spatial area, currenty: "snowcrab" or "SSE"
      areal_units_overlay = p$areal_units_overlay, # currently: "snowcrab_managementareas",  "groundfish_strata" .. additional polygon layers for subsequent analysis for now ..
      areal_units_resolution_km = p$areal_units_resolution_km, # km dim of lattice ~ 1 hr
      areal_units_proj4string_planar_km = p$areal_units_proj4string_planar_km,  # coord system to use for areal estimation and gridding for carstm
      inputdata_spatial_discretization_planar_km = p$inputdata_spatial_discretization_planar_km,  # 1 km .. some thinning .. requires 32 GB RAM and limit of speed -- controls resolution of data prior to modelling to reduce data set and speed up modelling
      auid = p$auid
    )
    return(pc)
  }


  # ---------------------



  if (DS=="parameters") {

    p$libs = unique( c( p$libs, project.library ( "carstm" ) ) )


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
        p$libs = unique( c( p$libs, project.library ("INLA" ) ) )
        p$carstm_modelcall = paste('
          inla(
            formula =', p$variabletomodel, ' ~ 1
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
            # num.threads=4,
            blas.num.threads=8,
            verbose=TRUE
          ) ' )
      }
      if ( grepl("glm", p$carstm_modelengine) ) {
        p$carstm_modelcall = paste('glm( formula =', p$variabletomodel, '~ 1 + StrataID,  family = gaussian(link="log"), data= M[ which(M$tag=="observations"), ], family=gaussian(link="log")  )  ' )  # for modelengine='glm'
      }
      if ( grepl("gam", p$carstm_modelengine) ) {
        p$libs = unique( c( p$libs, project.library ( "mgcv" ) ) )
        p$carstm_modelcall = paste('gam( formula =', p$variabletomodel, '~ 1 + StrataID,  family = gaussian(link="log"), data= M[ which(M$tag=="observations"), ], family=gaussian(link="log")  ) ' ) # for modelengine='gam'
      }
    }

    p = carstm_parameters( p=p, DS="basic" )

    return(p)
  }



  # ------------------



    if ( DS=="aggregated_data") {

      p = substrate_carstm(
        DS = "parameters",
        variabletomodel = p$variabletomodel,
        inputdata_spatial_discretization_planar_km=p$inputdata_spatial_discretization_planar_km
      )

      fn = file.path( p$datadir, paste( "substrate", "aggregated_data", p$inputdata_spatial_discretization_planar_km, "rdata", sep=".") )
      if (!redo)  {
        if (file.exists(fn)) {
          load( fn)
          return( M )
        }
      }

      M = substrate.db( p=p, DS="lonlat.highres" )
      M[,p$variabletomodel] = M$grainsize

      M$plon = round(M$plon / p$inputdata_spatial_discretization_planar_km + 1 ) * p$inputdata_spatial_discretization_planar_km
      M$plat = round(M$plat / p$inputdata_spatial_discretization_planar_km + 1 ) * p$inputdata_spatial_discretization_planar_km

      bb = as.data.frame( t( simplify2array(
        tapply( X=M[,p$variabletomodel], INDEX=list(paste(  M$plon, M$plat ) ),
          FUN = function(w) { c(
            mean(w, na.rm=TRUE),
            sd(w, na.rm=TRUE),
            length( which(is.finite(w)) )
          ) }, simplify=TRUE )
      )))
      M = NULL
      colnames(bb) = paste( p$variabletomodel, c("mean", "sd", "n"), sep=".")
      plonplat = matrix( as.numeric( unlist(strsplit( rownames(bb), " ", fixed=TRUE))), ncol=2, byrow=TRUE)

      bb$plon = plonplat[,1]
      bb$plat = plonplat[,2]
      plonplat = NULL

      M = bb[ which( is.finite( bb[, paste(p$variabletomodel, "mean", sep=".") ] )) ,]
      bb =NULL
      gc()
      M = planar2lonlat( M, p$aegis_proj4string_planar_km)
      save(M, file=fn, compress=TRUE)

      return( M )
    }


    # ---------------------------------------



  if ( DS=="carstm_inputs") {

    fn = file.path( p$modeldir, paste( "substrate", "carstm_inputs", p$auid,
      p$inputdata_spatial_discretization_planar_km,
      "rdata", sep=".") )

    if (!redo)  {
      if (file.exists(fn)) {
        load( fn)
        return( M )
      }
    }
    message( "Generating carstm_inputs ... ")

    # prediction surface
    sppoly = areal_units( p=p )  # will redo if not found

    # do this immediately to reduce storage for sppoly (before adding other variables)
    M = substrate_carstm ( p=p, DS="aggregated_data" )  # 16 GB in RAM just to store!
    names(M)[which(names(M)==paste(p$variabletomodel, "mean", sep=".") )] = p$variabletomodel

    # reduce size
    M = M[ which( M$lon > p$corners$lon[1] & M$lon < p$corners$lon[2]  & M$lat > p$corners$lat[1] & M$lat < p$corners$lat[2] ), ]
    # levelplot(substrate.grainsize.mean~plon+plat, data=M, aspect="iso")

    crs_lonlat = sp::CRS(projection_proj4string("lonlat_wgs84"))
    M$StrataID = over( SpatialPoints( M[, c("lon", "lat")], crs_lonlat ), spTransform(sppoly, crs_lonlat ) )$StrataID # match each datum to an area

    M$lon = NULL
    M$lat = NULL
    M$plon = NULL
    M$plat = NULL
    M = M[ which(is.finite(M$StrataID)),]
    M$tag = "observations"

    pB = aegis.bathymetry::bathymetry_carstm( p=p, DS="carstm_auid" ) # transcribes relevant parts of p to load bathymetry
    BI = bathymetry_carstm ( p=pB, DS="carstm_inputs" )  # unmodeled!
    jj = match( as.character( M$StrataID), as.character( BI$StrataID) )
    M$z = BI$z[jj]
    jj =NULL

    M = M[ which(is.finite(M[, pB$variabletomodel] )), ]

    BI = NULL

    sppoly_df = as.data.frame(sppoly)
    BM = carstm_model ( p=pB, DS="carstm_modelled" )  # modeled!
    kk = match( as.character(  sppoly_df$StrataID), as.character( BM$StrataID ) )
    sppoly_df[, pB$variabletomodel] = BM[ kk, paste(pB$variabletomodel, "predicted", sep=".") ]
    sppoly_df[,  p$variabletomodel] = NA
    BM = NULL
    sppoly_df$StrataID = as.character( sppoly_df$StrataID )
    sppoly_df$tag ="predictions"

    vn = c( p$variabletomodel, pB$variabletomodel, "tag", "StrataID")

    M = rbind( M[, vn], sppoly_df[, vn] )
    sppoly_df = NULL

    M$StrataID  = factor( as.character(M$StrataID), levels=levels( sppoly$StrataID ) ) # revert to factors
    M$strata  = as.numeric( M$StrataID)
    M$zi = discretize_data( M[, pB$variabletomodel], p$discretization[[pB$variabletomodel]] )
    M$iid_error = 1:nrow(M) # for inla indexing for set level variation

    save( M, file=fn, compress=TRUE )
    return( M )
  }

}
