

substrate_carstm = function( p=NULL, DS="aggregated_data", id=NULL, sppoly=NULL, redo=FALSE, ... ) {

  #\\ Note inverted convention: depths are positive valued
  #\\ i.e., negative valued for above sea level and positive valued for below sea level
  if ( is.null(p)) p = substrate_parameters(...)

  if ( !exists("project_name", p)) p$project_name = "substrate"
  if ( !exists("data_root", p) ) p$data_root = project.datadirectory( "aegis", p$project_name )
  if ( !exists("datadir", p) )   p$datadir  = file.path( p$data_root, "data" )
  if ( !exists("modeldir", p) )  p$modeldir = file.path( p$data_root, "modelled" )



  if (!exists("areal_units_strata_type", p )) p$areal_units_strata_type = "lattice" #
  if (!exists("areal_units_constraint", p )) p$areal_units_constraint = "none" #
  if (!exists("areal_units_overlay", p )) p$areal_units_overlay = "none" #
  if (!exists("areal_units_resolution_km", p )) stop( "areal_units_resolution_km should be defined ... " ) # km
  if (!exists("areal_units_proj4string_planar_km", p )) stop( "areal_units_proj4string_planar_km should be defined ... " ) # km
  if (!exists("timeperiod", p) )  p$timeperiod="default"



  if (is.null(id)) id = paste( p$spatial_domain, p$areal_units_overlay, p$areal_units_resolution_km, p$areal_units_strata_type, sep="_" )



  # -----------------------

  if ( DS=="aggregated_data") {

    fn = file.path( p$modeldir, paste( "bathymetry", "aggregated_data", id, "rdata", sep=".") )
    if (!redo)  {
      print( "Warning: aggregated_data is loading from a saved instance ... add redo=TRUE if data needs to be refresh" )
      if (file.exists(fn)) {
        load( fn)
        return( M )
      }
      print( "Warning: aggregated_data load from saved instance failed ... " )
    }

    print( "Warning: aggregated_data is being recreated ... " )
    print( "Warning: this needs a lot of RAM .. ~XX GB depending upon resolution of discretization .. a few hours " )

    M = substrate.db( p=p, DS="lonlat.highres" )
    M$substrate.grainsize = M$grainsize

    if (!exists("inputdata_spatial_discretization_planar_km", p) )  p$inputdata_spatial_discretization_planar_km = 1

    # thin data a bit ... remove potential duplicates and robustify
    M = lonlat2planar( M, proj.type=p$aegis_proj4string_planar_km )
    M$plon = round(M$plon / p$inputdata_spatial_discretization_planar_km + 1 ) * p$inputdata_spatial_discretization_planar_km
    M$plat = round(M$plat / p$inputdata_spatial_discretization_planar_km + 1 ) * p$inputdata_spatial_discretization_planar_km
    gc()

    bb = as.data.frame( t( simplify2array(
      tapply( X=M$z, INDEX=list(paste( paste( M$plon, M$plat) ) ),
        FUN = function(w) { c(
          mean(w, na.rm=TRUE),
          sd(w, na.rm=TRUE),
          length( which(is.finite(w)) )
        ) }, simplify=TRUE )
    )))
    M = NULL
    colnames(bb) = c("substrate.grainsize.mean", "substrate.grainsize.sd", "substrate.grainsize.n")
    plonplat = matrix( as.numeric( unlist(strsplit( rownames(bb), " ", fixed=TRUE))), ncol=2, byrow=TRUE)

    bb$plon = plonplat[,1]
    bb$plat = plonplat[,2]
    plonplat = NULL

    M = bb[ which( is.finite( bb$z.mean )) ,]
    bb =NULL
    gc()
    M = planar2lonlat( M, p$aegis_proj4string_planar_km)
    save(M, file=fn, compress=TRUE)

    return( M )
  }


  # ----------------------


  if ( DS=="carstm_inputs") {

    fn = file.path( p$modeldir, paste( "substrate", "carstm_inputs", id, "rdata", sep=".") )
    if (!redo)  {
      if (file.exists(fn)) {
        load( fn)
        return( S )
      }
    }

    # prediction surface
    sppoly = areal_units( p=p )  # will redo if not found
    sppoly = sppoly["StrataID"]

    crs_lonlat = sp::CRS(projection_proj4string("lonlat_wgs84"))

    # do this immediately to reduce storage for sppoly (before adding other variables)

    M = substrate.carstm ( p=p, DS="aggregated_data" )  # 16 GB in RAM just to store!

    # reduce size
    M = M[ which( M$lon > p$corners$lon[1] & M$lon < p$corners$lon[2]  & M$lat > p$corners$lat[1] & M$lat < p$corners$lat[2] ), ]
    # levelplot(z.mean~plon+plat, data=M, aspect="iso")

    M$StrataID = over( SpatialPoints( M[, c("lon", "lat")], crs_lonlat ), spTransform(sppoly, crs_lonlat ) )$StrataID # match each datum to an area
    M$lon = NULL
    M$lat = NULL
    M$plon = NULL
    M$plat = NULL
    M = M[ which(is.finite(M$StrataID)),]
    M$StrataID = as.character( M$StrataID )  # match each datum to an area
    M$tag = "observations"
    M$substrate.grainsize[!is.finite(M$substrate.grainsize)] = median(M$substrate.grainsize, na.rm=TRUE )  # missing data .. quick fix .. do something better


    B = bathymetry.carstm ( p=p, DS="aggregated_data" )  # unmodeled!
    # reduce size
    B = B[ which( B$lon > p$corners$lon[1] & B$lon < p$corners$lon[2]  & B$lat > p$corners$lat[1] & B$lat < p$corners$lat[2] ), ]
    locsmap = match(
      stmv::array_map( "xy->1", M[, c("plon","plat")], gridparams=p$gridparams ),
      stmv::array_map( "xy->1", B[, c("plon","plat")], gridparams=p$gridparams ) )
    M$z.mean = B$z.mean[locsmap]
    M$z.sd = B$z.sd[locsmap]
    M = M[ which(is.finite(M$z.mean)), ]

    B = NULL


    sppoly_df = as.data.frame(sppoly)
    sppoly_df$z = NA
    sppoly_df$StrataID = as.character( sppoly_df$StrataID )
    sppoly_df$tag ="predictions"

    vn = c("substrate.grainsize.mean", "z.mean", "tag", "StrataID")

    M = rbind( M[, vn], sppoly_df[, vn] )
    sppoly_df = NULL

    M$StrataID  = factor( as.character(M$StrataID), levels=levels( sppoly$StrataID ) ) # revert to factors
    sppoly = NULL
    M$strata  = as.numeric( M$StrataID)
    M$iid_error = 1:nrow(M) # for inla indexing for set level variation


    # prediction surface
    sppoly = areal_units( p=p )  # will redo if not found
    sppoly = sppoly["StrataID"]

    # M$StrataID  = as.character(M$StrataID)
    # M$tag = "observations"
    # M$Y = M$z

    ddepths = c(2.5, 5, 10, 20, 40, 80, 160, 320, 640 )  # depth cut points
    M$zi = as.numeric( as.character( cut( M$z, breaks=ddepths, labels=diff(ddepths)/2 + ddepths[-length(ddepths)], include.lowest=TRUE ) ))

    # do this immediately to reduce storage for sppoly (before adding other variables)
    M$StrataID = as.character( over( SpatialPoints( M[, c("lon", "lat")], crs_lonlat ), spTransform(sppoly, crs_lonlat ) )$StrataID) # match each datum to an area
    M$lon = NULL
    M$lat = NULL

    M$tag = "observations"

    sppoly = as.data.frame(sppoly)
    sppoly$z = NA
    sppoly$StrataID = as.character( sppoly$StrataID )
    sppoly$tag ="predictions"

    sppoly_depths = bathymetry.carstm( p=p, DS="" )

    M = rbind( M, sppoly[, names(M)] )

    M$StrataID  = factor( as.character(M$StrataID), levels=levels( sppoly$StrataID ) ) # revert to factors
    M$strata  = as.numeric( M$StrataID)
    M$iid_error = 1:nrow(M) # for inla indexing for set level variation

    gc()

    save( M, file=fn, compress=TRUE )
    return( M )
  }

  # -----------------------



  if ( DS %in% c("carstm_modelled", "carstm_modelled_fit") ) {

    fn = file.path( p$modeldir, paste( "substrate", "carstm_modelled", id, p$carstm_modelengine, p$carstm_family, "rdata", sep=".") )
    fn_fit = file.path( p$modeldir, paste( "substrate", "carstm_modelled_fit", id, p$carstm_modelengine, "rdata", sep=".") )

    if (!redo)  {
      print( "Warning: carstm_modelled is loading from a saved instance ... add redo=TRUE if data needs to be refresh" )
      if (DS=="carstm_modelled") {
        if (file.exists(fn)) {
          load( fn)
          return( sppoly )
        }
      }
      if (DS=="carstm_modelled_fit") {
        if (file.exists(fn_fit)) {
          load( fn_fit )
          return( fit )
        }
      }
      print( "Warning: carstm_modelled load from saved instance failed ... " )
    }

    print( "Warning: carstm_modelled is being recreated ... " )
    print( "Warning: this needs a lot of RAM .. ~XX GB depending upon resolution of discretization .. a few hours " )



    # prediction surface
    sppoly = areal_units( p=p )  # will redo if not found
#    sppoly = sppoly["StrataID"]

    M = substrate_carstm( p=p, DS="carstm_inputs" )  # will redo if not found
    fit  = NULL

    if ( grepl("glm", p$carstm_modelengine) ) {

      assign("fit", eval(parse(text=paste( "try(", p$carstm_modelcall, ")" ) ) ))
      if (is.null(fit)) error("model fit error")
      if ("try-error" %in% class(fit) ) error("model fit error")
      save( fit, file=fn_fit, compress=TRUE )

      # s = summary(fit)
      # AIC(fit)  # 104487274
      # reformat predictions into matrix form
      ii = which( M$tag=="predictions" & M$StrataID %in% M[ which(M$tag=="observations"), "StrataID"] )
      preds = predict( fit, newdata=M[ii,], type="link", na.action=na.omit, se.fit=TRUE )  # no/km2

      # out = reformat_to_matrix(
      #   input = preds$fit,
      #   matchfrom = list( StrataID=M$StrataID[ii] ),
      #   matchto   = list( StrataID=sppoly$StrataID  )
      # )
      # iy = match( as.character(sppoly$StrataID), aps$StrataID )
        sppoly@data[,"z.predicted"] = exp( preds$fit) - p$constant_offset
        sppoly@data[,"z.predicted_se"] = exp( preds$se.fit)
        sppoly@data[,"z.predicted_lb"] = exp( preds$fit - preds$se.fit ) - p$constant_offset
        sppoly@data[,"z.predicted_ub"] = exp( preds$fit + preds$se.fit ) - p$constant_offset
        save( sppoly, file=fn, compress=TRUE )
      }

      if ( grepl("gam", p$carstm_modelengine) ) {
        assign("fit", eval(parse(text=paste( "try(", p$carstm_modelcall, ")" ) ) ))
        if (is.null(fit)) error("model fit error")
        if ("try-error" %in% class(fit) ) error("model fit error")
        save( fit, file=fn_fit, compress=TRUE )

        s = summary(fit)
        AIC(fit)  # 104487274
        # reformat predictions into matrix form
        ii = which( M$tag=="predictions" & M$StrataID %in% M[ which(M$tag=="observations"), "StrataID"] )
        preds = predict( fit, newdata=M[ii,], type="link", na.action=na.omit, se.fit=TRUE )  # no/km2
        sppoly@data[,"z.predicted"] = exp( preds$fit) - p$constant_offset
        sppoly@data[,"z.predicted_se"] = exp( preds$se.fit)
        sppoly@data[,"z.predicted_lb"] = exp( preds$fit - preds$se.fit ) - p$constant_offset
        sppoly@data[,"z.predicted_ub"] = exp( preds$fit + preds$se.fit ) - p$constant_offset
        save( sppoly, file=fn, compress=TRUE )
    }


    if ( grepl("inla", p$carstm_modelengine) ) {

      H = carstm_hyperparameters( sd(log(M$z), na.rm=TRUE), alpha=0.5, median( log(M$z), na.rm=TRUE) )
      assign("fit", eval(parse(text=paste( "try(", p$carstm_modelcall, ")" ) ) ))
      if (is.null(fit)) error("model fit error")
      if ("try-error" %in% class(fit) ) error("model fit error")
      save( fit, file=fn_fit, compress=TRUE )

      s = summary(fit)
      s$dic$dic  # 31225
      s$dic$p.eff # 5200

      plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )

      # reformat predictions into matrix form
      ii = which(M$tag=="predictions")
      jj = match(M$StrataID[ii], sppoly$StrataID)
      sppoly@data$z.predicted = exp( fit$summary.fitted.values[ ii[jj], "mean" ]) - p$constant_offset
      sppoly@data$z.predicted_lb = exp( fit$summary.fitted.values[ ii[jj], "0.025quant" ]) - p$constant_offset
      sppoly@data$z.predicted_ub = exp( fit$summary.fitted.values[ ii[jj], "0.975quant" ]) - p$constant_offset
      sppoly@data$z.random_strata_nonspatial = exp( fit$summary.random$strata[ jj, "mean" ])
      sppoly@data$z.random_strata_spatial = exp( fit$summary.random$strata[ jj+max(jj), "mean" ])
      sppoly@data$z.random_sample_iid = exp( fit$summary.random$iid_error[ ii[jj], "mean" ])
      save( sppoly, file=fn, compress=TRUE )
    }


    if (map) {
      vn = "z.predicted"
      brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
      dev.new();  spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )
    }

    return( sppoly )

  }

}