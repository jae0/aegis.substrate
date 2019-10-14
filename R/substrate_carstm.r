

substrate_carstm = function( p=NULL, DS=NULL, redo=FALSE, ... ) {

  #\\ Note inverted convention: depths are positive valued
  #\\ i.e., negative valued for above sea level and positive valued for below sea level

  # ---------------------
  # deal with additional passed parameters
  if ( is.null(p) ) p=list()
  p_add = list(...)
  if (length(p_add) > 0 ) p = c(p, p_add)
  i = which(duplicated(names(p), fromLast = TRUE ))
  if ( length(i) > 0 ) p = p[-i] # give any passed parameters a higher priority, overwriting pre-existing variable



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
    M = substrate.db ( p=p, DS="aggregated_data" )  # 16 GB in RAM just to store!
    names(M)[which(names(M)=="substrate.grainsize.mean" )] = p$variabletomodel

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

    DS="unified"
    pb = aegis.bathymetry::bathymetry_parameters( p=p, project_class="carstm_auid" ) # transcribes relevant parts of p to load bathymetry
    BI = bathymetry_carstm ( p=pb, DS="carstm_inputs" )  # unmodeled!
    jj = match( as.character( M$StrataID), as.character( BI$StrataID) )
    M$z = BI$z[jj]
    jj =NULL

    M = M[ which(is.finite(M$z)), ]

    BI = NULL

    sppoly_df = as.data.frame(sppoly)
    BM = bathymetry_carstm ( p=pb, DS="carstm_modelled" )  # modeled!
    kk = match( as.character(  sppoly_df$StrataID), as.character( BM$StrataID ) )
    sppoly_df$z = BM$z.predicted[kk]
    BM = NULL
    sppoly_df[p$variabletomodel] = NA
    sppoly_df$StrataID = as.character( sppoly_df$StrataID )
    sppoly_df$tag ="predictions"

    vn = c(p$variabletomodel, "z", "tag", "StrataID")

    M = rbind( M[, vn], sppoly_df[, vn] )
    sppoly_df = NULL

    M$StrataID  = factor( as.character(M$StrataID), levels=levels( sppoly$StrataID ) ) # revert to factors

    save( M, file=fn, compress=TRUE )
    return( M )
  }

  # -----------------------



  if ( DS %in% c("carstm_modelled", "carstm_modelled_fit") ) {

    auids = paste(  p$auid, p$inputdata_spatial_discretization_planar_km, sep="_" )

    fn = file.path( p$modeldir, paste("substrate", "carstm_modelled", p$carstm_modelengine, auids, "rdata", sep="." ) )
    fn_fit = file.path( p$modeldir, paste( "substrate", "carstm_modelled_fit", p$carstm_modelengine, auids,  "rdata", sep=".") )

    if (!redo)  {
      if (DS=="carstm_modelled") {
        if (file.exists(fn)) {
          load( fn)
          return( res )
        }
      }
      if (DS=="carstm_modelled_fit") {
        if (file.exists(fn_fit)) {
          load( fn_fit )
          return( fit )
        }
      }
    }

    # prediction surface
    sppoly = areal_units( p=p )  # will redo if not found
    res = list(StrataID = sppoly[["StrataID"]])  # init results list
    res$strata = as.numeric(res$StrataID)

    M = substrate_carstm( p=p, DS="carstm_inputs" )  # will redo if not found
    M$strata  = as.numeric( M$StrataID)

    fit  = NULL

    if ( grepl("glm", p$carstm_modelengine) ) {

      assign("fit", eval(parse(text=paste( "try(", p$carstm_modelcall, ")" ) ) ))
      if (is.null(fit)) warning("model fit error")
      if ("try-error" %in% class(fit) ) warning("model fit error")
      save( fit, file=fn_fit, compress=TRUE )
      ii = which( M$tag=="predictions" & M$strata %in% M[ which(M$tag=="observations"), "strata"] )
      jj = match( M$strata[ii], res$strata)
      preds = predict( fit, newdata=M[ii,], type="link", na.action=na.omit, se.fit=TRUE )  # no/km2

      res[,paste(p$variabletomodel,"predicted", sep=".")] = exp( preds$fit[jj])
      res[,paste(p$variabletomodel,"predicted_se", sep=".")] = exp( preds$se.fit[jj])
      res[,paste(p$variabletomodel,"predicted_lb", sep=".")] = exp( preds$fit[jj] - preds$se.fit[jj] )
      res[,paste(p$variabletomodel,"predicted_ub", sep=".")] = exp( preds$fit[jj] + preds$se.fit[jj] )
      save( res, file=fn, compress=TRUE )
    }

    if ( grepl("gam", p$carstm_modelengine) ) {
      assign("fit", eval(parse(text=paste( "try(", p$carstm_modelcall, ")" ) ) ))
      if (is.null(fit)) warning("model fit error")
      if ("try-error" %in% class(fit) ) warning("model fit error")
      save( fit, file=fn_fit, compress=TRUE )
      ii = which( M$tag=="predictions" & M$strata %in% M[ which(M$tag=="observations"), "strata"] )
      jj = match( M$strata[ii], res$strata)
      preds = predict( fit, newdata=M[ii,], type="link", na.action=na.omit, se.fit=TRUE )  # no/km2
      res[,paste(p$variabletomodel,"predicted", sep=".")] = exp( preds$fit[jj] )
      res[,paste(p$variabletomodel,"predicted_se", sep=".")] = exp( preds$se.fit[jj])
      res[,paste(p$variabletomodel,"predicted_lb", sep=".")] = exp( preds$fit[jj] - preds$se.fit[jj] )
      res[,paste(p$variabletomodel,"predicted_ub", sep=".")] = exp( preds$fit[jj] + preds$se.fit[jj] )
      save( res, file=fn, compress=TRUE )
    }

    if ( grepl("inla", p$carstm_modelengine) ) {
      H = carstm_hyperparameters( sd(log(M[,p$variabletomodel]), na.rm=TRUE), alpha=0.5, median( log(M[,p$variabletomodel]), na.rm=TRUE) )
      M$zi = discretize_data( M$z, p$discretization$z )
      M$iid_error = 1:nrow(M) # for inla indexing for set level variation
      assign("fit", eval(parse(text=paste( "try(", p$carstm_modelcall, ")" ) ) ))
      if (is.null(fit)) warning("model fit error")
      if ("try-error" %in% class(fit) ) warning("model fit error")
      save( fit, file=fn_fit, compress=TRUE )
      ii = which(M$tag=="predictions")
      jj = match(M$strata[ii], res$strata)
      res[paste( p$variabletomodel, "predicted", sep=".")] = exp( fit$summary.fitted.values[ ii[jj], "mean" ])
      res[paste( p$variabletomodel, "predicted_lb", sep=".")] = exp( fit$summary.fitted.values[ ii[jj], "0.025quant" ])
      res[paste( p$variabletomodel, "predicted_ub", sep=".")] = exp( fit$summary.fitted.values[ ii[jj], "0.975quant" ])

      # simple spatial so just do the following here
      res[paste( p$variabletomodel, "random_strata_nonspatial", sep=".")] = exp( fit$summary.random$strata[ jj, "mean" ])
      res[paste( p$variabletomodel, "random_strata_spatial", sep=".")] = exp( fit$summary.random$strata[ jj+max(jj), "mean" ])
      res[paste( p$variabletomodel, "random_sample_iid", sep=".")] = exp( fit$summary.random$iid_error[ ii[jj], "mean" ])
      save( res, file=fn, compress=TRUE )
    }
    return( res )
  }

}
