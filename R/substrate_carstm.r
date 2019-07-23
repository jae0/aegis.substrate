
  substrate_carstm = function( p=NULL, DS="aggregated_data", id=NULL, sppoly=NULL, redo=FALSE, ... ) {

    #\\ Note inverted convention: depths are positive valued
    #\\ i.e., negative valued for above sea level and positive valued for below sea level
    if ( is.null(p)) p = substrate_parameters(...)

    if ( !exists("project.name", p)) p$project.name = "substrate"
    if ( !exists("data_root", p) ) p$data_root = project.datadirectory( "aegis", p$project.name )
    if ( !exists("datadir", p) )   p$datadir  = file.path( p$data_root, "data" )
    if ( !exists("modeldir", p) )  p$modeldir = file.path( p$data_root, "modelled" )

    if (is.null(id)) id = paste( p$resolution, p$spatial.domain, p$overlay, sep="_" )

    if ( DS=="aggregated_data") {

      fn = file.path( p$modeldir, paste( "substrate", "aggregated_data", id, "rdata", sep=".") )
      if (!redo) {
        print( "Warning: aggregated_data is loading from a saved instance ... add redo=TRUE if data needs to be refresh" )
        if (file.exists(fn)) {
          load( fn)
          return( sppoly )
        }
        print( "Warning: aggregated_data load from saved instance failed ... " )
      }

      print( "Warning: aggregated_data is being recreated ... " )
      M = substrate.db ( p=p, DS="lonlat.highres" )  # 16 GB in RAM just to store!
      crs_lonlat = sp::CRS("+proj=longlat +datum=WGS84")
      spp = sppoly["StrataID"]
      M$StrataID = over( SpatialPoints( M[, c("lon", "lat")], crs_lonlat ), spTransform(spp, crs_lonlat ) )$StrataID # match each datum to an area
      o = NULL
      M$lon = NULL
      M$lat = NULL
      gc()
      bb = as.data.frame( t( simplify2array(
        tapply( X=M$z, INDEX=list(paste( M$StrataID) ),
          FUN = function(w) { c(
            mean(w, na.rm=TRUE),
            sd(w, na.rm=TRUE),
            length( which(is.finite(w)) )
          ) }, simplify=TRUE )
      )))
      M = NULL
      colnames(bb) = c("substrate.grainsize.mean", "substrate.grainsize.sd", "substrate.grainsize.n")
      bb$StrataID = rownames(bb)
      sppoly$substrate.grainsize = NA
      sppoly$substrate.grainsize.sd = NA
      sppoly$substrate.grainsize.n = NA
      j = match( bb$StrataID, sppoly$StrataID )
      if (length(j) > 0)  {
        sppoly$substrate.grainsize[j] = bb$substrate.grainsize.mean
        sppoly$substrate.grainsize.sd[j] = bb$substrate.grainsize.sd
        sppoly$substrate.grainsize.n[j] = bb$substrate.grainsize.n
      }
      save( sppoly, file=fn, compress=TRUE )
      return( sppoly )
    }




    if ( DS=="carstm_modelled") {

      fn = file.path( p$modeldir, paste( "substrate", "carstm_modelled", id, "rdata", sep=".") )
      fn_fit = file.path( p$modeldir, paste( "substrate", "carstm_modelled_fit", id, "rdata", sep=".") )
      if (!redo) {
        print( "Warning: carstm_modelled is loading from a saved instance ... add redo=TRUE if data needs to be refresh" )
        if (file.exists(fn)) {
          load( fn)
          return( sppoly )
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
      M = substrate.db ( p=p, DS="lonlat.highres" )  # 16 GB in RAM just to store!
      crs_lonlat = sp::CRS("+proj=longlat +datum=WGS84")
      spp = sppoly["StrataID"]

      # do this immediately to reduce storage for spp (before adding other variables)
      M$StrataID = as.character( over( SpatialPoints( M[, c("lon", "lat")], crs_lonlat ), spTransform(spp, crs_lonlat ) )$StrataID) # match each datum to an area
      o = NULL
      M$lon = NULL
      M$lat = NULL
      # M$substrate.grainsize
      M$tag = "observations"

      spp = as.data.frame(spp)
      spp$substrate.grainsize = NA
      spp$StrataID = as.character( spp$StrataID )
      spp$tag ="predictions"

      M = rbind( M, spp[, names(M)] )
      M$StrataID  = factor( as.character(M$StrataID), levels=levels( sppoly$StrataID ) ) # revert to factors
      M$strata  = as.numeric( M$StrataID)
      M$iid_error = 1:nrow(M) # for inla indexing for set level variation

      spp = NULL
      gc()

      B = aegis.bathymetry::bathymetry_carstm( p=p, DS="carstm_modelled" )
      ddepths = c(2.5, 5, 10, 20, 40, 80, 160, 320, 640 )
      M$z = as.numeric( as.character( cut( M$z, breaks=ddepths, labels=diff(ddepths)/2 + ddepths[-length(ddepths)], include.lowest=TRUE ) ))

      mm = match( as.character(M$StrataID), as.character(B$StrataID) )
      tokeep = c("z", "z.sd")
      M = cbind(M, B[mm, tokeep])

      H = carstm_hyperparameters( sd(log(M$substrate.grainsize), na.rm=TRUE), alpha=0.5, median( log(M$substrate.grainsize), na.rm=TRUE) )
      fit = inla(
        formula =
          substrate.grainsize ~ 1 +
            + f(strata, model="bym2", graph=sppoly@nb, scale.model=TRUE, constr=TRUE, hyper=H$bym2)
            + f(z, model="rw2", hyper=H$rw2)
            + f(z.sd, model="rw2", hyper=H$rw2)
            + f(iid_error, model="iid", hyper=H$iid)
          ,
        family = "lognormal",
        data= M,
        control.compute=list(dic=TRUE, config=TRUE),
        control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
        control.predictor=list(compute=FALSE, link=1 ),
        # control.fixed=H$fixed,  # priors for fixed effects, generic is ok
        # control.inla=list(int.strategy="eb") ,# to get empirical Bayes results much faster.
        # control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
        num.threads=2,
        blas.num.threads=2,
        verbose=TRUE
      )
      s = summary(fit)
      s$dic$dic  # 31225
      s$dic$p.eff # 5200

      plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )

      # reformat predictions into matrix form
      ii = which(M$tag=="predictions")
      jj = match(M$StrataID[ii], sppoly$StrataID)
      sppoly@data$substrate.grainsize.predicted = exp( fit$summary.fitted.values[ ii[jj], "mean" ]) - 2500
      sppoly@data$substrate.grainsize.predicted_lb = exp( fit$summary.fitted.values[ ii[jj], "0.025quant" ]) - 2500
      sppoly@data$substrate.grainsize.predicted_ub = exp( fit$summary.fitted.values[ ii[jj], "0.975quant" ]) - 2500
      sppoly@data$substrate.grainsize.random_strata_nonspatial = exp( fit$summary.random$strata[ jj, "mean" ])
      sppoly@data$substrate.grainsize.random_strata_spatial = exp( fit$summary.random$strata[ jj+max(jj), "mean" ])
      sppoly@data$substrate.grainsize.random_sample_iid = exp( fit$summary.random$iid_error[ ii[jj], "mean" ])

      attr( spplot, "fit") = fit
      save( spplot, file=fn, compress=TRUE )

      vn = "substrate.grainsize.predicted"
      brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
      dev.new()
      spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )

      return( sppoly )

    }

  }
