
# construct basic parameter list defining the main characteristics of the study
# and some plotting parameters (bounding box, projection, bathymetry layout, coastline)

  p = aegis.substrate::substrate_parameters( project_class="carstm", areal_units_resolution_km=5  )  #10 km grid gives about 2000 au .. in the right range for optimal solutions that are not too slow, 5 km gives stable results .. 


    # adjust based upon RAM requirements and ncores
    # inla.setOption( num.threads="6:2" )
    # inla.setOption(blas.num.threads= 3 )

    if (0) {
      p$fraction_todrop = 1/4 # aggressiveness of solution finding ( fraction of counts to drop each iteration)
      p$fraction_cv = 1.0  #sd/mean no.
      p$fraction_good_bad = 0.9
      p$areal_units_constraint_ntarget = 500  # length(p$yrs)
      p$areal_units_constraint_nmin = 30  # length(p$yrs)
      p$nAU_min = 100
      # to recreate the underlying data
      # xydata=substrate_db(p=p, DS="areal_units_input", redo=TRUE)

      sppoly = areal_units( p=p , redo=T )  # this has already been done in aegis.polygons::01 polygons.R .. should nto have to redo
      plot( sppoly[ "AUID" ] )


      # prepare data

      sppoly = areal_units( p=p  )  # this has already been done in aegis.polygons::01 polygons.R .. should nto have to redo
   
      M = substrate_db( p=p, DS="carstm_inputs", sppoly=sppoly, redo=TRUE )  # will redo if not found
  
    }

  # run model and obtain predictions
    fit = carstm_model( 
      p=p, 
      sppoly=areal_units( p=p  ),
      data= substrate_db( p=p, DS="carstm_inputs", sppoly=sppoly), 
      num.threads="4:2",
      theta = c( 1.710, 3.729, 0.008, 5.591 ), 
      redo_fit = TRUE,  
      verbose=TRUE 
    ) 
    # fit = carstm_model( p=p, DS="carstm_modelled_fit" )  # extract currently saved model fit
    
      # extract results
      if (0) {
        fit = carstm_model( p=p, data=M, sppoly=sppoly ) # alt way of running
        # very large files .. slow 
        fit = carstm_model( p=p, DS="carstm_modelled_fit", sppoly=sppoly )  # extract currently saved model fit
    
        fit$summary$dic$dic
        fit$summary$dic$p.eff
        fit$dyear

        plot(fit)
        plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
      }

  # extract results and examine
    res = carstm_model( p=p, DS="carstm_modelled_summary", sppoly=sppoly  ) # to load currently saved results
    res$summary

    
    # maps of some of the results
    outputdir = file.path(p$data_root, "maps", p$carstm_model_label )
    if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

    tmout = carstm_map(  res=res, vn="predictions", 
        sppoly=sppoly,
        palette="viridis",
        title="Substrate grainsize (mm)", 
        plot_elements=c( "isobaths",  "compass", "scale_bar", "legend" ),
        tmap_zoom= c((p$lon0+p$lon1)/2 - 0.5, (p$lat0+p$lat1)/2 -0.8, 6.5)
    )  
    tmout

    outfilename= file.path( outputdir, paste("substrate_grain_size_carstm", "png", sep=".") )
    mapview::mapshot( tmap_leaflet(tmout), file=outfilename, vwidth = 1600, vheight = 1200 )
    
  # random effects  ..i.e.,  deviation from lognormal model
    tmout = carstm_map(  res=res, vn= c( "random", "space", "combined" ), 
        sppoly=sppoly,
        palette="viridis",
        title="Substrate grainsize spatial errors (mm)",
        plot_elements=c( "isobaths",  "compass", "scale_bar", "legend" ),
        tmap_zoom= c((p$lon0+p$lon1)/2-0.5, (p$lat0+p$lat1)/2 -0.8, 6.5)
    )  
    tmout
  
    outfilename= file.path( outputdir, paste("substrate_grain_size_spatialeffect_carstm", "png", sep=".") )
    mapview::mapshot( tmap_leaflet(tmout), file=outfilename, vwidth = 1600, vheight = 1200 )
    

  # end
