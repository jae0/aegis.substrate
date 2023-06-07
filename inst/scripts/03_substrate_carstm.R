
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

      sppoly = areal_units( p=p , hull_ratio=0.01, redo=TRUE )  # this has already been done in aegis.polygons::01 polygons.R .. should nto have to redo
      plot( sppoly[ "AUID" ] )


      # prepare data

      sppoly = areal_units( p=p  )  # this has already been done in aegis.polygons::01 polygons.R .. should nto have to redo
   
      M = substrate_db( p=p, DS="carstm_inputs", sppoly=sppoly, redo=TRUE )  # will redo if not found
  
    }

  # run model and obtain predictions
    res = carstm_model( 
      p=p, 
      sppoly=areal_units( p=p  ),
      data= substrate_db( p=p, DS="carstm_inputs", sppoly=sppoly), 
      space_id = sppoly$AUID,
      nposteriors=1000,
      redo_fit=TRUE, # to start optim from a solution close to the final in 2021 ... 
      # redo_fit=FALSE, # to start optim from a solution close to the final in 2021 ... 
      # debug = TRUE,
      theta = c( 1.710, 3.588, 0.008, 5.662 ) ,
      # control.mode = list( restart=FALSE, theta= c( 1.710, 3.588, 0.008, 5.662 ) ),  
      num.threads="4:2",
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

    
    # bbox = c(-71.5, 41, -52.5,  50.5 )
    additional_features = additional_features_tmap( 
        p=p, 
        isobaths=c( 10, 100, 200, 300, 500, 1000 ), 
        coastline =  c("canada"), 
        xlim=c(-80,-40), 
        ylim=c(38, 60) 
    )

    # maps of some of the results
    outputdir = file.path(p$data_root, "maps", p$carstm_model_label )
    if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

    outfilename= file.path( outputdir, paste("substrate_grain_size_carstm", "png", sep=".") )

    tmout = carstm_map(  res=res, vn="predictions", 
        sppoly=sppoly,
        palette="viridis",
        title="Substrate grainsize (mm)", 
        plot_elements=c(  "compass", "scale_bar", "legend" ),
        additional_features=additional_features,
        outfilename=outfilename
    )  
    tmout

 

  # random effects  ..i.e.,  deviation from lognormal model ( pure spatial effect )
    outfilename= file.path( outputdir, paste("substrate_grain_size_spatialeffect_carstm", "png", sep=".") )
    tmout = carstm_map(  res=res, vn= c( "random", "space", "combined" ), 
        sppoly=sppoly,
        palette="viridis",
        title="Substrate grainsize spatial errors (mm)",
        plot_elements=c(  "compass", "scale_bar", "legend" ),
        additional_features=additional_features,
        outfilename=outfilename
    )  
    tmout
  
 



  # end
