
# construct basic parameter list defining the main characteristics of the study
# and some plotting parameters (bounding box, projection, bathymetry layout, coastline)

  p = aegis.substrate::substrate_parameters( project_class="carstm", areal_units_resolution_km=10  )  #10 km grid gives about 2000 au .. in the right range for optimal solutions that are not too slow, 5 km gives unstable results .. 


    # adjust based upon RAM requirements and ncores
    inla.setOption( num.threads=6:2 )
    inla.setOption(blas.num.threads= 3 )

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
      M = substrate_db( p=p, DS="carstm_inputs", redo=TRUE )  # will redo if not found
  
    }

# run model and obtain predictions
  fit = carstm_model( 
    p=p, 
    data='substrate_db( p=p, DS="carstm_inputs" )', 
    num.threads="4:2",
    redo_fit = TRUE,  
    verbose=TRUE 
  ) 
  # fit = carstm_model( p=p, DS="carstm_modelled_fit" )  # extract currently saved model fit
  
    # extract results
    if (0) {
      fit = carstm_model( p=p, data=M ) # alt way of running
      # very large files .. slow 
      fit = carstm_model( p=p, DS="carstm_modelled_fit" )  # extract currently saved model fit
  
      fit$summary$dic$dic
      fit$summary$dic$p.eff
      fit$dyear

      plot(fit)
      plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
    }

# extract results and examine
  res = carstm_model( p=p, DS="carstm_modelled_summary"  ) # to load currently saved results
  res$summary

  
  # maps of some of the results
  outputdir = file.path( gsub( ".rdata", "", carstm_filenames(p, returntype="carstm_modelled_fit") ), "figures" )
  if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )


  tmout = carstm_map(  res=res, vn="predictions", 
      palette="viridis",
      title="Substrate grainsize (mm)", 
      plot_elements=c( "isobaths", "coastline", "compass", "scale_bar", "legend" ),
      outfilename= file.path( outputdir, paste("substrate_grain_size_carstm", "png", sep=".") ),
      tmap_zoom= c((p$lon0+p$lon1)/2 - 0.5, (p$lat0+p$lat1)/2 -0.8, 6.5)
  )  
  tmout

# random effects  ..i.e.,  deviation from lognormal model
  tmout = carstm_map(  res=res, vn= c( "random", "space", "combined" ), 
      palette="viridis",
      title="Substrate grainsize spatial errors (mm)",
      plot_elements=c( "isobaths", "coastline", "compass", "scale_bar", "legend" ),
      outfilename= file.path( outputdir, paste("substrate_grain_size_spatialeffect_carstm", "png", sep=".") ),
      tmap_zoom= c((p$lon0+p$lon1)/2-0.5, (p$lat0+p$lat1)/2 -0.8, 6.5)
  )  

  tmout

# end
