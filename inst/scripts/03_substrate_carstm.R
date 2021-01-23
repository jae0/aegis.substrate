
# construct basic parameter list defining the main characteristics of the study
# and some plotting parameters (bounding box, projection, bathymetry layout, coastline)

  p = aegis.substrate::substrate_parameters( project_class="carstm"  )


    # adjust based upon RAM requirements and ncores
    inla.setOption(num.threads= floor( parallel::detectCores() / 3 ) )
    inla.setOption(blas.num.threads= 3 )

    if (0) {
      p$fraction_todrop = 1/4 # aggressiveness of solution finding ( fraction of counts to drop each iteration)
      p$fraction_cv = 1.0  #sd/mean no.
      p$fraction_good_bad = 0.9
      p$areal_units_constraint_nmin = 500  # length(p$yrs)
      p$nAU_min = 100
      # to recreate the underlying data
      # xydata=substrate_db(p=p, DS="areal_units_input", redo=TRUE)

      sppoly = areal_units( p=p , redo=T )  # this has already been done in aegis.polygons::01 polygons.R .. should nto have to redo
      plot( sppoly[ "AUID" ] )


      # prepare data
      M = substrate_db( p=p, DS="carstm_inputs", redo=TRUE )  # will redo if not found
  
    }

# run model and obtain predictions
  fit = carstm_model( p=p, M='substrate_db( p=p, DS="carstm_inputs" )' ) 
  
    # extract results
    if (0) {
      fit = carstm_model( p=p, M=M ) # alt way of running
      # very large files .. slow 
      fit = carstm_model( p=p, DS="carstm_modelled_fit" )  # extract currently saved model fit
      plot(fit)
      plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
    }

# extract results and examine
  res = carstm_model( p=p, DS="carstm_modelled_summary"  ) # to load currently saved results
  
    res$summary$dic$dic
    res$summary$dic$p.eff
    res$dyear




  # maps of some of the results
  vn = paste(p$variabletomodel, "predicted", sep=".")
  outputdir = file.path( gsub( ".rdata", "", dirname(res$fn_res) ), "figures", vn )
  if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

  fn = file.path( outdir, paste("substrate_grain_size", "png", sep=".") )
  png( filename=fn, width=3072, height=2304, pointsize=40, res=300  )
    vn = paste(p$variabletomodel, "predicted", sep=".")
    print( carstm_map( res=res, vn=vn, at=p$discretization$z, sp.layout=p$coastLayout, col.regions=p$mypalette, main="Substrate grain size"  ) )
  dev.off()


  vn = paste(p$variabletomodel, "random_sample_iid", sep=".")
  carstm_map( res=res, vn=vn, at=p$discretization$z, sp.layout=p$coastLayout, col.regions=p$mypalette, main="Substrate sample iid"  )

  vn = paste(p$variabletomodel, "random_auid_nonspatial", sep=".")
  carstm_map( res=res, vn=vn, at=p$discretization$z, sp.layout=p$coastLayout, col.regions=p$mypalette, main="Substrate nonspatial"  )

  vn = paste(p$variabletomodel, "random_auid_spatial", sep="." )
  carstm_map( res=res, vn=vn, at=p$discretization$z, sp.layout=p$coastLayout, col.regions=p$mypalette, main="Substrate spatial effects"  )



# end
