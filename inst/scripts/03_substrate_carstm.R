
# construct basic parameter list defining the main characteristics of the study
# and some plotting parameters (bounding box, projection, bathymetry layout, coastline)

  p = aegis.substrate::substrate_parameters( project_class="carstm"  )
  p$inla_num.threads=  parallel::detectCores() 
  p$inla_blas.num.threads=1 


  sppoly = areal_units( p=p )  # this has already been done in aegis.polygons::01 polygons.R .. should nto have to redo
  plot(sppoly) # spplot( sppoly, "AUID", main="AUID", sp.layout=p$coastLayout )

  if (redo.inputs) {
    # prepare data
    sppoly = areal_units( p=p, redo=TRUE )  # this has already been done in aegis.polygons::01 polygons.R .. should nto have to redo
    M = substrate_db( p=p, DS="carstm_inputs", redo=TRUE )  # will redo if not found
  }

# run model and obtain predictions
  fit = carstm_model( p=p, M=M )

# extract results and examine
  fit =  carstm_model( p=p, DS="carstm_modelled_fit" )  # extract currently saved model fit
  res = carstm_summary( p=p  ) # to load currently saved results

  plot(fit)
  plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
  s = summary(fit)
  s$dic$dic
  s$dic$p.eff


  # maps of some of the results
  vn = paste(p$variabletomodel, "predicted", sep=".")
  carstm_plot( p=p, res=res, vn=vn )

  vn = paste(p$variabletomodel, "random_sample_iid", sep=".")
  carstm_plot( p=p, res=res, vn=vn )

  vn = paste(p$variabletomodel, "random_auid_nonspatial", sep=".")
  carstm_plot( p=p, res=res, vn=vn )

  vn = paste(p$variabletomodel, "random_auid_spatial", sep="._carstm  carstm_plot( p=p, res=res, vn=vn )


# end