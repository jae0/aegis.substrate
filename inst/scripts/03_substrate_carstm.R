
# construct basic parameter list defining the main characteristics of the study
# and some plotting parameters (bounding box, projection, bathymetry layout, coastline)

  p = aegis.substrate::substrate_parameters( project_class="carstm"  )


    # adjust based upon RAM requirements and ncores
    inla.setOption(num.threads= floor( parallel::detectCores() / 3 ) )
    inla.setOption(blas.num.threads= 3 )

if(0) {
      p$fraction_todrop = 1/4 # aggressiveness of solution finding ( fraction of counts to drop each iteration)
      p$fraction_cv = 1.0  #sd/mean no.
      p$fraction_good_bad = 0.9
      p$areal_units_constraint_nmin = 500  # length(p$yrs)
      p$nAU_min = 100
}
  # to recreate the underlying data
  # xydata=substrate_db(p=p, DS="areal_units_input", redo=TRUE)

  sppoly = areal_units( p=p , redo=T )  # this has already been done in aegis.polygons::01 polygons.R .. should nto have to redo
  plot( sppoly[ "AUID" ] )


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
