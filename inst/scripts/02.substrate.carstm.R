
# construct basic parameter list defining the main characteristics of the study
# and some plotting parameters (bounding box, projection, bathymetry layout, coastline)
  p = aegis.substrate::substrate_carstm(
    DS = "parameters",
    project_class = "carstm", # defines which parameter set to load
    inputdata_spatial_discretization_planar_km = 1,  # km controls resolution of data prior to modelling to reduce data set and speed up modelling
    spatial_domain = "SSE",  # defines spatial area, currenty: "snowcrab" or "SSE"
    variabletomodel ="substrate.grainsize",
    areal_units_resolution_km = 25, # km dim of lattice ~ 1 hr
    areal_units_proj4string_planar_km = projection_proj4string("utm20")  # coord system to use for areal estimation and gridding for carstm
  )

# prepare data
  # sppoly = areal_units( p=p, redo=TRUE )  # this has already been done in aegis.polygons::01 polygons.R .. should nto have to redo
  sppoly = areal_units( p=p )  # this has already been done in aegis.polygons::01 polygons.R .. should nto have to redo
  plot(sppoly)

  M = substrate_carstm( p=p, DS="aggregated_data", redo=TRUE )  # will redo if not found .. not used here but used for data matching/lookup in other aegis projects that use bathymetry
  M = substrate_carstm( p=p, DS="carstm_inputs", redo=TRUE )  # will redo if not found

# run model and obtain predictions
  res = carstm_model( p=p, M=M )

# extract results and examine
  res = carstm_model( p=p, DS="carstm_modelled" ) # to load currently saved results
  fit =  carstm_model( p=p, DS="carstm_modelled_fit" )  # extract currently saved model fit
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

  vn = paste(p$variabletomodel, "random_strata_nonspatial", sep=".")
  carstm_plot( p=p, res=res, vn=vn )

  vn = paste(p$variabletomodel, "random_strata_spatial", sep="._carstm  carstm_plot( p=p, res=res, vn=vn )


# end
