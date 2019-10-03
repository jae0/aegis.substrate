
# construct basic parameter list defining the main characteristics of the study
# and some plotting parameters (bounding box, projection, bathymetry layout, coastline)
p = aegis.substrate::substrate_parameters(
  project_class = "carstm", # defines which parameter set to load
  inputdata_spatial_discretization_planar_km = 1,  # km controls resolution of data prior to modelling to reduce data set and speed up modelling
  spatial_domain = "snowcrab",  # defines spatial area, currenty: "snowcrab" or "SSE"
  # spatial_domain = "SSE",  # defines spatial area, currenty: "snowcrab" or "SSE"
  areal_units_resolution_km = 25, # km dim of lattice ~ 1 hr
  areal_units_proj4string_planar_km = projection_proj4string("utm20")  # coord system to use for areal estimation and gridding for carstm
)


if (0) {
  # run model and obtain predictions
  sppoly = areal_units( p=p, redo=TRUE )  # this has already been done in aegis.polygons::01 polygons.R .. should nto have to redo
  M = substrate.db( p=p, DS="aggregated_data", redo=TRUE )  # will redo if not found .. not used here but used for data matching/lookup in other aegis projects that use bathymetry
  M = substrate_carstm( p=p, DS="carstm_inputs", redo=TRUE )  # will redo if not found
  sppoly = substrate_carstm( p=p, DS="carstm_modelled", redo=TRUE )


  sppoly = substrate_carstm( p=p, DS="carstm_modelled" ) # to load currently saved sppoly
  fit =  substrate_carstm( p=p, DS="carstm_modelled_fit" )  # extract currently saved model fit
  plot(fit)
  s = summary(fit)
  s$dic$dic
  s$dic$p.eff


  # maps of some of the results
  p$boundingbox = list( xlim=p$corners$lon, ylim=p$corners$lat) # bounding box for plots using spplot
  p$mypalette = RColorBrewer::brewer.pal(9, "YlOrRd")
  p = c(p, aegis.coastline::coastline_layout( p=p  ) )  # set up default map projection

  vn = "substrate.grainsize.predicted"
  dev.new();
  spplot( sppoly, vn, main=vn,
    col.regions=p$mypalette,
    at=interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile"),
    sp.layout=p$coastLayout,
    col="transparent"
  )

  vn = "substrate.grainsize.random_strata_nonspatial"
  dev.new();
  spplot( sppoly, vn, main=vn,
    col.regions=p$mypalette,
    at=interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile"),
    sp.layout=p$coastLayout,
    col="transparent"
  )

  vn = "substrate.grainsize.random_strata_spatial"
  dev.new();
  spplot( sppoly, vn, main=vn,
    col.regions=p$mypalette,
    at=interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile"),
    sp.layout=p$coastLayout,
    col="transparent"
  )

}



# end
