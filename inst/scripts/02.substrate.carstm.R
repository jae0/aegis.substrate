
# construct basic parameter list defining the main characteristics of the study
# and some plotting parameters (bounding box, projection, bathymetry layout, coastline)
p = aegis.substrate::substrate_parameters(
  project_class = "carstm", # defines which parameter set to load
  inputdata_spatial_discretization_planar_km = 1,  # km controls resolution of data prior to modelling to reduce data set and speed up modelling
  # spatial_domain = "snowcrab",  # defines spatial area, currenty: "snowcrab" or "SSE"
  spatial_domain = "SSE",  # defines spatial area, currenty: "snowcrab" or "SSE"
  areal_units_strata_type = "lattice", # "aegis_lattice" to use ageis fields instead of carstm fields ... note variables are not the same
  #areal_units_overlay = "snowcrab_managementareas", # currently: "snowcrab_managementareas",  "groundfish_strata" .. additional polygon layers for subsequent analysis for now ..
  areal_units_overlay = "groundfish_strata", #.. additional polygon layers for subsequent analysis for now ..
  # areal_units_resolution_km = 10, # km dim of lattice ~ 16 hrs
  # areal_units_resolution_km = 20, # km dim of lattice ~ 1 hr
  areal_units_resolution_km = 25, # km dim of lattice ~ 1 hr
  areal_units_proj4string_planar_km = projection_proj4string("utm20"),  # coord system to use for areal estimation and gridding for carstm
  # areal_units_proj4string_planar_km = "+proj=omerc +lat_0=44.0 +lonc=-63.0 +gamma=0.0 +k=1 +alpha=325 +x_0=0 +y_0=0 +ellps=WGS84 +units=km",  # oblique mercator, centred on Scotian Shelf rotated by 325 degrees
  carstm_modelengine = "inla.default",  # {model engine}.{label to use to store}
  carstm_modelcall = '
    inla(
      formula = substrate.grainsize ~ 1
        + f(zi, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
        + f(strata, model="bym2", graph=sppoly@nb, scale.model=TRUE, constr=TRUE, hyper=H$bym2)
        + f(iid_error, model="iid", hyper=H$iid),
      family = "lognormal", # "zeroinflatedpoisson0",
      data= M,
      control.compute=list(dic=TRUE, config=TRUE),
      control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
      control.predictor=list(compute=FALSE, link=1 ),
      control.fixed=H$fixed,  # priors for fixed effects, generic is ok
      control.inla=list(int.strategy="eb") ,# to get empirical Bayes results much faster.
      # control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
      num.threads=4,
      blas.num.threads=4,
      verbose=TRUE
    ) ',
  # carstm_modelcall = 'glm( formula = substrate.grainsize ~ 1 + StrataID + log(z),  family = gaussian(link="log"), data= M[ which(M$tag=="observations"), ] ) ',  # for modelengine='glm'
  # carstm_modelcall = 'gam( formula = substrate.grainsize ~ 1 + StrataID + s(log(z)),  family = gaussian(link="log"), data= M[ which(M$tag=="observations"), ] ) ',  # for modelengine='gam'
  libs = RLibrary ( "sp", "spdep", "rgeos", "spatialreg", "INLA", "raster", "aegis",  "aegis.polygons", "aegis.bathymetry", "aegis.substrate", "carstm" )
)


# run model and obtain predictions
sppoly = substrate_carstm( p=p, DS="carstm_modelled", redo=TRUE )

if (0) {
    sppoly = areal_units( p=p, redo=TRUE )  # this has already been done in aegis.polygons::01 polygons.R .. should nto have to redo
    M = substrate_carstm( p=p, DS="aggregated_data", redo=TRUE )  # will redo if not found .. not used here but used for data matching/lookup in other aegis projects that use bathymetry
    M = substrate_carstm( p=p, DS="carstm_inputs", redo=TRUE )  # will redo if not found


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
