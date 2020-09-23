
## NOTE:: substrate size is really only relevant for SSE/snowcrab domain right now as no
##        other data source has been found/identified
##        but working at the size of canada.east.highres for compatibility with bathymetry
## TODO:: add data collected by snow crab survey and any others for that matter
require(aegis.substrate)

substrate_db ( DS="substrate.initial.redo" ) # bring in Kostelev's data ... stored as a SpatialGridDataFrame
substrate_db ( DS="lonlat.highres.redo" ) # in future .. additional data would be added here

p = substrate_parameters(
  data_root = project.datadirectory( "aegis", "substrate" ),
  spatial_domain = "canada.east.highres" ,
  inputdata_spatial_discretization_planar_km = 0.5, # p$pres==0.5; controls resolution of data prior to modelling (km .. ie 20 linear units smaller than the final discretization pres)
  aegis_dimensionality="space"
)
M = substrate_db( p=p, DS="aggregated_data", redo=TRUE )  # will redo if not found .. not used here but used for data matching/lookup in other aegis projects that use bathymetry
