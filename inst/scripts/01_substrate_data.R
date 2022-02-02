
## NOTE:: substrate size is really only relevant for SSE/snowcrab domain right now as no
## NOTE:: other data source has been found/identified
## NOTE:: but working at the size of canada.east.highres for compatibility with bathymetry


## TODO:: add data collected by snow crab survey and any others for that matter
## TODO::
## TODO::


require(aegis.substrate)

substrate_db ( DS="substrate.initial.redo" ) # bring in Kostelev's data ... stored as a SpatialGridDataFrame
substrate_db ( DS="lonlat.highres.redo" ) # in future .. additional data would be added here

p = aegis.substrate::substrate_parameters()

M = substrate_db( p=p, DS="aggregated_data", redo=TRUE )  # will redo if not found .. not used here but used for data matching/lookup in other aegis projects that use bathymetry

