
## NOTE:: substrate size is really only relevant for SSE/snowcrab domain right now as no
##        other data source has been found/identified
##        but working at the size of canada.east.highres for compatibility with bathymetry
## TODO:: add data collected by snow crab survey and any others for that matter

substrate.db ( DS="substrate.initial.redo" ) # bring in Kostelev's data ... stored as a SpatialGridDataFrame
substrate.db ( DS="lonlat.highres.redo" ) # in future .. additional data would be added here
M = substrate_carstm( p=p, DS="aggregated_data", redo=TRUE )  # will redo if not found .. not used here but used for data matching/lookup in other aegis projects that use bathymetry
