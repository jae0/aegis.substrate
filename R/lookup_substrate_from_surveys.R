
lookup_substrate_from_surveys = function( p, locs, vnames="substrate.grainsize.mean" ) {

  B = substrate.db ( p=p, DS="aggregated_data" )
  locs = lonlat2planar( locs, proj.type=p$aegis_proj4string_planar_km )
  locs$plon = round(locs$plon / p$inputdata_spatial_discretization_planar_km + 1 ) * p$inputdata_spatial_discretization_planar_km
  locs$plat = round(locs$plat / p$inputdata_spatial_discretization_planar_km + 1 ) * p$inputdata_spatial_discretization_planar_km
  locs_map = paste(locs$plon, locs$plat, sep=".")
  domain_map = paste(B$plon, B$plat, sep=".")
  locs_index = match( locs_map, domain_map )
  out = B[locs_index, vnames]
  return(out)

}


