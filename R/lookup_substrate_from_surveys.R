
lookup_substrate_from_surveys = function( p, locs, vnames="substrate.grainsize.mean" ) {

  pS = spatial_parameters( spatial_domain=p$spatial_domain )
  if (!exists("inputdata_spatial_discretization_planar_km", pS)) {
    if (!exists("inputdata_spatial_discretization_planar_km", p)) {
      pS$inputdata_spatial_discretization_planar_km = p$inputdata_spatial_discretization_planar_km
    } else {
      pS$inputdata_spatial_discretization_planar_km = 1
    }
  }
  if (!exists("variabletomodel", pS)) pS$variabletomodel = "substrate.grainsize"

  B = substrate.db ( p=pS, DS="aggregated_data" )
  B = lonlat2planar( B, proj.type=pS$aegis_proj4string_planar_km )
  B$plon = round(B$plon / pS$inputdata_spatial_discretization_planar_km + 1 ) * pS$inputdata_spatial_discretization_planar_km
  B$plat = round(B$plat / pS$inputdata_spatial_discretization_planar_km + 1 ) * pS$inputdata_spatial_discretization_planar_km

  locs = lonlat2planar( locs, proj.type=pS$aegis_proj4string_planar_km )
  locs$plon = round(locs$plon / pS$inputdata_spatial_discretization_planar_km + 1 ) * pS$inputdata_spatial_discretization_planar_km
  locs$plat = round(locs$plat / pS$inputdata_spatial_discretization_planar_km + 1 ) * pS$inputdata_spatial_discretization_planar_km

  locs_map = paste(locs$plon, locs$plat, sep=".")
  domain_map = paste(B$plon, B$plat, sep=".")
  locs_index = match( locs_map, domain_map )
  out = B[locs_index, vnames]
  return(out)

}


