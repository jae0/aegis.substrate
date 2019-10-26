
lookup_substrate_from_stmv = function( p, locs, vnames="substrate.grainsize" ) {

  # only for stmv-based lookups
  varnames = unique( c( "plon", "plat", vnames ) )
  domain = substrate.db(p=p, DS="baseline", varnames=varnames )
  domain_map = stmv::array_map( "xy->1", domain[,c("plon","plat")], gridparams=p$gridparams )
  locs_map = stmv::array_map( "xy->1", locs, gridparams=p$gridparams )
  locs_index = match( locs_map, domain_map )
  out = domain[locs_index, vnames]
  return(out)

}


