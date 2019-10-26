
lookup_substrate_from_surveys = function( p, locs, vnames="substrate.grainsize.mean" ) {

  B = substrate.db ( p=p, DS="aggregated_data" )  # 16 GB in RAM just to store!
  varnames = unique( c( "plon", "plat", vnames ) )
  domain_map = stmv::array_map( "xy->1", B[,c("plon","plat")], gridparams=p$gridparams )
  locs_map = stmv::array_map( "xy->1", locs, gridparams=p$gridparams )
  locs_index = match( locs_map, domain_map )
  out = B[locs_index, vnames]
  return(out)

}


