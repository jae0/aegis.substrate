substrate_lookup = function( LOCS=NULL, spatial_domain=NULL, lookup_from="core", lookup_to="points", FUNC=mean,  vnames="substrate.grainsize", lookup_from_class="aggregated_data"  ) {
  # lookup from rawdata
 # substrate_lookup( LOCS=M[, c("lon", "lat")], spatial_domain=p$spatial_domain, lookup_from="core", lookup_to="points" , lookup_from_class="aggregated_data" ) # core=="rawdata"

  if (is.null(spatial_domain))  {
    pS = substrate_parameters(  project_class=lookup_from  )
  } else {
    pS = substrate_parameters( spatial_domain=spatial_domain, project_class=lookup_from  )
  }

  crs_lonlat =  st_crs(projection_proj4string("lonlat_wgs84"))


  if ( lookup_from %in% c("core") & lookup_to == "points" )  {
    # matching to point to point 
    # if any still missing then use stmv depths
    LU = substrate_db ( p=pB, DS=lookup_from_class )  # raw data
    LU = planar2lonlat(LU, pB$aegis_proj4string_planar_km)

    vn2 = "substrate.grainsize.mean"
    
    LOCS = lonlat2planar(LOCS, proj.type=pB$aegis_proj4string_planar_km) # get planar projections of lon/lat in km
    LOCS[, vnames] = LU[ match(
        array_map( "xy->1", LOCS[, c("plon","plat")], gridparams=pB$gridparams ),
        array_map( "xy->1", LU[,c("plon","plat")], gridparams=pB$gridparams )
    ), vn2 ]
    return( LOCS[,vnames] )
  }

  if ( lookup_from %in% c("core") & lookup_to == "areal_units" )  {
    # point -> areal unit
    LU = substrate_db ( p=pB, DS=lookup_from_class )  # raw data
    LU = planar2lonlat(LU, pB$aegis_proj4string_planar_km)
    
    LU = sf::st_as_sf( LU, coords=c("lon", "lat") )
    st_crs(LU) = st_crs( projection_proj4string("lonlat_wgs84") )
    LU = sf::st_transform( LU, crs=st_crs(LOCS) )
    vn2 = "substrate.grainsize.mean"
    LOCS[, vnames] = aggregate( LU[, vn2], LOCS, FUNC, na.rm=TRUE ) [[vn2]] [iAS]
    return( LOCS[,vnames] )
  }


  if ( lookup_from %in% c("stmv", "hybrid") & lookup_to == "points" )  {
    # matching to point to point 
    # if any still missing then use stmv depths
    pB = bathymetry_parameters( spatial_domain=pS$spatial_domain, project_class=lookup_from  )
    BA = bathymetry_db ( pB, DS="baseline", varnames=c("lon", "lat")  )
    LU = substrate_db ( pS, DS="complete", varnames="all" )  # raw data
    LU = cbind(LU, BA )
    LU = planar2lonlat(LU, proj.type=pS$aegis_proj4string_planar_km)
    
    LOCS = lonlat2planar(LOCS, proj.type=pB$aegis_proj4string_planar_km) # get planar projections of lon/lat in km
    LOCS[,vnames] = LU[ match(
        array_map( "xy->1", LOCS[, c("plon","plat")], gridparams=pS$gridparams ),
        array_map( "xy->1", LU[,c("plon","plat")], gridparams=pS$gridparams )
    ), vnames ]
    return( LOCS[,vnames] )
  }

  if ( lookup_from %in% c("stmv", "hybrid") & lookup_to == "areal_units" )  {
    # point -> areal unit
    pB = bathymetry_parameters( spatial_domain=pS$spatial_domain, project_class=lookup_from  )
    BA = bathymetry_db ( pB, DS="baseline", varnames=c("lon", "lat")  )
    LU = substrate_db ( pS, DS="complete", varnames="all" )  # raw data
    LU = cbind(LU, BA )
    LU = planar2lonlat(LU, pB$aegis_proj4string_planar_km)
    
    LU = sf::st_as_sf( LU, coords=c("lon", "lat") )
    st_crs(LU) = st_crs( projection_proj4string("lonlat_wgs84") )
    LU = sf::st_transform( LU, crs=st_crs(LOCS) )
    for (vn in vnames) {
      LOCS[, vn] = aggregate( LU[, vn], LOCS, FUNC, na.rm=TRUE ) [[vn]] [iAS]
    }
    return( LOCS[,vnames] )
  }



  if ( lookup_from %in% c("carstm" ) & lookup_to == "points" )  {
    # point to areal unit
    LU = carstm_model( p=pB, DS="carstm_modelled_summary" ) 
    if (is.null(LU)) stop("Carstm predicted fields not found")

    AU = areal_units( p=pB )  #  poly associated with LU
    bm = match( AU$AUID, LU$AUID )
    AU[[vnames]] = LU[,vnames][ bm ]
    LU = NULL
    # now rasterize and re-estimate

    LOCS = lonlat2planar(LOCS, proj.type=pB$aegis_proj4string_planar_km) # get planar projections of lon/lat in km
    raster_template = raster( LOCS, res=pB$areal_units_resolution_km, crs=st_crs( LOCS ) ) # +1 to increase the area
    for (vn in vnames) {
      LL = fasterize::fasterize( AU, raster_template, field=vn )
      LOCS[[vn]] = sp::over( LOCS, LL[, vn ], fn=FUNC, na.rm=TRUE )
    }
    return( LOCS[,vnames] )
  }


  if ( lookup_from %in% c("carstm") & lookup_to == "areal_units" )  {
    # areal unit to areal unit
    LU = carstm_model( p=pB, DS="carstm_modelled_summary" ) 
    if (is.null(LU)) stop("Carstm predicted fields not found")

    AU = areal_units( p=pB )  #  poly associated with LU
    bm = match( AU$AUID, LU$AUID )
    AU[[vnames]] = LU[,vnames][ bm ]
    LU = NULL

    # now rasterize and re-estimate
    raster_template = raster( LOCS, res=pB$areal_units_resolution_km, crs=st_crs( LOCS ) ) # +1 to increase the area
    for (vn in vnames) {
      LL = fasterize::fasterize( AU, raster_template, field=vn )
      LOCS[[vn]] = sp::over( LOCS, LL[, vn ], fn=FUNC, na.rm=TRUE )
    }
    return( LOCS[,vnames] )
  } 

}
