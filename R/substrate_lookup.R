substrate_lookup = function( LOCS=NULL, lookup_from="core", lookup_to="points", FUNC=mean,  vnames="substrate.grainsize", lookup_from_class="aggregated_data", pS=NULL  ) {
 
 # substrate_lookup( LOCS=M[, c("lon", "lat")], spatial_domain=p$spatial_domain, lookup_from="core", lookup_to="points" , lookup_from_class="aggregated_data" ) # core=="rawdata"
message("need to check::  [match( APS$AUID, as.character( sppoly$AUID ) )] ")

  if (is.null(pS)) pS = substrate_parameters(  project_class=lookup_from  )


  crs_lonlat =  st_crs(projection_proj4string("lonlat_wgs84"))


  if ( lookup_from %in% c("core") & lookup_to == "points" )  {
    # matching to point (LU) to point (LOCS) 
    vn = "substrate.grainsize"
    vn2 = paste( vn, "mean", sep="." )
 
    pB = bathymetry_parameters( spatial_domain=pS$spatial_domain, project_class=lookup_from  )
    BA = bathymetry_db ( p=pB, DS="baseline_prediction_locations", varnames=c("lon", "lat")  )
    BA = lonlat2planar(BA, pB$aegis_proj4string_planar_km)

    LU = substrate_db ( p=pS, DS=lookup_from_class )  # raw data
    names(LU)[ which(names(LU) == vn2 ) ] =  vn

    LOCS = lonlat2planar(LOCS, proj.type=pB$aegis_proj4string_planar_km) # get planar projections of lon/lat in km
    LOCS[, vn] = LU[ match(
        array_map( "xy->1", LOCS[, c("plon","plat")], gridparams=pB$gridparams ),
        array_map( "xy->1", BA[,c("plon","plat")], gridparams=pB$gridparams )
    ), vn ]
    return( LOCS[,vn] )
  }


  if ( lookup_from %in% c("core") & lookup_to == "areal_units" )  {
    # point (LU) -> areal unit (LOCS)
    vn = "substrate.grainsize"
    vn2 = paste( vn, "mean", sep="." )

    pB = bathymetry_parameters( spatial_domain=pS$spatial_domain, project_class=lookup_from  )
    BA = bathymetry_db ( p=pB, DS="baseline_prediction_locations", varnames=c("lon", "lat")  )
    BA = planar2lonlat(BA, pB$aegis_proj4string_planar_km)

    LU = substrate_db ( p=pB, DS=lookup_from_class )  # raw data
    LU = cbind( LU, BA )
    LU = sf::st_as_sf( LU, coords=c("lon", "lat") )
    st_crs(LU) = st_crs( projection_proj4string("lonlat_wgs84") )
    LU = sf::st_transform( LU, crs=st_crs(LOCS) )

    LOCS[, vn] = aggregate( LU[, vn2], LOCS, FUNC, na.rm=TRUE ) [[vn2]] 
    return( st_drop_geometry(LOCS)[,vn] )
  }


  if ( lookup_from %in% c("stmv", "hybrid") & lookup_to == "points" )  {
    # matching to point (LU) to point (LOCS)
    # if any still missing then use stmv depths
    pB = bathymetry_parameters( spatial_domain=pS$spatial_domain, project_class=lookup_from  )
    BA = bathymetry_db ( p=pB, DS="baseline", varnames=c("lon", "lat")  )
    BA = planar2lonlat(BA, proj.type=pS$aegis_proj4string_planar_km)

    LU = substrate_db ( p=pS, DS="complete", varnames="all" )  # raw data
    LU = cbind(LU, BA )
    
    LOCS = lonlat2planar(LOCS, proj.type=pB$aegis_proj4string_planar_km) # get planar projections of lon/lat in km
    LOCS[,vnames] = LU[ match(
        array_map( "xy->1", LOCS[, c("plon","plat")], gridparams=pB$gridparams ),
        array_map( "xy->1", LU[,c("plon","plat")], gridparams=pB$gridparams )
    ), vnames ]
    return( LOCS[,vnames] )
  }

  if ( lookup_from %in% c("stmv", "hybrid") & lookup_to == "areal_units" )  {
    # point (LU) -> areal unit (LOCS)
    pB = bathymetry_parameters( spatial_domain=pS$spatial_domain, project_class=lookup_from  )
    BA = bathymetry_db ( p=pB, DS="baseline", varnames=c("lon", "lat")  )
    BA = planar2lonlat(BA, pB$aegis_proj4string_planar_km)
    
    LU = substrate_db ( p=pS, DS="complete", varnames="all" )  # raw data
    LU = cbind(LU, BA )
    LU = sf::st_as_sf( LU, coords=c("lon", "lat") )
    st_crs(LU) = st_crs( projection_proj4string("lonlat_wgs84") )
    LU = sf::st_transform( LU, crs=st_crs(LOCS) )

    for (vn in vnames) {
      LOCS[, vn] = aggregate( LU[, vn], LOCS, FUNC, na.rm=TRUE ) [[vn]] 
    }
    return( st_drop_geometry(LOCS)[,vnames] )
  }



  if ( lookup_from %in% c("carstm" ) & lookup_to == "points" )  {
    #  areal unit (LU) to point (LOCS)
    LU_summ = carstm_model( p=pB, DS="carstm_modelled_summary" ) 
    if (is.null(LU_summ)) stop("Carstm predicted fields not found")

    LU = areal_units( p=pB )  #  poly associated with LU
    LU = sf::st_transform( LU, crs=st_crs(pB$aegis_proj4string_planar_km) )
    bm = match( LU$AUID, LU_summ$AUID )
    LU[[vnames]] = LU_summ[,vnames][ bm ]
    LU_summ = NULL

    LOCS = sf::st_as_sf( LOCS, coords=c("lon", "lat") )
    st_crs(LOCS) = st_crs( projection_proj4string("lonlat_wgs84") )
    LOCS = sf::st_transform( LOCS, crs=st_crs(pB$aegis_proj4string_planar_km) )

    raster_template = raster( LOCS, res=min(pB$gridparams$res), crs=st_crs( LOCS ) ) 
    for (vn in vnames) {
      LL = fasterize::fasterize( LU, raster_template, field=vn )
      o = sf::st_as_sf( as.data.frame( raster::rasterToPoints(LL)), coords=c("x", "y") )
      st_crs(o) = st_crs( LOCS )
    
      LOCS[,vn] = st_drop_geometry(o)[ match(
        array_map( "xy->1", st_coordinates(LOCS), gridparams=pB$gridparams ),
        array_map( "xy->1", st_coordinates(o), gridparams=pB$gridparams )
      ), "layer" ]
    }
    
    return( st_drop_geometry(LOCS)[,vnames] )
  }


  if ( lookup_from %in% c("carstm") & lookup_to == "areal_units" )  {
    # areal unit (LU) to areal unit (LOCS)
    LU_summ = carstm_model( p=pB, DS="carstm_modelled_summary" ) 
    if (is.null(LU_summ)) stop("Carstm predicted fields not found")

    LU = areal_units( p=pB )  #  poly associated with LU
    bm = match( LU$AUID, LU_summ$AUID )
    LU[[vnames]] = LU_summ[,vnames][ bm ]
    LU_summ = NULL

    # now rasterize and re-estimate
    raster_template = raster( LOCS, res=min(pB$gridparams$res), crs=st_crs( LOCS ) ) 
    for (vn in vnames) {
      LL = fasterize::fasterize( LU, raster_template, field=vn )
      o = sf::st_as_sf( as.data.frame( raster::rasterToPoints(LL)), coords=c("x", "y") )
      st_crs(o) = st_crs( LOCS )
      LOCS[, vn] = sf:::aggregate.sf( o, LOCS, FUNC, na.rm=TRUE ) [["layer"]]  
    }
    return( st_drop_geometry(LOCS)[,vnames] )
  } 

}
