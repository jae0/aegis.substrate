substrate_lookup_rawdata = function( M, spatial_domain=NULL, sppoly=NULL, lookup_mode="stmv" ) {
  # lookup from rawdata

  if (is.null(spatial_domain))  {
    pS = substrate_parameters(  project_class="core"  )
  } else {
    pS = substrate_parameters( spatial_domain=spatial_domain, project_class="core"  )
  }

  vnmod = pS$variabletomodel
  vnmod2 = paste(vnmod, "mean", sep="." )

  crs_lonlat =  st_crs(projection_proj4string("lonlat_wgs84"))

  LU = substrate_db ( p=pS, DS="aggregated_data" )  # raw data
  LU = LU[ which( LU$lon > pS$corners$lon[1] & LU$lon < pS$corners$lon[2]  & LU$lat > pS$corners$lat[1] & LU$lat < pS$corners$lat[2] ), ]
  LU = lonlat2planar(LU, proj.type=pS$aegis_proj4string_planar_km)

  M = lonlat2planar(M, proj.type=pS$aegis_proj4string_planar_km) # get planar projections of lon/lat in km
  M[,vnmod] = LU[ match(
    array_map( "xy->1", M[,c("plon","plat")], gridparams=pS$gridparams ),
    array_map( "xy->1", LU[,c("plon","plat")], gridparams=pS$gridparams )
  ), vnmod2 ]

  if (!is.null(sppoly)) {
    # if any still missing then use a mean depth by AUID
    ii = NULL
    ii =  which( !is.finite(M[,vnmod]))
    if (length(ii) > 0) {
      if (!exists("AUID", M)) {
        M_AUID = st_points_in_polygons(
          pts = st_as_sf( M[ii,], coords=c("lon","lat"), crs=crs_lonlat ),
          polys = sppoly[, "AUID"],
          varname = "AUID"
        )
        M_AUID = as.character( M_AUID )  # match each datum to an area
      }
      LU = lonlat2planar(LU, proj.type=pS$aegis_proj4string_planar_km)
      LU$AUID = st_points_in_polygons(
        pts = st_as_sf( LU, coords=c("lon","lat"), crs=crs_lonlat ),
        polys = sppoly[, "AUID"],
        varname="AUID"
      )
      LU = tapply( LU[, vnmod2], LU$AUID, FUN=median, na.rm=TRUE )
      jj = match( as.character( M_AUID ), as.character( names(LU )) )
      M[ii,vnmod] = LU[jj]
    }
  }

  if (lookup_mode %in% c("stmv",  "hybrid") ) {
      # if any still missing then use stmv depths
      pC = substrate_parameters( spatial_domain=pS$spatial_domain, project_class=lookup_mode  )
      ii = NULL
      ii =  which( !is.finite( M[,vnmod] ))
      if (length(ii) > 0) {
        LU = substrate_db ( pC, DS="complete", varnames="all" )  # raw data
        LU = planar2lonlat(LU, proj.type=pC$aegis_proj4string_planar_km)
        LU = LU[ which( LU$lon > pC$corners$lon[1] & LU$lon < pC$corners$lon[2]  & LU$lat > pC$corners$lat[1] & LU$lat < pC$corners$lat[2] ), ]
        M[ii,vnmod] = LU[ match(
          array_map( "xy->1", M[ii, c("plon","plat")], gridparams=pC$gridparams ),
          array_map( "xy->1", LU[,c("plon","plat")], gridparams=pC$gridparams )
        ), vnmod ]
      }
  }

  return( M[,vnmod] )

}
