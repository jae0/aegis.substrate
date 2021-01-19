substrate_lookup_rawdata = function( lonlat, spatial_domain=NULL, sppoly=NULL, lookup_mode="stmv" ) {
  # lookup from rawdata

  if (is.null(spatial_domain))  {
    pS = substrate_parameters(  project_class="core"  ) 
  } else {
    pS = substrate_parameters( spatial_domain=pS$spatial_domain, project_class="core"  )
  }

  vnmod = pS$variabletomodel
  crs_lonlat =  st_crs(projection_proj4string("lonlat_wgs84"))

  LU = substrate_db ( p=pS, DS="aggregated_data" )  # raw data
  LU = LU[ which( LU$lon > pS$corners$lon[1] & LU$lon < pS$corners$lon[2]  & LU$lat > pS$corners$lat[1] & LU$lat < pS$corners$lat[2] ), ]
  LU = lonlat2planar(LU, proj.type=pS$aegis_proj4string_planar_km)

  lonlat = lonlat2planar(lonlat, proj.type=pS$aegis_proj4string_planar_km) # get planar projections of lon/lat in km
  lonlat$Z = LU[ match( 
    array_map( "xy->1", lonlat[,c("plon","plat")], gridparams=pS$gridparams ), 
    array_map( "xy->1", LU[,c("plon","plat")], gridparams=pS$gridparams ) 
  ), paste(vnmod, "mean", sep=".")]

  if (!is.null(sppoly)) {
    # if any still missing then use a mean depth by AUID
    ii = NULL
    ii =  which( !is.finite(lonlat$Z))
    if (length(ii) > 0) {
      if (!exists("AUID", lonlat)) {
        lonlat$AUID = ""
        lonlat$AUID[ii] = st_points_in_polygons(
          pts = st_as_sf( lonlat[ii,], coords=c("lon","lat"), crs=crs_lonlat ),
          polys = sppoly[, "AUID"],
          varname = "AUID"
        )
        lonlat$AUID = as.character( lonlat$AUID )  # match each datum to an area
      }
      LU = lonlat2planar(LU, proj.type=pS$aegis_proj4string_planar_km)
      LU$AUID = st_points_in_polygons(
        pts = st_as_sf( LU, coords=c("lon","lat"), crs=crs_lonlat ),
        polys = sppoly[, "AUID"],
        varname="AUID"
      )
      LU = tapply( LU[, paste(vnmod, "mean", sep="." )], LU$AUID, FUN=median, na.rm=TRUE )
      jj = match( as.character( lonlat$AUID[ii]), as.character( names(LU )) )
      lonlat$Z[ii] = LU[jj]
    }
  }

  if (lookup_mode %in% c("stmv",  "hybrid") ) {
      # if any still missing then use stmv depths
      pC = substrate_parameters( spatial_domain=p$spatial_domain, project_class=lookup_mode  )
      ii = NULL
      ii =  which( !is.finite( lonlat$Z ))
      if (length(ii) > 0) {
        LU = substrate_db ( pC, DS="complete", varnames="all" )  # raw data
        LU = planar2lonlat(LU, proj.type=p$aegis_proj4string_planar_km)
        LU = LU[ which( LU$lon > p$corners$lon[1] & LU$lon < p$corners$lon[2]  & LU$lat > p$corners$lat[1] & LU$lat < p$corners$lat[2] ), ]
        lonlat$Z[ii] = LU[ match( 
          array_map( "xy->1", lonlat[ii, c("plon","plat")], gridparams=p$gridparams ), 
          array_map( "xy->1", LU[,c("plon","plat")], gridparams=p$gridparams ) 
        ), vnmod ]
      }
  }

  return( lonlat$Z )

}
