
substrate_lookup = function( p, locs, vnames="substrate.grainsize", output_data_class="points", source_data_class="aggregated_rawdata", locs_proj4string="lonlat" ) {

  # deprecated ... too slow when multiple data loads are required and needless transformations
  # if locs is points, then need to send info on projection as an attribute proj4string"

  require(aegis.substrate)

  if (p$project_name != "substrate") {
    p = substrate_parameters(p=parameters_reset(p), project_name="substrate" )
    warning( "Parameter list may be inconsistent")
  }


  # load input data or reformat it
   if (source_data_class=="rawdata") {

      B = substrate_db ( p=p, DS="lonlat.highres" )  # 16 GB in RAM just to store!
#      Bnames = c("lon", "lat", "grainsize", "plon", "plat"),

   } else if (source_data_class=="aggregated_rawdata") {

      B = substrate_db ( p=p, DS="aggregated_data" )
#       Bnames = c("substrate.grainsize.mean", "substrate.grainsize.sd",  "substrate.grainsize.n", "plon", "plat", "lon", "lat")
      B$substrate.grainsize = B$substrate.grainsize.mean
      B$substrate.grainsize.mean  = NULL

   } else if (source_data_class=="stmv") {

      B = substrate_db(p=p, DS="complete", varnames="all" )
    # Bnames = c( "plon", "plat", "substrate.grainsize", "substrate.grainsize.lb", "substrate.grainsize.ub",
    #   "s.sdTotal", "s.rsquared", "s.ndata", "s.sdSpatial", "s.sdObs", "s.phi", "s.nu", "s.localrange" )
      zname = "substrate.grainsize"

   } else if (source_data_class=="carstm") {

      Bcarstm = carstm_summary( p=p ) # to load currently saved sppoly
      B = areal_units( p=p )
      bm = match( B$AUID, Bcarstm$AUID )
      B$substrate.grainsize  = Bcarstm$substrate.grainsize.predicted[ bm ]
      B$substrate.grainsize.se = Bcarstm$substrate.grainsize.predicted_se[ bm ]
      Bcarstm = NULL
      zname = "substrate.grainsize"
  }

  Bnames = setdiff( names(B), c("AUID", "uid", "layer", "plon", "plat", "lon", "lat", "au_sa_km2",
    "cfanorth_surfacearea", "cfasouth_surfacearea", "cfa23_surfacearea",  "cfa24_surfacearea", "cfa4x_surfacearea" ) )


  if (output_data_class == "points ") {

    if ( source_data_class %in% c("rawdata", "aggregated_rawdata", "stmv" ) )  {

      if ( !is.null( attr( locs, "proj4string" )) ) locs_proj4string = attr( locs, "proj4string" )
      if ( locs_proj4string =="lonlat" ) {
        names( locs) = c("lon", "lat")
        locs = lonlat2planar( locs[, c("lon", "lat")], proj.type=p$aegis_proj4string_planar_km )
        locs_proj4string = p$aegis_proj4string_planar_km
      }
      if ( locs_proj4string != p$aegis_proj4string_planar_km ) {
        locs = planar2lonlat( locs[, c("plon", "plat")], proj.type=locs_proj4string )
        locs = lonlat2planar( locs[, c("lon", "lat")], proj.type=p$aegis_proj4string_planar_km )
        locs_proj4string = p$aegis_proj4string_planar_km
      }
      B_map = array_map( "xy->1", B[,c("plon","plat")], gridparams=p$gridparams )
      locs_map = array_map( "xy->1", locs[,c("plon","plat")], gridparams=p$gridparams )
      locs_index = match( locs_map, B_map )
      vnames = intersect( names(B), vnames )
      if ( length(vnames) ==0 ) vnames=names(B) # no match returns all
      return( B[locs_index, vnames] )
    }

    if ( source_data_class=="carstm") {
      # convert to raster then match
      require(raster)
      raster_template = raster(extent(locs))
      res(raster_template) = p$areal_units_resolution_km  # crs usually in meters, but aegis's crs is in km
      crs(raster_template) = projection(locs) # transfer the coordinate system to the raster

      locs = sf::st_as_sf( as.data.frame(locs), coords=c(1, 2) )
      st_crs(locs) = crs(B)
      for (vn in Bnames) {
        Bf = fasterize::fasterize( as(B, "sf"), raster_template, field=vn )
        vn2 = paste(vn, "sd", sep="." )
        locs[, vn ] = raster::extract( Bf, locs, fun=mean, na.rm=TRUE)
        locs[, vn2] = raster::extract( Bf, locs, fun=sd, na.rm=TRUE)
      }
      vnames = intersect( names(B), vnames )
      if ( length(vnames) ==0 ) vnames=names(B) # no match returns all
      return( as.matrix(locs[[vnames]]) )
    }
  }


  if ( output_data_class=="areal_units") {

    # expects loc to be a spatial polygon data frame

    if ( source_data_class %in% c("rawdata", "aggregated_rawdata", "stmv" ) ) {
      Bsf = sf::st_as_sf( B, coords=c("lon", "lat") )
      st_crs(Bsf) = CRS( projection_proj4string("lonlat_wgs84") )
      Bsf = sf::st_transform( Bsf, crs=CRS(proj4string(locs)) )
      for (vn in Bnames) {
        vn2 = paste(vn, "sd", sep="." )
        #Bf= ...
        locs[,vn] = st_polygons_in_polygons( locs, Bf[,vn], fn=mean, na.rm=TRUE )
        locs[,vn2] = st_polygons_in_polygons( locs, Bf[,vn], fn=sd, na.rm=TRUE )
      }
      vnames = intersect( names(B), vnames )
      if ( length(vnames) ==0 ) vnames=names(B) # no match returns all
      return(locs[, vnames] )
    }


    if ( source_data_class=="carstm") {
      # convert to raster then match
      require(raster)
      raster_template = raster(extent(locs)) # +1 to increase the area
      res(raster_template) = p$areal_units_resolution_km  # crs usually in meters, but aegis's crs is in km
      crs(raster_template) = projection(locs) # transfer the coordinate system to the raster
      Bsf = sf::st_transform( as(B, "sf"), crs=CRS(proj4string(locs)) )  # B is a carstm sppoly
      for (vn in Bnames) {
        Bf = fasterize::fasterize( Bsf, raster_template, field=vn )
        vn2 = paste(vn, "sd", sep="." )
        locs[,vn] = st_polygons_in_polygons( locs, Bf[,vn], fn=mean, na.rm=TRUE )
        locs[,vn2] = st_polygons_in_polygons( locs, Bf[,vn], fn=sd, na.rm=TRUE )
      }
      vnames = intersect( names(B), vnames )
      if ( length(vnames) ==0 ) vnames=names(B) # no match returns all
      return(locs[,vnames])
   }
  }

}


