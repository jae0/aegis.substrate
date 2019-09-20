

  substrate.db = function( p=substrate_parameters(), DS=NULL, varnames=NULL ) {

  #  if ( !exists("data_root", p) ) p$data_root = project.datadirectory( "aegis", p$project.name )
  #  if ( !exists("datadir", p) )   p$datadir  = file.path( p$data_root, "data" )
  #  if ( !exists("modeldir", p) )  p$modeldir = file.path( p$data_root, "modelled" )


    if ( DS %in% c("substrate.initial", "substrate.initial.redo") ) {
      # Read in the ArcInfo ascii grid file using library maptools and output a SpatialGridDataFrame
      # data provided by Kostelev:
      # Kostylev, V.E., and Hannah, C.G., 2007, Process-driven characterization and mapping of seabed habitats,
      # in Todd, B.J.,and Greene, H.G., eds., Mapping the Seafloor for Habitat Characterization:
      # Geological Association of Canada, Special Paper 47, p. 171-184.
      # Scotian shelf gridded grain size (mm).
      # NAD83 UTM zone 20 (I think)

      rawdata.file = file.path( p$datadir, "grainsize.txt" )
      filename = file.path( p$datadir, "substrate.asciigrid.rdata" )

      if (DS =="substrate.initial" ) {
        load( filename )
        return ( substrate )
      }
      proj4.params = "+proj=utm +zone=20 +ellps=GRS80  +datum=NAD83 +units=m" #resolution is 500m X 500m
      substrate = sp::read.asciigrid( rawdata.file, proj4string=CRS( proj4.params ), colname="grainsize" )  ## mm
      save( substrate, file=filename, compress=TRUE )
      return(filename)
    }


    # ---------------------------------------

    # lon - lat converted
    if (  DS %in% c("lonlat.highres", "lonlat.highres.redo") ) {
      filename = file.path( p$datadir, "substrate.lonlat.highres.rdata" )
      if (DS =="lonlat.highres" ) {
        load( filename)
        return( substrate)
      }
      # initial data stored in planar coords ... convert to lon/lats
      substrate = substrate.db( DS="substrate.initial" )
      substrate = as.data.frame( substrate )
      names(substrate) = c("grainsize", "plon", "plat" )
      substrate = substrate[,c("plon", "plat", "grainsize")]
      proj4.params = "+proj=utm +zone=20 +ellps=GRS80 +datum=NAD83 +units=m"  # original/raw data still in NAD83 geoid
      substrate= planar2lonlat ( substrate, proj4.params )
      substrate= substrate[ ,c("lon", "lat", "grainsize")]
      substrate= lonlat2planar ( substrate, p$internal.crs )
      save( substrate, file=filename, compress=TRUE   )
      return ( filename )
    }


    # ---------------------------------------

    if ( DS=="stmv.inputs") {

      varstokeep = unique( c( p$variables$Y, p$variables$LOCS, p$variables$COV ) )
      B = bathymetry.db( p=p, DS="baseline", varnames=varstokeep )

      # range checks
      if (exists("z", B))  B$z[ which( B$z < 0.5)]  = 0.5 # meters
      if (exists("dZ", B))  B$dZ[ which( B$dZ < 0.001)] = 0.001 # will be log transformed .. range check
      if (exists("ddZ", B))  B$ddZ[ which( B$ddZ < 0.001)] = 0.001
      if (exists("b.sdTotal", B))  B$b.sdTotal[ which(!is.finite(B$b.sdTotal))] = median( B$b.sdTotal, na.rm=TRUE )

      bid = stmv::array_map( "xy->1", B[,c("plon", "plat")], gridparams=p$gridparams )

      S = substrate.db( p=p, DS="lonlat.highres" )
      S = lonlat2planar( S,  proj.type=p$internal.crs )  # utm20, WGS84 (snowcrab geoid)
      S$substrate.grainsize = S$grainsize
      S = S[ ,c("plon", "plat", "substrate.grainsize" )]

      # discretize to speed up the rest
      S$plon = round(S$plon/p$inputdata_spatial_discretization_planar_km + 1) * p$inputdata_spatial_discretization_planar_km
      S$plat = round(S$plat/p$inputdata_spatial_discretization_planar_km + 1) * p$inputdata_spatial_discretization_planar_km

      gsrez = 0.001
      oo = paste( S$plon, S$plat, round(S$substrate.grainsize/gsrez +1 )*gsrez )

      S = S[!duplicated(oo),]  # drop all duplicated data
      rm(oo); gc()

      # clean up duplicated locations with variable data
      S = tapply( X=S$substrate.grainsize, INDEX=list(S$plon, S$plat),
          FUN = function(w) {median(w, na.rm=TRUE)},
          simplify=TRUE )
      S = as.data.frame( as.table (S) )
      S[,1] = as.numeric(as.character( S[,1] ))
      S[,2] = as.numeric(as.character( S[,2] ))
      S = S[ which( is.finite( S[,3] )) ,]
      names(S) = c("plon", "plat", "substrate.grainsize")  # geometric means ..

      # merge covars into S
      sid = stmv::array_map( "xy->1", S[,c("plon", "plat")], gridparams=p$gridparams )
      u = match( sid, bid )
      B_matched = B[u, ]
      B_matched$plon = B_matched$plat = NULL
      S = cbind(S, B_matched )
      S = S[ is.finite( rowSums(S) ), ]


      OUT  = list( LOCS=B[, p$variables$LOCS], COV=B[, p$variables$COV ] )

      return(  list( input=S, output=OUT ) )

    }

    #-------------------------


    if ( DS %in% c( "complete", "complete.redo" ) ) {
     #// substrate.db( DS="complete" .. ) returns the final form of the substrate data after
     #// regridding and selection to area of interest as specificied by girds.new=c("SSE", etc)

      fn = file.path( p$modeldir, paste( "substrate", "complete", p$spatial.domain, "rdata", sep=".") )

      if ( DS %in% c("complete") ) {

        defaultdir = project.datadirectory( "aegis", "substrate", "modelled" )
        if (!file.exists(fn)) fn = file.path( defaultdir, paste( "substrate", "complete", p$spatial.domain, "rdata", sep=".") )

        S = NULL
        if ( file.exists ( fn) ) load( fn)
        Snames = names(S)
        if (is.null(varnames)) {
          varnames=Snames
        } else {
          varnames = intersect( Snames, varnames )
        }
        if (length(varnames) == 0) varnames=Snames  # no match .. send all
        S = S[ , varnames]
        return( S )
      }

      Smean = stmv_db( p=p, DS="stmv.prediction", ret="mean" )
      # leave ranges high to reflect high uncertainty
      Slb = stmv_db( p=p, DS="stmv.prediction", ret="lb" )
      Sub = stmv_db( p=p, DS="stmv.prediction", ret="ub" )
      S = as.data.frame( cbind( Smean, Sub, Slb) )  # ub and lb swap due to log scaling and link space operation
      names(S) = c( "substrate.grainsize", "substrate.grainsize.lb", "substrate.grainsize.ub")
      rm (Smean, Slb, Sub); gc()


      if (0) {
        B = bathymetry.db(p=p, DS="baseline")
        levelplot( (S[,1]) ~ plon + plat, B, aspect="iso", labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
        levelplot( log(S[,1]) ~ plon + plat, B, aspect="iso", labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
      }

      # merge into statistics
      SS = stmv_db( p=p, DS="stmv.stats" )
      colnames(SS) = paste("s", colnames(SS), sep=".")
      S = cbind( S, SS )

      save (S, file=fn, compress=TRUE)

      p0 = p  # the originating parameters
      S0 = S
      L0 = bathymetry.db( p=p, DS="baseline" )
      L0i = array_map( "xy->2", L0, gridparams=p0$gridparams )

      varnames = setdiff( names(S0), c("plon","plat", "lon", "lat") )
      #using fields
      grids = setdiff( unique( p0$spatial.domain.subareas ), p0$spatial.domain )
      for (gr in grids ) {
        print(gr)
        p1 = spatial_parameters( spatial.domain=gr ) #target projection
        L1 = bathymetry.db( p=p1, DS="baseline" )
        L1i = array_map( "xy->2", L1[, c("plon", "plat")], gridparams=p1$gridparams )
        L1 = planar2lonlat( L1, proj.type=p1$internal.crs )
        S = L1
        L1$plon_1 = L1$plon # store original coords
        L1$plat_1 = L1$plat
        L1 = lonlat2planar( L1, proj.type=p0$internal.crs )
        p1$wght = fields::setup.image.smooth(
          nrow=p1$nplons, ncol=p1$nplats, dx=p1$pres, dy=p1$pres,
          theta=p1$pres/3, xwidth=4*p1$pres, ywidth=4*p1$pres )
        # theta=p1$pres/3 assume at pres most of variance is accounted ... correct if dense pre-intepolated matrices .. if not can be noisy

        for (vn in varnames) {
          S[,vn] = spatial_warp( S0[,vn], L0, L1, p0, p1, "fast", L0i, L1i )
        }
        S = S[, names(S0)]

        # range checks
        ii = which( S$substrate.grainsize < exp(-6) )
        if (length(ii) > 0 ) S$substrate.grainsize[ii] = exp(-6)

        ii = which( S$substrate.grainsize > exp(5) )
        if (length(ii) > 0 ) S$substrate.grainsize[ ii ] = exp(5)

        fn = file.path( p$modeldir, paste( "substrate", "complete", p1$spatial.domain, "rdata", sep=".") )
        save (S, file=fn, compress=TRUE)
      }

      if(0){
        datarange = quantile( S$substrate.grainsize, probs=c(0.01, 0.99), na.rm=TRUE )
        dr = log(seq( datarange[1], datarange[2], length.out=100))
        levelplot(log(substrate.grainsize) ~ plon+plat, S,  at=dr, col.regions=(color.code( "seis", dr)))

      }
      return ( "Completed subsets" )
    }


    # ---------------------------------------


    if (DS=="maps") {

      if (is.null(varnames)) varnames="substrate.grainsize"
      datarange=NULL
      logyvar=FALSE
      isodepths = c( 100, 300, 500 )

      b = substrate.db( p=p, DS="complete" )
      if (logyvar) {
        b = b[ b[,varnames]>0, ]
        b[,varnames] = log(b[,varnames])
      }

      # oc = landmask( db="worldHires", regions=c("Canada", "US"), return.value="not.land", tag="predictions",internal.crs=p$internal.crs )
      # if (length(oc) > 0) {
      #   b = b[oc,]
      # }

      if (is.null( datarange)) datarange = quantile( b[,varnames], probs=c(0.05, 0.95), na.rm=TRUE )
      dr = seq( datarange[1], datarange[2], length.out=100)

      print(
        levelplot( b[,varnames] ~ plon + plat, data=b[,], aspect="iso", main=NULL,
          at=dr, col.regions=rev(color.code( "seis", dr)) ,
          contour=FALSE, labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE),
          panel = function(x, y, subscripts, ...) {
            panel.levelplot (x, y, subscripts, aspect="iso", rez=c(1,1), ...)
            sp.lines( isobath.db( p=p, DS="isobath", depths=isodepths, crs=p$internal.crs ), col = "gray80", cex=0.1 )
            sp.lines( aegis.coastline::coastline.db( p=p, crs=p$internal.crs, DS="gshhg coastline highres" ), col = "steelblue", cex=0.1 )
          }
        )
      )
    }

  }
