

  substrate_db = function( p=NULL, DS=NULL, varnames=NULL, sppoly=NULL, redo=FALSE, ... ) {

    if ( is.null(p))  {
      p_add = list(...)
      if (length(p_add) > 0 ) {
        p = substrate_parameters(...)
      } else {
        p = substrate_parameters()
      }
    }


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
      message("FIX ME::: sp, maptools deprecated, use sf/stars")
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
      substrate = substrate_db( DS="substrate.initial" )
      substrate = as.data.frame( substrate )
      names(substrate) = c("grainsize", "plon", "plat" )
      substrate = substrate[,c("plon", "plat", "grainsize")]
      proj4.params = "+proj=utm +zone=20 +ellps=GRS80 +datum=NAD83 +units=m"  # original/raw data still in NAD83 geoid
      substrate= planar2lonlat ( substrate, proj4.params )
      substrate= substrate[ ,c("lon", "lat", "grainsize")]
      substrate= lonlat2planar ( substrate, p$aegis_proj4string_planar_km )
      save( substrate, file=filename, compress=TRUE   )
      return ( filename )
    }


    # ---------------------------------------



    if ( DS=="aggregated_data") {

      fn = file.path( p$datadir, paste( "substrate", "aggregated_data", round(p$inputdata_spatial_discretization_planar_km, 6) , "rdata", sep=".") )
      if (!redo)  {
        if (file.exists(fn)) {
          load( fn)
          return( M )
        }
      }

      M = substrate_db( p=p, DS="lonlat.highres" )
      setDT(M)
      
      # p$quantile_bounds = c(0.0005, 0.9995)
      if (exists("quantile_bounds", p)) {
        vn = "grainsize"
        TR = quantile(M[[vn]], probs=p$quantile_bounds, na.rm=TRUE )
        oo = which( M[[vn]] < TR[1])
        if (length(oo) > 0) M[[vn]][oo] = TR[1]
        oo = which( M[[vn]] > TR[2])
        if (length(oo) > 0) M[[vn]][oo] = TR[2]
      }

      M = lonlat2planar( M, p$aegis_proj4string_planar_km )  # ensure the use of correct projection

      M$plon = trunc(M$plon / p$inputdata_spatial_discretization_planar_km + 1 ) * p$inputdata_spatial_discretization_planar_km
      M$plat = trunc(M$plat / p$inputdata_spatial_discretization_planar_km + 1 ) * p$inputdata_spatial_discretization_planar_km

      M = M[, .(mean=mean(grainsize, na.rm=TRUE), sd=sd(grainsize, na.rm=TRUE), n=length(which(is.finite(grainsize))) ), by=list(plon, plat) ]

      colnames(M) = c( "plon", "plat", paste( p$variabletomodel, c("mean", "sd", "n"), sep=".") )
      M = planar2lonlat( M, p$aegis_proj4string_planar_km )

      setDF(M)
   
      attr( M, "proj4string_planar" ) =  p$aegis_proj4string_planar_km
      attr( M, "proj4string_lonlat" ) =  projection_proj4string("lonlat_wgs84")

      save(M, file=fn, compress=TRUE)

      return( M )
    }


   # ------------------------------

    if ( DS=="areal_units_input" ) {
      
      outdir = file.path( p$data_root, "modelled", p$carstm_model_label ) 
      fn = file.path( outdir, "areal_units_input.rdata"  )
      if ( !file.exists(outdir)) dir.create( outdir, recursive=TRUE, showWarnings=FALSE )

      xydata = NULL
      if (!redo)  {
        if (file.exists(fn)) {
          load( fn)
          return( xydata )
        }
      }
      xydata = substrate_db( p=p, DS="aggregated_data"   )  #
      names(xydata)[which(names(xydata)=="substrate.grainsize.mean" )] = "substrate.grainsize"
      xydata = xydata[ geo_subset( spatial_domain=p$spatial_domain, Z=xydata ) , ] # need to be careful with extrapolation ...  filter depths

      xydata = xydata[ , c("lon", "lat"  )]

      save(xydata, file=fn, compress=TRUE )
      return( xydata )
    }

    # -----------------------


    if ( DS=="carstm_inputs") {

      # prediction surface
      crs_lonlat = st_crs(projection_proj4string("lonlat_wgs84"))
      if (is.null(sppoly)) sppoly = areal_units( p=p )  # will redo if not found
      sppoly = st_transform(sppoly, crs=crs_lonlat )
      areal_units_fn = attributes(sppoly)[["areal_units_fn"]]

      fn = file.path( p$modeldir, p$carstm_model_label, paste("carstm_inputs", areal_units_fn, sep="_") )
      if (p$carstm_inputs_prefilter =="rawdata") {
        fn = file.path( p$modeldir, p$carstm_model_label, paste("carstm_inputs_rawdata", areal_units_fn, sep="_") )
      }

      # inputs are shared across various secneario using the same polys
      #.. store at the modeldir level as default
      outputdir = dirname( fn )
      if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

      if (!redo)  {
        if (file.exists(fn)) {
          load( fn)
          return( M )
        }
      }

      if (p$carstm_inputs_prefilter=="aggregated") {
        M = substrate_db ( p=p, DS="aggregated_data" )  # 16 GB in RAM just to store!
        names(M)[which(names(M)==paste(p$variabletomodel, "mean", sep=".") )] = p$variabletomodel
        setDT(M)

      } else  {

        M = substrate_db( p=p, DS="lonlat.highres" )
        setDT(M)
      
        names(M)[which(names(M)=="grainsize" )] = p$variabletomodel
        M = M[ which( !duplicated(M)), ]
        M = M[ which( M$lon > p$corners$lon[1] & M$lon < p$corners$lon[2]  & M$lat > p$corners$lat[1] & M$lat < p$corners$lat[2] ), ]
        M = lonlat2planar( M, proj.type=p$aegis_proj4string_planar_km)  # first ensure correct projection

      }
    
      if (p$carstm_inputs_prefilter =="sampled") {
        M$plon = trunc(M$plon / p$inputdata_spatial_discretization_planar_km + 1 ) * p$inputdata_spatial_discretization_planar_km
        M$plat = trunc(M$plat / p$inputdata_spatial_discretization_planar_km + 1 ) * p$inputdata_spatial_discretization_planar_km
        M = M[,.SD[sample(.N, min(.N, p$carstm_inputs_prefilter_n))], by =list(plon, plat) ]  # compact, might be slightly slower
        # M = M[ M[, sample(.N, min(.N, p$carstm_inputs_prefilter_n) ), by=list(plon, plat)], .SD[i.V1], on=list(plon, plat), by=.EACHI]  # faster .. just a bit
       }

      if (p$carstm_inputs_prefilter != "aggregated") {
        # aggregated data already have been truncated ... 
        if (exists("quantile_bounds", p)) {
          TR = quantile(M[[p$variabletomodel]], probs=p$quantile_bounds, na.rm=TRUE ) # this was -1.7, 21.8 in 2015
          keep = which( M[[p$variabletomodel]] >=  TR[1] & M[[p$variabletomodel]] <=  TR[2] )
          if (length(keep) > 0 ) M = M[ keep, ]
          # this was -1.7, 21.8 in 2015
        }
      }

      # levelplot(substrate.grainsize.mean~plon+plat, data=M, aspect="iso")
      M = carstm_prepare_inputdata( p=p, M=M, sppoly=sppoly )

      # remove duplicated data .. helps with variance estimation
#      i = duplicated(M)
      M$space = match( M$AUID, sppoly$AUID )  # for bym/car require index matching nb graph

#      if (length(i) >0) M = M[-i, ] 
      attr( M, "proj4string_planar" ) =  p$aegis_proj4string_planar_km
      attr( M, "proj4string_lonlat" ) =  projection_proj4string("lonlat_wgs84")

      save( M, file=fn, compress=TRUE )
      return( M )
    }

    # ---------------------------------------

    if ( DS=="stmv_inputs") {

      varstokeep = unique( c( p$stmv_variables$Y, p$stmv_variables$LOCS, p$stmv_variables$COV ) )
      pB = bathymetry_parameters( spatial_domain=p$spatial_domain, project_class="default"  )  # default is the "best" performing method
      B = bathymetry_db( p=pB, DS="baseline", varnames=varstokeep )

      # range checks
      if (exists("z", B))  B$z[ which( B$z < 0.5)]  = 0.5 # meters
      if (exists("dZ", B))  B$dZ[ which( B$dZ < 0.001)] = 0.001 # will be log transformed .. range check
      if (exists("ddZ", B))  B$ddZ[ which( B$ddZ < 0.001)] = 0.001
      if (exists("b.sdTotal", B))  B$b.sdTotal[ which(!is.finite(B$b.sdTotal))] = median( B$b.sdTotal, na.rm=TRUE )

      bid = array_map( "xy->1", B[,c("plon", "plat")], gridparams=p$gridparams )

      S = substrate_db( p=p, DS="aggregated_data" )
      names(S)[which(names(S) == paste(p$variabletomodel, "mean", sep="."))] = p$variabletomodel

      # merge covars into S
      sid = array_map( "xy->1", S[,c("plon", "plat")], gridparams=p$gridparams )
      u = match( sid, bid )
      B_matched = B[u, ]
      B_matched$plon = B_matched$plat = NULL

      S = cbind( S, B_matched )
      S = S[ is.finite( S$substrate.grainsize ), ]

      OUT  = list( LOCS=B[, p$stmv_variables$LOCS], COV=B[, p$stmv_variables$COV ] )

      return(  list( input=S, output=OUT ) )

    }

    #-------------------------


    if ( DS %in% c( "complete", "complete.redo" ) ) {
     #// substrate_db( DS="complete" .. ) returns the final form of the substrate data after
     #// regridding and selection to area of interest as specificied by girds.new=c("SSE", etc)

      outdir = file.path( p$modeldir, p$stmv_model_label, p$project_class, paste(  p$stmv_global_modelengine, p$stmv_local_modelengine, sep="_") )

      fn = file.path( outdir, paste( "substrate", "complete", p$spatial_domain, "rdata", sep=".") )

      if ( DS %in% c("complete") ) {

        defaultdir = project.datadirectory( "aegis", "substrate", "modelled" )
        if (!file.exists(fn)) fn = file.path( defaultdir, paste( "substrate", "complete", p$spatial_domain, "rdata", sep=".") )

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
      names(S) = paste( p$variabletomodel, c( "", ".lb", ".ub") , sep="")
      rm (Smean, Slb, Sub); gc()


      if (0) {
        pB = bathymetry_parameters( spatial_domain=p$spatial_domain, project_class="default"  )  # default is the best performing method
        B = bathymetry_db(p=pB, DS="baseline")
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
      p0 = bathymetry_parameters( spatial_domain=p0$spatial_domain, project_class="default"  )  # default is the best performing method
      L0 = bathymetry_db( p=p0, DS="baseline" )
      L0i = array_map( "xy->2", L0, gridparams=p0$gridparams )

      varnames = setdiff( names(S0), c("plon","plat", "lon", "lat") )
      #using fields
      grids = setdiff( unique( p0$spatial_domain_subareas ), p0$spatial_domain )
      for (gr in grids ) {
        print(gr)
        #target projection
        p1 = bathymetry_parameters( spatial_domain=gr, project_class="default"  )  # default is the best performing method
        L1 = bathymetry_db( p=p1, DS="baseline" )
        L1i = array_map( "xy->2", L1[, c("plon", "plat")], gridparams=p1$gridparams )
        L1 = planar2lonlat( L1, proj.type=p1$aegis_proj4string_planar_km )
        S = L1
        L1$plon_1 = L1$plon # store original coords
        L1$plat_1 = L1$plat
        L1 = lonlat2planar( L1, proj.type=p0$aegis_proj4string_planar_km )
        p1$wght = fields::setup.image.smooth(
          nrow=p1$nplons, ncol=p1$nplats, dx=p1$pres, dy=p1$pres,
          theta=p1$pres/3, xwidth=4*p1$pres, ywidth=4*p1$pres )
        # theta=p1$pres/3 assume at pres most of variance is accounted ... correct if dense pre-intepolated matrices .. if not can be noisy

        for (vn in varnames) {
          S[,vn] = spatial_warp( S0[,vn], L0, L1, p0, p1, "fast", L0i, L1i )
        }
        S = S[, names(S0)]

        # range checks
        ii = which( S[,p$variabletomodel] < exp(-6) )
        if (length(ii) > 0 ) S[,p$variabletomodel][ii] = exp(-6)

        ii = which( S[,p$variabletomodel] > exp(5) )
        if (length(ii) > 0 ) S[,p$variabletomodel][ ii ] = exp(5)

        fn = file.path( outdir, paste( "substrate", "complete", p1$spatial_domain, "rdata", sep=".") )
        save (S, file=fn, compress=TRUE)
      }

      if(0){
        datarange = quantile( S[,p$variabletomodel], probs=c(0.01, 0.99), na.rm=TRUE )
        dr = log(seq( datarange[1], datarange[2], length.out=100))
        locations = bathymetry_db(p=bathymetry_parameters( spatial_domain=p$spatial_domain, project_class="default"  ), DS="baseline")
        S$plon = locations$plon
        S$plat = locations$plat
        levelplot(log( eval(p$variabletomodel) ~ plon+plat, S,  at=dr, col.regions=(color.code( "seis", dr))))
      }
      return ( "Completed subsets" )
    }


    # ---------------------------------------


    if (DS=="maps") {

      if (is.null(varnames)) varnames=p$variabletomodel
      datarange=NULL
      logyvar=FALSE
      isodepths = c( 100, 300, 500 )

      b = substrate_db( p=p, DS="complete" )
      if (logyvar) {
        b = b[ b[,varnames]>0, ]
        b[,varnames] = log(b[,varnames])
      }

      # oc = landmask( db="worldHires", regions=c("Canada", "US"), return.value="not.land", tag="predictions", crs=p$aegis_proj4string_planar_km )
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
            sp.lines( as( isobath_db( DS="isobath", depths=isodepths, project_to=p$aegis_proj4string_planar_km ), "Spatial"), col = "gray80", cex=0.1 )
            sp.lines( as( aegis.coastline::coastline_db( p=p, project_to=p$aegis_proj4string_planar_km), "Spatial" ), , col = "steelblue", cex=0.1)
          }
        )
      )
    }

  }
