
substrate.figures = function( p=NULL, varnames="substrate.grainsize", datarange=NULL, logyvar=FALSE, isodepths = c( 100, 300, 500 ), savetofile="png", width=1365, height=1024, pointsize=12, res=96, quality=80 ) {


  #  if ( !exists("data_root", p) ) p$data_root = project.datadirectory( "aegis", p$project.name )
  #  if ( !exists("datadir", p) )   p$datadir  = file.path( p$data_root, "data" )
  #  if ( !exists("modeldir", p) )  p$modeldir = file.path( p$data_root, "modelled" )

  z = bathymetry.db(spatial.domain=p$spatial.domain, DS="baseline")
  b = substrate.db( p=p, DS="complete" )
#   lattice::levelplot( log(o$substrate.grainsize) ~ plon +plat, data=b, aspect="iso")

  for (vn in varnames) {
    print(vn)

    # oc = landmask( db="worldHires", regions=c("Canada", "US"), return.value="not.land", tag="predictions",internal.crs=p$internal.crs )
    # if (length(oc) > 0) {
    #   b = b[oc,]
    # }
    if (logyvar) b[,vn] = log( b[,vn] )
    datarange = quantile( b[,vn], probs=c(0.001, 0.999), na.rm=TRUE )
  #  if (is.null( datarange)) datarange = range( b[,vn], na.rm=TRUE )
    dr = seq( datarange[1], datarange[2], length.out=100)

    levplt = levelplot( b[,vn] ~ plon + plat, data=z, aspect="iso", main=NULL,
      at=dr, col.regions=rev(color.code( "seis", dr)) ,
      contour=FALSE, labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE),
      panel = function(x, y, subscripts, ...) {
        panel.levelplot (x, y, subscripts, aspect="iso", rez=c(1,1), ...)
        # sp.lines( isobath.db( p=p, DS="isobath", depths=isodepths, crs=p$internal.crs ), col = "gray80", cex=0.1 )
        sp.lines( coastline.db( p=p, crs=p$internal.crs ), col = "steelblue", cex=0.1 )
      }
    )

    if ( savetofile != "" ) {
      outdir = file.path( p$data_root, "maps", p$spatial.domain )
      for (i in 1:length(savetofile)){
        devtype = savetofile[i]
        if (devtype =="jpeg") devtype="jpg"
        fn = file.path( outdir, paste( "substrate", vn, p$spatial.domain, devtype, sep=".") )
        print(fn)
        if (devtype == "pdf" ) {
          pdf(file=fn, width=5, height=4, bg='white')
        } else if (devtype == "png" ) {
          png(filename=fn, width=width, height=height, pointsize=pointsize, res=res, bg='white' )
        } else if (devtype == "jpg" ) {
          jpeg(filename=fn, width=width, height=height, pointsize=pointsize, res=res, bg='white', quality=quality )
        } else {
          stop( "device not supported ... add it here")
        }
        print( levplt )
        dev.off()
      }
    } else {
      print( levplt )
    }

  }
  return(varnames)
}
