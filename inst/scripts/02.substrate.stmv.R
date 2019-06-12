
## NOTE:: substrate size is really only relevant for SSE/snowcrab domain right now as no
##        other data source has been found/identified
##        but working at the size of canada.east.highres for compatibility with bathymetry
## TODO:: add data collected by snow crab survey and any others for that matter


# about 1.5 hr
scale_ram_required_main_process = 1 # GB twostep / fft
scale_ram_required_per_process  = 1 # twostep / fft /fields vario ..  (mostly 0.5 GB, but up to 5 GB)
scale_ncpus = min( parallel::detectCores(), floor( (ram_local()- scale_ram_required_main_process) / scale_ram_required_per_process ) )

# nn hrs
interpolate_ram_required_main_process = 3.5 # GB twostep / fft
interpolate_ram_required_per_process  = 2  # twostep / fft /fields vario ..
interpolate_ncpus = min( parallel::detectCores(), floor( (ram_local()- interpolate_ram_required_main_process) / interpolate_ram_required_per_process ) )


p = aegis.substrate::substrate_parameters(
  project.mode="stmv",
  data_root = project.datadirectory( "aegis", "substrate" ),
  spatial.domain = "canada.east.highres" ,
  spatial.domain.subareas = c( "canada.east", "SSE", "snowcrab", "SSE.mpa" ),
  pres_discretization_substrate = 1 / 20, # 1==p$pres; controls resolution of data prior to modelling (km .. ie 20 linear units smaller than the final discretization pres)
  stmv_dimensionality="space",
  stmv_global_modelengine = "gam",
  stmv_global_modelformula = formula( paste(
    'substrate.grainsize ',
    ' ~ s( b.sdTotal, k=3, bs="ts") ',
    ' + s(log(z), k=3, bs="ts") + s(log(dZ), k=3, bs="ts") +s(log(ddZ), k=3, bs="ts") '
  ) ),
  stmv_global_family = gaussian(link="log"),
  # stmv_Y_transform =list(  # a log-normal works ok but a model of log-transformed data works better .. ie, working upon medians which is really OK
  #   transf = function(x) {log(x)} ,
  #   invers = function(x) {exp(x)}
  # ), # data range is from -1667 to 5467 m: make all positive valued
  stmv_local_modelengine="fft",  # currently the perferred approach
  stmv_lowpass_phi = 1*2, # p$res *2 = 1 *2:: FFT based method when operating gloablly
  stmv_lowpass_nu = 0.5, # this is exponential covar
  stmv_rsquared_threshold = 0.1, # lower threshold
  depth.filter = 0.1, # the depth covariate is input in m, so, choose stats locations with elevation > 0 m as being on land
  stmv_nmin = 1000, # stmv_nmin/stmv_nmax changes with resolution
  stmv_nmax = 4000, # numerical time/memory constraint -- anything larger takes too much time .. anything less .. errors
  stmv_distance_statsgrid = 4, # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
  stmv_distance_scale = c( 30, 40 ), # km ... approx guess of 95% AC range
  stmv_distance_prediction_fraction = 4/5, # i.e. 4/5 * 5 = 4 km
  stmv_clusters = list( scale=rep("localhost", scale_ncpus), interpolate=rep("localhost", interpolate_ncpus) )  # ncpus for each runmode
)

# Override stmv generics :
p$stmv_fft_filter = "lowpass_matern_tapered" #  act as a low pass filter first before matern .. depth has enough data for this. Otherwise, use:
p$stmv_range_correlation_fft_taper = 0.5  # in local smoothing convolutions occur of this correlation scale


# the scale /variogram is larger than 40, so not functional ...

stmv( p=p, runmode=c("globalmodel", "scale", "interpolate" ) ) # no global_model and force a clean restart


# as the interpolation process is so expensive, regrid based off the above run
substrate.db( p=p, DS="complete.redo" )


# quick map
b = bathymetry.db(p=p, DS="baseline")
o = substrate.db( p=p, DS="complete" )
lattice::levelplot( log(o$substrate.grainsize) ~ plon +plat, data=b, aspect="iso")


# or a cleaner map:
# p = aegis_parameters()
substrate.figures( p=p, varnames=c( "s.ndata" ), logyvar=FALSE, savetofile="png" )
substrate.figures( p=p, varnames=c( "substrate.grainsize", "s.nu", "s.range", "s.phi", "s.sdTotal", "s.sdSpatial", "s.sdObs"), logyvar=TRUE, savetofile="png" )


# to summarize just the global model
o = stmv_db( p=p, DS="global_model" )
summary(o)
plot(o)
AIC(o)  # [1] 813065

# Global model results:

Family: gaussian
Link function: log

Formula:
substrate.grainsize ~ s(b.sdTotal, k = 3, bs = "ts") + s(log(z),
    k = 3, bs = "ts") + s(log(dZ), k = 3, bs = "ts") + s(log(ddZ),
    k = 3, bs = "ts")

Parametric coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept)  -0.8787     0.0175   -50.3   <2e-16

Approximate significance of smooth terms:
             edf Ref.df      F p-value
s(b.sdTotal)   2      2 1595.6  <2e-16
s(log(z))      2      2 4663.4  <2e-16
s(log(dZ))     2      2   83.8  <2e-16
s(log(ddZ))    2      2  210.1  <2e-16

R-sq.(adj) =  0.157   Deviance explained = 15.4%
GCV = 5.5405  Scale est. = 5.5402    n = 179311
