
# construct basic parameter list defining the main characteristics of the study
# and some plotting parameters (bounding box, projection, bathymetry layout, coastline)

  p = aegis.substrate::substrate_parameters( project_class="carstm", areal_units_resolution_km=5  )  #10 km grid gives about 2000 au .. in the right range for optimal solutions that are not too slow, 5 km gives stable results .. 


    # adjust based upon RAM requirements and ncores
    # inla.setOption( num.threads="6:2" )
    # inla.setOption(blas.num.threads= 3 )

    if (0) {
      p$fraction_todrop = 1/4 # aggressiveness of solution finding ( fraction of counts to drop each iteration)
      p$fraction_cv = 1.0  #sd/mean no.
      p$fraction_good_bad = 0.9
      p$areal_units_constraint_ntarget = 500  # length(p$yrs)
      p$areal_units_constraint_nmin = 30  # length(p$yrs)
      p$nAU_min = 100
      # to recreate the underlying data
      # xydata=substrate_db(p=p, DS="areal_units_input", redo=TRUE)

      sppoly = areal_units( p=p, redo=TRUE )  # this has already been done in aegis.polygons::01 polygons.R .. should nto have to redo
      plot( sppoly[ "AUID" ] )

      # prepare data
      sppoly = areal_units( p=p  )  # this has already been done in aegis.polygons::01 polygons.R .. should nto have to redo
   
      M = substrate_db( p=p, DS="carstm_inputs", sppoly=sppoly, redo=TRUE )  # will redo if not found
   
    }

    sppoly = areal_units( p=p  )  # this has already been done in aegis.polygons::01 polygons.R .. should nto have to redo
   
    p$space_name = sppoly$AUID 
    p$space_id = 1:nrow(sppoly)  # numst match M$space
  
  # run model and obtain predictions
    carstm_model( 
      p=p, 
      sppoly=sppoly,
      data= substrate_db( p=p, DS="carstm_inputs"), 
      # nposteriors=1000,
      # redo_fit=FALSE, # to start optim from a solution close to the final in 2021 ... 
      # debug = TRUE,
      theta = c( 1.710, 3.588, 0.008, 5.662 ) ,
      toget = c("summary", "random_spatial", "predictions"),
      posterior_simulations_to_retain = c("predictions"),
      family = "lognormal",
      control.mode = list( restart=TRUE  ) ,
      # control.mode = list( restart=FALSE, theta= c( 1.710, 3.588, 0.008, 5.662 ) ),  
      num.threads="4:2",
      verbose=TRUE 
    ) 

    # fit = carstm_model( p=p, DS="modelled_fit" )  # extract currently saved model fit
    
      # extract results
      if (0) {
    
        fit = carstm_model( p=p, DS="modelled_fit", sppoly=sppoly )  # extract currently saved model fit
    
        fit$summary$dic$dic
        fit$summary$dic$p.eff
        fit$dyear

        plot(fit)
        plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
       
      }
 

    # posterior predictive check
    M = substrate_db( p=p, DS='carstm_inputs'  )
    carstm_posterior_predictive_check(p=p, M=M  )

    # EXAMINE POSTERIORS AND PRIORS
    res = carstm_model(  p=p, DS="carstm_summary" )  # parameters in p and summary

    names(res$hypers)
    for (i in 1:length(names(res$hypers)) ){
      o = carstm_prior_posterior_compare( hypers=res$hypers, all.hypers=res$all.hypers, vn=names(res$hypers)[i] )  
      dev.new(); print(o)
    }     

    
    # oeffdir = file.path(p$data_root, "figures")  # old ... delete files ..todo
    oeffdir = file.path(p$data_root, p$carstm_model_label, "figures") #new ...
    fn_root_prefix = p$variabletomodel
    carstm_plot_marginaleffects(  p=p, outputdir=oeffdir, fn_root_prefix=fn_root_prefix ) 

   
  # bbox = c(-71.5, 41, -52.5,  50.5 )
  additional_features = features_to_add( 
      p=p, 
      isobaths=c( 100, 200, 300, 400, 500  ), 
      xlim=c(-80,-40), 
      ylim=c(38, 60) 
  )

  # maps of some of the results
  outputdir = file.path(p$modeldir, p$carstm_model_label, "maps" )

  carstm_plot_map( p=p, outputdir=outputdir, additional_features=additional_features, 
    toplot="random_spatial", probs=c(0.025, 0.975), 
    colors=rev(RColorBrewer::brewer.pal(5, "RdYlBu")) ) 

  carstm_plot_map( p=p, outputdir=outputdir, additional_features=additional_features, 
    toplot="predictions", colors=rev(RColorBrewer::brewer.pal(5, "RdYlBu"))) 
    
    
  # more direct control over map
  # random effects  ..i.e.,  deviation from lognormal model ( pure spatial effect )
    res = carstm_model(  p=p, DS="carstm_randomeffects" )  

    outfilename= file.path( outputdir, paste("substrate_grain_size_spatialeffect_carstm", "png", sep=".") )
    plt = carstm_map(  res=res, vn= c(  "space", "re_total" ), 
        sppoly=sppoly,
        transformation=log10,
        title="Substrate grainsize spatial errors (mm)",
        colors=rev(RColorBrewer::brewer.pal(5, "RdYlBu")),
        additional_features=additional_features,
        outfilename=outfilename
    )  
    plt
  


  # end



    