#--script to test tmbRWM implementation
#----will create and test tmbRWM.cpp, 
#----so this is NOT a test of the package code
  require(TMB);
  if (TRUE){
    fn = "inst/data/tmbPIBKC.inp";
    compile=TRUE;
    verbose=1;
  }
  
  #--read input control/data file
  readVal<-function(str){return(as.numeric(strsplit(str,"#",fixed=TRUE)[[1]][1]));}
  res      = readLines(con=fn);
  #----control values
  styr     = readVal(res[1]);
  endyr    = readVal(res[2]);
  nobs     = readVal(res[3]);
  uncType  = readVal(res[4]);#--uncertainty type (0=cv, 1=sd)
  #----data (year, observed, uncertainty)
  mat      = matrix(data=NA,nrow=nobs,ncol=3);
  for (i in 1:nobs) mat[i,] = as.numeric(strsplit(res[i+5]," ",fixed=TRUE)[[1]][1:3]);
  colnames(mat) = c("year","obs","cv");
  #----create tibble for later use
  tbl_dat = tibble::as_tibble(mat,.name_repair="minimal") %>% 
              dplyr::mutate(cv=(uncType==0)*cv+(uncType==1)*cv/obs,
                            uci=exp(log(obs)+sqrt(log(1+cv*cv))),
                            lci=exp(log(obs)-sqrt(log(1+cv*cv))))
  srv_yrs=mat[,1];
  srv_obs=mat[,2];
  srv_cvs=mat[,3];
  
  nyrs = endyr-styr+1;
  mod_yrs = styr:endyr;
  o2m_yrs = -1+0*srv_yrs;#--mapping from srv_yrs to mod_yrs
  for (i in 1:nobs){
    if ((styr<=srv_yrs[i])&&(srv_yrs<=endyr)) o2m_yrs[i] = srv_yrs[i]-styr;
  }
  
  datlst = list(verbose=0,
              styr=styr,
              endyr=endyr,
              uncType=uncType,
              srv_yrs=srv_yrs,
              srv_obs=srv_obs,
              srv_cvs=srv_cvs,
              o2m_yrs=o2m_yrs);
  parlst = list(logSD=0,
                logPred=rep(1,nyrs)
                );
  
  DLL_PATH = "src/TMB";
  DLL_NAME = "tmbRWM";
  DLL = file.path(DLL_PATH,DLL_NAME);
  if (compile){
    objFile<-paste0(DLL,".o");
    if (file.exists(objFile)) file.remove(objFile);
    soFile<-paste0(DLL,".so");
    if (file.exists(soFile)) file.remove(soFile);
    if (TMB::compile(paste0(DLL,".cpp"))!=0){
      msg<-"Could not compile code.\n";
      stop(msg);
    }
  }
  dyn.load(TMB::dynlib(DLL));
  DLL = DLL_NAME;
  
  #create the objective function
  objFun<-TMB::MakeADFun(data=datlst,
                         parameters=parlst,
                         random=c("logPred"),
                         DLL=DLL,
                         checkParameterOrder=TRUE,
                         silent=FALSE);

  #optimize the parameters
  opt = TMBhelper::fit_tmb(objFun);
  
  #--check estimability by calculating the matrix of second-derivatives of the marginal likelihood 
  #--w.r.t. fixed effects, to see if any linear combinations are not estimable 
  #--(i.e. cannot be uniquely estimated conditional upon model structure and available data, e.g., 
  #--resulting in a likelihood ridge and singular, non-invertable Hessian matrix)
  chk = TMBhelper::check_estimability(objFun);
  
  #--get sdreport
  sdrep = TMB::sdreport(objFun,bias.correct=TRUE,
                        bias.correct.control=list(sd=TRUE,split=NULL,nsplit=NULL));
  sum_sdrep = tibble::as_tibble(summary(sdrep),rownames="name");
  
  #--extract and plot results
  require(magrittr); require(ggplot2);
  tbl = dplyr::bind_cols(sum_sdrep %>% dplyr::filter(name=="pred"),year=styr:endyr);
  tbl %<>% dplyr::mutate(uci=Estimate+`Std. Error`,
                         lci=Estimate-`Std. Error`,
                         uci_bc=`Est. (bias.correct)`+`Std. (bias.correct)`,
                         lci_bc=`Est. (bias.correct)`-`Std. (bias.correct)`);
  p1 =   ggplot(tbl,aes(x=year,y=`Est. (bias.correct)`,ymin=lci_bc,ymax=uci_bc)) + 
          geom_ribbon(alpha=0.8,colour="green",fill="green") + 
          geom_line(colour="green") + 
          geom_line(aes(x=year,y=`Estimate`),colour="dark blue",linetype=2) +
          geom_ribbon(aes(x=year,y=`Estimate`,ymin=lci,ymax=uci),fill="dark blue",alpha=0.3) +
          geom_point(data=tbl_dat,mapping=aes(x=year,y=obs),inherit.aes=FALSE) +
          geom_errorbar(data=tbl_dat,mapping=aes(x=year,ymin=lci,ymax=uci),inherit.aes=FALSE) +
          wtsPlots::getStdTheme();
  p2 =  ggplot(tbl,aes(x=year,y=`Est. (bias.correct)`,ymin=lci_bc,ymax=uci_bc)) + 
          geom_ribbon(alpha=0.8,colour="green",fill="green") + 
          geom_line(colour="green") + 
          geom_line(aes(x=year,y=`Estimate`),colour="dark blue",linetype=2) +
          geom_ribbon(aes(x=year,y=`Estimate`,ymin=lci,ymax=uci),fill="dark blue",alpha=0.3) +
          geom_point(data=tbl_dat,mapping=aes(x=year,y=obs),inherit.aes=FALSE) +
          geom_errorbar(data=tbl_dat,mapping=aes(x=year,ymin=lci,ymax=uci),inherit.aes=FALSE) +
          scale_y_log10() + 
          wtsPlots::getStdTheme();

  # --get ADMB results and add to plots
    res.RE<-PBSmodelling::readList('rwout.rep')
    #finish off the output
    res<-rPIBKC::calcCIs(res.RE$est,res.RE$cv,pdfType="lognormal",ci=0.66,verbose=FALSE);
    tbl_admb<-tibble::tibble(year=res.RE$yrs,type='RE',value=res.RE$est,lci=res$lci,uci=res$uci);
    p1 = p1 + geom_line(data=tbl_admb,mapping=aes(x=year,y=value),colour="red",inherit.aes=FALSE) + 
              geom_ribbon(data=tbl_admb,mapping=aes(x=year,ymin=lci,ymax=uci),colour="red",alpha=0.3,inherit.aes=FALSE)
    p2 = p2 + geom_line(data=tbl_admb,mapping=aes(x=year,y=value),colour="red",inherit.aes=FALSE) + 
              geom_ribbon(data=tbl_admb,mapping=aes(x=year,ymin=lci,ymax=uci),colour="red",alpha=0.3,inherit.aes=FALSE)
    print(p1); print(p2);
  
  #--compare ln-scale process error uncertainty
  tblSD = dplyr::bind_rows(
            dplyr::bind_cols(type="TMB",sum_sdrep %>% dplyr::filter(name %in% c("logSD","sd")) %>% dplyr::select(name,est=Estimate,se=`Std. Error`)),
            tibble::tibble_row(type="ADMB",name="logSD",est=res.RE$sdrepLogSdLam,se=res.RE$sdrepLogSdLam.sd)
          );
    ggplot(tblSD,aes(x=name,y=est,ymin=est-se,ymax=est+se,colour=type,fill=type)) +
      geom_col(position="dodge",alpha=0.5) +
      geom_linerange(position=position_dodge(width=1)) +
      labs(x="parameter",y="estimate") +
      wtsPlots::getStdTheme();
    
  if (is.loaded(DLL)) dyn.unload(DLL);
  