/// @file tmbRWM.hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

// name of function below **MUST** match filename
template <class Type>
Type tmbRWM(objective_function<Type>* obj) {
  //data inputs
  DATA_INTEGER(verbose);            //flag to print debugging info
  if (verbose) std::cout<<"1"<<std::endl;
  DATA_INTEGER(styr);               //start year for interpolation
  if (verbose) std::cout<<"2"<<std::endl;
  DATA_INTEGER(endyr);              //end year for interpolation
  if (verbose) std::cout<<"3"<<std::endl;
  DATA_INTEGER(uncType);            //uncertainty type (0=cv's, 1=arithmetic cv's)
  if (verbose) std::cout<<"5"<<std::endl;
  DATA_IVECTOR(srv_yrs);            //years with survey observations
  if (verbose) std::cout<<"6"<<std::endl;
  DATA_VECTOR(srv_obs);             //survey observations
  if (verbose) std::cout<<"7"<<std::endl;
  DATA_VECTOR(srv_cvs);             //survey observation cv's or sd's
  if (verbose) std::cout<<"8"<<std::endl;
  DATA_IVECTOR(o2m_yrs);            //map from survey year index to model year index
  if (verbose) std::cout<<"9"<<std::endl;
  //end of input data
  
  //parameters
  PARAMETER(logSD);             //ln-scale process standard deviation
  if (verbose) std::cout<<"9"<<std::endl;
  PARAMETER_VECTOR(logPred);   //ln-scale predictions
  if (verbose) std::cout<<"10"<<std::endl;
    
  //constants
  Type pi = Type(std::acos(-1.0));

  int nobs = srv_yrs.size();
  vector<Type> srv_sds(nobs);
  if (verbose) std::cout<<"1"<<std::endl;

  //convert uncertainty info
  //uncType==0 => cv's are given, so no need to do anything
  if (uncType==1){ //unc = arithmetic std devs's
      for (int i=0;i<nobs;i++) srv_cvs[i] = srv_cvs[i]/srv_obs[i];
  }
  for (int i=0;i<nobs;i++) srv_sds[i] = sqrt(log(1.0+(srv_cvs[i]*srv_cvs[i])));
  if (verbose) std::cout<<"2"<<std::endl;

  //compute process-related quantities
  int nyrs  = endyr-styr+1;  //number of process years
  if (verbose) std::cout<<"nyrs = "<<nyrs<<std::endl;
  int nyrsc = logPred.size();
  if (nyrs!=nyrsc){Rf_error("nyrs!=nyrsc");}
  Type sd = exp(logSD);     //process error sd
  if (verbose) std::cout<<"sd = "<<sd<<std::endl;
  vector<Type> zprc(nyrs-1);//z-scores for process error
  for (int i=0;i<(nyrs-1);i++) zprc[i] = (logPred[i+1]-logPred[i])/sd;
  if (verbose) std::cout<<"3"<<std::endl;
  
  //compute objective function
  Type nll = Type(0.0);
  //add in process error
  for (int i=0;i<(nyrs-1);i++){
    nll += 0.5*(log(2.0*pi*sd*sd)+zprc[i]*zprc[i]);
  }
  if (verbose) std::cout<<"1"<<std::endl;
  //add in observation error
  vector<Type> zscrs(nobs);
  zscrs = Type(0.0);
  for (int i=0;i<nobs;i++){
    int j = o2m_yrs[i];
    if (!(j<0)) {
      zscrs[i]= (logPred[j]-log(srv_obs[i]))/srv_sds[i];
      nll += 0.5*(log(2.0*pi*srv_sds[i]) + zscrs[i]*zscrs[i]);
    }
  }
  if (verbose) std::cout<<"2"<<std::endl;

  //--make report
  REPORT(nll);
  if (verbose) std::cout<<"0"<<std::endl;
  REPORT(nyrs);
  if (verbose) std::cout<<"1"<<std::endl;
  REPORT(zscrs);
  if (verbose) std::cout<<"4"<<std::endl;
  ADREPORT(sd);
  if (verbose) std::cout<<"2"<<std::endl;
  vector<Type> pred = exp(logPred);
  ADREPORT(pred);
  if (verbose) std::cout<<"3"<<std::endl;
  
  return nll;  
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this


