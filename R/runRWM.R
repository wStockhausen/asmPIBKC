runRWM<-function(){
  #--testing: require(TMB)
  path=".";
  #--testing: path = "src/TMB"
  TMB::compile(file.path(path,"tmbRWM.hpp")) # compile model
  dynload(TMB::dynlib("tmbRWM")) # link to R session
  
  # instantiate model object
  modRWM <- TMB::MakeADFun(data = data_list,
                         parameters = param_list,
                         DLL = "tmRWM")
}