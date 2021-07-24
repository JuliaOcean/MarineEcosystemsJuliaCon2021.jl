
logdir=joinpath(exps[iexp].folder,string(exps[iexp].ID),"log")
pardir=joinpath(logdir,"tracked_parameters")

# output files parameters

fil=joinpath(pardir,"data")
nml=read(fil,MITgcm_namelist())

nml.params[1][:useSingleCpuIO]=true

nml.params[3][:pChkptFreq] = 31104000.0
nml.params[3][:chkptFreq]  = 31104000.0
nml.params[3][:dumpFreq]   = 31104000.0
nml.params[3][:taveFreq]   = 31104000.0
nml.params[3][:taveFreq]   = 31104000.0
nml.params[3][:monitorFreq]= 86400.0

write(fil,nml)
git_log_fil(exps[iexp],fil,"update parameter file : "*split(fil,"/")[end])

# output files parameters

fil=joinpath(pardir,"data.diagnostics")
nml=read(fil,MITgcm_namelist())
nml.params[1][Symbol("frequency(1)")]=2592000.0
write(fil,nml)
git_log_fil(exps[iexp],fil,"update parameter file : "*split(fil,"/")[end])

# Set up Atmosphere CO2 concentration to pre-industrial level
#  (e.g. dic_pCO2 = 0.000278 for pre-industrial level
#   or   dic_pCO2 = 0.00035 for already changed climate)

fil=joinpath(pardir,"data.dic")
nml=read(fil,MITgcm_namelist())
nml.params[3][:dic_int1]=1
nml.params[3][:dic_pCO2]=0.000278
write(fil,nml)
git_log_fil(exps[iexp],fil,"update parameter file : "*split(fil,"/")[end])

# change model run duration 
#  (e.g. nTimeSteps = 720 for one 360-day year with 1/2 day time step
#    or  nTimeSteps = 2160 for three years)

fil=joinpath(pardir,"data")
nml=read(fil,MITgcm_namelist())
nml.params[3][:nTimeSteps] = 120
write(fil,nml)
git_log_fil(exps[iexp],fil,"update parameter file : "*split(fil,"/")[end])

# rerun model with updated parameters

clean(exps[iexp])
setup(exps[iexp])
launch(exps[iexp])
