
logdir=joinpath(exps[iexp].folder,string(exps[iexp].ID),"log")

#

mydats="data"
tmpfil=joinpath(logdir,"tracked_parameters",mydats)
nml=read(tmpfil,MITgcm_namelist())

nml.params[1][:useSingleCpuIO]=true

nml.params[3][:nTimeSteps] = 720

nml.params[3][:pChkptFreq] = 31104000.0
nml.params[3][:chkptFreq]  = 31104000.0
nml.params[3][:dumpFreq]   = 31104000.0
nml.params[3][:taveFreq]   = 31104000.0
nml.params[3][:taveFreq]   = 31104000.0
nml.params[3][:monitorFreq]= 86400.0

write(tmpfil,nml)
git_log_fil(exps[iexp],tmpfil,"update parameter file : "*mydats)

# 

mydats="data.diagnostics"
tmpfil=joinpath(logdir,"tracked_parameters",mydats)
nml=read(tmpfil,MITgcm_namelist())
nml.params[1][Symbol("frequency(1)")]=2592000.0
write(tmpfil,nml)
git_log_fil(exps[iexp],tmpfil,"update parameter file : "*mydats)

# 

mydats="data.dic"
tmpfil=joinpath(logdir,"tracked_parameters",mydats)
nml=read(tmpfil,MITgcm_namelist())
nml.params[3][:dic_int1]=1
nml.params[3][:dic_pCO2]=0.00035
write(tmpfil,nml)
git_log_fil(exps[iexp],tmpfil,"update parameter file : "*mydats)

#

clean(exps[iexp])
setup(exps[iexp])
launch(exps[iexp])
