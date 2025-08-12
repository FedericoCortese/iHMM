source("iHMM_Utils.R")
temp=sim_data_stud_t(seed=123,
                         TT=500,
                         P=1,
                         Ktrue=3,
                         mu=1.5,
                         rho=0,
                         nu=40,
                         pers=.95)

Y=temp$SimData[,1]
hypers=NULL
numb=100
nums=1000
numi=5
S0=NULL