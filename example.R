source("iHMM_Utils.R")
temp=sim_data_stud_t(seed=123,
                         TT=1000,
                         P=1,
                         Ktrue=3,
                         mu=2,
                         rho=0,
                         nu=40,
                         pers=0)
#
Y=temp$SimData[,1]

plot(Y,col=temp$mchain,pch=19)

hypers=NULL

hypers=list()
hypers$alpha0= 1;
hypers$gamma= 1;
# hypers$alpha0_a = 10;
# hypers$alpha0_b = 1;
# hypers$gamma_a = 5;
# hypers$gamma_b = 1;
hypers$sigma2 = 1.5;
hypers$mu_0 = 0.0;
hypers$sigma2_0 = 1.0;

numb=10
nums=10
numi=1000
S0=NULL
K0=4

fit=iHmmNormalSampleBeam(Y, hypers=hypers, numb=numb, nums=nums, numi=numi, S0=S0,K0=K0)
