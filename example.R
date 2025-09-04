source("iHMM_Utils.R")
temp=sim_data_stud_t(seed=123,
                         TT=1000,
                         P=1,
                         Ktrue=3,
                         mu=2,
                         rho=0,
                         nu=40,
                         pers=0)

Y=temp$SimData[,1]

plot(Y,col=temp$mchain,pch=19)

hypers=NULL

hypers=list()
hypers$alpha0_a = 4;
hypers$alpha0_b = 4;
hypers$gamma_a = 3;
hypers$gamma_b = 6;
# hypers$alpha0= 1;
# hypers$gamma= 2;
hypers$sigma2 = 1.5;
hypers$mu_0 = 0.0;
hypers$sigma2_0 = 1.0;

numb=100
nums=1000
numi=10
S0=NULL
K0=4

fit=iHmmNormalSampleBeam(Y, hypers=hypers, numb=100, nums=1000, numi=5, S0=S0,K0=K0)
