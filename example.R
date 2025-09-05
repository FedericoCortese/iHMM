source("iHMM_Utils.R")
temp=sim_data_stud_t(seed=123,
                         TT=1000,
                         P=1,
                         Ktrue=3,
                         mu=2,
                         rho=0,
                         nu=40,
                         pers=0.99)
#
Y=temp$SimData[,1]

plot(Y,col=temp$mchain,pch=19)

hypers=NULL

hypers=list()
# hypers$alpha0= 1;
# hypers$gamma= 1;
hypers$alpha0_a = 6;
hypers$alpha0_b = 1;
hypers$gamma_a = 4;
hypers$gamma_b = 1;
hypers$sigma2 = 1.5;
hypers$mu_0 = 0.0;
hypers$sigma2_0 = 1.0;

numb=10
nums=500
numi=5
S0=NULL
K0=4
niter <- numb + (nums - 1L) * numi;niter

fit=iHmmNormalSampleBeam(Y, hypers=hypers, numb=numb, nums=nums, numi=numi, S0=S0,K0=K0)
plot(fit$stats$K,type='l')
table(temp$mchain,fit$S[[nums]]$S)



# vonmises ----------------------------------------------------------------

source("iHMM_Utils.R")
TT=1000
mu_vec=c(0,pi)
kappa_vec=c(2,2)

y=sim_data_vonmises(seed = 123,
                                TT=TT,
                                Ktrue = 2,
                                mu_vec=mu_vec,
                                kappa_vec=kappa_vec,
                                pers = 0.95)

plot(y$Y,col=y$mchain,pch=19)

fitvm=iHmmvonMisesSampleBeam(y$Y, hypers=NULL, 
                                   numb=100, 
                                   nums=1000, 
                                   numi=5, 
                                   S0=NULL,
                                   K0=NULL)


