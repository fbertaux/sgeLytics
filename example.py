
# import module
import sgeLytics as sgel
from numpy import exp,log,sqrt

# define a model using Ton/Toff/EM
model = sgel.SgeModel ()
model.defineModel ( Ton=0.1 , Toff=2.6 , EM=17.0 , HLM=9.0 , HLP=27.0 , EP=1.0 )

# compute the CV and the mixing time (half-autocorrelation time)
# print 'CV = %f , Tau = %f' % ( model.giveCV() , model.giveMixingTime() )

# verify we can find back original parameters with our search functions !
# (EG,rg) = sgel.estimate_EG_RG_from_EM_rm_rp_CV_Tau (model.EM,model.rm,model.rp,model.giveCV(),model.giveMixingTime())
# print 'EG found = %f , EG expected = %f' % (EG , model.EG )
# print 'rg found = %f , rg expected = %f' % (rg , model.rg)

# test the function that parameterize a model from CV and Tau
# model = sgel.SgeModel ()
# model.defineModelFromCVTau ( 0.25 , 41.2 , 30.0 , 20.0 , 9.0 )

# print sgel.computeEG_from_rg_rm_rp_CV2PG ( 1. , 1. , 1. , 0.5 )


Tau_wanted = 5.0
CV_wanted = 0.25
HLP_wanted = 20.0

# (rg,EG,EM,HLM,best_Tau) = sgel.find_params_given_CV_Tau_HLP ( CV_wanted , Tau_wanted , HLP_wanted )
# print 'rg = %f , EG = %f , EM = %f , HLM = %f , best_Tau = %f' % (rg,EG,EM,HLM,best_Tau)
# print ' verified CV = %f , verified Tau = %f' % ( sgel.computeCV (rg,log(2.)/HLM,log(2.)/HLP_wanted,EG,EM) , sgel.estimateTau (rg,log(2.)/HLM,log(2.)/HLP_wanted,EG,EM) )

HLM_wanted = 4.0

(rg,EG,EM,HLP,best_Tau) = sgel.find_params_given_CV_Tau_HLM ( CV_wanted , Tau_wanted , HLM_wanted )
print 'rg = %f , EG = %f , EM = %f , HLP = %f , best_Tau = %f' % (rg,EG,EM,HLP,best_Tau)
print ' verified CV = %f , verified Tau = %f' % ( sgel.computeCV (rg,log(2.)/HLM_wanted,log(2.)/HLP,EG,EM) , sgel.estimateTau (rg,log(2.)/HLM_wanted,log(2.)/HLP,EG,EM) )



# (rg,EG,EM,best_Tau) = sgel.find_best_rg_EM_for_Tau_given_CV_rm_rp ( Tau_wanted , CV_wanted , log(2.)/HLM_wanted , log(2.)/HLP_wanted )
# print 'rg = %f , EG = %f , EM = %f , best_Tau = %f' % (rg,EG,EM,best_Tau)
# print ' verified CV = %f , verified Tau = %f' % ( sgel.computeCV (rg,log(2.)/HLM_wanted,log(2.)/HLP_wanted,EG,EM) , sgel.estimateTau (rg,log(2.)/HLM_wanted,log(2.)/HLP_wanted,EG,EM) )

