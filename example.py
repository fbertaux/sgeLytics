
# import module
import sgeLytics as sgel

# define a model using Ton/Toff/EM
model = sgel.SgeModel ()
model.defineModel ( Ton=0.1 , Toff=2.6 , EM=17.0 , HLM=9.0 , HLP=27.0 , EP=1.0 )

# compute the CV and the mixing time (half-autocorrelation time)
print 'CV = %f , Tau = %f' % ( model.giveCV() , model.giveMixingTime() )

# verify we can find back original parameters with our search functions !
(EG,rg) = sgel.estimate_EG_RG_from_EM_rm_rp_CV_Tau (model.EM,model.rm,model.rp,model.giveCV(),model.giveMixingTime())
print 'EG found = %f , EG expected = %f' % (EG , model.EG )
print 'rg found = %f , rg expected = %f' % (rg , model.rg)