
# import module
import sgeLytics as sgel
from numpy import exp,log,sqrt

# define a model using Ton/Toff/EM/HLM/HLP ('direct' parameterization)
model = sgel.SgeModel ()
model.defineModel ( Ton=0.1 , Toff=2.6 , EM=20.0 , HLM=10.0 , HLP=27.0 )

# compute the CV and the mixing time Tau (half-autocorrelation time)
print 'CV = %f , Tau = %f' % ( model.giveCV() , model.giveMixingTime() )

# Alternatively, a model can be parameterized by defining constraints on CV, Tau (and HLP)
model.defineModelFromCVTauHLP ( 0.25 , 42 , 27.0 )
print 'CV = %f , Tau = %f' % ( model.giveCV() , model.giveMixingTime() )

# Although an adequate solution will always be found when it exists, it is not unique
print 'Ton = %f , Toff = %f , EM = %f , HLM = %f , HLP = %f' % ( model.Ton , model.Toff , model.EM , model.HLM , model.HLP )

# The additional constraint on HLP can be replace by a constraint on HLM
model.defineModelFromCVTauHLM ( 0.25 , 42 , 10.0 )
print 'CV = %f , Tau = %f' % ( model.giveCV() , model.giveMixingTime() )
print 'Ton = %f , Toff = %f , EM = %f , HLM = %f , HLP = %f' % ( model.Ton , model.Toff , model.EM , model.HLM , model.HLP )

# Those two functions can be used with an additional argument if both HLP and HLM should be constrained
# However in that case, the second constraint might not be respected if not compatible with CV and Tau
model.defineModelFromCVTauHLP ( 0.25 , 42 , 27.0 , 7.0 )
print 'CV = %f , Tau = %f' % ( model.giveCV() , model.giveMixingTime() )
print 'Ton = %f , Toff = %f , EM = %f , HLM = %f , HLP = %f' % ( model.Ton , model.Toff , model.EM , model.HLM , model.HLP )
