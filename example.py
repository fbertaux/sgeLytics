# import module
import sgeLytics as sgel

# define a model using Ton/Toff/EM/HLM/HLP ('direct' parameterization)
model = sgel.SgeModel()
model.defineModel(Ton=0.1, Toff=2.6, EM=20.0, HLM=10.0, HLP=27.0)

# compute the CV and the mixing time Tau (half-autocorrelation time)
print('CV = %f, Tau = %f' % (model.giveCV(), model.giveMixingTime()))

print('')

# one can find a parameterization from its fluctuations properties
rg = model.rg
rm = model.rm
rp = model.rp
CV = 0.245
Tau = 42.718
model.defineModelFromCVTau_rg_rm_rp(CV, Tau, rg, rm, rp)
print('CV = %f, Tau = %f' % (model.giveCV(), model.giveMixingTime()))

print('')

# For this set of rg,rm,rp, there is no freedom in Tau
print('Tau_min = {}'.format(sgel.estimate_TauMin_from_rg_rm_rp(rg, rm, rp)))
print('Tau_max = {}'.format(sgel.estimate_TauMax_from_rg_rm_rp(rg, rm, rp)))

# when all timescales are comparable, there is more freedom in Tau
rg = rp
rm = rp
print('Tau_min = {}'.format(sgel.estimate_TauMin_from_rg_rm_rp(rg, rm, rp)))
print('Tau_max = {}'.format(sgel.estimate_TauMax_from_rg_rm_rp(rg, rm, rp)))

print('')

# Alternatively, it is possible to impose only one timescale in addition to CV, Tau
model.defineModelFromCVTauHLP(0.25, 42, 27.0)
print('CV = %f, Tau = %f' % (model.giveCV(), model.giveMixingTime()))

# But then, an adequate solution will always be found when it exists, it is not unique
print('Ton = %f, Toff = %f, EM = %f, HLM = %f, HLP = %f' %
      (model.Ton, model.Toff, model.EM, model.HLM, model.HLP))
