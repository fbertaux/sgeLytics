
#### imports
import os.path
from numpy import exp,log,sqrt
from scipy import optimize

### class to construct/parameterize a model instance
class SgeModel (object) :
	def __init__ (self) :
		pass
	def isDefined () :
		definingAttributes = [ "Ton" , "Toff" , "EM" , "HLM" , "HLP" ]
		areDefined = [ attrib for attrib in definingAttributes if self.hasattr (attrib) ]
		return ( len (areDefined) > 0 )
	def defineModel (self,Ton,Toff,EM,HLM,HLP,EP=1.) :
		self.Ton = Ton
		self.Toff = Toff
		self.EM = EM
		self.HLM = HLM
		self.HLP = HLP
		self.EP = EP
		self.rg = 1./Ton + 1./Toff
		self.EG = Ton / (Ton + Toff)
		self.rm = log(2.) / HLM 
		self.rp = log(2.) / HLP
	def defineModelFromCVTauHLP (self,CV,Tau,HLP,desired_HLM=None) :
		(rg,EG,EM,HLM,best_Tau) = find_params_given_CV_Tau_HLP (CV,Tau,HLP,desired_HLM)
		Toff = 1. / rg / EG
		Ton = Toff * EG / (1.-EG)
		self.defineModel (Ton,Toff,EM,HLM,HLP)
	def defineModelFromCVTauHLM (self,CV,Tau,HLM,desired_HLP=None) :
		(rg,EG,EM,HLP,best_Tau) = find_params_given_CV_Tau_HLM (CV,Tau,HLM,desired_HLP)
		Toff = 1. / rg / EG
		Ton = Toff * EG / (1.-EG)
		self.defineModel (Ton,Toff,EM,HLM,HLP)		
	def giveCV (self) : return computeCV (rg=self.rg,rm=self.rm,rp=self.rp,EG=self.EG,EM=self.EM)
	def giveAutocProt (self,t) : return computeAutocorrelationProt (rg=self.rg,rm=self.rm,rp=self.rp,EG=self.EG,EM=self.EM,t=t)
	def giveMixingTime (self,autoc=0.5) : return estimateTau (rg=self.rg,rm=self.rm,rp=self.rp,EG=self.EG,EM=self.EM,autoc=autoc)


### direct, explicit analytical expressions
def computeF (rg,rm,rp) : return rp * rm * (rp+rm+rg) / (rp+rm) / (rp+rg) / (rm+rg)
def computeCV_2 (rg,rm,rp,EG,EM) :
	return computeF (rg=rg,rm=rm,rp=rp) * (1.-EG)/EG + rp/(rm+rp)/EM
def computeCV (rg,rm,rp,EG,EM) : return sqrt (computeCV_2(rg=rg,rm=rm,rp=rp,EG=EG,EM=EM))
def computeCV2_M (rm,rp,EM) : return rp / (rm+rp) / EM
def computeHG (rg,rm,rp,t) :
	if ( (rp != rm) and (rp != rg) and (rg != rm) ) : 
		return ( (rm-rg)*exp(-rp*t) + (rp-rm)*exp(-rg*t) + (rg-rp)*exp(-rm*t) ) * rm * rp  / (rm-rg) / (rp-rg) / (rp-rm)
	else :
		if rp == rm :  
			if (rg == rp) : return  rp**2 * t**2 * exp(-rp*t) / 2.
			else : return rp**2. / (rp-rg)**2 * exp(-rp*t) * ( exp((rp-rg)*t) - 1. - (rp-rg)*t )
		else :
			if ( rp == rg) : return  rm*rp / (rp-rm)**2 * exp(-rp*t) * ( exp((rp-rm)*t) - 1. - (rp-rm)*t )
			else : return  rp*rm / (rm-rp)**2 * exp(-rm*t) * ( exp((rm-rp)*t) - 1. - (rm-rp)*t )
def computeHM (rm,rp,t) :
	if ( rm != rp ) : return rp / (rp-rm) * ( exp(-rm*t) - exp(-rp*t) )
	else : return rm*t*exp(-rm*t)
def computeAutocorrelationProt (rg,rm,rp,EG,EM,t) :
	return exp(-rp*t) + computeHM (rm,rp,t) + (1.-EG) / EG / computeCV_2(rg,rm,rp,EG,EM) * rp * rm / (rp+rg) / (rm+rg) * computeHG(rg,rm,rp,t)
def computeEG_from_F_CV2PG ( F , CV2PG ) : return F / ( F + CV2PG )
def computeEG_from_rg_rm_rp_CV2PG ( rg , rm , rp , CV2PG ) : return computeEG_from_F_CV2PG ( computeF(rg,rm,rp) , CV2PG )
def computeEM_from_rm_rp_CV2PM ( rm , rp , CV2PM ) : return rp / ( rp + rm ) / CV2PM

### analytical expressions under an implicit form
def costEstimationTau (t,rg,rm,rp,EG,EM,autoc) : return ( computeAutocorrelationProt (rg,rm,rp,EG,EM,t) - autoc) ** 2 
def estimateTau (rg,rm,rp,EG,EM,autoc=0.5) :
	result = optimize.minimize_scalar ( costEstimationTau , args=(rg,rm,rp,EG,EM,autoc), method='golden', tol=None, options=None )
	return result.x
def costEstimationRG_from_F_rm_rp (rg,F,rm,rp) : return ( F - computeF(rg,rm,rp) )**2. / F**2.
def estimateRG_from_F_rm_rp (F,rm,rp) :
	result = optimize.minimize_scalar ( costEstimationRG_from_F_rm_rp , args=(F,rm,rp) , method='golden' , tol=None , options=None )
	return result.x
def cost_find_best_rg_for_Tau_given_CV2PG_rm_rp_EM ( rg , Tau , CV2PG , rm , rp , EM ) :
	# find EG that matches CV2PG for this rg (and given rm and rp)
	EG = computeEG_from_rg_rm_rp_CV2PG ( rg , rm , rp , CV2PG )
	# compute the Tau associated
	this_Tau = estimateTau (rg,rm,rp,EG,EM)
	# return the relative square error on Tau
	return (this_Tau-Tau)**2. / Tau**2.
def find_best_rg_for_Tau_given_CV2PG_rm_rp_EM ( Tau , CV2PG , rm , rp , EM ) :
	result = optimize.minimize_scalar ( cost_find_best_rg_for_Tau_given_CV2PG_rm_rp_EM , args=(Tau,CV2PG,rm,rp,EM) , method='golden' , tol=None , options=None )
	rg = result.x
	EG = computeEG_from_rg_rm_rp_CV2PG ( rg , rm , rp , CV2PG )
	best_Tau = estimateTau (rg,rm,rp,EG,EM)
	return (rg,best_Tau)
def cost_find_rg_EM_for_Tau_given_CV_rm_rp ( CV2_gene_frac , Tau , CV , rm , rp ) :
	CV2PG = CV2_gene_frac * CV * CV
	# compute the EM that gives this gene frac
	EM = rp / ( rp + rm ) / ( CV*CV - CV2PG )
	# find the best rg for Tau
	(rg,best_Tau) = find_best_rg_for_Tau_given_CV2PG_rm_rp_EM ( Tau , CV2PG , rm , rp , EM )
	return ((best_Tau-Tau)/Tau)**2.
def find_best_rg_EM_for_Tau_given_CV_rm_rp ( Tau , CV , rm , rp ) :
	result = optimize.minimize_scalar ( cost_find_rg_EM_for_Tau_given_CV_rm_rp , args=(Tau,CV,rm,rp) , method='bounded' , bounds=(0.,1.) )
	CV2_gene_frac = result.x
	CV2PG = CV2_gene_frac * CV * CV
	EM = rp / ( rp + rm ) / ( CV*CV - CV2PG )
	(rg,best_tau) = find_best_rg_for_Tau_given_CV2PG_rm_rp_EM ( Tau , CV2PG , rm , rp , EM )
	EG = computeEG_from_rg_rm_rp_CV2PG ( rg , rm , rp , CV2PG )
	return ( rg , EG , EM , best_tau )
def find_params_given_CV_Tau_HLP ( CV , Tau , HLP , desired_HLM=None ) :
	if HLP > Tau : 
		raise Exception ( 'HLP cannot be longer than Tau for this model !' )
	if not desired_HLM : 
		desired_HLM = HLP/3.
	rp = log(2.) / HLP
	rm = log(2.) / desired_HLM
	(rg,EG,EM,best_Tau) = find_best_rg_EM_for_Tau_given_CV_rm_rp ( Tau , CV , rm , rp )
	while ( abs ( ( best_Tau - Tau ) / Tau ) ) > 1e-3 :
		rm = 2. * rm
		(rg,EG,EM,best_Tau) = find_best_rg_EM_for_Tau_given_CV_rm_rp ( Tau , CV , rm , rp )
	HLM = log(2.) / rm
	return (rg,EG,EM,HLM,best_Tau)
def find_params_given_CV_Tau_HLM ( CV , Tau , HLM , desired_HLP=None ) :
	if HLM > Tau : 
		raise Exception ( 'HLM cannot be longer than Tau for this model ! (did not try to prove it but seems true)' )
	if not desired_HLP : 
		desired_HLP = HLM*3.
	rp = log(2.) / desired_HLP
	rm = log(2.) / HLM
	(rg,EG,EM,best_Tau) = find_best_rg_EM_for_Tau_given_CV_rm_rp ( Tau , CV , rm , rp )
	while ( abs ( ( best_Tau - Tau ) / Tau ) ) > 1e-3 :
		rp = 2. * rp
		(rg,EG,EM,best_Tau) = find_best_rg_EM_for_Tau_given_CV_rm_rp ( Tau , CV , rm , rp )
	HLP = log(2.) / rp
	return (rg,EG,EM,HLP,best_Tau)
def cost_estimate_TauMin_from_rg_rm_rp ( var_rep_frac , rg , rm , rp ) :
	EM = computeEM_from_rm_rp_CV2PM(rm,rp,1.-var_rep_frac)
	EG = computeEG_from_rg_rm_rp_CV2PG(rg,rm,rp,var_rep_frac)
	return estimateTau (rg,rm,rp,EG,EM)
def cost_estimate_TauMax_from_rg_rm_rp ( var_rep_frac , rg , rm , rp ) :
	return ( - cost_estimate_TauMin_from_rg_rm_rp (var_rep_frac,rg,rm,rp) )
def estimate_TauMin_from_rg_rm_rp ( rg , rm , rp ) :
	result = optimize.minimize_scalar ( cost_estimate_TauMin_from_rg_rm_rp , args=(rg,rm,rp) , method='bounded' , bounds=(0.,1.) )
	return cost_estimate_TauMin_from_rg_rm_rp(result.x,rg,rm,rp)
def estimate_TauMax_from_rg_rm_rp ( rg , rm , rp ) :
	result = optimize.minimize_scalar ( cost_estimate_TauMax_from_rg_rm_rp , args=(rg,rm,rp) , method='bounded' , bounds=(0.,1.) )
	return cost_estimate_TauMin_from_rg_rm_rp(result.x,rg,rm,rp)
def cost_estimate_alpha_from_rg_rm_rp_Tau ( var_rep_frac , rg , rm , rp , Tau ) :
	EM = computeEM_from_rm_rp_CV2PM(rm,rp,1.-var_rep_frac)
	EG = computeEG_from_rg_rm_rp_CV2PG(rg,rm,rp,var_rep_frac)
	this_Tau = estimateTau (rg,rm,rp,EG,EM)
	return ((this_Tau-Tau)/Tau)**2.
def estimate_alpha_from_rg_rm_rp_Tau ( rg , rm , rp , Tau ) :
	result = optimize.minimize_scalar ( cost_estimate_alpha_from_rg_rm_rp_Tau , args=(rg,rm,rp,Tau) , method='bounded' , bounds=(0.,1.) )
	return result.x

# not 100% trustable
def estimateRG_from_rm_rp_CV_EG_EM (rm,rp,CV,EG,EM,error_if_impossible=True) :
	CV_2_G = CV*CV - computeCV2_M (rm=rm,rp=rp,EM=EM)
	print CV_2_G / CV / CV
	if (CV_2_G <= 0) : 
		if error_if_impossible : raise Exception ('Impossible constraints: variance part from gene superior to requested total variance')
		else : return 0. # why returning 0 here ?
	F = CV_2_G * EG / (1-EG)
	if (F >= 1) :
		if error_if_impossible : raise Exception ('Impossible constraints: Variance part from gene cannot be that big. Lower EG ?')
		else : return 0. # why returning 0 here ?
	return estimateRG_from_F_rm_rp (F,rm,rp) # I think this always works, i.e. rg can be found that matches F given rm/rp
def costEstimation_EG_RG_from_EM_rm_rp_CV_Tau (EG,EM,rm,rp,CV,Tau) :
	rg = estimateRG_from_rm_rp_CV_EG_EM (rm,rp,CV,EG,EM,error_if_impossible=False)
	this_tau = estimateTau (rg,rm,rp,EG,EM)
	return ((this_tau-Tau)/Tau)**2.
def estimate_EG_RG_from_EM_rm_rp_CV_Tau (EM,rm,rp,CV,Tau,error_if_impossible=True) :
	result = optimize.minimize_scalar ( costEstimation_EG_RG_from_EM_rm_rp_CV_Tau , args=(EM,rm,rp,CV,Tau) , method='bounded' , bounds=(0.,1.) )
	EG = result.x
	if result.fun > 1e-5 :
		if error_if_impossible : raise Exception ("No solution found")
		else : return None
	rg = estimateRG_from_rm_rp_CV_EG_EM (rm,rp,CV,EG,EM,error_if_impossible=error_if_impossible)
	return (EG,rg)

