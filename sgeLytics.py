
#### imports
import os.path
from numpy import exp,log,sqrt
from scipy import optimize

### class to construct a given model
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
		self.rg = 1/Ton + 1/Toff
		self.EG = Ton / (Ton + Toff)
		self.rm = log(2.) / HLM 
		self.rp = log(2.) / HLP
	def defineModelFromCVTau_Fixed_rp_EM_EG (self,rp,EM,EG) :
		pass
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

### analytical expressions under an implicit form
def costEstimationTau (t,rg,rm,rp,EG,EM,autoc) : return ( computeAutocorrelationProt (rg,rm,rp,EG,EM,t) - autoc) ** 2 
def estimateTau (rg,rm,rp,EG,EM,autoc=0.5) :
	result = optimize.minimize_scalar ( costEstimationTau , args=(rg,rm,rp,EG,EM,autoc), method='golden', tol=None, options=None )
	return result.x
def costEstimationRG_from_F_rm_rp (rg,F,rm,rp) : return ( F - computeF(rg,rm,rp) )**2. / F**2.
def estimateRG_from_F_rm_rp (F,rm,rp) :
	result = optimize.minimize_scalar ( costEstimationRG_from_F_rm_rp , args=(F,rm,rp) , method='golden' , tol=None , options=None )
	return result.x
def estimateRG_from_rm_rp_CV_EG_EM (rm,rp,CV,EG,EM,error_if_impossible=True) :
	CV_2_G = CV*CV - computeCV2_M (rm=rm,rp=rp,EM=EM)
	if (CV_2_G <= 0) : 
		if error_if_impossible : raise Exception ('Impossible constraints: variance part from MRNA superior to requested total variance')
		else : return 0.
	F = CV_2_G * EG / (1-EG)
	if (F >= 1) :
		if error_if_impossible : raise Exception ('Impossible constraints: Variance part of Prot cannot be that big. Lower EG ?')
		else : return 0.
	return estimateRG_from_F_rm_rp (F,rm,rp)
def costEstimation_EG_RG_from_EM_rm_rp_CV_Tau (EG,EM,rm,rp,CV,Tau) :
	rg = estimateRG_from_rm_rp_CV_EG_EM (rm,rp,CV,EG,EM,error_if_impossible=False)
	this_tau = estimateTau (rg,rm,rp,EG,EM)
	return (this_tau-Tau)**2. / Tau**2.
def estimate_EG_RG_from_EM_rm_rp_CV_Tau (EM,rm,rp,CV,Tau,error_if_impossible=True) :
	result = optimize.minimize_scalar ( costEstimation_EG_RG_from_EM_rm_rp_CV_Tau , args=(EM,rm,rp,CV,Tau) , method='bounded' , bounds=(0.,1.) )
	EG = result.x
	if result.fun > 1e-5 :
		if error_if_impossible : raise Exception ("No solution found")
		else : return None
	rg = estimateRG_from_rm_rp_CV_EG_EM (rm,rp,CV,EG,EM,error_if_impossible=error_if_impossible)
	return (EG,rg)
def costEstimation_EG_RG_from_EM_rm_rp_CV_Taus (EG,EM,rm,rp,CV,Taus) :
	rg = estimateRG_from_rm_rp_CV_EG_EM (rm,rp,CV,EG,EM,error_if_impossible=False)
	this_taus = [ (estimateTau (rg,rm,rp,EG,EM,autoc),Tau) for (Tau,autoc) in Taus ]
	errors = [ (this_tau-Tau)**2. / Tau**2. for (this_tau,Tau) in this_taus ]
	return sum (errors)
def estimate_EG_RG_from_EM_rm_rp_CV_Taus (EM,rm,rp,CV,Taus,error_if_impossible=False) :
	result = optimize.minimize_scalar ( costEstimation_EG_RG_from_EM_rm_rp_CV_Taus , args=(EM,rm,rp,CV,Taus) , method='bounded' , bounds=(0.,1.) )
	EG = result.x
	if result.fun > 1e-5 :
		if error_if_impossible : raise Exception ("No solution found")
	rg = estimateRG_from_rm_rp_CV_EG_EM (rm,rp,CV,EG,EM,error_if_impossible=error_if_impossible)
	return (EG,rg)