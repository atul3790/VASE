import numpy as np
from astropy import units as u
from astropy.constants import G,sigma_sb,k_B,m_e,m_p,b_wien,g0
########### Solar Values & Constants #####################################
L_sun=u.astrophys.Lsun.cgs
Lx_sun=u.CompositeUnit(scale=10**27.2,bases=[u.erg,u.s],powers=[1,-1])
Mv_sun=u.CompositeUnit(scale=4.82,bases=[u.mag],powers=[1])
BV_sun=u.CompositeUnit(0.65,bases=[u.mag],powers=[1])
Rhkp_sun=u.CompositeUnit(10**-4.899,bases=[],powers=[])
Z_sun=u.CompositeUnit(0.013,bases=[],powers=[])
V_sun=-u.CompositeUnit(-26.74,bases=[u.mag],powers=[1])
M_sun=u.astrophys.Msun.cgs
pc=u.pc.cgs
au=u.au.cgs
###########################################################################
def MtoL(M): # Mass to Luminosity MS stars. M in Msun. Returns L in Lsun
	if M<=0.85:
		return M**(-141.7*M**4+232.4*M**3-129.1*M**2+33.29*M+0.215)
	if M>0.85 and M<=2:
		return M**4
	if M>2 and M<=55:
		return 1.4*M**3.5
	if M>55:
		return 32000*M
	return 0	
def LtoM(L): #L to M. L in Lsun. Returns M in Msun
	M_array=np.arange(0.2,250,0.001)
	Ls=np.array([MtoL(i) for i in M_array])
	loc_u=np.where(Ls>L)[0][0]
	loc_d=np.where(Ls<L)[0][-1]
	return M_array[int(np.mean([loc_u,loc_d]))]			
class Spec_Prop:
	'''
	Define the class object as Sp=star_spec(B-V,T,M,R,Rotation_period,S,R_HK,R'_HK) 
	bv is compulsory. Do provide either S or R'_HK. Rotation period helps to compute convective 
	turnover time accurately. Provide Rotation period in days. But it is not a necessity as tau is a function of B-V
	'''
    
	def __init__(self,bv,T=np.nan,M=np.nan,R=np.nan,period=np.nan, S=np.nan,rhkp=np.nan):
		self.bv=bv
		self.T=T
		self.M=M
		self.R=R
		self.Cf=self.cf()
		self.Rhk_phot=self.rhk_phot()
		if np.isnan(rhkp) and np.isnan(S):
			print('Please provide either \n1. S or \n2. R\’_HK.')
			ch=int(input('Enter choice number (1 or 2):'))
			if ch == 1:
				S=float(input('Enter S index: '))
			else:
				rhkp=float(input('Enter R\'_HK: '))
		if ~np.isnan(rhkp):
			self.Rhkp=rhkp
			self.Rhk=self.Rhkp+self.Rhk_phot
			self.S=self.Sval()
		elif ~np.isnan(S):
			self.S=S
			self.Rhk=self.rhk_eval()
			self.Rhkp=self.Rhk-self.Rhk_phot
		self.R0=self.Ro()
		if ~np.isnan(period):
			self.period=period
			self.tau = period/self.R0
		else:
			self.tau=self.conv_time()
			self.period=self.R0*self.tau
		self.t_Myr=self.Age()
		if np.isnan(T):
			self.T=self.calcT()
		if np.isnan(M):
			self.M=self.calcM()
		self.Fhk=self.FHK()
		self.RHalp=self.RHA()
		self.Rx=self.calcRx()
		self.B_surf=self.calcB()

	def cf(self):
		return 10**(0.25*self.bv**3-1.33*self.bv**2+0.43*self.bv+0.24)

	def Sval(self):
		return self.Rhk/(1.34*10**-4*self.Cf)    

	def rhk_phot(self):
		return 10**(-4.898+1.918*self.bv**2-2.893*self.bv**3)
    
	def rhk_eval(self):
		return 1.34*10**-4*self.Cf*self.S        

	def Ro(self):
		y=np.log10(10**5*self.Rhkp)
		return 10**(0.324-0.4*y-0.283*y**2-1.325*y**3)
	
	def conv_time(self):
		if self.bv<0.8 > 0:
			return 5.2+53*(self.bv-0.5) 
		else:
			return 20.5 
	def Age(self):
		y=10**5*self.Rhkp
		return 10**-6*(10.725 - 1.334*y + 0.4085*y**2 - 0.0522*y**3)
	
	def FHK(self):
		return self.S*self.Cf*self.T**4*10**-14*1.29*10**6

	def RHA(self):
		return (2.447-1.787*np.log10(self.period))*10**-5

	def calcT(self):
		return 10**(3.908 - 0.234*self.bv)

	def calcM(self):
		return 10**(0.28-0.42*self.bv)

	def calcRx(self):
		return 10**(-3.71-1.46*self.R0)

	def calcB(self):
		return 10**(2.75+0.56*self.bv-0.76*self.R0)

	def prop_disp(self):
		print('Mass : ',self.M,' M\N{CIRCLED DOT OPERATOR}')
		print('T : ',self.T,' K')
		print('B - V : ',self.bv)
		print('Period : ',self.period,' days.')
		print('S : ',self.S)
		print('Flux in CaII H \& K lines at the star surface, F_HK: ',self.Fhk)
		print('R_HK : ',self.Rhk)
		print('R\'_HK : ',self.Rhkp)
		print('Log R\'_HK : ',np.log10(self.Rhkp))
		print('Rossby number R\N{SUBSCRIPT ZERO} : ',self.R0)
		print('Convective turnover time, \N{GREEK SMALL LETTER TAU} :', self.tau,' days.')
		print('Age (Myr) : ',self.t_Myr,' Myr \n(This is calculated using an old method, not reliable for B-V>0.6)')
		print('Log Rx : ',np.log10(self.Rx))
		print('RH\N{GREEK SMALL LETTER ALPHA} : ',self.RHalp)
		print('Surface Magnetic field, B: ',self.B_surf,' G')