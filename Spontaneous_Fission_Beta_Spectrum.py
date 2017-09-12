from __future__ import division
import math
import numpy as np,glob, scipy.constants as sciconst, scipy.special as ss,  matplotlib.pyplot as plt
from pylab import *
from scipy import integrate
from scipy.integrate import quad, simps, trapz, romberg
import re

# This code takes as input 1. Cumulative fission yields, in ENSDF format, for relevant
# fissile isotopes and 2. Beta decay scheme data files in the ENSDF format for the ab 
# initio calculation of neutrino spectra. Output is the reactor neutrino spectrum using
# the Daya Bay fuel fractions.

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define functions and relevant variables. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def nonblank_lines(f):
    for l in f:
        line = l.rstrip()
        if line:
            yield line

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def dn(E,Q,Z,A,forb):

	if 0 < E <= Q:
		a=sciconst.alpha
		PI=sciconst.pi
		r0=1.25E-15
		r=r0*np.power(A,1/3)
		S=np.sqrt(1-(a*Z)**2)
		P=np.sqrt((E+511)**2-511**2)
		eta=a*Z*(E+511)/P
		F= 4*np.power(2*P*r,-(2*(1-S)))*np.power(abs(ss.gamma(S+1j*eta)/ss.gamma(2*S+1)),2)*np.exp(PI*eta)
		
	# Define Forbiddenness correction
		if forb=='1U':
			forbiddenness = P**2+(Q-E)**2
		elif forb=='2U':
			forbiddenness=(Q-E)**4+(10/3)*P**2*(Q-E)**2+P**2
		elif forb=='3U':
			forbiddenness=(Q-E)**6+7*P**2*(Q-E)**4+7*P**4*(Q-E)**2+P**6
		else: 
			forbiddenness=1

	# Branch Spectrum
		return forbiddenness*F*(np.sqrt(np.power(E,2)+2*E*511)*np.power(Q-E,2)*(E+511))
	else:
		return 0
	
def convert_sci_not(s):
	try:
		if float(s)==None:
			raise TypeError
		else:
			return float(s)

	except ValueError:	
		if '-' in s:
				return float('{:.5e}'.format(float(s.split('-')[0])*math.pow(10,-int(s.split('-')[1]))))
		elif '+' in s:
			return float('{:.5e}'.format(float(s.split('+')[0])*math.pow(10,int(s.split('+')[1]))))

z_dict={'BE': 4, 'BA': 56, 'BH': 107, 'BI': 83, 'BK': 97, 'BR': 35, 'UUH': 116, 'RU': 44, 'RE': 75, 'RF': 104, 'LU': 71, 'RA': 88, 'RB': 37, 'RN': 86, 'RH': 45, 'H': 1, 'P': 15, 'GE': 32, 'GD': 64, 'GA': 31, 'UUT': 113, 'OS': 76, 'HS': 108, 'C': 6, 'HO': 67, 'HF': 72, 'HG': 80, 'HE': 2, 'PR': 59, 'PT': 78, 'PU': 94, 'UUO': 118, 'PB': 82, 'PA': 91, 'PD': 46, 'PO': 84, 'PM': 61, 'ZN': 30, 'K': 19, 'O': 8, 'S': 16, 'W': 74, 'EU': 63, 'ZR': 40, 'ER': 68, 'MD': 101, 'MG': 12, 'MO': 42, 'MN': 25, 'MT': 109, 'U': 92, 'FR': 87, 'FE': 26, 'FM': 100, 'NI': 28, 'NO': 102, 'NA': 11, 'NB': 41, 'ND': 60, 'NE': 10, 'RG': 111, 'ES': 99, 'NP': 93, 'B': 5, 'CO': 27, 'CM': 96, 'CL': 17, 'CA': 20, 'CF': 98, 'CE': 58, 'N': 7, 'V': 23, 'CS': 55, 'CR': 24, 'CP': 112, 'CU': 29, 'SR': 38, 'UUQ': 114, 'UUP': 115, 'UUS': 117, 'KR': 36, 'SI': 14, 'SN': 50, 'SM': 62, 'SC': 21, 'SB': 51, 'SG': 106, 'SE': 34, 'YB': 70, 'DB': 105, 'DY': 66, 'DS': 110, 'LA': 57, 'F': 9, 'LI': 3, 'TL': 81, 'TM': 69, 'LR': 103, 'TH': 90, 'TI': 22, 'TE': 52, 'TB': 65, 'TC': 43, 'TA': 73, 'AC': 89, 'AG': 47, 'I': 53, 'IR': 77, 'AM': 95, 'AL': 13, 'AS': 33, 'AR': 18, 'AU': 79, 'AT': 85, 'IN': 49, 'Y': 39, 'CD': 48, 'XE': 54}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read Fission Yields ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data = []
dat = [0]*6
Yields = []           
cum_yields = open('252CF_Fission.txt','r')

# Read 235U Fission Yields ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for line in nonblank_lines(cum_yields):
	line_no = int(line[76:80])
	if line_no >=3:
		dat[0] = convert_sci_not(line[1:11])
		dat[1] = convert_sci_not(line[12:22])
		dat[2] = convert_sci_not(line[23:33])
		dat[3] = convert_sci_not(line[34:44])
		dat[4] = convert_sci_not(line[45:55])
		dat[5] = convert_sci_not(line[56:66])
		data.extend(dat)

chunks = [data[x:x+4] for x in xrange(0, len(data), 4)]
for index in range(len(chunks)):
	print chunks[index]
	if chunks[index][0] is not None:
		#print chunks[index]
		Z = int(chunks[index][0]/1000)
		A = int(chunks[index][0]-int(chunks[index][0]/1000)*1000)
		isomer = int(chunks[index][1])
		fiss_yield = chunks[index][2]
		fiss_entry = [Z,A,isomer,fiss_yield]
		Yields.append(fiss_entry)
data = []
chunks = []

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Spectrum ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

E=np.linspace(0,14000,1000)
data = np.zeros(len(E))
Qval = []
dataout= open('/Users/Gabe/Research/252Cf_spontaneous_fission_beta.txt','w')
non_decimal = re.compile(r'[^\d.]+')
number=0
filelist=glob.glob('/Users/Gabe/Research/ENSDF/*') # retrieve your listing.
for f in filelist:
	f = open(f,'r')
	L = []
	In = []
	J = []
	N = []
	F=[]
	fb=[]
	fission=[]
	isomer = 0
	tot_yield = 0
	for line in f:
		# All ENSDF formatted files end in a blank line. If a decay from an
		# isomeric state exists, a new ENSDF entry is placed after the blank
		# line. 
		if line.isspace():
			spectrum = np.zeros(len(E))
			print 'LEVELS:', L
			print 'Q-VALUES:',Qval
			print 'INTENSITY:', In 
			print 'Jpi:', J
			print 'Isomer:',isomer,"\n"
			print 'Forbiddenness:',fb
			Q_eff = [0]*len(L)
			# If there is no beta scheme information, move on
			if L == [] or In == [] or J == []:
				continue
			# Calculate branch spectra, normalize, multiply by the branching ratio, then add and multiply entire spectrum by the appropriate fission yield
			for i in range(len(L)):
				branch = np.zeros(len(E))
				Q_eff[i] = Q - L[i]
				
				for j in range(len(E)):
					branch[j]+=dn(E[j],Q_eff[i],Z,A,fb[i])
				norm = simps(branch,E)
				if norm == 0:
					continue
				else:
					branch = [x/norm for x in branch]
				branch = [y*In[i] for y in branch]
				spectrum = [n+m for (n,m) in zip(branch,spectrum)]
			
			# Multiply by fission yield
			for i in range(len(Yields)):
				if Yields[i][0] == Z and Yields[i][1] == A and Yields[i][2] == isomer:
					tot_yield+= Yields[i][3]
			FY = tot_yield
			spectrum = [t*FY for t in spectrum]
			print "ISOTOPE SPECTRUM:",spectrum
			if FY != 0:
				data = [f+g for (f,g) in zip(spectrum,data)]
			print "TOTAL SPECTRUM:", data
			
			# Create plot of dN/dE versus KE(B)
			#if Q_eff[0] > 12000 and FY!=0:
			#	plt.semilogy(E,spectrum,label=p)
			
			#Reset lists and variables
			L = []
			In = []
			J = []
			tot_yield=0
			
			continue
		if line[6] != ' ':
			continue

		# Read nuclear information on parent isotope, endpoint, tag isomers  
		if line[7] == 'P' and line[6] == ' ' and line[5] == ' ':
			p = line[0:5].strip()
			P = line[3:5].strip()
			print 'Parent: ', line[0:5].strip()
			#print 'Halflife: ', line[39:49].strip()
			if line[8:19].strip() != '0.0' and line[8:19].strip() != '0':
				isomer = 1
			if is_number(line[8:19]):
				Q = float(line[63:73]) + float(line[8:19])
				print Q
			else:
				Q = float(line[63:73])
				print Q
			A = int(non_decimal.sub('',line[0:3]))
			Z=z_dict[P]
			print Q

		# Read energy levels, beta branching ratios, spin/pairity
		if line[7] == 'L' and line[6] == ' ' and line[5] == ' ':
			if is_number(line[9:19]):
				level = float(line[9:18])
			else:
				try:
					level = float(non_decimal.sub('',line[9:18]))
				except ValueError:
					continue
			jp = line[21:38].strip()
		if line[7] == 'B' and line[6] == ' ' and line[5] == ' ' and line[21:29].strip() != '':
			if line[9:18].strip()!='':
				q = float(line[9:18])
			L.append(level)
			J.append(jp)
			forb=line[77:79].strip()
			fb.append(forb)
	
		# Create array of normalized branch intensities
			I = round(float(line[21:29])/100,len(line[21:29])) #  *Read Pg 19 of ENSDF format manual for conversion to per 100 decays
			In.append(I)
		# Normalize branching ratios
		if line[7] == 'N' and line[6] == ' ' and line[5] == ' ' and line[41:49].strip()!='' and line[31:39].strip!='':
			In = [I * float(line[41:49])*float(line[31:39]) for I in In]

	# If there is no beta scheme information, move on
	if L ==[] or In == []:
		continue

	print 'DAUGHTER LEVELS:', L
	print 'Jpi:', J
	print 'BETA ENDPOINTS:',Qval
	print 'BRANCHING RATIOS:', In 
	print len(Qval)
	print len(In)
	print len(J)

print "Length of Spectrum Array: " + str(len(data))
print "Length of Energy Array: " + str(len(E))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Print data to file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#for x,y in zip(E,data):
#	dataout.write(str(x) + ' ' + str(y) + '\n')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot Calculated Spectrum ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plt.semilogy(E,data,'k')
plt.xlabel('Electron Kinetic Energy [keV]')
plt.ylabel('dN/dE [$keV^{-1}fission^{-1}$]')
plt.title('$^{252}$Cf Spontaneous Fission Beta Spectrum')
plt.rc('font', family='serif')
plt.rc('xtick', labelsize='x-small')
plt.rc('ytick', labelsize='x-small')
plt.grid(True)
#plt.legend()
plt.xlim([0,15000])
show()