import math
import numpy as np, scipy.constants as sciconst, scipy.special as ss,  matplotlib.pyplot as plt
from scipy.integrate import quad, simps, trapz, romberg
import re

def is_number(s):
	"""This function determines if a value is a number."""
	try:
		float(s)
		return True
	except ValueError:
		return False

def dn(E,Q,Z,A,forb):
	"""This function, dn(E,Q,Z,forb) Calculates the branch beta kinetic energy spectrum. It 
# takes the following variables: Q, the branch endpoint; Z, the number of protons in the 
# parent nucleus; forb, the forbiddenness. E the kinetic energy of the electron."""
	if 0 < E <= Q:
		a=sciconst.alpha
		PI=sciconst.pi
		r0=1.25E-15
		r=r0*np.power(A,1/3)
		S=np.sqrt(1-(a*Z)**2)
		P=np.sqrt((E+.511)**2-.511**2)
		eta=a*Z*(E+.511)/P
		F= 4*np.power(2*P*r,-(2*(1-S)))*np.power(abs(ss.gamma(S+1j*eta)/ss.gamma(2*S+1)),2)*np.exp(PI*eta)
		
	# Define Forbiddenness correction
		if forb=='1U':
			forbiddenness = P**2+(Q-E)**2
		elif forb=='2U':
			forbiddennes=(Q-E)**4+(10/3)*P**2*(Q-E)**2+P**2
		elif forb=='3U':
			forbiddenness=(Q-E)**6+7*P**2*(Q-E)**4+7*P**4*(Q-E)**2+P**6
		else: 
			forbiddenness=1

	# Branch Spectrum
		return forbiddenness*F*(np.sqrt(np.power(E,2)+2*E*.511)*np.power(Q-E,2)*(E+.511))
	else:
		return 0

def calc_spectrum(L,Q,In,Z,A,fb):
	""" This function calculates the total spectrum for the isotope. """
	global E
	normal=[]
	E=np.linspace(0,14,5000)
	spectrum = np.zeros(len(range(5000)))
	for i in range(len(L)):
		branch = np.zeros(len(range(5000)))
		Q_eff = Q - L[i]
		for j in range(len(E)):
			branch[j]+=dn(E[j],Q_eff,Z,A,fb[i])
		norm = simps(branch,E)
		normal.append(norm)
		branch = [x/norm for x in branch]
		branch = [y*In[i] for y in branch]
		spectrum = [n+m for (n,m) in zip(branch,spectrum)]
	print 'Normalizations: ',normal
	norm = simps(spectrum,E)
	print "Norm: ", norm
	spectrum = [x/norm for x in spectrum]
	#Optional: plot branch spectrum
	#create plot of dN/dE versus KE(B)
	#plt.errorbar(E,spectrum, yerr=None,ecolor='y',color='k', label = "Calculated")
	#plt.semilogy(E,spectrum)
	return spectrum

def read_ENSDF(infile):
	"""This function reads an ENSDF file, puts the relevant beta decay variables into arrays, and will plot the beta spectrum at the end of the file, with optional data output."""
	L = []
	Qval = []
	In = []
	J = []
	fb =[]
	isomer = 0
	z_dict={'BE': 4, 'BA': 56, 'BH': 107, 'BI': 83, 'BK': 97, 'BR': 35, 'UUH': 116, 'RU': 44, 'RE': 75, 'RF': 104, 'LU': 71, 'RA': 88, 'RB': 37, 'RN': 86, 'RH': 45, 'H': 1, 'P': 15, 'GE': 32, 'GD': 64, 'GA': 31, 'UUT': 113, 'OS': 76, 'HS': 108, 'C': 6, 'HO': 67, 'HF': 72, 'HG': 80, 'HE': 2, 'PR': 59, 'PT': 78, 'PU': 94, 'UUO': 118, 'PB': 82, 'PA': 91, 'PD': 46, 'PO': 84, 'PM': 61, 'ZN': 30, 'K': 19, 'O': 8, 'S': 16, 'W': 74, 'EU': 63, 'ZR': 40, 'ER': 68, 'MD': 101, 'MG': 12, 'MO': 42, 'MN': 25, 'MT': 109, 'U': 92, 'FR': 87, 'FE': 26, 'FM': 100, 'NI': 28, 'NO': 102, 'NA': 11, 'NB': 41, 'ND': 60, 'NE': 10, 'RG': 111, 'ES': 99, 'NP': 93, 'B': 5, 'CO': 27, 'CM': 96, 'CL': 17, 'CA': 20, 'CF': 98, 'CE': 58, 'N': 7, 'V': 23, 'CS': 55, 'CR': 24, 'CP': 112, 'CU': 29, 'SR': 38, 'UUQ': 114, 'UUP': 115, 'UUS': 117, 'KR': 36, 'SI': 14, 'SN': 50, 'SM': 62, 'SC': 21, 'SB': 51, 'SG': 106, 'SE': 34, 'YB': 70, 'DB': 105, 'DY': 66, 'DS': 110, 'LA': 57, 'F': 9, 'LI': 3, 'TL': 81, 'TM': 69, 'LR': 103, 'TH': 90, 'TI': 22, 'TE': 52, 'TB': 65, 'TC': 43, 'TA': 73, 'AC': 89, 'AG': 47, 'I': 53, 'IR': 77, 'AM': 95, 'AL': 13, 'AS': 33, 'AR': 18, 'AU': 79, 'AT': 85, 'IN': 49, 'Y': 39, 'CD': 48, 'XE': 54}
	non_decimal = re.compile(r'[^\d.]+')
	for line in infile:
		print line

		# Empty line indicates end of file or missing dataset
		if line[0:80].isspace():
			if L == [] or Qval == [] or In == [] or J == []:
				print "No beta scheme information available for this ENSDF file. \n"
				continue
			print "g.s. -> g.s. Q value:",Q
			print 'LEVELS:', L
			print 'Q-VALUES:',Qval
			print 'INTENSITY:', In 
			print 'Jpi:', J
			print 'Forbiddenness: ',fb
			
			spectrum=calc_spectrum(L,Q,In,Z,A,fb)
			
			# Data output
			user_prompt = raw_input('Print data to file (betaspectrum.txt)? ')
			if user_prompt.lower() == 'yes' or user_prompt.lower() == 'y':
				dataout= open('Output_Data/' + P + '_beta_spectrum.txt','w+')
				dataout.write( 'Energy(keV) dN/dE(keV^-1)\n')
				for x,y in zip(E,spectrum):
					dataout.write(str(x) + ' ' + str(y) + '\n')
				dataout.close()

			plt.plot(E,spectrum)
			plt.xlabel('Electron Kinetic Energy (MeV)')
			plt.ylabel('dN/dE ($MeV^{-1}$)')
			plt.title(str(A) + str(P) + ' Beta Spectrum')
			plt.grid(True)
			plt.xlim([0,14])
			plt.show()
			
			# clear parameter arrays
			L = []
			Qval = []
			In = []
			J = []

		#Get information on parent isotope and calculate Q value 
		if line[7] == 'P' and line[6] == ' ' and line[5] == ' ':
			p = line[0:5].strip()
			P = line[3:5].strip()
			print 'Parent: ', line[0:5].strip()
			#print 'Halflife: ', line[39:49].strip()
			if line[8:19].strip() != '0.0' and line[8:19].strip() != '0':
					isomer = 1
			Q = (float(line[64:73]) + float(non_decimal.sub('',line[8:19])))/1000
			A = int(non_decimal.sub('',line[0:3]))
			Z=z_dict[P]+1
			print Q
			if line[8:19].strip() != '0.0' and line[8:19].strip() != '0':
				isomer = 1
			#print Q
		# Calculate Q value for each branch
		if line[7] == 'L' and line[6] == ' ' and line[5] == ' ':
			if is_number(line[9:18]):
				level = float(line[9:18])/1000
			else: 
				continue
			jp = line[21:38].strip()
			q = round((Q - level),4)
		if line[7] == 'B' and line[6] == ' ' and line[5] == ' ' and line[21:29].strip() != '':
			if line[9:18].strip()!='':
				q = float(line[9:18])/1000
			L.append(level)
			Qval.append(q)
			J.append(jp)
			fb.append(line[77:79].strip())
			
		# Create array of normalized branch intensities
			I = round(float(line[21:29])/100,len(line[21:29])) #  *Read Pg 19 of ENSDF format manual for conversion to per 100 decays
			In.append(I)

		if line[7] == 'N' and line[6] == ' ' and line[5] == ' ' and line[41:49].strip()!='' and line[31:39].strip!='':
			In = [I * float(line[41:49])*float(line[31:39]) for I in In]

	infile.close()

ensdf_file = open('/Users/Gabe/Research/BetaData/9BE','r')
read_ENSDF(ensdf_file)


