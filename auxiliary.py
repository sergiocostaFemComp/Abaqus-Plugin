"""

E11=           136000
E22=           9150
E33=           7700
G12=           4400
G13=           3700
G23=           3018
nu12=          0.28
nu32=          0.43
nu13=          0.35
XT=            1787
XC=            631
YT=            29
tau0=        23	#13
SL=            78
ST=            34
GIc_ft=        67
GIc_fc=        103
GIc_mt=        0.149
GIIc_mt=       0.69
Gc_mm=         0.257
gammaLcr = 0.091
eps22cr  =  0.0171
eps33cr  = 0.0381
epsNcr   = 0.0032
volDist =0.4
MU=            0.4
PO=            60
exp=           -0.6
K=             3500
w=             0.2
theta_i_12=     0
theta_i_13=     3
# end props #
density=       1.509e-06
lt =       1.0
forient=       0
nplies=        10
fric=          0.3
NLGeom=        True
inter=         'Tie plies'
subroutine=    'C:\\SIMULIA\\EstProducts\\2021\\win_b64\\code\\python2.7\\lib\\abaqus_plugins\\SingleElementTest\\subroutine\\SHF_version2.for'
mc=False
mc3=False
fc=False
mt=False
ft=False
sh12=False
sh13=False
sh23=False
MC_Exp=False
mcf='C:\\SIMULIA\\CAE\\2016\\win_b64\\code\\python2.7\\lib\\abaqus_plugins\\Experiments\\mcExp.txt'
FC_Exp=False
fcf='C:\\SIMULIA\\CAE\\2016\\win_b64\\code\\python2.7\\lib\\abaqus_plugins\\Experiments\\fcExp.txt'
MT_Exp=False
mtf='C:\\SIMULIA\\CAE\\2016\\win_b64\\code\\python2.7\\lib\\abaqus_plugins\\Experiments\\mtExp.txt'
FT_Exp=False
ftf='C:\\SIMULIA\\CAE\\2016\\win_b64\\code\\python2.7\\lib\\abaqus_plugins\\Experiments\\ftExp.txt'
SH12_Exp=True
sh12f='C:\\SIMULIA\\EstProducts\\2021\\win_b64\\code\\python2.7\\lib\\abaqus_plugins\\SingleElementTest\\MatData\\sh12Exp.txt'
SH13_Exp=False
sh13f='C:\\SIMULIA\\CAE\\2016\\win_b64\\code\\python2.7\\lib\\abaqus_plugins\\Experiments\\sh13Exp.txt'
SH23_Exp=False
sh23f='C:\\SIMULIA\\CAE\\2016\\win_b64\\code\\python2.7\\lib\\abaqus_plugins\\Experiments\\sh23Exp.txt'
HFrame1 =   'Calibrate the Model' # 'Plot the results' # 'Create CAE Model'
currentDir=os.getcwd()
odb=currentDir+'\\SingleElemTest-.odb'
absDir='C:\\SIMULIA\\EstProducts\\2021\\win_b64\\code\\python2.7\\lib\\abaqus_plugins\\SingleElementTest\\'
eps2cr = 0.09
eps3cr = 0.09

BASE_DIR = os.getcwd()
PARENT_DIR = os.path.dirname(BASE_DIR)
matProps = 'C:\\SIMULIA\\EstProducts\\2021\\win_b64\\code\\python2.7\\lib\\abaqus_plugins\\SingleElementTest\\MatData\\SHF_version2.for'

with open('C:\\SIMULIA\\EstProducts\\2021\\win_b64\\code\\python2.7\\lib\\abaqus_plugins\\SingleElementTest\\MatData\\'+'MatProps.json','r') as f:
	d=json.load(f)
# 
"""

	# data={	"E11" 		:    d["E11"],
	# 		"E22" 		:    d["E22"],
	# 		"E33" 		:    d["E33"],
	# 		"G12" 		:    d["G12"],
	# 		"G13" 		:    d["G13"],
	# 		"G23" 		:    d["G23"],
	# 		"nu12"		:    d["nu12"],
	# 		"nu32"		:    d["nu32"],
	# 		"nu13"		:    d["nu13"],
	# 		"XT"		:    d["XT"],
	# 		"YT"		:    d["YT"],
	# 		"tau0"		:    d["tau0"],
	# 		"theta_i_12":    d["theta_i_12"],
	# 		"GIc_ft"	:    d["GIc_ft"],
	# 		"theta_i_13":    d["theta_i_13"],
	# 		"GIc_mt"	:    d["GIc_mt"],
	# 		"GIIc_mt"	:    d["GIIc_mt"],
	# 		"gammaLcr"	:    d["gammaLcr"],
	# 		"eps2cr"	:    d["eps2cr"],
	# 		"eps3cr"	:    d["eps3cr"],
	# 		"epsNcr"	:    d["epsNcr"],
	# 		"MU"		:    d["MU"],
	# 		"PO"		:    d["PO"],
	# 		"exp"		:    d["exp"],
	# 		"volDist"	:    d["volDist"],
	# 		"density"	:    d["density"],
	# 		"matProps"	:    d["matProps"],
	# 		"mc"		 :	d["mc"],
	# 		"fc"         :  d["fc"],
	# 		"mt"         :  d["mt"],
	# 		"ft"         :  d["ft"],
	# 		"sh12"       :  d["sh12"],
	# 		"sh13"       :  d["sh13"],
	# 		"sh23"       :  d["sh23"],
	# 		#"MC_Exp"     :  d["MC_Exp"],
	# 		"mcf"        :  d["mcf"],
	# 		#"FC_Exp"     :  d["FC_Exp"],
	# 		"fcf"        :  d["fcf"],
	# 		#"MT_Exp"     :  d["MT_Exp"],
	# 		"mtf"        :  d["mtf"],
	# 		#"FT_Exp"     :  d["FT_Exp"],
	# 		"ftf"        :  d["ftf"],
	# 		#"SH12_Exp"   :  d["SH12_Exp"],
	# 		"sh12f"      :  d["sh12f"],
	# 		#"SH13_Exp"   :  d["SH13_Exp"],
	# 		"sh13f"      :  d["sh13f"],
	# 		#"SH23_Exp"   :  d["SH23_Exp"],
	# 		"sh23f"      :  d["sh23f"],
	# 		}