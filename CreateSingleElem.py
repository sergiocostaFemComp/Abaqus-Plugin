# Code created by Sergio Costa - 2015-2024
#
from abaqus import *
from abaqusConstants import *
import __main__
import regionToolset
import xyPlot
#
import json
import numpy as np
import os, shutil
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

import numpy as np
import json


# Script directory #Only works when run from abaqus
absPath = os.path.abspath(__file__)  
absDir = os.path.dirname(absPath)
subDir = str(absDir)+ '\\subroutine\\'

# Abaqus working directory
currentDir=os.getcwd()



# if 2>1:
def CreateElms(E11,E22,E33,G12,G13,G23,nu12,nu32,nu13,		# Material Props
	XT,YT,thetai_12,thetai_13,GIc_ft,		# Damage props	
	gammaLcr,eps2cr,eps3cr,tau0,MU,PO,exp,volDist,density,	# Crash props
	XC,YC,SL, XC_eps,YC_eps,SL_eps,XT_eps,YT_eps,matProps,										# Element/Subroutine
 	sh12f,HFrame1):	         										
		
	lt = 1
	w  = 0.2
	
	theta_i = max(thetai_12, thetai_13)

	# Data to updated json file
	data={	"E11" 		:E11,
			"E22" 		:E22,
			"E33" 		:E33,
			"G12" 		:G12,
			"G13" 		:G13,
			"G23" 		:G23,
			"nu12"		:nu12,
			"nu32"		:nu32,
			"nu13"		:nu13,
	     ##########################
			"XT"		:XT,
			"YT"		:YT,
			"thetai_12":thetai_12,
			"thetai_13":thetai_13,
			"GIc_ft"	:GIc_ft,
			#"GIc_mt"	:GIc_mt,
			#"GIIc_mt"	:GIIc_mt,
		##########################
			"gammaLcr"	:gammaLcr,
			"eps2cr"	:eps2cr,
			"eps3cr"	:eps3cr,
			#"epsNcr"	:epsNcr,
			"tau0"		:tau0,
			"MU"		:MU,
			"PO"		:PO,
			"exp"		:exp,
			"volDist"	:volDist,
			"density"	:density,
			"matProps"	:matProps,
			"XC"		:XC,
			"YC"        :YC,
			"SL"        :SL,
			"XC_eps"	:XC_eps,
			"YC_eps"	:YC_eps,
			"SL_eps"	:SL_eps,
			"XT_eps"	:XT_eps,
			"YT_eps"	:YT_eps,

			"sh12f"     :sh12f,
			"fc10"		: - np.sign(theta_i)*1,   # 1 deg oposite direction
			"fc20"		:np.sign(theta_i)*1,
			"fc30"		:np.sign(theta_i)*2,
			"fc45"		:np.sign(theta_i)*3.5,
			}

	# Save all the parameters used to create the model
	with open(matProps,'w') as f:
			json.dump(data,f)		
	#
	if HFrame1=='Create single elements':
		print('CREATING THE CAE MODEL...')
		mlist=('mc','mc3','fc','mt','ft','sh12','sh13','sh23','fc10', 'fc20', 'fc30', 'fc45' ) 

		#
		mname='SingleElems'
		# 
		# Create profile
		mdb.Model(name=mname, modelType=STANDARD_EXPLICIT)
		mm=mdb.models[mname]
		#
		for i in range(0, len(mlist)):
		#
		# Define specific properties according to the name
		#	- Position  #Starts drawing the cube from top left to bot right - grows in z
		#   - Fibre orientation
		#   - Displacement value
		#   - Equation for shear cases
			ori = 0.0;	
			axisRot=AXIS_1 
			disp = lt/10.0
			shdisp=UNSET
			shlock=UNSET
			if mlist[i]== 'mc': 
				p1x = 0.0
				p1y = lt
				p2x = lt
				p2y = 0.0; 
				disp =  -lt/10.0
				
			if mlist[i] == 'fc': 
				axisRot=AXIS_3; 
				ori = 90.0
				p1x = 2*lt; p1y = lt;     p2x = 3*lt; p2y = 0.0; 
				disp =  -lt/10.0
				
			if mlist[i] == 'mt': 
				p1x = 4*lt; p1y = lt;     p2x = 5*lt; p2y = 0.0; 
				
			if mlist[i] == 'ft': 
				axisRot=AXIS_3;	
				ori  = 90.0
				p1x = 6*lt; p1y = lt;     p2x = 7*lt; p2y = 0.0; 
				
			# "1st row"
			if mlist[i] == 'fc10': 
				axisRot=AXIS_3;	
				ori  = 90.0 + data['fc10']  # there is a mat prop + this orientation
				p1x = 0.0
				p1y = 3*lt
				p2x = lt
				p2y = 2*lt; 
				disp =  -lt/10.0
				
			if mlist[i] == 'fc20': 
				axisRot=AXIS_3;	
				ori  = 90.0 + data['fc20'] 
				p1x = 2*lt
				p1y = 3*lt
				p2x = 3*lt
				p2y = 2*lt;  
				disp =  -lt/10.0
				
			if mlist[i] == 'fc30': 
				axisRot=AXIS_3;	
				ori  = 90.0 + data['fc30'] 
				p1x = 4*lt
				p1y = 3*lt
				p2x = 5*lt
				p2y = 2*lt; 
				disp =  -lt/10.0
				
			if mlist[i] == 'fc45': 
				axisRot=AXIS_3;	
				ori  = 90.0 + data['fc45'] 
				p1x = 6*lt
				p1y = 3*lt
				p2x = 7*lt
				p2y = 2*lt; 
				disp =  -lt/10.0
				
		# 3rd row
			if mlist[i]== 'mc3': 
				p1x = 0.0
				p1y = -2*lt
				p2x = lt
				p2y = -lt; 
				ori  = 90.0
				disp =  lt/10.0

			if mlist[i] == 'sh12':
				p1x = 2*lt
				p1y = -2*lt
				p2x = 3*lt
				p2y = -lt; 
				ori  = 0.0
				disp =  UNSET;	shdisp= lt/5.0; 		shlock=0.0
				
			if mlist[i] == 'sh13':
				p1x = 4*lt
				p1y = -2*lt
				p2x = 5*lt
				p2y = -lt
				ori  = 90.0
				disp =  UNSET;	shdisp= lt/5.0;		shlock=0.0
				
			if mlist[i] == 'sh23':
				p1x = 6*lt
				p1y = -2*lt
				p2x = 7*lt
				p2y = -lt
				axisRot=AXIS_2;	ori  = 90.0
				disp =  UNSET; 	shdisp= lt/5.0;		shlock=0.0
			
			s = mm.ConstrainedSketch(name='__profile__', sheetSize=200.0)
			g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
			s.setPrimaryObject(option=STANDALONE)
			s.rectangle(point1=(p1x, p1y), point2=(p2x, p2y))
			p = mm.Part(name=mlist[i], dimensionality=THREE_D, type=DEFORMABLE_BODY)
			p = mm.parts[mlist[i]]
			p.BaseSolidExtrude(sketch=s, depth=lt)
			s.unsetPrimaryObject()

			sv1=session.viewports['Viewport: 1']
			sv1.setValues(displayedObject=p)
			sv1.partDisplay.setValues(sectionAssignments=ON, engineeringFeatures=ON)
			sv1.partDisplay.geometryOptions.setValues(referenceRepresentation=OFF)
		#	
		# Material 
			mm.Material(name='NCF')
			mmat=mm.materials['NCF']
			mmat.Density(table=((density, ), ))
			mmat.UserMaterial(mechanicalConstants=(
				E11,E22,E33,G12,G13,G23,nu12,nu32,nu13,
				XT,YT,thetai_12,thetai_13,GIc_ft,
				gammaLcr,eps2cr,eps3cr,tau0,MU,PO,exp,volDist))
			mmat.Depvar(deleteVar=23, n=85)
		#
		# Assembly
			a = mm.rootAssembly
			sv1.setValues(displayedObject=a)
			sv1.assemblyDisplay.setValues(optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
			a.DatumCsysByDefault(CARTESIAN)

			a.Instance(name=mlist[i]+'-1', part=p, dependent=ON)
			ai=a.instances[mlist[i]+'-1']
			sv1.assemblyDisplay.setValues(adaptiveMeshConstraints=ON)
		#
		# Create step amplitude and time increment only once - Otherwise everything associated with the step will disappear
			if i==0: 
				mm.ExplicitDynamicsStep(name='ApplyLoad', previous='Initial', nlgeom=ON,maxIncrement=1e-07)
				mm.steps['ApplyLoad'].setValues(timePeriod=0.05)
				sv1.assemblyDisplay.setValues(loads=ON, bcs=ON, predefinedFields=ON, connectors=ON, adaptiveMeshConstraints=OFF)
				sv1.view.setValues(session.views['Front'])
			
			mm.rootAssembly.regenerate()
			f1 = a.instances[mlist[i]+'-1'].faces
		# 
		# Define geometry face sets for BCs
			BotSet=a.Set(faces=f1.findAt((( p2x-0.1,p2y,lt/2.0),),), name=mlist[i]+'BotSet')
			TopSet=a.Set(faces=f1.findAt((( p1x+0.1,p1y,lt/2.0),),), name=mlist[i]+'TopSet')
			EdgeSet =a.Set(edges=ai.edges.findAt((( p1x+0.1,p2y,lt),),), name=(mlist[i]+'EdgeSet'))
		#
		# Definition of Loading/Unloading and Pressure amplitudes # Creates as many amplitudes as asked
			mm.SmoothStepAmplitude(name='SmoothAmp', timeSpan=STEP, data=((0.0, 0.0), (0.05, 1.0)))
		#
		# Define BCs  # Only defines for last model
			mm.DisplacementBC(name=mlist[i]+'_Bc_z', createStepName='ApplyLoad', region=EdgeSet, u1=UNSET, u2=UNSET, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude='SmoothAmp', fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
			mm.DisplacementBC(name=mlist[i]+'_BC_Top', createStepName='ApplyLoad', region=TopSet, u1=shdisp, u2=disp, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude='SmoothAmp', fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
			mm.DisplacementBC(name=mlist[i]+'_BC_Bot', createStepName='ApplyLoad', region=BotSet, u1=shlock, u2=0.0, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude='SmoothAmp', fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
			
			sv1.setValues(displayedObject=mm.parts[mlist[i]])

			mm.HomogeneousSolidSection(name='Section-1', material='NCF', thickness=None)

			c = p.cells
			ElemRegi = c.getByBoundingBox(-100, -100,-100, 100, 100, 100)
			p.Set(cells=ElemRegi, name='Elm')
			region = regionToolset.Region(cells=ElemRegi)
			p.SectionAssignment(region=region, sectionName='Section-1', offset=0.0, 
					offsetType=MIDDLE_SURFACE, offsetField='', 
					thicknessAssignment=FROM_SECTION)
		#
		# Equation for shear tests
			if 'sh' in mlist[i]:
				nodeSet1 = ai.vertices.getByBoundingBox( p1x-0.1,p1y-0.1,lt-0.1, p2x+0.1,p1y,lt) # 2 nodes in the front
				nodeSet2 = ai.vertices.getByBoundingBox(p1x-0.1,p1y-0.1,-0.1, p1x+0.1,p1y,0.1)    # 1 nodes in the back
				nodeSet3 = ai.vertices.getByBoundingBox(p2x-0.1,p1y-0.1,-0.1, p2x+0.1,p1y,0.1)    # final nodes in the back
				a.Set(vertices=nodeSet1+nodeSet2, name=(mlist[i]+'nodes_3'))    
				a.Set(vertices=nodeSet3, name=(mlist[i]+'nodes_1'))
				mm.Equation(name=mlist[i]+'_eq1', terms=((1.0, mlist[i]+'nodes_3', 2), (-1.0, mlist[i]+'nodes_1', 2)))  # 2 for the y-dir
		#
		# Define fibre orientation
			orientation=None
			p.MaterialOrientation(region=region, 
				orientationType=SYSTEM, axis=axisRot, localCsys=orientation, 
				fieldName='', additionalRotationType=ROTATION_ANGLE, 
				additionalRotationField='', angle=ori, stackDirection=STACK_3)
		#
		# Mesh the cubes with more elements using bounding box method       
			Elem = p.edges.getByBoundingBox(-100, -100,-100, 100, 100, 100)
			p.seedEdgeBySize(edges=Elem , size=lt, deviationFactor=0.1, minSizeFactor=0.1, constraint=FINER)
			p.generateMesh()
		#
		# Declaring fieldOutputRequests SDV   
			fo=mm.fieldOutputRequests['F-Output-1']
			fo.setValues(variables=('S', 'LE', 'U','RF', 'SDV','STATUS','CSQUADSCRT','CSDMG',))
			fo.setValues(numIntervals=200)
		#
		#Creating the job with subroutine 
			mdb.Job(name=mname, model=mname, description='', 
			type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, 
			memory=90, memoryUnits=PERCENTAGE, explicitPrecision=DOUBLE_PLUS_PACK,  #
			nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF, 
			contactPrint=OFF, historyPrint=OFF, #userSubroutine="C:\\SIMULIA\\EstProducts\\2021\\win_b64\\code\\python2.7\\lib\\abaqus_plugins\\SingleElementTest\\subroutine\\SHF_version2.obj",
			scratch='', parallelizationMethodExplicit=DOMAIN, numDomains=1, 
			activateLoadBalancing=False, multiprocessingMode=DEFAULT, numCpus=1)
		# 
		# Regenerate the Assembly
			mm.rootAssembly.regenerate()

		print('#########################')
		print('The single elements have been created!')
		print('#########################')
		print('Copying pre-compiled subroutine')
		shutil.copy2(subDir + "explicitU-D.dll", currentDir)
		print('Please submit the job')
	###############################################################################
	# ###################       Plotting   ########################################           
	###############################################################################

	if not HFrame1 :
		getWarningReply('Please choose "Calibrate the Model", Create CAE model" or "Plot the results"!\n YES or CANCEL to continue', (YES,CANCEL))


	#####
	if HFrame1=='Calibrate the Model':	


	# Step 1 - Import calibration data
		with open(absDir+'\\MatData\\MatProps.json','r') as f:
			props=json.load(f)
		#
		# Step 1 - Create an array for strain
		epsIncs = np.linspace(0.001,0.15,num=100)
		tau12=np.zeros(len(epsIncs))
		dam=np.zeros(len(epsIncs))

		gamma0 = float(props["tau0"])/float(props["G12"])
		epsfKink = 0.6**props["exp"]

		for i in range(0, len(epsIncs)):
			
			# Step 2 - Define damage
			dam[i] = (abs(epsIncs[i])**props["exp"]-gamma0**props["exp"] ) / (epsfKink- gamma0**props["exp"])

			if epsIncs[i] > gammaLcr:
				dam[i] =  1

			dam[i] = min(max(dam[i], 0),1)

			# Step 3 - Friction 
			tauFric = props["MU"]*props["PO"]

			# Step 1 - Calculate the applied stress
			tau12[i] = props["G12"]*epsIncs[i]*(1-dam[i]) + tauFric*dam[i]

		# Step 2 - Calculate damage

		# Step 4 - Write stress back to the file
		np.savetxt(absDir+'\\MatData\\shearCurve.out', np.transpose([epsIncs, tau12]), delimiter=' ')  

		# Step 5 - Read Output file
			# Load data from file to create a XYData
		xy1 = xyPlot.XYDataFromFile(fileName=absDir+'\\MatData\\shearCurve.out', xField=1, yField=2, sourceDescription='Read from'+'shearCurve.out',) 
				#Make a plot: From XYDATA to XYPlots, it cannot be overwritten 
		try:
			session.xyDataObjects.changeKey(xy1.name,  'Calibrated shear' )
		except:
			del session.xyDataObjects['Calibrated shear']
			session.xyDataObjects.changeKey(xy1.name,  'Calibrated shear')
	#
		try:
			xyp = session.XYPlot('XYPlot-1')  	 # round brackets to create the plot
		except:
			xyp = session.xyPlots['XYPlot-1']     		# square brackets to access the plotSSSSSSSSSSS (with s)
	#
	#
	# Stylethe experimental curves	- First bring the XYData to the XYPlot
		xy1 = session.xyDataObjects['Calibrated shear']
		chartName = xyp.charts.keys()[0]
		chart = xyp.charts[chartName]
		c1 = session.Curve(xyData=xy1)
		chart.setValues(curvesToPlot=(c1, ), )
		session.viewports['Viewport: 1'].setValues(displayedObject=xyp)
		session.mdbData.summary()
	#
	# Stylethe experimental curves	
		charts=session.charts['Chart-1']
		charts.axes2[0].labelStyle.setValues(font='-*-verdana-bold-r-normal-*-*-180-*-*-p-*-*-*')
		charts.axes2[0].titleStyle.setValues(font='-*-verdana-bold-r-normal-*-*-180-*-*-p-*-*-*')
		charts.axes2[0].axisData.setValues(useSystemTitle=False, title='Stress (MPa)')
		charts.axes1[0].titleStyle.setValues(font='-*-verdana-bold-r-normal-*-*-180-*-*-p-*-*-*')
		charts.axes1[0].labelStyle.setValues(font='-*-verdana-bold-r-normal-*-*-180-*-*-p-*-*-*')
		charts.axes1[0].axisData.setValues(useSystemTitle=False, title='Strain (-)')
		session.curves['Calibrated shear'].lineStyle.setValues(style=SOLID)		
		session.curves['Calibrated shear'].lineStyle.setValues(thickness=1.2)
		session.curves['Calibrated shear'].lineStyle.setValues(color='#0037fc')
	# Axis
		charts.axes1[0].axisData.setValues(maxValue=0.1, maxAutoCompute=True)
		charts.axes1[0].axisData.setValues(minValue=0, minAutoCompute=True)
		charts.fitCurves(fitAxes1=True, fitAxes2=True)
	#########################################################
		# Shear EXPERIMENTAL curve
	###########################################################
				
		#Make a plot: From XYDATA to XYPlots, it cannot be overwritten 
		# Tries to read an experimental plot from file
		try:
			xy1 = xyPlot.XYDataFromFile(fileName=sh12f, xField=1, yField=2, sourceDescription='Read from'+ str(sh12f),) 
		except:
			getWarningReply('Please add a file with experimental shear values!\n YES or CANCEL to continue', (YES,CANCEL))

		try:
			session.xyDataObjects.changeKey(xy1.name, 'Experimental Shear' )
			xyp = session.XYPlot('XYPlot-1')  	 		# round brackets to create the plot
		except:
			xyp = session.xyPlots['XYPlot-1']     		# square brackets to access the plotSSSSSSSSSSS (with s
	# 
			xy1 = session.xyDataObjects['Experimental Shear']
			chartName = xyp.charts.keys()[0]	
			chart = xyp.charts[chartName]
			c2 = session.Curve(xyData=xy1)
			chart.setValues(curvesToPlot=(c1, c2, ), )	
			session.viewports['Viewport: 1'].setValues(displayedObject=xyp)
			session.mdbData.summary()

			# Read default properties from file
			print('#########################')
			print('The model has been calibrated!')
			print('#########################')

			
	# Initial parameters
	# Step 1 - Import calibration data

		G12 = float(props["G12"])
		tau_0 = float(props["tau0"])
		p = float(props["exp"])
		g_f = 0.6

		p_0 = float(props["PO"])
		mu = float(props["MU"])

		# Precompression
		s22 = -0.


		def solve_curve(G12, p, p_0, mu, tau_0):
			s22 = 0
			g_f = 0.6
			g_0 = tau_0 / G12
			g12 = np.linspace(0., g_f, 2000)
			s12 = []
			dmg = []
			for g_ in g12:
				if g_ < g_0:
					s12.append(g_ * G12)
					dmg.append(0)
				else:
					# Compute damage
					d = (g_ ** p - g_0 ** p) / (g_f ** p - g_0 ** p)
					# Friction
					s12.append((1. - d) * G12 * g_ - d * np.sign(g_) * mu * (s22 - p_0))
					dmg.append(d)

			return g12, s12, dmg


		fig, axes = plt.subplots(ncols=2, figsize=(8, 6), sharex=True)
		ax, ax2 = axes
		plt.subplots_adjust(left=0.10, bottom=0.35, top=0.95)
		g12, s12, dmg = solve_curve( G12,  p, p_0, mu, tau_0)
		# strain12 = np.tan(g12)
		l_, = ax.plot(g12, s12, 'b-', lw=2)
		#ax.plot(g12, s12, '-', color='grey', lw=0.5, label='HTS45/LY556 NCF')
		ax.set_xlabel(r'$\gamma \, (\%)$')
		ax.set_ylabel(r'$\tau \, (MPa)$')
		ax.set_xlim(0, 0.20)
		ax.set_ylim(0, 140)

		# Damage axis
		# ax2 = ax.twinx()
		ax2.set_ylim(0, 1)
		l_2, = ax2.plot(g12, dmg, 'b--', lw=2,  label='Damage evolution')
		# plt.sca(ax)

		axcolor = 'lightgoldenrodyellow'
		#ax_S22 = plt.axes([0.35, 0.02, 0.55, 0.03], facecolor='white')
		ax_p0 = plt.axes([0.35, 0.07, 0.55, 0.03], facecolor=axcolor)
		ax_p = plt.axes([0.35, 0.12, 0.55, 0.03], facecolor=axcolor)
		ax_tau0 = plt.axes([0.35, 0.17, 0.55, 0.03], facecolor=axcolor)
		#ax_gf = plt.axes([0.35, 0.22, 0.55, 0.03], facecolor=axcolor)

		#s_S22 = Slider(ax_S22, r'$\sigma_{22} \, (MPa)$', -100., 100.0, valinit=s22, valstep=1.0)
		s_p0 = Slider(ax_p0, r'$p_0 \, (MPa)$', 0.0, 200.0, valinit=p_0, valstep=1.0)
		s_p = Slider(ax_p, r'$p$', -1.0, 1.0, valinit=p)
		s_tau0 = Slider(ax_tau0, r'$\tau_0 \, (MPa)$', 1, 150, valinit=tau_0, valstep=1.0)
		#s_gf = Slider(ax_gf, r'$\gamma_f$', 0, 2.0, valinit=g_f, valstep=0.05)


		def update(val):
			# print(val)  # Value of the slider changed
			#gf_i = s_gf.val
			tau0_i = s_tau0.val
			p_i = s_p.val
			p0_i = s_p0.val
			#s22_i = s_S22.val

			mu_i = float(r_mu.value_selected)

			#g12, s12, dmg = solve_curve(gf_i, G12, s22_i, p_i, p0_i, mu_i, tau0_i)
			g12, s12, dmg = solve_curve(G12, p_i, p0_i, mu_i, tau0_i)

			#strain12 = np.tan(g12)
			l_.set_xdata(g12)
			l_.set_ydata(s12)
			l_2.set_xdata(g12)
			l_2.set_ydata(dmg)


		#s_S22.on_changed(update)
		s_p0.on_changed(update)
		s_p.on_changed(update)
		s_tau0.on_changed(update)
		#s_gf.on_changed(update)

		ax_mu = plt.axes([0.10, 0.02, 0.12, 0.23], facecolor=axcolor)
		r_mu = RadioButtons(ax_mu, ('0.6', '0.4', '0.2', '0.1', '0'), active=1, activecolor='grey')

		colors = {'0.6': 'red', '0.4': 'black', '0.2': 'blue', '0.1': 'yellow', '0': 'green'}

		def colorfunc(label):
			l_.set_color(colors[label])
			l_2.set_color(colors[label])
			# r_mu.activecolor = colors[label]
			update(0)
			fig.canvas.draw_idle()

		r_mu.on_clicked(colorfunc)

		colorfunc(r_mu.value_selected)

		data = np.loadtxt(sh12f, delimiter=',')
		gamma_ = data[:, 0]
		tau_ = data[:, 1]
		# Exp. plot
		ax.plot(gamma_, tau_, '-', color='cyan', label='HTS45/LY556 NCF')
		ax2.plot(g12, dmg, 'b--', lw=2,  label='Damage evolution')


		leg = ax.legend()
		#leg._draggable(True)

		plt.show()

###########################
###########################
	if HFrame1=='Plot the results':

		with open(absDir+'\\MatData\\MatProps.json','r') as f:
			props=json.load(f)

		PartInstance_label_component = {'FC-1':'S11','FC10-1':'S11','FC20-1':'S11','FC30-1':'S11','FC45-1':'S11','FT-1':'S11','MC-1':'S22','MC3-1':'S33','MT-1':'S22','SH12-1':'S12','SH13-1':'S13','SH23-1':'S23',
			}
		
		odb_path = session.viewports[session.currentViewportName].displayedObject.name  # Returns the full path with the name of the .odb otherwise XYPLOT
		odb = session.odbs[odb_path]
		step = odb.steps['ApplyLoad']

		# Initialize a list to store data from all frames
		stress_results = {'FC-1':[],'FC10-1':[],'FC20-1':[],'FC30-1':[],'FC45-1':[],'FT-1':[],'MC-1':[],'MC3-1':[],'MT-1':[],'SH12-1':[],'SH13-1':[],'SH23-1':[],}
		strain_results = {'FC-1':[],'FC10-1':[],'FC20-1':[],'FC30-1':[],'FC45-1':[],'FT-1':[],'MC-1':[],'MC3-1':[],'MT-1':[],'SH12-1':[],'SH13-1':[],'SH23-1':[],}

		engStrain = [0,-1.24064e-5,-0.98506e-5,-33.0083e-5,-77.6423e-5,-150.485e-5,-257.996e-5,-406.556e-5,-602.214e-5,-850.959e-5,-0.0115826,-0.0152967,-0.0197046,-0.0248567,-0.0308019,-0.0375867,-0.0452559,-0.0538523,-0.0634168,-0.0739886,-0.0856,-0.0982957,-0.112105,-0.127059,-0.143189,-0.160522,-0.179086,-0.198905,-0.220003,-0.242401,-0.266119,-0.291175,-0.317587,-0.34537,-0.374543,-0.405108,-0.43708,-0.470468,-0.505282,-0.541526,-0.579207,-0.618328,-0.658891,-0.700897,-0.744346,-0.789237,-0.835565,-0.883329,-0.932521,-0.983136,-1.03517,-1.0886,-1.14343,-1.19965,-1.25724,-1.31619,-1.37649,-1.43812,-1.50106,-1.5653,-1.63081,-1.69759,-1.76561,-1.83484,-1.90528,-1.97688,-2.04964,-2.12353,-2.19852,-2.27459,-2.35171,-2.42985,-2.50899,-2.58909,-2.67014,-2.75209,-2.83492,-2.9186,-3.0031,-3.08838,-3.17442,-3.26117,-3.34861,-3.43671,-3.52542,-3.61472,-3.70457,-3.79493,-3.88577,-3.97706,-4.06875,-4.16081,-4.25321,-4.34591,-4.43887,-4.53205,-4.62542,-4.71894,-4.81257,-4.90627,-5.00002,-5.09376,-5.18747,-5.2811,-5.37462,-5.46799,-5.56117,-5.65413,-5.74682,-5.83922,-5.93129,-6.02298,-6.11426,-6.20511,-6.29546,-6.38533,-6.47463,-6.56333,-6.65142,-6.73885,-6.8256,-6.91164,-6.99692,-7.08142,-7.16511,-7.24794,-7.32989,-7.41095,-7.49105,-7.57019,-7.64833,-7.72545,-7.80151,-7.8765,-7.95039,-8.02315,-8.09476,-8.16519,-8.23442,-8.30244,-8.36921,-8.43473,-8.49897,-8.56191,-8.62354,-8.68383,-8.74278,-8.80037,-8.85659,-8.91141,-8.96486,-9.01689,-9.06751,-9.1167,-9.16446,-9.21078,-9.25567,-9.29912,-9.34112,-9.38169,-9.4208,-9.45849,-9.49473,-9.52954,-9.56293,-9.5949,-9.62547,-9.65463,-9.68242,-9.70883,-9.73388,-9.7576,-9.78,-9.8011,-9.82092,-9.83948,-9.85681,-9.87294,-9.8879,-9.90171,-9.9144,-9.92602,-9.93659,-9.94615,-9.95475,-9.96242,-9.9692,-9.97515,-9.9803,-9.98471,-9.98842,-9.99149,-9.99398,-9.99593,-9.99742,-9.9985,-9.99922,-9.99967,-9.9999,-9.99999,-10]
		
		for label, component in PartInstance_label_component.items():
			variable ="S"

			xyList = xyPlot.xyDataListFromField(odb=odb, outputPosition=ELEMENT_CENTROID, 
				variable=((variable, INTEGRATION_POINT, ((COMPONENT, component), )), ), 
				elementLabels=((label, ('1', )), ))
			
			if label == "SH12-1": component, variable = "LE12", "LE"
			if label == "SH13-1": component, variable = "LE13", "LE"
			if label == "SH23-1": component, variable = "LE23", "LE"

			xyListStrain = xyPlot.xyDataListFromField(odb=odb, outputPosition=ELEMENT_CENTROID, 
					variable=((variable, INTEGRATION_POINT, ((COMPONENT, component), )), ), 
					elementLabels=((label, ('1', )), ))
						
			strain_results[label].append(xyListStrain[0].data)
			stress_results[label].append(xyList[0].data)
	
				
		# Extract time and stress data for plotting
		#time_fc1 = [data[1] for data in strain_results['FC-1'][0]]
		FC = [data[1] for data in stress_results['FC-1'][0]]
		FC10 = [data[1] for data in stress_results['FC10-1'][0]]
		FC20 = [data[1] for data in stress_results['FC20-1'][0]]
		FC30 = [data[1] for data in stress_results['FC30-1'][0]]
		FC45 = [data[1] for data in stress_results['FC45-1'][0]]
		FT = [data[1] for data in stress_results['FT-1'][0]]
		MC = [data[1] for data in stress_results['MC-1'][0]]
		MC3 = [data[1] for data in stress_results['MC3-1'][0]]
		MT = [data[1] for data in stress_results['MT-1'][0]]
		SH12 = [data[1] for data in stress_results['SH12-1'][0]]
		SH13 = [data[1] for data in stress_results['SH13-1'][0]]
		SH23 = [data[1] for data in stress_results['SH23-1'][0]]

		SH12_strain = [data[1] for data in strain_results['SH12-1'][0]]
		SH13_strain = [data[1] for data in strain_results['SH13-1'][0]]
		SH23_strain = [data[1] for data in strain_results['SH23-1'][0]]

		colors = ['#1f77b4', '#ff7f0e', 'grey', '#d62728',  'black']

		##### START PLOTTING ######
		# Tensile Loading 
		plt.figure(figsize=(10, 6))
		plt.plot(np.abs(engStrain[:len(FT)]), FT, label= "Fibre tension" , color=colors[0])
		plt.plot(np.abs(engStrain[:len(MT)]), MT, label= "Matrix tension" , color=colors[1])
		plt.scatter(props['XT_eps'],props['XT'], color=colors[0], marker='x', s=100, label='Longitudinal tensile strength')
		plt.scatter(props['YT_eps'],props['YT'], color=colors[1], marker='x', s=100, label='Transverse tensile strength')
		# Compressive loading 
		plt.plot(engStrain[:len(FC)], FC, label= "Fibre compression", color= colors[2])
		plt.plot(engStrain[:len(MC)], MC, label= "Matrix compression" , color= colors[3] )
		plt.scatter(props['XC_eps'], props['XC'], color=colors[2], marker='x', s=100, label='Longitudinal compressive strength')
		plt.scatter(props['YC_eps'],props['YC'], color=colors[3], marker='x', s=100, label='Transverse compressive strength')

		# Axis labels and legend
		plt.xlabel('Engineering strain (%)', fontsize=18)
		plt.ylabel('Stress (MPa)', fontsize=18)
		plt.legend(fontsize=13)
		# Grid and ticks
		plt.grid(True, linestyle='--', alpha=0.7)
		plt.xticks(fontsize=16)
		plt.yticks(fontsize=16)
		plt.show()

		###### Shear loading 
		plt.figure(figsize=(10, 6))
		half_index = len(SH12_strain) // 2

		plt.plot(np.abs(SH12_strain[:half_index])*100, np.abs(SH12[:half_index]), label= "In-plane shear", color= colors[0])
		plt.plot(np.abs(SH13_strain[:half_index])*100, np.abs(SH13[:half_index]), label= "Out-off-plane shear (1-3)", color= colors[1])
		plt.plot(np.abs(SH23_strain[:half_index])*100, np.abs(SH23[:half_index]), label= "Out-off-plane shear (2-3)", color= colors[2])

		plt.axhline(y=props['SL'], color=colors[0], linewidth=3, alpha=0.2, label='Experimental shear strength')

		# Axis labels and legend
		plt.xlabel('Engineering strain (%)', fontsize=18)
		plt.ylabel('Shear stress (MPa)', fontsize=18)
		plt.legend(fontsize=13)


		# Grid and ticks
		plt.grid(True, linestyle='--', alpha=0.7)
		plt.xticks(fontsize=16)
		plt.yticks(fontsize=16)
		plt.show()

		# Off-axis
		plt.figure(figsize=(10, 6))

		plt.plot(engStrain[:len(FC10)], FC10, label=r" $\theta_i$=" + str(theta_i + props['fc10'] ) + r"$^\circ$", color=colors[0], linewidth=2)
		plt.plot(engStrain[:len(FC)], FC, label=r" $\theta_i$="  + str(theta_i)  + r"$^\circ$", color=colors[2],alpha=0.6, linewidth=2)
		plt.plot(engStrain[:len(FC20)], FC20, label=r" $\theta_i$=" + str(theta_i + props['fc20'] ) + r"$^\circ$", color=colors[1], linewidth=2)
		plt.plot(engStrain[:len(FC30)], FC30, label=r" $\theta_i$=" + str(theta_i + props['fc30'] ) + r"$^\circ$", color=colors[4], linewidth=2)
		plt.plot(engStrain[:len(FC45)], FC45, label=r" $\theta_i$=" + str(theta_i + props['fc45'] ) + r"$^\circ$", color=colors[3], linewidth=2)

		plt.axhline(y=props['XC'], color=colors[4], linewidth=4, alpha=0.6, label='Experimental strength')

		# Axis labels and legend
		plt.xlabel('Engineering strain (%)', fontsize=18)
		plt.ylabel('Stress (MPa)', fontsize=18)
		plt.legend(fontsize=13)

		# Grid and ticks
		plt.grid(True, linestyle='--', alpha=0.7)
		plt.xticks(fontsize=16)
		plt.yticks(fontsize=16)
		plt.show()


		
"""


# Before starting the model type the following command into the Python editor box in CAE
# session.journalOptions.setValues(recoverGeometry=COORDINATE)

execfile('C:\\SIMULIA\\EstProducts\\2021\\win_b64\\code\\python2.7\\lib\\abaqus_plugins\\SingleElementTest\\CreateSingleElem.py', __main__.__dict__)

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
YC=120
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
thetai_12=     0
thetai_13=     2
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
HFrame1 =  'Plot the results' # 'Create single elements' #'Calibrate the Model' #
currentDir=os.getcwd()
odb=currentDir+'\\SingleElemTest-.odb'
absDir='C:\\SIMULIA\\EstProducts\\2021\\win_b64\\code\\python2.7\\lib\\abaqus_plugins\\SingleElementTest\\'
eps2cr = 0.09
eps3cr = 0.09

BASE_DIR = os.getcwd()
PARENT_DIR = os.path.dirname(BASE_DIR)
matProps = 'C:\\SIMULIA\\EstProducts\\2021\\win_b64\\code\\python2.7\\lib\\abaqus_plugins\\SingleElementTest\\MatData\\SHF_version2.for'

with open('C:\\SIMULIA\\EstProducts\\2021\\win_b64\\code\\python2.7\\lib\\abaqus_plugins\\SingleElementTest\\MatData\\'+'MatProps.json','r') as f:
	props=json.load(f)
# 
"""