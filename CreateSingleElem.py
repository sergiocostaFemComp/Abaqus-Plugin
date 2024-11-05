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
###############################################################################################################
	if HFrame1=='Create single elements':
###############################################################################################################
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

	if not HFrame1 :
		getWarningReply('Please choose "Calibrate the Model", Create CAE model" or "Plot the results"!\n YES or CANCEL to continue', (YES,CANCEL))

###############################################################################################################
	if HFrame1=='Calibrate the Model':	
###############################################################################################################

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
		theta_i = float(props["thetai_12"]) # Degrees
		p_0 = float(props["PO"])
		mu = float(props["MU"])
		g_f = 0.6
		sigY = 0
		tauXY = 0


		def solve_curve(g_f, G12, p, p_0, mu, tau_0, theta_i, sigY, tauXY):
			s22 = 0
			#g_f = 0.6
			mu = 0.4
			eps12m_0 = tau_0 / G12
			g12 = np.linspace(0., g_f, 2000)
			s12 = []
			dmg = []
			sig11 = []
			eps11 = []
			for eps12m in g12:
				if eps12m < eps12m_0:
					tau12m = eps12m * G12
					s12.append(tau12m)
					dmg.append(0)
				else:
					# Compute damage
					d = (eps12m ** p - eps12m_0 ** p) / (g_f ** p - eps12m_0 ** p)
					# Friction
					tau12m =((1. - d) * G12 * eps12m - d * np.sign(eps12m) * mu * (s22 - p_0))
					s12.append(tau12m)
					#dmg.append(d)

				# Simplified kinking
				theta = eps12m + theta_i * np.pi / 180
				sigX = sigY + 2*(tauXY*np.cos(2*theta) - tau12m )  / np.sin(2*theta)  
				sig11.append((sigX)) # Because sig11 will be negative since it is compressive
				#sig11.append (tau12m / (theta) )
				eps11.append ( eps12m * np.sin( theta ) )

			return g12, s12, eps11, sig11

		fig, axes = plt.subplots(ncols=2, figsize=(8, 6), sharex=True)
		ax, ax2 = axes
		plt.subplots_adjust(left=0.10, bottom=0.35, top=0.95)
		g12, s12, eps11, sig11 = solve_curve(g_f, G12, p, p_0, mu, tau_0, theta_i, sigY, tauXY)

		# Plot the curve with the initial values
		ax.plot(g12, s12, '-', color='blue', lw=2, label='Initial shear calibration')
		ax.set_xlabel(r'$\gamma_{12} \, (-)$', fontsize=14)
		ax.set_ylabel(r'$\tau_{12} \, (MPa)$', fontsize=14)
		ax.set_xlim(0, 0.10)
		ax.set_ylim(0, 100)

		# Kinking Axis (secondary axis)
		ax2.set_ylim(0, -1200)
		ax2.plot(eps11, sig11, '-', color='blue', lw=2, label='Resulting kinking stress')
		ax2.set_xlabel(r'$\epsilon_{11} \, (-)$', fontsize=14)

		# Add the legend for each axis
		ax.legend(loc='upper right')
		ax2.legend(loc='upper right')

		# Slider and Update Function (unchanged)
		axcolor = 'lightgoldenrodyellow'

		# Define sliders
		ax_p0 = plt.axes([0.2, 0.02, 0.7, 0.03], facecolor=axcolor)
		ax_p = plt.axes([0.2, 0.06, 0.7, 0.03], facecolor=axcolor)
		ax_tau0 = plt.axes([0.2, 0.10, 0.7, 0.03], facecolor=axcolor)
		ax_theta_i = plt.axes([0.2, 0.22, 0.7, 0.03], facecolor=axcolor)
		ax_S22 = plt.axes([0.2, 0.14, 0.7, 0.03], facecolor='white')
		ax_t22 = plt.axes([0.2, 0.18, 0.7, 0.03], facecolor='white')

		s_p0 = Slider(ax_p0, r'$p_0 \, (MPa)$', 0.0, 200.0, valinit=p_0, valstep=1.0)
		s_p = Slider(ax_p, r'$p \, (-)$', -1.0, 1.0, valinit=p)
		s_tau0 = Slider(ax_tau0, r'$\tau_0 \, (MPa)$', 1, 150, valinit=tau_0, valstep=1.0)
		t = Slider(ax_theta_i, r'$\theta_i \, (Deg.)$', 0.3, 45, valinit=theta_i)
		s_S22 = Slider(ax_S22, r'$\sigma_{y} \, (MPa)$', -100., 100.0, valinit=sigY, valstep=1.0)
		t_S22 = Slider(ax_t22, r'$\tau_{xy} \, (MPa)$', -100., 100.0, valinit=tauXY, valstep=1.0)


		# For the updated color
		l_, = ax.plot(g12, s12, 'b-', lw=2, color='red', label='Updated shear calibration')
		l_2, = ax2.plot(eps11, sig11, '--', lw=2, color='red', label='Updated kinking stress')


		def update(val):
			tau0_i = s_tau0.val
			p_i = s_p.val
			p0_i = s_p0.val
			theta_i_i = t.val
			sigY_i = s_S22.val
			tauXY_i = t_S22.val

			gf_i = 0.6
			mu_i = 0.4

			g12, s12, eps11, sig11 = solve_curve(gf_i, G12, p_i, p0_i, mu_i, tau0_i, theta_i_i, sigY_i, tauXY_i)

			l_.set_xdata(g12)
			l_.set_ydata(s12)
			l_2.set_xdata(eps11)
			l_2.set_ydata(sig11)

			fig.canvas.draw_idle()

		s_p0.on_changed(update)
		s_p.on_changed(update)
		s_tau0.on_changed(update)
		t.on_changed(update)
		s_S22.on_changed(update)
		t_S22.on_changed(update)

		Exp_shear = np.loadtxt(sh12f, delimiter=',')
		gamma_ = Exp_shear[:, 0]
		tau_ = Exp_shear[:, 1]

		# Exp. plot
		ax.plot(gamma_, tau_, '-', color='black', label='Exp: NCF HTS45/LY556')

		# Add the legend again to account for the experimental data
		ax.legend(loc='upper right')
		ax2.legend(loc='upper right')

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