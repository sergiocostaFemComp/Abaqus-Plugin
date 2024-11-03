import section
import regionToolset
import displayGroupMdbToolset as dgm
import job
import part
import material
import assembly
import optimization
import step
import interaction
import load
import mesh
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior
"""
execfile('C:\\Users\\soca\\abaqus_plugins\\SingleElementTest\\loadCurves.py', __main__.__dict__)
"""
mc=True
fc=True
mt=True
ft=True
sh12=True
sh13=True
sh23=True

mlist=[]
elmSets=[]
stressComp=[]
strainComp=[]
if mc==True: mlist.append('mc');		 elmSets.append('MC-1.ELM'); 	stressComp.append('S22'); 		strainComp.append('LE22')
if fc==True: mlist.append('fc');		 elmSets.append('FC-1.ELM');		stressComp.append('S11');		strainComp.append('LE11')
if mt==True: mlist.append('mt');		elmSets.append('MT-1.ELM');		stressComp.append('S22');	strainComp.append('LE22')
if ft==True: mlist.append('ft');		elmSets.append('FT-1.ELM');		stressComp.append('S11');	strainComp.append('LE11')
if sh12==True: mlist.append('sh12');	elmSets.append('SH12-1.ELM');	stressComp.append('S12');	strainComp.append('LE12')
if sh13==True: mlist.append('sh13');	elmSets.append('SH13-1.ELM');	stressComp.append('S13');	strainComp.append('LE13')		
if sh23==True: mlist.append('sh23');	elmSets.append('SH23-1.ELM');	stressComp.append('S23');	strainComp.append('LE23')



# Load data from file
xy1 = xyPlot.XYDataFromFile(fileName='C:/0_Abaqus/99_Colaboration/Engenuity/Results/shearResponse.txt', xField=1, yField=2, sourceDescription='Read from C:/0_Abaqus/99_Colaboration/Engenuity/Results/shearResponse.txt',) 

# find current viewport
currentViewport = session.viewports[session.currentViewportName]
# assign odb file from current viewport
odbFile = currentViewport.displayedObject
# get file name and path
odbFileNameFull = odbFile.path
odb = session.odbs[odbFileNameFull]

for i in range(0, len(mlist)):
	session.viewports[session.currentViewportName].setValues(displayedObject=session.openOdb(name=odbFileNameFull))															
	stress  = xyPlot.xyDataListFromField(odb=odb, outputPosition=INTEGRATION_POINT, variable=(('S', INTEGRATION_POINT, ((COMPONENT, stressComp[i]), )), ),  elementSets=(elmSets[i], ))
	strain  = xyPlot.xyDataListFromField(odb=odb, outputPosition=INTEGRATION_POINT, variable=(('LE', INTEGRATION_POINT, ((COMPONENT, strainComp[i]), )), ),  elementSets=(elmSets[i], ))

	plotStressName = str(stress).split("'")[1]
	plotStrainName = str(strain).split("'")[1]
	x1=session.xyDataObjects[plotStrainName]
	y1=session.xyDataObjects[plotStressName]
	xy2 = combine(x1, y1)
	session.xyDataObjects.changeKey(xy2.name, mlist[i])
