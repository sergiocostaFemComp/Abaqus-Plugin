
if HFrame1=='Plot the results':
	# Plot 


	print('*******')
	print("PLOTTING THE RESULTS!")
	print('*******')

	#
	# Open the ODB
	#myOdb = visualization.openOdb(path=odb)            # Open the odb
	#session.mdbData.summary()

	try:  # if the display of odb is active
		odb_path = session.viewports[session.currentViewportName].displayedObject.name  # Returns the full path with the name of the .odb otherwise XYPLOT
		odb = session.odbs[odb_path]
		odb_name = os.path.split(session.viewports[session.currentViewportName].displayedObject.path)[1].split(".")[0]  # Returns the name of the .odb only
	except:
		getWarningReply('The ODB with the single elements is not available in the display!\nPlease open the ODB or change the display', (YES,CANCEL))
		odb_name = odb.split("/")[-1].split(".")[0]

		#pass  # takes the odb in the script
		
	#session.viewports[session.currentViewportName].setValues(displayedObject=odb)

	mlist=[];		mlstFE=[];		elmSets=[]
	stressComp=[];	strainComp=[]
	expFile=[];		curveColor=[]
	#
	# For experimental models
	if MC_Exp==True: mlist.append('MC_Exp');			expFile.append(mcf);		curveColor.append('#000000')
	if FC_Exp==True: mlist.append('FC_Exp');			expFile.append(fcf);    	curveColor.append('#c93434')
	if MT_Exp==True: mlist.append('MT_Exp');			expFile.append(mtf);    	curveColor.append('#475757')
	if FT_Exp==True: mlist.append('FT_Exp');			expFile.append(ftf);    	curveColor.append('#ba6e6e')
	if SH12_Exp==True: mlist.append('SH12_Exp');		expFile.append(sh12f);  	curveColor.append('#0037fc')
	if SH13_Exp==True: mlist.append('SH13_Exp');		expFile.append(sh13f);		curveColor.append('#2dad8f')	
	if SH23_Exp==True: mlist.append('SH23_Exp');		expFile.append(sh23f);  	curveColor.append('#FFFF00')
	#
	mc3=False
	# For simulation models
	if mc==True: mlstFE.append('mc');		elmSets.append('MC-1.ELM'); 	stressComp.append('S22'); 	strainComp.append('LE22');	expFile.append(mcf);	curveColor.append('#000000')
	if mc3==True: mlstFE.append('mc3');		elmSets.append('MC3-1.ELM'); 	stressComp.append('S33'); 	strainComp.append('LE33');	expFile.append(mcf);  	curveColor.append('#666666')
	if fc==True: mlstFE.append('fc');		elmSets.append('FC-1.ELM');		stressComp.append('S11');	strainComp.append('LE11'); 	expFile.append(fcf);    	curveColor.append('#c93434')
	if mt==True: mlstFE.append('mt');		elmSets.append('MT-1.ELM');		stressComp.append('S22');	strainComp.append('LE22'); 	expFile.append(mtf);    	curveColor.append('#475757')
	if ft==True: mlstFE.append('ft');		elmSets.append('FT-1.ELM');		stressComp.append('S11');	strainComp.append('LE11'); 	expFile.append(ftf);    	curveColor.append('#ba6e6e')
	if sh12==True: mlstFE.append('sh12');	elmSets.append('SH12-1.ELM');	stressComp.append('S12');	strainComp.append('LE12'); 	expFile.append(sh12f);  	curveColor.append('#0037fc')
	if sh13==True: mlstFE.append('sh13');	elmSets.append('SH13-1.ELM');	stressComp.append('S13');	strainComp.append('LE13'); 	expFile.append(sh13f);	curveColor.append('#2dad8f')
	if sh23==True: mlstFE.append('sh23');	elmSets.append('SH23-1.ELM');	stressComp.append('S23');	strainComp.append('LE23'); 	expFile.append(sh23f);  	curveColor.append('#a82893')

	# Fix graphics according to precision
	# if precision==SINGLE: 
	# 	mystyle=SOLID
	# 	mythick = 0.5
	# else:
	# 	mystyle=DOTTED
	# 	mythick =1.2
	# FOr the experiments
	for i in range(0, len(mlist)):
	#
	# Load data from file to create a XYData
		xy1 = xyPlot.XYDataFromFile(fileName=str(expFile[i]), xField=1, yField=2, sourceDescription='Read from'+ str(expFile[i]),) 
	#
	#Make a plot: From XYDATA to XYPlots, it cannot be overwritten 
		try:
			expType=session.xyDataObjects.changeKey(xy1.name,  mlist[i] )
		except:
			usrInput=getWarningReply('The plot ' + mlist[i] +  ' already exist and it will be replaced!\nOkay to continue', (YES,NO))
			if usrInput=='No':break
			del session.xyDataObjects[mlist[i]]
			session.xyDataObjects.changeKey(xy1.name,  mlist[i])
		try:
			xyp = session.XYPlot('XYPlot-1')  	 # round brackets to create the plot
		except:
			xyp = session.xyPlots['XYPlot-1']     		# square brackets to access the plotSSSSSSSSSSS (with s)
	#
	# Style the experimental curves	- First bring the XYData to the XYPlot
		xy1 = session.xyDataObjects[ mlist[i]]
		chartName = xyp.charts.keys()[0]
		chart = xyp.charts[chartName]
		c1 = session.Curve(xyData=xy1)
		chart.setValues(curvesToPlot=(c1, ), )
		session.viewports['Viewport: 1'].setValues(displayedObject=xyp)
		session.mdbData.summary()
	#
	# Style the experimental curves	
		charts=session.charts['Chart-1']
		charts.axes2[0].labelStyle.setValues(font='-*-verdana-bold-r-normal-*-*-180-*-*-p-*-*-*')
		charts.axes2[0].titleStyle.setValues(font='-*-verdana-bold-r-normal-*-*-180-*-*-p-*-*-*')
		charts.axes2[0].axisData.setValues(useSystemTitle=False, title='Stress (MPa)')
		charts.axes1[0].titleStyle.setValues(font='-*-verdana-bold-r-normal-*-*-180-*-*-p-*-*-*')
		charts.axes1[0].labelStyle.setValues(font='-*-verdana-bold-r-normal-*-*-180-*-*-p-*-*-*')
		charts.axes1[0].axisData.setValues(useSystemTitle=False, title='Strain (-)')
		session.curves[ mlist[i]].lineStyle.setValues(style=SOLID)		
		session.curves[ mlist[i]].lineStyle.setValues(thickness=1.2)
		session.curves[ mlist[i]].lineStyle.setValues(color=curveColor[i])
	# Axis
		charts.axes1[0].axisData.setValues(maxValue=0.1, maxAutoCompute=False)
		charts.axes1[0].axisData.setValues(minValue=-0.1, minAutoCompute=False)
		charts.fitCurves(fitAxes1=True, fitAxes2=False)
	# 
	#################################################################
	##############################    FE     ########################
	#################################################################
	for i in range(0, len(mlstFE)):
	# Extract the respective strain and stress
		strainStress  = xyPlot.xyDataListFromField(odb=odb, outputPosition=INTEGRATION_POINT, variable=(
			('LE', INTEGRATION_POINT, ((COMPONENT, strainComp[i]), )),
			('S',   INTEGRATION_POINT, ((COMPONENT, stressComp[i]), ))),  
			elementSets=(elmSets[i], ))
	# #
	# # Combine the stress with the strain # strainStress[0].legendLabel
		x1 = session.xyDataObjects[strainStress[0].name]
		y1 = session.xyDataObjects[strainStress[1].name]
		xy1 = combine(x1, y1)
	# #
	# # Add the combined plot to the session
		try:
			session.xyDataObjects.changeKey(xy1.name, mlstFE[i]+'_FE-'+ odb_name)
		except:
			usrInput=getWarningReply('The plot ' + mlstFE[i]+'_FE-' + odb_name+ ' already exist and it will be replaced!\nOkay to continue', (YES,NO))
			if usrInput=='No':break
			del session.xyDataObjects[mlstFE[i]+'_FE-']
			session.xyDataObjects.changeKey(xy1.name,  mlstFE[i]+'_FE-'+odb_name)
		xyp = session.xyPlots['XYPlot-1']
		chartName = xyp.charts.keys()[0]
		chart = xyp.charts[chartName]
		xy1 = session.xyDataObjects[ mlstFE[i]+'_FE-'+odb_name]
		c1 = session.Curve(xyData=xy1)
		chart.setValues(curvesToPlot=(c1, ), )
		session.viewports['Viewport: 1'].setValues(displayedObject=xyp)
	#
	# Style the FE curves	
		charts=session.charts['Chart-1']
		charts.axes2[0].labelStyle.setValues(font='-*-verdana-bold-r-normal-*-*-180-*-*-p-*-*-*')
		charts.axes2[0].titleStyle.setValues(font='-*-verdana-bold-r-normal-*-*-180-*-*-p-*-*-*')
		charts.axes2[0].axisData.setValues(useSystemTitle=False, title='Stress (MPa)')
		charts.axes1[0].titleStyle.setValues(font='-*-verdana-bold-r-normal-*-*-180-*-*-p-*-*-*')
		charts.axes1[0].labelStyle.setValues(font='-*-verdana-bold-r-normal-*-*-180-*-*-p-*-*-*')
		charts.axes1[0].axisData.setValues(useSystemTitle=False, title='Strain (-)')
		session.curves[mlstFE[i]+'_FE-'+odb_name].lineStyle.setValues(style=SOLID)		
		session.curves[mlstFE[i]+'_FE-'+odb_name].lineStyle.setValues(thickness=mythick)
		session.curves[ mlstFE[i]+'_FE-'+odb_name].lineStyle.setValues(color=curveColor[i])
		session.charts['Chart-1'].legend.textStyle.setValues(color=curveColor[i])
		session.charts['Chart-1'].legend.titleStyle.setValues(color=curveColor[i])

	# Axis
		charts.axes1[0].axisData.setValues(maxValue=0.1, maxAutoCompute=False)
		charts.axes1[0].axisData.setValues(minValue=-0.1, minAutoCompute=False)
		charts.fitCurves(fitAxes1=True, fitAxes2=False)


		C:/SIMULIA/EstProducts/2021/win_b64/code/python2.7/lib/abaqus_plugins/SingleElementTest/MatData/sh12Exp.txt", "density": 1.594e-09, "G13": 3770, "eps3cr": 0.0381, "tau0": 30, "GIc_mt": 0.149, "PO": 40, "matProps": "C:\\SIMULIA\\EstProducts\\2021\\win_b64\\code\\python2.7\\lib\\abaqus_plugins\\SingleElementTest\\MatData\\MatProps.json",