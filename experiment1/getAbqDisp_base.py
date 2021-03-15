#set work directory first
#execfile('DisplacementExtraction.py')

##### NOTE TO SELF: NEED TO USE DATADOUBLE FOR EXPLICIT AND DATA FOR IMPLICIT

from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup
import numpy as np
executeOnCaeStartup()
import os
import time

#Attempt to define geometry set which auto defines element and node set for the necked region
#von Mises strain equation

###### CHECK FOR ABAQUS ERROR TO SEE IF ANY LICENSES ARE AVAILABLE


#==============================================================================
#  INPUT PARAMETERS
#==============================================================================
NeckNodeSet = NeckNodeSet_par
NeckElemSet = NeckElemSet_par
ForceSet = ForceSet_par
JobName = JobName_par
LoadDir = LoadDir_par
DataPrec = DataPrec_par
RawDir = RawDir_par
DispScaleX = float(DispScaleX_par)
DispScaleY = float(DispScaleY_par)

fieldList = ['COORD','U']#Nodal data to be read
fieldNameList = ['','u']#Append the axis to this string as file identifier
elemFieldList = ['COORD','LE']#Element data to be read
elemFieldNameList = ['coordinate','strain']#descriptions od element data. Used in header of output file
axisList = ['x','y','z']#List of axis, mostly here to control naming
tensorList = ['11','22','33','12','13','23'] #Need to double check this for order

#!!! STRAIN MEASURES NOT NECESSARILY VALID !!!


#Get working directory and seperate all of the .odb files out of it
worDir = os.getcwd()
fileDir = worDir + os.path.sep + RawDir
odb = session.openOdb(name = worDir + '/' + JobName + '.odb')

#Custom function to only search one level of directories
def walklevel(some_dir, level):
    some_dir = some_dir.rstrip(os.path.sep)
    assert os.path.isdir(some_dir)
    num_sep = some_dir.count(os.path.sep)
    for root, dirs, files in os.walk(some_dir):
        yield root, dirs, files
        num_sep_this = root.count(os.path.sep)
        if num_sep + level <= num_sep_this:
            del dirs[:]

def int2strPad(integer,padLength=3):
    outStr = ''
    if integer >= 10**padLength: print('Padding length may be insufficient to ensure ordering')
    if integer == 0: outStr = str('0'*padLength)
    else:
        for i in reversed(range(padLength+1)):
            if integer < 10**i and integer >= 10**(i-1):
               outStr = str('0'*(padLength-i) + str(integer))
    return outStr

def UpdateLogfile(message):#Consider passing urgency to each message so amount of info can be controlled
    path = worDir.replace(os.path.sep + 'temp' + os.path.sep + JobName,'')
    log = open(path + os.path.sep + JobName+'_log.txt', 'a+')
    log.write(message + ': ' + str(time.strftime("%c")) + '\n')
    log.close()

def DefineSubset(setName,setType):#Set type must be (n)ode or (e)lement
    if setType == 'n':
        try:
            tempSet = odb.rootAssembly.nodeSets[setName]
            #UpdateLogfile('Only points in node set '+ setName + ' will be tested')
        except KeyError:
            UpdateLogfile('No node subset found of name '+ setName + ', collecting data for whole database')
            tempSet = None
    if setType == 'e':
        try:
            tempSet = odb.rootAssembly.elementSets[setName]
            #UpdateLogfile('Only points in element set '+ setName + ' will be tested')
        except KeyError:
            UpdateLogfile('No element subset found of name '+ setName + ', collecting data for whole database')
            tempSet = None      
    return tempSet

def ReadLDCData(filename):
    tempDict = {}
    lineID = ['ID','force','extenX','extenY','extenZ','length','pos1x','pos1y','pos1z','pos2x','pos2y','pos2z']
    #Read from single file for now, likely will need to change
    with open(filename,'r') as file:
        for lCount,line in enumerate(file): 
            temp = line.split(',')
            temp = [x.strip() for x in temp]
            for i,val in enumerate(temp):
                if not val: temp[i] = 'nan'
            temp = [float(x) for x in temp[4:]]
            tempDict[lineID[lCount]] = temp
    return tempDict

def ReadLDCFile_txt(filename):
    # this is the implementation for a .txt file
    array1 = []
    array2 = []
    with open(filename, "r") as f:
        for line in f:
            temp1 = line.split('\t')
            temp2 = temp1[1].split('\n')[0]
#            if np.isnan(temp1[0]) == False and np.isnan(temp2) == False:
            array1.append(float(temp1[0]))
            array2.append(float(temp2))
        # convert to numpy arrays for better handeling
        disp = np.array(array1)
        force = np.array(array2)
        
        if OutputFlag > 2:
            print('Displacements: ' + str(disp))
            print('Forces: ' + str(force))
            
    return disp,force

def distance2D(point1,point2):
    dist = 0
    for i in range(2):
        dist = dist + np.square(point1[i] - point2[i])
    dist = sqrt(dist)
    return dist

# get the associated python index
if LoadDir == 'x':
    LoadDir = 0
elif LoadDir == 'y':
    LoadDir = 1
elif LoadDir == 'z':
    LoadDir = 2

#Look for a subset of the given name and use it if it exists
DispSet = DefineSubset(NeckNodeSet,'n')
frForceSet = DefineSubset(ForceSet,'n')
StrainSet = DefineSubset(NeckElemSet,'e')

#Read in the file that contains data for all frames for time, exten, and force
#LDCHolder = ReadLDCData('LDC_ref.csv')
#ExtenCoord = []
#for i in (1,2):
#    for j in ('x','y','z'):
#        ExtenCoord.append(abs(LDCHolder['pos'+str(i)+j][0]))
#LDC_ref = ReadLDCFile_txt('LDC_ref.txt')

padNum = len(str(len(odb.steps[odb.steps.keys()[-1]].frames)))

tensorList = ['11','22','33','12','13','23'] #Need to double check this for order

firstFrame = True #Used for flow control

#Begin reading the database
try:
	axisMaskList = [True] * len(axisList)#Keeps the program from reading inputs that haven't been recorded
	elemAxisMaskList = [True] * len(tensorList)
	for frInd, frame in enumerate(odb.steps[odb.steps.keys()[-1]].frames):
		#Get coordinate values from subset
		frameTime = frame.frameValue
		
		#Code to get the force and displacement for the frame
		frameForce = 0
		if frForceSet: readList = frame.fieldOutputs['RF'].getSubset(region=frForceSet)
		else: UpdateLogfile('No force set defined, Abaqus job failed')
		for node in readList.values:
			if DataPrec == 'double':
				frameForce = frameForce + node.dataDouble[LoadDir]
			else:
				frameForce = frameForce + node.data[LoadDir]
		if frForceSet: readList = frame.fieldOutputs['U'].getSubset(region=frForceSet)
		if DataPrec == 'double': frameDisp = readList.values[0].dataDouble[LoadDir]
		else: frameDisp = readList.values[0].data[LoadDir]

		dataOut = {}
		allData = {}
		elemOut = {}
		#Read in coordinates and displacement fields
		for field in fieldList:
			if DispSet: readList = frame.fieldOutputs[field].getSubset(region=DispSet)
			else: readList = frame.fieldOutputs[field]
			for axis in axisList:
				dataOut[field+axis] = []
			for node in readList.values:
				for i in range(3):
					try:
						if axisMaskList[i]:
							if DataPrec == 'double':
								dataOut[field+axisList[i]].append(node.dataDouble[i])
							else:
								dataOut[field+axisList[i]].append(node.data[i])
							#Store displacements for the origin
						else:
							dataOut[field+axisList[i]].append('nan')
					except IndexError:
						axisMaskList[i] = False
						dataOut[field+axisList[i]].append('nan')
			#Read in values for all nodes to get extensometer values
			for axis in axisList:
				allData[field+axis] = []
			for node in frame.fieldOutputs[field].values:
				for i in range(3):
					try:
						if axisMaskList[i]:
							if DataPrec == 'double':
								allData[field+axisList[i]].append(node.dataDouble[i])
							else:
								allData[field+axisList[i]].append(node.data[i])
							#Store displacements for the origin
						else:
							allData[field+axisList[i]].append('nan')
					except IndexError:
						axisMaskList[i] = False
						allData[field+axisList[i]].append('nan')
		#End of field loop
		
		#Adjust displacements by their proper factors
		dataOut['Ux'] = [val*DispScaleX for val in dataOut['Ux']]
		dataOut['Uy'] = [val*DispScaleY for val in dataOut['Uy']]
		
		#Define the extensometer points
	#        xyzList = list(zip(allData['COORDx'],allData['COORDy'],allData['COORDz']))
	#        if firstFrame:
	#            dist1, dist2 = [], []
	#            for entry in xyzList:
	#                dist1.append(distance2D(ExtenCoord[:2],entry))
	#                dist2.append(distance2D(ExtenCoord[3:],entry))
	#            dist1Ind = dist1.index(min(dist1))
	#            dist2Ind = dist2.index(min(dist2))
	#            firstFrame = False
	#        frameDisp = allData['Uy'][dist1.index(min(dist1))]
	#        print(frameDisp)
	#        frameExten = allData['COORDy'][dist1.index(min(dist1))]
	#        print(frameExten)
	#        frameExten = frameExten + allData['COORDy'][dist2.index(min(dist2))]
	#        print(frameExten)
	#        print('')
		
		#Collect element-based field data
		for field in elemFieldList:
			if StrainSet: readList = frame.fieldOutputs[field].getSubset(region=StrainSet)
			else: readList = frame.fieldOutputs[field]
			#initialize elemOut dict
			if field == 'LE':
				for tens in tensorList:
					elemOut[field+tens] = []
			if field == 'COORD':
				for i in range(1,4):
					elemOut[field+str(i)] = []
			else:
				elemOut[field] = []
			
			for element in readList.values:
				if field == 'LE':
					for i in range(len(tensorList)):
						try:
							if elemAxisMaskList[i]:
								if DataPrec == 'double':
										elemOut[field+tensorList[i]].append(element.dataDouble[i])
								else:
										elemOut[field+tensorList[i]].append(element.data[i])
							else:
								elemOut[field+tensorList[i]].append('nan')
						except IndexError:
							elemAxisMaskList[i] = False
							elemOut[field+tensorList[i]].append('nan')
				elif field == 'COORD':
					#loop over directions
					for i in range(3):
						try:
							if elemAxisMaskList[i]:
								if DataPrec == 'double':
										elemOut[field+str(i+1)].append(element.dataDouble[i])
								else:
										elemOut[field+str(i+1)].append(element.data[i])
							else:
								elemOut[field+str(i+1)].append('nan')
						except IndexError:
							elemAxisMaskList[i] = False
							elemOut[field+str(i+1)].append('nan')	
				elif 'SDV' in field:
					try:
						if elemAxisMaskList[i]:
							if DataPrec == 'double':
								elemOut[field].append(element.dataDouble)
							else:
								elemOut[field].append(element.data)
						else:
							elemOut[field].append('nan')
					except IndexError:
						elemAxisMaskList[i] = False
						elemOut[field].append('nan')

		#Write nodal data into file, including load-displacement
		numStr = int2strPad(frInd,padNum)
		fileName = JobName + '_ABQFrame' + numStr + '_NodeData.txt'
		file = open(os.path.join(fileDir + os.path.sep + fileName),'w+')
		#Any once per frame data goes here
	#        file.write('Time\t'+str(frameTime)+'\tForce\t'+str(frameForce)+'\tExtens\t'+str(frameExten)+'\tDisp\t'+str(frameDisp)+'\n')
		# ALEX: remove extensometer data
		file.write('Time\t'+str(frameTime)+'\tForce\t'+str(frameForce)+'\tDisp\t'+str(frameDisp)+'\n')
		for name in fieldNameList:
			for axis in axisList:
				file.write(name+axis+'\t')
		file.write('\n')
		#Any once per point data goes here
		for k in range(np.size(dataOut[fieldList[0]+axisList[0]])):
			for name in fieldList:
				for axis in axisList:
					file.write(str(dataOut[name+axis][k])+'\t')
			file.write('\n')
		file.close()
		
		#Write element-wise data into file
		fileName = JobName + '_ABQFrame' + numStr + '_ElementData.txt'
		file = open(os.path.join(fileDir + os.path.sep + fileName),'w+')
		#Any once per frame data goes here
	#        file.write('Time\t'+str(frameTime)+'\tForce\t'+str(frameForce)+'\tExtens\t'+str(frameExten)+'\n')
		# ALEX: not necessary anymore
		file.write('Time\t'+str(frameTime)+'\tForce\t'+str(frameForce)+'\n')
		for i,field in enumerate(elemFieldList):
			name = elemFieldNameList[i]
			if field == 'LE':
				for tens in tensorList:
					file.write(name+tens+'\t')
			elif field == 'COORD':
				for name in axisList:
					file.write('IntPoint_'+name+'\t')
			else:
				file.write(name+'\t')
		file.write('\n')
		#Any once per point data goes here
		#loop over elements
		for k in range(np.size(elemOut[elemFieldList[0]+'1'])):
			for field in elemFieldList:
				if field == 'LE':
					for tens in tensorList:
						file.write(str(elemOut[field+tens][k])+'\t')
				elif field == 'COORD':
					for i in range(1,4):
						file.write(str(elemOut[field+str(i)][k])+'\t')
				else:
					file.write(str(elemOut[field][k])+'\t')
			file.write('\n')
		file.close()
        
except Exception as ex:
	UpdateLogfile('An error was encountered in the post-processing script, try changing the DataPrec')
	template = "An exception of type {0} occurred. Arguments: {1!r}"
	message = template.format(type(ex).__name__, ex.args)
	UpdateLogfile(message)
	print(message)
#End of frame loop
odb.close()
sys.exit()