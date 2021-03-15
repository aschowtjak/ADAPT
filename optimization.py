# LDC

# -*- coding: utf-8 -*-
from scipy.optimize import fmin
import scipy.optimize as sopt
import time
from scipy import interpolate
from scipy import special
from numpy import *
from math import sqrt 
import os as os
import sys
from collections import defaultdict
from imp import reload
from matplotlib import cm as cm
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import numpy as np
import shutil
from itertools import islice
import subprocess
from sklearn.neighbors import NearestNeighbors
from colorama import Fore, Back, Style #colored terminal outputs

global weighting
global EnableAllGraphs
global firstRun
global fLDC0, fDIC0
global CreateInterpDistPlots
global feIn

import matplotlib.path as mpltPath
from random import uniform

#==============================================================================
#  NOTES FOR DOCUMENTATION
#==============================================================================
#   -Node sets must be created on the assembly level instead of on the part
FortDir = None          # needs to be initialized for if FortDir: ...
UserMat = None          # needs to be initialized for if UserMat: ...
DataPrec = ''

#!!TEMPORARY ONLY
# Define the maximum number of iterations to be performed
IterNum = 1000

EnableAllGraphs = False #USE FOR TESTING ONLY - massively slows program down
CreateInterpDistPlots = False
    
################################# optimization ################################
#
# optimizes a set of parameters and plots the result
#
###############################################################################

# advantages of Python script:
# - no licenses needed
# - implementation of aramis data possible
# - possible to parallelize for gradient based optimization
# - smoothening of the input data (simulation or experimental data)

#==============================================================================
# TODO
#==============================================================================
#   -Create ouput functions that can be optionally used
#   -Examine efficiency of the post-processing script
#   -Certain Outputs should only be generated for the first iteration
#   -automatically detect elastic regime and neglect in objective function
#   -use python module 'pickle' to save data of optimization run

#==============================================================================
#  POSSIBLE FUNCTIONALITY
#==============================================================================
#   -Choose input and output directories besides the cwd (breaks a few path dependencies)

#==============================================================================
# INPUT PARAMETERS
#==============================================================================
JobNames = ['r15']
FrameLabelsDIC   = ['r15_t2p5_00'] #labels to identify Aramis frame data. Need to be in same order as jobnames!
LabelsDMG = ['damageField.csv']
LabelsLDC = ['LDC_ref.txt']
DICheaderLines = 4 #number of header lines in DIC-files
quadrantDIC = {'r15': 3} #specify the qudrant of the DIC-data to use.
PostProcessingScript = 'getAbqDisp'# post processing script (without file extension)
NeckNodeSet = 'NSETDISPFIELD'#Name of node set in Abaqus to be analyzed, None for all
ForceSet = 'NSETFORCE'#Name of the node set in Abaqus for the force
NeckElemSet = 'ESETSTRAIN' #Element set that should likely corrsepond to the node set in NeckNodeSet
#DataPrec = 'double'#Comment this out to read in double mode
LoadDirFE = 'y'
LoadDirDIC = 'y'
NumCPUs = 1

#UserMat = 'Lemaitre_3D'
UserMat = None # This has to be specified if no user material is used
AbqMode = '3D' #Let the program know if the Abaqus test will be 2D or 3D

weighting = {'LDC': 1.0, 'strain': 1.0}#,'displacement': 0.5} #weigthing factors
#weighting = {'LDC': 0.5,'strain': 0.5} #weigthing factors
normFlag = False #devide DIC data by frame displacement and LDC data by frame force
firstRunNorm = False #activate to normalize objective function values by vlues of first run

validationFlag = False # set to true when using FE data as input for validation pruposes
autodelInnerNoodes = True # automatically delete FE-nodes that are not on the top surface area

objectiveField = 'displacement'
#objectiveField = 'strain'


fixParString=[]
xFix=[ ]

#ParString = ['x_K_a_par','x_K_b_par','x_K_n_par','Big_S_par','Small_s_par','Beta_par','sigma_y_0_11_par']
# initial guess
#x0 = [ 1238.1,1.153e-4,1.38e-1,1.0e5,1.0,0.0]
x0 = [ 1200.0, 0.001, 0.13]
# bounds must be any positive number, need at minimum as many tuples as parameters
bnds = ((0., 999999.0),(0., 999999.0),(0., 999999.0),(0., 999999.0),(0., 999999.0),(0., 999999.0),(0., 999999.0),(0., 999999.0))


#global criterion method vector
idealVector = {'LDC': 907.302196,'strain': 1.879976772186068e-05}

# restart a aborted simulation by using the iter-file as a restart file
# does not work with evolutionary optimization algorithms
#restartFile = 'restart.txt'
restartFile = None


# symmetry data - This is for flat dogbone w/ all symm planes, should eventually have multiple cases
ForceScale          = 4
DispScale           = 1


#==============================================================================
# CREATE LOGFILE
#==============================================================================
#get run name
RunName = os.path.split(os.getcwd())[-1]
if os.path.exists(RunName + '_log.txt'):
    os.remove(RunName + '_log.txt')
log = open(RunName+'_log.txt', 'w')
log.write(str(time.strftime("%c")) + ':  # O P T I M I Z A T I O N   S T A R T E D\n')
StartTime = time.clock()
log.close()


#==============================================================================
# INPUT PARAMETERS (only for advanced users)
#==============================================================================
# directories for starting abaqus:
#AbaqusCmd = 'abq6142'       # version 6.14
AbaqusCmd = 'abq2016'       # version 2016

#disable call of abaqus (just for debugging purposes):
noAbaqus = False
oneIter = False

# Depending on the computer used, fortran needs to be called explicitly
FortDir = '\"C:\Program Files (x86)\Intel\Composer XE 2013 SP1'+os.path.sep+'bin\ifortvars_intel64.bat\"'
if FortDir:
    FortDir = FortDir + ' && '
else:
    FortDir = ''
#UserMat = 'user.for'
if UserMat:
    UserMatABQ = ' user=' + UserMat + '.for'
else:
    UserMatABQ = ''

# name of the file containing the material parameters to be identified
ParameterFileName = UserMat
ParameterFileNameExt = '.for'

# interpolation data
EquidistantSamples = False
NumPoints = 50
DynamicWeight = True

#displacement interval considered in LDC-objective-function
#used to neglect elastic regime
global x_start
global x_end
x_start = 0.00
x_end   = 3.0


RawDir = 'RawFrameData'#Place where raw data for each frame is temporarily saved
ExpFolder = 'Experimental Results'#Name of the folder the DIC data is in
    #The last iteration stays in this folder
##### NOTE TO SELF: NEED TO USE DATADOUBLE FOR EXPLICIT AND DATA FOR IMPLICIT

# set up everything for storage of temporary files
global TempDir, cwd
cwd = os.getcwd()               # current work directory
TempFolder = os.path.sep+'temp' # Folder for temporary files
TempDir    = cwd + TempFolder   # Directory to TempFolder
ResDir = cwd+os.path.sep+'Optimization Results' # Directory to optimization results
GraphDir = ResDir+os.path.sep+'FrameGraphs'
LdcDir = ResDir+os.path.sep+'LDCplots'

# delete old optimization results
if os.path.exists(ResDir):
    shutil.rmtree(ResDir)
# delete old temp folder
if not noAbaqus:
    if os.path.exists(TempDir):
        shutil.rmtree(TempDir)
# delete old iter file
if os.path.exists(RunName + '_iter.txt'):
    os.remove(RunName + '_iter.txt')

DispScaleX = DispScale
DispScaleY = DispScale

global DICPropList#List of the properties as a string expected from the DIC input
if validationFlag:
    if objectiveField == 'damage':
        DICPropList = ['x','y','z','eps11','eps22','eps33','eps12','eps13','eps23','damage']
    else:
        DICPropList = ['x','y','z','ux','uy','uz']
elif LoadDirFE == LoadDirDIC:
    if objectiveField == 'displacement' or objectiveField == 'strain':
        DICPropList = ['ID','x','y','z','ux','uy','uz','eps11','eps22']
    elif objectiveField == 'damage':
        DICPropList = ['ID','x','y','z','damage']
else:
    if objectiveField == 'displacement' or objectiveField == 'strain':
        DICPropList = ['ID','y','x','z','ux','uy','uz','eps11','eps22']
    elif objectiveField == 'damage':
        DICPropList = ['ID','y','x','z','damage']
        

# output generation (user and output options)
global OutputFlag
global LogOutputFlag
global SaveInitGuessData
global SaveLastIterationData

OutputFlag = 1 #0 is for minimal output, 1 gives developer output, 2 enables command lines
LogOutputFlag = 1 #Same as above but no option for 2
SaveInitGuessData = 0
SaveLastIterationData = 0

# perturbation parameters for gradient based optimization methods
eps = 1e-6
epsilon = np.asarray(x0)*eps

# helper to determine first run and dicts to save obj. fun. value of 1st iteration
firstRun = True
fDIC0, fLDC0 = {}, {}

terminalWidth = shutil.get_terminal_size().columns
print('\n\n')
message = 'S T A R T I N G   O P T I M I Z A T I O N   R O U T I N E\n'
print(message.center(terminalWidth))




#==============================================================================
# FUNCTION DEFINITIONS
#==============================================================================
def AdjustData(Data1,Data2,EquidistantSamples=False,NumPoints=100):
    # Rearrange Data1 and Data2, such that both have the same format
    # For the evaluation of the objective function, the Load-Displacement-Curve
    # needs to be evaluated at the same sample points. Therefore the data
    # with higher resolution (= more sample points) is interpolated on the 
    # sample points with lower resolution.
    # If EquidistantMesh is True, [NumPoints] data samples with equivalent
    # distance are generated. If non-equidistant sample points are used,
    # it would mean that non-uniform weights are being used.
    
    # TODO: daten anpassen: wenn ein array eine größere verschiebung aufweist,
    #       dann muss das andere array entsprechend angepasst werden
    
    # INFO in case of errors: Scale factors or the adaption
    # of the PostProcessingScript might be incorrect

    # rearrange data
    x1 = abs(Data1[0])  # Simulation
    y1 = abs(Data1[1])  # Simulation
    x2 = abs(Data2[0])  # Experiment
    y2 = abs(Data2[1])  # Experiment

#    import shit

#    # ALEX
#    print('\nAdjustData After Fine to Coarse:')
#    print('Length of Data1 (LDC): ' + str(len(x1)))
#    print('Length of Data2 (LDC_ref): ' + str(len(x2)))
#    print('Maximum x Data1: ' + str(x1[-1]))
#    print('Maximum x Data2: ' + str(x2[-1]))
#    print(x1[-1]-x2[-1])
#    print('')
    
    # if simulation data shows higher maximum displacement than experimental data
    if x1[-1] > x2[-1]:
        if len(fvec) == 0:
            message = 'The simulation shows larger displacements than the experiment. Hence, the displacements exceeding the experimental values are ignored.'
            sendWarning(message,True)
        searching = True
        for eCt,entry in enumerate(x1):
            #find first entry of sim data that is higher than max exp displacement
            if entry > x2[-1] and searching:
                # if simulation is finer (interpolation from simulation to experiment)
                if len(x1[:eCt]) > len(x2[:eCt]): 
                    x1, y1 = x1[:eCt+1], y1[:eCt+1]   # TODO: this needs to be tested
                else: 
                    x1, y1 = x1[:eCt], y1[:eCt]
                searching = False
    
#        print('\nAdjustData After Large Simulation Displacements Filtered:')
#        print('Length of Data1 (LDC): ' + str(len(x1)))
#        print('Length of Data2 (LDC_ref): ' + str(len(x2)))
#        print('Maximum x Data1: ' + str(x1[-1]))
#        print('Maximum x Data2: ' + str(x2[-1]))
#        print(x1[-1]-x2[-1])
#        print('')
    
    # if experimental data shows higher maximum displacement than simulation data
    elif x2[-1] > x1[-1]:
        
        if len(fvec) == 0:
            message = 'The simulation shows lower displacements than the experiment. Hence, the experimental displacements exceeding the simulation values are ignored.'
            sendWarning(message,True)
        searching = True
        for eCt,entry in enumerate(x2):
            if entry > x1[-1] and searching:
                # if simulation is finer (interpolation from simulation to experiment)
                if len(x2[:eCt]) > len(x1[:eCt]): 
                    x2, y2 = x2[:eCt+1], y2[:eCt+1]   # TODO: this needs to be tested
                else: 
                    x2, y2 = x2[:eCt], y2[:eCt]
                searching = False
                
#        print('\nAdjustData After Large Simulation Displacements Filtered:')
#        print('Length of Data1 (LDC): ' + str(len(x1)))
#        print('Length of Data2 (LDC_ref): ' + str(len(x2)))
#        print('Maximum x Data1: ' + str(x1[-1]))
#        print('Maximum x Data2: ' + str(x2[-1]))
#        print(x1[-1]-x2[-1])
#        print('')   

    
    # define interpolation function
    try:
        f1 = interpolate.interp1d(x1,y1,fill_value='extrapolate')
        f2 = interpolate.interp1d(x2,y2,fill_value='extrapolate')
    except ValueError:
        plt.figure(figsize=(8,0.5))
        plt.title('Sample points plot')
        plt.vlines(Data1[0],0,0.66,'red',alpha=0.6,linewidth=1,label='Simulation sample points(Data1)')
        plt.vlines(Data2[0],0.33,1,'blue',alpha=0.6,linewidth=1,label='Experimental sample points(Data2)')
        plt.gca().axes.get_yaxis().set_visible(False)
        plt.xlabel('displacement in mm')
        plt.legend(loc=9, bbox_to_anchor=(0.5, -1), ncol=2)
        
        message = 'ERROR: LDC data interpolation failed! Check sample points plot. Each line represents a sample point. Sample points should overlap in the region of interest.'        
        sendError(message)
    
    if EquidistantSamples == True:
        # create equidistant x vector
        x1 = np.linspace(x1[0],x1[-1],num=NumPoints,endpoint=True)
        x2 = x1
    else:
        # interpolate from fine to coarse resolution
        if len(x1) > len(x2):
            x1 = x2
#            print('interpolating from simulation to experimental data')
        elif len(x2) > len(x1):
            x2 = x1
#            print('interpolating from experimental to simulation data')
        else:       # if both have the same length, but different sample points
            x2 = x1
    
#    print('\nAdjustData After Fine to Coarse:')
#    print('Length of Data1 (LDC): ' + str(len(x1)))
#    print('Length of Data2 (LDC_ref): ' + str(len(x2)))
#    print('Maximum x Data1: ' + str(x1[-1]))
#    print('Maximum x Data2: ' + str(x2[-1]))
#    print('')
    
    # interpolate data
    y1 = f1(x1)
    y2 = f2(x2)
    
#    print('\nAdjustData Final:')
#    print('Maximum y Data1: ' + str(max(y1)))
#    print('Maximum y Data2: ' + str(max(y2)))
#    print('Maximum x Data1: ' + str(x1[-1]))
#    print('Maximum x Data2: ' + str(x2[-1]))
#    print('')
#    xfrac = ReadElementDeletion('ElementDeletion.txt')
#    xfrac = 2*xfrac     # FIXME: introduce SymDispMultiplier 
#    for i,e in enumerate(y1):
#        if x1[i] > xfrac:
#            y1[i:] = np.zeros(len(y1)-i)
#            break
    
#    SaveLDC_InitialGuess([x1,y1],[x2,y2])

    SaveLDC_Optimum([x1,y1])
    return [x1,y1],[x2,y2]


def call_abaqus_single(job):
    # this version of call_abaqus is more robust but works only for one single job
    # job cannot be termiated upon user interrupt of PI-script
    
    os.chdir(job)
    
    try:
       os.remove('abaqus.rpy')
    except OSError:
       pass
    # all file extensions to be deleted
    DelString = ['.dat','.com','.odb','.msg','.prt','.sta','.log','.sim',
                 '.lck','.fil','.par','.pes','.pmg','.tec','.txt']
    
    UpdateLogfile('Running Abaqus job ...')
    print('Running Abaqus job ...')

    # delete old job-files
    for ext in DelString:
        try:
           os.remove(job + ext)
        except OSError:
           pass
        
    startupinfo = subprocess.STARTUPINFO()
    if OutputFlag <= 2:
        startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
    ErrorCode = subprocess.call('run_'+job+'.bat', startupinfo=startupinfo, shell=False)


    if ErrorCode == 0:
        error_bool = False
        if LogOutputFlag > 0: UpdateLogfile('Abaqus job has finished. Call post processing script ...')
        if OutputFlag > 0: print('Abaqus job has finished. Call post processing script ...')
        
        # call post processing script
        startupinfo = subprocess.STARTUPINFO()
        if OutputFlag <= 2:
            startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
        command = AbaqusCmd + ' cae noGUI=' + PostProcessingScript
        PPerror = subprocess.call(command, startupinfo=startupinfo, shell=True)
        if PPerror != 0:
            error_bool = True
            if OutputFlag > 0: print('Potential error in the post processing script')
            if LogOutputFlag > 0: UpdateLogfile('Potential error in the post processing script')
        elif PPerror == 0:
            if OutputFlag > 0: print('Post processing script has finished')
            if LogOutputFlag > 0: UpdateLogfile('Post processing script has finished')
        
        
    else:
        error_bool = True
        print('WARNING: Abaqus job did not finish properly.\n')
        UpdateLogfile('WARNING: Abaqus job did not finish properly')
        

    # change back to temp directory
    os.chdir('..')        
    return error_bool
        
def call_abaqus_multi(JobList):
    # this version of call_abaqus is less robust but works for multiple jobs
    # call all jobs and post processing scripts
    
    #command to call post processing script
    command = AbaqusCmd + ' cae noGUI=' + PostProcessingScript
    
    for job in JobList:
        # change to directory of current job
        os.chdir(job)
        try:
           os.remove('abaqus.rpy')
        except OSError:
           pass
        # all file extensions to be deleted
        DelString = ['.dat','.com','.odb','.msg','.prt','.sta','.log','.sim',
                     '.lck','.fil','.par','.pes','.pmg','.tec','.txt','.abq',
                     '.mdl','.pac','.res','.sel','.stt']
    
        # delete old job-files
        for ext in DelString:
            try:
               os.remove(job + ext)
            except OSError:
               pass

        UpdateLogfile('Call Abaqus job '+job+' ...')
        print('Call Abaqus job '+job+' ...')
        # call abaqus
        startupinfo = subprocess.STARTUPINFO()
        if OutputFlag <= 2:
            startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
        subprocess.call('run_'+job+'.bat', startupinfo=startupinfo, shell=False)
        # change back to temp directory
        os.chdir('..')

    
    # monitor all jobs
    
    # create bool array to save job status
    jobCompleted  = np.full(len(JobList),False, dtype=bool)

    # set boolean to false
    error_bool = False
    PPerror = 0
    error_string = 'Abaqus/Analysis exited with errors' # string if sim. failed
    time.sleep(5) # initial pause
    while (not all(jobCompleted) and error_bool == False and PPerror == 0):
        pp_called = False # marker to check if post processing script was called
        # loop over jobs
        for ind,job in enumerate(JobList):
            # check only jobs that are not completed yet
            if not jobCompleted[ind]:
                compl_string = 'Abaqus JOB '+ job +' COMPLETED'      # string if sim. completed
                os.chdir(job)
                with open(job+'.log','r') as fidc:  # check log-file for strings
                    for i,line in enumerate(fidc):
                        if line.startswith(compl_string):
                            jobCompleted[ind] = True
                            if LogOutputFlag > 0: UpdateLogfile('Abaqus job ' + job + ' has finished. Call post processing script ...')
                            if OutputFlag > 0: print('Abaqus job ' + job + ' has finished. Call post processing script ...')
                            
                            # call post processing script
                            startupinfo = subprocess.STARTUPINFO()
                            if OutputFlag <= 2:
                                startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
                            PPerror = subprocess.call(command, startupinfo=startupinfo, shell=True)
                            if PPerror != 0:
                                if OutputFlag > 0: print('Potential error in the post processing script')
                                if LogOutputFlag > 0: UpdateLogfile('Potential error in the post processing script')
                            elif PPerror == 0:
                                if OutputFlag > 0: print('Post processing script for job '+job+' has finished')
                                if LogOutputFlag > 0: UpdateLogfile('Post processing script for job '+job+' has finished')
                            pp_called = True
                            
                        elif line.startswith(error_string):
                            error_bool = True
                            message = 'Abaqus job ' + job + ' did not finish properly.'
                            sendWarning(message)
                            # kill all remaining jobs
                            for job2kill in JobList:
                                if not job2kill == job:
                                    os.chdir('..')
                                    os.chdir(job2kill)
                                    if LogOutputFlag > 0: UpdateLogfile('Terminate '+job2kill)
                                    if OutputFlag > 0: print('Terminate job '+job2kill)
                                    startupinfo = subprocess.STARTUPINFO()
                                    if OutputFlag <= 2:
                                        startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
                                    subprocess.call('kill_'+job2kill+'.bat', startupinfo=startupinfo, shell=False)
                fidc.close()
                os.chdir('..')

        # pause for each time, only if pp script was not called
        if not pp_called:
            time.sleep(5)
    
    if error_bool == False and PPerror == 0:
        print('All jobs have finished successfully')
        UpdateLogfile('All jobs have finished successfully')
    return error_bool

def CopyFilesToTempDir(FileList,labelDIC,labelLDC,labelDMG,destination):
    # copy all files to temporary directory
    for e in FileList:
        try:
            shutil.copy(e,os.path.join(TempDir,destination))
        except OSError:
            pass
    for roots, dirs, files in os.walk(ExpFolder):
        padLen = len(str(len(files)))
        fcount = 0
        for frame in files:
            # copy DIC frame files
            if labelDIC in frame:
                numStr = int2strPad(fcount,padLen)
                fcount = fcount+1
                RawFile = destination + '_DICFrame' + numStr + '.txt'
                shutil.copy(os.path.join(ExpFolder,frame),os.path.join(TempDir,destination,RawDir,RawFile))
            # copy LDC and DMG file
            elif labelLDC in frame or labelDMG in frame:
                shutil.copy(os.path.join(ExpFolder,frame),os.path.join(TempDir,destination,RawDir,frame))




def distance2D(point1,point2):
    dist = 0
    for i in range(2):
        dist = dist + np.square(point1[i] - point2[i])
    dist = sqrt(dist)
    return dist

def ExtractLDCData(LDCHolder):
    #Adjust the force and displacement measurments
    if LDCHolder['force'][0] != 0:
        message = 'Non-zero initial force of value %.3f has been shifted to 0' %LDCHolder['force'][0]
        sendWarning(message)
    LDCHolder['force'] = [x - LDCHolder['force'][0] for x in LDCHolder['force']]
    dispYTemp = [x - LDCHolder['extenY'][0] for x in LDCHolder['extenY']]
    LDC_ref = (np.array(dispYTemp),np.array(LDCHolder['force']))
    LDC_ref = RemoveCorruptLDCData(LDC_ref)
    
    return LDC_ref


def findSetPoints(setName,inputfile):
    reading = False
    tag = None
    pointList = []
    with open(inputfile,'r') as file:
        for line in file:
            #find beginning of set keyword
            if setName in line.upper() and 'PART-1-1' in line.upper():
                reading = True
                tag = line[-9:]
            #fing end of set keyword
            elif reading and '*' in line:
                reading = False
            #read lines of set keyword
            elif reading:
                temp = line.split(',')
                temp = [x.strip() for x in temp]
                for x in temp:
                    if not x == '':pointList.append(int(x))
    #check for generate-option that specifies a set by a range
    if tag == 'generate\n':
        pointList = list(range(pointList[0],pointList[1]+1,pointList[2]))
        
    if len(pointList) == 0:
        message = 'No nodes/elements found in the set ' + setName +'! Check set names in line 80 to 82.'
        sendWarning(message)
    else:
        message = str(len(pointList)) + ' nodes/elements found in the set ' + setName + ' in input file ' + os.path.split(inputfile)[1]
        if OutputFlag > 0: print(message)
        if LogOutputFlag > 0: UpdateLogfile(message)
    return pointList



def getDirs(JobNames):
    # get names of directories containing abaqus input files
    dirNames = {}
    for root, dirs, files in os.walk("."):
        for name in files:
            for job in JobNames:
                if job in name and not 'temp' in root and not 'Experimental Results' in root :
                    dirNames[job] = (os.path.split(root)[-1])
    return dirNames


def getObjectiveFunctionValueDIC(DICdata,FEdata,normFlag=False):
    # Calculate objective function of the form
    
    if normFlag:
        Weight = np.copy(DICdata)
        for i in range(len(Weight)):
            if abs(DICdata[i]) > 0.001:
                Weight[i] = abs(1/Weight[i])
            else:
                Weight[i] = 0
    else:
        Weight=np.ones(len(DICdata))
    
    #Check assumptions about length
    tempSum = [abs(x) - abs(y) for x,y in zip(DICdata,FEdata)]
    
    #number of data points that are used in objective function / dont have zero weight
    #exeption frame 1: weight of all datapoints is zero --> set ndp to 1 to avoid devision by zero
    ndp = 1
    if np.count_nonzero(Weight)>0: ndp = np.count_nonzero(Weight)
    
    f = sum((tempSum*Weight)**2)/ndp
#    f = np.linalg.norm(np.array(tempSum)*Weight)/len(DICdata)
    
    return f



def getObjectiveFunctionValueLDC(Data1,Data2,DynamicWeight=True,Weight=None):
    # Calculate objective function of the form
    # f = sqrt( sum(Data1[i]-Data2[i]) ) / len(Data1)
    # The "Weight" argument is optional. If its not specified, all weights
    # are set to ones. It has to be of the same length as Data1 and Data2.
    global x_start
    global x_end
    
    # if no weight is specified, set to ones
    if Weight==None:
        Weight=np.ones(len(Data1[0]))
        
    ErrFlag = 0
    
    if DynamicWeight == True:
        DWeight = []
        for i in range(len(Data1[0])):
            if np.isnan(Data1[0][i]) == True or np.isnan(Data1[1][i]) == True \
            or np.isnan(Data2[0][i]) == True or np.isnan(Data1[0][i]) == True:
                message = 'Corrupt simulation results! (LDC) Some data is nan.'
                sendWarning(message)
                ErrFlag = 1
                break
            
            # for first element only
            if i == 0:
                w = Data1[0][i+1]/2
            # for last element only
            elif i == len(Data1[0])-1:
                w = (Data1[0][i]-Data1[0][i-1])/2
            # for all elements between the first and last ones
            else:
                w = (Data1[0][i+1]-Data1[0][i-1])/2
            
            DWeight.append(w)

        # convert to numpy array
        DWeight = np.array(DWeight/np.linalg.norm(DWeight))
    else: DWeight=np.ones(len(Data1[0]))
    
    # set weights to zero for all points with greater displacement than xcrit
    for i in range(len(Data1[0])):
        if Data1[0][i] < x_start or Data1[0][i] > x_end:
            Weight[i] = 0
    
    # error catcher, if weight is not given correctly
    if len(Weight) != len(Data1[0]):
        Weight=np.ones(len(Data1[0]))
        print('Error! Dimension mismatch in getObjectiveFunction(...)')
        print('Continuing computation with Weight = ones ...')

    if ErrFlag == 0:
#        plt.figure()
#        plotvec=((abs(Data1[1][:len(Weight)])-abs(Data2[1][:len(Weight)]))*Weight*DWeight)**2
#        plt.plot(plotvec[1:],'x-')
#        plt.title(job)
#        plt.xlabel('frame no.')
#        plt.ylabel('objective function value LDC')
#        plt.grid()

#        f = np.linalg.norm((abs(Data1[1][:len(Weight)])-abs(Data2[1][:len(Weight)]))*Weight*DWeight)/np.count_nonzero(Weight)
        f = np.sum(((abs(Data1[1][:len(Weight)])-abs(Data2[1][:len(Weight)]))*Weight*DWeight)**2)/np.count_nonzero(Weight)
        
    else:
        f = 100000000

#    xfrac = ReadElementDeletion('ElementDeletion.txt')
#    xfrac = 2*xfrac     # FIXME: introduce SymDispMultiplier 
#    global xfrac_ref
#    f = abs(xfrac_ref-xfrac)
    return f



def int2strPad(integer,padLength=3):
    outStr = ''
    if integer >= 10**padLength: print('Padding length may be insufficient to ensure ordering')
    if integer == 0: outStr = str('0'*padLength)
    else:
        for i in reversed(range(padLength+1)):
            if integer < 10**i and integer >= 10**(i-1):
               outStr = str('0'*(padLength-i) + str(integer))
    return outStr


#inter_err_bool,eps11Out2,eps22Out2,distance2 = interpolate2DSpace(feDat,dicDat2,'strain')
def interpolate2DSpace (mainList,interList,objectiveField): #Normally FE is main and DIC inter
    NumNeigh = 7            # number of neighbors
    Radius = 10              # radius for the nearest neighbor search
    neigh = NearestNeighbors(n_neighbors=NumNeigh,radius=Radius)
    # extract coordinates from input
    inter_data = [i[0:2] for i in interList]
    neigh.fit(inter_data)
    main_data = [i[0:2] for i in mainList]
    neighbors = neigh.kneighbors(main_data, NumNeigh, return_distance=False)
#    print("--- %s seconds ---" % (time.time() - start_time))
    if objectiveField == 'displacement' or objectiveField == 'strain':
        d1List, d2List, d3List, sumDist = [], [], [], []
    elif objectiveField == 'damage':
        dList, sumDist = [], []
    polygon = [([],[]),([],[]),([],[])]
    x1,y1 = [0,0,0],[0,0,0]
    errorFlag = False
    

    # assemble comninations for choosing neighbors
    m=0
    combinations = np.zeros((int(special.binom(NumNeigh-1,2)),3))
    for k in range(1,NumNeigh):
        for l in range(k+1,NumNeigh):
            combinations[m] = [0,k,l]
            m = m+1
    combinations = combinations.astype(int)
    
    #Loop through the list to be interpolated to
    for i,item in enumerate(mainList):
        if objectiveField == 'displacement':
            xneigh,yneigh,Dux,Duy,Duz = [],[],[],[],[]
        elif objectiveField == 'strain':
            xneigh,yneigh,Deps11,Deps22 = [],[],[],[]
        elif objectiveField == 'damage':
            xneigh,yneigh,DD = [],[],[]
        #Determine where the NumNeigh closest points in each list are
        interInd = neighbors[i]
        #Define the variables needed for the interpolation equation
        x,y = item[0],item[1]
        
        if validationFlag:
            #nearest neighbor interpolation
            if objectiveField == 'displacement':
                projCalcx = interList[interInd[0]][-3]
                projCalcy = interList[interInd[0]][-2]
                projCalcz = interList[interInd[0]][-1]
            elif objectiveField == 'strain':
                projCalcx = interList[interInd[0]][-2]
                projCalcy = interList[interInd[0]][-1]
            elif objectiveField == 'damage':
                projCalcD = interList[interInd[0]][-1]
            distHold = distance2D((interList[interInd[0]][0],interList[interInd[0]][1]),(x,y))
        else:
            # get coordinates of NumNeigh neighbors
            for j in range(NumNeigh):
                xneigh.append(interList[interInd[j]][0])
                yneigh.append(interList[interInd[j]][1])
                
            # make sure query point is inside a triangle (=> no extrapolation)
            # this also guarantees that the three support points are not on one line
            activeNeigh = combinations[0]
            isInside = False
            itr = 0
            while not isInside:
                # get polygon of the 3 active neighbors
                for k,n in enumerate(activeNeigh):
                    polygon[k] = (xneigh[n],yneigh[n])
                if not RayTracingMethod(x,y,polygon):
#                    if itr > 5:
#                        xneighplot = [0,0,0]
#                        yneighplot = [0,0,0]
#                        colorarray = ['blue','green','yellow','red','black','purple','orange']
#                        plt.figure()
#                        for l in range(itr+1):
#                            activeNeigh = [0,1,2+l]
#                            plt.axis('equal')
#                            for k,n in enumerate(activeNeigh):
#                                xneighplot[k] = xneigh[n]
#                                yneighplot[k] = yneigh[n]
#                            xneighplot.append(xneighplot[0])
#                            yneighplot.append(yneighplot[0])
#                            plt.plot(xneighplot,yneighplot,c=colorarray[l],linewidth=1)
#                            plt.scatter(xneigh,yneigh,c='blue')
#                            plt.annotate(l,(xneigh[l],yneigh[l]))
#                            plt.scatter(x,y,c='red',marker="x")
                    # choose new neighbors
                    itr = itr + 1
                    if itr > len(combinations)-1:
                        errorFlag = True
                        break
                    activeNeigh = combinations[itr]
                else:
                    isInside = True
                    
            #save the total distance of the chosen neighbors
            #and get displacements and coordinates of chosen neighbors
            distHold = 0.0
            for j,n in enumerate(activeNeigh):
                distHold = distHold + distance2D((xneigh[n],yneigh[n]),(x,y))  
                x1[j] = xneigh[n]
                y1[j] = yneigh[n]
                if objectiveField == 'displacement':
                    Dux.append(interList[interInd[n]][-3])
                    Duy.append(interList[interInd[n]][-2])
                    Duz.append(interList[interInd[n]][-1])
                elif objectiveField == 'strain':
                    Deps11.append(interList[interInd[n]][-2])
                    Deps22.append(interList[interInd[n]][-1])
                elif objectiveField == 'damage':
                    DD.append(interList[interInd[n]][-1])    
            
            #Calcualte u and v for both lists
            uTop = (y1[2]-y1[0])*(x-x1[0]) - (x1[2]-x1[0])*(y-y1[0])
            vTop = (y1[0]-y1[1])*(x-x1[0]) + (x1[1]-x1[0])*(y-y1[0])
            bot = (y1[2]-y1[0])*(x1[1]-x1[0]) - (x1[2]-x1[0])*(y1[1]-y1[0])
    
            #Added to prevent random zero occurences
    #        if bot < 1e-10:
            if bot == 0:
                bot = bot + 1e-10
                message = 'Denominator of interpolation equation evaluated to zero'
                sendWarning(message)
            u = uTop/bot
            v = vTop/bot
            
            #Apply U and V to interpolate one point to the other
            if objectiveField == 'displacement':
                projCalcx = (1-u-v)*Dux[0] + u*Dux[1] + v*Dux[2]
                projCalcy = (1-u-v)*Duy[0] + u*Duy[1] + v*Duy[2]
                projCalcz = (1-u-v)*Duz[0] + u*Duz[1] + v*Duz[2]
            elif objectiveField == 'strain':
                projCalcx = (1-u-v)*Deps11[0] + u*Deps11[1] + v*Deps11[2]
                projCalcy = (1-u-v)*Deps22[0] + u*Deps22[1] + v*Deps22[2]
            elif objectiveField == 'damage':
                projCalcD = (1-u-v)*DD[0] + u*DD[1] + v*DD[2]

        #Append values of the current point to list
        sumDist.append(distHold)
        if objectiveField == 'displacement':
            d1List.append(projCalcx)
            d2List.append(projCalcy)
            d3List.append(projCalcz)
        elif objectiveField == 'strain':
            d1List.append(projCalcx)
            d2List.append(projCalcy)
        elif objectiveField == 'damage':
            dList.append(projCalcD)

    if objectiveField == 'displacement':
        return errorFlag, d1List, d2List, d3List, sumDist
    elif objectiveField == 'strain':
        return errorFlag, d1List, d2List, sumDist
    elif objectiveField == 'damage':
        return errorFlag, dList, sumDist
    
def kill_all_abaqus(JobList):
    for job in JobList:
        os.chdir(TempDir + os.path.sep + job)
        UpdateLogfile('Send terminate to Abaqus job '+job+' ...')
        print('Send terminate to Abaqus job '+job+' ...')
        
        startupinfo = subprocess.STARTUPINFO()
        if OutputFlag <= 2:
            startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
        subprocess.call('kill_'+job+'.bat', startupinfo=startupinfo, shell=False)


def line2listDIC(line):
    out = []
    line = line.strip('\n')
    if validationFlag:
        change = line.split('\t')
    else:
        change = line.split(',')
    for item in change:
        out.append(item.strip())
    return out



def line2listFE(line):
    line = line.strip('\n')
    change = line.split('\t')
    try:
        change.remove('')
    except ValueError:
        pass
    return change



def multiMinIndex(inList, numInd): #Returns a specified number of index minimums
    temporaryList=inList[:]#Avoid changing the input list
    minList = []
    for i in range(numInd):
        tempMin = temporaryList.index(min(temporaryList))
        minList.append(tempMin)
        temporaryList[tempMin] = 99999999 #Arbitrarily large number, tried with nan but caused errors
    return minList



def optimfunc(x,*bounds):
    global fvec
    global DICPropList
#    global LDCHolder
    global LDC_ref
    global weighting
    global fLDC0
    global fDIC0
    global firstRun
    global CreateInterpDistPlots
    global feIn
    global dicIn
    global Weight
    # initialize function value
    f = 0
    
    # dict, containing the overall objective function value for each job
    f_dict= {}
    # dict, containing the objective function value for each QOI and each job
    fQOI = defaultdict(dict)

    #Initialize storage for frame data
    feIn,dicIn,feInElem = {},{},{}
    for job in JobNames:
        feIn[job],dicIn[job],feInElem[job] = [],[],[]
    
    #initialize error boolean
    err_bool = False
    
    #False to speed up, true to get more data about frame
    if EnableAllGraphs: CreateInterpDistPlots = True

        
    terminalWidth = shutil.get_terminal_size().columns
    print('\n')
    message = 'O B J E C T I V E   F U N C T I O N   C A L L  ' + str(len(fvec)+1)
    UpdateLogfile(message)
    print(message.center(terminalWidth))
    print(RunName.center(terminalWidth))
    print('\n')
    print('current parameter set: ')
    print('[',end='')
    for i in x:
        print('% 10.8E' %i +' ',end='')
    print(']')
    
    
    if 'fListRestart' in globals() and len(fListRestart)>0:
        f = fListRestart.pop(0)
        message = 'Using data from restart file'
        if OutputFlag > 0: print(message)
        if LogOutputFlag > 0: UpdateLogfile(message)

        if firstRun:
            fQOI = fQOI_ini
            firstRun = False
        else:
            for qoi in weighting.keys():
                for jo in JobNames:
                    fQOI[qoi][job] = 1
        if OutputFlag > 0:
            print('Global objective function value = ' + str(f) + ' (from restart file)\n')
    else:
        
        SetUpFlowCurve(x)
#        sigma_y0 = x[0] * x[1] ** x[2]
#        y = ([x[0],x[1],x[2],x[3],x[4],x[5],sigma_y0])
    
    #    sigma_y0 = x[0] * xFix[0] ** xFix[1]
    #    yFix = xFix.copy()
    #    yFix.append(sigma_y0)
#        SetUpParameters(y,ParString,ParameterFileName,ParameterFileNameExt,xFix,fixParString)
        
        # penalty, if bounds are violated
        for i,e in enumerate(x):
            if e < bnds[i][0] or e > bnds[i][1]:
                sendWarning('Bounds violated! Objective function is penalized and next iteration started...')
                # do not evaluate abaqus model if bounds are violated
                err_bool = True
        
        # start simulation
        if not err_bool:
            if noAbaqus:
                err_bool = False
                if len(fvec) == 0:
                    message = 'Abaqus is not called! This is a debugging run. Switch noAbaqus to false for normal run.'
                    sendWarning(message,True)
            else:
                if len(JobNames)>1:
                    err_bool = call_abaqus_multi(JobNames)
                elif len(JobNames)==1:
                    err_bool = call_abaqus_single(JobNames[0])
#                    err_bool = False
                    
            
        # computation of the error function only if all simulations were finished properly
        # and bounds are not violated
        if (err_bool == False):
            for job in JobNames:
                #go to directory of current job
                os.chdir(job)
    
    #           --SUBSECTION DATA SETUP--
                #Read data from all of the files according to its type
                #loop over frames/files
                for root, dirs, files in os.walk(RawDir):
                    for specfile in files:
                        if 'ABQ' in specfile and 'Node' in specfile:
                            feIn[job].append(ReadFEData(os.path.join(RawDir,specfile)))
                        elif 'ABQ' in specfile and 'Elem' in specfile:
                            feInElem[job].append(ReadFEDataElem(os.path.join(RawDir,specfile)))
                        elif 'DIC' in specfile:
                            dicIn[job].append(ReadDICData(os.path.join(RawDir,specfile),DICPropList,DICheaderLines))
                #merge nodal and element fe data to one struct
                for frameInd,frame in enumerate(feIn[job]):
                    frame.update(feInElem[job][frameInd])
                #Remove the empty lists if filtering out
    #            feIn[job][:] = [li for li in feIn if li != None]
    #            dicIn[job][:] = [li for li in dicIn if li != None]

                        
                if len(feIn[job])!=len(dicIn[job]) and len(fvec) == 0:
                    message = 'FE and DIC data of job '+job+' have different numbers of frames.'
                    sendWarning(message,True)
                    
                # delete inner nodes
                if autodelInnerNoodes :
                    # initial thickness
                    t = max(feIn[job][0]['z'])
                    # list of nodes to delete
                    delList = []
                    for ind, z in enumerate(feIn[job][0]['z']):
                        if abs(z-t) > 1e-8:
                            delList.append(ind)
                    
                    for frameind, frame in enumerate(feIn[job]):
                        for ind in sorted(delList, reverse=True):
                            for key in ['x','y','z','ux','uy','uz']:
                                frame[key].pop(ind)
        
                

                if ('displacement' in weighting and weighting['displacement'] > 0) or ('strain' in weighting and weighting['strain'] > 0):
                    #note:following 3 steps could be done in same loop
                    #Add extensometer data to the dicIn values
                    for frCount,frame in enumerate(dicIn[job]):
                        frame['frameExten'] = LDC_ref[job][0][frCount]
     
                    #adjust sign of displacements
                    for frInd, frame in enumerate(dicIn[job]):
                        for i, elem in enumerate(frame['uy']):
                            dicIn[job][frInd]['uy'][i] = abs(elem)
                            dicIn[job][frInd]['ux'][i] = -abs(frame['ux'][i])
                
                    #choose the specified quadrant
                    if quadrantDIC[job] > 1:
                        for frInd, frame in enumerate(dicIn[job]):
                            for i, elem in enumerate(frame['x']):
                                if quadrantDIC[job] == 2 or quadrantDIC[job] == 3:
                                    #mirror x-axis
                                    dicIn[job][frInd]['x'][i] = -elem
                                if quadrantDIC[job] == 3 or quadrantDIC[job] == 4:
                                    #mirror y-axis
                                    dicIn[job][frInd]['y'][i] = -frame['y'][i]
    
    
        #        --SUBSECTION FRAME LOOP--
                #Work on each fe frame as the basis, assumes FE coarser than DIC
                print('Beginning data interpolation for job '+job)
                fFrame, LDC1, LDC2 = [], [], []
                numInterp = 0
                totNumInterp = len(feIn[job])*2
                for frCount, frame in enumerate(feIn[job]):
                    #Read LDC stats for this frame
                    LDC1.append(float(frame['frameDisp'])*DispScaleY)
                    LDC2.append(float(frame['frameForce'])*ForceScale)
                    
                    if ('displacement' in weighting and weighting['displacement'] > 0) or ('strain' in weighting and weighting['strain'] > 0):
                        # get the 2 LDC entries closest to current FE-frame
                        tempList = []
                        for dicFrame in dicIn[job]:
                            tempList.append(abs(dicFrame['frameExten'] - frame['frameDisp']))
                        interpFrame = multiMinIndex(tempList,2)
                        
                        # S P A T I A L   I N T E R P O L A T I O N
                        if 'displacement' in weighting:
                            propListFe  = ['x','y','ux','uy','uz']
                            propListDic = ['x','y','ux','uy','uz']
                            
                        elif 'strain' in weighting:
#                            propListFe  = ['IntPoint_x','IntPoint_y','IntPoint_z','strain11','strain22']
#                            propListDic = ['x','y','z','eps11','eps22']
                            propListFe  = ['IntPoint_x','IntPoint_y','strain11','strain22']
                            propListDic = ['x','y','eps11','eps22']
                            
                            
                        feDat = list(zip(*[frame[prop] for prop in propListFe]))  #the * unpacks the zip-arguemnts from list

                        tmpFrame = dicIn[job][interpFrame[0]]
                        dicDat1 = list(zip(*[tmpFrame[prop] for prop in propListDic]))
                        tmpFrame = dicIn[job][interpFrame[1]]
                        dicDat2 = list(zip(*[tmpFrame[prop] for prop in propListDic]))
                        if 'displacement' in weighting:
                            inter_err_bool,uxOut1,uyOut1,uzOut1,distance1 = interpolate2DSpace(feDat,dicDat1,'displacement')
                            numInterp = numInterp+1
                            inter_err_bool,uxOut2,uyOut2,uzOut2,distance2 = interpolate2DSpace(feDat,dicDat2,'displacement')
                            numInterp = numInterp+1      
#                            if frCount == len(feIn[job])-1:
#                                print('\nBLABLABLA DISP')
#                                print(feDat[0])
#                                print('bla')
                        elif 'strain' in weighting:
                            inter_err_bool,eps11Out1,eps22Out1,distance1 = interpolate2DSpace(feDat,dicDat1,'strain')
                            numInterp = numInterp+1
                            inter_err_bool,eps11Out2,eps22Out2,distance2 = interpolate2DSpace(feDat,dicDat2,'strain')
                            numInterp = numInterp+1       
#                            if frCount == len(feIn[job])-1:
#                                print('\nBLABLABLA STRAIN')
#                                print(feDat[0])
#                                print('bla')
#                                print(zip(eps11Out2))
                            
                        msgNum = numInterp/totNumInterp*100
                        if OutputFlag > 0: sys.stdout.write('\r%.2f'  % msgNum + '% done with interpolation')
            
                        #Geometry Matching Plots
                        if CreateInterpDistPlots:
                            f, (ax1, ax2) = plt.subplots(1,2, sharey=True)
                            ax1.scatter(dicIn[job][interpFrame[0]]['x'],dicIn[job][interpFrame[0]]['y'],c='blue',s=2)
                            ax1.scatter(frame['x'],frame['y'],c='red',s=2)
                            ax1.set_title('Interp 1')
                            ax2.scatter(dicIn[job][interpFrame[1]]['x'],dicIn[job][interpFrame[1]]['y'],c='blue',s=2)
                            ax2.scatter(frame['x'],frame['y'],c='red',s=2)
                            ax2.set_title('Interp 2')
                            f.suptitle('Geometry Comparison Frame '+str(frCount))
                            plt.savefig(ResDir+os.path.sep+'GeometryComparisonFrame'+str(frCount)+'.png')
                            if OutputFlag > 1: plt.show()
                            plt.close()
                            
                        
                        if CreateInterpDistPlots:
                            plt.figure()
                            plt.scatter(frame['x'],frame['y'],c=distance1)
                            sm = plt.cm.ScalarMappable(cmap=cm.jet)#, norm=plt.Normalize(vmin=minHeat, vmax=maxHeat))
                            # Create an empty array for the colorbar
                            sm._A = []
                            plt.colorbar(sm)
                            tText = 'Distance1_Scatter_Frame'+str(frCount)
                            plt.suptitle(tText)
                            plt.savefig(ResDir+os.path.sep+tText)
                            plt.close()
                        
                        
                        if CreateInterpDistPlots:
                            plt.figure()
                            plt.scatter(frame['x'],frame['y'],c=distance2)
            #                print('distance:')
            #                print(distance2)
            
                            sm = plt.cm.ScalarMappable(cmap=cm.jet)#, norm=plt.Normalize(vmin=minHeat, vmax=maxHeat))
                            # Create an empty array for the colorbar
                            sm._A = []
                            plt.colorbar(sm)
                            tText = 'Distance2_Scatter_Frame'+str(frCount)
                            plt.suptitle(tText)
                            plt.savefig(ResDir+os.path.sep+tText)
                            plt.close()
    
    
                        #  T I M E / E X T E N S O M E T E R   W I S E   I N T E R P O L A T I O N
                        #Calculate weights
                        T1 = (float(dicIn[job][interpFrame[1]]['frameExten']) - float(frame['frameDisp'])) / (float(dicIn[job][interpFrame[1]]['frameExten']) - float(dicIn[job][interpFrame[0]]['frameExten']))
                        T2 = (float(frame['frameDisp']) - float(dicIn[job][interpFrame[0]]['frameExten'])) / (float(dicIn[job][interpFrame[1]]['frameExten']) - float(dicIn[job][interpFrame[0]]['frameExten']))
#                        print('\n')
#                        print(T1)
#                        print(T2)
                        #interpolate
                        if 'displacement' in weighting:
                            uxOut = [ux1*T1 + ux2*T2 for ux1,ux2 in zip(uxOut1,uxOut2)]
                            uyOut = [uy1*T1 + uy2*T2 for uy1,uy2 in zip(uyOut1,uyOut2)]
#                            uzOut = [uz1*T1 + uz2*T2 for uz1,uz2 in zip(uzOut1,uzOut2)]
                        elif 'strain' in weighting:
                            eps11Out = [eps1*T1 + eps2*T2 for eps1,eps2 in zip(eps11Out1,eps11Out2)]
                            eps22Out = [eps1*T1 + eps2*T2 for eps1,eps2 in zip(eps22Out1,eps22Out2)]
#                            if frCount == len(feIn[job])-1:
#                                print('\nTIME INTERPOLATION:')
#                                
#                                print(eps22Out)
#                                import shit


                            
                        #scatter plot for last frame
                        if 'displacement' in weighting:
                            direction = ['ux','uy']
                            if frCount == len(feIn[job])-1:
                                plt.figure(figsize=(6,3))
                                for i,item in enumerate([uxOut,uyOut]):
                                    len(dicIn[job][frCount][direction[i]])
                                    len(item)
                                    plt.subplot(1,2,i+1)
                                    vmin = min(dicIn[job][frCount][direction[i]] + item)
                                    vmax = max(dicIn[job][frCount][direction[i]] + item)
                                    plt.scatter(dicIn[job][frCount]['x'],dicIn[job][frCount]['y'],3,dicIn[job][frCount][direction[i]],cmap='jet',vmin=vmin,vmax=vmax,alpha=1)
                                    plt.scatter(frame['x'],frame['y'],10,item,cmap='jet',vmin=vmin,vmax=vmax)                                            
                                    plt.title(direction[i])
                                    plt.xlabel('x-coordinate')
                                    plt.ylabel('y-coordinate')
                                    plt.axis([-0.5,max(dicIn[job][frCount]['x']),-0.5,max(dicIn[job][frCount]['y'])])
                                    plt.colorbar()
                                    # make sure colorbar and plots do not overlap:
                                    plt.tight_layout() 
                        elif 'damage' in weighting:
                            direction = ['damage']
                            if frCount == len(feIn[job])-1:
                                plt.figure(figsize=(6,3))
                                for i,item in enumerate([damageOut]):
                                    plt.subplot(1,2,i+1)
                                    vmin = min(dicIn[job][frCount][direction[i]] + item)
                                    vmax = max(dicIn[job][frCount][direction[i]] + item)
                                    plt.scatter(dicIn[job][frCount]['x'],dicIn[job][frCount]['y'],3,dicIn[job][frCount][direction[i]],cmap='jet',vmin=vmin,vmax=vmax,alpha=1)
                                    plt.scatter(frame['IntPoint_x'],frame['IntPoint_y'],10,item,cmap='jet',vmin=vmin,vmax=vmax)                                            
                                    plt.title(direction[i])
                                    plt.xlabel('x-coordinate')
                                    plt.ylabel('y-coordinate')
                                    plt.axis([-0.5,max(dicIn[job][frCount]['x']),-0.5,max(dicIn[job][frCount]['y'])])
                                    plt.colorbar()
                                    # make sure colorbar and plots do not overlap:
                                    plt.tight_layout() 
                        elif 'strain' in weighting:
                            direction = ['eps11','eps22']
                            if frCount == len(feIn[job])-1:
                                plt.figure(figsize=(6,3))
                                for i,item in enumerate([eps11Out,eps22Out]):
                                    if i == 0:
                                        plt.subplot(1,2,i+1)
                                        vmin = min(dicIn[job][frCount]['eps22'] + item)
                                        vmax = max(dicIn[job][frCount]['eps22'] + item)
    #                                    plt.scatter(dicIn[job][frCount]['x'],dicIn[job][frCount]['y'],3,dicIn[job][frCount][direction[i]],cmap='jet',vmin=vmin,vmax=vmax,alpha=1)
                                        plt.scatter(frame['IntPoint_x'],frame['IntPoint_y'],10,eps22Out,cmap='jet')#,vmin=vmin,vmax=vmax) 
#                                        plt.scatter(frame['IntPoint_x'],frame['IntPoint_y'],10,frame['strain22'],cmap='jet',vmin=vmin,vmax=vmax)                                            
                                        plt.title('DIC data')
                                        plt.xlabel('x-coordinate')
                                        plt.ylabel('y-coordinate')
                                        plt.axis([-0.5,max(dicIn[job][frCount]['x']),-0.5,max(dicIn[job][frCount]['y'])])
                                        plt.colorbar()
                                        # make sure colorbar and plots do not overlap:
                                        plt.tight_layout() 
                                    elif i == 1:
                                        plt.subplot(1,2,i+1)
                                        vmin = min(frame['strain22'])
                                        vmax = max(frame['strain22'])
    #                                    plt.scatter(dicIn[job][frCount]['x'],dicIn[job][frCount]['y'],3,dicIn[job][frCount][direction[i]],cmap='jet',vmin=vmin,vmax=vmax,alpha=1)
#                                        plt.scatter(frame['IntPoint_x'],frame['IntPoint_y'],10,item,cmap='jet')#,vmin=vmin,vmax=vmax) 
                                        plt.scatter(frame['IntPoint_x'],frame['IntPoint_y'],10,frame['strain22'],cmap='jet',vmin=vmin,vmax=vmax)                                            
                                        plt.title('FE DATA')
                                        plt.xlabel('x-coordinate')
                                        plt.ylabel('y-coordinate')
                                        plt.axis([-0.5,max(dicIn[job][frCount]['x']),-0.5,max(dicIn[job][frCount]['y'])])
                                        plt.colorbar()
                                        # make sure colorbar and plots do not overlap:
                                        plt.tight_layout() 
                                    
                        #plot differences
                        if 'displacement' in weighting:
                            direction = ['ux','uy']
                            if frCount == len(feIn[job])-1:
                                plt.figure(figsize=(6,3))
                                for i,item in enumerate([uxOut,uyOut]):
                                    plt.subplot(1,2,i+1)
                                    diff = np.array(item) - np.array(frame[direction[i]])
                                    plt.scatter(frame['x'],frame['y'],20,diff,cmap='jet')
                                    plt.title(direction[i] + ' difference')
                                    plt.xlabel('x-coordinate')
                                    plt.ylabel('y-coordinate')
                                    plt.axis('tight')
                                    plt.axis('equal')
                                    plt.colorbar()
                                    # make sure colorbar and plots do not overlap:
                                    plt.tight_layout() 
                        if 'strain' in weighting:
                            direction = ['strain11','strain22']
                            if frCount == len(feIn[job])-1:
                                plt.figure(figsize=(6,3))
                                for i,item in enumerate([eps11Out,eps22Out]):
                                    plt.subplot(1,2,i+1)
                                    diff = np.array(item) - np.array(frame[direction[i]])
                                    plt.scatter(frame['IntPoint_x'],frame['IntPoint_y'],20,diff,cmap='jet')
                                    plt.title(direction[i] + ' difference')
                                    plt.xlabel('x-coordinate')
                                    plt.ylabel('y-coordinate')
                                    plt.axis('tight')
                                    plt.axis('equal')
                                    plt.colorbar()
                                    # make sure colorbar and plots do not overlap:
                                    plt.tight_layout() 

                                    
                        # save plot
                        plt.savefig(GraphDir+os.path.sep+job+'_'+objectiveField+'_iter'+str(len(fvec))+'.pdf', bbox_inches='tight')
    
                        if 'displacement' in weighting:
                            #concatenate lists of x- and y-displacements
                            ExpObjectiveInterList = np.array(uxOut+uyOut)
                            FeObjectiveList = np.array(frame['ux']+frame['uy'])
                            
                        if 'strain' in weighting:
                            #concatenate lists of x- and y-strains
                            ExpObjectiveInterList = np.array(eps11Out+eps22Out)
                            FeObjectiveList = np.array(frame['strain11']+frame['strain22'])

    
                        #claculate field-related objective function value for current frame
                        fu = getObjectiveFunctionValueDIC(ExpObjectiveInterList,FeObjectiveList,normFlag)
                            
                        fFrame.append(fu)        
                                         
    
            #            --SUBSECTION LAST ITERATION--
                        if SaveLastIterationData == 1 or frCount > 0: #All this stuff is only done at the end
                            #Save interpolated outputs as a file
                            frameNum = int2strPad(frCount,len(str(len(feIn[job]))))
                            with open(ResDir+os.path.sep+'TempInterpDataFile_' + job + '_Frame' + frameNum + '.txt','w') as file:
                                #Write general information
                                if 'displacement' in weighting:
                                    tempZip = list(zip(frame['x'],frame['y'],frame['ux'],frame['uy'],uxOut,uyOut))
                                    file.write('x\ty\tuxFE\tuyFE\tuxDIC_interpolated\tuyDIC_interpolated\t\n')
                                if 'strain' in weighting:
                                    tempZip = list(zip(frame['IntPoint_x'],frame['IntPoint_y'],frame['strain11'],frame['strain22'],eps11Out,eps22Out))
                                    file.write('x_intPoint\ty_intPoint\teps11FE\teps22FE\teps11DIC_interpolated\teps22DIC_interpolated\t\n')
                                for point in tempZip:
                                    for pt in point:
                                        file.write(str(pt) + '\t')
                                    file.write('\n')
                                    
                                
                        if inter_err_bool:
                            message = 'Interpolation failed: FE-node might be outside DIC-data. Check geometry matching plot!'
                            sendWarning(message)
                            break
                #End of frame loop
                
                
    
                sys.stdout.write('\rInterpolation completed            \n')
        
        #        --SUBSECTION RESULTS--  
                
                #assemble FE-LDC:
                LDC = (np.array(LDC1),np.array(LDC2))
                
                
                if 'displacement' in weighting and weighting['displacement'] > 0:
                    if not inter_err_bool:
                        fQOI['displacement'][job] = sum(fFrame)/np.count_nonzero(fFrame)
                    else:
                        fQOI['displacement'][job] = 1000000000
                        #Geometry Matching Plots
                        plt.figure
                        plt.scatter(dicIn[job][interpFrame[0]]['x'],dicIn[job][interpFrame[0]]['y'],c='blue',s=2,label='DIC-data')
                        plt.scatter(frame['x'],frame['y'],c='red',s=2,label='FE-data')
                        plt.suptitle('Geometry matching plot')
                        plt.legend()
                        plt.axis('equal')
                        plt.show()
                        plt.close()
                elif 'displacemnt' not in weighting or weighting['displacement'] == 0:
                    #dummy for unused DIC objective function value
                    fQOI['displacement'][job] = -1
                  
                if 'strain' in weighting and weighting['strain'] > 0:
                    if not inter_err_bool:
                        fQOI['strain'][job] = sum(fFrame)/np.count_nonzero(fFrame)
#                    else:
#                        TO DO ....                        
                elif 'strain' not in weighting or weighting['strain'] == 0:
                    #dummy for unused DIC objective function value
                    fQOI['strain'][job] = -1
                    
                    
                if 'LDC' in weighting and weighting['LDC'] > 0:
                    #Create LDC temp file
                    with open('LoadDispCurve.txt','w') as file:
                        for entry in range(len(LDC[1])):
                            file.write(str(LDC[0][entry]) + '\t' + str(LDC[1][entry]) + '\n')
                            
         
                    # remove corrupt data from LDC
                    LDC     = RemoveCorruptLDCData(LDC)
                    
                    # ALEX: size changes during iterations, this is fixed only temporarily
                    # needs to be done soon
                    LDC_ref_backup = LDC_ref[job]
                    
                    # interpolate LDC data if necessary
                    [LDC,LDC_ref_backup] = AdjustData(LDC,LDC_ref_backup)  
                elif 'LDC' not in weighting or weighting['LDC'] == 0:
                    #dummy for unused LDC objective function value
                    fQOI['LDC'][job] = -1
                    
                # G E T   O B J E C T I V E    F U N C T I O N   V A L U E   W R T   D A M A G E
                if 'damage' in weighting and weighting['damage'] > 0:
                    # interpolate damage field data. Only for last frame of simulation. Make sure the specimen is unloaded in simulation
                    feDat = list(zip(*[frame[prop] for prop in ['IntPoint_x','IntPoint_y','damage']]))  #the * unpacks the zip-arguemnts from list
                    DMGDat = list(zip(*[DMG_ref[prop] for prop in ['x','y','damage']]))
                    inter_err_bool,damageOut,distance1 = interpolate2DSpace(feDat,DMGDat,'damage')
                    
                    fQOI['damage'][job] = getObjectiveFunctionValueDIC(np.array(damageOut),np.array(frame['damage']),normFlag)

                    
                PlotOrigLDC(LDC,LDC_ref_backup,FileName=job+'_LDC_iter_'+str(len(fvec)))
                if normFlag:
                    fQOI['LDC'][job] = getObjectiveFunctionValueLDC(LDC,LDC_ref_backup,DynamicWeight,1/LDC_ref_backup[1])
                
                else:
                    fQOI['LDC'][job] = getObjectiveFunctionValueLDC(LDC,LDC_ref_backup,DynamicWeight)

                # write header to console
                print('\nStats for job ' + job + ':')
                  
                f_dict[job] = 0
                if firstRunNorm:
                    # save objective function values of first iteration
                    if firstRun:
                        f_dict[job] = 1
                    else:
                        # sum up all individual objective fun. falues with associated weighting factors
                        for qoi in weighting.keys():
                            if weighting[qoi] > 0:
                                val = fQOI[qoi][job]/fQOIvec[0][qoi][job]
                                if OutputFlag > 0:
                                     print('normed f_' + qoi + ': %.5f' %val)
                                f_dict[job] = f_dict[job] + weighting[qoi]*val
                                
                elif not idealVector == {}:
                    tmpSum =  0
                    for qoi in idealVector.keys():
                        tmpSum = tmpSum + ((fQOI[qoi][job] - idealVector[qoi]) / idealVector[qoi])**2
                    f_dict[job] = np.sqrt(tmpSum)
                else:
                    f_dict[job] = 0
                    # sum up all individual objective fun. falues with associated weighting factors
                    for qoi in weighting.keys():
                        if weighting[qoi] > 0:
                            f_dict[job] = f_dict[job] + weighting[qoi]*fQOI[qoi][job]
                            
                #calculate overall objective function value
                f = f + f_dict[job]
                
                os.chdir('..')   
            #End of job loop
            
            f = f/len(JobNames)
            
            if firstRun:
                firstRun = False
            
#            #print results to console
            for job in JobNames:
                if OutputFlag > 0:

                    for qoi in weighting.keys():
                        print('f_' + qoi + ': %.5f' %fQOI[qoi][job])                    
                
        elif err_bool == True: #jobs did not finish properly
            f = 1000000000
            for qoi in weighting.keys():
                fQOI[qoi] = dict.fromkeys(JobNames,1000000000)


    
        if OutputFlag > 0:
            print('\nGlobal stats:')
            print('Global objective function value = ' + str(f) + '\n')
        
        
    # append current objective function values to vector
    fvec.append(f)
    fQOIvec.append(fQOI)
    
    UpdateIterfile(x,f,fQOI,len(fvec))
#    WriteObjFunctionData(f,fDIC,fLDC)
    if len(fvec) > 5:
        PlotObjFunction(fvec)
    if oneIter:
        sys.exit()
    
    return f



def patchFillSort(pointsList):
    if len(pointsList) == 4:
        sortPointsList = [pointsList[0],pointsList[1],pointsList[3],pointsList[2]]
    elif len(pointsList) == 3:
        sortPointsList = pointsList
    else:
        print('Unsupported element shape')
    
    return(sortPointsList)



def PatchPlots(inDict, nameStr):
#    plotTypeList = ['strainX','strainY','ux','uy','uMag','uxInt','uyInt','uMagInt']
#    plotTypeList = ['uy','uyInt','ux','uxInt']
    plotTypeList = ['uy','uyInt']
    for plotType in plotTypeList:
        if plotType =='uy':
            name = nameStr + 'uy FE-data'
        elif plotType=='uyInt':
            name = nameStr + 'uy DIC-data mapped to FE-nodes'
        elif plotType=='ux':
            name = nameStr + 'ux x FE-data'
        elif plotType=='uxInt':
            name = nameStr + 'ux DIC-data mapped to FE-nodes'
        fig, ax = plt.subplots()
        patchHold = []
        patchVal = []
        tX,tY = [],[]
        for i,patch in enumerate(inDict):
            tX.append([x[0] for x in patch['NodeDat']])
            tY.append([y[1] for y in patch['NodeDat']])
            patchXY = np.zeros([len(patch['NodeDat']),2])
            patchVal.append(patch[plotType])
            for j,point in enumerate(patch['NodeDat']):
                patchXY[j,:] = [point[0],point[1]]
            #Sort patchXY so the patches are square
            if AbqMode == '3D': patchXY = patchFillSort(patchXY)
            patchHold.append(Polygon(patchXY, True))
        p = PatchCollection(patchHold, cmap=cm.jet,norm=plt.Normalize(vmin=min(patchVal),vmax=max(patchVal)))
        p.set_array(np.array(patchVal))
        ax.add_collection(p)
        xMin,xMax = min(min(tX)),max(max(tX))
        yMin,yMax = min(min(tY)),max(max(tY))
        ax.set_xlim([xMin,xMax])
        ax.set_ylim([yMin,yMax])
        plt.suptitle(name)
        fig.colorbar(p,ax=ax)
        plt.savefig(GraphDir+os.path.sep+name+'.png')
        if OutputFlag > 1: plt.show()
        plt.close()
    
    # ALEX: plot differences in y-displacements
    name = nameStr + ' y-displacement Difference'
    fig, ax = plt.subplots()
    patchHold = []
    patchVal = []
    tX,tY = [],[]
    for i,patch in enumerate(inDict):
        tX.append([x[0] for x in patch['NodeDat']])
        tY.append([y[1] for y in patch['NodeDat']])
        patchXY = np.zeros([len(patch['NodeDat']),2])
        patchVal.append(patch['uy']-patch['uyInt'])
        for j,point in enumerate(patch['NodeDat']):
            patchXY[j,:] = [point[0],point[1]]
        #Sort patchXY so the patches are square
        if AbqMode == '3D': patchXY = patchFillSort(patchXY)
        patchHold.append(Polygon(patchXY, True))
    p = PatchCollection(patchHold, cmap=cm.jet,norm=plt.Normalize(vmin=min(patchVal),vmax=max(patchVal)))
    p.set_array(np.array(patchVal))
    ax.add_collection(p)
    xMin,xMax = min(min(tX)),max(max(tX))
    yMin,yMax = min(min(tY)),max(max(tY))
    ax.set_xlim([xMin,xMax])
    ax.set_ylim([yMin,yMax])
    plt.suptitle(name)
    fig.colorbar(p,ax=ax)
    plt.savefig(GraphDir+os.path.sep+name+'.png')
    if OutputFlag > 1: plt.show()
    plt.close()
    
    # ALEX: plot differences in x-displacements
    name = nameStr + ' x-displacement Difference'
    fig, ax = plt.subplots()
    patchHold = []
    patchVal = []
    tX,tY = [],[]
    for i,patch in enumerate(inDict):
        tX.append([x[0] for x in patch['NodeDat']])
        tY.append([y[1] for y in patch['NodeDat']])
        patchXY = np.zeros([len(patch['NodeDat']),2])
        patchVal.append(patch['ux']-patch['uxInt'])
        for j,point in enumerate(patch['NodeDat']):
            patchXY[j,:] = [point[0],point[1]]
        #Sort patchXY so the patches are square
        if AbqMode == '3D': patchXY = patchFillSort(patchXY)
        patchHold.append(Polygon(patchXY, True))
    p = PatchCollection(patchHold, cmap=cm.jet,norm=plt.Normalize(vmin=min(patchVal),vmax=max(patchVal)))
    p.set_array(np.array(patchVal))
    ax.add_collection(p)
    xMin,xMax = min(min(tX)),max(max(tX))
    yMin,yMax = min(min(tY)),max(max(tY))
    ax.set_xlim([xMin,xMax])
    ax.set_ylim([yMin,yMax])
    plt.suptitle(name)
    fig.colorbar(p,ax=ax)
    plt.savefig(GraphDir+os.path.sep+name+'.png')
    if OutputFlag > 1: plt.show()
    plt.close()
    
    
def PlotOrigLDC(Data1,Data2,FileName='test'):
    # Plot original Load Displacement Curve without any manipulations
    # (only absolute values)

    # rearrange data
    x1 = abs(Data1[0])
    y1 = abs(Data1[1])
    x2 = abs(Data2[0])
    y2 = abs(Data2[1])

    plt.figure()
    plt.plot(x2,y2,linewidth=2.0)
    plt.plot(x1,y1,linewidth=2.0)
    plt.grid(True)
    #plot region used for optimization
    ax = plt.gca()
    ylim = ax.get_ylim()
    xlim = ax.get_xlim()
    # draw lines and rectangles
    plt.plot((x_start,x_start),ylim,linewidth=0.5,color='k')
    ax.add_patch(Rectangle((xlim[0],ylim[0]),x_start-xlim[0],ylim[1]-ylim[0],color='grey',alpha=0.5,zorder=10))
    if x_end < xlim[1]:
        plt.plot((x_end,x_end),ylim,linewidth=0.5,color='k')
        ax.add_patch(Rectangle((x_end,ylim[0]),xlim[1]-x_end,ylim[1]-ylim[0],color='grey',alpha=0.5,zorder=10))
        
    #reset old axis limits
    ax.set_ylim(ylim) 
    ax.set_xlim(xlim) 
       
    plt.xlabel('Displacement in mm')
    plt.ylabel('Force in N')
    plt.legend(['Experimental Data','Simulation Results'],loc='lower center')
#    plt.savefig(FileName+'.eps')
    plt.savefig(LdcDir + os.path.sep + FileName+'.pdf', bbox_inches='tight')
    if OutputFlag >= 1: plt.show()
    plt.close()
#    plt.figure()
#    plt.plot(x2,y2,linewidth=1.5)
#    plt.plot(x1,y1,linewidth=1.5)    
#    plt.grid(True)
#    plt.xlabel('x')
#    plt.ylabel('y')
#    plt.legend(['Experimental Data','Simulation Results'],loc='upper left')
#    plt.savefig(FileName+'_Latex.eps')
#    plt.savefig(FileName+'_Latex.png')


def RayTracingMethod(x,y,poly):
# checks wether query point is in specified polygon
# returns False if query point is on one of the edges
    n = len(poly)
    inside = False
    p1x,p1y = poly[0]
    for i in range(n+1):
        p2x,p2y = poly[i % n]
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xints = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xints:
                        inside = not inside
        p1x,p1y = p2x,p2y

    return inside


def ReadConnectivity(fName):
    with open(fName) as file:
        ReadIn = 0 #0 for none, 1 for Node, 2 for Element, 3 for Done
        connectList, coordList = [], []
        for line in file:
            if ReadIn == 0:
                if '*Node' in line and not 'Output' in line: ReadIn = 1
            elif ReadIn == 1:
                if '*Element,' in line: ReadIn = 2
                else:
                    line = line.strip()
                    temp = line.split(',')
                    temp = [x.strip() for x in temp]
                    temp = [float(x) for x in temp]
                    coordList.append(temp[1:])
            elif ReadIn == 2:
                if '*Nset' in line or '*Elset' in line: ReadIn = 3
                else:
                    line = line.strip()
                    temp = line.split(',')
                    temp = [x.strip() for x in temp]
                    temp = [int(x) for x in temp]
                    connectList.append(temp[1:])
    return connectList, coordList

def ReadDMGFile(filename):
    outDict = {}
    for key in ['x','y','damage']: outDict[key] = []

    with open(filename) as file:
        for lCount, line in enumerate(file):
            temp = line.split(',')
            temp = [x.strip() for x in temp]
            for i,key in enumerate(['x','y','damage']):
                outDict[key].append(float(temp[i]))
    file.close()
    return outDict
    
def ReadDICData(filename,properties,headerlines):
    frameSave = {}
    headerErr = False
    for prop in properties: frameSave[prop] = []
    with open(filename) as file:
        for lCount, line in enumerate(file):
            if validationFlag:
                headerlines = 1

            #check header of input file for tags
            if lCount == headerlines-1:
                if 'displacement' in weighting:
                    if not any([item in line for item in ['disp','ux']]):
                        headerErr = True
                if 'strain' in weighting:
                    if not any([item in line for item in ['epsilon','strain','eps','LE']]):
                        headerErr = True
                if headerErr:
                    message = 'Current quantity-of-interest was not found in header of DIC input file! Check DIC files and number of header lines'
                    sendWarning(message)
            if lCount > headerlines-1:
                temp = line2listDIC(line)
                for kCount,keys in enumerate(properties):
                    frameSave[keys].append(float(temp[kCount]))

    file.close()            
    #Empty values for overwriting            
#    frameSave['frameForce'] = 'nan'
#    frameSave['frameExten'] = 'nan'
    #Check to see if this is a frame of zero data, in which case it shouldn't be saved
#    if sum([abs(y) for y in frameSave['uy']]) != 0: UNCOMMENT AND TAB NEXT LINE TO FILTER ZERO FRAMES
    return frameSave



def ReadFEData(filename):
    frameSave = {}
    with open(filename) as file:
        for lCount, line in enumerate(file):
            if not lCount == 0 and not lCount == 1:
                temp = line2listFE(line)
                for kCount,keys in enumerate(enumPropList):
                    frameSave[keys].append(float(temp[kCount]))
            elif lCount == 0:
                temp = line2listFE(line)
#                frameKeyList = ['frameTime','frameForce','frameExten','frameDisp']
                frameKeyList = ['frameTime','frameForce','frameDisp'] # ALEX: extensometer data not necessary anymore
                for kCount, keys in enumerate(frameKeyList):
                    frameSave[keys] = float(temp[2*kCount+1])
            else:
                temp = line2listFE(line)
                enumPropList = temp
                for item in temp:
                    frameSave[item]=[]
    file.close()
#    if sum([abs(y) for y in frameSave['uy']]) != 0: UNCOMMENT AND TAB NEXT LINE TO FILTER ZERO FRAMES
    return frameSave



def ReadFEDataElem(filename):
    frameSave = {}
    with open(filename) as file:
        for lCount, line in enumerate(file):
            temp = line2listFE(line)
            if not lCount == 0 and not lCount == 1:
                for kCount,keys in enumerate(enumPropList):
                    frameSave[keys].append(float(temp[kCount]))
            elif lCount == 1:
                enumPropList = temp
                for item in temp:
                    frameSave[item]=[]
    file.close()
    return frameSave



def ReadLDCFile(filename):
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
            temp1 = line.split(',')
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

def ReadRestartFile(filename):
    with open(filename) as file:
        f_list = []
        f_ini_dict = defaultdict(dict)
        
        # read initial guess
        temp = file.readline().split('[')
        temp = temp[-1].split(']')[0]
        x0 = [float(x) for x in temp.split(',')]
        # read weighting factor and firstRunNorm-flag
        temp = file.readline().split('norm to 1st run value:')
        firstRunNorm = bool(temp[-1])
        temp = temp[0].split('Weighting: ')[-1]
        exec('weighting =' + temp)

        #skip header line
        file.readline()
        
        if firstRunNorm:
            # read all objective function values of first iteration
            temp = file.readline().split(']')[-1]
            temp = temp.split()
            f_list.append(float(temp.pop(0)))
            for qoi in weighting.keys():
                for job in JobNames:
                    f_ini_dict[qoi][job] = float(temp.pop(0))
        # read all objective function values
        for line in file:
            temp = line.split(']')[-1]
            f_list.append(float(temp.split()[0]))

    file.close()
    
    return x0, f_list, f_ini_dict, weighting, firstRunNorm

#Adapt to my data        
def RemoveCorruptLDCData(Data,DelTol=100.0):
    # delete corrupt data from array
    # delete all displacement data, that are greater than chosen tol. (DelTol)
    Data1 = Data[0]
    Data2 = Data[1]

    for i,e in enumerate(Data1):
        if i > len(Data1):
            break
        if Data1[i] > DelTol:
            Data1 = np.delete(Data1,i)
            Data2 = np.delete(Data2,i)
            i = i-1                  
    return Data1,Data2



def SaveLDC_InitialGuess(Data1,Data2):
    # save Load displacement curve for initial guess
    global SaveInitGuessData#Need to declare global to read it in
    if SaveInitGuessData == 0:
        shutil.copy('LoadDispCurve.txt','LoadDispCurve_InitialGuess.txt')
        UpdateLogfile('\nData for initial guess has been saved.')
        SaveInitGuessData = 1
        # save load displacement curve for reference data
        # for interpolated data points
        shutil.copy('LoadDispCurve.txt','LoadDispCurve_ref_Interp.txt')



def SaveLDC_Optimum(Data):
    # save load displacement curve for optimized parameter set
    global SaveLastIterationData
    if SaveLastIterationData == 1:
        shutil.copy('LoadDispCurve.txt','LoadDispCurve_Optimum.txt')
        UpdateLogfile('\nData for optimized parameters has been saved.')

def sendError(message='An error has occured'):
    print('\n')
    print(Fore.RED + 'ERROR:  ' + message)
    print('\n')
    UpdateLogfile('ERROR:  ' + message)
    sys.exit()

def sendWarning(message='A warning has occured',suppress = False):
    print(Fore.GREEN + 'WARNING:  ' + Style.RESET_ALL + message + suppress*'(This warning will be suppressed from now on)')
    UpdateLogfile('WARNING:  ' + message + suppress*'(This warning will be suppressed from now on)')

def SetUpFlowCurve(x,ext='.inp'):
    # Set up a flow curve file containing the flow curve as tabulated data
    # to be used in the Abaqus input file
    eps_min = 0.0
    eps_max = 1.0
    n = 1000
    for job in JobNames:
        os.chdir(job)
        file = open('FlowCurve'+ext,'w+')
        
        for i in range(n+1):
            eps = (eps_max-eps_min)*i/n
            sig = x[0] * (eps + x[1]) ** x[2]
            file.write(str(sig) + ',' + str(eps) + '\n')
    
        file.close()
        os.chdir('..')
    

#Look for more elegant solution to this function    
def SetUpParameters_old(x,ParString,FileName,ext='.for'):
    # Set up all parameters (x) in the input or user material files
    # by replacing the strings (ParString) with the associated parameter (x).
    
    # x: parameter vector
    # ParString: list containing all parameter names as strings
    # FileName: name of the file to be manipulated
    #           it has to be named FileName_base
    # ext: file extension, can be either '.inp' or '.for'
    
    # if more than one parameter needs to be identified
    for job in JobNames:
        os.chdir(job)
        if len(ParString) > 1:
            for i,e in enumerate(ParString):
                # for the first element only
                if i == 0:
                    with open(FileName+'_v'+str(i)+ext,'wt') as fout:
                        with open(FileName+'_base'+ext,'rt') as fin:
                            for line in fin:
                                fout.write(line.replace(e, str(x[i])))
                    fin.close()
                    fout.close()
    
                # for any element between the first and the last
                elif i > 0 and i < len(ParString)-1:
                    with open(FileName+'_v'+str(i)+ext,'wt') as fout:
                        with open(FileName+'_v'+str(i-1)+ext,'rt') as fin:
                            for line in fin:
                                fout.write(line.replace(e, str(x[i])))
                    fin.close()
                    fout.close()
    
                # for the last element only
                elif i == len(ParString)-1: 
                    with open(FileName+ext,'wt') as fout:
                        with open(FileName+'_v'+str(i-1)+ext,'rt') as fin:
                            for line in fin:
                                fout.write(line.replace(e, str(x[i])))
                    fin.close()
                    fout.close()
    
            # delete all redundant files (if they exist)
            for i in range(len(ParString)-1):
                try:
                   os.remove(FileName+'_v'+str(i)+ext)
                except OSError:
                   pass
    
        # if only one parameter needs to be identified
        else: 
            with open(FileName+ext,'wt') as fout:
                with open(FileName+'_base'+ext,'rt') as fin:
                    for line in fin:
                        fout.write(line.replace(ParString[0], str(x[0])))
            fin.close()
            fout.close()
        
        if LogOutputFlag > 0:
            UpdateLogfile('Parameters have been set up (' + str(x) + ')')
        if OutputFlag > 0:
            sys.stdout.write('Parameters have been set up (' + str(x) + ')\n')
    
        os.chdir('..')
        
def SetUpParameters(x,ParString,FileName,ext='.for',xFix=[],fixParString=[]):
    # Set up all parameters (x) in the input or user material files
    # by replacing the strings (ParString) with the associated parameter (x).
    
    # x: parameter vector
    # ParString: list containing all parameter names as strings
    # FileName: name of the file to be manipulated
    #           it has to be named FileName_base
    # ext: file extension, can be either '.inp' or '.for'
    
    # if more than one parameter needs to be identified
    x = np.concatenate((x,xFix),0)
    ParString = ParString + fixParString
    for job in JobNames:
        os.chdir(job)

        with open(FileName + ext,'wt') as fout:
            with open(FileName + '_base' + ext,'rt') as fin:
                for line in fin:
                    found = False
                    for i,par in enumerate(ParString):
                            if par in line:
                                found = True
                                fout.write(line.replace(par, str(x[i])))
                    if not found:
                        fout.write(line)
                            
        fin.close()
        fout.close()

        if LogOutputFlag > 1:
            UpdateLogfile('Parameters have been set up (' + str(x) + ')')
        if OutputFlag > 1:
            sys.stdout.write('Parameters have been set up (' + str(x) + ')\n')
    
        os.chdir('..')

#Write a file to pass variables to Abaqus
#CURRENT IMPLEMENTATION RELIES ON THE SET NAMES BEING MOSTLY GLOBAL SINCE DEFINED IN THE MAIN BLOCK
#THIS IS NOT AWESOME PROGRAMMING AND IT WOULD BE SLIGHTLY BETTER IF EACH OF THE PARAMETERS WAS PASSED TO THE FUNCTION
def SetUpPPScript(PostProcessingScript,JobName,RunName):
    replaceList = [['JobName_par',JobName],['RunName_par',RunName],['NeckNodeSet_par',NeckNodeSet],['ForceSet_par',ForceSet],
                   ['LoadDir_par',LoadDirFE],['DataPrec_par',DataPrec],['RawDir_par',RawDir],
                   ['NeckElemSet_par',NeckElemSet],['DispScaleX_par',DispScaleX],['DispScaleY_par',DispScaleY]]
    if OutputFlag > 0:
        print('\nThe following parameters are being updated in the post processing script:')
        [print('\t-' + x[0]) for x in replaceList]
        print('')
    with open(PostProcessingScript+'.py','w') as fout, open(PostProcessingScript+'_base.py','r') as fin:
        for line in fin:
            changed = False
            for par in replaceList:
                if type(par[1]) != str: par[1] = str(par[1])
                if par[0] in line:
                    tempLine = line.replace(par[0], '\'' + par[1] + '\'')
                    changed = True
                elif not changed:
                    tempLine = line
            fout.write(tempLine)
    fin.close()
    fout.close()




def UpdateIterfile(x,f,xfQOI,niter):
    # create iter file
    if not os.path.exists(cwd + os.path.sep + RunName+'_iter.txt'):
        Iter = open(cwd + os.path.sep + RunName+'_iter.txt', 'w')
        Iter.write('# Starting vector: ' + str(x0) + '\n')
        Iter.write('# Weighting: ' + str(weighting) + '\tnorm to 1st run value: ' + str(firstRunNorm) + '\n')
        Iter.write('# Iteration \t[Parameter vector]\tObjective function value')
        for qoi in weighting.keys():
            for job in JobNames:
                Iter.write('\tf_'+qoi+'_'+job)
        Iter.write('\tCPU time [hhhh:mm:ss]\n')
        Iter.close()
    
    # Update log-file: parameter vector, objective function value, CPU time
    Iter = open(cwd + os.path.sep + RunName+'_iter.txt', 'a+')
    t = (time.clock()-StartTime)
    h = int(t/60/60)                # hours
    m = int(t/60 - h*60)            # minutes
    s = int(t - h*60*60 - m*60)     # seconds
    
    #write iteration number
    Iter.write(str(niter) + '\t')
    #write parameter vector by loop to avoid new line (as done by print(x))
    Iter.write('[')
    for i in x:
        Iter.write('% 10.8E' %i +' ')
    Iter.write(']')
    #write objective function values
    Iter.write('\t' + '% 10.6f' %f + '\t')
    for qoi in weighting.keys():
        for job in JobNames:
            Iter.write('% 10.8E' %xfQOI[qoi][job] +' ')
    #write timestamp
    Iter.write('\t' + '%04d' %h + ':' + '%02d' %m + ':' + '%02d' %s + '\n')
    Iter.close()
    
    
    
def UpdateLogfile(message):#Consider passing urgency to each message so amount of info can be controlled
	log = open(cwd + os.path.sep + RunName+'_log.txt', 'a+')
	log.write( str(time.strftime("%c")) + ':  ' + message + '\n')
	log.close()
 

def WriteObjFunctionData(f,fDIC,fLDC):
    
    # append current f,fDIC,fLDC to file
    padLen = 4
    count = len(fvec)-1
    ObjFunFile = open(ResDir+os.path.sep+RunName+'_ObjectiveFunctionValues.txt', 'a')
    ObjFunFile.write('Iteration ' + int2strPad(count,padLen) + ':\t' + str(f) + '\n')
    ObjFunFile.close()
    
    ObjFunFile = open(ResDir+os.path.sep+RunName+'_ObjectiveFunctionValuesDIC.txt', 'a')
    if count == 0:
        ObjFunFile.write('Iteration')
        for job in JobNames: ObjFunFile.write('\t'+job)
    line = '\n' + int2strPad(count,padLen)
    for job in JobNames: line = line + '\t' + str(fDIC[job])
    ObjFunFile.write(line)
    ObjFunFile.close()

    ObjFunFile = open(ResDir+os.path.sep+RunName+'_ObjectiveFunctionValuesLDC.txt', 'a')
    if count == 0:
        ObjFunFile.write('Iteration')
        for job in JobNames: ObjFunFile.write('\t'+job)
    line = '\n' + int2strPad(count,padLen)
    for job in JobNames: line = line + '\t' + str(fLDC[job])
    ObjFunFile.write(line)
    ObjFunFile.close()   
    

def PlotObjFunction(vect,log=True):
    #plots the objective function value over objective function calls
    
    #remove aborted runs
    vect = [x for x in vect if x != 1000000000]
    plt.figure()
    plt.plot(vect)
    plt.grid(True)
    if log: plt.semilogy()
    plt.ylabel(log*'log10(' + 'Objective function value' + log*')')
    plt.xlabel('Objective function calls')
    plt.suptitle('Objective Function Value over Calls')
    plt.savefig(ResDir+os.path.sep+RunName+'_ObjectiveFunctionValues_Graph.png')
    plt.show()
    plt.savefig('ObjectiveFunction.eps')
    
def PlotObjFunctionDict(vect,ftype,log=True):
    #plots specific objective function values of different jobs over function calls
    
    #remove aborted runs
    plotdict = {}
    for job in JobNames: plotdict[job] = []
    for element in vect:
        if not any([x == 1000000000 for x in element.values()]):
            for job in JobNames:
                plotdict[job].append(element[job])
    plt.figure()
    for job in JobNames:
        plt.plot(plotdict[job],label=job)
    if log: plt.semilogy()
    plt.ylabel(log*'log10(' + 'Objective function value' + log*')')
    plt.xlabel('Objective function calls')
    plt.suptitle('Objective Function Value over Calls '+ftype)
    plt.savefig(ResDir+os.path.sep+RunName+'_ObjectiveFunctionValues'+ftype+'_Graph.png')
    plt.legend()
    plt.show()
    plt.savefig('ObjectiveFunction'+ftype+'.eps')

def WriteResultsFile(x0,xopt,ParString):
    # write all parameters to be identified, the initial guess
    # and the optimized parameter set into a file
    
    results = open(ResDir+os.path.sep+RunName+'_results.txt', 'w')
    results.write('Parameters: ' + str(ParString))
    results.write('\nInitial guess: ' + str(x0))
    results.write('\nOptimized parameterset: ')
    results.write('[')
    for i in xopt:
        results.write('% 10.8E' %i +' ')
    results.write(']')
    results.close()    
    



#==============================================================================
# MAIN PROGRAM
#==============================================================================
try:
    start_time = time.time()
    global fvec
    global LDC_ref
    global f_ini_dict
    
    # vector containing the total objective function value of all iterations
    fvec = []
    # vector containing the objective function values of specific qunatities of interest
    fQOIvec = []

    
    
    # CHECK FOR RESTART FILE
    if not restartFile == None:
        try:
            x0, fListRestart, fQOI_ini, weighting, firstRunNorm = ReadRestartFile(restartFile)
            message = 'Restart file has been read. Using ' + str(len(fListRestart)) + ' iterations from restart file.'
            if OutputFlag > 0: print(message)
            if LogOutputFlag > 0: UpdateLogfile(message)
                
        except FileNotFoundError:
            message = 'Specified restart file not found! Beginning optimisation, starting from inital guess.'
            sendWarning(message)
        


    # create folder structure
    directories = (getDirs(JobNames))
    # create folder, if it does not exist already
    makeDirList = [TempDir,ResDir,GraphDir,LdcDir]
    for job in JobNames:
        makeDirList.append(os.path.join(os.path.join(TempDir,job),RawDir))
    
    for dirName in makeDirList:
        if not os.path.exists(dirName):
            os.makedirs(dirName)
    
 
    # copy files to temp directory
    for i,job in enumerate(JobNames):
        os.chdir(directories[job])
        
        # check wether DIC data is available
        foundDICData = False
        for roots, dirs, files in os.walk(ExpFolder):
            padLen = len(str(len(files)))
            fcount = 0
            for frames in files:
                if FrameLabelsDIC[i] in frames:
                    foundDICData = True
        if not foundDICData and not weighting == 0:
            message = 'No DIC data found for job ' + str(job) + '! Check DIC-labels in line 81 and experimental results folder or set weighting to zero.'
            sendError(message)
            
        if UserMat:
    
            FileList = [job+'.inp','getAbqDisp_base.py',
                    'LoadDispCurve_ref.txt',job+'.bat',UserMat+'_base.for']
        else:
            FileList = [job+'.inp','getAbqDisp_base.py',
                    'LoadDispCurve_ref.txt',job+'.bat']   
             
        CopyFilesToTempDir(FileList,FrameLabelsDIC[i],LabelsLDC[i],LabelsDMG[i],job)
        os.chdir('..')
    
    
    
    # get connectivity lists and sets for all jobs and set up postprocessing scripts and batch files
    connectivity = {}
    NeckNodeSetList = {}
    ForceSetList = {}
    NeckElemSetList = {}
    LDC_ref = {}
    NeckNodeSet = NeckNodeSet.upper()
    ForceSet = ForceSet.upper()
    NeckElemSet = NeckElemSet.upper()
    for i,job in enumerate(JobNames):
        file = os.path.join(TempDir,job,job+'.inp')
        connectivity[job]    = ReadConnectivity(file)[0]
        ForceSetList[job]    = findSetPoints(ForceSet,file)
        if 'displacement' in weighting and weighting['displacement'] > 0:
            NeckNodeSetList[job] = findSetPoints(NeckNodeSet,file)
            NeckElemSetList[job] = findSetPoints(NeckElemSet,file)
        
        # change to directory of each job in temp folder
        os.chdir(os.path.join(TempDir,job))
        
        # get LDC  data from experiments
        LDC_ref[job] = ReadLDCFile_txt(RawDir+os.path.sep+LabelsLDC[i])
        
        # get damage field data from experiments
        DMG_ref = ReadDMGFile(RawDir+os.path.sep+LabelsDMG[i])

        # set up post processing script
        SetUpPPScript(PostProcessingScript,job,RunName)
        
        # write batch files for each job
        StartJob = AbaqusCmd + ' job=' + job + UserMatABQ + ' cpus=' + str(NumCPUs)
        batchfile_run=open('run_'+job+'.bat','w+')
        batchfile_run.write('@ECHO off\n')
        if len(JobNames) > 1:
            batchfile_run.write('CALL ' + StartJob + ' double=both')
        elif len(JobNames) == 1:
            batchfile_run.write('CALL ' + StartJob + ' double=both interactive ask_delete=OFF')
        batchfile_run.close()

        batchfile_kill=open('kill_'+job+'.bat','w+')
        batchfile_kill.write('@ECHO off\n')
        batchfile_kill.write('CALL ' + AbaqusCmd + ' job=' + job +' terminate')
        batchfile_kill.close()
        
        os.chdir('..')
        
    # work from TempDir
    os.chdir(TempDir)
    # optimization process     
    if EnableAllGraphs: SaveLastIterationData = 1
    
    #generate Datapoints for Pareto minimization
    
#    x0 = [ 1.27136985E+03,  1.28480709E-05,  1.46127748E-01 ]
#    
#    delta = 1
#    
#    p1_min = x0[0]*(1+delta/100)
#    p1_max = x0[0]*(1-delta/100)
#
#    p2_min = x0[1]*(1+delta/100)
#    p2_max = x0[1]*(1-delta/100)
# 
#    p3_min = x0[2]*(1+delta/100)
#    p3_max = x0[2]*(1-delta/100)
#    
#    numiter = 10
#
#    p1 = np.linspace(p1_min,p1_max,numiter).tolist()
#    p2 = np.linspace(p2_min,p2_max,numiter).tolist()
#    p3 = np.linspace(p3_min,p3_max,numiter).tolist()
#
#    p = np.array([p1,p2,p3],dtype=object)
#    p.shape
#  
#     
#    for i in p[0]:
#        for j in p[1]:
#            for k in p[2]:
#                paramset = np.array([i,j,k]).tolist()
#                print(paramset)
#                print('test')
#                optimfunc(paramset)
#                #(x0, f_list, f_ini_dict, weighting, firstRunNorm) = ReadRestartFile("E:\Arbeit_IUL\Parameteridentifikation_08102020\create_Pareto_Output\R15 and R5\R15 and R5_iter.txt")
#                #print(f_list)
#                #print(fQOIvec)
    
    #result = sopt.differential_evolution(optimfunc,bnds,args=(bnds),maxiter=IterNum, popsize=15, tol=0.01, mutation=(0.5, 1), recombination=0.7, seed=None, callback=None); xopt = result.x
    (xopt, fopt, iters, funcalls, warnflag, allvecs) = sopt.fmin(optimfunc,x0,args=(bnds),xtol=0.0001,ftol=0.0001,maxiter=IterNum,maxfun=None,full_output=1,disp=1,retall=1,callback=None)
    #(xopt, fopt, funcalls, gradcalls, warnflag, allvecs) = sopt.fmin_cg(optimfunc,x0,args=(bnds),gtol=1e-8,epsilon=epsilon,maxiter=2000,full_output=1,disp=1,retall=1,callback=None)
    #(xopt, fopt, direc, iters, funcalls, warnflag, allvecs) = sopt.fmin_powell(optimfunc,x0,args=(bnds),xtol=0.0001,ftol=0.0001,maxiter=2000,maxfun=None,full_output=1,disp=1,retall=1,callback=None)
    #fvec.append(fopt)
    
    
    UpdateLogfile('\nOptimization finished successfully.')
    
    #==============================================================================
    # POST PROCESSING
    #==============================================================================
    SaveLastIterationData = 1
    
    # run simulation with optimized parameterset and save all requested outputs
#    optimfunc(xopt)
    
    
#    WriteResultsFile(x0,xopt,ParString)
    
    PlotObjFunction(fvec)
#    for qoi in weighting.keys():
#        PlotObjFunctionDict(fQOIvec[qoi])
    
    # change to original working directory
    os.chdir(cwd)
except KeyboardInterrupt:
    message = 'Parameter identification script interrupted by user.'
    print('\n'+message+'\n')
    UpdateLogfile(message)
    kill_all_abaqus(JobNames)
    sys.exit()