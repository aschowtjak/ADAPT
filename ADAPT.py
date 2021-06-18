# Authors: Alexander Schowtjak, Robin Schulte, Robin Gitschel
# e-mail:   alexander.schowtjak@tu-dortmund.de
#           robin.schulte@tu-dortmund.de
# Creative Commons Attribution 3.0 Unported License.

# You should have received a copy of the license along with this
# work.  If not, see <http://creativecommons.org/licenses/by/3.0/>.

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
import platform
import matplotlib.path as mpltPath
from random import uniform

# define global variables
global fLDC0, fDIC0
global feIn
global colorterm
global weighting
global DisableAllGraphs
global firstRun
global x_start
global x_end
global TempDir
global cwd
global DICPropList 


#==============================================================================
############################ START OF USER INPUT ##############################
#==============================================================================

#==============================================================================
# GENERAL SETTINGS
#==============================================================================

#FortDir = None          # needs to be initialized for if FortDir: ...
AbaqusCmd = ''           # if not specified, 'abq2016' is used for Windows and '/opt/abaqus/2018/SIMULIA/Commands/abq*' for Linux

# restart an aborted simulation by using the name of the iter-file as a restart file
# does not work with evolutionary optimization algorithms
restartFile = None

# specify the amount of CPUs to use for FE-simulations
NumCPUs = 2


#==============================================================================
# DEBUGGING SETTINGS
#==============================================================================

noAbaqus = False # set to true in order to not compute the abaqus job
oneIter = False # if only one iteration should be accomplished
validationFlag = False # set to true when using FE-data as input for validation pruposes


#==============================================================================
# SETTINGS FOR THE PARAMETER IDENTIFICATION
#==============================================================================

# optimisation algorithm
OptimAlg = 1 # 1: simplex, 2: steepest decent, 3: differential evolution
# define the maximum number of iterations to be performed
IterNum = 1000
# perturbation parameters for gradient based optimization methods
eps = 1e-6

# parameters to be optimized
ParString = ['A','eps0','n']
# initial guess
x0 = [ 800.0, 0.001, 0.1]
# bounds must be any positive number, need at minimum as many tuples as parameters (there can be more bounds than parameters)
bnds = ((0., 999999.0),(0., 999999.0),(0., 999999.0),(0., 999999.0),(0., 999999.0),(0., 999999.0),(0., 999999.0),(0., 999999.0))

# specify the type of data to be used
objectiveField = 'strain' # possible options are 'displacement' or 'strain'

# specify the weighting for load and displacement field contributions to the objective function
weighting = {'LDC': 0.5, 'strain': 0.5} #, 'displacement': 0.5} # weigthing factors (only the relative size matters and not the absolute values)
normFlag = False # divide DIC data by frame displacement and LDC data by frame force (this helps to better interpret the objective function value)

# further settings
firstRunNorm = False # activate to normalize objective function values by values of first run (this is useful if the minima of the objective function w.r.t. each parameter haven't been found)

# global criterion method vector (values have been identified by optimising w.r.t. to each quantity separately)
idealVector = {'LDC': 907.302196,'strain': 1.879976772186068e-05}


#==============================================================================
# SETTINGS REGARDING THE FE-SIMULATION
#==============================================================================

DataPrec = ''            # set to "double", if double precision is used within Abaqus

# specify your list of abaqus job-names, the postprocessing script (if different from standard) and the required node sets
JobNames = ['r15'] #,'r5'] # uncomment this if two experiments should be used
PostProcessingScript = 'getAbqDisp' # post processing script (without file extension)
NeckNodeSet = 'NSETDISPFIELD' # name of node set in Abaqus to be analyzed, None for all
ForceSet = 'NSETFORCE' # name of the node set in Abaqus for the force
NeckElemSet = 'ESETSTRAIN' # element set that should likely corrsepond to the node set in NeckNodeSet

# name your user material and whether the material parameters are specified in the input file or user material - automatic replacement by the algorithm with the current node set
UserMat = None # e.g. 'C02_UMAT_clean_v1d'
LocMatParInp = True # set to true if material parameters are defined in the input file - False if in the UMAT
YieldLimitCorrection = True # set to true if yield limit needs to be corrected within the UMAT
AbqMode = '3D' # let the program know if the Abaqus test will be 2D or 3D

autodelInnerNoodes = True # automatically delete FE-nodes that are not on the top surface area

# symmetry assumptions
ForceScale          = 4
DispScale           = 1


#==============================================================================
# SETTINGS FOR EXPERIMENTAL DATA
#==============================================================================

# name of the folder the experimental data is in
ExpFolder = 'Experimental Results'

# specify the names of your DIC frame data, the number of lines of the header and the quadrant of the DIC-data to use
FrameLabelsDIC   = ['r15_t2p5_00'] #,'r5_t2p5_6'] # labels to identify Aramis frame data. Need to be in same order as jobnames!
DICheaderLines = 4 # number of header lines in DIC-files
quadrantDIC = {'r15': 3,'r5': 3} # specify the qudrant of the DIC-data to use.

# specify your files containing the load-displacement-data
LabelsLDC = ['LDC_ref.txt'] #,'LDC_ref.txt'] # uncomment this if two experiments should be used
LDCFileDelimiter = ',' # delimiter between the columns in the LDC-file, e.g. ',' or '\t'

# specify the loading directions for the FE and DIC-data
LoadDirFE = 'y'
LoadDirDIC = 'y'

# displacement interval considered in LDC-objective-function (used to neglect elastic regime)
x_start = -np.inf
x_end   =  np.inf

# domain to compare DIC and FE-displacements
xDIC_min = -np.inf
xDIC_max =  np.inf
yDIC_min = -np.inf
yDIC_max =  np.inf

# List of the properties as a string expected from the DIC input
if validationFlag:
    DICPropList = ['x','y','z','ux','uy','uz']
elif LoadDirFE == LoadDirDIC:
    DICPropList = ['ID','x','y','z','ux','uy','uz','eps11','eps22']
else:
    DICPropList = ['ID','y','x','z','ux','uy','uz','eps11','eps22']


#==============================================================================
# INTERPOLATION SETTINGS
#==============================================================================

# interpolation data
EquidistantSamples = False
NumPoints = 50
DynamicWeight = True
# interpolation once at the beginning
interonce = True


#==============================================================================
# OUTPUT SETTINGS
#==============================================================================

# colored terminal output - should only sometimes (e.g. some clusters) lead to problems
colorterm = True

# disable all graphical output - only sometimes necessary (e.g. some clusters)
DisableAllGraphs = False

# create distance plots for the interpolation
CreateInterpDistPlots = False

# select frequency of output
OutputFlag = 1 # 0 is for minimal output, 1 gives developer output, 2 enables command lines
LogOutputFlag = 1 # same as above but no option for 2


#==============================================================================
############################# END OF USER INPUT ###############################
#==============================================================================


#==============================================================================
# INITIALISATION
#==============================================================================

RawDir = 'RawFrameData'# place where raw data for each frame is temporarily saved

# variable to determine first run
firstRun = True

if colorterm: from colorama import Fore, Back, Style #colored terminal outputs

# create log-file
RunName = os.path.split(os.getcwd())[-1]
if os.path.exists(RunName + '_log.txt'):
    os.remove(RunName + '_log.txt')
log = open(RunName+'_log.txt', 'w')
log.write(str(time.strftime("%c")) + ':  # O P T I M I Z A T I O N   S T A R T E D\n')
StartTime = time.time()
log.close()

# create iteration-file
if os.path.exists(RunName + '_iter.txt'):
    os.remove(RunName + '_iter.txt')
Iter = open(RunName+'_iter.txt', 'w')
Iter.write('# Starting vector: ' + str(x0) + '\n')
Iter.write('# Weighting: ' + str(weighting) + '\tnorm to 1st run value: ' + str(firstRunNorm) + '\n')
Iter.write('# Iteration \t[Parameter vector]\tObjective function value')
for job in JobNames:
    Iter.write('\tf_LDC_'+job)
for job in JobNames:
    Iter.write('\tf_DIC_'+job)
Iter.write('\tCPU time [hhhh:mm:ss]\n')
Iter.close()

# check operating platform
Platform = platform.system()
if Platform == 'Windows':
    WinPlat = True
elif Platform == 'Linux':
    WinPlat = False
else:
    SendError('Operating system not properly detected or not compatible with routine. Use Windows or Linux.')

# set commands for calling Abaqus
if not AbaqusCmd:
    if WinPlat:
        AbaqusCmd = 'abq2016'       # version 2016
    else:
        AbaqusCmd = '/opt/abaqus/2018/SIMULIA/Commands/abq*'

# depending on the computer used, fortran needs to be called explicitly
if WinPlat:
    FortDir = '\"C:\Program Files (x86)\Intel\Composer XE 2013 SP1'+os.path.sep+'bin\ifortvars_intel64.bat\"'
else:
    FortDir = None
if FortDir:
    FortDir = FortDir + ' && '
else:
    FortDir = ''

# set up user material
if UserMat:
    if WinPlat:
        UserMatABQ = ' user=' + UserMat + '.for'
    else:
        UserMatABQ = ' user=' + UserMat + '.f'
else:
    UserMatABQ = ''

# name of the file containing the material parameters to be identified
if LocMatParInp:
    ParameterFileName = JobNames[0]
    ParameterFileNameExt = '.inp'
else:
    ParameterFileName = UserMat
    if WinPlat:
        ParameterFileNameExt = '.for'
    else:
        ParameterFileNameExt = '.f'

# set up everything for storage of temporary files
cwd = os.getcwd()               # current work directory
TempFolder = os.path.sep+'temp' # Folder for temporary files
TempDir    = cwd + TempFolder   # Directory to TempFolder
ResDir = cwd+os.path.sep+'Optimization_Results' # Directory to optimization results

# delete old job files and corresponding results while avoiding deletion of other job results
for idel,ijob in enumerate(JobNames):
    cJDir = ResDir + os.path.sep + ijob
    icount = 0
    if os.path.exists(cJDir):
        while os.path.exists(cJDir):
            icount += 1
            sicount = str(icount)
            if os.path.exists(cJDir):
                cJDir = ResDir+os.path.sep+ijob+os.path.sep+sicount
                GraphDir = cJDir+os.path.sep+'FrameGraphs'
                LdcDir = cJDir+os.path.sep+'LDCplots'
    else:
        cJDir = cJDir+os.path.sep+str(icount)
        GraphDir = cJDir+os.path.sep+'FrameGraphs'
        LdcDir = cJDir+os.path.sep+'LDCplots'
    ResDir = cJDir
if not noAbaqus:
    if os.path.exists(TempDir):
        for idel,ijob in enumerate(JobNames):
            cJDir = TempDir + os.path.sep + ijob
            if os.path.exists(cJDir):
                shutil.rmtree(cJDir)

DispScaleX = DispScale # scale displacements in x-direction
DispScaleY = DispScale # scale displacements in y-direction

# for postprocessing purposes
SaveInitGuessData = 0
SaveLastIterationData = 0

# perturbation parameters for gradient based optimization methods
epsilon = np.asarray(x0)*eps

# helper to determine dicts to save obj. fun. value of 1st iteration
fDIC0, fLDC0 = {}, {}

terminalWidth = shutil.get_terminal_size().columns
print('\n\n')
message = 'S T A R T I N G   O P T I M I Z A T I O N   R O U T I N E\n'
print(message.center(terminalWidth))


#==============================================================================
# FUNCTION DEFINITIONS
#==============================================================================

def AdjustData(Data1,Data2,EquidistantSamples=False,NumPoints=100):
    '''
    Rearrange Data1 and Data2, such that both have the same format for the 
    evaluation of the objective function, the Load-Displacement-Curve needs to 
    be evaluated at the same sample points. Therefore the data with higher 
    resolution (more sample points) is interpolated on the sample points with 
    lower resolution.
    
    Info in case of errors: Scale factors or the adaption of the 
    PostProcessingScript might be incorrect.
    
    Parameters
    ----------
    Data1 : numpy.array
        Numerical data
    Data2 : numpy.array
        Experimental data
    EquidistantSamples : Boolean
        True: data samples with equivalent distance are generated
        False: original and probably non-equidistant samples are used,
        resulting in non-uniform weights        
    NumPoints : integer
        Number of points, default set to 100
        
    Returns
    -------
    [x1,y1] : numpy.array
        Numerical input data that has been adjusted and interpolated to the 
        same data points than the experimental data.
    [x2,y2] : numpy.array
        Experimental input data that has been adjusted and interpolated to the 
        same data points than the numerical data.
    '''

    # rearrange data
    x1 = abs(Data1[0])  # simulation
    y1 = abs(Data1[1])  # simulation
    x2 = abs(Data2[0])  # experiment
    y2 = abs(Data2[1])  # experiment

    # if simulation data shows higher maximum displacement than experimental data
    if x1[-1] > x2[-1]:
        if len(fvec) == 0:
            message = 'The simulation shows larger displacements than the experiment. Hence, the displacements exceeding the experimental values are ignored.'
            SendWarning(message,True)
        searching = True
        for eCt,entry in enumerate(x1):
            #find first entry of sim data that is higher than max exp displacement
            if entry > x2[-1] and searching:
                # if simulation is finer (interpolation from simulation to experiment)
                if len(x1[:eCt]) > len(x2[:eCt]): 
                    x1, y1 = x1[:eCt+1], y1[:eCt+1]
                else: 
                    x1, y1 = x1[:eCt], y1[:eCt]
                searching = False
    
    # if experimental data shows higher maximum displacement than simulation data
    elif x2[-1] > x1[-1]:
        
        if len(fvec) == 0:
            message = 'The simulation shows lower displacements than the experiment. Hence, the experimental displacements exceeding the simulation values are ignored.'
            SendWarning(message,True)
        searching = True
        for eCt,entry in enumerate(x2):
            if entry > x1[-1] and searching:
                # if simulation is finer (interpolation from simulation to experiment)
                if len(x2[:eCt]) > len(x1[:eCt]): 
                    x2, y2 = x2[:eCt+1], y2[:eCt+1]
                else: 
                    x2, y2 = x2[:eCt], y2[:eCt]
                searching = False
                
    # define interpolation function
    try:
        f1 = interpolate.interp1d(x1,y1,fill_value='extrapolate')
        f2 = interpolate.interp1d(x2,y2,fill_value='extrapolate')
    except ValueError:
        if DisableAllGraphs == False:
            plt.figure(figsize=(8,0.5))
            plt.title('Sample points plot')
            plt.vlines(Data1[0],0,0.66,'red',alpha=0.6,linewidth=1,label='Simulation sample points(Data1)')
            plt.vlines(Data2[0],0.33,1,'blue',alpha=0.6,linewidth=1,label='Experimental sample points(Data2)')
            plt.gca().axes.get_yaxis().set_visible(False)
            plt.xlabel('displacement in mm')
            plt.legend(loc=9, bbox_to_anchor=(0.5, -1), ncol=2)
        
        message = 'ERROR: LDC data interpolation failed! Check sample points plot. Each line represents a sample point. Sample points should overlap in the region of interest.'        
        SendError(message)
    
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
      
    # interpolate data
    y1 = f1(x1)
    y2 = f2(x2)

    SaveLDC_Optimum([x1,y1])
    return [x1,y1],[x2,y2]


def Call_abaqus_single(job):
    '''
    This function starts and evaluates a single Abaqus job. First, all 
    associated files like the ODB and DAT files are deleted. Subsequently, 
    the Abaqus job is submitted using the batch or shell script. Here, the 
    interactive keyword guarantees, that the computation is finished. Finally, 
    when the job has finished and no error occured, the output data base is 
    evaluated using the postprocessing script.
    Using Windows, the job is started as a subprocess to avoid the terminal 
    popping up regularly with each job submission.
    
    Info: The job cannot be terminated upon user interrupt of the script.
    
    Parameters
    ----------
    job : string
        name of the Abaqus input file
    '''    
    
    # command = AbaqusCmd + ' cae noGUI=' + PostProcessingScript # command to call post processing script
    
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
    
    # submit Abaqus job
    if WinPlat:
        startupinfo = subprocess.STARTUPINFO()
        if OutputFlag <= 2:
            startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
        ErrorCode = subprocess.call('run_'+job+'.bat', startupinfo=startupinfo, shell=False)
    else:
        ErrorCode = os.system('./run_'+job+'.sh > /dev/null')

    if ErrorCode == 0:
        error_bool = False
        if LogOutputFlag >= 0: UpdateLogfile('Abaqus job has finished. Call post processing script ...')
        if OutputFlag >= 0: print('Abaqus job has finished. Call post processing script ...')
        
        # call post processing script
        command = AbaqusCmd + ' cae noGUI=' + PostProcessingScript
        if WinPlat:
            startupinfo = subprocess.STARTUPINFO()
            if OutputFlag <= 2:
                startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
            PPerror = subprocess.call(command, startupinfo=startupinfo, shell=True)
        else:
            PPerror = os.system(command + ' > /dev/null')
        if PPerror != 0:
            error_bool = True
            if OutputFlag >= 0: print('Potential error in the post processing script')
            if LogOutputFlag >= 0: UpdateLogfile('Potential error in the post processing script')
        elif PPerror == 0:
            if OutputFlag > 0: print('Post processing script has finished')
            if LogOutputFlag >= 0: UpdateLogfile('Post processing script has finished')
        
        
    else:
        error_bool = True
        print('WARNING: Abaqus job did not finish properly.\n')
        UpdateLogfile('WARNING: Abaqus job did not finish properly')
        

    # change back to temp directory
    os.chdir('..')        
     
    return error_bool
        
def Call_abaqus_multi(JobList):
    '''
    This function starts and evaluates all specified Abaqus jobs
    simulateneously to ensure maximum efficiency due to parallelisation. First,
    all associated files like the ODB and DAT files are deleted. Subsequently, 
    the Abaqus job is submitted using the batch or shell script. Using multiple
    jobs simulateneously, the interactive keyword cannot be used and the status
    of the jobs is monitored by reading the status files. Finally, when the job
    has finished and no error occured, the output data base is evaluated using
    the postprocessing script.    
    Using Windows, the job is started as a subprocess to avoid the terminal 
    popping up regularly with each job submission.

    Parameters
    ----------
    job : list of strings
        names of all Abaqus input files
    '''   
    
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
        if WinPlat:
            startupinfo = subprocess.STARTUPINFO()
            if OutputFlag <= 2:
                startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
            subprocess.call('run_'+job+'.bat', startupinfo=startupinfo, shell=False)
        else:
            os.system('./run_'+job+'.sh > /dev/null')
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
                            if WinPlat:
                                startupinfo = subprocess.STARTUPINFO()
                                if OutputFlag <= 2:
                                    startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
                                PPerror = subprocess.call(command, startupinfo=startupinfo, shell=True)
                                
                            else:
                                PPerror = os.system(command + ' > /dev/null')
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
                            SendWarning(message)
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

def CopyFilesToTempDir(FileList,labelDIC,labelLDC,destination):
    '''
    All specified files are copied to the temporary directory to ensure that 
    the original data is unchanged while allowing to keep track of all files 
    that are being used.

    Parameters
    ----------
    FileList : list of strings

    labelDIC : string
    
    labelLDC : string
    ''' 
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
                numStr = Int2strPad(fcount,padLen)
                fcount = fcount+1
                RawFile = destination + '_DICFrame' + numStr + '.txt'
                shutil.copy(os.path.join(ExpFolder,frame),os.path.join(TempDir,destination,RawDir,RawFile))
            # copy LDC file
            elif labelLDC in frame:
                shutil.copy(os.path.join(ExpFolder,frame),os.path.join(TempDir,destination,RawDir,frame))

def Distance2D(point1,point2):
    '''
    Compute the distance between the points point1 and point2.

    Parameters
    ----------
    point1 : list
        Contains all coordinates of point1
    point2 : list
        Contains all coordinates of point2
        
    Returns
    -------
    dist : double
        Scalard valued distance
    ''' 
    dist = 0
    for i in range(2):
        dist = dist + np.square(point1[i] - point2[i])
    dist = sqrt(dist)
    return dist


def GetDirs(JobNames):
    '''
    Get names of directories containing Abaqus input files.

    Parameters
    ----------
    JobNames : list of strings
        Contains all jobnames
        
    Returns
    -------
    dirNames : list of strings
        List of all directory names
    ''' 

    dirNames = {}
    for root, dirs, files in os.walk("."):
        for name in files:
            for job in JobNames:
                if job in name and not 'temp' in root and not 'Experimental Results' in root and not 'Optimization_Results' in root:
                    dirNames[job] = (os.path.split(root)[-1])
    return dirNames


def GetObjectiveFunctionValueDIC(DICdata,FEdata,normFlag=False):
    '''
    Get objective function value f for the DIC/field data that is computed 
    based on a mean squared error function normalised by the number of data 
    points.

    Parameters
    ----------
    DICdata : np.array
        Contains all experimental data points
    FEdata : np.array
        Contains all numerical data points
    normFlag : Boolean
        True: divide DIC-data by extensometer displacement
        False: no weights are being used
        
    Returns
    -------
    f : float
        Value of the objective function
    ''' 
    
    # divide DIC-data by frame displacement
    if normFlag:
        Weight = np.copy(DICdata)
        for i in range(len(Weight)):
            if abs(DICdata[i]) > 0.001:
                Weight[i] = abs(1/Weight[i])
            else:
                Weight[i] = 0
    else:
        Weight=np.ones(len(DICdata))
    
    # check assumptions about length
    tempSum = [abs(x) - abs(y) for x,y in zip(DICdata,FEdata)]
    
    # number of data points that are used in objective function / dont have zero weight
    # exeption frame 1: weight of all datapoints is zero --> set ndp to 1 to avoid devision by zero
    ndp = 1
    if np.count_nonzero(Weight)>0: ndp = np.count_nonzero(Weight)
    
    f = sum((tempSum*Weight)**2)/ndp
    
    return f

def GetObjectiveFunctionValueLDC(Data1,Data2,DynamicWeight=True,Weight=None):
    '''
    Get objective function value f for the load-displacement-data that is 
    computed based on a mean squared error function normalised by the number of
    data points. If DynamicWeight is true, the weights for all points are 
    computed based on its distance to their neighbours to prevent intrinsic 
    biases for non-equidistant sample point distribution. If weights are 
    specified, they are multiplied with the associated values of the data 
    points.
    
    In case of an error, the objective function value is returned as 100000000.

    Parameters
    ----------
    Data1 : np.array
        Contains all numerical data points
    Data2 : np.array
        Contains all experimental data points
    DynamicWeight : Boolean (optional)
        True: use weights based on distance to each point's neighbours
        False: no weights are being used
    Weight : np.array (optional)
        
    Returns
    -------
    f : float
        Value of the objective function
    ''' 

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
                SendWarning(message)
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
        f = np.sum(((abs(Data1[1][:len(Weight)])-abs(Data2[1][:len(Weight)]))*Weight*DWeight)**2)/np.count_nonzero(Weight)
    else:
        f = 100000000

    return f


def Int2strPad(integer,padLength=3):
    '''
    Rename frame number in order to ensure correct ordering by beginning with 
    x leading zeros. Additionally, the frame number is returned as string.
    
    Parameters
    ----------
    integer: integer
        Current frame number
    padLength: integer
        Number of digits, default set to 3
    
    Returns
    -------
    outStr: string
        Properly renamed frame number
    '''
    outStr = ''
    if integer >= 10**padLength: print('Padding length may be insufficient to ensure ordering')
    if integer == 0: outStr = str('0'*padLength)
    else:
        for i in reversed(range(padLength+1)):
            if integer < 10**i and integer >= 10**(i-1):
               outStr = str('0'*(padLength-i) + str(integer))
    return outStr


def Interpolate2DSpace (mainList,interList,objectiveField): # normally FE is main and DIC inter
    '''
    Two-dimensional spatial interpolation according to the ADAPT publication.
    To minimise the interpolation error, the support points for the 
    interpolation are chosen such that they comprise a trigangle surrounding 
    the query point with minimum distance to the vertices. This is done with a 
    nearest neighbour search. The ray casting algorithm ensures that the query 
    point lies within the triangle formed by the support points.
    
    In case of an error, the objective function value is returned as 100000000.

    Parameters
    ----------
    mainList : list
        List of points that is being interpolated
    interList : list
        List of points that is being interpolated to
    objectiveField : Boolean (optional)
        True: use weights based on distance to each point's neighbours
        False: no weights are being used
        
    Returns
    -------
    errorFlag : Boolean
        True: Error within the interpolation
        False: No error within the interpolation
    d1List : list
        X-coordinate of the interpolated data set
    d2List : list
        X-coordinate of the interpolated data set
    d3List : list
        X-coordinate of the interpolated data set
    ''' 
                       
    NumNeigh = 7 # number of neighbors
    Radius = 10 # radius for the nearest neighbor search
    neigh = NearestNeighbors(n_neighbors=NumNeigh,radius=Radius)
    # extract coordinates from input
    inter_data = [i[0:2] for i in interList]
    neigh.fit(inter_data)
    main_data = [i[0:2] for i in mainList]
    neighbors = neigh.kneighbors(main_data, NumNeigh, return_distance=False)
#    print("--- %s seconds ---" % (time.time() - start_time))
    if objectiveField == 'displacement' or objectiveField == 'strain':
        d1List, d2List, d3List, sumDist = [], [], [], []
    polygon = [([],[]),([],[]),([],[])]
    x1,y1 = [0,0,0],[0,0,0]
    errorFlag = False
    
    # assemble combinations for choosing neighbors
    m=0
    combinations = np.zeros((int(special.binom(NumNeigh-1,2)),3))
    for k in range(1,NumNeigh):
        for l in range(k+1,NumNeigh):
            combinations[m] = [0,k,l]
            m = m+1
    combinations = combinations.astype(int)
    
    # loop through the list to be interpolated to
    for i,item in enumerate(mainList):
        if objectiveField == 'displacement':
            xneigh,yneigh,Dux,Duy,Duz = [],[],[],[],[]
        elif objectiveField == 'strain':
            xneigh,yneigh,Deps11,Deps22 = [],[],[],[]

        # determine where the NumNeigh closest points in each list are
        interInd = neighbors[i]
        # define the variables needed for the interpolation equation
        x,y = item[0],item[1]
        
        if validationFlag:
            # nearest neighbor interpolation
            if objectiveField == 'displacement':
                projCalcx = interList[interInd[0]][-3]
                projCalcy = interList[interInd[0]][-2]
                projCalcz = interList[interInd[0]][-1]
            elif objectiveField == 'strain':
                projCalcx = interList[interInd[0]][-2]
                projCalcy = interList[interInd[0]][-1]
            distHold = Distance2D((interList[interInd[0]][0],interList[interInd[0]][1]),(x,y))
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

                    # choose new neighbors
                    itr = itr + 1
                    if itr > len(combinations)-1:
                        errorFlag = True
                        break
                    activeNeigh = combinations[itr]
                else:
                    isInside = True
                    
            # save the total distance of the chosen neighbors
            # and get displacements and coordinates of chosen neighbors
            distHold = 0.0
            for j,n in enumerate(activeNeigh):
                distHold = distHold + Distance2D((xneigh[n],yneigh[n]),(x,y))  
                x1[j] = xneigh[n]
                y1[j] = yneigh[n]
                if objectiveField == 'displacement':
                    Dux.append(interList[interInd[n]][-3])
                    Duy.append(interList[interInd[n]][-2])
                    Duz.append(interList[interInd[n]][-1])
                elif objectiveField == 'strain':
                    Deps11.append(interList[interInd[n]][-2])
                    Deps22.append(interList[interInd[n]][-1])
            
            # calcualte u and v for both lists
            uTop = (y1[2]-y1[0])*(x-x1[0]) - (x1[2]-x1[0])*(y-y1[0])
            vTop = (y1[0]-y1[1])*(x-x1[0]) + (x1[1]-x1[0])*(y-y1[0])
            bot = (y1[2]-y1[0])*(x1[1]-x1[0]) - (x1[2]-x1[0])*(y1[1]-y1[0])
    
            # added to prevent random zero occurences
            if bot == 0:
                bot = bot + 1e-10
                message = 'Denominator of interpolation equation evaluated to zero'
                SendWarning(message)
            u = uTop/bot
            v = vTop/bot
            
            # apply U and V to interpolate one point to the other
            if objectiveField == 'displacement':
                projCalcx = (1-u-v)*Dux[0] + u*Dux[1] + v*Dux[2]
                projCalcy = (1-u-v)*Duy[0] + u*Duy[1] + v*Duy[2]
                projCalcz = (1-u-v)*Duz[0] + u*Duz[1] + v*Duz[2]
            elif objectiveField == 'strain':
                projCalcx = (1-u-v)*Deps11[0] + u*Deps11[1] + v*Deps11[2]
                projCalcy = (1-u-v)*Deps22[0] + u*Deps22[1] + v*Deps22[2]

        # append values of the current point to list
        sumDist.append(distHold)
        if objectiveField == 'displacement':
            d1List.append(projCalcx)
            d2List.append(projCalcy)
            d3List.append(projCalcz)
        elif objectiveField == 'strain':
            d1List.append(projCalcx)
            d2List.append(projCalcy)

    if objectiveField == 'displacement':
        return errorFlag, d1List, d2List, d3List, sumDist
    elif objectiveField == 'strain':
        return errorFlag, d1List, d2List, sumDist
    
def Kill_all_abaqus(JobList):
    '''
    All specified Abaqus jobs are being terminated.

    Parameters
    ----------
    JobList : list
        List of all jobs to be terminated
    ''' 
    for job in JobList:
        os.chdir(TempDir + os.path.sep + job)
        UpdateLogfile('Send terminate to Abaqus job '+job+' ...')
        print('Send terminate to Abaqus job '+job+' ...')
        
        if WinPlat == True:
            startupinfo = subprocess.STARTUPINFO()
            if OutputFlag <= 2:
                startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
            subprocess.call('kill_'+job+'.bat', startupinfo=startupinfo, shell=False)
        else:
            os.system('./kill_'+job+'.sh')


def Line2listDIC(line):
    '''
    Restructure the DIC data.

    Parameters
    ----------
    line : string
    ''' 
    out = []
    line = line.strip('\n')
    if validationFlag:
        change = line.split('\t')
    else:
        change = line.split(',')
    for item in change:
        out.append(item.strip())
    return out


def Line2listFE(line):
    '''
    Restructure the FE data.

    Parameters
    ----------
    line : string
    ''' 
    line = line.strip('\n')
    change = line.split('\t')
    try:
        change.remove('')
    except ValueError:
        pass
    return change


def MultiMinIndex(inList, numInd = 2):
    '''
    Get the closest frames of the experimental data associated to the 
    current frame required for the interpolation.
    
    Parameters
    ----------
    inList: list
        List of frames
    numInd: integer
        Number of closest frames, default set to 2
    
    Returns
    -------
    minList: list
        List of the associated frames
    '''
    
    temporaryList=inList[:]# avoid changing the input list
    minList = []
    for i in range(numInd):
        tempMin = temporaryList.index(min(temporaryList))
        minList.append(tempMin)
        temporaryList[tempMin] = 99999999 # arbitrarily large number, tried with nan but caused errors
    return minList


def Optimfunc(x,*bounds):
    '''
    Computes the objective function for a specific parameter set x. This 
    function is evaluated within the optimisation algorithm. Here, the material
    parameters are being set up and the simulation is computed. After 
    extraction and evaluation of the output data, the experimental and 
    numerical data are being processed and interpolated onto one another. This 
    is used to compute the objective function value based on a mean squared 
    error function.
    
    Parameters
    ----------
    x : list
        List containing all parameters
        
    Returns
    -------
    f : float
        Objective function value
    '''     
    
    global fvec
    global DICPropList
    global LDC_ref
    global weighting
    global fLDC0
    global fDIC0
    global firstRun
    global CreateInterpDistPlots
    global feIn
    global dicIn
    global Weight
    global uxOutlist
    global uyOutlist
    global uzOutlist
    global eps11Outlist
    global eps22Outlist
    global idList
    global inter_err_bool
    
    # initialize function value
    f = 0
    
    # dict, containing the overall objective function value for each job
    f_dict= {}
    # dict, containing the objective function value for each QOI and each job
    fQOI = defaultdict(dict)

    # initialize storage for frame data
    if interonce and firstRun: uxOutlist, uyOutlist, uzOutlist, eps11Outlist, eps22Outlist = [], [], [], [], []
    if interonce == True and firstRun == False:
        feIn,feInElem = {},{}
        for job in JobNames:
            feIn[job],feInElem[job] = [],[]
    else:
        feIn,dicIn,feInElem = {},{},{}
        for job in JobNames:
            feIn[job],dicIn[job],feInElem[job] = [],[],[]    
    
    # initialize error boolean
    err_bool = False
      
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
        if UserMat == None:
            SetUpFlowCurve(x)
        elif YieldLimitCorrection == True: # if yield limits needs to be corrected within the UMAT
            sigma_y0 = x[0] * x[1] ** x[2]
            ParVec = []
            for e in x:
                ParVec.append(e)
            ParVec.append(sigma_y0)
            SetUpParameters(ParVec,ParString,ParameterFileName,ParameterFileNameExt)
        else:
            SetUpParameters(x,ParString,ParameterFileName,ParameterFileNameExt)
        
        # penalty, if bounds are violated
        for i,e in enumerate(x):
            if e < bnds[i][0] or e > bnds[i][1]:
                SendWarning('Bounds violated! Objective function is penalized and next iteration started...')
                # do not evaluate abaqus model if bounds are violated
                err_bool = True
        
        # start simulation
        if not err_bool:
            if noAbaqus:
                err_bool = False
                if len(fvec) == 0:
                    message = 'Abaqus is not called! This is a debugging run. Switch noAbaqus to false for normal run.'
                    SendWarning(message,True)
            else:
                if len(JobNames)>1:
                    err_bool = Call_abaqus_multi(JobNames)
                elif len(JobNames)==1:
                    err_bool = Call_abaqus_single(JobNames[0])
#                    err_bool = False
                    
            
        # computation of the error function only if all simulations were finished properly
        # and bounds are not violated
        if (err_bool == False):
            for job in JobNames:
                #go to directory of current job
                os.chdir(job)
    
    #           --SUBSECTION DATA SETUP--
                # read data from all of the files according to its type
                # loop over frames/files
                for root, dirs, files in os.walk(RawDir):
                    for specfile in sorted(files):
                        if 'ABQ' in specfile and 'Node' in specfile:
                            feIn[job].append(ReadFEData(os.path.join(RawDir,specfile)))
                        elif 'ABQ' in specfile and 'Elem' in specfile:
                            feInElem[job].append(ReadFEDataElem(os.path.join(RawDir,specfile)))
                        elif 'DIC' in specfile and (interonce == False or (interonce == True and firstRun == True)):
                            dicIn[job].append(ReadDICData(os.path.join(RawDir,specfile),DICPropList,DICheaderLines))
                
                # merge nodal and element fe data to one struct
                for frameInd,frame in enumerate(feIn[job]):
                    frame.update(feInElem[job][frameInd])

                        
                if len(feIn[job])!=len(dicIn[job]) and len(fvec) == 0 and (interonce == False or (interonce == True and firstRun == True)):
                    message = 'FE and DIC data of job '+job+' have different numbers of frames.'
                    SendWarning(message,True)
                    
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
                

                if ('displacement' in weighting and weighting['displacement'] > 0 and (interonce == False or (interonce == True and firstRun == True))) or ('strain' in weighting and weighting['strain'] > 0 and (interonce == False or (interonce == True and firstRun == True))):
                    # add extensometer data to the dicIn values
                    for frCount,frame in enumerate(dicIn[job]):
                        frame['frameExten'] = LDC_ref[job][0][frCount]
     
                    # adjust sign of displacements
                    for frInd, frame in enumerate(dicIn[job]):
                        for i, elem in enumerate(frame['uy']):
                            dicIn[job][frInd]['uy'][i] = abs(elem)
                            dicIn[job][frInd]['ux'][i] = -abs(frame['ux'][i])
                
                    # choose the specified quadrant
                    if quadrantDIC[job] > 1:
                        for frInd, frame in enumerate(dicIn[job]):
                            for i, elem in enumerate(frame['x']):
                                if quadrantDIC[job] == 2 or quadrantDIC[job] == 3:
                                    #mirror x-axis
                                    dicIn[job][frInd]['x'][i] = -elem
                                if quadrantDIC[job] == 3 or quadrantDIC[job] == 4:
                                    #mirror y-axis
                                    dicIn[job][frInd]['y'][i] = -frame['y'][i]
    
                # loop over all frames
                # work on each fe frame as the basis, assumes FE coarser than DIC
                if (interonce == False or (interonce == True and firstRun == True)): print('Beginning data interpolation for job '+job)
                
                fFrame, LDC1, LDC2 = [], [], []
                numInterp = 0
                totNumInterp = len(feIn[job])*2
                for frCount, frame in enumerate(feIn[job]):
                    # read LDC stats for this frame
                    LDC1.append(float(frame['frameDisp'])*DispScaleY)
                    LDC2.append(float(frame['frameForce'])*ForceScale)
                    
                    if ('displacement' in weighting and weighting['displacement'] > 0) or ('strain' in weighting and weighting['strain'] > 0):
                        # get the 2 LDC entries closest to current FE-frame
                        if (interonce == False or (interonce == True and firstRun == True)):
                            tempList = []
                            for dicFrame in dicIn[job]:
                                tempList.append(abs(dicFrame['frameExten'] - frame['frameDisp']))
                            interpFrame = MultiMinIndex(tempList,2)
                        
                        # S P A T I A L   I N T E R P O L A T I O N
                        if 'displacement' in weighting:
                            propListFe  = ['x','y','ux','uy','uz']
                            propListDic = ['x','y','ux','uy','uz']
                            
                        elif 'strain' in weighting:
#                            propListFe  = ['IntPoint_x','IntPoint_y','IntPoint_z','strain11','strain22']
#                            propListDic = ['x','y','z','eps11','eps22']
                            propListFe  = ['IntPoint_x','IntPoint_y','strain11','strain22']
                            propListDic = ['x','y','eps11','eps22']
                        else:
                            print('Invalid objective field tag. Check line 104.')
                            
                        feDat = list(zip(*[frame[prop] for prop in propListFe]))  #the * unpacks the zip-arguemnts from list
                        if firstRun and frCount == 0: 
                            idList = []
                            for ife,fedata in enumerate(feDat):
                                if fedata[0] < xDIC_min or fedata[0] > xDIC_max or fedata[1] < yDIC_min or fedata[1] > yDIC_max: idList.append(ife)
                        for ife in idList[::-1]: del feDat[ife]
                        feDatx, feDaty, feDat1, feDat2 = [], [], [], []
                        for fedata in feDat:
                            feDatx.append(fedata[0])
                            feDaty.append(fedata[1])
                            feDat1.append(fedata[2])
                            feDat2.append(fedata[3])
                        
                        if (interonce == False or (interonce == True and firstRun == True)):
                            tmpFrame = dicIn[job][interpFrame[0]]
                            dicDat1 = list(zip(*[tmpFrame[prop] for prop in propListDic]))
                            tmpFrame = dicIn[job][interpFrame[1]]
                            dicDat2 = list(zip(*[tmpFrame[prop] for prop in propListDic]))
                            if 'displacement' in weighting:
                                inter_err_bool,uxOut1,uyOut1,uzOut1,distance1 = Interpolate2DSpace(feDat,dicDat1,'displacement')
                                numInterp = numInterp+1
                                inter_err_bool,uxOut2,uyOut2,uzOut2,distance2 = Interpolate2DSpace(feDat,dicDat2,'displacement')
                                numInterp = numInterp+1                        
                            elif 'strain' in weighting:
                                inter_err_bool,eps11Out1,eps22Out1,distance1 = Interpolate2DSpace(feDat,dicDat1,'strain')
                                numInterp = numInterp+1
                                inter_err_bool,eps11Out2,eps22Out2,distance2 = Interpolate2DSpace(feDat,dicDat2,'strain')
                                numInterp = numInterp+1  
                    
                                
                            msgNum = numInterp/totNumInterp*100
                            if OutputFlag > 0: sys.stdout.write('\r%.2f'  % msgNum + '% done with interpolation')
            
                        # geometry Matching Plots
                        if CreateInterpDistPlots and (interonce == False or (interonce == True and firstRun == True)):
                            fpl, (ax1, ax2) = plt.subplots(1,2, sharey=True)
                            ax1.scatter(dicIn[job][interpFrame[0]]['x'],dicIn[job][interpFrame[0]]['y'],c='blue',s=2)
                            ax1.scatter(frame['x'],frame['y'],c='red',s=2)
                            ax1.set_title('Interp 1')
                            ax2.scatter(dicIn[job][interpFrame[1]]['x'],dicIn[job][interpFrame[1]]['y'],c='blue',s=2)
                            ax2.scatter(frame['x'],frame['y'],c='red',s=2)
                            ax2.set_title('Interp 2')
                            fpl.suptitle('Geometry Comparison Frame '+str(frCount))
                            plt.savefig(ResDir+os.path.sep+'GeometryComparisonFrame'+str(frCount)+'.png')
                            if OutputFlag > 1: plt.show()
                            plt.close()
                        
                        if CreateInterpDistPlots and (interonce == False or (interonce == True and firstRun == True)):
                            plt.figure()
                            plt.scatter(feDatx,feDaty,c=distance1)
                            sm = plt.cm.ScalarMappable(cmap=cm.jet)#, norm=plt.Normalize(vmin=minHeat, vmax=maxHeat))
                            # Create an empty array for the colorbar
                            sm._A = []
                            plt.colorbar(sm)
                            tText = 'Distance1_Scatter_Frame'+str(frCount)
                            plt.suptitle(tText)
                            plt.savefig(ResDir+os.path.sep+tText)
                            plt.close()
                        
                        
                        if CreateInterpDistPlots and (interonce == False or (interonce == True and firstRun == True)):
                            plt.figure()
                            plt.scatter(feDatx,feDaty,c=distance2)
            
                            sm = plt.cm.ScalarMappable(cmap=cm.jet)#, norm=plt.Normalize(vmin=minHeat, vmax=maxHeat))
                            # Create an empty array for the colorbar
                            sm._A = []
                            plt.colorbar(sm)
                            tText = 'Distance2_Scatter_Frame'+str(frCount)
                            plt.suptitle(tText)
                            plt.savefig(ResDir+os.path.sep+tText)
                            plt.close()
    
    
                        #  T I M E / E X T E N S O M E T E R   W I S E   I N T E R P O L A T I O N
                        # calculate weights
                        if (interonce == False or (interonce == True and firstRun == True)):
                            T1 = (float(dicIn[job][interpFrame[1]]['frameExten']) - float(frame['frameDisp'])) / (float(dicIn[job][interpFrame[1]]['frameExten']) - float(dicIn[job][interpFrame[0]]['frameExten']))
                            T2 = (float(frame['frameDisp']) - float(dicIn[job][interpFrame[0]]['frameExten'])) / (float(dicIn[job][interpFrame[1]]['frameExten']) - float(dicIn[job][interpFrame[0]]['frameExten']))
                            # interpolate
                            if 'displacement' in weighting:
                                uxOut = [ux1*T1 + ux2*T2 for ux1,ux2 in zip(uxOut1,uxOut2)]
                                uyOut = [uy1*T1 + uy2*T2 for uy1,uy2 in zip(uyOut1,uyOut2)]
                                uzOut = [uz1*T1 + uz2*T2 for uz1,uz2 in zip(uzOut1,uzOut2)]
                                if interonce:
                                    uxOutlist.append(uxOut)
                                    uyOutlist.append(uyOut)
                                    uzOutlist.append(uzOut)
                            elif 'strain' in weighting:
                                eps11Out = [eps1*T1 + eps2*T2 for eps1,eps2 in zip(eps11Out1,eps11Out2)]
                                eps22Out = [eps1*T1 + eps2*T2 for eps1,eps2 in zip(eps22Out1,eps22Out2)]
                                if interonce:
                                    eps11Outlist.append(eps11Out)
                                    eps22Outlist.append(eps22Out)
                        else:
                            if 'displacement' in weighting:
                                uxOut = uxOutlist[frCount]
                                uyOut = uyOutlist[frCount]
                            elif 'strain' in weighting:
                                eps11Out = eps11Outlist[frCount]
                                eps22Out = eps22Outlist[frCount]
                            
                        # scatter plot for last frame
                        if 'displacement' in weighting and DisableAllGraphs == False:
                            direction = ['ux','uy']
                            if frCount == len(feIn[job])-1:
                                plt.figure(figsize=(6,3))
                                for i,item in enumerate([uxOut,uyOut]):
                                    plt.subplot(1,2,i+1)
                                    vmin = min(dicIn[job][frCount][direction[i]] + item)
                                    vmax = max(dicIn[job][frCount][direction[i]] + item)
                                    plt.scatter(dicIn[job][frCount]['x'],dicIn[job][frCount]['y'],3,dicIn[job][frCount][direction[i]],cmap='jet',vmin=vmin,vmax=vmax,alpha=1)
                                    plt.scatter(feDatx,feDaty,10,item,cmap='jet',vmin=vmin,vmax=vmax)                                            
                                    plt.title(direction[i])
                                    plt.xlabel('x-coordinate')
                                    plt.ylabel('y-coordinate')
                                    plt.axis([-0.5,max(dicIn[job][frCount]['x']),-0.5,max(dicIn[job][frCount]['y'])])
                                    plt.colorbar()
                                    # make sure colorbar and plots do not overlap:
                                    plt.tight_layout()
                        elif 'strain' in weighting and DisableAllGraphs == False:
                            direction = ['eps11','eps22']
                            if frCount == len(feIn[job])-1:
                                plt.figure(figsize=(6,3))
                                for i,item in enumerate([eps11Out,eps22Out]):
                                    if i == 0:
                                        plt.subplot(1,2,i+1)
                                        vmin = min(dicIn[job][frCount]['eps22'] + item)
                                        vmax = max(dicIn[job][frCount]['eps22'] + item)
                                        plt.scatter(feDatx,feDaty,10,eps22Out,cmap='jet')#,vmin=vmin,vmax=vmax) 
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
                                        plt.scatter(feDatx,feDaty,10,feDat2,cmap='jet',vmin=vmin,vmax=vmax)                                            
                                        plt.title('FE DATA')
                                        plt.xlabel('x-coordinate')
                                        plt.ylabel('y-coordinate')
                                        plt.axis([-0.5,max(dicIn[job][frCount]['x']),-0.5,max(dicIn[job][frCount]['y'])])
                                        plt.colorbar()
                                        # make sure colorbar and plots do not overlap:
                                        plt.tight_layout() 
                                    
                        # plot differences
                        if 'displacement' in weighting and DisableAllGraphs == False:
                            directionStr = ['ux','uy']
                            direction = [feDat1, feDat2]
                            if frCount == len(feIn[job])-1:
                                plt.figure(figsize=(6,3))
                                for i,item in enumerate([uxOut,uyOut]):
                                    plt.subplot(1,2,i+1)
                                    diff = np.array(item) - np.array(direction[i])
                                    plt.scatter(feDatx,feDaty,20,diff,cmap='jet')
                                    plt.title(directionStr[i] + ' difference')
                                    plt.xlabel('x-coordinate')
                                    plt.ylabel('y-coordinate')
                                    plt.axis('tight')
                                    plt.axis('equal')
                                    plt.colorbar()
                                    # make sure colorbar and plots do not overlap:
                                    plt.tight_layout() 
                        if 'strain' in weighting and DisableAllGraphs == False:
                            directionStr = ['strain11','strain22']
                            direction = [feDat1, feDat2]
                            if frCount == len(feIn[job])-1:
                                plt.figure(figsize=(6,3))
                                for i,item in enumerate([eps11Out,eps22Out]):
                                    plt.subplot(1,2,i+1)
                                    diff = np.array(item) - np.array(direction[i])
                                    plt.scatter(feDatx,feDaty,20,diff,cmap='jet')
                                    plt.title(directionStr[i] + ' difference')
                                    plt.xlabel('x-coordinate')
                                    plt.ylabel('y-coordinate')
                                    plt.axis('tight')
                                    plt.axis('equal')
                                    plt.colorbar()
                                    # make sure colorbar and plots do not overlap:
                                    plt.tight_layout()
                                    
                        # save plot
                        if  DisableAllGraphs == False:
                            plt.savefig(GraphDir+os.path.sep+job+'_'+objectiveField+'_iter'+str(len(fvec))+'.pdf', bbox_inches='tight')
    
                        if 'displacement' in weighting:
                            # concatenate lists of x- and y-displacements
                            ExpObjectiveInterList = np.array(uxOut+uyOut)
                            FeObjectiveList = np.array(frame['ux']+frame['uy'])
                            
                        if 'strain' in weighting:
                            # concatenate lists of x- and y-strains
                            ExpObjectiveInterList = np.array(eps11Out+eps22Out)
                            FeObjectiveList = np.array(feDat1+feDat2)
    
                        # calculate field-related objective function value for current frame
                        fu = GetObjectiveFunctionValueDIC(ExpObjectiveInterList,FeObjectiveList,normFlag)
                            
                        fFrame.append(fu)        
                                         
                        # last iteration only
                        if SaveLastIterationData == 1 or frCount > 0: # all this stuff is only done at the end
                            # save interpolated outputs as a file
                            frameNum = Int2strPad(frCount,len(str(len(feIn[job]))))
                            with open(ResDir+os.path.sep+'TempInterpDataFile_' + job + '_Frame' + frameNum + '.txt','w') as file:
                                # write general information
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
                                
                        if (interonce == False or (interonce == True and firstRun == True)):    
                            if inter_err_bool:
                                message = 'Interpolation failed: FE-node might be outside DIC-data. Check geometry matching plot!'
                                SendWarning(message)
                                break
                #End of frame loop
                
                
    
                if (interonce == False or (interonce == True and firstRun == True)): sys.stdout.write('\rInterpolation completed            \n')
        
                # results                
                # assemble FE-LDC:
                LDC = (np.array(LDC1),np.array(LDC2))                
                
                if 'displacement' in weighting and weighting['displacement'] > 0:
                    if not inter_err_bool:
                        fQOI['displacement'][job] = sum(fFrame)/np.count_nonzero(fFrame)
                    else:
                        fQOI['displacement'][job] = 1000000000
                        #Geometry Matching Plots
                        if  DisableAllGraphs == False:
                            plt.figure
                            plt.scatter(dicIn[job][interpFrame[0]]['x'],dicIn[job][interpFrame[0]]['y'],c='blue',s=2,label='DIC-data')
                            plt.scatter(frame['x'],frame['y'],c='red',s=2,label='FE-data')
                            plt.suptitle('Geometry matching plot')
                            plt.legend()
                            plt.axis('equal')
                            plt.show()
                            plt.close()
                elif 'displacemnt' not in weighting or weighting['displacement'] == 0:
                    # dummy for unused DIC objective function value
                    fQOI['displacement'][job] = -1
                  
                if 'strain' in weighting and weighting['strain'] > 0:
                    if not inter_err_bool:
                        fQOI['strain'][job] = sum(fFrame)/np.count_nonzero(fFrame)
                 
                elif 'strain' not in weighting or weighting['strain'] == 0:
                    # dummy for unused DIC objective function value
                    fQOI['strain'][job] = -1
                    
                    
                if 'LDC' in weighting and weighting['LDC'] > 0:
                    # create LDC temp file
                    with open('LoadDispCurve.txt','w') as file:
                        for entry in range(len(LDC[1])):
                            file.write(str(LDC[0][entry]) + '\t' + str(LDC[1][entry]) + '\n')                            
         
                    # remove corrupt data from LDC
                    LDC     = RemoveCorruptLDCData(LDC)
                    
                    LDC_ref_backup = LDC_ref[job]
                    
                    # interpolate LDC data if necessary
                    [LDC,LDC_ref_backup] = AdjustData(LDC,LDC_ref_backup)  
                elif 'LDC' not in weighting or weighting['LDC'] == 0:
                    # dummy for unused LDC objective function value
                    fQOI['LDC'][job] = -1
                    
                PlotOrigLDC(LDC,LDC_ref_backup,FileName=job+'_LDC_iter_'+str(len(fvec)))
                if normFlag:
                    fQOI['LDC'][job] = GetObjectiveFunctionValueLDC(LDC,LDC_ref_backup,DynamicWeight,1/LDC_ref_backup[1])
                
                else:
                    fQOI['LDC'][job] = GetObjectiveFunctionValueLDC(LDC,LDC_ref_backup,DynamicWeight)

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
                                if OutputFlag >= 0:
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
            # end of job loop
            
            f = f/len(JobNames)
            
            if firstRun:
                firstRun = False
            
            # print results to console
            for job in JobNames:
                if OutputFlag >= 0:

                    for qoi in weighting.keys():
                        print('f_' + qoi + ': %.5f' %fQOI[qoi][job])                    
                
        elif err_bool == True: # jobs did not finish properly
            f = 1000000000
            for qoi in weighting.keys():
                fQOI[qoi] = dict.fromkeys(JobNames,1000000000)


    
        if OutputFlag >= 0:
            print('\nGlobal stats:')
            print('Global objective function value = ' + str(f) + '\n')
        
        
    # append current objective function values to vector
    fvec.append(f)
    fQOIvec.append(fQOI)
    
    UpdateIterfile(x,f,fQOI,len(fvec))
    if len(fvec) > 5 and DisableAllGraphs == False:
        PlotObjFunction(fvec)
    if oneIter:
        print('\nUser-requested exit due to \'oneIter\' being true.\n')
        sys.exit()
    
    return f

   
def PlotOrigLDC(Data1,Data2,FileName='temp'):
    '''
    This function plots the original, unprocessed load displacement data and 
    saves it with the given filename.
    
        Parameters
    ----------
    Data1 : numpy.array
        Numerical data
    Data2 : numpy.array
        Experimental data
    FileName : string
        Name of the file
    '''   

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
    plt.savefig(LdcDir + os.path.sep + FileName+'.pdf', bbox_inches='tight')
    if OutputFlag >= 1: plt.show()
    plt.close()


def RayTracingMethod(x,y,poly):
    '''
    This function checks whether the query point lies within the specified 
    polygon using the ray tracing / ray casting algorithm. Returns false if 
    the query point lies on one of the edges.
    
    Parameters
    ----------
    x : float
        x-coordinate
    y : float
        y-coordinate
    poly : list
        list of x and y-coordinates of the points building the polygon
        
    Returns
    -------
    inside : boolean
        True: query point lies within the triangle
        False: query point lies within the triangle
    '''

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


def ReadDICData(filename,properties,headerlines):
    '''
    The experimental DIC data such as coordinates and displacements are being 
    read from the specified text file.
    
    Parameters
    ----------
    filename : string 
        Name of the file being read
            
    Returns
    -------
    frameSave : dict
        Dict containing all relevant data
    '''
    
    frameSave = {}
    headerErr = False
    for prop in properties: frameSave[prop] = []
    with open(filename) as file:
        for lCount, line in enumerate(file):
            if validationFlag:
                headerlines = 1

            # check header of input file for tags
            if lCount == headerlines-1:
                if 'displacement' in weighting:
                    if not any([item in line for item in ['disp','ux']]):
                        headerErr = True
                if 'strain' in weighting:
                    if not any([item in line for item in ['epsilon','strain','eps','LE']]):
                        headerErr = True
                if headerErr:
                    message = 'Current quantity-of-interest was not found in header of DIC input file! Check DIC files and number of header lines'
                    SendWarning(message)
            if lCount > headerlines-1:
                temp = Line2listDIC(line)
                for kCount,keys in enumerate(properties):
                    frameSave[keys].append(float(temp[kCount]))

    file.close()            
    return frameSave



def ReadFEData(filename):
    '''
    The relevant data from the requested nodes in the finite element simulation 
    such as nodal coordinates and displacements are being read from the 
    specified text file.
    
    Parameters
    ----------
    filename : string 
        Name of the file being read
            
    Returns
    -------
    frameSave : dict
        Dict containing all relevant data
    '''
    
    frameSave = {}
    with open(filename) as file:
        for lCount, line in enumerate(file):
            if not lCount == 0 and not lCount == 1:
                temp = Line2listFE(line)
                for kCount,keys in enumerate(enumPropList):
                    frameSave[keys].append(float(temp[kCount]))
            elif lCount == 0:
                temp = Line2listFE(line)
#                frameKeyList = ['frameTime','frameForce','frameExten','frameDisp']
                frameKeyList = ['frameTime','frameForce','frameDisp'] 
                for kCount, keys in enumerate(frameKeyList):
                    frameSave[keys] = float(temp[2*kCount+1])
            else:
                temp = Line2listFE(line)
                enumPropList = temp
                for item in temp:
                    frameSave[item]=[]
    file.close()
#    if sum([abs(y) for y in frameSave['uy']]) != 0: UNCOMMENT AND TAB NEXT LINE TO FILTER ZERO FRAMES
    return frameSave



def ReadFEDataElem(filename):
    '''
    The relevant data from the requested elements in the finite element 
    simulation such as integration point coordinates and stresses/strains are 
    being read from the specified text file.
    
    Parameters
    ----------
    filename : string 
        Name of the file being read
            
    Returns
    -------
    frameSave : dict
        Dict containing all relevant data
    '''
    
    frameSave = {}
    with open(filename) as file:
        for lCount, line in enumerate(file):
            temp = Line2listFE(line)
            if not lCount == 0 and not lCount == 1:
                for kCount,keys in enumerate(enumPropList):
                    frameSave[keys].append(float(temp[kCount]))
            elif lCount == 1:
                enumPropList = temp
                for item in temp:
                    frameSave[item]=[]
    file.close()
    return frameSave


def ReadLDCFile(filename, LDCFileDelimiter = ','):
    '''
    The load-displacement-curve is being read from the specified file.
    
    Parameters
    ----------
    filename : string 
        Name of the file being read
    LDCFileDelimiter: string
        Delimiter between the columns of the LDC-file, e.g. ',' or '\t'
            
    Returns
    -------
    frameSave : dict
        Dict containing all relevant data
    '''
    
    # this is the implementation for a .txt file
    array1 = []
    array2 = []
    with open(filename, "r") as f:
        for line in f:
            temp1 = line.split(LDCFileDelimiter)
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
    '''
    Restarts the script in case of errors or abortions by skipping the 
    computation of the FE-simulations and instead using the objective function 
    values from the iter file. This works well for the simplex algorithm as 
    well as the steepest decent method. However, this does not work for 
    random-based algorithms like differential evolution.
    
    Parameters
    ----------
    filename : string 
        Name of the iter file
    '''
    
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

        # skip header line
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

       
def RemoveCorruptLDCData(Data,DelTol=100.0):
    '''
    Remove all displacement data that is greater than the chosen tolerance.
    
    Parameters
    ----------
    Data : list 
        Load-displacement-curve
    DelTol: float
        Tolerance for the displacement data to be deleted. Default is 100
            
    Returns
    -------
    Data1, Data2 : list
        Adjusted input data
    '''

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
    '''
    Save load-displacement-curve from the initial guess for visualisation 
    purposes.
    
    Parameters
    ----------
    Data1, Data2 : list 
        Load-displacement-curve
    '''
    
    # save Load displacement curve for initial guess
    global SaveInitGuessData
    if SaveInitGuessData == 0:
        shutil.copy('LoadDispCurve.txt','LoadDispCurve_InitialGuess.txt')
        UpdateLogfile('\nData for initial guess has been saved.')
        SaveInitGuessData = 1
        # save load displacement curve for reference data
        # for interpolated data points
        shutil.copy('LoadDispCurve.txt','LoadDispCurve_ref_Interp.txt')


def SaveLDC_Optimum(Data):
    '''
    Save load-displacement-curve for the optimised parameter set for 
    visualisation purposes.
    
    Parameters
    ----------
    Data : list 
        Load-displacement-curve
    '''
    
    # save load displacement curve for optimized parameter set
    global SaveLastIterationData
    if SaveLastIterationData == 1:
        shutil.copy('LoadDispCurve.txt','LoadDispCurve_Optimum.txt')
        UpdateLogfile('\nData for optimized parameters has been saved.')


def SendError(message='An error has occured'):
    '''
    Save load-displacement-curve for the optimised parameter set for 
    visualisation purposes.
    
    Parameters
    ----------
    Data : list 
        Load-displacement-curve
    '''
    
    print('\n')
    if colorterm:
        print(Fore.RED + 'ERROR:  ' + message)
    else:
        print('ERROR:  ' + message)
    print('\n')
    UpdateLogfile('ERROR:  ' + message)
    sys.exit()


def SendWarning(message='A warning has occured',suppress = False):
    '''
    Send a warning message to the output window and the log file.
    
    Parameters
    ----------
    message : string 
        Text for the warning message
    supress : boolean (optional)
        True: warning message will be displayed only once
        False: warning message will be displayed each iteration and evaluation
    '''
    
    if colorterm:
        print(Fore.GREEN + 'WARNING:  ' + Style.RESET_ALL + message + suppress*'(This warning will be suppressed from now on)')
    else:
        print('WARNING:  ' + message + suppress*'(This warning will be suppressed from now on)')
    UpdateLogfile('WARNING:  ' + message + suppress*'(This warning will be suppressed from now on)')


def SetUpFlowCurve(x,ext='.inp'):
    '''
    Compute a flow curve according to Swift's law (A*(e0+e)^n) and write it to 
    a text file named "FlowCurve" as tabulatead data that should be linked to 
    the FE-model.
    
    Parameters
    ----------
    x : list 
        List containing three material parameters
    '''

    eps_min = 0.0 # lower bound for the strains
    eps_max = 1.0 # upper bound for the strains
    n = 1000 # number of data points
    for job in JobNames:
        os.chdir(job)
        file = open('FlowCurve'+ext,'w+')
        
        for i in range(n+1):
            eps = (eps_max-eps_min)*i/n
            sig = x[0] * (eps + x[1]) ** x[2]
            file.write(str(sig) + ',' + str(eps) + '\n')
    
        file.close()
        os.chdir('..')
    
        
def SetUpParameters(x,ParString,FileName,ext='.for'):
    '''
    Set up all parameters (x) in the specified file that is typically the 
    input or user material file by replacing the strings (ParString) with the 
    associated parameters.
    
    Parameters
    ----------
    x : list 
        List containing all material parameters
    ParString : list of strings
        List containing all placeholders that are to be replaced with the 
        parameter values
    FileName : string
        Name of the file (it has to be FileName_base)
    ext : string
        File name extension
    '''

    # if more than one parameter needs to be identified
    for job in JobNames:
        os.chdir(job)

        with open(FileName + ext,'wt') as fout:
            with open(FileName + '_base' + ext,'rt') as fin:
                for line in fin:
                    for i,par in enumerate(ParString):
                            if par in line:
                                line = line.replace(par, str(x[i]))                    
                    fout.write(line)

        fin.close()
        fout.close()

        if LogOutputFlag > 1:
            UpdateLogfile('Parameters have been set up (' + str(x) + ')')
        if OutputFlag > 1:
            sys.stdout.write('Parameters have been set up (' + str(x) + ')\n')
    
        os.chdir('..')


def SetUpPPScript(PostProcessingScript,JobName,RunName):
    '''
    Pass all variables to the postprocessing script by replacing placeholders
    with the specified strings.
    
    Parameters
    ----------
    PostProcessingScript : string 
        Name of the postprocessing-script
    JobName : string
        Name of the job
    RunName : string
        Name of the folder that is generated automatically
    '''
    
    replaceList = [['JobName_par',JobName],['RunName_par',RunName],['NeckNodeSet_par',NeckNodeSet],['ForceSet_par',ForceSet],
                   ['LoadDir_par',LoadDirFE],['DataPrec_par',DataPrec],['RawDir_par',RawDir],
                   ['NeckElemSet_par',NeckElemSet],['DispScaleX_par',DispScaleX],['DispScaleY_par',DispScaleY]]
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
    '''
    Create iteration file if not existing and then updates it by appending the
    current parameter set and each objective function value for every 
    iteration.
    
    Parameters
    ----------
    x : float 
        Parameter set
    f : float
        Global objective function value
    xfQOI : float
        Local objective function value (forces, strains, displacements)
    niter : integer
        Iteration number
    '''
    
    # create iter file if not existing already
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
    
    # update log-file: parameter vector, objective function value, CPU time
    Iter = open(cwd + os.path.sep + RunName+'_iter.txt', 'a+')
    t = (time.time()-StartTime)
    h = int(t/60/60)                # hours
    m = int(t/60 - h*60)            # minutes
    s = int(t - h*60*60 - m*60)     # seconds
    
    # write iteration number
    Iter.write(str(niter) + '\t')
    # write parameter vector by loop to avoid new line (as done by print(x))
    Iter.write('[')
    for i in x:
        Iter.write('% 10.8E' %i +' ')
    Iter.write(']')
    # write objective function values
    Iter.write('\t' + '% 10.6f' %f + '\t')
    for qoi in weighting.keys():
        for job in JobNames:
            Iter.write('% 10.8E' %xfQOI[qoi][job] +' ')
    # write timestamp
    Iter.write('\t' + '%04d' %h + ':' + '%02d' %m + ':' + '%02d' %s + '\n')
    Iter.close()
        
    
def UpdateLogfile(message):
    '''
    Write message into the log file.
    
    Parameters
    ----------
    message : string 
        Message to be written into the log file
    '''
    
    log = open(cwd + os.path.sep + RunName+'_log.txt', 'a+')
    log.write( str(time.strftime("%c")) + ':  ' + message + '\n')
    log.close()
 

def WriteObjFunctionData(f,fDIC,fLDC):
    '''
    Write objective function values (f, fDIC and fLDC) into the associated text
    file.
    
    Parameters
    ----------
    f : float 
        Global objective function value
    fDIC : float 
        Force-related objective function value 
    fLDC : float 
        DIC-data-related objective function value
    '''
    
    padLen = 4
    count = len(fvec)-1
    ObjFunFile = open(ResDir+os.path.sep+RunName+'_ObjectiveFunctionValues.txt', 'a')
    ObjFunFile.write('Iteration ' + Int2strPad(count,padLen) + ':\t' + str(f) + '\n')
    ObjFunFile.close()
    
    ObjFunFile = open(ResDir+os.path.sep+RunName+'_ObjectiveFunctionValuesDIC.txt', 'a')
    if count == 0:
        ObjFunFile.write('Iteration')
        for job in JobNames: ObjFunFile.write('\t'+job)
    line = '\n' + Int2strPad(count,padLen)
    for job in JobNames: line = line + '\t' + str(fDIC[job])
    ObjFunFile.write(line)
    ObjFunFile.close()

    ObjFunFile = open(ResDir+os.path.sep+RunName+'_ObjectiveFunctionValuesLDC.txt', 'a')
    if count == 0:
        ObjFunFile.write('Iteration')
        for job in JobNames: ObjFunFile.write('\t'+job)
    line = '\n' + Int2strPad(count,padLen)
    for job in JobNames: line = line + '\t' + str(fLDC[job])
    ObjFunFile.write(line)
    ObjFunFile.close()   
    

def PlotObjFunction(vect,log=True):
    '''
    Create figure to show convergence by illustrating objective function value 
    over the objective function calls. Only sucessful iterations are shown 
    here.
    
    Parameters
    ----------
    vect : list of floats 
        Vector containing all objective function values
    log : boolean 
        True: show logarithmic axes
        False: show linear axes
    '''

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


def WriteResultsFile(x0,xopt,ParString):
    '''
    Write initial guess as well as optimal solution to the results file. This 
    function is called after successful optmisation.
    
    Parameters
    ----------
    x0 : list of floats
        Vector containing the initial parameter set
    xopt : list of floats
        Vector containing the optimised parameter set
    ParString : list of strings
        Names of all parameters
    '''
    
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
if __name__ == '__main__':
    try:
        start_time = time.time()
        global fvec
        global LDC_ref
        global f_ini_dict
        
        # vector containing the total objective function value of all iterations
        fvec = []
        # vector containing the objective function values of specific qunatities of interest
        fQOIvec = []    
        
        # check for restart file
        if not restartFile == None:
            try:
                x0, fListRestart, fQOI_ini, weighting, firstRunNorm = ReadRestartFile(restartFile)
                message = 'Restart file has been read. Using ' + str(len(fListRestart)) + ' iterations from restart file.'
                if OutputFlag > 0: print(message)
                if LogOutputFlag > 0: UpdateLogfile(message)
                    
            except FileNotFoundError:
                message = 'Specified restart file not found! Beginning optimisation, starting from inital guess.'
                SendWarning(message)
        
        # create folder structure
        directories = (GetDirs(JobNames))
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
                SendError(message)
                
            if UserMat:
        
                if WinPlat:
                    FileList = [job+'.inp','getAbqDisp_base.py',
                        'LoadDispCurve_ref.txt',job+'.bat',UserMat+'_base.for']
                else:
                    FileList = [job+'.inp', job+'_base.inp','getAbqDisp_base.py',
                            'LoadDispCurve_ref.txt',job+'.sh',UserMat+'.f']
                
            else:
                if WinPlat:
                    FileList = [job+'.inp','getAbqDisp_base.py',
                        'LoadDispCurve_ref.txt',job+'.bat']   
                else:
                    FileList = [job+'.inp', job+'_base.inp','getAbqDisp_base.py',
                        'LoadDispCurve_ref.txt',job+'.sh']
                 
            CopyFilesToTempDir(FileList,FrameLabelsDIC[i],LabelsLDC[i],job)
            os.chdir('..')    
        
        LDC_ref = {}
        for i,job in enumerate(JobNames):
            file = os.path.join(TempDir,job,job+'.inp')
            
            # change to directory of each job in temp folder
            os.chdir(os.path.join(TempDir,job))
            
            # get LDC  data from experiments
            LDC_ref[job] = ReadLDCFile(RawDir+os.path.sep+LabelsLDC[i], LDCFileDelimiter)
    
            # set up post processing script
            SetUpPPScript(PostProcessingScript,job,RunName)
            
            # write batch files for each job
            StartJob = AbaqusCmd + ' job=' + job + UserMatABQ + ' cpus=' + str(NumCPUs)
            if WinPlat:
                batchfile_run=open('run_'+job+'.bat','w+')
                batchfile_run.write('@ECHO off\n')
                if len(JobNames) > 1:
                    batchfile_run.write('CALL ' + StartJob + ' double=both')
                elif len(JobNames) == 1:
                    batchfile_run.write('CALL ' + StartJob + ' double=both interactive ask_delete=OFF')
            else:
                batchfile_run=open(os.open('run_'+job+'.sh', os.O_WRONLY | os.O_CREAT, 0o777),"w+")
                if len(JobNames) > 1:
                    batchfile_run.write(StartJob + '\n')
                elif len(JobNames) == 1:
                    batchfile_run.write(StartJob+ ' interactive' + '\n')
            batchfile_run.close()
    
            if WinPlat:
                batchfile_kill=open('kill_'+job+'.bat','w+')
                batchfile_kill.write('@ECHO off\n')
                batchfile_kill.write('CALL ' + AbaqusCmd + ' job=' + job +' terminate')
            else:
                batchfile_kill=open(os.open('kill_'+job+'.sh', os.O_WRONLY | os.O_CREAT, 0o777),"w+")
                batchfile_kill.write(AbaqusCmd + ' job=' + job +' terminate')
            batchfile_kill.close()
            
            os.chdir('..')
        
        # work from TempDir
        os.chdir(TempDir)
        
        # perform optimisation according to the specified algorithm
        if OptimAlg == 1: # simplex algorithm
            message = 'Optimisation is done using the simplex algorithm.'    
            print(message)
            UpdateLogfile(message)
            (xopt, fopt, iters, funcalls, warnflag, allvecs) = sopt.fmin(Optimfunc,x0,args=(bnds),xtol=0.0001,ftol=0.0001,maxiter=IterNum,maxfun=None,full_output=1,disp=1,retall=1,callback=None)
        elif OptimAlg == 2: # steepest decent algorithm
            message = 'Optimisation is done using the steepest decent algorithm.'    
            print(message)
            UpdateLogfile(message)
            (xopt, fopt, funcalls, gradcalls, warnflag, allvecs) = sopt.fmin_cg(Optimfunc,x0,args=(bnds),gtol=1e-8,epsilon=epsilon,maxiter=IterNum,full_output=1,disp=1,retall=1,callback=None)
        elif OptimAlg == 3: # differential evoltuion algorithm
            message = 'Optimisation is done using the differential evolution algorithm.'    
            print(message)
            UpdateLogfile(message)
            result = sopt.differential_evolution(Optimfunc,bnds,args=(bnds),maxiter=IterNum, popsize=15, tol=0.01, mutation=(0.5, 1), recombination=0.7, seed=None, callback=None); xopt = result.x
        else:
            message = 'Optimisation algorithm has not been selected correctly.'
            print(message)
            SendWarning(message,True)
            
        UpdateLogfile('\nOptimization finished successfully.')
        
        # run simulation with optimized parameterset and save all requested outputs
        Optimfunc(xopt)
        
        WriteResultsFile(x0,xopt,ParString)
        
        # copy _log.txt and _iter.txt to corresponding result folder to save them
        shutil.copy(os.path.join(cwd,RunName+'_iter.txt'),ResDir)
        shutil.copy(os.path.join(cwd,RunName+'_log.txt'),ResDir)
        
        if DisableAllGraphs == False:
            PlotObjFunction(fvec)
        
        # change to original working directory
        os.chdir(cwd)
        
    except KeyboardInterrupt:
        message = 'Parameter identification script interrupted by user.'
        print('\n'+message+'\n')
        UpdateLogfile(message)
        Kill_all_abaqus(JobNames)
        sys.exit()
