#!/usr/bin/env python

import sys
import os
# use the following statement if mantidplotnightly is required.
sys.path.insert(0,"/opt/mantidnightly/bin")
sys.path.insert(0,"/opt/mantidnightly/lib")
import imp
import stat
import shutil as shutil

#sys.path.insert(0,"/opt/mantid37/bin")

from XMLparser import XMLparser
from mantid import *
from mantid.simpleapi import *
from MaskAngleold import *
from numpy import *
import mantid.api


def GetPathFromRunNumber(instrument,run):
    if instrument in ['ARCS','CNCS','SEQUOIA','HYSPEC']:
        #First try find_data
        finder=FileFinder.Instance()
        print(instrument, run)    
        try:
            #path = finder.findRuns(str(run))[0] - ORIGINAL
            path = finder.findRuns(instrument+str(run))[0]
            return path
        except:
            pass

        instrumentsn = config.getFacility().instrument(instrument).shortName()
        f = os.popen('findnexus --event -i %s %s' % (instrumentsn,run), 'r')
        path = f.readline().strip()

        #check if there is an "ERROR" in the output string
        if path.find('ERROR') == -1 :
            return path
        else:
            raise ValueError("Event Nexus file not found")
    else:
        raise ValueError("Instrument not yet implemented")



class dgsreduction(object):
    def __init__(self, XMLfile=None):
        self.instrument=None
        self.filterbadpulses=None
        self.vanruns=None
        self.units=None
        self.vanmin=None
        self.vanmax=None
        self.processedfilename=None
        self.maskfilename=None
        self.mask=None
        self.normalizedcalibration=None
        self.ipts=None
        self.runs=None
        self.efixed=None
        self.t0=None
        self.calce=None
        self.eimon1=None
        self.eimon2=None
        self.emin=None
        self.emax=None
        self.loadmon=True
        self.ebin=None
        self.qstep=None
        self.qmin=None
        self.qmax=None
        self.kiokf=None
        self.tibg=None
        self.tibgstart=None
        self.tibgstop=None
        self.grouping=None
        self.powderanglestep=None
        self.powderanglemin=None
        self.powderanglemax=None        
        self.goniometermotor=None
        self.goniometermotoroffset=None
        self.goniometermotoraxis=None
        self.goniometermotordirection=None
        self.save=None
        self.friendlyname=None
        self.calibrationtext=''
        self.datatext=''
        self.scantype = 'single'
        self.logvalue=''
        self.logvaluestep = None
        self.logvaluemin = None
        self.logvaluemax = None
        self.logconststep =True
        self.friendlynamelogs = None
        self.vanpath = ''
        self.datapath = ''
        self.lattice = ''
        self.ub = ''
        self.filternames = None
        self.filtermin = None
        self.filtermax = None
        self.datadate = None
        self.vanmine = None
        self.vanmaxe = None
        self.vanminangle = None
        self.vanmaxangle = None
        self.vanmass = None
        self.samplemass = None
        self.sample_rmm = None
        self.scalefactor = None


        if XMLfile!=None:
            if not(os.path.isfile(XMLfile)):
                raise IOError ("data text file "+ XMLfile+ " not found")
            self.RunFromXML(XMLfile)    


    def RunFromXML(self,filename):      
        #TODO check the xml against a schema
        parsed=XMLparser(filename)
        if parsed.root.tag !='dgsreduction':
            raise RuntimeError("This is not a dgsreduction xml file")
 
        if 'instrument' in list(parsed.calibdict.keys()):
            self.LoadInstrumentSettings(parsed.calibdict['instrument'])
        else:
            raise RuntimeError("No instrument defined in the XML file")
        self.LoadParameters(parsed.calibdict)

        if len(self.vanpath) > 1:
            config.appendDataSearchDir(self.vanpath)


        self.PerformCalibration()
        if parsed.datadicts!=[]:    
            for d in range(len(parsed.datadicts)-1):
                dgsi=dgsreduction()
                dgsi.LoadInstrumentSettings(parsed.calibdict['instrument'])
                dgsi.LoadParameters(parsed.calibdict)
                dgsi.PerformCalibration() 
                dgsi.LoadParameters(parsed.datadicts[d])
                if len(self.datapath) > 1:
                    config.appendDataSearchDir(self.datapath)
                dgsi.Execute()

            self.LoadParameters(parsed.datadicts[-1])
            self.Execute()


    #Where the data are actually reduced.
    def Execute(self):
        self.loadmon = self.calce
        self.resetEnergyToNone=False            
        if self.efixed==None:
            self.resetEnergyToNone=True
        if self.scantype == 'single':
            #load and filter bad pulses
            #load and add the runs together.
            #get the path for the first file
            path = GetPathFromRunNumber(self.instrument,self.runs[0])
            #load the file.
            Load(Filename=path,OutputWorkspace = 'data')
            print(("Datafile "+path+" loaded."))

            #Fix all of the time series log values to start at the same time as the proton_charge
            r=mtd['data'].getRun()
            for x in list(r.keys()):
                if x not in ['duration','proton_charge','start_time','run_title','run_start','run_number','gd_prtn_chrg','end_time']:
                    try:
                        ShiftTime('data',x)
                    except:
                        pass

            if self.loadmon:            
                LoadNexusMonitors(Filename=path,OutputWorkspace = 'monitorws')

            if self.filterbadpulses:
                FilterBadPulses(InputWorkspace = 'data', OutputWorkspace = 'data')
            
            #2/1/12 added code to filter by additional log values.
            if self.filternames != None:
                cntr = 0
                for part in self.filternames:
                    FilterByLogValue(InputWorkspace = 'data', OutputWorkspace = 'data', LogName=part, MinimumValue=self.filtermin[cntr],MaximumValue=self.filtermax[cntr],TimeTolerance=0,LogBoundary='Left')
                    self.datatext += "Data filtered by "+part+" between "+str(self.filtermin[cntr])+" and "+str(self.filtermax[cntr])+".\n"
                    print(("Data filtered by "+part+" between "+str(self.filtermin[cntr])+" and "+str(self.filtermax[cntr])+".\n"))
                    cntr += 1
                
            self.datatext += "Loaded data run from "+path +"\n"

            if len(self.runs) > 1:
                for i in range(1,len(self.runs)):
                    path = GetPathFromRunNumber(self.instrument,self.runs[i])
                    Load(Filename=path,OutputWorkspace = 'datatemp')
                    print(("Datafile "+path+" loaded."))

                    #Fix all of the time series log values to start at the same time as the proton_charge
                    r=mtd['datatemp'].getRun()
                    for x in list(r.keys()):
                        if x not in ['duration','proton_charge','start_time','run_title','run_start','run_number','gd_prtn_chrg','end_time']:
                            try:
                                ShiftTime('datatemp',x)
                            except:
                                pass

                    if self.loadmon:                    
                        LoadNexusMonitors(Filename=path,OutputWorkspace = 'monitortemp')
                    if self.filterbadpulses:
                        FilterBadPulses(InputWorkspace = 'datatemp', OutputWorkspace = 'datatemp')

                    #2/1/12 added code to filter by additional log values.
                    if self.filternames != None:
                     cntr = 0
                     for part in self.filternames:
                        FilterByLogValue(InputWorkspace = 'data', OutputWorkspace = 'data', LogName=part, MinimumValue=self.filtermin[cntr],MaximumValue=self.filtermax[cntr],TimeTolerance=0,LogBoundary='Left')
                        self.datatext += "Data filtered by "+part+" between "+str(self.filtermin[cntr])+" and "+str(self.filtermax[cntr])+".\n"
                        print(("Data filtered by "+part+" between "+str(self.filtermin[cntr])+" and "+str(self.filtermax[cntr])+".\n"))
                        cntr += 1

                    Plus(LHSWorkspace='data', RHSWorkspace = 'datatemp', OutputWorkspace='data')
                    if self.loadmon:                    
                        Plus(LHSWorkspace='monitorws', RHSWorkspace = 'monitortemp', OutputWorkspace='monitorws')
                    self.datatext += "Added data run from "+path +"\n"

            if self.filterbadpulses:
                self.datatext += "Bad pulses have been filterd from the data file(s).\n"

            #This is where the reduction is done.
            self.ProcessWorkspace('data')


        if self.scantype == 'step':
            #load each file and process individually, ONE summary file.
            for run in self.runs:
                #get the path for the first file
                path = GetPathFromRunNumber(self.instrument,run)
                #load the file.
                Load(Filename=path,OutputWorkspace = 'data')
                print(("Datafile "+path+" loaded."))

                #Fix all of the time series log values to start at the same time as the proton_charge
                r=mtd['data'].getRun()
                for x in list(r.keys()):
                    if x not in ['duration','proton_charge','start_time','run_title','run_start','run_number','gd_prtn_chrg','end_time']:
                        try:
                            ShiftTime('data',x)
                        except:
                            pass

                if self.loadmon:
                    LoadNexusMonitors(Filename=path,OutputWorkspace = 'monitorws')
                if self.filterbadpulses:
                    FilterBadPulses(InputWorkspace = 'data', OutputWorkspace = 'data')

                #2/1/12 added code to filter by additional log values.
                if self.filternames != None:
                    cntr = 0
                    for part in self.filternames:
                        FilterByLogValue(InputWorkspace = 'data', OutputWorkspace = 'data', LogName=part, MinimumValue=self.filtermin[cntr],MaximumValue=self.filtermax[cntr],TimeTolerance=0,LogBoundary='Left')
                        self.datatext += "Data filtered by "+part+" between "+str(self.filtermin[cntr])+" and "+str(self.filtermax[cntr])+".\n"
                        print(("Data filtered by "+part+" between "+str(self.filtermin[cntr])+" and "+str(self.filtermax[cntr])+".\n"))
                        cntr += 1

                self.datatext += "Loaded data run from "+path +"\n"
                if self.filterbadpulses:
                    self.datatext += "Bad pulses have been filtered from the data file(s).\n"

                #This is where the reduction is done.
                self.ProcessWorkspace('data')
                if self.resetEnergyToNone:
                    self.efixed=None

 
        if self.scantype == 'sweep':
            #load and filter bad pulses
            #load and add the runs together.
            #get the path for the first file
            path = GetPathFromRunNumber(self.instrument,self.runs[0])
            #load the file.
            Load(Filename=path,OutputWorkspace = 'data')
            print(("Datafile "+path+" loaded."))

            #Fix all of the time series log values to start at the same time as the proton_charge
            r=mtd['data'].getRun()
            for x in list(r.keys()):
                if x not in ['duration','proton_charge','start_time','run_title','run_start','run_number','gd_prtn_chrg','end_time']:
                    try:
                        ShiftTime('data',x)
                    except:
                        pass

            if self.loadmon:
                LoadNexusMonitors(Filename=path,OutputWorkspace = 'monitorws')
            if self.filterbadpulses:
                FilterBadPulses(InputWorkspace = 'data', OutputWorkspace = 'data')
            #2/1/12 added code to filter by additional log values.
            if self.filternames != None:
                cntr = 0
                for part in self.filternames:
                    FilterByLogValue(InputWorkspace = 'data', OutputWorkspace = 'data', LogName=part, MinimumValue=self.filtermin[cntr],MaximumValue=self.filtermax[cntr],TimeTolerance=0,LogBoundary='Left')
                    self.datatext += "Data filtered by "+part+" between "+str(self.filtermin[cntr])+" and "+str(self.filtermax[cntr])+".\n"
                    print(("Data filtered by "+part+" between "+str(self.filtermin[cntr])+" and "+str(self.filtermax[cntr])+".\n"))
                    cntr += 1

            self.datatext += "Loaded data run from "+path +"\n"

            if len(self.runs) > 1:
                for i in range(1,len(self.runs)):
                    path = GetPathFromRunNumber(self.instrument,self.runs[i])
                    Load(Filename=path,OutputWorkspace = 'datatemp')
                    print(("Datafile "+path+" loaded."))

                    #Fix all of the time series log values to start at the same time as the proton_charge
                    r=mtd['datatemp'].getRun()
                    for x in list(r.keys()):
                        if x not in ['duration','proton_charge','start_time','run_title','run_start','run_number','gd_prtn_chrg','end_time']:
                            try:
                                ShiftTime('datatemp',x)
                            except:
                                pass

                    if self.loadmon:                    
                        LoadNexusMonitors(Filename=path,OutputWorkspace = 'monitortemp')
                    if self.filterbadpulses:
                        FilterBadPulses(InputWorkspace = 'datatemp', OutputWorkspace = 'datatemp')
                    #2/1/12 added code to filter by additional log values.
                    if self.filternames != None:
                        cntr = 0
                        for part in self.filternames:
                            FilterByLogValue(InputWorkspace = 'data', OutputWorkspace = 'data', LogName=part, MinimumValue=self.filtermin[cntr],MaximumValue=self.filtermax[cntr],TimeTolerance=0,LogBoundary='Left')
                            self.datatext += "Data filtered by "+part+" between "+str(self.filtermin[cntr])+" and "+str(self.filtermax[cntr])+".\n"
                            print(("Data filtered by "+part+" between "+str(self.filtermin[cntr])+" and "+str(self.filtermax[cntr])+".\n"))
                            cntr += 1

                    Plus(LHSWorkspace='data', RHSWorkspace = 'datatemp', OutputWorkspace='data')
                    if self.loadmon:                    
                        Plus(LHSWorkspace='monitorws', RHSWorkspace = 'monitortemp', OutputWorkspace='monitorws')
                    self.datatext += "Added data run from "+path +"\n"

            if self.filterbadpulses:
                self.datatext += "Bad pulses have been filterd from the data file(s).\n"
            
            wsrun = mtd['data'].run()
            #Check if the workspace has the variable of interest
            
            if (self.logvalue == None or wsrun.hasProperty(self.logvalue)==False):
                raise ValueError("No log value given OR the given log value was not found in the file.")

            #need to split the data by an independt variable , some log value.
            #Create the array of logvalue BOUNDARIES
            if self.logvaluemin == None:
                self.logvaluemin= array(wsrun.getProperty(self.logvalue).value).min()
            if self.logvaluemax == None:
                self.logvaluemax= array(wsrun.getProperty(self.logvalue).value).max()
            if self.logvaluestep == None:
                self.logvaluestep = self.logvaluemax - self.logvaluemin

            bounds = arange(self.logvaluemin, self.logvaluemax+self.logvaluestep, self.logvaluestep)

            #Get the time correlation correct if you set the time correlation keyword.
            #To first approximation, set the time to zero for the first.
            for i in range(len(bounds)-1):
                FilterByLogValue(InputWorkspace="data",OutputWorkspace = 'dataslice', LogName= self.logvalue,MinimumValue=float(bounds[i]) ,MaximumValue = float(bounds[i+1]))
                dataslice=mtd['dataslice']
                if dataslice.getNumberEvents()>0:
                    values=array(dataslice.run().getProperty(self.logvalue).value)
                    self.datatext+= "Processing data for "+self.logvalue+" between "+str(bounds[i])+" and "+str(bounds[i+1])+", mean="+str(values.mean())+" std="+str(values.std())+"\n"
                    self.ProcessWorkspace('dataslice')                
                    if self.resetEnergyToNone:
                        self.efixed=None


    def ProcessWorkspace(self,datawsname):
        if self.efixed == None:
            if mtd[datawsname].run().hasProperty('EnergyRequest'):
                self.efixed =  float(mean(mtd[datawsname].run().getProperty('EnergyRequest').value))
            else:
                raise ValueError("no Efixed has been set, and not found in file.")

	    #get Ei, or use Ei
        if self.calce == False:
            efixed = self.efixed
            if self.t0==None:
                if self.instrument in ['HYSPEC','CNCS']:
                    t0 = self.t0fromei(efixed,self.instrument)
                else:
                    t0 = 0.0
            else:
                t0 = self.t0
            self.datatext += "User set value of incident energy, Ei="+str(efixed)+" meV, and t0="+str(t0)+" micro-seconds.\n"
       
	    #now deal with calculating the incident energy
        else:
            #check that the monitors are in memory
            try:
                mtd['monitorws']
            except:
                raise RuntimeError("monitor workspace not found")

            self.datatext += "Incident energy is calculated from monitor data.\n"

            mon1spec=mtd['monitorws'].getInstrument().getNumberParameter("ei-mon1-spec")[0]
            mon2spec=mtd['monitorws'].getInstrument().getNumberParameter("ei-mon2-spec")[0]

            alg=GetEiT0atSNS(MonitorWorkspace="monitorws",IncidentEnergyGuess=self.efixed)	        

            efixed = float(alg[0])
            t0     = float(alg[1])		
           
            #alg=GetEiT0atSNS(InputWorkspace="monitorws",Monitor1Spec=int(mon1spec),Monitor2Spec=int(mon2spec),EnergyEstimate=self.efixed)
            
            #alg=GetEiT0atSNS(InputWorkspace="monitorws",EnergyEstimate=self.efixed)	 	        

#            efixed = float(alg[0])
#            t0     = float(alg[3])		
	
            self.datatext += "Ei ="+str(efixed)+" meV, t0 ="+str(t0)+" microseconds\n"

        #Now adjust the data by the value of t0
        ChangeBinOffset(InputWorkspace=datawsname,OutputWorkspace=datawsname,Offset=-t0)

        #define Erange
        if self.emin == None:
            emin = -0.5*efixed
        else:
            emin = self.emin

        if self.emax == None:
            emax = efixed
        else:
            emax = self.emax

        if self.ebin == None:
            ebin = (emax-emin)/100.0
        else:
            ebin = self.ebin

        Erange = str(emin)+","+str(ebin)+","+str(emax)    


	    #Time-ind-bg subtraction.
        if (self.tibg):
            #check if tibmin and tibmax have been defined.
            tibmin = self.tibgstart
            tibmax = self.tibgstop
            if tibmin== None or tibmax == None:
                raise ValueError("Time independent background subtraction selected, but no limits set.")
            tibstep=tibmax-tibmin
            tibpar=str(tibmin)+","+str(tibstep)+","+str(tibmax)

            Rebin(InputWorkspace=datawsname,OutputWorkspace="background_origin_ws",Params=tibpar,PreserveEvents=False)
            ConvertUnits(InputWorkspace=datawsname,OutputWorkspace=datawsname,Target="DeltaE",EMode="Direct",Efixed=efixed)

	        #Do the Binning into energy bins
            Rebin(InputWorkspace=datawsname,OutputWorkspace=datawsname,Params=Erange,PreserveEvents=False)
            ConvertUnits(InputWorkspace=datawsname,OutputWorkspace=datawsname,Target="TOF",EMode="Direct",Efixed=efixed)

            ConvertToDistribution(Workspace=datawsname)
            CalculateFlatBackground(InputWorkspace="background_origin_ws",OutputWorkspace="background_ws",StartX=tibmin,EndX=tibmax,Mode="Mean",OutputMode="Return Background")
            ConvertToDistribution(Workspace="background_ws")
            Minus(LHSWorkspace=datawsname,RHSWorkspace="background_ws",OutputWorkspace=datawsname)
            ConvertFromDistribution(Workspace=datawsname)
            self.datatext  += "Time-independent background between "+str(tibmin)+" and "+str(tibmax)+" microseconds was subtracted.\n"
        else:
            self.datatext += "No time-independent background subtraction performed.\n"


        #normalize by charge
        NormaliseByCurrent(InputWorkspace=datawsname,OutputWorkspace=datawsname)
        w=mtd[datawsname]
        totalmuAhr = w.run().getProtonCharge()
        totalcoul  = totalmuAhr/1000*3.6
        self.datatext += "Data normalized by proton charge ("+str(totalmuAhr) + " micro-Ah), ("+str(totalcoul)+" C).\n"

        #detector wavelength sensitivity
        ConvertUnits(InputWorkspace=datawsname,OutputWorkspace=datawsname,Target="Wavelength",EMode="Direct",EFixed=efixed)
        He3TubeEfficiency(InputWorkspace=datawsname,OutputWorkspace=datawsname)												
        self.datatext += "Data corrected for He3 Tube Efficiency.\n"

        #Convert the data to units of energy
        ConvertUnits(InputWorkspace=datawsname,OutputWorkspace=datawsname,Target="DeltaE",EMode="Direct",EFixed=efixed)
        
	    #ki/kf
        if self.kiokf == True:
            CorrectKiKf(InputWorkspace=datawsname,OutputWorkspace=datawsname)
            self.datatext += "ki/kf factor has been applied to the data.\n" 

        #Rebinning the data if it was not already done.
        Rebin(InputWorkspace=datawsname,OutputWorkspace=datawsname,Params=Erange,PreserveEvents=False)
        self.datatext += "Data binned with emin="+str(emin)+", emax="+str(emax)+", ebin=" + str(ebin)+ " meV.\n"

        #normalize data by 1/(bin)
        #Convert to differential cross section by dividing by the energy bin width
        ConvertToDistribution(Workspace=datawsname)
        self.datatext += "Data converted to differential cross section by dividing by the energy bin width.\n"

	    #do the masking and calibration
        MaskDetectors(Workspace=datawsname,MaskedWorkspace="calibration")
        self.datatext += "Data have been masked by the calibration file.\n"


#        Divide(LHSWorkspace=datawsname,RHSWorkspace="calibration",OutputWorkspace=datawsname)
#        self.datatext += "Data have been normalized and masked by the calibration file.\n"

	    #deal with angles
        [psiangle, angletext]= definegoniometer(self.goniometermotor, self.goniometermotoroffset, self.goniometermotordirection, self.goniometermotoraxis, datawsname)
        self.datatext += angletext

        #deal with UB OR lattice definintion
        if self.lattice != None:
            #check if there are 12 elements.
            #remove any white space
            splitlattice = self.lattice.replace(' ','').split(',')
            if len(splitlattice) == 12:
                [a,b,c,alpha, beta, gamma, ux, uy, uz, vx, vy, vz] = splitlattice
            else:
                [a,b,c,alpha, beta, gamma, ux, uy, uz, vx, vy, vz] = ['1','1','1','90','90','90','1','0','0','0','1','0']
            u=ux+','+uy+','+uz
            v=vx+','+vy+','+vz
        if self.ub !=None:
            if len(self.ub.split(','))!=9:
                self.ub="0,0,0,0,0,0,0,0,0"
            else:
                self.ub.replace(' ','')
        else:
            self.ub="0,0,0,0,0,0,0,0,0"
        SetUB(Workspace=datawsname,a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma,u=u,v=v,UB=self.ub)
        self.datatext+="Set UB matrix with a="+a+",b="+b+",c="+c+",alpha="+alpha+",beta="+beta+",gamma="+gamma+",u="+u+",v="+v+",UB="+self.ub+"\n"

	    #grouping
        #in general grouping files will be stored in Mantid/instrument/Grouping
        #check if the requested grouping file is present
        #Two types of grouping, powder and pixel
        if ((self.grouping == 'powder') and (self.powderanglestep != None)):
            #do the powder work
            outputdir = os.path.abspath(os.curdir)
            GenerateGroupingPowder(InputWorkspace=datawsname,GroupingFilename=outputdir+"/powdergroup.xml",AngleStep=self.powderanglestep)
            changepermissions(outputdir+"powdergroup.xml")
            changepermissions(outputdir+"powdergroup.par")
            GroupDetectors(InputWorkspace=datawsname,OutputWorkspace=datawsname,MapFile=outputdir+"/powdergroup.xml",Behaviour="Sum")
            #GroupDetectors(InputWorkspace=datawsname,OutputWorkspace=datawsname,MapFile=outputdir+"/powdergroup.xml")
            #SolidAngle(InputWorkspace=datawsname,OutputWorkspace="sa")
            #Divide(LHSWorkspace=datawsname,RHSWorkspace="sa",OutputWorkspace=datawsname)
            #DeleteWorkspace(Workspace = "sa")
            self.datatext += "Detectors grouped by angle with a step of "+str(self.powderanglestep)+ ".\n"

            #Group the calibration file as a powder also
            GroupDetectors(InputWorkspace="calibration",OutputWorkspace="groupcalibration",MapFile=outputdir+"/powdergroup.xml",Behaviour="Sum")
            Divide(LHSWorkspace=datawsname,RHSWorkspace="groupcalibration",OutputWorkspace=datawsname)
            self.datatext += "Data have been normalized by the grouped calibration file.\n"






        #case of powder grouping and NO valid angle step
        elif ((self.grouping == 'powder') and (self.powderanglestep == None)):
            raise ValueError("Powder grouping chosen, but anglestep is invalid.")

        elif (self.grouping != None and len(self.grouping)>0):
            #check that the self.grouping value has the correct form "NXM" with the correct integers.
            twoelem = self.grouping.split('X')
            grouppath = ""
            if len(twoelem) == 2:
                #good, now check that they are the correct value
                try:
                    pixely = int(twoelem[0])
                    pixelx = int(twoelem[1])
                except:
                    pixelx = 0
                    pixely = 0
                if (pixely in [1,2,4,8,16,32,64,128]) and (pixelx in [1,2,4,8]):
                    #generate the grouping file in the current directory (if it doesn't already exist)
                    grouppath = os.path.abspath(os.curdir)+'/'+self.instrument+"_"+self.grouping+"_grouping.xml"
                    if not (os.path.isfile(grouppath)):
                        #if it is not there, then make it
                        GenerateGroupingSNSInelastic(AlongTubes=str(pixely),AcrossTubes=str(pixelx),Instrument=self.instrument,Filename=grouppath)
            if len(grouppath) < 1:
                grouppath = self.grouping #Filename may be specified by user for any grouping they want, without using the standard "NXM" schema.
            if (os.path.isfile(grouppath)):
                GroupDetectors(InputWorkspace=datawsname,OutputWorkspace=datawsname,MapFile=grouppath,Behaviour="Sum")
                GroupDetectors(InputWorkspace="calibration",OutputWorkspace="groupcalibration",MapFile=grouppath,Behaviour="Sum")
                Divide(LHSWorkspace=datawsname,RHSWorkspace="groupcalibration",OutputWorkspace=datawsname)
	
#        Divide(LHSWorkspace=datawsname,RHSWorkspace="calibration",OutputWorkspace=datawsname)
#        self.datatext += "Data have been normalized and masked by the calibration file.\n"


                self.datatext += "Detectors grouped by averaging using the "+grouppath+" file.\n"
                self.datatext += "Data have been normalized by the grouped calibration file.\n"
            else:
                raise ValueError("Grouping file "+grouppath+" NOT FOUND.")






        #Now scale the data if the scale factor is listed, default operation is multiply.
        if self.scalefactor != None:
            Scale(InputWorkspace=datawsname,OutputWorkspace=datawsname,Factor=self.scalefactor)
            self.datatext += "Dave have been scaled by a factor of "+str(self.scalefactor)+".\n"

 
         #Now deal with saving files.
        if self.save != None:
            #create a friendly name
            friendlynamebase = self.CreateFriendlyFilename(datawsname)

            #Parse the filetypes to save
            if 'nxspe' in self.save:
                #if there is a powder maping file in the current directory, then use it with the .nxspe file
                if self.grouping == 'powder':
                    SaveNXSPE(Filename=friendlynamebase+".nxspe", InputWorkspace=datawsname, Efixed=str(efixed),Psi=str(psiangle), KiOverKfScaling=self.kiokf, ParFile=os.path.abspath(os.curdir) +"/powdergroup.par")
                else:
                    SaveNXSPE(Filename=friendlynamebase+".nxspe", InputWorkspace=datawsname, Efixed=str(efixed),Psi=str(psiangle), KiOverKfScaling=self.kiokf)
                self.datatext += "Data have been saved as a .nxspe file, FILENAME="+friendlynamebase+".nxspe.\n"
                #change permissions of the directory and file
                changepermissions(friendlynamebase+".nxspe")
 
            if 'nxs' in self.save:
                #save the nxs
                SaveNexus(Filename=friendlynamebase+".nxs", InputWorkspace=datawsname)
                self.datatext += "Data have been saved as a .nxs file, FILENAME="+friendlynamebase+".nxs.\n"
                #change permissions of the directory and file
                changepermissions(friendlynamebase+".nxs")

            #Save the diffraction pattern as I(Q) three column text format
            if 'iofq' in self.save:
                #get the qmin and q max
                if self.qmin == None:
                    [qmin, tempqmax] = calqrangefromworkspace(datawsname)
                else:
                    qmin = self.qmin

                if self.qmax == None:
                    [tempqmin, qmax] = calqrangefromworkspace(datawsname)
                else:
                    qmax = self.qmax

                #check that qbining has been set
                if self.qstep == None:
                    #if not set then set it to 0.1 inv angstroms
                    qstep = (qmax-qmin)/150.0
                else:
                    qstep = self.qstep

                qbinparams = str(qmin)+','+str(qstep)+','+str(qmax)

                ws = mtd[datawsname]
                #get the energy binning using the binning set by emin, and emax with a binsize of the full range.
                efullbin = (emax-emin)*1.01
                ebinparams = str(emin)+","+str(efullbin)+","+","+str(emax)
                SofQW3(InputWorkspace=datawsname,OutputWorkspace='SofQWdata',QAxisBinning=qbinparams,Emode="Direct",Efixed=efixed)
                Transpose(InputWorkspace='SofQWdata',OutputWorkspace='SofQWdata')
                wsofq = Rebin2D(InputWorkspace='SofQWdata',Axis1Binning=qbinparams,Axis2Binning=ebinparams,UseFractionalArea=True)
                SaveAscii(Filename=friendlynamebase+"_iofq.dat",InputWorkspace='wsofq')
                self.datatext += "Data have been saved as a iofq.dat file, FILENAME="+friendlynamebase+"_iofq.dat.\n" 
                #change permissions of the directory and file
                changepermissions(friendlynamebase+"_iofq.dat")

            if 'mdnxs' in self.save:
                #convert to MD events
                #calculate Ki
                ki = sqrt(efixed/2.072) #units of Angstroms inverse
                #calculate Kf
                kf = sqrt((efixed-emin)/2.072) # units of Angstroms inverse

                Qvecmax = ki + kf # Angstroms inverse
                #get the lattice parameters from setting the UB
                ol = mtd[datawsname].sample().getOrientedLattice()
                hmax = Qvecmax*ol.a()/2/pi
                kmax = Qvecmax*ol.b()/2/pi
                lmax = Qvecmax*ol.c()/2/pi
                minval = str(-hmax)+','+str(-kmax)+','+str(-lmax)+','+str(emin)
                maxval = str(hmax)+','+str(kmax)+','+str(lmax)+','+str(emax)

                ConvertToMD(InputWorkspace = datawsname, OutputWorkspace = datawsname+"MD",QDimensions = "Q3D", QConversionScales="HKL", dEAnalysisMode="Direct", MinValues=minval , MaxValues=maxval,MaxRecursionDepth='1' )
                #save the Md workspace.
                SaveMD(InputWorkspace=datawsname+"MD", Filename=friendlynamebase+"_MD.nxs")
                self.datatext += "Data have been saved as a MD nexus file, FILENAME="+friendlynamebase+"_MD.nxs.\n" 
                #change permissions of the directory and file
                changepermissions(friendlynamebase+"_MD.nxs")


            if 'par' in self.save:
                if self.grouping == 'powder':
                    shutil.copy(os.path.abspath(os.curdir) +"/powdergroup.par",friendlynamebase+".par")
                else:
                    SavePAR(Filename=friendlynamebase+".par", InputWorkspace=datawsname)
                self.datatext += "PAR file has been saved as FILENAME="+friendlynamebase+".par.\n"
                #change permissions of the directory and file
                changepermissions(friendlynamebase+".par")

            if 'phx' in self.save:
                SavePHX(Filename=friendlynamebase+".phx",InputWorkspace=datawsname)
                self.datatext += "PHX file has been saved as FILENAME="+friendlynamebase+".phx.\n"

            if 'spe' in self.save:
                SaveSPE(Filename=friendlynamebase+".spe",InputWorkspace=datawsname)
                self.datatext += "SPE file has been saved as FILENAME="+friendlynamebase+".spe.\n"

            if 'sqw' in self.save:
                #get the qmin and q max
                if self.qmin == None:
                    [qmin, tempqmax] = calqrangefromworkspace(datawsname)
                else:
                    qmin = self.qmin

                if self.qmax == None:
                    [tempqmin, qmax] = calqrangefromworkspace(datawsname)
                else:
                    qmax = self.qmax

                #check that qbining has been set
                if self.qstep == None:
                    #if not set then set it to 150 qbins for the full range
                    qstep = (qmax-qmin)/150.0
                else:
                    qstep = self.qstep

                qbinparams = str(qmin)+","+str(qstep)+","+str(qmax)

                SofQW3(InputWorkspace=datawsname,OutputWorkspace='SofQWdata',QAxisBinning=qbinparams,Emode="Direct",Efixed=efixed)
                #save the nxs
                SaveNexus(Filename=friendlynamebase+"_sqw3.nxs", InputWorkspace='SofQWdata')
                self.datatext += "Sqw3 Data have been saved as a .nxs file, FILENAME="+friendlynamebase+".nxs.\n"
                #change permissions of the directory and file
                changepermissions(friendlynamebase+"_sqw3.nxs")


            if 'iofe' in self.save:
                #get the q binning.
                if self.qmin == None:
                    [qmin, tempqmax] = calqrangefromworkspace(datawsname)
                else:
                    qmin = self.qmin

                if self.qmax == None:
                    [tempqmin, qmax] = calqrangefromworkspace(datawsname)
                else:
                    qmax = self.qmax

                #check that qbining has been set
                if self.qstep == None:
                    #if not set then set it to 150 bins for the full range
                    qstep = (qmax-qmin)/150.0
                else:
                    qstep = self.qstep

                qbinparams = str(qmin)+","+str(qstep)+","+str(qmax)
                qfullstep = (qmax-qmin)*1.01
                qfullstepbinparams = str(qmin)+","+str(qfullstep)+","+str(qmax)

                SofQW3(InputWorkspace=datawsname,OutputWorkspace='SofQWdata',QAxisBinning=qbinparams,Emode="Direct",Efixed=efixed)
                Transpose(InputWorkspace='SofQWdata',OutputWorkspace='SofQWdata')
                wsofe = Rebin2D(InputWorkspace='SofQWdata',Axis1Binning=qfullstepbinparams,Axis2Binning=Erange,UseFractionalArea=True,Transpose=True)
                SaveAscii(Filename=friendlynamebase+"_iofe.dat",InputWorkspace='wsofe')
                self.datatext += "Data have been saved as a iofe.dat file, FILENAME="+friendlynamebase+"_iofe.dat.\n" 
                #change permissions of the directory and file
                changepermissions(friendlynamebase+"_iofe.dat")
 

            if 'iofphiecolumn' in self.save:
                #convert to two-theta energy space
                ConvertSpectrumAxis(InputWorkspace=datawsname,OutputWorkspace="iofphiwdata",Target='theta',Emode="Direct", Efixed=efixed)
                #rebin the data in energy and phi space
                ebinparams = str(emin)+","+str(ebin)+","+","+str(emax)

                #Check if powderanglemin and max were set.  if they were not set then need to determine them from the data
                w=mtd['iofphiwdata']
                if self.powderanglemin == None:
                    angleaxis = w.getAxis(1)
                    powderanglemin = angleaxis.getMin()
                else:
                    powderanglemin = self.powderanglemin

                if self.powderanglemax == None:
                    angleaxis = w.getAxis(1)
                    powderanglemax = angleaxis.getMax()
                else:
                    powderanglemax = self.powderanglemax

                anglebins = str(powderanglemin)+","+str(self.powderanglestep)+","+str(powderanglemax)
                Rebin2D(InputWorkspace='iofphiwdata',OutputWorkspace='iofphiwdataRebin',Axis1Binning=ebinparams,Axis2Binning=anglebins)

                w1 = mtd['iofphiwdataRebin']

                #convert data from the bin edges to bin centeres
                w2 = ConvertToPointData('iofphiwdataRebin')
                energy = w2.extractX()
                intensity = w2.extractY()
                error = w2.extractE()

                #Need to deal with the angle dimension seperately (go back to the w1 data)
                a = w1.getAxis(1)
                anglebins = a.extractValues()
                #get the bin centers, not the bin edges
                anglecent = (anglebins[1:] + anglebins[:-1])*0.5
                #get an array to match the dimensions of the energy, intensity and error arrays
                anglecentarr = transpose(tile(anglecent,(w2.blocksize(),1)))

                #fltatten all 4 arrays.
                energy = energy.reshape(-1)
                intensity = intensity.reshape(-1)
                error = error.reshape(-1)
                anglecentarr = anglecentarr.reshape(-1)
                
                iofphiecolumn = transpose(array([anglecentarr,energy,intensity,error]))

                savetxt(friendlynamebase+"_iofphiecolumn.dat",iofphiecolumn)
                self.datatext += "Data have been saved as a iofphiecolumn.dat file, FILENAME="+friendlynamebase+"_iofphiecolumn.dat.\n" 
                changepermissions(friendlynamebase+"_iofphiecolumn.dat")
                

            if 'iofphiearray' in self.save:
                #convert to two-theta energy space
                ConvertSpectrumAxis(InputWorkspace=datawsname,OutputWorkspace="iofphiwdata",Target='theta',Emode="Direct", Efixed=efixed)
                #rebin the data in energy and phi space
                ebinparams = str(emin)+","+str(ebin)+","+","+str(emax)

                #Check if powderanglemin and max were set.  if they were not set then need to determine them from the data
                w=mtd['iofphiwdata']
                if self.powderanglemin == None:
                    angleaxis = w.getAxis(1)
                    powderanglemin = angleaxis.getMin()
                else:
                    powderanglemin = self.powderanglemin

                if self.powderanglemax == None:
                    angleaxis = w.getAxis(1)
                    powderanglemax = angleaxis.getMax()
                else:
                    powderanglemax = self.powderanglemax

                anglebins = str(powderanglemin)+","+str(self.powderanglestep)+","+str(powderanglemax)
                Rebin2D(InputWorkspace='iofphiwdata',OutputWorkspace='iofphiwdataRebin',Axis1Binning=ebinparams,Axis2Binning=anglebins)

                w1 = mtd['iofphiwdataRebin']

                #convert data from the bin edges to bin centeres
                w2 = ConvertToPointData('iofphiwdataRebin')
                energy = w2.extractX()
                intensity = transpose(w2.extractY())
                error = transpose(w2.extractE())

                #Need to deal with the angle dimension seperately (go back to the w1 data)
                a = w1.getAxis(1)
                anglebins = a.extractValues()
                #get the bin centers, not the bin edges
                anglecent = (anglebins[1:] + anglebins[:-1])*0.5
                #get an array to match the dimensions of the energy, intensity and error arrays
                anglecentarr = transpose(tile(anglecent,(w2.blocksize(),1)))

                evals = w2.getAxis(0)
                evalsbincenters = evals.extractValues()

                phiearrayfile = open(friendlynamebase+"_iofphiearray.dat",'w')
                phiearrayfile.write(str(intensity.shape[0])+" "+str(intensity.shape[1])+'\n')
                phiearrayfile.write('\n')
                for val in anglecent:
                    phiearrayfile.write(str(val)+'\n')
                phiearrayfile.write('\n')
                for val in evalsbincenters:
                    phiearrayfile.write(str(val)+'\n')
                phiearrayfile.write('\n')
                for line in intensity:
                    phiearrayfile.write(" ".join(str(elem) for elem in line) + "\n")
                phiearrayfile.write('\n')
                for line in error:
                    phiearrayfile.write(" ".join(str(elem) for elem in line) + "\n")
                phiearrayfile.close()
                
                self.datatext += "Data have been saved as a iofphiearray.dat file, FILENAME="+friendlynamebase+"_iofphiearray.dat.\n" 
                changepermissions(friendlynamebase+"_iofphiearray.dat")
             

            if 'iofqecolumn' in self.save:
                #get the q binning.
                if self.qmin == None:
                    [qmin, tempqmax] = calqrangefromworkspace(datawsname)
                else:
                    qmin = self.qmin

                if self.qmax == None:
                    [tempqmin, qmax] = calqrangefromworkspace(datawsname)
                else:
                    qmax = self.qmax

                #check that qbining has been set
                if self.qstep == None:
                    #if not set then set it to 150 bins for the full range
                    qstep = (qmax-qmin)/150.0
                else:
                    qstep = self.qstep

                qbinparams = str(qmin)+","+str(qstep)+","+str(qmax)
                qfullstep = (qmax-qmin)*1.01
                qfullstepbinparams = str(qmin)+","+str(qfullstep)+","+str(qmax)

                SofQW3(InputWorkspace=datawsname,OutputWorkspace='SofQWdata',QAxisBinning=qbinparams,Emode="Direct",Efixed=efixed)

                w1 = mtd['SofQWdata']

                #convert data from the bin edges to bin centeres
                w2 = ConvertToPointData('SofQWdata')
                energy = w2.extractX()
                intensity = w2.extractY()
                error = w2.extractE()

                #Need to deal with the Q dimension seperately (go back to the w1 data)
                a = w1.getAxis(1)
                qbins = a.extractValues()
                #get the bin centers, not the bin edges
                qcent = (qbins[1:] + qbins[:-1])*0.5
                #get an array to match the dimensions of the energy, intensity and error arrays
                qcentarr = transpose(tile(qcent,(w2.blocksize(),1)))

                #fltatten all 4 arrays.
                energy = energy.reshape(-1)
                intensity = intensity.reshape(-1)
                error = error.reshape(-1)
                qcentarr = qcentarr.reshape(-1)
                
                iofqecolumn = transpose(array([qcentarr,energy,intensity,error]))

                savetxt(friendlynamebase+"_iofqecolumn.dat",iofqecolumn)
                self.datatext += "Data have been saved as a iofqecolumn.dat file, FILENAME="+friendlynamebase+"_iofqecolumn.dat.\n" 
                changepermissions(friendlynamebase+"_iofqecolumn.dat")
                

            if 'iofqearray' in self.save:
                #get the q binning.
                if self.qmin == None:
                    [qmin, tempqmax] = calqrangefromworkspace(datawsname)
                else:
                    qmin = self.qmin

                if self.qmax == None:
                    [tempqmin, qmax] = calqrangefromworkspace(datawsname)
                else:
                    qmax = self.qmax

                #check that qbining has been set
                if self.qstep == None:
                    #if not set then set it to 150 bins for the full range
                    qstep = (qmax-qmin)/150.0
                else:
                    qstep = self.qstep

                qbinparams = str(qmin)+","+str(qstep)+","+str(qmax)
                qfullstep = (qmax-qmin)*1.01
                qfullstepbinparams = str(qmin)+","+str(qfullstep)+","+str(qmax)

                SofQW3(InputWorkspace=datawsname,OutputWorkspace='SofQWdata',QAxisBinning=qbinparams,Emode="Direct",Efixed=efixed)

                w1 = mtd['SofQWdata']

                #convert data from the bin edges to bin centeres
                w2 = ConvertToPointData('SofQWdata')
                energy = w2.extractX()
                intensity = transpose(w2.extractY())
                error = transpose(w2.extractE())

                #Need to deal with the Q dimension seperately (go back to the w1 data)
                a = w1.getAxis(1)
                qbins = a.extractValues()
                #get the bin centers, not the bin edges
                qcent = (qbins[1:] + qbins[:-1])*0.5
                #get an array to match the dimensions of the energy, intensity and error arrays
                qcentarr = transpose(tile(qcent,(w2.blocksize(),1)))

                evals = w2.getAxis(0)
                evalsbincenters = evals.extractValues()


                phiearrayfile = open(friendlynamebase+"_iofqearray.dat",'w')
                phiearrayfile.write(str(intensity.shape[0])+" "+str(intensity.shape[1])+'\n')
                phiearrayfile.write('\n')
                for val in qcent:
                    phiearrayfile.write(str(val)+'\n')
                phiearrayfile.write('\n')
                for val in evalsbincenters:
                    phiearrayfile.write(str(val)+'\n')
                phiearrayfile.write('\n')
                for line in intensity:
                    phiearrayfile.write(" ".join(str(elem) for elem in line) + "\n")
                phiearrayfile.write('\n')
                for line in error:
                    phiearrayfile.write(" ".join(str(elem) for elem in line) + "\n")
                phiearrayfile.close()
                
                self.datatext += "Data have been saved as a iofqearray.dat file, FILENAME="+friendlynamebase+"_iofqearray.dat.\n" 
                changepermissions(friendlynamebase+"_iofqearray.dat")
                

            if 'vannorm' in self.save:               
                if self.vanmine == None:
                    vanmine = -40.0
                else:
                    vanmine = self.vanmine
                if self.vanmaxe == None:
                    vanmaxe =  40.0
                else:
                    vanmaxe = self.vanmaxe
                if vanmaxe > emax:
                    vanmaxe = emax
                if vanmine < emin:
                    vanmine = emin
                vanebinparams = str(vanmine)+","+str(ebin)+","+","+str(vanmaxe)


                #get angle1 and angle2 (these are the min and max scattering angles).
                angle1 = 125.0 #max angle
                angle2 = 4.0   #min angle

 
                #Check if vanminangle and max were set.  if they were not set then need to determine them from the data
                if self.vanminangle == None:
                    vanminangle = angle2
                else:
                    vanminangle = self.vanminangle
                if self.vanmaxangle == None:
                    vanmaxangle = angle1
                else:
                    vanmaxangle = self.vanmaxangle
                
                #NO weighting, but include the zero intensities, use the the error in zero counts
                #is changed to be the minimum non-zero error.
  
                #Integrate over the energy axis.
                Integration(InputWorkspace=datawsname, OutputWorkspace='van_int', RangeLower=vanmine, RangeUpper=vanmaxe,IncludePartialBins='1')

                ### A significant portion of the following code between the dollar signs is lifted from the inelastic scripts of isis on github.com
                #$
                van_ws = mtd['van_int']
                nhist = van_ws.getNumberHistograms()

                signal1_sum = 0.0
                weight1_sum = 0.0
                err1_sum =    0.0

                signal1_mean = 0.0
                signal1_meanofsqr  = 0.0
                
                ic = 0
                izerc=0
                ioutsideangle=0
                imask = 0


                cntrlist = []
                intlist = []
                errlist = []
                scatanglelist = []
                effectivecountlist = []
                effectivescalefactorlist = []
                #summaryfilename = friendlynamebase+"_summary.txt"
                #parentdir = os.path.dirname(summaryfilename)    

                #determine the smallest non-zero error value.  This will correspond to the
                #error associated with the smallest number of counts
                errlistprelim=[]
                for i in range(nhist):
                    errlistprelim.append(van_ws.readE(i)[0])
                errlist2 = [a for a in errlistprelim if a!= 0]
                minerr = min(errlist2)
              

                for i in range(nhist):
                    try:
                        det=van_ws.getDetector(i)
                    except Exception:
                        continue
                    if det.isMasked():
                        imask += 1
                        continue

                    origin = van_ws.getInstrument().getSample().getPos()
                    scatteringangle = van_ws.getDetector(i).getTwoTheta(origin,V3D(0,0,1))*180/pi
                    if scatteringangle > vanmaxangle or scatteringangle < vanminangle:
                        ioutsideangle += 1
                        continue

                    signal = van_ws.readY(i)[0]
                    error  = van_ws.readE(i)[0]
                    if error==0:
                        error = minerr

                    if signal != signal:                #ignore NaNs
                        continue

                    cntrlist.append(i)
                    intlist.append(signal)
                    errlist.append(error)
                    scatanglelist.append(scatteringangle)

                    effectivecount = 0
                    if (signal!=0) and (error != 0):
                        effectivecount = signal*signal/error/error
                    effectivecountlist.append(effectivecount)

                    effectivescalefactor = 0
                    if (signal!=0) and (error != 0):
                        effectivescalefactor = signal/effectivecount
                    effectivescalefactorlist.append(effectivescalefactor)

                    weight = 1.0
                    signal1_sum += signal
                    weight1_sum += weight

                    signal1_mean += signal
                    signal1_meanofsqr  += signal*signal

                    ic += 1

                integral_monovanLibISIS=signal1_sum/weight1_sum

                err1_sum = sqrt(ic)/weight1_sum

                signal1_mean /= ic
                signal1_meanofsqr  /= ic
 

                van_rmm = 50.942
                if self.vanmass == None:
                    vanmass = 1.0
                else:
                    vanmass = self.vanmass
                #integral_monovan=signal_sum /(wbVan_sum)
                van_multiplier = (float(van_rmm)/float(vanmass))
                absnorm_factorLibISIS = integral_monovanLibISIS * van_multiplier
                err1_sum= err1_sum*van_multiplier


                if efixed >= 210.0:
                    xsection = 421.0  # vanadium cross-section in mBarn/sR (402 mBarn/Sr) (!!!modified to fit high energy limit?!!!)
                else: # old textbook cross-section for vanadium for ei=20mEv
                    xsection = 400.0 + (efixed/10.0)  

                absnorm_factorLibISIS /= xsection
                err1_sum /=xsection

                if self.samplemass == None:
                    samplemass = 1.0
                else:
                    samplemass = self.samplemass

                if self.samplermm == None:
                    samplermm = 1.0
                else:
                    samplermm = self.samplermm

                sample_multiplier = (float(samplemass)/float(samplermm))

                absnorm_factorLibISIS= absnorm_factorLibISIS *sample_multiplier
                err1_sum = err1_sum*sample_multiplier

                inverr1_sum = err1_sum/absnorm_factorLibISIS/absnorm_factorLibISIS             

                    #$
                self.datatext += "\n"                
                self.datatext += "--------------ABSOLUTE UNITS FROM MONOCHROMATIC VANADIUM----INCLUDE zeros, Straight AVG----------\n"
                self.datatext += "--------------ERROR IN ZEROS CHANGED TO BE THE MINIMUM ERROR VALUE -------------------------------\n"
                self.datatext += "Monochormatic vanadium energy integration range:  min="+str(vanmine)+" max="+str(vanmaxe)+" meV.\n"
                self.datatext += "Monochromatic vanadium two-theta angle integration range:  min="+str(vanminangle) +" max="+str(vanmaxangle)+" degrees.\n"
                self.datatext += "Summed "+str(ic)+" spectra with total value "+str(signal1_sum)+" and total weight "+str(weight1_sum)+".\n"
                self.datatext += "Dropped "+str(izerc)+" empty spectra.\n"
                self.datatext += "Dropped "+str(ioutsideangle)+" spectra outside of the angular range chosen.\n"
                self.datatext += "Dropped "+str(imask)+" spectra that were masked.\n"
                self.datatext += "mean of the signal <I>:  "+str(signal1_mean)+".\n"
                self.datatext += "mean of the square of signal <I^2>:  "+str(signal1_meanofsqr)+".\n"
                self.datatext += "\n"
                self.datatext += "Atomic mass of vanadium:  "+str(van_rmm)+" grams per mole.\n"
                self.datatext += "Mass of vanadium standard:  "+str(vanmass)+" grams.\n"
                self.datatext += "Van multiplier:  "+str(van_multiplier)+" 1/moles of vanadium.\n"
                self.datatext += "Moles of vanadium:  "+str(1/van_multiplier)+".\n"
                self.datatext += "\n"
                self.datatext += "Formula weight of sample:  "+str(samplermm)+" grams per mole.\n"
                self.datatext += "Mass of sample:  "+str(samplemass)+" grams.\n"
                self.datatext += "sample multiplier:  "+str(sample_multiplier)+" moles of sample.\n"
                self.datatext += "\n"
                self.datatext += "Xsection used for vanadium:  "+str(xsection)+" mBarn/sR.\n"
                self.datatext += "\n"
                self.datatext += "Absolute normalization factor:     "+str(1.0/absnorm_factorLibISIS)+" +/- "+str(inverr1_sum)+"\n"
                self.datatext += "\n"
                self.datatext += "In order to obtain absolute units of mBarn/meV/Sr/formula unit, one multiplies the data by the above number.\n"
                self.datatext += "\n"
                self.datatext += "NOTE: errors in normalization factors are small compared to other errors not accounted for in this calculation:\n"
                self.datatext += "     shape anisotropy, guide effects, neutron absorption effects, vanadium composition, \n"
                self.datatext += "     error in sample and vanadium mass, etc.\n"
                self.datatext += "-----------------------------------------------------------------------------\n"

                #save the list of int values to a file also.
                intlistfilename = friendlynamebase+"_vannorm_intlist.txt"
                parentdir = os.path.dirname(intlistfilename)    
                if os.path.isdir(parentdir) == False:
                    os.mkdir( parentdir, 0o755 )
                f = open(intlistfilename,"w")
                f.write("Counter"+"\t"+"Intensity"+"\t"+"Error"+"\t"+"twotheta"+"\t"+"EffectiveCounts"+"\t"+"EffectiveScaleFactor"+"\n")
                for i,j,k,l,m,n in zip(cntrlist,intlist,errlist,scatanglelist,effectivecountlist,effectivescalefactorlist):
                    f.write(str(i)+'\t'+str(j)+str(k)+'\t'+str(l)+'\t'+str(m)+'\t'+str(n)+'\n')
                # COMMAND FOR WRITING At SINGLE VALUEf.write("\n".join(map(lambda x: str(x), intlist)))
                f.close
                

            #writing the summary file
            if 'summary' or 'vannorm' in self.save:
                #check if the friendlyname directory exists.
                summaryfilename = friendlynamebase+"_summary.txt"
                parentdir = os.path.dirname(summaryfilename)    
                if os.path.isdir(parentdir) == False:
                    os.mkdir(parentdir,0o755)
                sumfile = open(summaryfilename, 'w')
                sumfile.write("-----------VANADIUM CALIBRATION AND MASKING-----------\n")
                sumfile.write(self.calibrationtext)
                sumfile.write("\n-----------DATA REDUCTION-----------------------------\n")
                sumfile.write(self.datatext)
                sumfile.close()
                changepermissions(summaryfilename)


    def SetInstrument(self,instrument):
        if instrument not in ['ARCS','CNCS','HYSPEC','SEQUOIA']:
            raise ValueError("Instrument not defined")
        else:
            self.instrument=instrument
            self.LoadInstrumentSettings(instrument)

    def SetFilterBadPulses(self,value):
        self.fileterbadpulses=value

    def SetMask(self,algorithm=None,**kwargs):
        pass

    def LoadInstrumentSettings(self,instrumentname):
        """
        instrumentname is a string
        """
        modulename = instrumentname.lower()+'default'
        #We have the instrument name.
        #Try to import the default file for that instrument
        #First try to import from the /SNS/"instrumentname"/shared
        try:
            path = "/SNS/"+instrumentname+"/shared/"+modulename+'.py'
            if (os.path.isfile(path)):
                m = imp.load_source(modulename,path)
            else:
                m = imp.load_source(modulename,os.path.abspath(os.curdir)+"/"+modulename+'.py')
        except:
            raise RuntimeError("Could not find instrument definition module"+os.path.abspath(os.curdir)+"/"+modulename+'.py')
        params=m.instrumentparameters()
        self.LoadParameters(params)


    def t0fromei(self,ei,instrumentname):
        modulename = instrumentname.lower()+'default'
        try:
            path = "/SNS/"+instrumentname+"/shared/"+modulename+'.py'
            if (os.path.isfile(path)):
                m = imp.load_source(modulename,path)
            else:
                m = imp.load_source(modulename,os.path.abspath(os.curdir)+"/"+modulename+'.py')
        except:
            raise RuntimeError("Could not find instrument definition module"+os.path.abspath(os.curdir)+"/"+modulename+'.py')
        return m.t0fromei(ei)


    def LoadParameters(self,params):
        """
        Load parameters from a dictionary
        """
        if 'instrument' in list(params.keys()):
            if self.instrument==None:
                self.instrument=params['instrument']
            else:
                if self.instrument!=params['instrument']:
                    raise RuntimeError("instruments not compatible")
        if 'filterbadpulses' in params:
            self.filterbadpulses=params['filterbadpulses']
        if 'vanruns' in params:
            self.vanruns=params['vanruns']
        if 'units' in params:
            self.units=params['units']
        if 'vanmin' in params:
            self.vanmin=params['vanmin']
        if 'vanmax' in params:
            self.vanmax=params['vanmax']
        if 'processedfilename' in params:
            self.processedfilename=params['processedfilename']
        if 'maskfilename' in params:
            self.maskfilename=params['maskfilename']
        if 'mask' in params:
            self.mask=params['mask']
        if 'normalizedcalibration' in params:
            self.normalizedcalibration=params['normalizedcalibration']
        if 'ipts' in params:
            self.ipts=params['ipts']
        if 'runs' in params:
            self.runs=params['runs']
        if 'efixed' in params:
            self.efixed=params['efixed']
        if 't0' in params:
            self.t0=params['t0']
        if 'calce' in params:
            self.calce=params['calce']
        if 'ei-mon1-spec' in params:
            self.eimon1=params['ei-mon1-spec']
        if 'ei-mon2-spec' in params:
            self.eimon2=params['ei-mon2-spec']
        if 'emin' in params:
            self.emin=params['emin']
        if 'emax' in params:
            self.emax=params['emax']
        if 'ebin' in params:
            self.ebin=params['ebin']
        if 'qstep' in params:
            self.qstep=params['qstep']
        if 'qmax' in params:
            self.qmax=params['qmax']
        if 'qmin' in params:
            self.qmin=params['qmin']
        if 'kiokf' in params:
            self.kiokf=params['kiokf']
        if 'tibg' in params:
            self.tibg=params['tibg']
        if 'tibgstart' in params:
            self.tibgstart=params['tibgstart']
        if 'tibgstop' in params:
            self.tibgstop=params['tibgstop']
        if 'grouping' in params:
            self.grouping=params['grouping']
        if 'powderanglestep' in params:
            self.powderanglestep=params['powderanglestep']
        if 'powderanglemin' in params:
            self.powderanglemin=params['powderanglemin']      
        if 'powderanglemax' in params:
            self.powderanglemax=params['powderanglemax']      
        if 'goniometermotor' in params:
            self.goniometermotor=params['goniometermotor']
        if 'goniometermotoroffset' in params:
            self.goniometermotoroffset=params['goniometermotoroffset']
        if 'goniometermotoraxis' in params:
            self.goniometermotoraxis=params['goniometermotoraxis']
        if 'goniometermotordirection' in params:
            self.goniometermotordirection=params['goniometermotordirection']
        if 'save' in params:
            self.save=params['save']
        if 'friendlyname' in params:
            self.friendlyname=params['friendlyname']
        if 'scantype' in params:
            self.scantype=params['scantype']
        if 'logvalue' in params:
            self.logvalue=params['logvalue']
        if 'logvaluestep' in params:
            self.logvaluestep=params['logvaluestep']
        if 'logvaluemin' in params:
            self.logvaluemin=params['logvaluemin']
        if 'logvaluemax' in params:
            self.logvaluemax=params['logvaluemax']
        if 'friendlynamelogs' in params:
            self.friendlynamelogs=params['friendlynamelogs']
        if 'vanpath' in params:
            self.vanpath= params['vanpath']
        if 'datapath' in params:
            self.datapath = params['datapath']
        if 'lattice' in params:
            self.lattice = params['lattice']
        if 'ub' in params:
            self.ub = params['ub']
        if 'filternames' in params:
            self.filternames = params['filternames']
        if 'filtermin' in params:
            self.filtermin = params['filtermin']
        if 'filtermax' in params:
            self.filtermax = params['filtermax']
        if 'datadate' in params:
            self.datadate = params['datadate']
        if 'vanmine' in params:
            self.vanmine = params['vanmine']
        if 'vanmaxe' in params:
            self.vanmaxe = params['vanmaxe']
        if 'vanminangle' in params:
            self.vanminangle = params['vanminangle']
        if 'vanmaxangle' in params:
            self.vanmaxangle = params['vanmaxangle']
        if 'vanmass' in params:
            self.vanmass = params['vanmass']
        if 'samplemass' in params:
            self.samplemass = params['samplemass']
        if 'samplermm' in params:
            self.samplermm = params['samplermm']
        if 'scalefactor' in params:
            self.scalefactor = params['scalefactor']


    def PerformCalibration(self):
        #first check if there is currently a file ready to use
        #if the processedfilename exists in the current directory,
        #then load it and exit.
        if (self.processedfilename!=None and (os.path.isfile(self.processedfilename))):
            path=os.getcwd()+'/'+self.processedfilename
            Load(Filename = path, OutputWorkspace = 'calibration')
            self.calibrationtext += "Calibration loaded from "+path+"\n"
            print(("Calibration file "+path+" loaded."))
        #otherwise, this file does not exist, and it must be made.
        else:
            self.calibrationtext+="Performing calibration:\n"
            if self.vanruns == None:
                #calibration=LoadEmptyInstrument(Filename=config['instrumentDefinition.directory']+self.instrument+'_Definition_20121011-.xml',
                #                DetectorValue='1')
                #2/8/2013 changed to automatic loading of instrument definition file
                #if no date given use most current.
                if self.datadate == None:
                    calibration = LoadEmptyInstrument(Filename=mantid.api.ExperimentInfo.getInstrumentFilename(self.instrument),DetectorValue='1')
                    configfilestring = mantid.api.ExperimentInfo.getInstrumentFilename(self.instrument)


                #if date has been given use the appropriate version by date.
                else:
                    calibration = LoadEmptyInstrument(Filename=mantid.api.ExperimentInfo.getInstrumentFilename(self.instrument, self.datadate),DetectorValue='1')
                    configfilestring = mantid.api.ExperimentInfo.getInstrumentFilename(self.instrument, self.datadate)

                self.calibrationtext += "Empty instrument file loaded:  "+configfilestring+"\n"

                #if there is any error in loading the empty instrument, then issue an error.
                if calibration == None:
                    raise RuntimeError("Instrument definition could not be loaded")

                #remove any beam monitors from the empty instrument calibration
                i=0
                while(calibration.getDetector(i).isMonitor()):
                    i += 1
                    #i is the index of the first true detector
                #now, crop the workspace of the monitors
                calibration = CropWorkspace(calibration,StartWorkspaceIndex=i)                
    
                self.calibrationtext += "Loaded empty instrument file.   Calibration set to 1.0 \n"
                print("No calibration file loaded, calibration set to 1.0.")
            else:
                path=GetPathFromRunNumber(self.instrument,self.vanruns[0])
                calibration=Load(Filename=path)
                if self.filterbadpulses:
                    FilterBadPulses(InputWorkspace = 'calibration', OutputWorkspace = 'calibration')
                self.calibrationtext += "Loaded vanadium run from "+path +"\n"
                print(("Calibration file "+path+" loaded."))
                if len(self.vanruns) > 1:
                    for i in range(1,len(self.vanruns)):
                        path = GetPathFromRunNumber(self.instrument,self.vanruns[i])
                        calibrationtemp=Load(Filename=path)
                        if self.filterbadpulses:
                            FilterBadPulses(InputWorkspace = 'calibrationtemp', OutputWorkspace = 'calibrationtemp')
                        calibration=calibration+calibrationtemp
                        self.calibrationtext += "Added vanadium run from "+path +"\n"
                if self.filterbadpulses:
                    self.calibrationtext+= "Bad pulses have been filtered from the vanadium calibration.\n"

                #If the file is not being normalized to one than normalise by current.    
                #Normalize the intensity by charge in milliAmphours (1mA hr = 3.6 coul)
                if self.normalizedcalibration!=True:
                    NormaliseByCurrent(InputWorkspace='calibration', OutputWorkspace='calibration')
                    self.calibrationtext += "Normalized vanadium to proton charge.\n"

                #get the total current and store it in a variable.
                totalmuAhr = mtd['calibration'].run().getProtonCharge()
                totalcoul  = totalmuAhr/1000*3.6
                self.calibrationtext += "Proton charge ("+str(totalmuAhr) + " micro-Ah), ("+str(totalcoul)+" C).\n"


                #Change the units of the calibration workspace
                if self.units != None:
                    ConvertUnits(InputWorkspace = 'calibration', OutputWorkspace = 'calibration', Target = self.units, EMode='Elastic')
                #Need to integrate the data from vanmin to vanmax
                #check the limits.
                if self.vanmin == None:
                    vanmin = mtd['calibration'].getTofMin()
                else:
                    vanmin = self.vanmin
                if self.vanmax == None:
                    vanmax = mtd['calibration'].getTofMax()
                else:
                    vanmax = self.vanmax
                #integrate the data.
                Rebin(InputWorkspace = 'calibration', OutputWorkspace = 'calibration', Params=str(vanmin)+','+str(vanmax-vanmin)+','+str(vanmax),PreserveEvents=False)
                calibration=mtd['calibration']
                self.calibrationtext += "Vanadium integrated between "+ calibration.getAxis(0).getUnit().unitID() + " " +str(vanmin)+", "+str(vanmax)+ " " +calibration.getAxis(0).getUnit().label()+"\n"
            #Do the Masking.
            #loop over the self.mask array
            for elem in self.mask:
                self.calibrationtext += "Masking vanadium with "+str(elem) +"\n"
                if elem['algorithm'].lower()=='banktubepixel':
                    myparser=XMLparser(None)
                    #get the bank tubes and pixels
                    if 'bank' in elem:
                        bank = myparser.parseto('bank',elem['bank'])
                    else:
                        bank = None
                    if 'tube' in elem:
                        #Masking of banks in maskBTP is an integer list, but masking of tubes
                        #is based upon the input being a string.  
                        #tube = myparser.parseto('tube',elem['tube'])
                        tube = elem['tube']
                    else:
                        tube = None
                    if 'pixel' in elem:
                        #Masking of banks in maskBTP is an integer list, but masking of pixels
                        #is based upon the input being a string.  
                        #pixel = myparser.parseto('pixel',elem['pixel'])
                        pixel = elem['pixel'] 
                    else:
                        pixel = None
                    #7/1/2019 - removed the maskbtpold code - switched to standard MaskBTP


                    #MaskBTP(**btpdict)
                    MaskBTP(Workspace='calibration',Bank=bank,Tube=tube,Pixel=pixel)

                if elem['algorithm'].lower()=='angle':
                    #get the ttmin and max
                    if 'twothetamin' in elem:
                        ttmin = float(elem['twothetamin'])
                    else:
                        ttmin = None
                    if 'twothetamax' in elem:
                        ttmax = float(elem['twothetamax'])
                    else:
                        ttmax = None
                    #mask the angle
                    print((ttmin,ttmax))
                    MaskAngleold(Workspace='calibration',twothetamin=ttmin,twothetamax=ttmax)
#                    angledict = dict(elem)
#                    angledict['Workspace'] = 'calibration'
#                    MaskAngle(**angledict)

                if elem['algorithm'].lower() == 'mediandetectortest':
                    detectordiagdict = dict(elem)
                    detectordiagdict['InputWorkspace'] = 'calibration'
                    detectordiagdict['OutputWorkspace']= 'maskdetectordiag'
                    #Need to pop off the 'algorithm' keyvalue from the dict.
                    detectordiagdict.pop('algorithm')
                    MedianDetectorTest(**detectordiagdict)
                    MaskDetectors(Workspace='calibration',MaskedWorkspace='maskdetectordiag')
                    DeleteWorkspace(Workspace='maskdetectordiag')

                if elem['algorithm'].lower() == 'finddetectorsoutsidelimits':
                    detectordiagdict = dict(elem)
                    detectordiagdict['InputWorkspace'] = 'calibration'
                    detectordiagdict['OutputWorkspace']= 'maskdetectordiag'
                    #Need to pop off the 'algorithm' keyvalue from the dict.
                    detectordiagdict.pop('algorithm')
                    FindDetectorsOutsideLimits(**detectordiagdict)
                    MaskDetectors(Workspace='calibration',MaskedWorkspace='maskdetectordiag')
                    DeleteWorkspace(Workspace='maskdetectordiag')


            #Merging a prior mask with the current mask.
            if self.maskfilename != None:
                #check if the file exists
                try:
                    Load(Filename = self.maskfilename, OutputWorkspace = 'loadedmask')
                    MaskDetectors(Workspace='calibration',MaskedWorkspace='loadedmask')
                    DeleteWorkspace(Workspace='loadedmask')
                    self.calibrationtext += "Mask vanadium merged with file " + self.maskfilename+'\n'
                except:
                    raise RuntimeError("Could not load preexisting maskfile " + self.maskfilename)


            #This is the section for normalizing the vanadium calibration to fluctuate about 1.0
            if self.normalizedcalibration==True:
                #normalize the calibration file about 1.0
                #get the mean calibration intensity for non-masked points
                datay = mtd['calibration'].extractY()
                meanval = float(datay[datay>0].mean())
                CreateSingleValuedWorkspace(OutputWorkspace='meanval',DataValue=meanval)
                Divide(LHSWorkspace='calibration',RHSWorkspace='meanval',OutputWorkspace='calibration')
                DeleteWorkspace(Workspace='meanval')
                self.calibrationtext += "Calibration normalized to fluctuate about 1.0\n"      

            #save the calibration file as a loadable workspace.            
            SaveNexus(Filename = os.getcwd()+'/'+self.processedfilename, InputWorkspace = 'calibration')
            changepermissions(os.getcwd()+'/'+self.processedfilename)
            self.calibrationtext += "Saving calibration and mask to "+self.processedfilename
            
            sumfile = open(os.getcwd()+'/'+self.processedfilename+'.calibration_summary.txt', 'w')
            sumfile.write(self.calibrationtext)
            sumfile.close()
            changepermissions(os.getcwd()+'/'+self.processedfilename+'.calibration_summary.txt')


#        print self.calibrationtext


    def CreateFriendlyFilename(self,wsname):

        #access the sample log to generate the friendlyname
        friendlyfilename = os.path.abspath(os.curdir)+'/'+self.friendlyname+'/'+self.friendlyname
        if self.friendlynamelogs != None:
            #get the handle to the run
            run = mtd[wsname].run()
            for part in self.friendlynamelogs if not isinstance(self.friendlynamelogs,str) else [self.friendlynamelogs]:
                if run.hasProperty(part):
                    value = run.getProperty(part).value
                    try:
                        friendlyfilename += "_"+part+"_"+value
                    except:
                        #splitting up the value in case of deciaml points.  Typically, decimal points in filenames
                        #causes issues with Horace.
                        value=array(value).mean()

                        roundedvalue = "%.2f" % value
                        valuestringwithoutdot = str(roundedvalue).replace('.', 'p')
                        friendlyfilename += "_"+valuestringwithoutdot

#                        intpart = int(value)
#                        decpart = int(round(abs(value-intpart)*100))
#                        friendlyfilename += "_"+str(intpart)+"p"+str(decpart)
        return friendlyfilename



def definegoniometer(names, offsets, directions, axes, workspace):
    """
    Function for defining the goniometer
    Returns the psi value which is written to the nxspe file.
    workspace is a string
    returns [psi, text string of goniometer settings]
    """


    #get a handle to the workspace
    ws = mtd[workspace]

    outputstring = ""
    psivalue     = 0


    #Transform all inputs to LISTS.
    if str(type(offsets)) != "<class 'list'>" and str(type(offsets)) != "<class 'NoneType'>":
        offsets = [offsets]

    if str(type(names)) != "<class 'list'>" and str(type(names)) != "<class 'NoneType'>":
        names = [names]

    if str(type(directions)) != "<class 'list'>" and str(type(directions)) != "<class 'NoneType'>":
        directions = [directions]

    if str(type(axes)) != "<class 'list'>" and str(type(axes)) != "<class 'NoneType'>":
        axes = [axes]


    #first case, no motor name given (i.e. None)
    if names == None:
        if offsets == None:
            #no motor name, and no offset, just exit
            outputstring = "No goniometer set.\n"
            return [psivalue, outputstring]
        else:
           #check the number of offsets = N axes = Ndirections
            if (len(offsets)==len(directions) and len(directions)==len(axes) and len(offsets)<7):
                #because no motor names given, we need to make a fake log.
                #for loop over the angles
                #list of 6 empty strings that will be filled in
                anglelist = ["","","","","",""]
                try:
                    for i in range(len(offsets)):
                        AddSampleLog(Workspace=workspace,LogName="angle"+str(i),
                                    LogText=str(offsets[i]),LogType="Number Series")
                        anglelist[i] = "angle"+str(i)+","+axes[i]+","+str(directions[i])
                        #print anglelist[i]
                except:
                    raise RuntimeError("Could not find goniometer axis(axes)")
                SetGoniometer(Workspace=workspace,Axis0=anglelist[0],Axis1=anglelist[1],Axis2=anglelist[2],
                              Axis3=anglelist[3],Axis4=anglelist[4],Axis5=anglelist[5])
                outputstring = "The following axes have been set:\n"

                for i in range(len(offsets)):
                    tempstr = "CCW"
                    if directions[i] == -1:
                        tempstr = "CW"
                    outputstring += "   Axis"+str(i)+" along the "+ axes[i] +" direction, "+tempstr+", rotation angle="+str(offsets[i])+"\n"

                psivalue = offsets[0]
                return [psivalue, outputstring]
            else:
                raise ValueError("Number of angle offsets, directions and axes do not match.")

    #other big case, Motor name is given
    else:
        #Are there any offsets listed
        print(offsets,"****----***")
        if offsets==None:
            #create offsets = 0 for ALL motor names listed
            offsets = zeros(len(names))
        #check if the noffsets = Naxes = Ndirectiosn = Nnames
        if (len(offsets)==len(directions) and len(directions)==len(axes) and len(axes)==len(names) and len(offsets)<7):
            #everything is ready
            anglelist = ["","","","","",""]
            anglevalues = []
            try:
                for i in range(len(offsets)):
                    #get the correct log from the workspace
                    #print(offsets[i],"****-----****")
                    angle = mean(ws.run().get(names[i]).value) + offsets[i]
                    anglevalues.append(angle)
                    AddSampleLog(Workspace=workspace,LogName="angle"+str(i),
                                LogText=str(angle),LogType="Number Series")
                    anglelist[i] = "angle"+str(i)+","+axes[i]+","+str(directions[i])
            except:
                raise RuntimeError("Could not find goniometer axis "+names[i])
            SetGoniometer(Workspace=workspace,Axis0=anglelist[0],Axis1=anglelist[1],Axis2=anglelist[2],
                              Axis3=anglelist[3],Axis4=anglelist[4],Axis5=anglelist[5])
            outputstring = "The following axes have been set:\n"

            for i in range(len(offsets)):
                tempstr = "CCW"
                if directions[i] == -1:
                    tempstr = "CW"
                outputstring += "   Axis"+str(i)+" along the "+ axes[i] +" direction, "+tempstr+", rotation angle="+names[i]+"+" +str(offsets[i])+"="+str(anglevalues[i])+"\n"

            return [anglevalues[0],outputstring]
        else:
            raise ValueError("Number of angle names, offsets, directions and axes do not match.")



#def createanglelist(ws,astep):
#    """
#    Function to create a map of detectors corresponding to angles in a certain range.
#    ws is a string refering to the workspace that contains an instrument
#    astep is the step size for the angular binning (degrees).
#    """
#
#    amin = 0.0
#    amax = 180.0
#    bin_angles=arange(amin+astep*0.5,amax+astep*0.5,astep)
#    a=[[] for i in range(len(bin_angles))] #list of list with detector IDs
#    w=mtd[ws]
#    origin = w.getInstrument().getSample().getPos()
#    for i in range(w.getNumberHistograms()):
#        ang=w.getDetector(i).getTwoTheta(origin,V3D(0,0,1))*180/pi
#        index=int((ang-amin)/astep)
#        if (index>=0) and (index<len(a)) and ((w.getDetector(i).getID())>0):
#            a[index].append(w.getDetector(i).getID())
#    #create lists with angles and detector ID only for bins where there are detectors 
#    ang_list=[]
#    detIDlist=[]
#    for elem,ang in zip(a,bin_angles):
#        if len(elem)>0:
#            detIDlist.append(elem)
#            ang_list.append(ang)
#	# file with grouping information, saved to current directory
#    outputdir = os.path.abspath(os.curdir)
#    f=open(outputdir+"/powdergroup.map",'w')
#    print >>f,len(ang_list)
#    for i in range(len(ang_list)):
#        print >>f,i
#        print >>f,len(detIDlist[i])
#        mystring=str(detIDlist[i]).strip(']').strip('[')
#        mystring=mystring.replace(',','')
#        print >>f,mystring
#    f.close()
#    #print ang_list
#    # par file
#    f=open(outputdir+"/powdergroup.par",'w')
#    print >>f,len(ang_list)
#    for i in range(len(ang_list)):
#        print >>f,3.5,ang_list[i],0.0,1.0,1.0,1
#    f.close()
#    return [ang_list,detIDlist]


def calqrangefromworkspace(workspace):
    """ Direct geometry calculator for determining the max and min
    q-values that will be measured for a given workspace, ws.  The
    workspace must have a log value for incident energy called Ei.
    ws must be a 2D workspace with the x-axis energy.
    THIS IS NOT CHECKED.
    The function uses the incident energy (Ei),
    corresponding to zero energy transfer at scattering angle1
    and hwmin energy transfer at scattering angle2"""

    ws = mtd[workspace]
    
    #Get the Ei from the log
    ei = float(ws.run().get('Ei').value)

    #get angle1 and angle2 (these are the min and max scattering angles).
    angle1 = 180.0
    angle2 = 0.0

    #loop over all the detectors and find the
    for i in range(ws.getNumberHistograms()):
        #get the detectorIDs
        ids = ws.getSpectrum(i).getDetectorIDs()
        for j in ids:
            det = ws.getInstrument().getDetector(j)
            if not det.isMasked():
                angle =  degrees(det.getTwotheta(V3D(0,0,0),V3D(0,0,1)))
                if angle < angle1:
                    angle1 = angle
                if angle > angle2:
                    angle2 = angle
 
    #need the minimum energy transfer the data have been binned to.
    hwmin = ws.readX(0)[0]

    hbarsqrd_over_2m = 2.072

    Qvecmin = sqrt((1.0/hbarsqrd_over_2m)*(2*ei-0.0-2.0*sqrt(ei*(ei-0.0))*cos(angle1*pi/180.0)))
    Qvecmax = sqrt((1.0/hbarsqrd_over_2m)*(2*ei-hwmin-2.0*sqrt(ei*(ei-hwmin))*cos(angle2*pi/180.0)))
        
    return array([Qvecmin, Qvecmax])



def changepermissions(filename):
    """
    change permissions of the directory and file of filename
    to read, write for everyone
    directory also allows execute.
    """
    parentdir = os.path.dirname(filename)
    try:
        os.chmod(parentdir, stat.S_IRWXG | stat.S_IRWXO | stat.S_IRWXU)
    except:
        print(("Not able to change permissions of " + parentdir))

    try:
        os.chmod(filename, stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH |stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH )
    except:
       print(("Not able to change permissions of " + filename))


def parsedbl(val):
	"""
	####PARSEDBL##########################################
	#This function parses the string and converts to a double floating
	#point value
	#
	# example:
	#  print parsedbl("123.7")
	#	  123.7
	#
	######################################################
	"""
	return float(val)

def parsebool(val):
	"""
	####PARSEBOOL##########################################
	#This function parses text to see if it is a true statement.
	#It will return true for "yes", "true", "t", "1", "tru", "tr", "y".
	#can be single or double quotes.
	#
	#Returns True for all cases listed here, independent of case (upper/lower/mixed)
	#Returns False for all other cases.
	#
	# example:
	#  print parsebool("y")
	#	 True
	#  print parsebool("test")
	#	 False
	#
	#######################################################
	"""
	return val.lower() in ("yes", "true", "t", "1", "tru", "tr", "y")



def ShiftTime(WName,lg_name):
	"""
	shift the time in a given log (lg_name) of a workspace (WName) to match the time in the proton charge log"
	"""
	H_IN = mtd[WName]
	PC =  H_IN.getRun()['proton_charge'].firstTime()
	#print "P="+str(PC)+"\n"
	P =  H_IN.getRun()[lg_name].firstTime()
	#print "P="+str(P)+"\n"
	Tdiff = PC-P
	Tdiff_num = Tdiff.total_milliseconds()*1E-3
	#print "Tdiff="+str(Tdiff_num)+"\n"
	ChangeLogTime(InputWorkspace=WName, OutputWorkspace = WName, LogName = lg_name, TimeOffset = Tdiff_num)




if __name__ == "__main__":
    #check number of arguments
    if (len(sys.argv) != 2): 
        print("reduction code requires a datatext file")
        sys.exit()
    if not(os.path.isfile(sys.argv[1])):
        print(("data text file ", sys.argv[1], " not found"))
        sys.exit()
    dgsreduction(XMLfile=sys.argv[1])
