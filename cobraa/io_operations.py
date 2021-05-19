from .load import *
from .globals import *
from ROOT import TChain, gROOT, TFile


# The purpose of this class is to handle the input/ouput operations of
# Watchcore (WC). This includes creating directories and files for the different
# operations of WC, such as creating macros, jobs and so forth.

def testCreateDirectory(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
    else:
        rmtree(directory)
        os.makedirs(directory)

def deleteDirectory(directory):
    if os.path.exists(directory):
        rmtree(directory)

def testCreateDirectoryIfNotExist(directory):

    if os.path.exists(directory):
        print('''There is already a directory here. %s
No new directory has been made. Bad idea. Consider saving current files
and using --force.        \n'''%(directory))
    if not os.path.exists(directory):
        os.makedirs(directory)

def generateMacros():

    # Writes the macros for each element (_element) of each physical process (_p)
    # and location (_loc) to be simulated
    
    # Calls macroGenerator(), which provides the macro content 
   
    testCreateDirectoryIfNotExist("mac")

    # write the macros for each type of event to be generated
    # and its location in the detector
    for _p in proc:
        for _loc in proc[_p]:
            for _element in d[_p]:
                if 'singles' not in _p:
                    header,processors,generator,detectorvolume = macroGenerator(_loc,_element,_p,nruns)
                    outfile = open("mac/phys_%s.mac"%(_element).replace(" ",""),"w+")
                    outfile.writelines(generator)
                    outfile.close
                    outfile = open("mac/geo_%s.mac"%(_loc),"w+")
                    outfile.writelines(detectorvolume)
                    outfile.close
                    print(_p,_loc,_element)

    # write the macros for the detector geometry and rat processors
    header,processors,generator,detectorvolume = macroGenerator("","","",0)
    outfile = open(f"mac/detector_{detectorStr}.mac","w+")
    outfile.writelines(header)
    outfile.close
    outfile = open("mac/process.mac","w+")
    outfile.writelines(processors)
    outfile.close

    # write the macros for the number of events to be simulated and
    # expected total event rates in the detector (before detector effects)
    for _k in rates:
        _events = int(float(arguments['-e'])*rates[_k][1])
        if 'singles' in _k:
            print("\n\n\n Warning - only %f days of singles events will be simulated!!!!!\n\n\n"%(_events*nsetSingles*nruns/float(singlespersec*86400)))
            outfile = open(f"mac/evts_singles.mac","w+")
            outfile.writelines(f"/run/beamOn {_events}")
        elif 'pn_ibd' in _k or 'A_Z' in _k or 'fast' in _k:
            outfile = open(f"mac/rates_{_k}.mac","w+")
            outfile.writelines(f"/generator/rate/set {rates[_k][0]}")
            outfile.close
            outfile = open(f"mac/evts_{_k}.mac","w+")
            outfile.writelines(f"/run/beamOn {int(_events)}")
            outfile.close
        else:
            outfile = open(f"mac/rates_{_k}.mac","w+")
            outfile.writelines(f"/generator/rate/set {rates[_k][0]}")
            outfile.close



def generateJobs():


    # Create the directory trees
    for _p in proc:
        for _loc in proc[_p]:
            for _element in d[_p]:
                if 'pn_ibd' in _p or 'A_Z' in _p or 'fast' in _p or 'singles' in _p:
                
                    dir = "root_files%s/%s_%s_%s"%(additionalString,_element,_loc,_p)
                    dir = dir.replace(" ","")
                    if arguments['--force']:
                        print('Using force to recreate dir:',dir)
                        testCreateDirectory(dir)
                    else:
                        testCreateDirectoryIfNotExist(dir)
                    dir = "fred_root_files%s/%s_%s_%s"%(additionalString,_element,_loc,_p)
                    dir = dir.replace(" ","")
                    if arguments['--force']:
                        print('Using force to recreate dir:',dir)
                        testCreateDirectory(dir)
                    else:
                       testCreateDirectoryIfNotExist(dir)
                    dir = "log%s/%s_%s_%s"%(additionalString,_element,_loc,_p)
                    dir = dir.replace(" ","")
                    if arguments['--force']:
                        print('Using force to recreate dir:',dir)
                        testCreateDirectory(dir)
                    else:
                        testCreateDirectoryIfNotExist(dir)

    ratDir      = os.environ['RATROOT']
    nameJob     = "nameJob"
    timeJob     = arguments["--jobTime"]
    outFile     = "out_file.log"
    errFile     = "err_file.log"

    directory   = os.getcwd()
    dir =  "%s/job"%(directory)
    if arguments['--force']:
        print('Using force to recreate dir:',dir)
        testCreateDirectory(dir)
    else:
        testCreateDirectoryIfNotExist(dir)

    # write out the scripts to be run
    # and the jobs which will run the scripts the correct number of times
    # jobs are set up to run locally by default but can be tailored to 
    # job submission on a cluster (currently --lassen and --sheff)
    singlesscript = f"{dir}/script{additionalString}_singles.sh".replace(" ","")
    outfile_singlesscript = open(singlesscript, "w+")
    outfile_singlesscript.writelines(f"""#!/bin/sh
source {ratDir+'/../../env.sh'} && TMPNAME=$(date +%s%N)  && rat mac/detector_{detectorStr}.mac mac/process.mac """)
    for _p in proc:
        for _loc in proc[_p]:
            for _element in d[_p]:
                _element = _element.replace(" ","")
                if 'NA' in _p or 'STEEL' in _p or 'RADIOGENIC' in _p: 
                 outfile_singlesscript.writelines(f" mac/phys_{_element}.mac mac/geo_{_loc}.mac mac/rates_{_element}_{_loc}_{_p}.mac") 
                elif 'pn_ibd' in _p or 'A_Z' in _p or 'FASTNEUTRONS' in _p:
                    script = f"{dir}/script{additionalString}_{_element}_{_loc}_{_p}.sh".replace(" ","")
                    outfile_script = open(script,"w+")
                    outfile_script.writelines(f"""#!/bin/sh
source {ratDir+'/../../env.sh'} && TMPNAME=$(date +%s%N)  && rat mac/detector_{detectorStr}.mac mac/process.mac mac/phys_{_element}.mac mac/geo_{_loc}.mac mac/rates_{_element}_{_loc}_{_p}.mac mac/evts_{_element}_{_loc}_{_p}.mac -o root_files{additionalString}/{_element}_{_loc}_{_p}/run$TMPNAME.root -l log{additionalString}/{_element}_{_loc}_{_p}/run$TMPNAME.log""")
                    outfile_script.close
                    os.chmod(script,S_IRWXU)
                    file = f"{dir}/job{additionalString}_{_element}_{_loc}_{_p}.sh".replace(" ","") 
                    outfile_jobs = open(file,"w+")
                    jobheader = jobSubmissionCommands(_element,timeJob,file,outFile,errFile,script,arguments)
                    outfile_jobs.writelines(jobheader)

                    outfile_jobs.close
                elif _element in d['singles']:
                    for i in range(nsetSingles):
                        file = f"{dir}/job{additionalString}_{_element}_{_loc}_{_p}_{i}.sh".replace(" ","")
                        outfile_jobs = open(file,"w+")
                        jobheader = jobSubmissionCommands(_element,timeJob,file,outFile,errFile,singlesscript,arguments)
                        outfile_jobs.writelines(jobheader)
    
    outfile_singlesscript.writelines(f" mac/evts_singles.mac -o root_files{additionalString}/singles_SINGLES_singles/run$TMPNAME.root -l log{additionalString}/singles_SINGLES_singles/run$TMPNAME.log")
    outfile_singlesscript.close()
    os.chmod(singlesscript,S_IRWXU)


def deleteAllWorkDirectories():

    # does what it says on the tin

    dir = "log"
    if os.path.exists(dir):
        rmtree(dir)

    dir = "jobs"
    if os.path.exists(dir):
        rmtree(dir)

    dir = "mac"
    if os.path.exists(dir):
        rmtree(dir)

    if os.path.exists('fit_param.dat'):
        os.remove('fit_param.dat')

    if os.path.exists('like.bin'):
        os.remove('like.bin')

    if os.path.exists('job'):
        os.remove('job')


def mergeRootFiles():

    # merges the CoRe root files
    # use --mergeRATFiles to also merge the raw root files

    for _p in proc:
        for _loc in proc[_p]:
            for _element in d[_p]:
                _p = _p.replace(" ","")
                print("Generating jobs:",_p,_loc,_element)
                outfile = "root_files%s/merged_%s_%s_%s_%s.root"%(additionalString,_element,_loc,_p)
                outfile = outfile.replace(" ","")
                files = "root_files%s/%s_%s_%s_%s/run*.root"%(additionalString,_element,_loc,_p)
                files = files.replace(" ","")
                # merge the raw root files if required
                if arguments['--mergeRATFiles']:
                    os.system(f'hadd -f -k -v 0 {outfile} {files}')
                #otherwise merge the bonsai root files
                else:
                    os.system(f'hadd -f -k -v 0 core_{outfile} core_{files}')


def macroGenerator(location,element,_dict,nruns):

    # Generates the text to go inside the rat macros
    # and is called by generateMacros()

    # First the header and processor macros;
    # these macros are the same for all jobs

    dir = os.getcwd()
    detectorOption = ""
    detectorvolume = ''

    AO = additionalMacOpt.splitlines()
    for _str in AO:
        if "material" in _str or "optic" in _str or "pmt_model" in _str or "light_cone" in _str or "mu_metal" in _str:
            detectorOption = detectorOption + _str + "\n"
    depth = float(arguments["--depth"])
    rate = 1.0
    header = f"""
/glg4debug/glg4param omit_muon_processes  0.0
/glg4debug/glg4param omit_hadronic_processes  0.0

/rat/db/set DETECTOR experiment "{detectorStr}"
/rat/db/set DETECTOR geo_file "{detectorStr}/{detectorStr}.geo"
{additionalMacOpt}


/run/initialize
"""
    processors=f"""# BEGIN EVENT LOOP
/rat/proc lesssimpledaq
# /rat/proc fitbonsai
#/rat/proc fitcentroid
#/rat/proc fitpath
/rat/proc count
/rat/procset update 200

# Use IO.default_output_filename
/rat/proclast outroot
#END EVENT LOOP
""" 

    # Then the generator (phys) and location (geo) macros;
    # these set the generator, generator conditions and location for a given event type

    if element in d['CHAIN_238U_NA'] or element in d['CHAIN_232Th_NA'] or element in d['40K_NA'] or element in d['STEEL_ACTIVITY'] or element in d['CHAIN_235U_NA']:
        if location == 'PMT':
            generator = f'''
/generator/add decaychain {element}:regexfill:poisson
'''
            detectorvolume = f'''/generator/pos/set inner_pmts_body_phys+
'''
        elif location == "PSUP":
            generator = f"""
/generator/add decaychain {element}:regexfill:poisson
"""
            detectorvolume = f"""
/generator/pos/set   psup
"""
        else:
            locat = location.lower()
            if locat == 'liquid' or locat == 'gd':
                locat = 'detector'
                xTimes = 2
            else:
                xTimes = 10
            generator = f'''
/generator/add decaychain {element}:regexfill:poisson
'''
            detectorvolume = f'''
/generator/pos/set {locat}+
'''

    elif element in d['pn_ibd']:
        generator = f'''
/generator/add combo ibd:regexfill:poisson
/generator/ibd/spectrum {element}
/generator/vtx/set IBD 1 0 0
'''
        detectorvolume = f'''/generator/vtx/set {element} 1 0 0
/generator/pos/set detector+
'''

    elif element in d['A_Z']:
        #A =  int(int(element)/1000)
        #Z = int(element) - A*1000
        generator = f'''
/generator/add {element}:regexfill:poisson

'''
        detectorvolume = f'''/generator/pos/set detector+
'''

    elif element in d['RADIOGENIC']:
        locat = location.lower()
        generator =f'''
/generator/add combo spectrum:regexfill:poisson
/generator/vtx/set neutron {element}
'''
        detectorvolume = f'''/generator/pos/set {locat}
'''

    elif element in d['FASTNEUTRONS']:
        locat = location.lower()
        generator =f'''
/generator/add combo fastneutron:regexfill:poisson
/generator/vtx/set 0 0 0
/generator/fastneutron/depth 2805.0
/generator/fastneutron/enthresh 10.0
/generator/fastneutron/sidewalls 1.0 
'''
        detectorvolume = f'''/generator/pos/set {locat}
'''


    else:
        print('Could not find ',element,location.lower())
        generator = ''
        detectorvolume = ''

    return header,processors,generator,detectorvolume


# Specify the header for the job submission script
def jobSubmissionCommands(_element,timeJob,file,outFile,errFile,script,arguments):
    jobheader=""
    if arguments['--cluster']=='lassen':
        jobheader = f"""#!/bin/sh
#BSUB -nnodes 1  
#BSUB -J job_{_element}    #name of job
#BSUB -W {timeJob}      #time in minute
#BSUB -G ait         # sets bank account
#BSUB -q pbatch         #pool
#BSUB -o {file+outFile}
#BSUB -e {file+errFile}
#BSUB                     # no more psub commands

jsrun -p{nruns} {script}

""" 

    elif arguments['--cluster']=='sheffield':
        jobheader = f"""#!/bin/sh

executable = {script}
output     = {file+outFile}
error      = {file+errFile}
getenv     = True

queue {nruns} {script}

"""

    elif arguments['--cluster']=='glasgow':
        jobheader = f"""#!/bin/sh

for i in `seq {nruns}`; do source {script}; done

    """ 

    elif arguments['--cluster']=='--edinburgh':
        jobheader = f"""#!/bin/sh

for i in `seq {nruns}`; do source {script}; done

    """ 

    else:
        jobheader = f"""#!/bin/sh

for i in `seq {nruns}`; do source {script}; done

    """ 

    return jobheader
