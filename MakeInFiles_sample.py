"""
Script to generate a set of job folders, each containing some number of job scripts.
Make sure it can be run at SLAC, since our NERSC ROOT does not have xml enabled.
We put the input files into a hierarchy rooted at pwd.

ToDO:
Extract the number of entries in the processed run file from metadata.
Verify that the processed file I see is truly available -- don't want to crash later for missing file.
"""
from __future__ import with_statement
import os
import glob
import ROOT
ROOT.gROOT.SetBatch()
ROOT.gSystem.Load("libEXOUtilities")

DenoisedOutDir = "/scratch1/scratchdirs/claytond/LightOnly"
ProcsPerJob = [40]*1000 # 40 simultaneous jobs is about all that we can safely handle.

NoiseFileBase = "/global/u1/c/claytond/NoiseCorrFiles"
RunWindows = [(2401, 2423), (2424, 2690), (2691, 2852), (2853, 2891), (2892, 3117), (3118, 3326), (3327, 3700),
              (3701, 3949), (3950, 4140), (4141, 4579), (4580, 4779), (4780, 5197), (5198, 5590), (5591, 5892)]
def GetNoiseFile(runNo):
    noiseFile = ""
    for noiseWindow in RunWindows:
        if runNo >= noiseWindow[0] and runNo <= noiseWindow[1]:
            return NoiseFileBase + "/%i_to_%i.dat" % noiseWindow
    # Currently I refuse to pick a noise file by any heuristic here.
    raise ValueError("No noise file for run %i." % runNo)

# RunTypes are "Data-Source calibration", "Data-Physics", etc.
ProcDataset = ROOT.EXORunInfoManager.GetDataSet("Data/Processed/masked",
                                                "run>=2464&&run<=5367&&runType==\"Data-Source calibration\"")

def JobForProc(procFile, runNo):
    FilePartsFull = procFile.split('/')
    FileParts = FilePartsFull[-6:]
    if FilePartsFull[-3] == "processed":
        # Check for the existence of the masked file.  (We may not have permission to check for the proc file.)
        FileToVerify_Parts = list(FilePartsFull)
        FileToVerify_Parts[-3] = "masked"
        FileToVerify_Parts[-1] = "masked" + FileToVerify_Parts[-1][-17:]
        FileToVerify = '/' + '/'.join(FileToVerify_Parts)
    else:
        FileToVerify = procFile
    if not os.path.isfile(FileToVerify):
        raise OSError("File %s does not exist. (Processing must have failed.)" % FileToVerify)
    FileBase = FileParts[-1][-17:]
    xrootd_procfile = '/' + '/'.join(FileParts)
    FileParts[-3] = 'root'
    FileParts[-1] = 'run' + FileBase
    xrootd_rawfile = '/' + '/'.join(FileParts)
    return ("%s\n%s\n%s\n%s\n%i\n%i\n%f\n%f\n" %
            (xrootd_procfile, # Processed file
             xrootd_rawfile, # Raw file
             "%s/%i/denoised%s" % (DenoisedOutDir, runNo, FileBase), # Out file
             GetNoiseFile(runNo), # noise file
             0, -1, 0.1, 1.6)) # Run parameters

OutRunList = []
ProcList = []
for runInfo in ProcDataset:
    try:
        ProcsToInsert = []
        for runFile in runInfo.GetRunFiles():
            Job = JobForProc(runFile.GetFileLocation(), runInfo.GetRunNumber())
            numEntries = runFile.FindMetaData("eventCount").AsInt()
            ProcsToInsert.append((numEntries, Job))
    except (ValueError, OSError), exc:
        print "Failed to add jobs for run %i: %s." % (runInfo.GetRunNumber(), exc.message)
    else:
        ProcList += ProcsToInsert
        OutRunList.append(os.path.join(DenoisedOutDir, str(runInfo.GetRunNumber())))
        print "Added jobs for run %i." % runInfo.GetRunNumber()

# Try to group processes of similar length together, so they finish in similar times.
# Also make the shortest jobs go into the smallest group, since this should also reduce waste.
ProcList.sort(key = lambda x: x[0])
JobIndex = 0
while len(ProcList) > 0:
    ProcsInThisJob = []
    NumProcsInThisJob = ProcsPerJob.pop(0)
    for i in xrange(NumProcsInThisJob):
        if len(ProcList) > 0: ProcsInThisJob.append(ProcList.pop())
    try:
        os.mkdir('Job%04i' % JobIndex)
    except:
        for oldfile in glob.glob('Job%04i' % JobIndex): os.remove(oldfile)
    with open('Job%04i/LengthOfJob.txt' % JobIndex, 'w') as LengthFile:
        LengthFile.write("%i entries in the longest job." % ProcsInThisJob[0][0])
    for i in xrange(len(ProcsInThisJob)):
        with open('Job%04i/infile%04i.txt' % (JobIndex, i), 'w') as infile:
            infile.write(ProcsInThisJob[i][1])
    JobIndex += 1

# We run this script at SLAC, and the necessary output folders may not exist -- print the command to run.
print 'To create the necessary output folders, run \'bash MakeOutDirs.sh\''
with open('MakeOutDirs.sh', 'w') as make_out_dirs:
    make_out_dirs.write('#!/bin/bash\n')
    make_out_dirs.write('mkdir -p ' + ' '.join(OutRunList))
