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
ROOT.gSystem.Load("libEXOUtilities")

ProcsPerJob = 100
DenoisedOutDir = "/scratch1/scratchdirs/claytond/LightOnly"

NoiseFileBase = "/global/u1/c/claytond/NoiseCorrFiles"
RunWindows = [(2464, 2699), (2700, 2852), (2853, 2891), (2892, 3117), (3118, 3329), (3330, 3699),
              (3700, 3949), (3950, 4149), (4150, 4579), (4580, 4779), (4780, 5197), (5198, 5367)]
def GetNoiseFile(runNo):
    noiseFile = ""
    for noiseWindow in RunWindows:
        if runNo >= noiseWindow[0] and runNo <= noiseWindow[1]:
            return NoiseFileBase + "/%i_to_%i.dat" % noiseWindow
    # Currently I refuse to pick a noise file by any heuristic here.
    raise ValueError("No noise file for run %i." % runNo)

ProcDataset = ROOT.EXORunInfoManager.GetDataSet("Data/Processed/masked",
                                                "run>=2464&&run<=5367")

def JobForProc(procFile, runNo):
    FileParts = procFile.split('/')[-6:]
    FileBase = FileParts[-1][-17:]
    RawFileParts = FileParts
    RawFileParts[3] = 'root'
    RawFileParts[-1] = 'run' + FileBase
    return ("%s\n%s\n%s\n%s\n%i\n%i\n%f\n" %
            ('/' + '/'.join(FileParts), # Processed file
             '/' + '/'.join(RawFileParts), # Raw file
             "%s/%i/denoised%s" % (DenoisedOutDir, runNo, FileBase), # Out file
             GetNoiseFile(runNo), # noise file
             0, -1, 0.1)) # Run parameters

OutRunList = []
ProcList = []
for runInfo in ProcDataset:
    try:
        ProcsToInsert = []
        for runFile in runInfo.GetRunFiles():
            Job = JobForProc(runFile.GetFileLocation(), runInfo.GetRunNumber())
            numEntries = runFile.FindMetaData("eventCount").AsInt()
            ProcsToInsert.append((numEntries, Job))
        ProcList += ProcsToInsert
        print "Added jobs for run %i." % runInfo.GetRunNumber()
    except ValueError, exc:
        print "Failed to add jobs for run %i: %s." % (runInfo.GetRunNumber(), exc.message)

# Try to group processes of similar length together, so they finish in similar times.
# Also make the shortest jobs go into the smallest group, since this should also reduce waste.
ProcList.sort(key = lambda x: x[0])
JobIndex = 0
while len(ProcList) > 0:
    ProcsInThisJob = []
    for i in xrange(ProcsPerJob):
        if len(ProcList) > 0: ProcsInThisJob.append(ProcList.pop()[1])
    try:
        os.mkdir('Job%04i' % JobIndex)
    except:
        for oldfile in glob.glob('Job%04i' % JobIndex): os.remove(oldfile)
    for i in xrange(len(ProcsInThisJob)):
        with open('Job%04i/infile%04i.txt' % (JobIndex, i), 'w') as infile:
            infile.write(ProcsInThisJob[i])
    JobIndex += 1

