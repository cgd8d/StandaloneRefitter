
import glob, os

InOutDir = "/scratch1/scratchdirs/claytond/LightOnly/InOutFiles"
DenoisedOutDir = "/scratch1/scratchdirs/claytond/LightOnly"

NoiseFileBase = "/global/u1/c/claytond/NoiseCorrFiles"
RunWindows = [(2464, 2699), (2700, 2852), (2853, 2891), (2892, 3117), (3118, 3329), (3330, 3699),
              (3700, 3949), (3950, 4149), (4150, 4579), (4580, 4779), (4780, 5197), (5198, 5367)]

import glob
rawFiles = glob.glob("/project/projectdirs/exo200/cal_data/raw/*/run*.root")
procFiles = glob.glob("/project/projectdirs/exo200/cal_data/processed/*/proc*.root")

RawProcNoiseGroups = []
for rawFile in rawFiles:
    tail = rawFile[-17:]
    runNo = int(tail[:8])

    noiseFile = ""
    for noiseWindow in RunWindows:
        if runNo >= noiseWindow[0] and runNo <= noiseWindow[1]:
            noiseFile = NoiseFileBase + "/%i_to_%i.dat" % noiseWindow
            break
    if noiseFile == "":
        if runNo < RunWindows[0][0]: noiseFile = NoiseFileBase + "/%i_to_%i.dat" % RunWindows[0]
        elif runNo > RunWindows[-1][1]: noiseFile = NoiseFileBase + "/%i_to_%i.dat" % RunWindows[-1]
    if noiseFile == "":
        print "Run %i is in a noise gap." % runNo
        continue

    for procFile in procFiles:
        if procFile[-17:] == tail:
            RawProcNoiseGroups.append((rawFile, procFile, noiseFile))
            break

# Sort by size of processed file; this helps us pair short file segments together, which will finish faster.
RawProcNoiseGroups.sort(key = lambda x: os.path.getsize(x[1]))

NodeNum = -1
while len(RawProcNoiseGroups) > 0:
    NodeNum += 1
    proc1 = RawProcNoiseGroups.pop()
    if len(RawProcNoiseGroups) > 0: proc2 = RawProcNoiseGroups.pop()
    else: proc2 = None
    try:
        os.mkdir("%s/Node%03i" % (InOutDir, NodeNum))
    except OSError: pass
    with open("%s/Node%03i/infile0000.txt" % (InOutDir, NodeNum), 'w') as infile1:
        infile1.write("%s\n%s\n%s\n%s\n%i\n%i\n%f\n" %
                      (proc1[1], proc1[0],
                       DenoisedOutDir + "/denoised" + proc1[1][-17:],
                       proc1[2],
                       0, 1000000, 0.1))
    if proc2 == None: continue
    with open("%s/Node%03i/infile0001.txt" % (InOutDir, NodeNum), 'w') as infile2:
        infile2.write("%s\n%s\n%s\n%s\n%i\n%i\n%f\n" %
                      (proc2[1], proc2[0],
                       DenoisedOutDir + "/denoised" + proc2[1][-17:],
                       proc2[2],
                       0, 1000000, 0.1))

