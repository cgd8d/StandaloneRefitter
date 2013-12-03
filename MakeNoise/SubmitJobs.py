#!/usr/local/bin/python

EntriesPerRun = 10
Queue = 'medium'

# Inclusive run ranges.
# This set of windows is meant to only include Run2ab for now.
# Boundaries which are explicitly commented are confirmed exactly;
# others are not independently confirmed by me with environmental correlations.
# (I still trust them, though.)
RunWindows = [(2424, 2699),
              (2700, 2852), # Ebox 1 fan installed
              (2853, 2891), # Ebox 2 fan installed
              (2892, 3117), # Power outage here.
              (3118, 3326), # APD board swap
              (3327, 3700), # Utility power swap
              (3701, 3949),
              (3950, 4140), # Ralph's diode box installed
              (4141, 4579),
              (4580, 4779),
              (4780, 5197), # LC filters removed from FECs
              (5198, 5892)] # Run2c ends.

import subprocess
import ROOT
ROOT.gSystem.Load("libEXOUtilities")

for runWindow in RunWindows:
    ds = ROOT.EXORunInfoManager.GetDataSet('Data/Processed/masked', '%i<=run&&run<=%i&&quality==\"GOLDEN\"&&runType==\"Data-Physics\"' % runWindow)
    runList = [str(ri.GetRunNumber()) for ri in ds]
    subprocess.call(['bsub', '-q', Queue, '-R', 'rhel50', '-o', '%i_to_%i.log' % runWindow,
                     './MakeNoiseFile', '%i_to_%i.dat' % runWindow, str(EntriesPerRun)] + runList)

# We also generate a noise window for runs 2401-2423 (09-28-11 APD biases).
# But, lacking physics data, we use a noise run there.
subprocess.call(['bsub', '-q', Queue, '-R', 'rhel50', '-o', '2401_to_2423.log',
                 './MakeNoiseFile', '2401_to_2423.dat', '10000', '2401'])

