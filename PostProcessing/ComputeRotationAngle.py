"""
For a given theta, find the resolution; then optimize resolution by varying theta.
This script is adapted by Clayton from a script by Manuel,
at http://java.freehep.org/svn/repos/exo/list/EXOCalibration/trunk/?revision=6540.
It is only intended to run on data, not MC.  (Different thresholds are not taken into account.)
It has its own rudimentary peak-finder, so it's intended to be a little robust against changes to the detector.
(That's not to say it can't fail.  Please report any unsound output it may produce.
 Note that in addition to the output values, it produces two pdf files with fits and so forth.
 This can help to pinpoint a problem.)
"""

import RotationAngle_common

def Run(prefix, **kwargs):
    RotationAngle_common.prefix = prefix
    RotationAngle_common.WhichScint = ROOT.kScintRawEnergy

    return RotationAngle_common.Run(**kwargs)

if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2:
        print 'We need the run number to be passed in as an argument:'
        print 'python ComputeRotationAngle.py <fileglob>'
        sys.exit()
    import ROOT
    ROOT.gROOT.SetBatch()
    if ROOT.gSystem.Load("libEXOUtilities") < 0: sys.exit('Failed to load EXOUtilities.')
    if ROOT.gSystem.Load("/nfs/slac/g/exo_data4/users/cgd8d/rhel5-64/EXO_Fitting/EXO_Fitting/lib/libEXOFitting") < 0: sys.exit('Failed to load EXOFitting.')
    EventTree = ROOT.TChain('tree')
    EventTree.Add(sys.argv[1])
    AllMaskedFiles = glob.glob(sys.argv[1])
    print "Using files:", AllMaskedFiles
    AllMaskedFiles.sort(reverse = True)
    LastFile = ROOT.TFile(AllMaskedFiles[0])
    LastTree = LastFile.Get('tree')
    ControlRecordList = LastTree.GetUserInfo().At(1)
    result = Run('RotationAngle', EventTree=EventTree, ControlRecordList=ControlRecordList)
    print result
