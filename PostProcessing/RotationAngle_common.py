"""
For a given theta, find the resolution; then optimize resolution by varying theta.
This script is adapted by Clayton from a script by Manuel,
at http://java.freehep.org/svn/repos/exo/list/EXOCalibration/trunk/?revision=6540.
It is only intended to run on data, not MC.  (Different thresholds are not taken into account.)
It has its own rudimentary peak-finder, so it's intended to be a little robust against changes to the detector.
(That's not to say it can't fail.  Please report any unsound output it may produce.
 Note that in addition to the output values, it produces two pdf files with fits and so forth.
 This can help to pinpoint a problem.)
This module is meant to be used by others, with slightly different options; this allows easy comparison
of the resolution produced by slightly tweaking options.
For examples, see the ComputeRotationAngle*.py scripts.
Be sure to set the module attributes prefix and WhichScint.
"""

import numpy
import ROOT
import math
import sys
import os

MaxRotatedEnergy = 10000.

canvas = None
FileNames = { }

class PeakFinderError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

def IsFiducial(cluster):
    if abs(cluster.fZ) < 20*ROOT.CLHEP.millimeter or abs(cluster.fZ) > 172*ROOT.CLHEP.millimeter: return False
    if pow(cluster.fX, 2) + pow(cluster.fY, 2) > pow(163*ROOT.CLHEP.millimeter, 2): return False
    return True

def IsSingleSite(scint):
    if scint.GetNumChargeClusters() != 1: return False
    cluster = scint.GetChargeClusterAt(0)
    UniqueU = set([cluster.GetUWireSignalAt(i).fChannel for i in range(cluster.GetNumUWireSignals())])
    if len(UniqueU) > 2: return False
    else: return True

def DivideWithErrors(x, y):
    # If x and y have uncorrelated errors, return their quotient (with errors).
    return (x[0]/y[0],
            math.sqrt(pow(x[1]/y[0], 2) +
                      pow(x[0]*y[1] / pow(y[0], 2), 2)))

def FindPeak(dataSet):
    # Given an unbinned RooDataSet, locate the fit range for the photopeak.
    # OK -- have to find peaks near the end.  First bin it coursely.
    hist = ROOT.TH1F('PeakFinder', 'PeakFinder', int(MaxRotatedEnergy/200), 0., MaxRotatedEnergy)
    if dataSet.fillHistogram(hist, ROOT.RooArgList(dataSet.get())) == None:
        raise PeakFinderError('Failed to bin data for peak-finding.')
    for i in range(hist.GetNbinsX()-1, 0, -1): # Step from right to left.
        # We want the chance of an accidental peak to be very small.
        # Here we demand it's more than 3 sigmas.
        # (4 occasionally misses real scint-only peaks, which are hard because of the low resolution.)
        if hist.GetBinContent(i+1) <= 0: continue
        sigma = math.sqrt(hist.GetBinContent(i) + hist.GetBinContent(i+1))
        if (hist.GetBinContent(i+1) - hist.GetBinContent(i))/sigma >= 3:
            Peak = hist.GetBinCenter(i+1)
            break
    try:
        if i == 1: raise PeakFinderError('Failed to find a peak.')
        if Peak < 300: raise PeakFinderError('Peak found, but the energy is only ' + str(Peak) + '.')
        return Peak
    except PeakFinderError:
        canvas = ROOT.TCanvas()
        hist.Draw()
        canvas.Print(prefix + 'FailedToFindPeak.pdf')
        raise

def DoFit(dataSet, PeakPosGuess, SigmaGuess, ErfcFracGuess, ID = None):
    # Using input guesses, pick a fit window and refine those guesses.
    # If ID is given, use it to save a plot.
    RotatedEnergy = dataSet.get()['RotatedEnergy']

    # Note, the bounds used here imply resolution can never be worse than 30% at the peak.
    mean = ROOT.RooRealVar('Mean', 'Mean', PeakPosGuess, PeakPosGuess - SigmaGuess, PeakPosGuess + SigmaGuess)
    sigma = ROOT.RooRealVar('Sigma', 'Sigma', SigmaGuess, 0., 0.3*PeakPosGuess)
    erfc_coeff = ROOT.RooRealVar('coeff', 'Erfc coeff', ErfcFracGuess, 0., 1.)
    GaussianComponent = ROOT.RooGaussian('GaussianComponent', 'GaussianComponent',
                                         RotatedEnergy, mean, sigma)
    ErfcComponent = ROOT.RooGenericPdf('ErfcComponent', 'ErfcComponent',
                                       'TMath::Erfc((@0 - @1) / (TMath::Sqrt(2)*@2))',
                                       ROOT.RooArgList(RotatedEnergy, mean, sigma))
    FitFunction = ROOT.RooAddPdf('FitFunction', 'FitFunction',
                                 GaussianComponent, ErfcComponent, erfc_coeff)
    FitFunction.fitTo(dataSet, ROOT.RooFit.Range(PeakPosGuess - 3*SigmaGuess, PeakPosGuess + 2.5*SigmaGuess))

    if ID:
        # ID[0] should identify the fit as single-site or multi-site.
        # ID[1] should be a formatted identifier of theta.

        # The plot should be wider than the fit window, so we can see what wasn't fit.
        frame = RotatedEnergy.frame(ROOT.RooFit.Bins(200),
                                    ROOT.RooFit.Range(PeakPosGuess - 6*SigmaGuess,
                                                      PeakPosGuess + 6*SigmaGuess))
        dataSet.plotOn(frame,
                       ROOT.RooFit.DrawOption('ZP'),
                       ROOT.RooFit.MarkerStyle(20),
                       ROOT.RooFit.MarkerSize(0.5))
        FitFunction.plotOn(frame)
        FitFunction.plotOn(frame,
                           ROOT.RooFit.Components("ErfcComponent"),
                           ROOT.RooFit.LineWidth(2), ROOT.RooFit.LineStyle(2))
        FitFunction.plotOn(frame,
                           ROOT.RooFit.Components("GaussianComponent"),
                           ROOT.RooFit.LineWidth(2), ROOT.RooFit.LineStyle(2))
        FitFunction.paramOn(frame,
                            ROOT.RooFit.Layout(0.6))
        frame.SetTitle(ID[0] + ' rotated spectrum (Theta = ' + ID[1] + ' rad)')
        frame.SetXTitle('Rotated Energy (uncalibrated)')
        frame.Draw()
        canvas.Print(FileNames[ID[0][:2]])

    return ((mean.getVal(), mean.getError()),
            (sigma.getVal(), sigma.getError()),
            (erfc_coeff.getVal(), sigma.getError()))

def GetResolutionForTheta(ID, EnergyList2D, Theta):
    # Build rotated energy dataset for this theta
    RotatedEnergy = ROOT.RooRealVar('RotatedEnergy', 'RotatedEnergy', 0., MaxRotatedEnergy)
    ArgSet = ROOT.RooArgSet(RotatedEnergy)
    dataSet = ROOT.RooDataSet('dataSet', 'dataset', ArgSet)
    for Energy2D in EnergyList2D:
        RotatedEnergy.setVal(Energy2D[0]*math.cos(Theta) + Energy2D[1]*math.sin(Theta))
        dataSet.add(ArgSet)

    # Find the peak.
    InitPeakPos = FindPeak(dataSet)

    # Initialize guesses.  (Uncertainties are bogus.)
    PeakPos = (InitPeakPos, 0.1*InitPeakPos)
    Sigma = (0.03*InitPeakPos, 0.01*InitPeakPos)
    ErfcFrac = (0.1, 0.01)

    # Do a series of intermediate fits, to improve guesses.
    # This turns out to be necessary to get smooth results in some cases.
    for NumFit in range(2):
        PeakPos, Sigma, ErfcFrac = DoFit(dataSet, PeakPos[0], Sigma[0], ErfcFrac[0])

    # Final fit.  Save the results.
    PeakPos, Sigma, ErfcFrac = DoFit(dataSet,
                                     PeakPos[0],
                                     Sigma[0],
                                     ErfcFrac[0],
                                     ID = (ID, '%.5f' % Theta))

    Resolution = DivideWithErrors(Sigma, PeakPos)
    return Resolution, PeakPos

def LoadRun(data_ss, data_ms, inchain):
    # Load the DataSet, making use of EXOFitter to ensure consistency.

    # Create a fitter object with a workspace.
    fitter = ROOT.EXOFitter()
    fitter.BuildWorkspace()

    # Set basic options.
    fitter.fUseBetaScale = False
    fitter.SetApplyVetoes(False)
    fitter.SetRandomTrigger(True)
    fitter.fWhichScintVariable = WhichScint
    ROOT.EXOFiducialVolume.SetUserHexCut( 153., 15., 182 )

    # Build a dataset, and save it into a temporary file.
    print("Reading data into a condensed tree.")
    fitter.LoadData(inchain)
    fitter.BuildDataTree(prefix + '_temptreefile.root')

    # Create the RooDataSet from the tree in the temporary file.
    print("Reading condensed tree into a RooDataSet.")
    e_charge = ROOT.RooRealVar("e_charge", "e_charge", 0., 10000., "keV")
    e_scint = ROOT.RooRealVar("e_scint", "e_scint", 0., 10000., "keV")
    multiplicity = ROOT.RooRealVar("multiplicity", "multiplicity", 0, 10)
    observables = ROOT.RooArgSet("obs")
    observables.add(e_charge)
    observables.add(e_scint)
    observables.add(multiplicity)
    fitter.BuildDataSetFromTree(prefix + '_temptreefile.root', observables)
    alldata = fitter.GetWorkspace().data("EXOData_run2")
    if alldata == None: print "Failed to retrieve EXOData_run2 from the workspace."

    # Fill data_ss and data_ms from the RooDataSet.
    print("Reading RooDataSet into python lists.")
    for i in range(alldata.numEntries()):
        row = alldata.get(i)
        mult = row.find("multiplicity").getVal()
        if 0.5 < mult and mult < 1.5: # SS
            data_ss.append( (row.find("e_charge").getVal(), row.find("e_scint").getVal()) )
        elif 1.5 <= mult: # MS
            data_ms.append( (row.find("e_charge").getVal(), row.find("e_scint").getVal()) )

    # Delete the temporary file which was created.
    os.remove(prefix + '_temptreefile.root')

def Run(**kwargs):

    # Check that you've added all needed attributes.
    # Possibly this code design, where attributes may or may not be defined,
    # may be bad code design; but it is efficient, and
    # just reinforces the idea that this module is not meant to be used directly.
    if 'prefix' not in globals():
        sys.exit('You need to define the \'prefix\' attribute.')
    if 'WhichScint' not in globals():
        sys.exit('You need to define the \'WhichScint\' attribute.')

    # Only process thorium source runs.
    beginRecord = kwargs['ControlRecordList'].GetNextRecord('EXOBeginRecord')()
    if not isinstance(beginRecord, ROOT.EXOBeginSourceCalibrationRunRecord):
        print "This is not a source run; skipping."
        return
    if 'Th' not in beginRecord.GetSourceTypeString():
        print "We only handle thorium runs here; skipping."
        return

    # ROOT settings.
    ROOT.gROOT.SetStyle("Plain")
    ROOT.gStyle.SetPalette(1)
    ROOT.gStyle.SetOptFit(11)
    ROOT.TH1.AddDirectory(False)

    # Initialize canvases.
    global canvas
    global FileNames
    FileNames['SS'] = prefix + '_ss.pdf'
    FileNames['MS'] = prefix + '_ms.pdf'
    canvas = ROOT.TCanvas()
    for name in FileNames.itervalues(): canvas.Print(name + '[')

    # We'll need the drift velocity.
    ROOT.gSystem.Load("libEXOCalibUtilities")
    CalibManager = ROOT.EXOCalibManager.GetCalibManager()
    print "There are " + str(kwargs['EventTree'].GetEntries()) + " entries."

    # Create lists of charge and scint energy, using EXOFitter functions to do the heavy lifting.
    EnergyList2D_ss = []
    EnergyList2D_ms = []
    LoadRun(EnergyList2D_ss, EnergyList2D_ms, kwargs['EventTree'])

    ThetaToTry = [0.05 + 0.01*i for i in range(60)]
    graph_ss = ROOT.TGraphErrors()
    graph_ms = ROOT.TGraphErrors()

    BestTheta_ss = None
    BestRes_ss = None
    BestTheta_ms = None
    BestRes_ms = None

    try:
        ssScintResolution, ssScintMean = GetResolutionForTheta('SS Scint', EnergyList2D_ss, ROOT.TMath.Pi()/2)
        ssIonizResolution, ssIonizMean = GetResolutionForTheta('SS Ioniz', EnergyList2D_ss, 0.)
        msScintResolution, msScintMean = GetResolutionForTheta('MS Scint', EnergyList2D_ms, ROOT.TMath.Pi()/2)
        msIonizResolution, msIonizMean = GetResolutionForTheta('MS Ioniz', EnergyList2D_ms, 0.)

        for TestTheta in ThetaToTry:
            print "TestTheta ", TestTheta
            ssResolution, _ = GetResolutionForTheta('SS', EnergyList2D_ss, TestTheta)
            msResolution, _ = GetResolutionForTheta('MS', EnergyList2D_ms, TestTheta)

            graph_ss.Set(graph_ss.GetN()+1)
            graph_ss.SetPoint(graph_ss.GetN()-1, TestTheta, ssResolution[0])
            graph_ss.SetPointError(graph_ss.GetN()-1, 0., ssResolution[1])
            graph_ms.Set(graph_ms.GetN()+1)
            graph_ms.SetPoint(graph_ms.GetN()-1, TestTheta, msResolution[0])
            graph_ms.SetPointError(graph_ms.GetN()-1, 0., msResolution[1])

            if BestRes_ss == None or ssResolution[0] < BestRes_ss:
                BestRes_ss = ssResolution[0]
                BestTheta_ss = TestTheta
            if BestRes_ms == None or msResolution[0] < BestRes_ms:
                BestRes_ms = msResolution[0]
                BestTheta_ms = TestTheta

        # ToDO:  Consider trying to pick a fit window more intelligently.
        ssPolynomialFunc = ROOT.TF1("ssFit", "[0] + [1]*(x-[2])**2", BestTheta_ss - 0.07, BestTheta_ss + 0.07)
        ssPolynomialFunc.SetParameters(BestRes_ss, 4.0, BestTheta_ss)
        graph_ss.Fit(ssPolynomialFunc, "ER")
        msPolynomialFunc = ROOT.TF1("msFit", "[0] + [1]*(x-[2])**2", BestTheta_ms - 0.07, BestTheta_ms + 0.07)
        msPolynomialFunc.SetParameters(BestRes_ms, 4.0, BestTheta_ms)
        graph_ms.Fit(msPolynomialFunc, "ER")

        ROOT.gStyle.SetOptFit(11)
        graph_ss.SetTitle('SS Resolution vs Theta')
        graph_ss.GetXaxis().SetTitle('Theta (radians)')
        graph_ss.GetYaxis().SetTitle('Resolution (at thorium peak)')
        graph_ss.Draw("A*")
        canvas.Modified()
        canvas.Update()
        canvas.SaveAs(FileNames['SS'])
        graph_ms.SetTitle('MS Resolution vs Theta')
        graph_ms.GetXaxis().SetTitle('Theta (radians)')
        graph_ms.GetYaxis().SetTitle('Resolution (at thorium peak)')
        graph_ms.Draw("A*")
        canvas.Modified()
        canvas.Update()
        canvas.SaveAs(FileNames['MS'])

        ssResolution, ssMean = GetResolutionForTheta('SS', EnergyList2D_ss, ssPolynomialFunc.GetParameter(2))
        msResolution, msMean = GetResolutionForTheta('MS', EnergyList2D_ms, msPolynomialFunc.GetParameter(2))

        for name in FileNames.itervalues(): canvas.Print(name + ']')

        return { 'Theta_ss' : (ssPolynomialFunc.GetParameter(2), ssPolynomialFunc.GetParError(2)),
                 'PeakPos_ss' : ssMean,
                 'Resolution_ss' : ssResolution,
                 'Theta_ms' : (msPolynomialFunc.GetParameter(2), msPolynomialFunc.GetParError(2)),
                 'PeakPos_ms' : msMean,
                 'Resolution_ms' : msResolution,
                 'ssScintPeak' : ssScintMean,
                 'ssScintRes' : ssScintResolution,
                 'ssIonizPeak' : ssIonizMean,
                 'ssIonizRes' : ssIonizResolution,
                 'msScintPeak' : msScintMean,
                 'msScintRes' : msScintResolution,
                 'msIonizPeak' : msIonizMean,
                 'msIonizRes' : msIonizResolution }

    except PeakFinderError, err:
        # Clean up and terminate -- hopefully a person will look at this, and fix it.
        for name in FileNames.itervalues(): canvas.Print(name + ']')
        print 'Peak-finding failed: ', err
        print 'Terminating this script!!!  (Someone please figure out what went wrong.)'
        return

