#ifndef EXORefitSignals_hh
#define EXORefitSignals_hh

#include "EXOAnalysisManager/EXOAnalysisModule.hh"
#include "EXOUtilities/EXOTemplWaveform.hh"
#include "TStopwatch.h"
#include <string>
#include <vector>
#include <map>

class EXOUWireSignal;
class EXOTalkToManager;
class EXOTransferFunction;
class EXOEventData;
class EXOWaveformFT;
class TH3D;
class TGraph;

class EXORefitSignals : public EXOAnalysisModule
{
 public:
  EXORefitSignals() : fLightmapFilename("data/lightmap/LightMaps.root"),
                      fRThreshold(0.1),
                      fThoriumEnergy_keV(2615),
                      fMinF(1),
                      fMaxF(1024),
                      fNumEntriesSolved(0),
                      fTotalNumberOfIterationsDone(0),
                      fTotalIterationsForWires(0),
                      fTotalIterationsForAPDs(0) {}
  int TalkTo(EXOTalkToManager *tm);
  int Initialize();
  EXOAnalysisModule::EventStatus ProcessEvent(EXOEventData *ED);
  int ShutDown();

  void SetNoiseFilename(std::string name) { fNoiseFilename = name; }
  void SetLightmapFilename(std::string name) { fLightmapFilename = name; }
  void SetRThreshold(double threshold) { fRThreshold = threshold; }

 protected:

  // fNoiseCorrelations[f-fMinF] stores the matrix of noise correlations at frequency f.
  // This matrix is a single contiguous array, ordered like:
  // <N^R_0 N^R_0> <N^I_0 N^R_0> <N^R_1 N^R_0> ... <N^R_0 N^I_0> ...
  // (Ie. in column-major format -- this facilitates the use of GEMM if a BLAS library is available.)
  // The number corresponds to the *index* of the various channels, of course.
  std::string fNoiseFilename;
  std::vector<std::vector<double> > fNoiseCorrelations;
  void FillNoiseCorrelations(const EXOEventData& ED);

  std::string fLightmapFilename;
  std::map<unsigned char, TH3D*> fLightMaps;
  std::map<unsigned char, TGraph*> fGainMaps;
  std::vector<unsigned char> fChannels;
  size_t fFirstAPDChannelIndex;

  double GetGain(unsigned char channel) const;

  double fRThreshold;

  // Various stopwatches, to understand the fraction of time spent actually solving the matrix.
  TStopwatch fWatch_GetNoise;
  TStopwatch fWatch_ProcessEvent;
  TStopwatch fWatch_InitialGuess;
  TStopwatch fWatch_Solve;
  TStopwatch fWatch_NoiseMul;
  TStopwatch fWatch_RestMul;

  const double fThoriumEnergy_keV;
  const size_t fMinF;
  const size_t fMaxF;

  // Wire digitization.
  EXODoubleWaveform fWireDeposit;
  EXODoubleWaveform fWireInduction;
  std::vector<double> MakeWireModel(EXODoubleWaveform& in,
                                    const EXOTransferFunction& transfer,
                                    const double Gain,
                                    const double Time) const;

  // Collect statistics on the number of iterations required.
  unsigned long int fNumEntriesSolved;
  unsigned long int fTotalNumberOfIterationsDone;
  unsigned long int fTotalIterationsForWires;
  unsigned long int fTotalIterationsForAPDs;

  // Block-BiCGSTAB solver stuff.
  void BiCGSTAB_iteration(std::vector<double>& X,
                          std::vector<double>& R,
                          std::vector<double>& P,
                          const std::vector<double>& r0hat) const;
  bool DoBiCGSTAB(std::vector<double>& X,
                  double Threshold);


  struct EventHandler {
    // So we can grab the event again when we're done.
    Int_t fRunNumber;
    Int_t fEventNumber;
    double fUnixTimeOfEvent;
    size_t fColumnLength;

    // fWireModel keeps, for each u-wire signal we're fitting:
    //   the index of the u-wire signal within the event.
    //   model waveforms for channels which it affects.
    // The map keys are software channels.
    // The model waveforms are normalized so that the shaped deposition has a peak-baseline of 1 ADC.
    // The ordering in the matrix is defined by the vector index, of course.
    std::vector<std::pair<size_t,
                          std::map<unsigned char, std::vector<double> > > > fWireModel;

    // APD model information.
    std::map<unsigned char, double> fExpectedYieldPerGang; // Expected magnitudes of Th gamma line (ADC).
    std::vector<double> fmodel_realimag;
    double fExpectedEnergy_keV; // For appropriate handling of Poisson noise.

    // Information on the current status and data of the solver.
    // We need enough information so that when a matrix multiplication with noise finishes,
    // we can pick up the pieces.
    // We can re-enter in the setup phase, just after computing V, or just after computing T.
    // So, identify the phase based on which vectors have a size of zero.
    std::vector<double> fX;
    std::vector<double> fR;
    std::vector<double> fP;
    std::vector<double> fR0hat;
    std::vector<double> fV; // only needed within an iteration.
    std::vector<double> fAlpha; // only needed within an iteration.
    std::vector<double> fR0hat_V_Inv; // only needed within an iteration; not strictly needed, but convenient.

    // Where in the result matrix can we expect to find the required result?
    size_t fResultIndex;
  };
  void DoBlBiCGSTAB_iteration(EventHandler& event);


  // Functions to multiply by the noise matrix.
  std::vector<double> fNoiseMulQueue;
  std::vector<double> fNoiseMulResult;
  size_t fNumVectorsInQueue;
  void DoNoiseMultiplication();

  // Function to perform the rest of matrix multiplication.
  void DoRestOfMultiplication(const std::vector<double>& in,
                              std::vector<double>& out,
                              EventHandler& event);

  EXOWaveformFT GetModelForTime(double time) const;

 DEFINE_EXO_ANALYSIS_MODULE(EXORefitSignals)
};
#endif
