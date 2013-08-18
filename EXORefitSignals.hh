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
  size_t fColumnLength;

  double GetGain(unsigned char channel) const;

  double fRThreshold;

  // Various stopwatches, to understand the fraction of time spent actually solving the matrix.
  TStopwatch fWatch_GetNoise;
  TStopwatch fWatch_ProcessEvent;
  TStopwatch fWatch_InitialGuess;
  TStopwatch fWatch_Solve;
  mutable TStopwatch fWatch_MatrixMul;
  mutable TStopwatch fWatch_MatrixMul_NoiseTerms;

  std::map<unsigned char, double> fExpectedYieldPerGang; // Expected magnitudes of Th gamma line (ADC).
  std::vector<double> fmodel_realimag;
  double fUnixTimeOfEvent;
  double fExpectedEnergy_keV;
  const double fThoriumEnergy_keV;
  const size_t fMinF;
  const size_t fMaxF;

  // Wire digitization.
  // fWireModel keeps, for each u-wire signal we're fitting, a pointer back to that signal plus
  // model waveforms for channels which it affects.
  // The map keys are software channels.
  // The model waveforms are normalized so that the shaped deposition has a peak-baseline of 1 ADC.
  // The ordering in the matrix is defined by the vector index, of course.
  EXODoubleWaveform fWireDeposit;
  EXODoubleWaveform fWireInduction;
  std::vector<std::pair<EXOUWireSignal*,
                        std::map<unsigned char, std::vector<double> > > > fWireModel;
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

  std::vector<double> MatrixTimesVector(const std::vector<double>& in) const;

  EXOWaveformFT GetModelForTime(double time) const;

 DEFINE_EXO_ANALYSIS_MODULE(EXORefitSignals)
};
#endif
