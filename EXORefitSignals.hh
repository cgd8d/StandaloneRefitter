#ifndef EXORefitSignals_hh
#define EXORefitSignals_hh

#include "EXOUtilities/EXOTemplWaveform.hh"
#include "TStopwatch.h"
#include <string>
#include <vector>
#include <map>
#include <list>

class EXOUWireSignal;
class EXOTransferFunction;
class EXOEventData;
class EXOWaveformFT;
class EXOTreeInputModule;
class EXOTreeOutputModule;
class TH3D;
class TGraph;
class TTree;

class EXORefitSignals
{
 public:
  // Functions specified in the order they should be called.
  EXORefitSignals(EXOTreeInputModule& inputModule,
                  TTree& wfTree,
                  EXOTreeOutputModule& outputModule);

  void SetNoiseFilename(std::string name) { fNoiseFilename = name; }
  void SetLightmapFilename(std::string name) { fLightmapFilename = name; }
  void SetRThreshold(double threshold) { fRThreshold = threshold; }
  bool fUseWireAPDCorrelations; // For now, this isn't higher performance -- just for testing.

  int Initialize();
  void AcceptEvent(EXOEventData* ED, Long64_t entryNum);
  void FlushEvents();

  ~EXORefitSignals();

 protected:

  // Handles so that we can read events in and save events as necessary.
  EXOTreeInputModule& fInputModule;
  TTree& fWFTree;
  EXOEventData* fWFEvent;
  EXOTreeOutputModule& fOutputModule;

  struct EventHandler {
    // So we can grab the event again when we're done.
    Long64_t fEntryNumber;
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
    std::vector<double> fprecon_tmp; // For storing the right-preconditioned version of a vector.

    // Preconditioner stuff.
    std::vector<double> fLPrecon; // K1_inv
    std::vector<double> fRPrecon; // K2_inv
    void DoInvLPrecon(std::vector<double>& in);
    void DoInvRPrecon(std::vector<double>& in);
    void DoLPrecon(std::vector<double>& in);
    void DoRPrecon(std::vector<double>& in);

    // Where in the result matrix can we expect to find the required result?
    size_t fResultIndex;
  };

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

  double GetGain(unsigned char channel, EventHandler& event) const;

  double fRThreshold;

  // Various stopwatches, to understand the fraction of time spent actually solving the matrix.
  TStopwatch fWatch_BiCGSTAB;
  TStopwatch fWatch_BiCGSTAB_part1;
  TStopwatch fWatch_BiCGSTAB_part2;
  TStopwatch fWatch_NoiseMul;
  TStopwatch fWatch_RestMul;
  size_t fNumEventsHandled;
  size_t fNumSignalsHandled;
  size_t fTotalIterationsDone;
  double fInitialNormWires;
  double fInitialNormAPDs;

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

  // Interact with files.
  std::list<EventHandler*> fEventHandlerList;
  void FinishEvent(EventHandler* event);

  // Block BiCGSTAB algorithm.
  bool DoBlBiCGSTAB(EventHandler& event);

  // Functions to multiply by the noise matrix.
  size_t fNoiseColumnLength;
  std::vector<double> fNoiseMulQueue;
  std::vector<double> fNoiseMulResult;
  size_t fNumVectorsInQueue;
  void DoNoiseMultiplication();
  size_t RequestNoiseMul(std::vector<double>& vec,
                         size_t ColLength);
  void FillFromNoise(std::vector<double>& vec,
                     size_t NumCols,
                     size_t ColLength,
                     size_t ResultIndex);

  // Function to perform the rest of matrix multiplication.
  void DoRestOfMultiplication(const std::vector<double>& in,
                              std::vector<double>& out,
                              EventHandler& event);

  // Produce the light model, used on all gangs.
  EXOWaveformFT GetModelForTime(double time) const;
};
#endif
