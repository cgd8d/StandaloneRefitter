#ifndef EXORefitSignals_hh
#define EXORefitSignals_hh

#include "EXOUtilities/EXOTemplWaveform.hh"
#include "TStopwatch.h"
#include "mkl_cblas.h"
#include "mkl_lapacke.h"
#include <string>
#include <vector>
#include <map>
#include <list>
#include <cassert>

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
  bool fOnlyEstimateConditionNumbers; // Use for studying the quality of preconditioners.

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
    std::vector<double> fR0hat_V_factors;
    std::vector<lapack_int> fR0hat_V_pivot;
    std::vector<double> fprecon_tmp; // For storing the right-preconditioned version of a vector.

    // Preconditioner stuff.
    // We approximate approx(A) = {{D L} {trans(L) 0}}, where D is diagonal.
    // Then approx(A) = {{D^(0.5) 0} {trans(L)D^(-0.5) -trans(X)}}.{{D^(0.5) D^(-0.5)L}{0 X}},
    // for some X.
    // X can be obtained by doing Cholesky factoring with trans(X)X = trans(L) D^(-1) L,
    // which is done using LAPACK.
    // We precondition with K1 as the first, and K2 as the second;
    // both are easy to invert.
    std::vector<double> fDiag; // Diagonal entries of noise matrix.
    std::vector<double> fPreconX; // X, which is upper-triangular (unpacked).

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
  std::vector<double> fNoiseDiag;
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
  TStopwatch fWatch_FillNoise;
  TStopwatch fWatch_TotalTime;
  size_t fNumEventsHandled;
  size_t fNumSignalsHandled;
  size_t fTotalIterationsDone;

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
  void DoPoissonMultiplication(const std::vector<double>& in,
                               std::vector<double>& out,
                               EventHandler& event);
  template<char WHICH>
  void DoLagrangeAndConstraintMul(const std::vector<double>& in,
                                  std::vector<double>& out,
                                  EventHandler& event);

  // Preconditioner functions.
  std::vector<double> DoInvLPrecon(std::vector<double>& in, EventHandler& event);
  std::vector<double> DoInvRPrecon(std::vector<double>& in, EventHandler& event);
  std::vector<double> DoLPrecon(std::vector<double>& in, EventHandler& event);
  std::vector<double> DoRPrecon(std::vector<double>& in, EventHandler& event);
  std::pair<double, double> EstimateConditionNumber(EventHandler& event);

  // Produce the light model, used on all gangs.
  EXOWaveformFT GetModelForTime(double time) const;
};

template<char WHICH>
void EXORefitSignals::DoLagrangeAndConstraintMul(const std::vector<double>& in,
                                                 std::vector<double>& out,
                                                 EventHandler& event)
{
  // Do multiplication by Lagrange and constraint terms.
  // If WHICH = 'L', only handle Lagrange terms.
  // If WHICH = 'C', only handle constraint terms.
  // If WHICH = 'A', do both.
  // All other values of WHICH result in an error.  (This is done to force compile-time optimization.)
  // In and out should only be equal if WHICH != 'A'.
  // We also package the "preconditioning" matrix, so this actually does:
  // (0                D^(-1/2)L)
  // (trans(L)D^(-1/2) 0        )
  assert(WHICH == 'L' or WHICH == 'C' or WHICH == 'A');
  assert(WHICH == 'L' or WHICH == 'C' or &in[0] != &out[0]);
  bool Lagrange = (WHICH == 'L' or WHICH == 'A');
  bool Constraint = (WHICH == 'C' or WHICH == 'A');

  // First loop through wire signals.
  for(size_t m = 0; m < event.fWireModel.size(); m++) {
    const std::map<unsigned char, std::vector<double> >& models = event.fWireModel[m].second;
    for(std::map<unsigned char, std::vector<double> >::const_iterator it = models.begin();
        it != models.end();
        it++) {
      unsigned char ChannelWithWireSignal = it->first;
      size_t channel_index = 0;
      while(fChannels[channel_index] != ChannelWithWireSignal) {
        channel_index++;
        if(channel_index >= fChannels.size()) LogEXOMsg("Index exceeded -- why can this happen?", EEAlert);
      }
      const std::vector<double>& modelWF = it->second;
      for(size_t f = 0; f <= fMaxF - fMinF; f++) {
        size_t Index1 = event.fColumnLength - (event.fWireModel.size()+1) + m;
        size_t Index2 = 2*fChannels.size()*f + channel_index*(f < fMaxF-fMinF ? 2 : 1);
        for(size_t n = 0; n <= event.fWireModel.size(); n++) {
          if(Lagrange) out[Index2] += modelWF[2*f]*fNoiseDiag[Index2 % event.fColumnLength]*in[Index1];
          if(Constraint) out[Index1] += modelWF[2*f]*fNoiseDiag[Index2 % event.fColumnLength]*in[Index2];
          if(f < fMaxF-fMinF) {
            if(Lagrange) out[Index2+1] += modelWF[2*f+1]*fNoiseDiag[(Index2+1)%event.fColumnLength]*in[Index1];
            if(Constraint) out[Index1] += modelWF[2*f+1]*fNoiseDiag[(Index2+1)%event.fColumnLength]*in[Index2+1];
          }
          Index1 += event.fColumnLength;
          Index2 += event.fColumnLength;
        }
      }
    }
  } // Done with Lagrange and constraint terms for wires.
  // Now, Lagrange and constraint terms for APDs.
  for(size_t k = fFirstAPDChannelIndex; k < fChannels.size(); k++) {
    double ExpectedYieldOnGang = event.fExpectedYieldPerGang.at(fChannels[k]);
    for(size_t f = 0; f <= fMaxF - fMinF; f++) {
      size_t Index1 = 2*fChannels.size()*f + k*(f < fMaxF-fMinF ? 2 : 1);
      size_t Index2 = event.fColumnLength - 1;
      for(size_t n = 0; n <= event.fWireModel.size(); n++) {
        if(Constraint) out[Index2] += event.fmodel_realimag[2*f]*ExpectedYieldOnGang*fNoiseDiag[Index1%event.fColumnLength]*in[Index1];
        if(Lagrange) out[Index1] += event.fmodel_realimag[2*f]*ExpectedYieldOnGang*fNoiseDiag[Index1%event.fColumnLength]*in[Index2];
        if(f < fMaxF-fMinF) {
          if(Constraint) out[Index2] += event.fmodel_realimag[2*f+1]*ExpectedYieldOnGang*fNoiseDiag[(Index1+1)%event.fColumnLength]*in[Index1+1];
          if(Lagrange) out[Index1+1] += event.fmodel_realimag[2*f+1]*ExpectedYieldOnGang*fNoiseDiag[(Index1+1)%event.fColumnLength]*in[Index2];
        }
        Index1 += event.fColumnLength;
        Index2 += event.fColumnLength;
      }
    }
  }
}
#endif
