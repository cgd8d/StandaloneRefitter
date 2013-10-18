#ifndef EXORefitSignals_hh
#define EXORefitSignals_hh

#include "SafeStopwatch.hh"
#include "EXOUtilities/EXOWaveformData.hh"
#include "Rtypes.h"
#include "mkl_cblas.h"
#include "mkl_lapacke.h"
#include "mkl_vml_functions.h"
#include <string>
#include <vector>
#include <set>
#include <map>
#include <cassert>

#ifdef USE_THREADS
#include <boost/atomic.hpp>
#include <boost/thread/mutex.hpp>

// Some mutexes need to be accessible in multiple translation units -- declare them here.
extern boost::mutex ProcessedFileMutex;
#endif

class EXOEventData;
class EXOWaveformFT;
class EXOTreeInputModule;
class EXOTreeOutputModule;
class TH3D;
class TGraph;
class TTree;

#ifdef ENABLE_CHARGE
#include "EXOUtilities/EXOTemplWaveform.hh"
class EXOUWireSignal;
class EXOTransferFunction;
#endif

struct EventHandler {
  // So we can grab the event again when we're done.
  Long64_t fEntryNumber;
  Int_t fRunNumber;
  Int_t fEventNumber;
  double fUnixTimeOfEvent;
  size_t fColumnLength;
  size_t fNumSignals;
  size_t fNumIterSinceReset;
  size_t fNumIterations; // Count the number of times we've tried to terminate.

#ifdef ENABLE_CHARGE
  // fWireModel keeps, for each u-wire signal we're fitting:
  //   the index of the u-wire signal within the event.
  //   model waveforms for channels which it affects.
  // The map keys are software channels.
  // The model waveforms are normalized so that the shaped deposition has a peak-baseline of 1 ADC.
  // The ordering in the matrix is defined by the vector index, of course.
  std::vector<std::pair<size_t,
                        std::map<unsigned char, std::vector<double> > > > fWireModel;
#endif
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
  std::vector<double> fPreconX; // X, which is upper-triangular (unpacked).

  // Where in the result matrix can we expect to find the required result?
  size_t fResultIndex;
};

// We want to sort event handlers by file position; makes file reading much more efficient.
struct CompareEventHandlerPtrs {
  bool operator()(const EventHandler* const &evt1, const EventHandler* const &evt2) {
    return (evt1->fRunNumber < evt2->fRunNumber) or
           (evt1->fRunNumber == evt2->fRunNumber and evt1->fEventNumber < evt2->fEventNumber);
  }
};

#ifdef USE_LOCKFREE
// Use a lock-free queue so that multiple threads can push and pop events to be handled without a manager.
#include <boost/lockfree/queue.hpp>
typedef boost::lockfree::queue<EventHandler*> queue_type;
#else
// If we're threading, we'll have to use regular locks around a regular queue.
// If we're not threading, of course just use a queue anyway.
#include <queue>
typedef std::queue<EventHandler*> queue_type;
#endif

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
#ifdef ENABLE_CHARGE
  bool fAPDsOnly; // Do not denoise wire signals; and do not use u-wires to denoise APDs.
  bool fUseWireAPDCorrelations; // For now, this isn't higher performance -- just for testing.
#endif
  bool fVerbose;
  size_t fDoRestarts; // 0 if we never restart; else, value indicates number of iterations before a restart.
  size_t fNumMulsToAccumulate;

  int Initialize();
  void AcceptEvent(EXOEventData* ED, Long64_t entryNum);
  void FlushEvents();

#ifdef USE_THREADS
  boost::atomic<bool> fProcessingIsDone; // Permit notification that no more events will be finished.
#endif
  void FinishEventThread();
  size_t GetFinishEventQueueLength();

  ~EXORefitSignals();

 protected:

  // Handles so that we can read events in and save events as necessary.
  EXOTreeInputModule& fInputModule;
  TTree& fWFTree;
  EXOWaveformData fWFData;
  EXOTreeOutputModule& fOutputModule;

  // fNoiseCorrelations[f-fMinF] stores the matrix of noise correlations at frequency f.
  // This matrix is a single contiguous array, ordered like:
  // <N^R_0 N^R_0> <N^I_0 N^R_0> <N^R_1 N^R_0> ... <N^R_0 N^I_0> ...
  // (Ie. in column-major format -- this facilitates the use of GEMM if a BLAS library is available.)
  // The number corresponds to the *index* of the various channels, of course.
  std::string fNoiseFilename;
  std::vector<std::vector<double> > fNoiseCorrelations;
  std::vector<double> fNoiseDiag;
  std::vector<double> fInvSqrtNoiseDiag;
  void FillNoiseCorrelations(const EXOEventData& ED);

  std::string fLightmapFilename;
  std::map<unsigned char, TH3D*> fLightMaps;
  std::map<unsigned char, TGraph*> fGainMaps;
  std::vector<unsigned char> fChannels;
  size_t fFirstAPDChannelIndex;

  double GetGain(unsigned char channel, EventHandler& event) const;

  double fRThreshold;
  bool CanTerminate(EventHandler& event);

  // Various counters.
  size_t fNumEventsHandled;
  size_t fNumSignalsHandled;
  size_t fTotalIterationsDone;

  const double fThoriumEnergy_keV;
  const size_t fMinF;
  const size_t fMaxF;

#ifdef ENABLE_CHARGE
  // Wire digitization.
  EXODoubleWaveform fWireDeposit;
  EXODoubleWaveform fWireInduction;
  std::vector<double> MakeWireModel(EXODoubleWaveform& in,
                                    const EXOTransferFunction& transfer,
                                    const double Gain,
                                    const double Time) const;
#endif

  // Interact with files.
  queue_type fEventHandlerQueue;
  queue_type fEventHandlerResults;
  void HandleEventsInThread();
  void DoPassThroughEvents();
  EventHandler* PopAnEvent();
  void PushAnEvent(EventHandler* evt);

  std::set<EventHandler*, CompareEventHandlerPtrs> fEventsToFinish;
  void PushFinishedEvent(EventHandler* event);
  void FinishEvent(EventHandler* event);

  // Block BiCGSTAB algorithm.
  bool DoBlBiCGSTAB(EventHandler& event);
  void DoRestart(EventHandler& event);

  // Functions to multiply by the noise matrix.
  size_t fNoiseColumnLength;
  std::vector<double> fNoiseMulQueue;
  std::vector<double> fNoiseMulResult;
  size_t fNumVectorsInQueue;
  void DoNoiseMultiplication();
  void DoNoiseMultiplication_Range(size_t flo, size_t fhi);
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
  template<char WHICH, bool Add>
  void DoLagrangeAndConstraintMul(const std::vector<double>& in,
                                  std::vector<double>& out,
                                  EventHandler& event);

  // Preconditioner functions.
  void DoInvLPrecon(std::vector<double>& in, EventHandler& event);
  void DoInvRPrecon(std::vector<double>& in, EventHandler& event);
  void DoLPrecon(std::vector<double>& in, EventHandler& event);
  void DoRPrecon(std::vector<double>& in, EventHandler& event);

  // Produce the light model, used on all gangs.
  EXOWaveformFT GetModelForTime(double time) const;
};

template<char WHICH, bool Add>
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
  // If Add is true (default), we add this to out; otherwise, we subtract it.
  assert(WHICH == 'L' or WHICH == 'C' or WHICH == 'A');
  assert(WHICH == 'L' or WHICH == 'C' or &in[0] != &out[0]);
  bool Lagrange = (WHICH == 'L' or WHICH == 'A');
  bool Constraint = (WHICH == 'C' or WHICH == 'A');

#ifdef ENABLE_CHARGE
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
        if(channel_index >= fChannels.size()) {
          std::cout<<"On entry "<<event.fEntryNumber<<", in DoLagrangeAndConstraintMul, "<<
                     "Index exceeded -- why can this happen?"<<std::endl;
          std::exit(1);
        }
      }
      const std::vector<double>& modelWF = it->second;
      for(size_t f = 0; f <= fMaxF - fMinF; f++) {
        size_t Index1 = event.fColumnLength - event.fNumSignals + m;
        size_t Index2 = 2*fChannels.size()*f + channel_index;
        const size_t DiagIndex = Index2;
        for(size_t n = 0; n < event.fNumSignals; n++) {
          if(Lagrange) out[Index2] += (Add ? 1 : -1)*modelWF[2*f]*fInvSqrtNoiseDiag[DiagIndex]*in[Index1];
          if(Constraint) out[Index1] += (Add ? 1 : -1)*modelWF[2*f]*fInvSqrtNoiseDiag[DiagIndex]*in[Index2];
          if(f < fMaxF-fMinF) {
            if(Lagrange) out[Index2+fChannels.size()] += (Add ? 1 : -1)*modelWF[2*f+1]*fInvSqrtNoiseDiag[DiagIndex+fChannels.size()]*in[Index1];
            if(Constraint) out[Index1] += (Add ? 1 : -1)*modelWF[2*f+1]*fInvSqrtNoiseDiag[DiagIndex+fChannels.size()]*in[Index2+fChannels.size()];
          }
          Index1 += event.fColumnLength;
          Index2 += event.fColumnLength;
        }
      }
    }
  } // Done with Lagrange and constraint terms for wires.
#endif
  // Now, Lagrange and constraint terms for APDs.
  // This is where most of the time is spent, because we need to loop through all APD channels.
  // Furthermore, the APD channels are grouped together, making this portion of the code more matrix-friendly.
  // So, we exploit MKL as much as possible.
  std::vector<double> ExpectedYields;
  for(size_t k = fFirstAPDChannelIndex; k < fChannels.size(); k++) {
    ExpectedYields.push_back(event.fExpectedYieldPerGang.at(fChannels[k]));
  }
  std::vector<double> Workspace(ExpectedYields.size(), 0);
  for(size_t f = 0; f <= fMaxF - fMinF; f++) {
    size_t StartIndex = 2*fChannels.size()*f + fFirstAPDChannelIndex;

    // Start with real blocks.
    vdMul(ExpectedYields.size(), &ExpectedYields[0], &fInvSqrtNoiseDiag[StartIndex], &Workspace[0]);
    if(Constraint) {
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                  1, event.fNumSignals, ExpectedYields.size(),
                  (Add ? 1 : -1)*event.fmodel_realimag[2*f], &Workspace[0], 1,
                  &in[StartIndex], event.fColumnLength,
                  1, &out[event.fColumnLength-1], event.fColumnLength);
    }
    if(Lagrange) {
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                  ExpectedYields.size(), event.fNumSignals, 1,
                  (Add ? 1 : -1)*event.fmodel_realimag[2*f], &Workspace[0], ExpectedYields.size(),
                  &in[event.fColumnLength-1], event.fColumnLength,
                  1, &out[StartIndex], event.fColumnLength);
    }

    // Now do imaginary blocks.
    if(f == fMaxF - fMinF) continue;
    StartIndex += fChannels.size();
    vdMul(ExpectedYields.size(), &ExpectedYields[0], &fInvSqrtNoiseDiag[StartIndex], &Workspace[0]);
    if(Constraint) {
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                  1, event.fNumSignals, ExpectedYields.size(),
                  (Add ? 1 : -1)*event.fmodel_realimag[2*f+1], &Workspace[0], 1,
                  &in[StartIndex], event.fColumnLength,
                  1, &out[event.fColumnLength-1], event.fColumnLength);
    }
    if(Lagrange) {
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                  ExpectedYields.size(), event.fNumSignals, 1,
                  (Add ? 1 : -1)*event.fmodel_realimag[2*f+1], &Workspace[0], ExpectedYields.size(),
                  &in[event.fColumnLength-1], event.fColumnLength,
                  1, &out[StartIndex], event.fColumnLength);
    }
  }
}
#endif
