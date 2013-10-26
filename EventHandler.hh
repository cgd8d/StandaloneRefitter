#ifndef EventHandler_hh
#define EventHandler_hh

#include "Rtypes.h"
#include "mkl_lapacke.h"
#include <vector>
#include <map>
#include <utility>

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
  std::map<unsigned char, double> fAPDGainMapEval;

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

#endif
