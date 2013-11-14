#ifndef EventHandler_hh
#define EventHandler_hh

#include "ModelManager.hh"
#include "Rtypes.h"
#include "mkl_lapacke.h"
#include <vector>
#include <map>
#include <utility>

#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>

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
  std::vector<unsigned char> fChannels;

#ifdef ENABLE_CHARGE
  std::vector<ModelManager> fWireModel; // U-wire model information.
#endif
  std::vector<ModelManager> fAPDModel; // APD model information.
  double fExpectedEnergy_keV; // For appropriate handling of Poisson noise.
  std::map<unsigned char, double> fAPDGainMapEval;

  std::vector<ModelManager*> fModels; // Pointers to APD and u-wire models.

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

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive & ar, unsigned int) {
    // Only serialize/unserialize the data structures we need in EventFinisher.
    ar & fEntryNumber;
    ar & fRunNumber;
    ar & fEventNumber;
    ar & fColumnLength;
    ar & fNumSignals;
    ar & fChannels;
#ifdef USE_CHARGE
    ar & fWireModel;
#endif
    ar & fAPDModel;
    ar & fX;
  }
};

// We want to sort event handlers by file position; makes file reading much more efficient.
struct CompareEventHandlerPtrs {
  bool operator()(const EventHandler* const &evt1, const EventHandler* const &evt2) {
    return (evt1->fRunNumber < evt2->fRunNumber) or
           (evt1->fRunNumber == evt2->fRunNumber and evt1->fEventNumber < evt2->fEventNumber);
  }
};

#endif
