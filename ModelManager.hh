#ifndef ModelManager_hh
#define ModelManager_hh
/*
For an individual signal, this object holds a mapping from channel number to expected signal.
It is general enough to be used for charge or light signals.
*/

#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <vector>
#include <map>
#include <set>
#include <cstddef>

struct ModelManager
{
  ModelManager() : fSignalNumber(-1) {} // For the benefit of deserialization.

  ModelManager(size_t ModelSize, size_t NumChannels)
  : fSignalNumber(-1),
    fNumChannels(NumChannels),
    fModel(ModelSize, 0)
  {
    assert(ModelSize % NumChannels == 0);
  }

  // Insert a Channel hit, replacing previously-held zero values.
  // Do not attempt to insert models for the same channel more than once.
  void AddChannelHit(unsigned char channel, const std::vector<double>& model) {
    assert(channel < fNumChannels);
    assert(model.size()*fNumChannels == fModel.size());
    assert(fHitChannels.count(channel) == 0);
    fHitChannels.insert(channel);
    for(size_t i = 0; i < model.size(); i++) fModel[fNumChannels*i + channel] = model[i];
  }

  // Pre-multiply the models by fInvSqrtNoiseDiag, for efficiency.
  // Also locate contiguous runs of hit channels.
  // This function should be called when you are done adding hit channels.
  void Finalize(const std::vector<double>& InvSqrtNoiseDiag) {
    assert(fHitChannels.size() != 0);
    assert(InvSqrtNoiseDiag.size() == fModel.size());
    for(size_t i = 0; i < fModel.size(); i++) fModel[i] *= InvSqrtNoiseDiag[i];

    std::set<unsigned char>::iterator it = fHitChannels.begin();
    while(it != fHitChannels.end()) {
      unsigned char startRange = *it;
      unsigned char endRange = startRange;
      do {
        it++;
        endRange++;
      } while(it != fHitChannels.end() and *it == endRange);
      fContiguousChannels.insert(std::make_pair(startRange, endRange));
    }
  }

  // The index number of this signal in EXOEventData, eg GetScintillationCluster(fSignalNumber).
  // Needed only so that when it is time to write the denoised magnitude,
  // we know where to put it.
  size_t fSignalNumber;

  // The number of channels -- this is known elsewhere, but permits utility functions defined here.
  // This is equivalent to the stride, if scanning through 
  size_t fNumChannels;

  // Model, with indexing to match other vectors.
  std::vector<double> fModel;

  // Which channels have non-zero models (by channel index, not software channel).
  std::set<unsigned char> fHitChannels;

  // For cases where it is convenient to use contiguous channels as a block, store ranges.
  // These ranges, like any C++ iterator, are inclusive on the lower end and exclusive on the higher end.
  std::set<std::pair<unsigned char, unsigned char> > fContiguousChannels;

  // For APDs only -- expected yield (ADC counts for a 2615 keV deposit) on each gang.
  std::map<unsigned char, double> fExpectedYieldPerGang;

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive & ar, unsigned int) {
    // Only serialize/unserialize the data structures we need in EventFinisher.
    ar & fSignalNumber;
  }
};
#endif
