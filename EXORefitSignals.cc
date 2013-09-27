//___________________________________________________
//
// Fit for scintillation and u-wire magnitudes in a noise-tolerant fashion.
// For more information, please see the slides and derivation from the Energy Meeting
// on July 22, 2013:
// https://confluence.slac.stanford.edu/display/exo/Energy+Meeting+22+July+2013
// (A full note will be written up in the near future.)
//
// Limitations:
// * Although this algorithm can be extended to handle events with more than one scintillation cluster,
//   for simplicity it does not currently do so.
// * The APD gains may not be right, which results in a suboptimal energy resolution.
// * In principle we could extract cluster-by-cluster light yield down the road.
// * It should be fairly cheap to also extract the estimated fit error due to electronic and Poisson noise.
//

#include "EXORefitSignals.hh"
#include "EXOAnalysisManager/EXOTreeInputModule.hh"
#include "EXOAnalysisManager/EXOTreeOutputModule.hh"
#include "EXOCalibUtilities/EXOChannelMapManager.hh"
#include "EXOCalibUtilities/EXOElectronicsShapers.hh"
#include "EXOCalibUtilities/EXOUWireGains.hh"
#include "EXOUtilities/EXOEventData.hh"
#include "EXOUtilities/EXODigitizeWires.hh"
#include "EXOUtilities/EXOWaveform.hh"
#include "EXOUtilities/EXOWaveformFT.hh"
#include "EXOUtilities/EXOFastFourierTransformFFTW.hh"
#include "EXOUtilities/EXOTransferFunction.hh"
#include "EXOUtilities/EXOMiscUtil.hh"
#include "EXOUtilities/EXOUWireSignal.hh"
#include "TFile.h"
#include "TTree.h"
#include "TArrayI.h"
#include "TH3D.h"
#include "TGraph.h"
#include "MKL_adaptors.hh"
#include "mkl_trans.h"
#include <iomanip>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include <fstream>

#ifdef USE_THREADS
// We currently use the boost::threads library, if threading is enabled.
#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread/detail/thread_group.hpp>

// Create basic mutexes for two multiple-writer situations.
boost::mutex FinishEventMutex; // tinput and toutput are not thread-safe.
boost::mutex RequestNoiseMulMutex; // Currently events are appended to the request queue.

#ifndef USE_LOCKFREE
// We'll need additional mutexes to protect our event queues.
boost::mutex EventQueueMutex;
boost::mutex EventResultsMutex;
#endif

#endif

EXORefitSignals::EXORefitSignals(EXOTreeInputModule& inputModule,
                                 TTree& wfTree,
                                 EXOTreeOutputModule& outputModule)
: fAPDsOnly(false),
  fUseWireAPDCorrelations(true),
  fInputModule(inputModule),
  fWFTree(wfTree),
  fOutputModule(outputModule),
  fLightmapFilename("data/lightmap/LightMaps.root"),
  fRThreshold(0.1),
  fThoriumEnergy_keV(2615),
  fMinF(1),
  fMaxF(1024),
#ifdef USE_LOCKFREE
  fEventHandlerQueue(0), // The default lockfree constructor is not allowed.
  fEventHandlerResults(0), // http://boost.2283326.n4.nabble.com/lockfree-Faulty-static-assert-td4635029.html
#endif
  fNumVectorsInQueue(0)
{
  fWFEvent = NULL;
  fWFTree.SetBranchAddress("EventBranch", &fWFEvent);
  fWatches["Time from construction to destruction"].Start();
}

void EXORefitSignals::FillNoiseCorrelations(const EXOEventData& ED)
{
  // Reorganize the noise matrix entries for fast use in matrix-vector multiplication.
  // This depends on the set of available waveforms; but we assume that this set doesn't
  // change much.  Compare to the static set in this function.
  // Note that we don't really expect it to change within a run; but the cost of protecting
  // against that is small.
  // For the future, might consider whether we can rearrange EXONoiseCorrelations to be
  // more convenient for this.
  fWatches["FillNoiseCorrelations"].Start(false);

  // Get the channel map.
  const EXOChannelMap& ChannelMap = GetChanMapForHeader(ED.fEventHeader);

  // Construct the set of channels to keep.
  std::vector<unsigned char> ChannelsToUse;
  for(unsigned char i = 0; i < NUMBER_READOUT_CHANNELS; i++) {
    if(EXOMiscUtil::TypeOfChannel(i) == EXOMiscUtil::kVWire) continue; // No v wires for now.
    if(fAPDsOnly and EXOMiscUtil::TypeOfChannel(i) != EXOMiscUtil::kAPDGang) continue;
    if(ChannelMap.channel_suppressed_by_daq(i) or not ChannelMap.good_channel(i)) continue;
    ChannelsToUse.push_back(i);
  }

  // If the channel mapping is unchanged, do nothing.
  if(ChannelsToUse == fChannels) {
    fWatches["FillNoiseCorrelations"].Stop();
    return;
  }

  // Else, we'll need to extract the noise information to match the new ordering.
  // Start by flushing all currently-held events, since the noise information will change.
  FlushEvents();

  // Build a map from channel index to where it can be found in the noise file.
  std::map<size_t, size_t> ChannelIndexMap;
  for(size_t i = 0; i < ChannelsToUse.size(); i++) {
    unsigned char software_channel = ChannelsToUse[i];
    if(software_channel < NCHANNEL_PER_WIREPLANE) ChannelIndexMap[i] = software_channel;
    else if(software_channel < 3*NCHANNEL_PER_WIREPLANE) ChannelIndexMap[i] = software_channel - NCHANNEL_PER_WIREPLANE;
    else ChannelIndexMap[i] = software_channel - 2*NCHANNEL_PER_WIREPLANE;
  }

  // For convenience, pre-store useful parameters.
  fChannels = ChannelsToUse;
  for(size_t i = 0; i < fChannels.size(); i++) {
    if(EXOMiscUtil::TypeOfChannel(fChannels[i]) == EXOMiscUtil::kAPDGang) {
      fFirstAPDChannelIndex = i;
      break;
    }
  }
  fNoiseColumnLength = fChannels.size() * (2*(fMaxF-fMinF) + 1);
  size_t FileNumChannels = NUMBER_READOUT_CHANNELS - 2*NCHANNEL_PER_WIREPLANE;

  // Then fill fNoiseCorrelations.
  // Note that we store the same-frequency blocks in column-major order,
  // to simplify GEMM calls (if BLAS is provided).
  // We also extract the diagonal entries, for the purpose of preconditioning.
  fNoiseCorrelations.resize(fMaxF - fMinF + 1);

  std::filebuf NoiseFile;
  fNoiseDiag.clear();
  fNoiseDiag.reserve(fNoiseColumnLength);
  assert(NoiseFile.open(fNoiseFilename.c_str(), std::ios_base::in | std::ios_base::binary) != NULL);
  assert(NoiseFile.pubseekoff(0, std::ios_base::beg, std::ios_base::in) == 0);
  assert(NoiseFile.pubseekoff(0, std::ios_base::end, std::ios_base::in) ==
         FileNumChannels*FileNumChannels*(4*1023+1)*sizeof(double));
  assert(sizeof(double) == 8);
  assert(fMinF == 1); // Else I'll need to generalize the code that produces these files.
  for(size_t f = fMinF; f <= fMaxF; f++) {
    bool IsFullBlock = (f != fMaxF);
    std::vector<double>& block = fNoiseCorrelations[f-fMinF];
    block.clear();
    block.reserve(fChannels.size()*fChannels.size()*(IsFullBlock ? 4 : 1));

    size_t FreqFilePos = (f-fMinF)*4*sizeof(double)*FileNumChannels*FileNumChannels;

    // Go ahead and fetch the entire block from the file.
    // A few rows/columns aren't needed (suppressed or bad channels),
    // but this prevents individual small queries which don't scale well.
    std::vector<double> FileBlock(FileNumChannels*FileNumChannels*(IsFullBlock ? 4 : 1), 0);
    assert(NoiseFile.pubseekpos(FreqFilePos, std::ios_base::in) == FreqFilePos);
    if(NoiseFile.sgetn((char*)&FileBlock[0], FileBlock.size()*sizeof(double)) !=
       FileBlock.size()*sizeof(double)) {
      std::cout<<"Failed to read out block corresponding to f = "<<f<<std::endl;
      LogEXOMsg("sgetn failed", EEAlert);
    }

    // Iterate through columns (real).
    for(size_t index1 = 0; index1 < fChannels.size(); index1++) {
      size_t NoiseFileIndex1 = ChannelIndexMap[index1];
      size_t ColumnPos = (IsFullBlock ? 2 : 1) * FileNumChannels * NoiseFileIndex1;

      // Start with the real rows.
      for(size_t index2 = 0; index2 < fChannels.size(); index2++) {
        if(not fUseWireAPDCorrelations) {
          if(EXOMiscUtil::TypeOfChannel(fChannels[index1]) !=
             EXOMiscUtil::TypeOfChannel(fChannels[index2])) {
            block.push_back(0);
            continue;
          }
        }
        size_t NoiseFileIndex2 = ChannelIndexMap[index2];
        size_t RowPos = ColumnPos + NoiseFileIndex2;
        block.push_back(FileBlock[RowPos]);
        if(index1 == index2) fNoiseDiag.push_back(FileBlock[RowPos]);
      } // Done with real rows of real column.

      // Now imaginary rows.
      if(not IsFullBlock) continue;
      ColumnPos += FileNumChannels;
      for(size_t index2 = 0; index2 < fChannels.size(); index2++) {
        if(not fUseWireAPDCorrelations) {
          if(EXOMiscUtil::TypeOfChannel(fChannels[index1]) !=
             EXOMiscUtil::TypeOfChannel(fChannels[index2])) {
            block.push_back(0);
            continue;
          }
        }
        size_t NoiseFileIndex2 = ChannelIndexMap[index2];
        size_t RowPos = ColumnPos + NoiseFileIndex2;
        block.push_back(FileBlock[RowPos]);
      } // Done with imaginary rows of real column.
    } // Done with real columns.

    // Now the imaginary columns.
    if(not IsFullBlock) continue;
    for(size_t index1 = 0; index1 < fChannels.size(); index1++) {
      size_t NoiseFileIndex1 = ChannelIndexMap[index1];
      size_t ColumnPos = 2*FileNumChannels*(FileNumChannels+NoiseFileIndex1);

      // Start with the real rows.
      for(size_t index2 = 0; index2 < fChannels.size(); index2++) {
        if(not fUseWireAPDCorrelations) {
          if(EXOMiscUtil::TypeOfChannel(fChannels[index1]) !=
             EXOMiscUtil::TypeOfChannel(fChannels[index2])) {
            block.push_back(0);
            continue;
          }
        }
        size_t NoiseFileIndex2 = ChannelIndexMap[index2];
        size_t RowPos = ColumnPos + NoiseFileIndex2;
        block.push_back(FileBlock[RowPos]);
      } // Done with real rows of imag column.

      // Now imaginary rows.
      if(not IsFullBlock) continue;
      ColumnPos += FileNumChannels;
      for(size_t index2 = 0; index2 < fChannels.size(); index2++) {
        if(not fUseWireAPDCorrelations) {
          if(EXOMiscUtil::TypeOfChannel(fChannels[index1]) !=
             EXOMiscUtil::TypeOfChannel(fChannels[index2])) {
            block.push_back(0);
            continue;
          }
        }
        size_t NoiseFileIndex2 = ChannelIndexMap[index2];
        size_t RowPos = ColumnPos + NoiseFileIndex2;
        block.push_back(FileBlock[RowPos]);
        if(index1 == index2) fNoiseDiag.push_back(FileBlock[RowPos]);
      } // Done with imaginary rows of imag column.
    } // End loop over imag columns (index1).
  } // End loop over frequencies.
  assert(fNoiseDiag.size() == fNoiseColumnLength);
  NoiseFile.close();

  // Precompute useful transformations of the noise diagonal.
  fInvSqrtNoiseDiag.resize(fNoiseDiag.size());
  for(size_t i = 0; i < fNoiseDiag.size(); i++) fInvSqrtNoiseDiag[i] = double(1)/std::sqrt(fNoiseDiag[i]);

  // Go ahead and precondition the noise matrices.  This should improve the accuracy of multiplications.
  // So, N -> D^(-1/2) N D^(-1/2).
  for(size_t f = fMinF; f <= fMaxF; f++) {
    size_t BlockSize = fChannels.size()*(f < fMaxF ? 2 : 1);
    size_t DiagIndex = 2*fChannels.size()*(f-fMinF);
    for(size_t row = 0; row < BlockSize; row++) {
      for(size_t col = 0; col < BlockSize; col++) {
        fNoiseCorrelations[f-fMinF][row + BlockSize*col] *= fInvSqrtNoiseDiag[DiagIndex + row];
        fNoiseCorrelations[f-fMinF][row + BlockSize*col] *= fInvSqrtNoiseDiag[DiagIndex + col];
      }
    }
  }

  fWatches["FillNoiseCorrelations"].Stop();
}

int EXORefitSignals::Initialize()
{
  // Open the lightmap file, and extract their information.
  // Create unshaped wire drift waveforms.
  // Also initialize our various timers.

  // Initialize counters.
  fNumEventsHandled = 0;
  fNumSignalsHandled = 0;
  fTotalIterationsDone = 0;

  std::string FullLightmapFilename = EXOMiscUtil::SearchForFile(fLightmapFilename);
  if(FullLightmapFilename == "") LogEXOMsg("Failed to find lightmap file.", EEAlert);
  TFile* LightmapFile = TFile::Open(FullLightmapFilename.c_str());

  // Get the list of active APDs.
  TArrayI* APDs = (TArrayI*)LightmapFile->GetObjectUnchecked("APDs");
  for(Int_t i = 0; i < APDs->GetSize(); i++) {
    Int_t gang = APDs->At(i);
    // Get the lightmaps.
    std::ostringstream lightmapname;
    lightmapname << "lightmap_" << std::setw(3) << std::setfill('0') << gang;
    std::string old_lightmap = lightmapname.str();
    std::string new_lightmap = old_lightmap + "_clone";
    fLightMaps[gang] = (TH3D*)LightmapFile->Get(old_lightmap.c_str())->Clone(new_lightmap.c_str());

    // Get the gainmaps.
    std::ostringstream gainmapname;
    gainmapname << "gainmap_" << std::setw(3) << std::setfill('0') << gang;
    std::string old_gainmap = gainmapname.str();
    std::string new_gainmap = old_gainmap + "_clone";
    fGainMaps[gang] = (TGraph*)LightmapFile->Get(old_gainmap.c_str())->Clone(new_gainmap.c_str());
  }
  delete LightmapFile;

  // If we're doing APDs only, we don't need to generate wire models.
  if(fAPDsOnly) return 0;

  // Create wire drift waveforms.
  // The results will be high-bandpass waveforms, with the deposit occurring at 256us.
  // The whole waveform has a length of 512us.
  EXODigitizeWires dig;
  dig.set_drift_velocity(0.171*CLHEP::cm/CLHEP::microsecond);
  dig.set_collection_drift_velocity(0.2041*CLHEP::cm/CLHEP::microsecond);
  dig.set_trigger_time(256*CLHEP::microsecond); // Shift later so deposit is at 256us.
  dig.set_nsample(512);
  EXOCoordinates coord(EXOMiscUtil::kUVCoordinates, 0.25*CHANNEL_WIDTH, 0, 100*CLHEP::mm, 0);
  EXOMCPixelatedChargeDeposit pixel;
  pixel.SetCoordinates(coord);
  pixel.fTotalEnergy = pixel.fTotalIonizationEnergy = 1; // 1 MeV -- we'll rescale this.
  dig.GenerateUnshapedSignals(&pixel);
  fWireDeposit = dig.GetDdata(NCHANNEL_PER_WIREPLANE/2);
  fWireInduction = dig.GetDdata(NCHANNEL_PER_WIREPLANE/2-1);
  fWireInduction += dig.GetDdata(NCHANNEL_PER_WIREPLANE/2+1);
  fWireInduction /= 2;

  // Shift to move deposit time to 256us.
  size_t ShiftNeeded = size_t(pixel.fWireHitTime/SAMPLE_TIME_HIGH_BANDWIDTH) - 256*BANDWIDTH_FACTOR;
  for(size_t i = 0; i < fWireDeposit.GetLength() - ShiftNeeded; i++) {
    fWireDeposit[i] = fWireDeposit[i+ShiftNeeded];
    fWireInduction[i] = fWireInduction[i+ShiftNeeded];
  }

  // Normalize the model waveforms -- not strictly necessary since they're unshaped, but
  // it makes it easier to modularize code later.
  double MaxVal = fWireDeposit.GetMaxValue();
  fWireDeposit /= MaxVal;
  fWireInduction /= MaxVal;

  return 0;
}

EXORefitSignals::~EXORefitSignals()
{
  // Print statistics and timing information.

  // Stop any watches which are still going.
  for(std::map<std::string, TStopwatch>::iterator it = fWatches.begin(); it != fWatches.end(); it++) {
    it->second.Stop();
  }

  std::cout<<fNumEventsHandled<<" events were handled by signal refitting."<<std::endl;
  std::cout<<"Those events contained a total of "<<fNumSignalsHandled<<" signals to refit."<<std::endl;
  std::cout<<fTotalIterationsDone<<" iterations were required."<<std::endl;
  std::cout<<"Timings follow."<<std::endl;

  // Print all watch information.
  for(std::map<std::string, TStopwatch>::iterator it = fWatches.begin(); it != fWatches.end(); it++) {
    std::cout<<it->first<<":"<<std::endl;
    it->second.Print();
  }
}

EXOWaveformFT EXORefitSignals::GetModelForTime(double time) const
{
  // Return an EXOWaveformFT corresponding to a scintillation signal at time T.
  // The magnitude should be normalized so peak-baseline = 1.
  // The baseline itself is zero.
  // Note that at the moment, this assumes a waveform of length 2048 is required.
  // time is in ns.
  // No accounting for APD-by-APD shaping time variations is currently made.
  //
  // It might seem reasonable to do this just once, and apply a time shift in fourier space.
  // However, generating it in real space allows us to deal with signals near the end of
  // the trace, where periodicity is violated.
  // There is still some potential for caching the time-domain waveform, though, if needed.

  EXODoubleWaveform timeModel_fine;
  int refinedFactor = 5;
  timeModel_fine.SetLength(2048*refinedFactor);
  timeModel_fine.SetSamplingFreq(refinedFactor*CLHEP::megahertz);
  timeModel_fine.Zero();
  size_t NonzeroIndex = size_t(time/(CLHEP::microsecond/refinedFactor));
  for(size_t i = NonzeroIndex; i < timeModel_fine.GetLength(); i++) timeModel_fine[i] = 1;

  EXOTransferFunction tf;
  tf.AddIntegStageWithTime(3.*CLHEP::microsecond);
  tf.AddIntegStageWithTime(3.*CLHEP::microsecond);
  tf.AddDiffStageWithTime(10.*CLHEP::microsecond);
  tf.AddDiffStageWithTime(10.*CLHEP::microsecond);
  tf.AddDiffStageWithTime(300.*CLHEP::microsecond);

  tf.Transform(&timeModel_fine);
  timeModel_fine /= tf.GetGain();

  EXODoubleWaveform timeModel;
  timeModel.SetLength(2048);
  for(size_t i = 0; i < timeModel.GetLength(); i++) timeModel[i] = timeModel_fine[i*refinedFactor];

  EXOWaveformFT fwf;
  EXOFastFourierTransformFFTW::GetFFT(timeModel.GetLength()).PerformFFT(timeModel, fwf);
  assert(fwf.GetLength() == 1025); // Just to make sure I'm reasoning properly.
  return fwf;
}

double EXORefitSignals::GetGain(unsigned char channel, EventHandler& event) const
{
  // Return the gain of an apd channel.  This is the conversion factor from number
  // of photons incident on the APD to number of ADC counts (peak-baseline) in the
  // digitized signal.
  // It's a rough estimate, since this number isn't well-known, but we only need it to
  // set the scale for how important Poisson noise is, relative to electronic noise.
  // We currently use laser data from run 4540, and extract time-dependence from the
  // gainmap (the time-dependence of the lightmap).  It would be interesting to
  // do a fit of laser and gainmap data for the times when they overlap,
  // and get a higher-quality set of values.

  double Gain = 1.9; // 1.9 electron-hole pairs per photon, on average.

  // APD gains from the laser run 4540.
  switch(channel) {
    case 152: Gain *= 201.230438146; break;
    case 153: Gain *= 178.750438779; break;
    case 154: Gain *= 194.228589338; break;
    case 155: Gain *= 183.33801615; break;
    case 156: Gain *= 218.485999976; break;
    case 157: Gain *= 222.139259152; break;
    case 158: Gain *= 169.982559736; break;
    case 159: Gain *= 140.385120552; break;
    case 160: Gain *= 137.602725389; break;
    case 161: Gain *= 197.78183714; break;
    case 162: Gain *= 155.478773762; break;
    // case 163: Gain *= 0; // Bad channel, omitted.
    case 164: Gain *= 175.875067527; break;
    case 165: Gain *= 160.014408865; break;
    case 166: Gain *= 183.408055613; break;
    case 167: Gain *= 189.600819126; break;
    case 168: Gain *= 160.339214431; break;
    case 169: Gain *= 168.547991045; break;
    case 170: Gain *= 182.670039836; break;
    case 171: Gain *= 205.567802982; break;
    case 172: Gain *= 195.87450621; break;
    case 173: Gain *= 224.956647122; break;
    case 174: Gain *= 232.062359991; break;
    case 175: Gain *= 241.822881767; break;
    case 176: Gain *= 194.740435753; break;
    case 177: Gain *= 189.867775084; break;
    // case 178: Gain *= 0; // Bad channel, omitted.
    case 179: Gain *= 206.755206938; break;
    case 180: Gain *= 207.822617603; break;
    case 181: Gain *= 207.501985741; break;
    case 182: Gain *= 218.213137769; break;
    case 183: Gain *= 234.369354843; break;
    case 184: Gain *= 99.908111992; break;
    case 185: Gain *= 238.381809313; break;
    case 186: Gain *= 225.118270743; break;
    case 187: Gain *= 199.078450518; break;
    case 188: Gain *= 221.863823239; break;
    case 189: Gain *= 177.032783679; break;
    case 190: Gain *= 196.787332164; break;
    // case 191: Gain *= 0; // Bad channel, omitted.
    case 192: Gain *= 194.923448865; break;
    case 193: Gain *= 197.027984846; break;
    case 194: Gain *= 202.757086104; break;
    case 195: Gain *= 194.432937658; break;
    case 196: Gain *= 208.992809367; break;
    case 197: Gain *= 224.762562055; break;
    case 198: Gain *= 217.696006443; break;
    case 199: Gain *= 222.380158829; break;
    case 200: Gain *= 218.358804472; break;
    case 201: Gain *= 209.573057132; break;
    case 202: Gain *= 194.684536629; break;
    case 203: Gain *= 182.543842783; break;
    case 204: Gain *= 193.469930111; break;
    // case 205: Gain *= 0; // Bad channel, omitted.
    case 206: Gain *= 193.627191472; break;
    case 207: Gain *= 196.073150574; break;
    case 208: Gain *= 189.597962521; break;
    case 209: Gain *= 198.824317108; break;
    case 210: Gain *= 222.747770671; break;
    case 211: Gain *= 216.928470825; break;
    case 212: Gain *= 223.437239807; break;
    case 213: Gain *= 224.316404923; break;
    case 214: Gain *= 216.26783603; break;
    case 215: Gain *= 209.612423384; break;
    case 216: Gain *= 223.041660884; break;
    case 217: Gain *= 202.642254512; break;
    case 218: Gain *= 213.904993632; break;
    case 219: Gain *= 221.988942321; break;
    case 220: Gain *= 201.427174798; break;
    case 221: Gain *= 196.689200146; break;
    case 222: Gain *= 191.457656123; break;
    case 223: Gain *= 186.183873541; break;
    case 224: Gain *= 217.033080346; break;
    case 225: Gain *= 205.858374653; break;
    default: Gain *= 0; // Bad or non-existent channel.
  }
  // Time-dependence from the gainmap.
  const TGraph* GainGraph = fGainMaps.at(channel);
  Gain *= GainGraph->Eval(event.fUnixTimeOfEvent)/GainGraph->Eval(1355409118.254096);

  Gain *= 32.e-9; // Convert from electrons to volts in the preamp. Roughly 1/(5 pF) gain.
  Gain *= 12.10; // Gain from shapers (amplification factor, and gain from transfer function.
  Gain *= 4096./2.5; // Conversion from volts to ADC counts -- full-scale is 2.5 volts.

  return Gain;
}

std::vector<double> EXORefitSignals::MakeWireModel(EXODoubleWaveform& in,
                                                   const EXOTransferFunction& transfer,
                                                   const double Gain,
                                                   const double Time) const
{
  // Helper function for dealing with shaping and FFT of wire models.
  EXODoubleWaveform shapedIn;
  transfer.Transform(&in, &shapedIn);
  shapedIn /= Gain;

  EXODoubleWaveform wf;
  wf.SetLength(2048);
  wf.Zero();
  for(size_t i = 0; i < 2048; i++) {
    double RelTime = SAMPLE_TIME*i - Time;
    int HighBandwidthIndex = int((256.*CLHEP::microsecond + RelTime)/SAMPLE_TIME_HIGH_BANDWIDTH);
    if(HighBandwidthIndex >= 0 and HighBandwidthIndex < int(shapedIn.GetLength())) {
      wf[i] = shapedIn[HighBandwidthIndex];
    }
  }

  EXOWaveformFT fwf;
  EXOFastFourierTransformFFTW::GetFFT(2048).PerformFFT(wf, fwf);

  std::vector<double> out;
  out.resize(2*1024-1);
  for(size_t f = 1; f <= 1024; f++) {
    out[2*(f-1)] = fwf[f].real();
    if(f != 1024) out[2*(f-1)+1] = fwf[f].imag();
  }

  return out;
}

void EXORefitSignals::AcceptEvent(EXOEventData* ED, Long64_t entryNum)
{
  // Push in one more event to handle; return whatever events are able to finish.
  // Do whatever matrix multiplications are now warranted.

  EventHandler* event = new EventHandler;
  event->fEntryNumber = entryNum;

  // If we don't have previously-established scintillation times, we can't do anything -- skip.
  if(ED->GetNumScintillationClusters() == 0) {
    FinishEvent(event);
    return;
  }

  // If the waveforms aren't full-length, skip for now (although we should be able to handle them later).
  if(ED->fEventHeader.fSampleCount != 2047) {
    FinishEvent(event);
    return;
  }

  // For now, we also only deal with events containing *exactly* one scintillation cluster.
  // There's nothing theoretical that warrants this; it's just easier to code up.
  if(ED->GetNumScintillationClusters() != 1) {
    FinishEvent(event);
    return;
  }
  EXOScintillationCluster* scint = ED->GetScintillationCluster(0);

  // If there are no fully-reconstructed clusters, then we can't do anything -- so, skip them too.
  // Otherwise, extract a list of clusters for future convenience.
  std::vector<EXOChargeCluster*> FullClusters;
  for(size_t i = 0; i < scint->GetNumChargeClusters(); i++) {
    EXOChargeCluster* clu = scint->GetChargeClusterAt(i);
    if(std::abs(clu->fX) > 200 or std::abs(clu->fY) > 200 or std::abs(clu->fZ) > 200) continue;
    if(clu->fPurityCorrectedEnergy < 1) continue;
    FullClusters.push_back(clu);
  }
  if(FullClusters.empty()) {
    FinishEvent(event);
    return;
  }

  // If necessary, extract the noise correlations object with a proper ordering.
  FillNoiseCorrelations(*ED);

  // Save the unix time of the event (as a double, since ROOT will convert it anyway).
  event->fUnixTimeOfEvent = double(ED->fEventHeader.fTriggerSeconds);
  event->fUnixTimeOfEvent += double(ED->fEventHeader.fTriggerMicroSeconds)/1e6;

  // Given the positions of the clusters, estimate how the light should be distributed among gangs.
  // ExpectedYieldPerGang will be the expected peak-baseline (ADC counts) of a 2615 keV event.
  event->fExpectedYieldPerGang.clear();
  event->fExpectedEnergy_keV = 0;
  for(size_t i = fFirstAPDChannelIndex; i < fChannels.size(); i++) {
    event->fExpectedYieldPerGang[fChannels[i]] = 0;
  }
  for(size_t i = 0; i < FullClusters.size(); i++) {
    EXOChargeCluster* clu = FullClusters[i];
    event->fExpectedEnergy_keV += clu->fPurityCorrectedEnergy;
    for(size_t j = fFirstAPDChannelIndex; j < fChannels.size(); j++) {
      unsigned char gang = fChannels[j];
      Double_t GainFuncVal = fGainMaps[gang]->Eval(event->fUnixTimeOfEvent);

      // Make sure cluster is in the proper range for interpolation -- else return 0.
      Double_t LightMapVal;
      TAxis* Xaxis = fLightMaps[gang]->GetXaxis();
      TAxis* Yaxis = fLightMaps[gang]->GetYaxis();
      TAxis* Zaxis = fLightMaps[gang]->GetZaxis();
      if(Xaxis->GetBinCenter(1) <= clu->fX and clu->fX < Xaxis->GetBinCenter(Xaxis->GetNbins()) and
         Yaxis->GetBinCenter(1) <= clu->fY and clu->fY < Yaxis->GetBinCenter(Yaxis->GetNbins()) and
         Zaxis->GetBinCenter(1) <= clu->fZ and clu->fZ < Zaxis->GetBinCenter(Zaxis->GetNbins())) {
        LightMapVal = fLightMaps[gang]->Interpolate(clu->fX, clu->fY, clu->fZ);
      }
      else {
        // Interpolate would return 0, and I'm actually OK with that -- but I still want to kill the warning.
        LightMapVal = 0;
      }

      event->fExpectedYieldPerGang[gang] += LightMapVal*GainFuncVal*clu->fPurityCorrectedEnergy;
    }
  }
  // We just want to weight the clusters appropriately when we guess where light should be collected.
  // Divide out to ensure that at the end, a result of 1 corresponds to a 2615 keV event (roughly).
  for(size_t i = fFirstAPDChannelIndex; i < fChannels.size(); i++) {
    event->fExpectedYieldPerGang[fChannels[i]] /= event->fExpectedEnergy_keV;
  }

  // If we don't expect any yield, then clearly there will be a degenerate matrix.
  // So, instead drop such events.
  // (Specifically, if a 2615 keV event would produce less than 1ADC on every gang, drop it.)
  bool HasYield = false;
  for(size_t i = fFirstAPDChannelIndex; i < fChannels.size(); i++) {
    if(event->fExpectedYieldPerGang[fChannels[i]] > 1) HasYield = true;
  }
  if(not HasYield) {
    FinishEvent(event);
    return;
  }

  fWatches["GenerateExpectedSignals"].Start(false);
  // Generate the expected light signal shape (normalized), given the time of the scintillation.
  // Alternate between real and imaginary parts, mimicking the variable ordering we use throughout.
  // Also drop the zero-frequency component (which isn't used)
  // and the last imaginary component (which is identically zero).
  EXOWaveformFT modelFT = GetModelForTime(scint->fTime);
  event->fmodel_realimag.resize(2*modelFT.GetLength() - 3);
  for(size_t i = 1; i < modelFT.GetLength(); i++) {
    event->fmodel_realimag[2*i - 2] = modelFT[i].real();
    if(i != modelFT.GetLength()-1) event->fmodel_realimag[2*i - 1] = modelFT[i].imag();
  }

  // Now produce the expected wire signals.
  event->fWireModel.clear();
  if(not fAPDsOnly) {
    std::set<size_t> UWireSignals;
    EXOElectronicsShapers* electronicsShapers = GetCalibrationFor(EXOElectronicsShapers,
                                                                  EXOElectronicsShapersHandler,
                                                                  "timevartau",
                                                                  ED->fEventHeader);
    const EXOUWireGains* GainsFromDatabase = GetCalibrationFor(EXOUWireGains,
                                                               EXOUWireGainsHandler,
                                                               "source_calibration",
                                                               ED->fEventHeader);
    if(not electronicsShapers or not GainsFromDatabase) {
      LogEXOMsg("Unable to get necessary information from DB", EEAlert);
    }
    for(size_t i = 0; i < ED->GetNumUWireSignals(); i++) {
      EXOUWireSignal* sig = ED->GetUWireSignal(i);
      if(sig->fIsInduction) continue;
      UWireSignals.insert(i);
    }
    for(std::set<size_t>::iterator sigIt = UWireSignals.begin();
        sigIt != UWireSignals.end();
        sigIt++) {
      EXOUWireSignal* sig = ED->GetUWireSignal(*sigIt);

      std::map<unsigned char, std::vector<double> > ModelForThisSignal;

      // Deposit channel.
      const EXOTransferFunction& transferDep =
        electronicsShapers->GetTransferFunctionForChannel(sig->fChannel);
      double Gain = transferDep.GetGain();
      double DepChanGain = GainsFromDatabase->GetGainOnChannel(sig->fChannel);
      ModelForThisSignal[sig->fChannel] = MakeWireModel(fWireDeposit,
                                                        transferDep,
                                                        Gain,
                                                        sig->fTime);

      if(EXOMiscUtil::TypeOfChannel(sig->fChannel-1) == EXOMiscUtil::kUWire) {
        const EXOTransferFunction& transferInd =
          electronicsShapers->GetTransferFunctionForChannel(sig->fChannel-1);
        double ThisChanGain = Gain * GainsFromDatabase->GetGainOnChannel(sig->fChannel-1)/DepChanGain;
        ModelForThisSignal[sig->fChannel-1] = MakeWireModel(fWireInduction,
                                                            transferInd,
                                                            ThisChanGain,
                                                            sig->fTime);
      }

      if(EXOMiscUtil::TypeOfChannel(sig->fChannel+1) == EXOMiscUtil::kUWire) {
        const EXOTransferFunction& transferInd =
          electronicsShapers->GetTransferFunctionForChannel(sig->fChannel+1);
        double ThisChanGain = Gain * GainsFromDatabase->GetGainOnChannel(sig->fChannel+1)/DepChanGain;
        ModelForThisSignal[sig->fChannel+1] = MakeWireModel(fWireInduction,
                                                            transferInd,
                                                            ThisChanGain,
                                                            sig->fTime);
      }

      event->fWireModel.push_back(std::make_pair(*sigIt, ModelForThisSignal));
    } // End loop over u-wire signals.
  } // (which we only did if we're handling wire signals.
  fWatches["GenerateExpectedSignals"].Stop();

  // For convenience, store the column length we'll be dealing with.
  event->fColumnLength = 2*fChannels.size()*(fMaxF-fMinF) + fChannels.size() + event->fWireModel.size() + 1;

  // We can find the appropriate preconditioner here.
  // This is a pretty good preconditioner, obtained by approximating A ~ {{D L} {trans(L) 0}},
  // where D is diagonal.  Thus, the approximation comes from ignoring noise cross-terms and
  // Poisson noise terms.
  // Haven't decided yet whether it's important to include Poisson terms on the diagonal; currently I don't.
  // Find X using trans(X)X = trans(L) D^(-1) L.
  std::vector<double> Temp1(event->fColumnLength * (event->fWireModel.size()+1), 0);
  std::vector<double> Temp2(event->fColumnLength * (event->fWireModel.size()+1), 0);
  for(size_t i = 0; i < event->fWireModel.size()+1; i++) {
    Temp1[i*event->fColumnLength + fNoiseColumnLength + i] = 1;
  } // Temp1 = {{0}{I}}
  DoLagrangeAndConstraintMul<'L'>(Temp1, Temp2, *event); // Temp2 = {{D^(-1/2)L} {0}}
  Temp1.assign(event->fColumnLength * (event->fWireModel.size()+1), 0);
  DoLagrangeAndConstraintMul<'C'>(Temp2, Temp1, *event); // Temp1 = {{0} {trans(L)D^(-1)L}}
  // Here I could produce a packed version; but I don't think the performance boost will be significant.
  // Remember X is a small matrix.
  // Produce the Cholesky decomposition, and store it in event->fPreconX.
  lapack_int ret;
  ret = LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'U', event->fWireModel.size()+1,
                       &Temp1[fNoiseColumnLength], event->fColumnLength);
  if(ret != 0) LogEXOMsg("Failed factorization", EEAlert);
  event->fPreconX.assign((event->fWireModel.size()+1)*(event->fWireModel.size()+1), 0);
  for(size_t i = 0; i < event->fWireModel.size()+1; i++) {
    for(size_t j = 0; j <= i; j++) {
      event->fPreconX[i*(event->fWireModel.size()+1) + j] =
        Temp1[i*event->fColumnLength + fNoiseColumnLength + j];
    }
  }

  // Give a really nice initial guess for X, obtained by solving exactly
  // with the approximate version of the matrix used for preconditioning.
  // Since the RHS only has non-zero entries in the lower square, simplifications are used.
  event->fX.assign(event->fColumnLength * (event->fWireModel.size()+1), 0);
  for(size_t i = 0; i <= event->fWireModel.size(); i++) {
    size_t Index = (i+1)*event->fColumnLength; // Next column; then subtract.
    Index -= event->fWireModel.size() + 1; // Step backward.
    Index += i; // Go forward to the right entry.
    event->fX[Index] = 1; // All models are normalized to 1.
  }
  event->fX = DoInvLPrecon(event->fX, *event);
  event->fprecon_tmp = DoInvRPrecon(event->fX, *event); // Beginning of multiplying by matrix.

  // Request a matrix multiplication of X.
  event->fResultIndex = RequestNoiseMul(event->fprecon_tmp, event->fColumnLength);
  fNumEventsHandled++; // One more event that will be actually handled.
  fNumSignalsHandled += event->fWireModel.size() + 1;

  // Push event onto the list of event handlers.
#ifdef USE_LOCKFREE
  assert(fEventHandlerQueue.unsynchronized_push(event));
#else
  fEventHandlerQueue.push(event);
#endif

#ifdef USE_LOCKFREE
  // Ensure that fEventHandlerResults has enough nodes reserved to accept all of the queued events.
  // Unfortunately, the only way I know to do this is really awkward.
  assert(fEventHandlerResults.empty());
  while(not fEventHandlerQueue.empty()) {
    EventHandler* evt = NULL;
    assert(fEventHandlerQueue.unsynchronized_pop(evt));
    assert(fEventHandlerResults.unsynchronized_push(evt));
  }
  while(not fEventHandlerResults.empty()) {
    EventHandler* evt = NULL;
    assert(fEventHandlerResults.unsynchronized_pop(evt));
    assert(fEventHandlerQueue.unsynchronized_push(evt));
  }
#endif

  // Now, while there are enough requests in the queue, satisfy those requests.
  while(fNumVectorsInQueue > 40) DoPassThroughEvents();
}

void EXORefitSignals::FlushEvents()
{
  // Finish processing for all events in the event handler list,
  // regardless of how many pending multiplication requests are queued.
  while(not fEventHandlerQueue.empty()) DoPassThroughEvents();
}

void EXORefitSignals::FinishEvent(EventHandler* event)
{
  // Compute and fill denoised signals, as appropriate.
  // Then pass the filled event to the output module.

#ifdef USE_THREADS
  FinishEventMutex.lock(); // Wait until we get a lock.
#endif

  EXOEventData* ED = fInputModule.GetEvent(event->fEntryNumber);

  // We need to clear out the denoised information here, since we just freshly read the event from file.
  for(size_t i = 0; i < ED->GetNumScintillationClusters(); i++) {
    ED->GetScintillationCluster(i)->fDenoisedEnergy = 0;
  }
  for(size_t i = 0; i < ED->GetNumUWireSignals(); i++) {
    ED->GetUWireSignal(i)->fDenoisedEnergy = 0;
  }
  for(size_t i = 0; i < ED->GetNumChargeClusters(); i++) {
    ED->GetChargeCluster(i)->fDenoisedEnergy = 0;
  }

  if(not event->fX.empty()) {
    // Undo preconditioning of X.
    event->fX = DoInvRPrecon(event->fX, *event);
    for(size_t i = 0; i < event->fX.size(); i++) {
      size_t imod = i % event->fColumnLength;
      if(imod < fNoiseColumnLength) event->fX[i] *= fInvSqrtNoiseDiag[imod];
    }

    // We need to compute denoised signals.
    fWFTree.GetEntryWithIndex(ED->fRunNumber, ED->fEventNumber);
    fWFEvent->GetWaveformData()->Decompress();

    // Collect the fourier-transformed waveforms.  Save them split into real and complex parts.
    std::vector<EXODoubleWaveform> WF_real, WF_imag;
    for(size_t i = 0; i < fChannels.size(); i++) {
      const EXOWaveform* wf = fWFEvent->GetWaveformData()->GetWaveformWithChannel(fChannels[i]);

      // Take the Fourier transform.
      EXODoubleWaveform dwf = wf->Convert<Double_t>();
      EXOWaveformFT fwf;
      EXOFastFourierTransformFFTW::GetFFT(dwf.GetLength()).PerformFFT(dwf, fwf);

      // Extract the real part.
      EXODoubleWaveform rwf;
      rwf.SetLength(fwf.GetLength());
      for(size_t f = 0; f < fwf.GetLength(); f++) rwf[f] = fwf[f].real();
      WF_real.push_back(rwf);

      // Extract the imaginary part.
      // Note that the first and last components are strictly real (though we get them anyway).
      EXODoubleWaveform iwf;
      iwf.SetLength(fwf.GetLength());
      for(size_t f = 0; f < fwf.GetLength(); f++) iwf[f] = fwf[f].imag();
      WF_imag.push_back(iwf);
    }

    // Produce estimates of the signals.
    std::vector<double> Results(event->fWireModel.size()+1, 0);
    for(size_t i = 0; i < Results.size(); i++) {
      for(size_t f = 0; f <= fMaxF-fMinF; f++) {
        for(size_t chan_index = 0; chan_index < fChannels.size(); chan_index++) {
          size_t XIndex = event->fColumnLength*i + 2*fChannels.size()*f + chan_index;
          Results[i] += event->fX[XIndex]*WF_real[chan_index][f+fMinF];
        }
        if(f == fMaxF-fMinF) continue;
        for(size_t chan_index = 0; chan_index < fChannels.size(); chan_index++) {
          size_t XIndex = event->fColumnLength*i + 2*fChannels.size()*f + fChannels.size() + chan_index;
          Results[i] += event->fX[XIndex]*WF_imag[chan_index][f+fMinF];
        }
      }
    }

    // Translate signal magnitudes into corresponding objects.
    if(not fAPDsOnly) {
      for(size_t i = 0; i < event->fWireModel.size(); i++) {
        size_t sigIndex = event->fWireModel[i].first;
        EXOUWireSignal* sig = ED->GetUWireSignal(sigIndex);
        double UWireScalingFactor = ADC_FULL_SCALE_ELECTRONS_WIRE * W_VALUE_LXE_EV_PER_ELECTRON /
                                    (CLHEP::keV * ADC_BITS);
        sig->fDenoisedEnergy = Results[i]*UWireScalingFactor;
      }
    }
    ED->GetScintillationCluster(0)->fDenoisedEnergy = Results.back()*fThoriumEnergy_keV;
  } // End setting of denoised energy signals.

  fOutputModule.ProcessEvent(ED);

#ifdef USE_THREADS
  FinishEventMutex.unlock(); // Release access to tinput and toutput for other threads.
#endif

  delete event;
}

bool EXORefitSignals::DoBlBiCGSTAB(EventHandler& event)
{
  // Pick up wherever we left off.
  // This function gets called when a noise matrix multiplication just happened.
  // So, we have to first identify where we are, then proceed as far as we can until:
  //   Another matrix multiplication needs to be done, or
  //   The solver has terminated.
  // If the solver terminates, return true; otherwise, return false.
  // For more information on the Block-BiCGSTAB algorithm, see:
  // Electronic Transactions on Numerical Analysis, vol 16, 129-142 (2003).
  // "A BLOCK VERSION OF BICGSTAB FOR LINEAR SYSTEMS WITH MULTIPLE RIGHT-HAND SIDES"
  // A. EL GUENNOUNI, K. JBILOU, AND H. SADOK.
  lapack_int ret;

  if(event.fR.size() == 0) {
    // We're still in the setup phase.
    // Start by copying result into R.
    FillFromNoise(event.fR, event.fWireModel.size()+1, event.fColumnLength, event.fResultIndex);
    // Now need to finish multiplying by A, accounting for the other terms.
    DoRestOfMultiplication(event.fprecon_tmp, event.fR, event);
    // Now, R <-- B - R = B - AX.
    for(size_t i = 0; i < event.fR.size(); i++) event.fR[i] = -event.fR[i];
    for(size_t i = 0; i < event.fWireModel.size()+1; i++) {
      event.fR[i*event.fColumnLength + fNoiseColumnLength + i] += 1;
    }
    // Now precondition R appropriately.
    event.fR = DoInvLPrecon(event.fR, event);
    // Might as well check if we can terminate right off the bat.  Not impossible!
    if(CanTerminate(event)) return true;
    // Set up other pieces of the handler.
    event.fP = event.fR;
    event.fR0hat = event.fR;
    // We'll take K1_inv A K2_inv P, so compute K2_inv P.
    event.fprecon_tmp = DoInvRPrecon(event.fP, event);
    // Now we want V <- AP, so request a multiplication by P.
    event.fResultIndex = RequestNoiseMul(event.fprecon_tmp, event.fColumnLength);
    return false;
  }
  else if(event.fV.size() == 0) {
    // At the beginning of the iteration, we just computed V = AP.
    fTotalIterationsDone++;
    FillFromNoise(event.fV, event.fWireModel.size()+1, event.fColumnLength, event.fResultIndex);
    // Now need to finish multiplying by A, accounting for the other terms.
    DoRestOfMultiplication(event.fprecon_tmp, event.fV, event);
    event.fV = DoInvLPrecon(event.fV, event);
    // Factorize fR0hat*V, so that we can solve equations using it twice.
    event.fR0hat_V_factors.assign((event.fWireModel.size()+1)*(event.fWireModel.size()+1), 0);
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                event.fWireModel.size() + 1, event.fWireModel.size() + 1, event.fColumnLength,
                1, &event.fR0hat[0], event.fColumnLength, &event.fV[0], event.fColumnLength,
                0, &event.fR0hat_V_factors[0], event.fWireModel.size() + 1);
    event.fR0hat_V_pivot.resize(event.fWireModel.size()+1);
    ret = LAPACKE_dgetrf(LAPACK_COL_MAJOR, event.fWireModel.size()+1, event.fWireModel.size()+1,
                         &event.fR0hat_V_factors[0], event.fWireModel.size()+1,
                         &event.fR0hat_V_pivot[0]);
    if(ret != 0) LogEXOMsg("Factorization failed", EEAlert);
    // Now compute alpha.
    event.fAlpha.assign((event.fWireModel.size()+1)*(event.fWireModel.size()+1), 0);
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                event.fWireModel.size() + 1, event.fWireModel.size() + 1, event.fColumnLength,
                1, &event.fR0hat[0], event.fColumnLength, &event.fR[0], event.fColumnLength,
                0, &event.fAlpha[0], event.fWireModel.size() + 1);
    ret = LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N',
                         event.fWireModel.size()+1, event.fWireModel.size()+1,
                         &event.fR0hat_V_factors[0], event.fWireModel.size()+1,
                         &event.fR0hat_V_pivot[0],
                         &event.fAlpha[0], event.fWireModel.size()+1);
    if(ret != 0) LogEXOMsg("Solving failed", EEAlert);
    // Update R <-- R - V*alpha.
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                event.fColumnLength, event.fWireModel.size() + 1, event.fWireModel.size() + 1,
                -1, &event.fV[0], event.fColumnLength, &event.fAlpha[0], event.fWireModel.size() + 1,
                1, &event.fR[0], event.fColumnLength);
    // Now we desire T = AR (AS in paper).  Request a matrix multiplication, and return.
    // Remember to apply preconditioner here too.
    event.fprecon_tmp = DoInvRPrecon(event.fR, event);
    event.fResultIndex = RequestNoiseMul(event.fprecon_tmp, event.fColumnLength);
    return false;
  }
  else {
    // We're in the second half of the iteration, where T was just computed.
    std::vector<double> T;
    FillFromNoise(T, event.fWireModel.size()+1, event.fColumnLength, event.fResultIndex);
    // Now need to finish multiplying by A, accounting for the other terms.
    DoRestOfMultiplication(event.fprecon_tmp, T, event);
    // Finish preconditioner.
    T = DoInvLPrecon(T, event);
    // Compute omega.
    double omega = std::inner_product(T.begin(), T.end(), event.fR.begin(), double(0)) /
                   std::inner_product(T.begin(), T.end(), T.begin(), double(0));
    // Modify X.
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                event.fColumnLength, event.fWireModel.size() + 1, event.fWireModel.size() + 1,
                1, &event.fP[0], event.fColumnLength, &event.fAlpha[0], event.fWireModel.size() + 1,
                1, &event.fX[0], event.fColumnLength);
    for(size_t i = 0; i < event.fX.size(); i++) event.fX[i] += omega*event.fR[i];
    // Modify R.
    for(size_t i = 0; i < event.fR.size(); i++) event.fR[i] -= omega*T[i];
    // Check if we should conclude here.
    if(CanTerminate(event)) return true;
    // Now compute beta, solving R0hat_V beta = -R0hat_T
    std::vector<double> Beta((event.fWireModel.size()+1)*(event.fWireModel.size()+1), 0);
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                event.fWireModel.size() + 1, event.fWireModel.size() + 1, event.fColumnLength,
                -1, &event.fR0hat[0], event.fColumnLength, &T[0], event.fColumnLength,
                0, &Beta[0], event.fWireModel.size() + 1);
    ret = LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N',
                         event.fWireModel.size()+1, event.fWireModel.size()+1,
                         &event.fR0hat_V_factors[0], event.fWireModel.size()+1,
                         &event.fR0hat_V_pivot[0],
                         &Beta[0], event.fWireModel.size()+1);
    if(ret != 0) LogEXOMsg("Solving failed", EEAlert);
    // Update P.  Overwrite T for temporary work.
    T = event.fR;
    for(size_t i = 0; i < event.fP.size(); i++) event.fP[i] -= omega*event.fV[i];
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                event.fColumnLength, event.fWireModel.size() + 1, event.fWireModel.size() + 1,
                1, &event.fP[0], event.fColumnLength, &Beta[0], event.fWireModel.size() + 1,
                1, &T[0], event.fColumnLength);
    std::swap(T, event.fP);
    // Clear vectors in event which are no longer needed -- this helps us keep track of where we are.
    event.fV.clear();
    event.fAlpha.clear();
    event.fR0hat_V_factors.clear();
    event.fR0hat_V_pivot.clear();
    // And request AP.  Remember preconditioner.
    event.fprecon_tmp = DoInvRPrecon(event.fP, event);
    event.fResultIndex = RequestNoiseMul(event.fprecon_tmp, event.fColumnLength);
    return false;
  }
}

void EXORefitSignals::DoRestOfMultiplication(const std::vector<double>& in,
                                             std::vector<double>& out,
                                             EventHandler& event)
{
  // After noise terms have already been handled, deal with all of the others.
  // This should not be the bottleneck.
  DoPoissonMultiplication(in, out, event);
  DoLagrangeAndConstraintMul<'A'>(in, out, event);
}

void EXORefitSignals::DoPoissonMultiplication(const std::vector<double>& in,
                                              std::vector<double>& out,
                                              EventHandler& event)
{
  // Poisson terms for APD channels.
  for(size_t k = fFirstAPDChannelIndex; k < fChannels.size(); k++) { // APD gangs
    double ChannelFactors = event.fExpectedEnergy_keV/fThoriumEnergy_keV;
    ChannelFactors *= GetGain(fChannels[k], event);
    ChannelFactors *= event.fExpectedYieldPerGang.at(fChannels[k]);

    for(size_t n = 0; n <= event.fWireModel.size(); n++) { // signals
      // Compute the factors common to all frequencies.
      double CommonFactor = 0;
      for(size_t g = 0; g <= fMaxF - fMinF; g++) {
        size_t DiagIndex = 2*fChannels.size()*g + k;
        size_t InIndex = n*event.fColumnLength + DiagIndex;
        CommonFactor += event.fmodel_realimag[2*g] * in[InIndex] *
                        fInvSqrtNoiseDiag[DiagIndex];
        if(g < fMaxF-fMinF) CommonFactor += event.fmodel_realimag[2*g+1] * in[InIndex+fChannels.size()] *
                                            fInvSqrtNoiseDiag[DiagIndex+fChannels.size()];
      }
      CommonFactor *= ChannelFactors;

      // Now actually transfer the changes to the out vector.
      for(size_t f = 0; f <= fMaxF - fMinF; f++) {
        size_t DiagIndex = 2*fChannels.size()*f + k;
        size_t OutIndex = n*event.fColumnLength + DiagIndex;
        out[OutIndex] += CommonFactor * event.fmodel_realimag[2*f] *
                         fInvSqrtNoiseDiag[DiagIndex];
        if(f < fMaxF-fMinF) out[OutIndex+fChannels.size()] += CommonFactor * event.fmodel_realimag[2*f+1] *
                                               fInvSqrtNoiseDiag[DiagIndex+fChannels.size()];
      }
    }
  } // End Poisson terms.
}

void EXORefitSignals::DoNoiseMultiplication()
{
  // Multiply everything in fNoiseMulQueue by the noise.
  // Note that we expect columns in the input to contain only the noise portion, not the constraint rows;
  // otherwise, the vector lengths would not match.
  // The result is placed in fNoiseMulResult.
  assert(fNoiseMulQueue.size() == fNoiseColumnLength * fNumVectorsInQueue);
  fNoiseMulResult.assign(fNoiseMulQueue.size(), 0); // Probably don't need to fill with 0.

  // Do the multiplication -- one call for every frequency.
  fWatches["DoNoiseMultiplication"].Start(false); // Don't count vector allocation.
#ifdef USE_THREADS
  boost::thread_group threads;

  // frequencies_per_thread is the number of frequencies to distribute, rounding down.
  // threads_with_extra is the number of threads we'll actually give one extra frequency to.
  size_t frequencies_per_thread = (fMaxF-fMinF+1)/(NUM_THREADS); // round down
  size_t threads_with_extra = (fMaxF-fMinF+1) - frequencies_per_thread*(NUM_THREADS);
  assert(frequencies_per_thread*NUM_THREADS <= fMaxF-fMinF+1);
  assert(threads_with_extra < NUM_THREADS);

  for(size_t i = 0; i < threads_with_extra; i++) {
    size_t flo = i*(frequencies_per_thread+1);
    threads.add_thread(new boost::thread(&EXORefitSignals::DoNoiseMultiplication_Range,
                                         this,
                                         flo,
                                         size_t(flo + frequencies_per_thread + 1)));
  }
  for(size_t i = threads_with_extra; i < (NUM_THREADS)-1; i++) {
    size_t flo = threads_with_extra + i*frequencies_per_thread;
    threads.add_thread(new boost::thread(&EXORefitSignals::DoNoiseMultiplication_Range,
                                         this,
                                         flo,
                                         flo + frequencies_per_thread));
  }
  // Don't create the last thread; *this* is the last thread.
  size_t flo = threads_with_extra + ((NUM_THREADS)-1)*frequencies_per_thread;
  DoNoiseMultiplication_Range(flo, (fMaxF-fMinF+1)); // to handle [0, fMaxF-fMinF], inclusive.
  threads.join_all();
#else
  // The sequential version.
  DoNoiseMultiplication_Range(0, (fMaxF-fMinF+1)); // to handle [0, fMaxF-fMinF], inclusive.
#endif
  fWatches["DoNoiseMultiplication"].Stop();

  // Clean up, to be ready for the next call.
  fNoiseMulQueue.clear(); // Hopefully doesn't free memory, since I'll need it again.
  fNumVectorsInQueue = 0;
}

void EXORefitSignals::DoNoiseMultiplication_Range(size_t flo, size_t fhi)
{
  // Do noise multiplications in the range [flo, fhi).
  // In case we are parallelizing the noise multiplication, this lets us
  // break up DoNoiseMultiplication into pieces easily.
  // (If we are not parallelizing, the penalty is probably zero.)
  assert(fNoiseMulQueue.size() == fNoiseColumnLength*fNumVectorsInQueue);
  assert(fNoiseMulResult.size() == fNoiseColumnLength*fNumVectorsInQueue);
  assert(fNoiseCorrelations.size() >= fhi);
  for(size_t f = flo; f < fhi; f++) {
    size_t StartIndex = 2*fChannels.size()*f;
    size_t BlockSize = fChannels.size() * (f < fMaxF - fMinF ? 2 : 1);
    assert(fNoiseCorrelations[f].size() == BlockSize*BlockSize);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                BlockSize, fNumVectorsInQueue, BlockSize,
                1, &fNoiseCorrelations[f][0], BlockSize, &fNoiseMulQueue[StartIndex], fNoiseColumnLength,
                0, &fNoiseMulResult[StartIndex], fNoiseColumnLength);
  }
}

std::vector<double> EXORefitSignals::DoInvLPrecon(std::vector<double>& in, EventHandler& event)
{
  // Multiply by K1_inv.
  std::vector<double> out(in.size(), 0);
  ddiamm(fNoiseColumnLength, event.fWireModel.size()+1, fNoiseColumnLength,
         &fInvSqrtNoiseDiag[0], &in[0], event.fColumnLength,
         &out[0], event.fColumnLength); // out = {{D^(-1/2) v1} {0}}
  mkl_domatcopy('C', 'N', event.fWireModel.size()+1, event.fWireModel.size()+1,
                -1, &in[fNoiseColumnLength], event.fColumnLength,
                &out[fNoiseColumnLength], event.fColumnLength); // out = {{D^(-1/2) v1} {-v2}}
  DoLagrangeAndConstraintMul<'C'>(out, out, event); // out = {{D^(-1/2)v1} {trans(L)D^(-1/2)v1 - v2}}
  mkl_domatcopy('C', 'N', fNoiseColumnLength, event.fWireModel.size()+1,
                1, &in[0], event.fColumnLength,
                &out[0], event.fColumnLength); // out = {{v1} {trans(L)D^(-1/2)v1 - v2}}
  cblas_dtrsm(CblasColMajor, CblasLeft, CblasUpper, CblasTrans, CblasNonUnit,
              event.fWireModel.size()+1, event.fWireModel.size()+1,
              1, &event.fPreconX[0], event.fWireModel.size()+1,
              &out[fNoiseColumnLength], event.fColumnLength);
  // out = {{v1} {Inv(trans(X))(trans(L)D^(-1)v1 - v2)}}.
  return out;
}

std::vector<double> EXORefitSignals::DoInvRPrecon(std::vector<double>& in, EventHandler& event)
{
  // Multiply by K2_inv.
  std::vector<double> out(in.size(), 0);
  mkl_domatcopy('C', 'N', event.fWireModel.size()+1, event.fWireModel.size()+1,
                1, &in[fNoiseColumnLength], event.fColumnLength,
                &out[fNoiseColumnLength], event.fColumnLength); // out = {{0} {v2}}
  cblas_dtrsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit,
              event.fWireModel.size()+1, event.fWireModel.size()+1,
              1, &event.fPreconX[0], event.fWireModel.size()+1,
              &out[fNoiseColumnLength], event.fColumnLength); // out = {{0} {X^(-1)v2}}
  DoLagrangeAndConstraintMul<'L'>(out, out, event); // out = {{LX^(-1)v2} {X^(-1)v2}}

  // Although there is MKL functionality to do the following, it requires a fairly large buffer
  // to perform out-of-place arithmetic.
  // It is probably more effective to just do the looping ourselves in an intelligent manner.
  for(size_t i = 0; i < fNoiseColumnLength; i++) {
    double InvSqrtNoise = fInvSqrtNoiseDiag[i];
    for(size_t j = 0; j < event.fWireModel.size()+1; j++) {
      size_t Index = i + j*event.fColumnLength;
      out[Index] = in[Index] - InvSqrtNoise*out[Index];
    }
  } // out = {{v1 - D^(-1/2)LX^(-1)v2} {X^(-1)v2}}
  return out;
}

std::vector<double> EXORefitSignals::DoLPrecon(std::vector<double>& in, EventHandler& event)
{
  // Multiply by K1.
  std::vector<double> out(in.size(), 0);
  mkl_domatcopy('C', 'N', event.fWireModel.size()+1, event.fWireModel.size()+1,
                1, &in[fNoiseColumnLength], event.fColumnLength,
                &out[fNoiseColumnLength], event.fColumnLength); // out = {{0} {v2}}
  cblas_dtrmm(CblasColMajor, CblasLeft, CblasUpper, CblasTrans, CblasNonUnit,
              event.fWireModel.size()+1, event.fWireModel.size()+1,
              -1, &event.fPreconX[0], event.fWireModel.size()+1,
              &out[fNoiseColumnLength], event.fColumnLength); // out = {{0} {-trans(X)v2}}
  ddiamm(fNoiseColumnLength, event.fWireModel.size()+1, fNoiseColumnLength,
         &fInvSqrtNoiseDiag[0], &in[0], event.fColumnLength,
         &out[0], event.fColumnLength); // out = {{D^(-1/2) v1} {-trans(X)v2}}
  DoLagrangeAndConstraintMul<'C'>(out, out, event); // out = {{D^(-1/2)v1} {trans(L)D^(-1/2)v1 - trans(X)v2}}
  mkl_domatcopy('C', 'N', fNoiseColumnLength, event.fWireModel.size()+1,
                1, &in[0], event.fColumnLength,
                &out[0], event.fColumnLength); // out = {{v1} {trans(L)D^(-1/2)v1 - trans(X)v2}}
  return out;
}

std::vector<double> EXORefitSignals::DoRPrecon(std::vector<double>& in, EventHandler& event)
{
  // Multiply by K2.
  std::vector<double> out(in.size(), 0);
  mkl_domatcopy('C', 'N', event.fWireModel.size()+1, event.fWireModel.size()+1,
                1, &in[fNoiseColumnLength], event.fColumnLength,
                &out[fNoiseColumnLength], event.fColumnLength); // out = {{0} {v2}}
  DoLagrangeAndConstraintMul<'L'>(out, out, event); // out = {{Lv2} {v2}}
  cblas_dtrmm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit,
              event.fWireModel.size()+1, event.fWireModel.size()+1,
              1, &event.fPreconX[0], event.fWireModel.size()+1,
              &out[fNoiseColumnLength], event.fColumnLength); // out = {{Lv2} {Xv2}}

  // Although there is MKL functionality to do the following, it requires a fairly large buffer
  // to perform out-of-place arithmetic.
  // It is probably more effective to just do the looping ourselves in an intelligent manner.
  for(size_t i = 0; i < fNoiseColumnLength; i++) {
    double InvSqrtNoise = fInvSqrtNoiseDiag[i];
    for(size_t j = 0; j < event.fWireModel.size()+1; j++) {
      size_t Index = i + j*event.fColumnLength;
      out[Index] = in[Index] + InvSqrtNoise*out[Index];
    }
  } // out = {{v1 + D^(-1/2)Lv2} {Xv2}}
  return out;
}

size_t EXORefitSignals::RequestNoiseMul(std::vector<double>& vec,
                                        size_t ColLength)
{
  // Request a noise multiplication on vec.
  // Return value indicates where to retrieve results.
  assert(fNoiseColumnLength <= ColLength);
  assert(vec.size() % ColLength == 0);
  size_t NumCols = vec.size() / ColLength;

#ifdef USE_THREADS
  RequestNoiseMulMutex.lock(); // Lock the shared noise multiplication queue for this thread.
#endif

  size_t InitSize = fNoiseMulQueue.size();
  assert(InitSize == fNumVectorsInQueue*fNoiseColumnLength);
  fNoiseMulQueue.reserve(InitSize + NumCols*fNoiseColumnLength);
  for(size_t i = 0; i < NumCols; i++) {
    fNoiseMulQueue.insert(fNoiseMulQueue.end(),
                          vec.begin() + i*ColLength,
                          vec.begin() + i*ColLength + fNoiseColumnLength);
  }
  fNumVectorsInQueue += NumCols;

#ifdef USE_THREADS
  RequestNoiseMulMutex.unlock();
#endif

  return InitSize;
}

void EXORefitSignals::FillFromNoise(std::vector<double>& vec,
                                    size_t NumCols,
                                    size_t ColLength,
                                    size_t ResultIndex)
{
  // Overwrite in with the results from noise multiplication.
  // Leave zeros for rows which are not subject to noise multiplication terms.
  assert(fNoiseColumnLength <= ColLength);
  assert(ResultIndex % fNoiseColumnLength == 0);
  assert(ResultIndex + NumCols*fNoiseColumnLength <= fNoiseMulResult.size());

  vec.assign(NumCols*ColLength, 0);
  for(size_t i = 0; i < NumCols; i++) {
    std::copy(fNoiseMulResult.begin() + ResultIndex +  i   *fNoiseColumnLength,
              fNoiseMulResult.begin() + ResultIndex + (i+1)*fNoiseColumnLength,
              vec.begin() + i*ColLength);
  }
}

EventHandler* EXORefitSignals::PopAnEvent()
{
  // Return an event to process, or NULL if there isn't one.
  // Make sure this is thread-safe.
  EventHandler* evt = NULL;
#ifdef USE_LOCKFREE
  if(not fEventHandlerQueue.pop(evt)) evt = NULL;
#else
#ifdef USE_THREADS
  EventQueueMutex.lock();
#endif
  if(fEventHandlerQueue.empty()) evt = NULL;
  else {
    evt = fEventHandlerQueue.front();
    fEventHandlerQueue.pop();
  }
#ifdef USE_THREADS
  EventQueueMutex.unlock();
#endif
#endif
  return evt;
}

void EXORefitSignals::PushAnEvent(EventHandler* evt)
{
  // Push an event requiring further processing.
  // Make sure this is thread-safe.
#ifdef USE_LOCKFREE
  // We should have already allocated sufficient space.  (Else, it's not really thread-safe.)
  assert(fEventHandlerResults.bounded_push(evt));
#else
#ifdef USE_THREADS
  EventResultsMutex.lock();
#endif
  fEventHandlerResults.push(evt);
#ifdef USE_THREADS
  EventResultsMutex.unlock();
#endif
#endif
}

void EXORefitSignals::HandleEventsInThread()
{
  // Function for a thread to keep grabbing events to handle until no more exist.
  EventHandler* evt = NULL;
  while(evt = PopAnEvent()) {
    // We just drew an event pointer from the queue; handle it.
    bool Result = DoBlBiCGSTAB(*evt);
    if(Result) {
      // This event is done.
      FinishEvent(evt); // Deletes the EventHandler object.
    }
    else {
      // This event is not done, so it must have requested a noise multiplication.
      PushAnEvent(evt);
    }
  }
}

void EXORefitSignals::DoPassThroughEvents()
{
  // Do a pass through the event handlers we've accumulated;
  // do a round of noise multiplication followed by a round of BiCGSTAB.
  DoNoiseMultiplication();
  assert(fEventHandlerResults.empty());
  fWatches["HandleEvents"].Start(false);
#ifdef USE_THREADS
  boost::thread_group threads;
  for(size_t i = 0; i < (NUM_THREADS)-1; i++) {
    threads.add_thread(new boost::thread(&EXORefitSignals::HandleEventsInThread, this));
  }
#endif
  // This thread is the last one.
  // If we're in unthreaded code, this thread is the only one.
  HandleEventsInThread();
#ifdef USE_THREADS
  threads.join_all();
#endif
  fWatches["HandleEvents"].Stop();
  assert(fEventHandlerQueue.empty());

  // Transfer entries from results into queue.
  while(not fEventHandlerResults.empty()) {
#ifdef USE_LOCKFREE
    EventHandler* evt = NULL;
    assert(fEventHandlerResults.unsynchronized_pop(evt));
    assert(fEventHandlerQueue.unsynchronized_push(evt));
#else
    fEventHandlerQueue.push(fEventHandlerResults.front());
    fEventHandlerResults.pop();
#endif
  }
}

bool EXORefitSignals::CanTerminate(EventHandler& event)
{
  // Test whether the residual matrix R indicates we can terminate yet.
  // For now we test for termination against the *unpreconditioned* residual matrix.
  // This should be compared to the alternative of terminating against the preconditioned residual matrix.
  std::vector<double> R_unprec = DoLPrecon(event.fR, event);

  for(size_t col = 0; col <= event.fWireModel.size(); col++) {
    size_t ColIndex = col*event.fColumnLength;
    size_t NextCol = ColIndex + event.fColumnLength;
    double Norm = 0;
    for(size_t i = 0; i < fNoiseColumnLength; i++) {
      Norm += R_unprec[ColIndex + i]*R_unprec[ColIndex+i]*fNoiseDiag[i];
    }
    for(size_t i = fNoiseColumnLength; i < event.fColumnLength; i++) {
      Norm += R_unprec[ColIndex + i]*R_unprec[ColIndex+i];
    }
    if(Norm > fRThreshold*fRThreshold) return false;
  }

  return true;
}
