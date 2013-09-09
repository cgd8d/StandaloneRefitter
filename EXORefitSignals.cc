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
#include "EXOAnalysisManager/EXOGridCorrectionModule.hh"
#include "EXOAnalysisManager/EXOTreeInputModule.hh"
#include "EXOAnalysisManager/EXOTreeOutputModule.hh"
#include "EXOCalibUtilities/EXOChannelMapManager.hh"
#include "EXOCalibUtilities/EXOElectronicsShapers.hh"
#include "EXOCalibUtilities/EXOUWireGains.hh"
#include "EXOCalibUtilities/EXOLifetimeCalib.hh"
#include "EXOCalibUtilities/EXOGridCorrectionCalib.hh"
#include "EXOUtilities/EXONoiseCorrelations.hh"
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
#include <iomanip>
#include <iostream>
#include <sstream>
#include <cstdlib>

EXORefitSignals::EXORefitSignals(EXOTreeInputModule& inputModule,
                                 TTree& wfTree,
                                 EXOTreeOutputModule& outputModule)
: fUseWireAPDCorrelations(true),
  fOnlyEstimateConditionNumbers(false),
  fInputModule(inputModule),
  fWFTree(wfTree),
  fOutputModule(outputModule),
  fLightmapFilename("data/lightmap/LightMaps.root"),
  fRThreshold(0.1),
  fThoriumEnergy_keV(2615),
  fMinF(1),
  fMaxF(1024),
  fNumVectorsInQueue(0)
{
  fWFEvent = NULL;
  fWFTree.SetBranchAddress("EventBranch", &fWFEvent);
  fWatch_TotalTime.Start();
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
  fWatch_FillNoise.Start(false);

  // Get the channel map.
  const EXOChannelMap& ChannelMap = GetChanMapForHeader(ED.fEventHeader);

  // Construct the set of channels to keep.
  std::vector<unsigned char> ChannelsToUse;
  for(unsigned char i = 0; i < NUMBER_READOUT_CHANNELS; i++) {
    if(EXOMiscUtil::TypeOfChannel(i) == EXOMiscUtil::kVWire) continue; // No v wires for now.
    if(ChannelMap.channel_suppressed_by_daq(i) or not ChannelMap.good_channel(i)) continue;
    ChannelsToUse.push_back(i);
  }

  // If the channel mapping is unchanged, do nothing.
  if(ChannelsToUse == fChannels) {
    fWatch_FillNoise.Stop();
    return;
  }

  // Else, we'll need to extract the noise information to match the new ordering.
  // Start by flushing all currently-held events, since the noise information will change.
  FlushEvents();

  // Then fill fNoiseCorrelations.
  // Note that we store the same-frequency blocks in column-major order,
  // to simplify GEMM calls (if BLAS is provided).
  fChannels = ChannelsToUse;
  fNoiseCorrelations.resize(fMaxF - fMinF + 1);
  TFile* NoiseFile = TFile::Open(fNoiseFilename.c_str());
  EXONoiseCorrelations* NoiseCorr =
    dynamic_cast<EXONoiseCorrelations*>(NoiseFile->Get("EXONoiseCorrelations"));
  for(size_t f = fMinF; f <= fMaxF; f++) {
    bool IsFullBlock = (f != fMaxF);
    std::vector<double>& block = fNoiseCorrelations[f-fMinF];
    block.assign(fChannels.size()*fChannels.size() * (IsFullBlock ? 4 : 1), 0);

    // Iterate through column pairs.
    for(size_t index1 = 0; index1 < fChannels.size(); index1++) {
      unsigned char noiseIndex1 = NoiseCorr->GetIndexOfChannel(fChannels[index1]);
      size_t ColPos = index1*fChannels.size()*(IsFullBlock ? 4 : 1);

      // Start with the real column.
      for(size_t index2 = 0; index2 < fChannels.size(); index2++) {
        if(not fUseWireAPDCorrelations) {
          if(EXOMiscUtil::TypeOfChannel(fChannels[index1]) !=
             EXOMiscUtil::TypeOfChannel(fChannels[index2])) continue;
        }
        unsigned char noiseIndex2 = NoiseCorr->GetIndexOfChannel(fChannels[index2]);
        size_t RowPos = ColPos + index2*(IsFullBlock ? 2 : 1);

        // real row.
        block[RowPos] = NoiseCorr->GetRR(f)[noiseIndex2][noiseIndex1];

        // imag row.
        if(IsFullBlock) block[RowPos+1] = NoiseCorr->GetRI(f)[noiseIndex1][noiseIndex2];
      }

      // Now the imag column.
      if(IsFullBlock) {
        ColPos += fChannels.size()*(IsFullBlock ? 2 : 1);
        for(size_t index2 = 0; index2 < fChannels.size(); index2++) {
          unsigned char noiseIndex2 = NoiseCorr->GetIndexOfChannel(fChannels[index2]);
          size_t RowPos = ColPos + index2*(IsFullBlock ? 2 : 1);

          // real row.
          block[RowPos] = NoiseCorr->GetRI(f)[noiseIndex2][noiseIndex1];

          // imag row.
          if(IsFullBlock) block[RowPos+1] = NoiseCorr->GetII(f)[noiseIndex2][noiseIndex1];
        }
      }
    } // End loop over column pairs (index1).
  } // End loop over frequencies.  fNoiseCorrelations is initialized.

  // Cleanup -- this should do it.
  delete NoiseFile;

  // For convenience, pre-store useful parameters.
  for(size_t i = 0; i < fChannels.size(); i++) {
    if(EXOMiscUtil::TypeOfChannel(fChannels[i]) == EXOMiscUtil::kAPDGang) {
      fFirstAPDChannelIndex = i;
      break;
    }
  }
  fNoiseColumnLength = fChannels.size() * (2*(fMaxF-fMinF) + 1);

  fWatch_FillNoise.Stop();
}

int EXORefitSignals::Initialize()
{
  // Open the lightmap file, and extract their information.
  // Create unshaped wire drift waveforms.
  // Also initialize our various timers.

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

  // Initialize stopwatches too.
  fWatch_BiCGSTAB.Reset();
  fWatch_NoiseMul.Reset();
  fWatch_RestMul.Reset();
  fWatch_BiCGSTAB_part1.Reset();
  fWatch_BiCGSTAB_part2.Reset();
  fWatch_FillNoise.Reset();
  fNumEventsHandled = 0;
  fNumSignalsHandled = 0;
  fTotalIterationsDone = 0;
  return 0;
}

EXORefitSignals::~EXORefitSignals()
{
  // Print statistics and timing information.
  fWatch_TotalTime.Stop();

  std::cout<<fNumEventsHandled<<" events were handled by signal refitting."<<std::endl;
  std::cout<<"Those events contained a total of "<<fNumSignalsHandled<<" signals to refit."<<std::endl;
  std::cout<<fTotalIterationsDone<<" iterations were required."<<std::endl;
  std::cout<<"Multiplying by noise blocks:"<<std::endl;
  fWatch_NoiseMul.Print();
  std::cout<<"Total time spent in BiCGSTAB iterations (excluding noise blocks):"<<std::endl;
  fWatch_BiCGSTAB.Print();
  std::cout<<"Multiplying by the rest of the matrix entries:"<<std::endl;
  fWatch_RestMul.Print();
  std::cout<<"Part one of iterations (exploiting V = AP):"<<std::endl;
  fWatch_BiCGSTAB_part1.Print();
  std::cout<<"Part two of iterations (exploiting T = AS):"<<std::endl;
  fWatch_BiCGSTAB_part2.Print();
  std::cout<<"Setting up noise correlations:"<<std::endl;
  fWatch_FillNoise.Print();
  std::cout<<"Total time between construction and destruction:"<<std::endl;
  fWatch_TotalTime.Print();
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

  // For convenience, store the column length we'll be dealing with.
  event->fColumnLength = 2*fChannels.size()*(fMaxF-fMinF) + fChannels.size() + event->fWireModel.size() + 1;

  // Set up a simple, but not quite crazy, initial guess for X.
  event->fX.assign(event->fColumnLength * (event->fWireModel.size() + 1), 0);

  // Start with the wires.
  for(size_t i = 0; i < event->fWireModel.size(); i++) {
    size_t ColIndex = i*event->fColumnLength;
    double Normalization = 0;
    for(std::map<unsigned char, std::vector<double> >::iterator it = event->fWireModel[i].second.begin();
        it != event->fWireModel[i].second.end();
        it++) {
      std::vector<double>& model = it->second;
      size_t channel_index = 0;
      while(fChannels[channel_index] != it->first) {
        channel_index++;
        if(channel_index >= fChannels.size()) LogEXOMsg("Index exceeded -- why can this happen?", EEAlert);
      }
      for(size_t f = fMinF; f <= fMaxF; f++) {
        size_t step = (f < fMaxF ? 2 : 1);
        size_t ColIndex = step*step*fChannels.size()*channel_index;
        double RNoiseVal = fNoiseCorrelations[f-fMinF][ColIndex + step*channel_index];
        Normalization += model[2*(f-fMinF)]*model[2*(f-fMinF)]/RNoiseVal;
        if(step == 2) {
          ColIndex += step*fChannels.size();
          double INoiseVal = fNoiseCorrelations[f-fMinF][ColIndex + step*channel_index + 1];
          Normalization += model[2*(f-fMinF)+1]*model[2*(f-fMinF)+1]/INoiseVal;
        }
      }
    } // Have overall normalization for this signal.
    for(std::map<unsigned char, std::vector<double> >::iterator it = event->fWireModel[i].second.begin();
        it != event->fWireModel[i].second.end();
        it++) {
      std::vector<double>& model = it->second;
      size_t channel_index = 0;
      while(fChannels[channel_index] != it->first) {
        channel_index++;
        if(channel_index >= fChannels.size()) LogEXOMsg("Index exceeded -- why can this happen?", EEAlert);
      }
      for(size_t f = fMinF; f <= fMaxF; f++) {
        size_t step = (f < fMaxF ? 2 : 1);
        size_t RowIndex = ColIndex + 2*fChannels.size()*(f-fMinF);
        RowIndex += channel_index*step;
        size_t NoiseColIndex = step*step*fChannels.size()*channel_index;

        double RNoise = fNoiseCorrelations[f-fMinF][NoiseColIndex + step*channel_index];
        event->fX[RowIndex] = model[2*(f-fMinF)]/(RNoise*Normalization);
        if(f != fMaxF) {
          NoiseColIndex += step*fChannels.size();
          double INoise = fNoiseCorrelations[f-fMinF][NoiseColIndex + step*channel_index + 1];
          event->fX[RowIndex+1] = model[2*(f-fMinF) + 1]/(INoise*Normalization);
        }
      }
    }
  }
  // And then do the one light signal.
  double norm_APDmodel = std::inner_product(event->fmodel_realimag.begin(), event->fmodel_realimag.end(),
                                            event->fmodel_realimag.begin(), double(0));
  double SumSqYieldExpected = 0;
  for(size_t i = fFirstAPDChannelIndex; i < fChannels.size(); i++) {
    SumSqYieldExpected += std::pow(event->fExpectedYieldPerGang[fChannels[i]], 2);
  }
  for(size_t i = fFirstAPDChannelIndex; i < fChannels.size(); i++) {
    double ExpectedYieldOnChannel = event->fExpectedYieldPerGang[fChannels[i]];
    double LeadingFactor = ExpectedYieldOnChannel/(SumSqYieldExpected*norm_APDmodel);

    size_t ColIndex = event->fWireModel.size()*event->fColumnLength;
    for(size_t f = fMinF; f <= fMaxF; f++) {
      size_t RowIndex = ColIndex + 2*fChannels.size()*(f-fMinF);
      RowIndex += i*(f != fMaxF ? 2 : 1);

      event->fX[RowIndex] = LeadingFactor*event->fmodel_realimag[2*(f-fMinF)];
      if(f != fMaxF) event->fX[RowIndex+1] = LeadingFactor*event->fmodel_realimag[2*(f-fMinF)+1];
    }
  }
  // Defer guesses for Lagrange multipliers; we'll use AX to produce these.

  // We can also find the appropriate preconditioner here.
  // This is a pretty good preconditioner, obtained by approximating A ~ {{D L} {trans(L) 0}},
  // where D is diagonal.  Thus, the approximation comes from ignoring noise cross-terms and
  // Poisson noise terms.
  // Haven't decided yet whether it's important to include Poisson terms on the diagonal; currently I don't.
  event->fDiag.resize(fNoiseColumnLength);
  for(size_t i = 0; i < fNoiseColumnLength; i++) {
    size_t f = i/(2*fChannels.size());
    size_t c = i%(2*fChannels.size());
    size_t block = fChannels.size()*(f < fMaxF-fMinF ? 2 : 1);
    event->fDiag[i] = fNoiseCorrelations[f][block*c + c];
  }
  // Find X using trans(X)X = trans(L) D^(-1) L.
  std::vector<double> Temp1(event->fColumnLength * (event->fWireModel.size()+1), 0);
  std::vector<double> Temp2(event->fColumnLength * (event->fWireModel.size()+1), 0);
  for(size_t i = 0; i < event->fWireModel.size()+1; i++) {
    Temp1[i*event->fColumnLength + fNoiseColumnLength + i] = 1;
  } // Temp1 = {{0}{I}}
  DoLagrangeAndConstraintMul<'L'>(Temp1, Temp2, *event); // Temp2 = {{L} {0}}
  for(size_t i = 0; i < Temp2.size(); i++) {
    size_t modi = i % event->fColumnLength;
    if(modi >= fNoiseColumnLength) continue;
    Temp2[i] /= event->fDiag[modi];
  } // Now Temp2 = {{ D^(-1)L} {0}}
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

  // If our objective is to report condition numbers, proceed to do that.
  if(fOnlyEstimateConditionNumbers) {
    std::cout<<"For entry "<<event->fEntryNumber<<":"<<std::endl;
    std::pair<double, double> MinMax = EstimateConditionNumber(*event);
    std::cout<<"\tCondition number >= "<<MinMax.second/MinMax.first<<std::endl;
    for(size_t i = 1; i < 10; i++) {
      std::pair<double, double> TempMinMax = EstimateConditionNumber(*event);
      if(TempMinMax.first < MinMax.first) MinMax.first = TempMinMax.first;
      if(TempMinMax.second > MinMax.second) MinMax.second = TempMinMax.second;
      std::cout<<"\tCondition number >= "<<MinMax.second/MinMax.first<<std::endl;
    }
    delete event;
    return;
  }

  // Push event onto the list of event handlers.
  fEventHandlerList.push_back(event);

  // Request a matrix multiplication of X.
  event->fResultIndex = RequestNoiseMul(event->fX, event->fColumnLength);
  fNumEventsHandled++; // One more event that will be actually handled.
  fNumSignalsHandled += event->fWireModel.size() + 1;

  // Now, while there are enough requests in the queue, satisfy those requests.
  while(fNumVectorsInQueue > 300) {
    DoNoiseMultiplication();
    std::list<EventHandler*>::iterator it = fEventHandlerList.begin();
    while(it != fEventHandlerList.end()) {
      bool Result = DoBlBiCGSTAB(**it);
      if(Result) {
        // This event is done.
        FinishEvent(*it); // Deletes the EventHandler object.
        it = fEventHandlerList.erase(it); // Removes it from the list.
      }
      else {
        it++;
      }
    }
    if(fNumVectorsInQueue == 0 xor fEventHandlerList.empty()) { // Sanity check.
      LogEXOMsg("Inconsistent state between fNumVectorsInQueue and fEventHandlerList", EEAlert);
    }
  }
}

void EXORefitSignals::FlushEvents()
{
  // Finish processing for all events in the event handler list,
  // regardless of how many pending multiplication requests are queued.

  while(not fEventHandlerList.empty()) {
    DoNoiseMultiplication();
    std::list<EventHandler*>::iterator it = fEventHandlerList.begin();
    while(it != fEventHandlerList.end()) {
      bool Result = DoBlBiCGSTAB(**it);
      if(Result) {
        // This event is done.
        FinishEvent(*it); // Deletes the EventHandler object.
        it = fEventHandlerList.erase(it); // Removes it from the list.
      }
      else {
        it++;
      }
    }
    if(fNumVectorsInQueue == 0 xor fEventHandlerList.empty()) { // Sanity check.
      LogEXOMsg("Inconsistent state between fNumVectorsInQueue and fEventHandlerList", EEAlert);
    }
  }
}

void EXORefitSignals::FinishEvent(EventHandler* event)
{
  // Compute and fill denoised signals, as appropriate.
  // Then pass the filled event to the output module.
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
          size_t XIndex = event->fColumnLength*i + 2*fChannels.size()*f + chan_index*(f < fMaxF-fMinF ? 2 : 1);
          Results[i] += event->fX[XIndex]*WF_real[chan_index][f+fMinF];
          if(f < fMaxF-fMinF) Results[i] += event->fX[XIndex+1]*WF_imag[chan_index][f+fMinF];
        }
      }
    }

    // Translate signal magnitudes into corresponding objects.
    const EXOUWireGains* GainsFromDatabase = GetCalibrationFor(EXOUWireGains,
                                                               EXOUWireGainsHandler,
                                                               "source_calibration",
                                                               ED->fEventHeader);
    for(size_t i = 0; i < event->fWireModel.size(); i++) {
      size_t sigIndex = event->fWireModel[i].first;
      EXOUWireSignal* sig = ED->GetUWireSignal(sigIndex);
      double UWireScalingFactor = ADC_FULL_SCALE_ELECTRONS_WIRE * W_VALUE_LXE_EV_PER_ELECTRON /
                                  (CLHEP::keV * ADC_BITS);
      double GainCorrection = GainsFromDatabase->GetGainOnChannel(sig->fChannel)/300.0;
      sig->fDenoisedEnergy = Results[i]*GainCorrection*UWireScalingFactor;
    }
    ED->GetScintillationCluster(0)->fDenoisedEnergy = Results.back()*fThoriumEnergy_keV;
  } // End setting of denoised energy signals.

  fOutputModule.ProcessEvent(ED);
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
  fWatch_BiCGSTAB.Start(false);
  lapack_int ret;

  if(event.fR.size() == 0) {
    // We're still in the setup phase.
    // Start by copying result into R.
    FillFromNoise(event.fR, event.fWireModel.size()+1, event.fColumnLength, event.fResultIndex);
    // But wait, there's more: use this partial AX to guess lagrange multipliers.
    std::vector<double> Y = event.fR;
    DoRestOfMultiplication(event.fX, Y, event); // To make sure we factor in poisson terms.
    std::vector<double> B(event.fColumnLength * (event.fWireModel.size()+1), 0);
    for(size_t i = 0; i <= event.fWireModel.size(); i++) {
      size_t Index = (i+1)*event.fColumnLength; // Next column; then subtract.
      Index -= event.fWireModel.size() + 1; // Step backward.
      Index += i; // Go forward to the right entry.
      B[Index] = 1; // All models are normalized to 1.
    }
    // Use B to extract the Lagrange multiplier terms of the matrix.
    std::vector<double> Lterms(event.fColumnLength * (event.fWireModel.size()+1), 0);
    DoRestOfMultiplication(B, Lterms, event); // A little wasteful, but that's OK for now.
    // Solve for a left-inverse of L.
    // The idea here is, if A = { {M L} {C 0} } (blocks) and our initial guess is {{X}{K}},
    // we want MX + LK ~ 0.  Now we've computed MX; so we can solve for a value of K
    // which makes this smallish.
    // First, compute L_trans L.
    std::vector<double> Ltrans_L((event.fWireModel.size()+1)*(event.fWireModel.size()+1),0);
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                event.fWireModel.size()+1, event.fWireModel.size()+1, fNoiseColumnLength,
                1, &Lterms[0], event.fColumnLength, &Lterms[0], event.fColumnLength,
                0, &Ltrans_L[0], event.fWireModel.size()+1);
    // Now find the inverse of Ltrans_L, and put it back into Ltrans_L.
    lapack_int* ipiv = new lapack_int[event.fWireModel.size()+1];
    ret = LAPACKE_dgetrf(LAPACK_COL_MAJOR, event.fWireModel.size()+1, event.fWireModel.size()+1,
                         &Ltrans_L[0], event.fWireModel.size()+1,
                         ipiv);
    if(ret != 0) LogEXOMsg("Failed Factorization", EEAlert);
    ret = LAPACKE_dgetri(LAPACK_COL_MAJOR, event.fWireModel.size()+1,
                         &Ltrans_L[0], event.fWireModel.size()+1,
                         ipiv);
    if(ret != 0) LogEXOMsg("Failed Inversion", EEAlert);
    delete ipiv;
    // Multiply on the right by transpose-L; put the left inverse of L in LInv.
    std::vector<double> LInv(event.fColumnLength * (event.fWireModel.size()+1), 0);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                event.fWireModel.size()+1, event.fColumnLength, event.fWireModel.size()+1,
                1, &Ltrans_L[0], event.fWireModel.size()+1, &Lterms[0], event.fColumnLength,
                0, &LInv[0], event.fWireModel.size()+1);
    // Add an approximation of the lagrange terms to X; then we can finish multiplying R <-- AX.
    // The terms being modified would not interact with the noise anyway, so it's OK.
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                event.fWireModel.size()+1, event.fWireModel.size()+1, fNoiseColumnLength,
                -1, &LInv[0], event.fWireModel.size()+1, &Y[0], event.fColumnLength,
                0, &event.fX[fNoiseColumnLength], event.fColumnLength);
    // Now need to finish multiplying by A, accounting for the other terms.
    DoRestOfMultiplication(event.fX, event.fR, event);
    // Now, R <-- B - R = B - AX.
    for(size_t i = 0; i < event.fR.size(); i++) {
      event.fR[i] = B[i] - event.fR[i];
    }
    // Now precondition R and X0 appropriately.
    event.fX = DoRPrecon(event.fX, event);
    event.fR = DoInvLPrecon(event.fR, event);
    // Set up other pieces of the handler.
    event.fP = event.fR;
    event.fR0hat = event.fR;
    // We'll take K1_inv A K2_inv P, so compute K2_inv P.
    event.fprecon_tmp = DoInvRPrecon(event.fP, event);
    // Now we want V <- AP, so request a multiplication by P.
    event.fResultIndex = RequestNoiseMul(event.fprecon_tmp, event.fColumnLength);
    fWatch_BiCGSTAB.Stop();
    return false;
  }
  else if(event.fV.size() == 0) {
    // At the beginning of the iteration, we just computed V = AP.
    fTotalIterationsDone++;
    fWatch_BiCGSTAB_part1.Start(false);
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
    fWatch_BiCGSTAB_part1.Stop();
    fWatch_BiCGSTAB.Stop();
    return false;
  }
  else {
    // We're in the second half of the iteration, where T was just computed.
    fWatch_BiCGSTAB_part2.Start(false);
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
    // Check if we should conclude here -- R is the residual matrix.
    // Note: for the benefit of preconditioning comparison studies, I will currently terminate based on
    // the unpreconditioned residual.  That may be the wrong thing to do down the road.
    double WorstNorm = 0;
    std::vector<double> R_unprec = DoLPrecon(event.fR, event);
    for(size_t col = 0; col <= event.fWireModel.size(); col++) {
      size_t ColIndex = col*event.fColumnLength;
      size_t NextCol = ColIndex + event.fColumnLength;
      double Norm = std::inner_product(R_unprec.begin() + ColIndex,
                                       R_unprec.begin() + NextCol,
                                       R_unprec.begin() + ColIndex,
                                       double(0));
      std::cout<<"Column "<<col<<" has norm "<<Norm<<std::endl;
      if(Norm > WorstNorm) WorstNorm = Norm;
    }
    if(WorstNorm < fRThreshold*fRThreshold) {
      fWatch_BiCGSTAB_part2.Stop();
      fWatch_BiCGSTAB.Stop();
      return true;
    }
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
    fWatch_BiCGSTAB_part2.Stop();
    fWatch_BiCGSTAB.Stop();
    return false;
  }
}

void EXORefitSignals::DoRestOfMultiplication(const std::vector<double>& in,
                                             std::vector<double>& out,
                                             EventHandler& event)
{
  // After noise terms have already been handled, deal with all of the others.
  // This should not be the bottleneck.
  fWatch_RestMul.Start(false);
  DoPoissonMultiplication(in, out, event);
  DoLagrangeAndConstraintMul<'A'>(in, out, event);
  fWatch_RestMul.Stop();
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
        size_t InIndex = n*event.fColumnLength + 2*fChannels.size()*g + k*(g < fMaxF-fMinF ? 2 : 1);
        CommonFactor += event.fmodel_realimag[2*g] * in[InIndex];
        if(g < fMaxF-fMinF) CommonFactor += event.fmodel_realimag[2*g+1] * in[InIndex+1];
      }
      CommonFactor *= ChannelFactors;

      // Now actually transfer the changes to the out vector.
      for(size_t f = 0; f <= fMaxF - fMinF; f++) {
        size_t OutIndex = n*event.fColumnLength + 2*fChannels.size()*f + k*(f < fMaxF-fMinF ? 2 : 1);
        out[OutIndex] += CommonFactor * event.fmodel_realimag[2*f];
        if(f < fMaxF-fMinF) out[OutIndex+1] += CommonFactor * event.fmodel_realimag[2*f+1];
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
  fWatch_NoiseMul.Start(false); // Don't count vector allocation.
  for(size_t f = 0; f <= fMaxF - fMinF; f++) {
    size_t StartIndex = 2*fChannels.size()*f;
    size_t BlockSize = fChannels.size() * (f < fMaxF - fMinF ? 2 : 1);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                BlockSize, fNumVectorsInQueue, BlockSize,
                1, &fNoiseCorrelations[f][0], BlockSize, &fNoiseMulQueue[StartIndex], fNoiseColumnLength,
                0, &fNoiseMulResult[StartIndex], fNoiseColumnLength);
  }
  fWatch_NoiseMul.Stop();

  // Clean up, to be ready for the next call.
  fNoiseMulQueue.clear(); // Hopefully doesn't free memory, since I'll need it again.
  fNumVectorsInQueue = 0;
}

std::vector<double> EXORefitSignals::DoInvLPrecon(std::vector<double>& in, EventHandler& event)
{
  // Multiply by K1_inv.
  std::vector<double> out(in.size(), 0);
  for(size_t i = 0; i < out.size(); i++) {
    size_t imod = i % event.fColumnLength;
    if(imod < fNoiseColumnLength) out[i] = in[i] / event.fDiag[imod];
  } // out = {{D^(-1) v1} {0}}
  DoLagrangeAndConstraintMul<'C'>(out, out, event); // out = {{D^(-1)v1} {trans(L)D^(-1)v1}}
  for(size_t i = 0; i < in.size(); i++) {
    size_t imod = i % event.fColumnLength;
    if(imod < fNoiseColumnLength) out[i] = in[i]/std::sqrt(event.fDiag[imod]);
    else out[i] -= in[i];
  } // out = {{D^(-1/2)v1} {trans(L)D^(-1)v1 - v2}}
  cblas_dtrsm(CblasColMajor, CblasLeft, CblasUpper, CblasTrans, CblasNonUnit,
              event.fWireModel.size()+1, event.fWireModel.size()+1,
              1, &event.fPreconX[0], event.fWireModel.size()+1,
              &out[fNoiseColumnLength], event.fColumnLength);
  // out = {{D^(-1/2)v1} {Inv(trans(X))(trans(L)D^(-1)v1 - v2)}}.
  return out;
}

std::vector<double> EXORefitSignals::DoInvRPrecon(std::vector<double>& in, EventHandler& event)
{
  // Multiply by K2_inv.
  std::vector<double> out(in.size(), 0);
  for(size_t i = 0; i < in.size(); i++) {
    size_t imod = i % event.fColumnLength;
    if(imod >= fNoiseColumnLength) out[i] = in[i];
  }
  cblas_dtrsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit,
              event.fWireModel.size()+1, event.fWireModel.size()+1,
              1, &event.fPreconX[0], event.fWireModel.size()+1,
              &out[fNoiseColumnLength], event.fColumnLength); // out = {{0} {X^(-1)v2}}
  DoLagrangeAndConstraintMul<'L'>(out, out, event); // out = {{LX^(-1)v2} {X^(-1)v2}}
  for(size_t i = 0; i < out.size(); i++) {
    size_t imod = i % event.fColumnLength;
    if(imod < fNoiseColumnLength) out[i] = in[i]/sqrt(event.fDiag[imod]) - out[i]/event.fDiag[imod];
  } // out = {{D^(-1/2)v1 - D^(-1)LX^(-1)v2} {X^(-1)v2}}
  return out;
}

std::vector<double> EXORefitSignals::DoLPrecon(std::vector<double>& in, EventHandler& event)
{
  // Multiply by K1.
  std::vector<double> out = in;
  cblas_dtrmm(CblasColMajor, CblasLeft, CblasUpper, CblasTrans, CblasNonUnit,
              event.fWireModel.size()+1, event.fWireModel.size()+1,
              -1, &event.fPreconX[0], event.fWireModel.size()+1,
              &out[fNoiseColumnLength], event.fColumnLength); // out = {{v1} {-trans(X)v2}}
  for(size_t i = 0; i < out.size(); i++) {
    size_t imod = i % event.fColumnLength;
    if(imod < fNoiseColumnLength) out[i] /= std::sqrt(event.fDiag[imod]);
  } // out = {{D^(-1/2)v1} {-trans(X)v2}}
  DoLagrangeAndConstraintMul<'C'>(out, out, event); // out = {{D^(-1/2)v1} {trans(L)D^(-1/2)v1 - trans(X)v2}}
  for(size_t i = 0; i < out.size(); i++) {
    size_t imod = i % event.fColumnLength;
    if(imod < fNoiseColumnLength) out[i] = std::sqrt(event.fDiag[imod])*in[i];
  } // out = {{D^(1/2)v1} {trans(L)D^(-1/2)v1 - trans(X)v2}}
  return out;
}

std::vector<double> EXORefitSignals::DoRPrecon(std::vector<double>& in, EventHandler& event)
{
  // Multiply by K2.
  std::vector<double> out(in.size(), 0);
  for(size_t i = 0; i < out.size(); i++) {
    size_t imod = i % event.fColumnLength;
    if(imod >= fNoiseColumnLength) out[i] = in[i];
  } // out = {{0} {v2}}  
  DoLagrangeAndConstraintMul<'L'>(out, out, event); // out = {{Lv2} {v2}}
  cblas_dtrmm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit,
              event.fWireModel.size()+1, event.fWireModel.size()+1,
              1, &event.fPreconX[0], event.fWireModel.size()+1,
              &out[fNoiseColumnLength], event.fColumnLength); // out = {{Lv2} {Xv2}}
  for(size_t i = 0; i < out.size(); i++) {
    size_t imod = i % event.fColumnLength;
    if(imod < fNoiseColumnLength) out[i] = std::sqrt(event.fDiag[imod])*in[i] +
                                           out[i]/std::sqrt(event.fDiag[imod]);
  } // out = {{D^(1/2)v1 + D^(-1/2)Lv2} {Xv2}}
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
  size_t InitSize = fNoiseMulQueue.size();
  fNoiseMulQueue.reserve(InitSize + NumCols*fNoiseColumnLength);
  for(size_t i = 0; i < NumCols; i++) {
    fNoiseMulQueue.insert(fNoiseMulQueue.end(),
                          vec.begin() + i*ColLength,
                          vec.begin() + i*ColLength + fNoiseColumnLength);
  }
  fNumVectorsInQueue += NumCols;
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
    std::copy(fNoiseMulResult.begin() +  i   *fNoiseColumnLength,
              fNoiseMulResult.begin() + (i+1)*fNoiseColumnLength,
              vec.begin() + i*ColLength);
  }
}

std::pair<double, double> EXORefitSignals::EstimateConditionNumber(EventHandler& event)
{
  // Stupid algorithm for finding the condition number -- just keep guessing vectors,
  // and tracking the maximum and minimum value of |Ax|/|x|.
  // The condition number is defined as max/min of that expression.
  // The point here is to avoid solving Ax=b, since that's the whole challenge to begin with.
  // Better might be to do a gradient-based search, instead of randomized brute force,
  // but at the moment I don't know how many false optima there are in |Ax|/|x|.
  assert(fNumVectorsInQueue == 0);
  const size_t NumColumns = 100;

  std::vector<double> RandomX;
  RandomX.reserve(event.fColumnLength*NumColumns); // Avoid value initialization.
  for(size_t i = 0; i < event.fColumnLength*NumColumns; i++) RandomX.push_back(double(rand())/RAND_MAX - 0.5);

  std::vector<double> PreconRandomX = DoInvRPrecon(RandomX, event);
  size_t RequestIndex = RequestNoiseMul(PreconRandomX, event.fColumnLength);
  DoNoiseMultiplication();
  std::vector<double> UnpreconRandomAX;
  FillFromNoise(UnpreconRandomAX, NumColumns, event.fColumnLength, RequestIndex);
  DoRestOfMultiplication(PreconRandomX, UnpreconRandomAX, event);
  std::vector<double> PreconRandomAX = DoInvLPrecon(UnpreconRandomAX, event);

  double Max, Min;
  for(size_t i = 0; i < NumColumns; i++) {
    double Norm_AX = std::inner_product(PreconRandomAX.begin() + i*event.fColumnLength,
                                        PreconRandomAX.begin() + (i+1)*event.fColumnLength,
                                        PreconRandomAX.begin() + i*event.fColumnLength,
                                        double(0));
    double Norm_X = std::inner_product(RandomX.begin() + i*event.fColumnLength,
                                       RandomX.begin() + (i+1)*event.fColumnLength,
                                       RandomX.begin() + i*event.fColumnLength,
                                       double(0));
    double Ratio = Norm_AX/Norm_X;
    if(i == 0) Max = Min = Ratio;
    else {
      if(Ratio > Max) Max = Ratio;
      if(Ratio < Min) Min = Ratio;
    }
  }
  return std::make_pair(Min, Max);
}
