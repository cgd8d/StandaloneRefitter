/*
For whatever reason, the noise I'm using isn't correct.
To try to clarify this, write one program which outputs noise correlations in exactly the format I ultimately need.  It's not as nice as Mike's stuff, but I think I need to give up on using ROOT files due to the memory concerns with reading them.

Program format: the first argument should be the output file;
the second should give the (maximum) number of events to draw from each run;
all subsequent arguments should be run numbers to use.
So, eg. we could call:
./MakeNoiseCorrelationObject OutFile.dat 1000 4001 4002 4003 ...
For simplicity in the interface, I assume every LB run has exactly one file segment.

Be sure to compile with optimization enabled!  It makes a big difference.
*/

#include "EXOUtilities/EXOEventData.hh"
#include "EXOUtilities/EXOEventHeader.hh"
#include "EXOUtilities/EXOCoincidences.hh"
#include "EXOUtilities/EXOTemplWaveform.hh"
#include "EXOUtilities/EXOWaveformFT.hh"
#include "EXOUtilities/EXOWaveform.hh"
#include "EXOUtilities/EXOFastFourierTransformFFTW.hh"
#include "TFile.h"
#include "TTree.h"
#include <cassert>
#include <cstdlib>
#include <sstream>
#include <iomanip>
#include <fstream>

EXOCoincidences coinc;

bool IsEventAcceptable(const EXOEventData* event)
{
  // Check that this event is an acceptable example of pure electronic noise.
  const EXOEventHeader& header = event->fEventHeader;

  // Only use solicited triggers.
  if(header.fIndividualTriggerRequest or header.fSumTriggerRequest) return false;

  // Require full-length events.
  if(header.fSampleCount != 2047) return false; // 0 = 1 sample; 2047 = 2048 samples.

  // Don't use events with abnormal tags.
  if(event->fHasSaturatedChannel or header.fTaggedAsNoise or header.fSirenActiveInCR) return false;

  // Reconstruction has the lowest threshold available -- if it found a signal, skip the event.
  // Ignore v-wire signals.
  if(event->GetNumAPDSignals() or event->GetNumUWireSignals()) return false;

  // Skip events which are vetoed by any sort of coincidence.
  if(coinc.IsVetoed(*event)) return false;

  // OK :)
  return true;
}

void MakeProductWFs(const EXOWaveformFT& wf1, const EXOWaveformFT& wf2,
                    EXODoubleWaveform& rrOut, EXODoubleWaveform& iiOut,
                    EXODoubleWaveform& riOut, EXODoubleWaveform& irOut)
{
  const size_t Length = wf1.GetLength();
  assert(wf2.GetLength() == Length and
         rrOut.GetLength() == Length and
         iiOut.GetLength() == Length and
         riOut.GetLength() == Length and
         irOut.GetLength() == Length);
  for(size_t i = 0; i < Length; i++) {
    rrOut[i] += wf1[i].real()*wf2[i].real();
    iiOut[i] += wf1[i].imag()*wf2[i].imag();
    riOut[i] += wf1[i].real()*wf2[i].imag();
    irOut[i] += wf1[i].imag()*wf2[i].real();
  }
}

void MakeProductWFs(const EXOWaveformFT& wfIn,
                    EXODoubleWaveform& rrOut, EXODoubleWaveform& iiOut,
                    EXODoubleWaveform& riOut)
{
  // Product of waveform with itself.
  // Faster, but also avoids overlapping outputs.
  const size_t Length = wfIn.GetLength();
  assert(rrOut.GetLength() == Length and
         iiOut.GetLength() == Length and
         riOut.GetLength() == Length);
  for(size_t i = 0; i < Length; i++) {
    rrOut[i] += wfIn[i].real()*wfIn[i].real();
    iiOut[i] += wfIn[i].imag()*wfIn[i].imag();
    riOut[i] += wfIn[i].real()*wfIn[i].imag();
  }
}

int main(int argc, char** argv)
{
  assert(sizeof(double) == 8);
  assert(argc > 3);
  const int EntriesPerRun = std::atoi(argv[2]);

  std::vector<int> Runs;
  for(int i = 3; i < argc; i++) Runs.push_back(std::atoi(argv[i]));

  // Number of channels we're going to be dealing with.
  size_t ChannelsToHold = 2*NCHANNEL_PER_WIREPLANE + 2*NUMBER_APD_CHANNELS_PER_PLANE;

  // Hold the pairwise products of FFT waveforms.
  std::map<std::pair<size_t, size_t>, EXODoubleWaveform> RRProducts;
  std::map<std::pair<size_t, size_t>, EXODoubleWaveform> IIProducts;
  std::map<std::pair<size_t, size_t>, EXODoubleWaveform> RIProducts;
  for(size_t i = 0; i < ChannelsToHold; i++) {
    for(size_t j = 0; j < ChannelsToHold; j++) {
      // RI -- not symmetric in ij, so always do it.
      RIProducts[std::make_pair(i, j)].SetLength(1025);
      RIProducts[std::make_pair(i, j)].Zero();

      // The others are symmetric.
      if(i <= j) {
        RRProducts[std::make_pair(i, j)].SetLength(1025);
        RRProducts[std::make_pair(i, j)].Zero();
        IIProducts[std::make_pair(i, j)].SetLength(1025);
        IIProducts[std::make_pair(i, j)].Zero();
      }
    }
  }

  size_t NumEntriesAccepted = 0;

  for(size_t runIndex = 0; runIndex < Runs.size(); runIndex++) {
    std::cout<<"Starting on run "<<Runs[runIndex]<<std::endl;

    std::ostringstream RawRunName;
    RawRunName << "root://exo-rdr.slac.stanford.edu//exo_data/data/WIPP/root/";
    RawRunName << Runs[runIndex] << "/run";
    RawRunName << std::setw(8) << std::setfill('0') << Runs[runIndex];
    RawRunName << "-000.root";

    std::ostringstream ProcRunName;
    ProcRunName << "/nfs/slac/g/exo_data3/exo_data/data/WIPP/masked/";
    ProcRunName << Runs[runIndex] << "/masked";
    ProcRunName << std::setw(8) << std::setfill('0') << Runs[runIndex];
    ProcRunName << "-000.root";

    TFile* ProcFile = TFile::Open(ProcRunName.str().c_str());
    TTree* ProcTree = (TTree*)ProcFile->Get("tree");
    EXOEventData* ProcEvent = NULL;
    ProcTree->SetBranchAddress("EventBranch", &ProcEvent);

    TFile* RawFile = TFile::Open(RawRunName.str().c_str());
    TTree* RawTree = (TTree*)RawFile->Get("tree");
    EXOEventData* RawEvent = NULL;
    RawTree->SetBranchAddress("EventBranch", &RawEvent);
    assert(RawTree->BuildIndex("fRunNumber", "fEventNumber") >= 0);

    coinc.Load(ProcRunName.str()); // Clears any previously-loaded file.

    int NumAcceptedFromThisRun = 0;
    Long64_t EntryNum = 0;
    while(NumAcceptedFromThisRun < EntriesPerRun and EntryNum < ProcTree->GetEntries()) {
      ProcTree->GetEntry(EntryNum);
      EntryNum++;
      if(not IsEventAcceptable(ProcEvent)) continue;
      std::cout<<"\tAccepted event "<<ProcEvent->fEventNumber<<std::endl;

      // Get the raw entry by index, since entry numbers won't generally match (due to masking).
      assert(RawTree->GetEntryWithIndex(ProcEvent->fRunNumber, ProcEvent->fEventNumber));

      // Convert EXOIntWaveforms to EXODoubleWaveforms.
      // Also establish the channel ordering.
      std::vector<EXOWaveformFT> FourierWaveforms;
      for(Int_t channel = 0; channel < NUMBER_READOUT_CHANNELS; channel++) {
        // For now, we get u-wires and APDs, but not v-wires.
        if(EXOMiscUtil::TypeOfChannel(channel) == EXOMiscUtil::kVWire) continue;
        EXOWaveform* rawWF = RawEvent->GetWaveformData()->GetWaveformWithChannelToEdit(channel);
        if(rawWF) {
          rawWF->Decompress();
          EXODoubleWaveform dblWF = *rawWF;
          EXOWaveformFT ftWF;
          assert(dblWF.GetLength() == 2048);
          EXOFastFourierTransformFFTW::GetFFT(2048).PerformFFT(dblWF, ftWF);
          assert(ftWF.GetLength() == 1025); // Just to make sure I'm not confused.
          FourierWaveforms.push_back(ftWF);
        }
        else {
          // This channel doesn't exist in this run -- but don't leave an empty space,
          // since that introduces the nuisance of including a mapping from index to software channel.
          // Instead, just fill with zeros.
          EXOWaveformFT tempFT;
          tempFT.SetLength(1025);
          tempFT.Zero();
          FourierWaveforms.push_back(tempFT);
        }
      }
      assert(FourierWaveforms.size() == ChannelsToHold);

      // Take products of the vectors, in the format most convenient for now.
      for(size_t i = 0; i < ChannelsToHold; i++) {
        // Along the diagonal.
        MakeProductWFs(FourierWaveforms[i],
                       RRProducts[std::make_pair(i, i)],
                       IIProducts[std::make_pair(i, i)],
                       RIProducts[std::make_pair(i, i)]);
        // And off-diagonal.
        for(size_t j = i+1; j < ChannelsToHold; j++) {
          MakeProductWFs(FourierWaveforms[i], FourierWaveforms[j],
                         RRProducts[std::make_pair(i, j)],
                         IIProducts[std::make_pair(i, j)],
                         RIProducts[std::make_pair(i, j)],
                         RIProducts[std::make_pair(j, i)]);
        }
      }

      NumAcceptedFromThisRun++;
      NumEntriesAccepted++;
    } // Finish loop over entries in a particular run.
    delete RawFile;
    delete ProcFile;
  } // Finish loop over runs.

  // RRProducts[i, j][f] = < Re(X_i[f]) Re(X_j[f]) >; only entries with i <= j are stored.
  // IIProducts[i, j][f] = < Im(X_i[f]) Im(X_j[f]) >; only entries with i <= j are stored.
  // RIProducts[i, j][f] = < Re(X_i[f]) Im(X_j[f]) >.

  // Initialize the very large structure to hold rearranged noise correlations.
  // We'll drop the baseline component, leaving only 1024 frequency components.
  std::cout<<"Allocate NoiseCorrelations."<<std::endl;
  std::vector<std::vector<double> > NoiseCorrelations;
  NoiseCorrelations.resize(1024);
  for(size_t i = 0; i < NoiseCorrelations.size(); i++) {
    size_t Multiplier = (i == 1023 ? 1 : 4);
    NoiseCorrelations[i].assign(Multiplier*ChannelsToHold*ChannelsToHold, 0);
  }
  std::cout<<"Done allocating NoiseCorrelations."<<std::endl;

  // RR entries.
  std::cout<<"Fill RR entries."<<std::endl;
  for(size_t i = 0; i < ChannelsToHold; i++) {
    for(size_t j = i; j < ChannelsToHold; j++) {
      const EXODoubleWaveform& wf = RRProducts[std::make_pair(i, j)];
      for(size_t f = 1; f < 1024; f++) {
        NoiseCorrelations[f-1][2*ChannelsToHold*j + i] = wf[f]/NumEntriesAccepted;
        NoiseCorrelations[f-1][2*ChannelsToHold*i + j] = wf[f]/NumEntriesAccepted;
      }
      // f=1024 is differently sized, so handle specially.
      NoiseCorrelations[1023][ChannelsToHold*j + i] = wf[1024]/NumEntriesAccepted;
      NoiseCorrelations[1023][ChannelsToHold*i + j] = wf[1024]/NumEntriesAccepted;
    }
  }
  std::cout<<"Done filling RR entries."<<std::endl;

  // II entries.
  std::cout<<"Fill II entries."<<std::endl;
  for(size_t i = 0; i < ChannelsToHold; i++) {
    for(size_t j = i; j < ChannelsToHold; j++) {
      const EXODoubleWaveform& wf = IIProducts[std::make_pair(i, j)];
      for(size_t f = 1; f < 1024; f++) {
        NoiseCorrelations[f-1][2*ChannelsToHold*(ChannelsToHold+j) + ChannelsToHold + i] =
          wf[f]/NumEntriesAccepted;
        NoiseCorrelations[f-1][2*ChannelsToHold*(ChannelsToHold+i) + ChannelsToHold + j] =
          wf[f]/NumEntriesAccepted;
      }
      // For f=1024, imaginary entries are identically zero.
    }
  }
  std::cout<<"Done filling II entries."<<std::endl;

  // RI entries.
  std::cout<<"Fill RI entries."<<std::endl;
  for(size_t i = 0; i < ChannelsToHold; i++) {
    for(size_t j = 0; j < ChannelsToHold; j++) {
      const EXODoubleWaveform& wf = RIProducts[std::make_pair(i, j)];
      for(size_t f = 1; f < 1024; f++) {
        NoiseCorrelations[f-1][2*ChannelsToHold*i + ChannelsToHold + j] = wf[f]/NumEntriesAccepted;
        NoiseCorrelations[f-1][2*ChannelsToHold*(ChannelsToHold+j) + i] = wf[f]/NumEntriesAccepted;
      }
    }
  }
  std::cout<<"Done filling RI entries."<<std::endl;

  // Write out to file.
  std::filebuf outfile;
  outfile.open(argv[1],
               std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
  assert(outfile.pubseekoff(0, std::ios_base::cur, std::ios_base::out) == std::streampos(0));
  for(size_t f = 0; f < NoiseCorrelations.size(); f++) {
    size_t numChars = sizeof(double)*NoiseCorrelations[f].size();
    assert(outfile.sputn((char*)&NoiseCorrelations[f][0], numChars) == numChars);
  }
  std::cout<<"Done; closing file as we terminate."<<std::endl;
}
