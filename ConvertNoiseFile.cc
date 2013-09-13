#include "EXOUtilities/EXONoiseCorrelations.hh"
#include "EXOUtilities/EXOMiscUtil.hh"
#include "EXOUtilities/EXODimensions.hh"
#include "TFile.h"
#include <cassert>
#include <fstream>

int main()
{
  assert(sizeof(double) == 8);
  TFile* infile = TFile::Open("/global/u1/c/claytond/TestEXOAnalysisBuild/noise_manyruns_withuwires_100000.root");
  EXONoiseCorrelations* noise = dynamic_cast<EXONoiseCorrelations*>(infile->Get("EXONoiseCorrelations"));

  std::vector<unsigned char> ChannelsToUse;
  for(unsigned char i = 0; i < NUMBER_READOUT_CHANNELS; i++) {
    if(EXOMiscUtil::TypeOfChannel(i) == EXOMiscUtil::kVWire) continue; // No v wires for now.
    assert(noise->HasChannel(i));
    ChannelsToUse.push_back(i);
  }
  assert(ChannelsToUse.size() == NUMBER_READOUT_CHANNELS - 2*NCHANNEL_PER_WIREPLANE);

  std::filebuf outfile;
  outfile.open("/global/u1/c/claytond/TestEXOAnalysisBuild/noise_manyruns_withuwires_100000.dat",
               std::ios_base::out | std::ios_base::binary);
  assert(outfile.pubseekoff(0, std::ios_base::cur, std::ios_base::out) == 0);

  std::vector<Double_t> Buffer; // Collect batches of noise information.

  for(size_t f = 1; f <= 1024; f++) {
    Buffer.clear();

    // First do real columns.
    for(size_t index1 = 0; index1 < ChannelsToUse.size(); index1++) {
      unsigned char noiseIndex1 = noise->GetIndexOfChannel(ChannelsToUse[index1]);

      // First do real rows.
      for(size_t index2 = 0; index2 < ChannelsToUse.size(); index2++) {
        unsigned char noiseIndex2 = noise->GetIndexOfChannel(ChannelsToUse[index2]);
        Buffer.push_back(noise->GetRR(f)[noiseIndex2][noiseIndex1]);
      }

      // Now do imag rows.
      if(f == 1024) continue;
      for(size_t index2 = 0; index2 < ChannelsToUse.size(); index2++) {
        unsigned char noiseIndex2 = noise->GetIndexOfChannel(ChannelsToUse[index2]);
        Buffer.push_back(noise->GetRI(f)[noiseIndex1][noiseIndex2]);
      }
    }

    // Then do imag columns.
    if(f < 1024) {
    for(size_t index1 = 0; index1 < ChannelsToUse.size(); index1++) {
      unsigned char noiseIndex1 = noise->GetIndexOfChannel(ChannelsToUse[index1]);

      // First do real rows.
      for(size_t index2 = 0; index2 < ChannelsToUse.size(); index2++) {
        unsigned char noiseIndex2 = noise->GetIndexOfChannel(ChannelsToUse[index2]);
        Buffer.push_back(noise->GetRI(f)[noiseIndex2][noiseIndex1]);
      }

      // Now do imag rows.
      if(f == 1024) continue;
      for(size_t index2 = 0; index2 < ChannelsToUse.size(); index2++) {
        unsigned char noiseIndex2 = noise->GetIndexOfChannel(ChannelsToUse[index2]);
        Buffer.push_back(noise->GetII(f)[noiseIndex1][noiseIndex2]);
      }
    }
    }

    // Push the buffer into the output file.
    size_t numChars = 8*Buffer.size();
    assert(Buffer.size() == ChannelsToUse.size()*ChannelsToUse.size()*(f == 1024 ? 1 : 4));
    assert(outfile.sputn((char*)&Buffer[0], numChars) == numChars);
  }
  assert(outfile.pubseekoff(0, std::ios_base::cur, std::ios_base::out) ==
         8*ChannelsToUse.size()*ChannelsToUse.size()*(4*1023 + 1));
  assert(outfile.pubseekoff(0, std::ios_base::end, std::ios_base::out) ==
         8*ChannelsToUse.size()*ChannelsToUse.size()*(4*1023 + 1));
  assert(outfile.close() != NULL);
}
