/*
Refit-signals was originally written as an EXOAnalysis module; however, it quickly became clear that
the code could work much more efficiently if we collected many events together to handle in parallel.
The bottleneck to the code is multiplying vectors by the noise blocks of A

Note that this module DOES rearrange the ordering of events.  If you're not referring to events by
run/event number, well, you should.  (We index our trees by run/event number, for your convenience.)

The usability is less than EXOAnalysis, with fewer options settable on the command line.
It's possible I'll address this, but no promises.

Should be called like:
./Refit <InputProcessedFile> <InputWaveformFile> <OutputFile>
*/


#include "EXORefitSignals.hh"
#include "EXOUtilities/EXOEventData.hh"
#include "EXOAnalysisManager/EXOTreeInputModule.hh"
#include "EXOAnalysisManager/EXOTreeOutputModule.hh"
#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <cstdlib>


int main(int argc, char** argv)
{
  if(argc < 4) {
    std::cout<<"Usage: ./Refit <InputProcessedFile> <InputWaveformFile> <OutputFile>"<<std::endl;
    std::exit(1);
  }
  Long64_t StartEntry = 0;
  if(argc >= 5) StartEntry = std::atol(argv[4]);
  Long64_t NumEntries = 40;
  if(argc >= 6) NumEntries = std::atol(argv[5]);

  EXOTreeInputModule InputModule;
  InputModule.SetFilename(argv[1]);
  TFile WaveformFile(argv[2]);
  TTree* WaveformTree = dynamic_cast<TTree*>(WaveformFile.Get("tree"));

  EXOTreeOutputModule OutputModule;
  OutputModule.SetOutputFilename(argv[3]);
  OutputModule.Initialize();
  OutputModule.BeginOfRun(NULL); // OK, fine -- shortcut here, I assume input has only one run.

  EXORefitSignals RefitSig(InputModule, *WaveformTree, OutputModule);
  RefitSig.SetNoiseFilename("/global/u1/c/claytond/TestEXOAnalysisBuild/noise_manyruns_withuwires_100000.dat");
  RefitSig.SetRThreshold(0.1);
  RefitSig.Initialize();

  for(Long64_t entryNum = StartEntry; entryNum < StartEntry + NumEntries; entryNum++) {
    EXOEventData* ED = InputModule.GetEvent(entryNum);
    if(ED == NULL) break;
    RefitSig.AcceptEvent(ED, entryNum);
  }

  RefitSig.FlushEvents();
}
