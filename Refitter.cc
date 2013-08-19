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
#include <iostream>


int main(int argc, char** argv)
{
  if(argc != 4) {
    std::cout<<"Usage: ./Refit <InputProcessedFile> <InputWaveformFile> <OutputFile>"<<std::endl;
    std::exit(1);
  }

  EXOTreeInputModule InputModule;
  InputModule.SetFilename(argv[1]);
  TFile WaveformFile(argv[2]);
  TTree* WaveformTree = WaveformFile.Get("tree");

  EXOTreeOutputModule OutputModule;
  OutputModule.SetFilename(argv[3]);
  OutputModule.BeginOfRun(NULL); // OK, fine -- shortcut here, I assume input has only one run.

  int MaxEvents = 5; // Hard-code it for now.
  EXORefitSignals RefitSig(InputModule, *WaveformTree, OutputModule);

  while((EXOEventData* ED = InputModule.GetNext())) {

    // Start by resetting old values.
    for(size_t i = 0; i < ED->GetNumScintillationClusters(); i++) {
      ED->GetScintillationCluster(i)->fDenoisedEnergy = 0;
    }
    for(size_t i = 0; i < ED->GetNumUWireSignals(); i++) {
      ED->GetUWireSignal(i)->fDenoisedEnergy = 0;
    }
    for(size_t i = 0; i < ED->GetNumChargeClusters(); i++) {
      ED->GetChargeCluster(i)->fDenosiedEnergy = 0;
    }

    // Pass event into EXORefitSignals.
    std::list<EventHandler*> doneEvents = RefitSig.Process(ED);

    // Handle finished events.
    for(std::list<EventHandler*>::iterator it = doneEvents.begin(); it != doneEvents.end(); it++) {
      EventHandler* result = *it;
      if(result->fWasHandled) {
        // Apply result to actually compute the various signals.

        // Collect the fourier-transformed waveforms.  Save them split into real and complex parts.
        // Skip channels which aren't included in our noise or lightmap models, but warn.
        std::vector<EXODoubleWaveform> WF_real, WF_imag;
        for(size_t i = 0; i < RefitSig.fChannels.size(); i++) {
          const EXOWaveform* wf = ED->GetWaveformData()->GetWaveformWithChannel(fChannels[i]);
    if(not wf) LogEXOMsg("A waveform disappeared!", EEAlert);

    // Take the Fourier transform.
    EXODoubleWaveform dwf = wf->Convert<Double_t>();
    EXOWaveformFT fwf;
    EXOFastFourierTransformFFTW::GetFFT(dwf.GetLength()).PerformFFT(dwf, fwf);





