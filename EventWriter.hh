#ifndef EventWriter_hh
#define EventWriter_hh

/*
This class receives ready-to-write events.
It reads from the processed file at SLAC, fills in the appropriate values from denoising,
and writes out a new entry.
It will enforce entry ordering; so, it holds on to EventHandlers until their entry number is the
next one we want to write.
The event handlers should be automatically cleaned up when they're finished here; otherwise,
ownership has somehow been retained elsewhere (which would be bad).

Don't let this run in parallel with the raw-event reader -- this leads to issues with ROOT.
*/

#include "EventHandler.hh"
#include "Constants.hh"
#include "EXOUtilities/EXODimensions.hh"
#include "EXOUtilities/EXOEventData.hh"
#include "EXOUtilities/EXOScintillationCluster.hh"
#ifdef ENABLE_CHARGE
#include "EXOUtilities/EXOChargeCluster.hh"
#include "EXOUtilities/EXOUWireSignal.hh"
#endif
#include "EXOAnalysisManager/EXOTreeInputModule.hh"
#include "EXOAnalysisManager/EXOTreeOutputModule.hh"
#include "Rtypes.h"
#include <set>
#include <memory>
#include <cassert>

class EventWriter
{
 private:
  // Here we order events by entry number, since we want to duplicate that entry ordering.
  struct CompareEventHandlerPtrs_byentry {
    bool operator()(const std::shared_ptr<EventHandler> &evt1,
                    const std::shared_ptr<EventHandler> &evt2) {
      return evt1->fEntryNumber < evt2->fEntryNumber;
    }
  };
  std::set<std::shared_ptr<EventHandler>, CompareEventHandlerPtrs_byentry> fEventsToWrite;
  Long64_t fNextEntryNumber;
  EXOTreeInputModule& fInputModule;
  EXOTreeOutputModule fOutputModule;

  void WriteEvent(const std::shared_ptr<EventHandler> event) {
    // Grab processed entry; fill in denoised information; and write out.
    // Verify that the run/event numbers match too, as an end-to-end check.
    std::cout<<"\tWriteEvent for entry "<<event->fEntryNumber<<std::endl;
    EXOEventData* ED = fInputModule.GetEvent(event->fEntryNumber);
    assert(ED and ED->fRunNumber == event->fRunNumber and ED->fEventNumber == event->fEventNumber);

    // We need to clear out the denoised information here, since we just freshly read the event from file.
    for(size_t i = 0; i < ED->GetNumScintillationClusters(); i++) {
      ED->GetScintillationCluster(i)->fEnergy = ED->GetScintillationCluster(i)->fRawEnergy;
      ED->GetScintillationCluster(i)->fRawEnergy = 0;
      ED->GetScintillationCluster(i)->fDenoisedEnergy = 0;
      ED->GetScintillationCluster(i)->fDenoisingInternalCode = event->fStatusCode;
    }
#ifdef ENABLE_CHARGE
    for(size_t i = 0; i < ED->GetNumUWireSignals(); i++) {
      ED->GetUWireSignal(i)->fDenoisedEnergy = 0;
    }
    for(size_t i = 0; i < ED->GetNumChargeClusters(); i++) {
      ED->GetChargeCluster(i)->fDenoisedEnergy = 0;
    }
#endif

    // Fill denoised information that is available.
    if(not event->fResults.empty()) {
      // Translate signal magnitudes into corresponding objects.
#ifdef ENABLE_CHARGE
      for(size_t i = 0; i < event->fWireModel.size(); i++) {
        size_t sigIndex = event->fWireModel[i].fSignalNumber;
        EXOUWireSignal* sig = ED->GetUWireSignal(sigIndex);
        double UWireScalingFactor = ADC_FULL_SCALE_ELECTRONS_WIRE * W_VALUE_LXE_EV_PER_ELECTRON /
                                    (CLHEP::keV * ADC_BITS);
        sig->fDenoisedEnergy = event->fResults[event->fAPDModel.size() + i]*UWireScalingFactor;
      }
#endif
      for(size_t i = 0; i < event->fAPDModel.size(); i++) {
        size_t sigIndex = event->fAPDModel[i].fSignalNumber;
        ED->GetScintillationCluster(sigIndex)->fDenoisedEnergy = event->fResults[i]*THORIUM_ENERGY_KEV;
        ED->GetScintillationCluster(sigIndex)->fRawEnergy = ED->GetScintillationCluster(sigIndex)->fDenoisedEnergy;
      }
    }

    // Finish.
    fOutputModule.ProcessEvent(ED);
    std::cout<<"\tDone with entry "<<event->fEntryNumber<<std::endl;
    fNextEntryNumber = event->fEntryNumber+1;
  }

 public:
  EventWriter(EXOTreeInputModule& inputModule,
              std::string OutFileName)
  : fNextEntryNumber(0),
    fInputModule(inputModule)
  {
    fOutputModule.SetOutputFilename(OutFileName);
    fOutputModule.Initialize();
    fOutputModule.BeginOfRun(NULL);
  }

  void AcceptEvent(const std::shared_ptr<EventHandler> event) {
    // Accept an event.
    assert(fNextEntryNumber <= event->fEntryNumber);
    if(event->fEntryNumber == fNextEntryNumber) {
      WriteEvent(event);
      auto it = fEventsToWrite.begin();
      while(it != fEventsToWrite.end() and fNextEntryNumber == (*it)->fEntryNumber) {
        WriteEvent(*it);
        fEventsToWrite.erase(it);
        it = fEventsToWrite.begin();
      }
    }
    else fEventsToWrite.insert(event);
  }

  void Finish() {
    // Verify that we really did write everything -- if not, somehow an entry got lost.
    assert(fEventsToWrite.empty());
    fOutputModule.ShutDown();
  }
};
#endif
