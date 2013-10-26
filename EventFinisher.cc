#include "EventFinisher.hh"
#include "SafeStopwatch.hh"
#include "EXOUtilities/EXOWaveform.hh"
#include "EXOUtilities/EXOFastFourierTransformFFTW.hh"

#ifdef USE_THREADS
// So we can lower the priority of the FinishEvent thread.
#include <unistd.h>
#include <sys/syscall.h>
#include <sys/resource.h>
#endif

EventFinisher& EventFinisher::Get(EXOTreeInputModule& inputModule,
                                  std::string RawFileName,
                                  std::string OutFileName)
{
  static EventFinisher gEventFinisher(inputModule, RawFileName, OutFileName);
  return gEventFinisher;
}

EventFinisher::EventFinisher(EXOTreeInputModule& inputModule, std::string RawFileName, std::string OutFileName)
: fVerbose(true),
  fInputModule(inputModule),
  fWaveformFile(RawFileName.c_str()),
#ifdef USE_THREADS
  fDesiredQueueLength(2000),
  fProcessingIsDone(false),
#endif
  fThoriumEnergy_keV(2615)
{
  fWaveformTree = dynamic_cast<TTree*>(fWaveformFile.Get("tree"));
  fWaveformTree->GetBranch("fWaveformData")->SetAddress(&fWFData);
  fOutputModule.SetOutputFilename(OutFileName);
  fOutputModule.Initialize();
  fOutputModule.BeginOfRun(NULL);
}

EventFinisher::~EventFinisher()
{
  fOutputModule.ShutDown();
}

void EventFinisher::QueueEvent(EventHandler* eventHandler)
{
  // Take ownership of an event; don't actually finish it yet,
  // but insert it into out set of events to finish.
#ifdef USE_THREADS
  boost::mutex::scoped_lock sL(fFinisherMutex);
#endif
  fEventsToFinish.insert(eventHandler);
#ifdef USE_THREADS
  assert(not fProcessingIsDone); // Certainly shouldn't be done, if we got here.
  if(fEventsToFinish.size() >= fDesiredQueueLength) fFinisherCondition.notify_one();
#endif
}

void EventFinisher::FinishProcessedEvent(EventHandler* event, const std::vector<double>& Results)
{
  // Grab the processed event from input; apply denoised results as necessary; and write to output.
  // We will also delete event here.
  // This function is particularly useful when called directly for events with no denoising to do;
  // in that case, it allows us to complete handling of an event without unnecessary lock management.
  // Prerequisite: RootInterfaceMutex must be locked by the caller.
  if(fVerbose) std::cout<<"\tFinishProcessedEvent for entry "<<event->fEntryNumber<<std::endl;

  // If this function is called by AcceptEvent, trees are fast at re-returning an already-gotten entry.
  EXOEventData* ED = fInputModule.GetEvent(event->fEntryNumber);

  // We need to clear out the denoised information here, since we just freshly read the event from file.
  for(size_t i = 0; i < ED->GetNumScintillationClusters(); i++) {
    ED->GetScintillationCluster(i)->fEnergy = ED->GetScintillationCluster(i)->fRawEnergy;
    ED->GetScintillationCluster(i)->fRawEnergy = 0;
    ED->GetScintillationCluster(i)->fDenoisedEnergy = 0;
  }
#ifdef ENABLE_CHARGE
  for(size_t i = 0; i < ED->GetNumUWireSignals(); i++) {
    ED->GetUWireSignal(i)->fDenoisedEnergy = 0;
  }
  for(size_t i = 0; i < ED->GetNumChargeClusters(); i++) {
    ED->GetChargeCluster(i)->fDenoisedEnergy = 0;
  }
#endif

  if(not Results.empty()) {
    // Translate signal magnitudes into corresponding objects.
#ifdef ENABLE_CHARGE
    if(not fAPDsOnly) {
      for(size_t i = 0; i < event->fWireModel.size(); i++) {
        size_t sigIndex = event->fWireModel[i].first;
        EXOUWireSignal* sig = ED->GetUWireSignal(sigIndex);
        double UWireScalingFactor = ADC_FULL_SCALE_ELECTRONS_WIRE * W_VALUE_LXE_EV_PER_ELECTRON /
                                    (CLHEP::keV * ADC_BITS);
        sig->fDenoisedEnergy = Results[i]*UWireScalingFactor;
      }
    }
#endif
    ED->GetScintillationCluster(0)->fDenoisedEnergy = Results.back()*fThoriumEnergy_keV;
    ED->GetScintillationCluster(0)->fRawEnergy = ED->GetScintillationCluster(0)->fDenoisedEnergy;
  }

  fOutputModule.ProcessEvent(ED);
  if(fVerbose) std::cout<<"\tDone with entry "<<event->fEntryNumber<<std::endl;
  delete event;
}

void EventFinisher::FinishEvent(EventHandler* event)
{
  // Compute and fill denoised signals, as appropriate.
  // Then pass the filled event to the output module.
  if(fVerbose) std::cout<<"Finishing entry "<<event->fEntryNumber<<std::endl;

  std::vector<double> Results;
  if(not event->fX.empty()) {
    if(fVerbose) std::cout<<"\tThis entry has denoised results to compute."<<std::endl;
    Results.assign(event->fNumSignals, 0);

    // We need to compute denoised signals.
    static SafeStopwatch GetRawWatch("FinishEvent::GetRawEntry (threaded, mostly)");
    SafeStopwatch::tag GetRawTag = GetRawWatch.Start();
    Long64_t RawEntryNum = fWaveformTree->GetEntryNumberWithIndex(event->fRunNumber, event->fEventNumber);
    fWaveformTree->GetBranch("fWaveformData")->GetEntry(RawEntryNum);
    GetRawWatch.Stop(GetRawTag);
    fWFData.Decompress();

    // Collect the fourier-transformed waveforms.  Save them split into real and complex parts.
    std::vector<EXODoubleWaveform> WF_real, WF_imag;
#ifdef USE_THREADS
    FFTWMutex.lock();
#endif
    for(size_t i = 0; i < event->fChannels.size(); i++) {
      const EXOWaveform* wf = fWFData.GetWaveformWithChannel(event->fChannels[i]);

      // Take the Fourier transform.
      EXODoubleWaveform dwf = wf->Convert<Double_t>();
      EXOWaveformFT fwf;
      // Note: not thread-safe (we use a buffer), but covered by RootInterfaceMutex.
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
#ifdef USE_THREADS
    FFTWMutex.unlock();
#endif

    // Produce estimates of the signals.
    for(size_t i = 0; i < Results.size(); i++) {
      for(size_t f = 0; f <= event->fMaxF - event->fMinF; f++) {
        for(size_t chan_index = 0; chan_index < event->fChannels.size(); chan_index++) {
          size_t XIndex = event->fColumnLength*i + 2*event->fChannels.size()*f + chan_index;
          Results[i] += event->fX[XIndex]*WF_real[chan_index][f + event->fMinF];
        }
        if(f == event->fMaxF - event->fMinF) continue;
        for(size_t chan_index = 0; chan_index < event->fChannels.size(); chan_index++) {
          size_t XIndex = event->fColumnLength*i + 2*event->fChannels.size()*f + event->fChannels.size() + chan_index;
          Results[i] += event->fX[XIndex]*WF_imag[chan_index][f + event->fMinF];
        }
      }
    }
  } // End setting of denoised energy signals.

#ifdef USE_THREADS
  boost::mutex::scoped_lock sL(RootInterfaceMutex);
#endif
  FinishProcessedEvent(event, Results); // Deletes event.
}

void EventFinisher::Run()
{
  // Call FinishEvent repeatedly until the queue is empty.
  // When the queue is empty, either return or sleep until more events are available.

#ifdef USE_THREADS
  // Reduce thread priority -- we'd like for the threads which insert new events to get priority on the lock.
  pid_t tid = syscall(SYS_gettid);
  int main_priority = getpriority(PRIO_PROCESS, tid);
  std::cout<<"FinishEventThread has priority "<<main_priority<<"; reducing to "<<main_priority+5<<std::endl;
  setpriority(PRIO_PROCESS, tid, main_priority+5);
#endif

  while(true) {
#ifdef USE_THREADS
    boost::mutex::scoped_lock sL(fFinisherMutex);

    // Force ourselves to wait until the finish-events queue has a certain length (or processing is done).
    while(not fProcessingIsDone and
          fEventsToFinish.size() < fDesiredQueueLength) {
      fFinisherCondition.wait(sL);
    }
#endif
    std::set<EventHandler*, CompareEventHandlerPtrs>::iterator it = fEventsToFinish.begin();
    if(it == fEventsToFinish.end()) {
#ifdef USE_THREADS
      // We assume only one io thread; in that case,
      // if we weren't woken because there were events to finish,
      // then we'd better have gotten here because processing is done.
      assert(fProcessingIsDone);
#endif
      // Both sequential and threaded code should only get here because they're done.
      return;
    }

    EventHandler* event = *it;
    fEventsToFinish.erase(it);

#ifdef USE_THREADS
    sL.unlock();
#endif
    FinishEvent(event);
  } // while(true)
}

#ifdef USE_THREADS
void EventFinisher::SetProcessingIsFinished()
{
  // Processing is completed, set and notify
  boost::mutex::scoped_lock sL(fFinisherMutex);
  fProcessingIsDone = true;
  fFinisherCondition.notify_all();
}
#endif

size_t EventFinisher::GetFinishEventQueueLength()
{
  // Number of events currently queued to be finished.
#ifdef USE_THREADS
  boost::mutex::scoped_lock sL(fFinisherMutex);
#endif
  return fEventsToFinish.size();
}
