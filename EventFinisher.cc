#include "EventFinisher.hh"
#include "SafeStopwatch.hh"
#include "Constants.hh"
#include "EXOUtilities/EXOWaveform.hh"
#include "EXOUtilities/EXOFastFourierTransformFFTW.hh"

#ifdef USE_PROCESSES
#include <boost/date_time/posix_time/posix_time_types.hpp>
#include <boost/thread/thread.hpp>
#endif

#ifdef USE_PROCESSES
#include <boost/mpi/communicator.hpp>
static boost::mpi::communicator gMPIComm;
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
  fWaveformFile(RawFileName.c_str())
#if defined(USE_THREADS) || defined(USE_PROCESSES)
  ,
  fDesiredQueueLength(2000),
  fProcessingIsDone(false)
#endif
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
  if(fVerbose) std::cout<<"Queued an entry; queue length is "<<fEventsToFinish.size()<<std::endl;
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
    ED->GetScintillationCluster(0)->fDenoisedEnergy = Results.back()*THORIUM_ENERGY_KEV;
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
#if defined(USE_THREADS) && !defined(USE_PROCESSES)
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
#if defined(USE_THREADS) && !defined(USE_PROCESSES)
    FFTWMutex.unlock();
#endif

    // Produce estimates of the signals.
    for(size_t i = 0; i < Results.size(); i++) {
      for(size_t f = 0; f <= MAX_F - MIN_F; f++) {
        for(size_t chan_index = 0; chan_index < event->fChannels.size(); chan_index++) {
          size_t XIndex = event->fColumnLength*i + 2*event->fChannels.size()*f + chan_index;
          Results[i] += event->fX[XIndex]*WF_real[chan_index][f + MIN_F];
        }
        if(f == MAX_F - MIN_F) continue;
        for(size_t chan_index = 0; chan_index < event->fChannels.size(); chan_index++) {
          size_t XIndex = event->fColumnLength*i + 2*event->fChannels.size()*f + event->fChannels.size() + chan_index;
          Results[i] += event->fX[XIndex]*WF_imag[chan_index][f + MIN_F];
        }
      }
    }
  } // End setting of denoised energy signals.

#if defined(USE_THREADS) && !defined(USE_PROCESSES)
  boost::mutex::scoped_lock sL(RootInterfaceMutex);
#endif
  FinishProcessedEvent(event, Results); // Deletes event.
}

void EventFinisher::Run()
{
  // Call FinishEvent repeatedly until the queue is empty.
  // When the queue is empty, either return or sleep until more events are available.

#if defined(USE_THREADS) && defined(USE_PROCESSES)
  boost::thread pull_data(&EventFinisher::ListenForArrivingEvents, this);
#endif

  while(true) {
#ifdef USE_THREADS
    boost::mutex::scoped_lock sL(fFinisherMutex);
#endif

#if defined(USE_THREADS) || defined(USE_PROCESSES)
    // Force ourselves to wait until the finish-events queue has a certain length (or processing is done).
    while(not fProcessingIsDone and
          fEventsToFinish.size() < fDesiredQueueLength) {
#ifdef USE_THREADS
      fFinisherCondition.wait(sL);
#else // USE_PROCESSES but not USE_THREADS -- we have to ask for another event within this thread.
      ListenForArrivingEvents();
#endif
    }
#endif
    std::set<EventHandler*, CompareEventHandlerPtrs>::iterator it = fEventsToFinish.begin();
    if(it == fEventsToFinish.end()) {
#if defined(USE_THREADS) || defined(USE_PROCESSES)
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

#ifdef USE_PROCESSES
void EventFinisher::ListenForArrivingEvents()
{
#ifdef USE_THREADS
  // If we're using threads, we'll busy-wait until done.
  while (1)
#endif // If we're not using threads, just receive once.
  {
    EventHandler* eh = new EventHandler;
    boost::mpi::request req = gMPIComm.irecv( gMPIComm.rank() - 1, boost::mpi::any_tag, *eh );
    boost::optional<boost::mpi::status> s; 
    while ( ! (s = req.test()) ) {
      boost::this_thread::sleep(boost::posix_time::millisec(1));
    }
    if(s->tag()) {
#ifdef USE_THREADS
      SetProcessingIsFinished();
#else
      ProcessingIsDone = true;
#endif
      delete eh;
      break;
    }
    QueueEvent(eh);
  }
}
#endif
