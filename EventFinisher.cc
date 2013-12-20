#include "EventFinisher.hh"
#include "SafeStopwatch.hh"
#include "Constants.hh"
#include "EXOUtilities/EXOWaveform.hh"
#include "EXOUtilities/EXOFastFourierTransformFFTW.hh"

#include <boost/thread/thread.hpp>
#include <boost/date_time/posix_time/posix_time_types.hpp>

#include <boost/mpi/communicator.hpp>
#include <memory>

EventFinisher& EventFinisher::Get(EXOTreeInputModule& inputModule,
                                  std::string RawFileName,
                                  std::string OutFileName)
{
  static EventFinisher gEventFinisher(inputModule, RawFileName, OutFileName);
  return gEventFinisher;
}

EventFinisher::EventFinisher(EXOTreeInputModule& inputModule, std::string RawFileName, std::string OutFileName)
: fVerbose(true),
  fWaveformFile(RawFileName.c_str()),
  fDesiredQueueLength(2000),
  fProcessingIsDone(false),
  fHasAskedForPause(false),
  fWriter(inputModule, OutFileName)
{
  fWaveformTree = dynamic_cast<TTree*>(fWaveformFile.Get("tree"));
  fWaveformTree->GetBranch("fWaveformData")->SetAddress(&fWFData);
}

void EventFinisher::QueueEvent(EventHandler* eventHandler)
{
  // Take ownership of an event; don't actually finish it yet,
  // but insert it into out set of events to finish.
  boost::mutex::scoped_lock sL(fEventsToFinishMutex);
  fEventsToFinish.insert(eventHandler);
  if(fVerbose) std::cout<<"Queued an entry; queue length is "<<fEventsToFinish.size()<<std::endl;

  // If we've queued too many events, ask the compute process to give us a chance to catch up.
  if(not fHasAskedForPause and fEventsToFinish.size() > 5000) {
    sL.unlock();
    if(fVerbose) std::cout<<"Asking the compute process to pause for a bit."<<std::endl;
    static boost::mpi::communicator gMPIComm;
    gMPIComm.send(gMPIComm.rank() - 1, 1); // Please be merciful, and pause for a bit.
    fHasAskedForPause = true;
  }
}

void EventFinisher::FinishEvent(EventHandler* event)
{
  // Compute and fill denoised signals, as appropriate.
  // Then pass the filled event to the output module.
  if(fVerbose) std::cout<<"Finishing entry "<<event->fEntryNumber<<std::endl;

  if(not event->fX.empty()) {
    if(fVerbose) std::cout<<"\tThis entry has denoised results to compute."<<std::endl;
    event->fResults.assign(event->fNumSignals, 0);

    // We need to compute denoised signals.
    static SafeStopwatch GetRawWatch("FinishEvent::GetRawEntry (threaded, mostly)");
    SafeStopwatch::tag GetRawTag = GetRawWatch.Start();
    Long64_t RawEntryNum = fWaveformTree->GetEntryNumberWithIndex(event->fRunNumber, event->fEventNumber);
    fWaveformTree->GetBranch("fWaveformData")->GetEntry(RawEntryNum);
    GetRawWatch.Stop(GetRawTag);

    static SafeStopwatch RestOfRawWatch("FinishEvent::RestOfRaw");
    SafeStopwatch::tag RestOfRawTag = RestOfRawWatch.Start();
    fWFData.Decompress();

    // Collect the fourier-transformed waveforms.  Save them split into real and complex parts.
    std::vector<EXODoubleWaveform> WF_real, WF_imag;
    for(size_t i = 0; i < event->fChannels.size(); i++) {
      const EXOWaveform* wf = fWFData.GetWaveformWithChannel(event->fChannels[i]);

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
    for(size_t i = 0; i < event->fResults.size(); i++) {
      for(size_t f = 0; f <= MAX_F - MIN_F; f++) {
        for(size_t chan_index = 0; chan_index < event->fChannels.size(); chan_index++) {
          size_t XIndex = event->fColumnLength*i + 2*event->fChannels.size()*f + chan_index;
          event->fResults[i] += event->fX[XIndex]*WF_real[chan_index][f + MIN_F];
        }
        if(f == MAX_F - MIN_F) continue;
        for(size_t chan_index = 0; chan_index < event->fChannels.size(); chan_index++) {
          size_t XIndex = event->fColumnLength*i + 2*event->fChannels.size()*f + event->fChannels.size() + chan_index;
          event->fResults[i] += event->fX[XIndex]*WF_imag[chan_index][f + MIN_F];
        }
      }
    }
    RestOfRawWatch.Stop(RestOfRawTag);
    std::vector<double>().swap(event->fX); // Done with fX, so force-clear it.
  } // End setting of denoised energy signals.

  static SafeStopwatch FinishProcessedWatch("FinishProcessedWatch");
  SafeStopwatch::tag FinishProcessedTag = FinishProcessedWatch.Start();
  fWriter.AcceptEvent(std::shared_ptr<EventHandler>(event));
  FinishProcessedWatch.Stop(FinishProcessedTag);
}

void EventFinisher::Run()
{
  // Call FinishEvent repeatedly until the queue is empty.
  // When the queue is empty, either return or sleep until more events are available.

#ifdef USE_THREADS
  boost::thread finish_data(&EventFinisher::FinishReceivedEvents, this);
  ListenForArrivingEvents(); // This needs to be the in the main thread, for MPI to be "funneled".
  finish_data.join();
#else
  while(not (fProcessingIsDone and fEventsToFinish.empty())) {
    if(not fProcessingIsDone) ListenForArrivingEvents();
    FinishReceivedEvents();
  }
#endif
}

void EventFinisher::FinishReceivedEvents()
{
#ifdef USE_THREADS
  while(true) {

    boost::mutex::scoped_lock sL(fEventsToFinishMutex);

    // Sleep until we have enough events queued to start.
    static SafeStopwatch watch("Waiting to process more events");
    SafeStopwatch::tag tag = watch.Start();
    while(fEventsToFinish.size() < fDesiredQueueLength and not fProcessingIsDone) {
      sL.unlock();
      boost::this_thread::sleep_for(boost::chrono::seconds(1));
      sL.lock();
    }
    watch.Stop(tag);

    // As long as we have the lock, might as well check for being totally done.
    if(fEventsToFinish.empty()) {
      assert(fProcessingIsDone);
      fWriter.Finish();
      return;
    }
#endif

    assert(not fEventsToFinish.empty()); // If threaded, we slept; if not, listener guarantees this.
    EventHandler* event = *fEventsToFinish.begin();
    fEventsToFinish.erase(fEventsToFinish.begin());
#ifdef USE_THREADS
    sL.unlock();
#endif
    FinishEvent(event);
#ifdef USE_THREADS
  } // while(true)
#endif
}

void EventFinisher::ListenForArrivingEvents()
{
  static boost::mpi::communicator gMPIComm;
  static EventHandler* eh = new EventHandler;
  static boost::mpi::request req = gMPIComm.irecv( gMPIComm.rank() - 1, boost::mpi::any_tag, *eh );

  while (1)
  {
    static SafeStopwatch MPITestWatch("MPI_Test");
    SafeStopwatch::tag MPITestTag = MPITestWatch.Start();
    boost::optional<boost::mpi::status> s_opt = req.test();
    MPITestWatch.Stop(MPITestTag);
    if(s_opt.is_initialized()) {
      // We received a message.
      boost::mpi::status& s = s_opt.get();
      if(s.tag()) {
        // This was the signal that we're done listening.
        fProcessingIsDone = true;
        delete eh;
        return;
      }
      else {
        // Receive an event, then get ready for the next one.
        static SafeStopwatch watch("QueueEvent in listener");
        SafeStopwatch::tag tag = watch.Start();
        QueueEvent(eh);
        eh = new EventHandler;
        req = gMPIComm.irecv( gMPIComm.rank() - 1, boost::mpi::any_tag, *eh );
        watch.Stop(tag);
        continue;
      }
    }
    else {
      // We haven't received a message.

      {
        // If we're not receiving messages, it could be because we requested a pause.
        // Check to see if we've caught up.
        // Note: it is important to only interact with MPI from a single thread.
        // boost::mpi only started supporting multithreaded MPI in version 1.55,
        // and we'll need to explicitly enable it (with some performance penalty) if we want it.
        static SafeStopwatch watch("Getting lock in listener");
        SafeStopwatch::tag tag = watch.Start();
        boost::mutex::scoped_lock sL(fEventsToFinishMutex);
        if(fHasAskedForPause and fEventsToFinish.size() < 4000) {
          if(fVerbose) std::cout<<"Letting the compute process know that it can continue."<<std::endl;
          gMPIComm.send(gMPIComm.rank() - 1, 0); // OK for compute process to continue.
          fHasAskedForPause = false;
        }
        watch.Stop(tag);
      }

// If we're running with just one thread, then we only want to wait if
// there's no IO we could be doing.
// Otherwise, this thread does nothing but listen -- so wait for sure.
#ifndef USE_THREADS
      // If we're not using threads, we need to decide whether to return or wait here.
      // If we could be writing events, go back to doing that; else, wait.
      if(fEventsToFinish.size() < fDesiredQueueLength) {
#endif
        // Wait.
        static SafeStopwatch watch("Waiting in listener");
        SafeStopwatch::tag tag = watch.Start();
        // We have no choice but to wait.
        boost::this_thread::yield();
        watch.Stop(tag);

// On the other hand, if we're not using threads and we could be doing IO, return.
#ifndef USE_THREADS
      }
      else {
        // We should go back to io.
        return;
      }
#endif
    }
  }
}
