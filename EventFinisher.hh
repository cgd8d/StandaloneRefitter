#ifndef EventFinisher_hh
#define EventFinisher_hh

#include "EventHandler.hh"
#include "EXOAnalysisManager/EXOTreeInputModule.hh"
#include "EXOAnalysisManager/EXOTreeOutputModule.hh"
#include "EXOUtilities/EXOWaveformData.hh"
#include "TXNetFile.h"
#include "TTree.h"
#include <string>
#include <set>

#ifdef USE_THREADS
#include <boost/thread/thread.hpp>
#ifndef USE_PROCESSES
extern boost::mutex RootInterfaceMutex;
extern boost::mutex FFTWMutex;
#endif
#endif

class EventFinisher
{
 public:
  static EventFinisher& Get(EXOTreeInputModule& inputModule, std::string RawFileName, std::string OutFileName);
  ~EventFinisher();

  void QueueEvent(EventHandler* eventHandler);
  void FinishProcessedEvent(EventHandler* event, const std::vector<double>& Results = std::vector<double>());

  void Run();
  size_t GetFinishEventQueueLength();

#ifdef USE_PROCESSES
  void ListenForArrivingEvents();
#endif

#ifdef USE_THREADS
  void SetProcessingIsFinished();
  boost::mutex fFinisherMutex;
  boost::condition_variable fFinisherCondition;
#endif

  bool fVerbose;
 private:
  void FinishEvent(EventHandler* event);

  EventFinisher(EXOTreeInputModule& inputModule, std::string RawFileName, std::string OutFileName);
  EXOTreeInputModule& fInputModule;
  EXOTreeOutputModule fOutputModule;
  TXNetFile fWaveformFile;
  TTree* fWaveformTree;
  EXOWaveformData fWFData;
  std::set<EventHandler*, CompareEventHandlerPtrs> fEventsToFinish;
#if defined(USE_THREADS) || defined(USE_PROCESSES)
  size_t fDesiredQueueLength;
  bool fProcessingIsDone;
#endif
};
#endif
