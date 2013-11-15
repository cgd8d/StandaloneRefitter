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

class EventFinisher
{
 public:
  static EventFinisher& Get(EXOTreeInputModule& inputModule, std::string RawFileName, std::string OutFileName);

  void QueueEvent(EventHandler* eventHandler);
  void FinishProcessedEvent(EventHandler* event, const std::vector<double>& Results = std::vector<double>());

  void Run();

  void ListenForArrivingEvents();

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
  size_t fDesiredQueueLength;
  bool fProcessingIsDone;
  bool fHasAskedForPause;
};
#endif
