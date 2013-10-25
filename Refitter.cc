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
#include "EXOCalibUtilities/EXOCalibManager.hh"
#include "EXOAnalysisManager/EXOTreeInputModule.hh"
#include "EXOAnalysisManager/EXOTreeOutputModule.hh"
#include "TFile.h"
#include "TXNetFile.h"
#include "TTree.h"
#include <iostream>

#ifdef USE_THREADS
#include <boost/thread/thread.hpp>
#include <boost/chrono/duration.hpp>
#endif

#ifdef USE_MPI
#include <fstream>
#include <iomanip>
#include <mpi.h>
// To control order of destruction -- we'd like for all function-static objects to be destroyed
// before MPI_Finalize is called, so their destructors have a chance to execute side-effects.
// Also redirect output properly at this stage.
struct mpi_handler
{
  int rank;
  std::string RankString;
  std::ofstream Output;
  mpi_handler(int argc, char** argv) {
    assert(argc >= 2);
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::ostringstream ProcessRankString;
    ProcessRankString << std::setw(4) << std::setfill('0') << rank << ".txt";
    RankString = ProcessRankString.str();
    Output.open(std::string(argv[1]) + "/outfile" + RankString);
    std::cout.rdbuf(Output.rdbuf()); // Redirect all output from this process to file.
    std::cerr.rdbuf(Output.rdbuf());
  }
  ~mpi_handler() {
    MPI_Finalize();
  }
};
#else
#include <cstdlib>
#endif

int main(int argc, char** argv)
{
#ifdef USE_MPI
  // This should be the very first thing created, and the very last destroyed (as nearly as possible).
  static mpi_handler mpi(argc, argv);
#endif
  std::cout<<"Entered program."<<std::endl;
  static SafeStopwatch WholeProgramWatch("Whole program (sequential)");
  SafeStopwatch::tag WholeProgramTag = WholeProgramWatch.Start();
  std::string ProcessedFileName;
  std::string RawFileName;
  std::string OutFileName;
  std::string NoiseFileName;
  Long64_t StartEntry = 0;
  Long64_t NumEntries = 100;
  double Threshold = 10;

#ifdef USE_MPI
  std::ifstream OptionFile(std::string(argv[1]) + "/infile" + mpi.RankString);
  OptionFile >> ProcessedFileName
             >> RawFileName
             >> OutFileName
             >> NoiseFileName
             >> StartEntry
             >> NumEntries
             >> Threshold;
  // On NERSC, we always use xrootd.  Need the IP address of the MOM node.
  assert(argc == 3);
  std::string mom_ip = argv[2]; // Should also include port number used.
  ProcessedFileName = "root://cgd8d@" + mom_ip + "/" + ProcessedFileName + "?cachesz=300000000";
  RawFileName = "root://cgd8d@" + mom_ip + "/" + RawFileName;
#else
  assert(argc >= 4);
  ProcessedFileName = argv[1];
  RawFileName = argv[2];
  OutFileName = argv[3];
  NoiseFileName = argv[4];
  if(argc >= 6) StartEntry = std::atol(argv[5]);
  if(argc >= 7) NumEntries = std::atol(argv[6]);
  if(argc >= 8) Threshold = std::atof(argv[7]);
#endif

  std::cout<<"Input processed file: "<<ProcessedFileName<<std::endl;
  std::cout<<"Input raw file: "<<RawFileName<<std::endl;
  std::cout<<"Output file: "<<OutFileName<<std::endl;
  std::cout<<"Noise file: "<<NoiseFileName<<std::endl;
  std::cout<<"Starting at entry "<<StartEntry<<std::endl;
  std::cout<<"Handle "<<NumEntries<<" entries."<<std::endl;

  EXOTreeInputModule InputModule;
  std::cout<<"About to set filename."<<std::endl;
  InputModule.SetFilename(ProcessedFileName);
  std::cout<<"Successfully set filename."<<std::endl;
  TXNetFile WaveformFile(RawFileName.c_str());
  TTree* WaveformTree = dynamic_cast<TTree*>(WaveformFile.Get("tree"));

  EXOTreeOutputModule OutputModule;
  OutputModule.SetOutputFilename(OutFileName);
  OutputModule.Initialize();
  OutputModule.BeginOfRun(NULL); // OK, fine -- shortcut here, I assume input has only one run.

  EXORefitSignals RefitSig(InputModule, *WaveformTree, OutputModule);
  EXOCalibManager::GetCalibManager().SetMetadataAccessType("text");
  RefitSig.SetNoiseFilename(NoiseFileName);
  RefitSig.SetRThreshold(Threshold);
  RefitSig.fVerbose = true;
  RefitSig.Initialize();

#ifdef USE_THREADS
  std::cout<<"Using "<<NUM_THREADS<<" threads."<<std::endl;
#else
  std::cout<<"Sequential code."<<std::endl;
#endif
#ifdef USE_LOCKFREE
  std::cout<<"Using the boost::lockfree library."<<std::endl;
#endif

#ifdef USE_THREADS
  // Start a finish-up thread.
  boost::thread finishThread(&EXORefitSignals::FinishEventThread, &RefitSig);
#endif

  for(Long64_t entryNum = StartEntry;
      NumEntries == -1 or entryNum < StartEntry + NumEntries;
      entryNum++) {
    if(entryNum % 10 == 0) std::cout << "Grabbing entry " << entryNum << std::endl;
    if(entryNum % 100 == 0) {
#ifdef USE_THREADS
      while(RefitSig.GetFinishEventQueueLength() > 5000) {
        // If we're getting ahead of FinishEvent, wait a bit for it to catch up.
        // This is purely a memory concern of how many events we can save on the queue.
        std::cout<<"Stalling the computation threads for FinishEvent to catch up."<<std::endl;
        boost::this_thread::sleep_for(boost::chrono::seconds(5));
      }
#else
      RefitSig.FinishEventThread();
#endif
    }
#ifdef USE_THREADS
    static SafeStopwatch GetLockWatch("Get lock in main (sequential)");
    SafeStopwatch::tag GetLockTag = GetLockWatch.Start();
    boost::mutex::scoped_lock locLock(RootInterfaceMutex);
    GetLockWatch.Stop(GetLockTag);
#endif
    static SafeStopwatch InputModuleWatch("InputModule in main (sequential)");
    SafeStopwatch::tag InputModuleTag = InputModuleWatch.Start();
    EXOEventData* ED = InputModule.GetEvent(entryNum);
    InputModuleWatch.Stop(InputModuleTag);
    if(ED == NULL) break;
    static SafeStopwatch AcceptEventWatch("AcceptEvent (sequential)");
    SafeStopwatch::tag AcceptEventTag = AcceptEventWatch.Start();
    RefitSig.AcceptEvent(ED, entryNum
#ifdef USE_THREADS
      , locLock
#endif
    );
    AcceptEventWatch.Stop(AcceptEventTag);
  }

  static SafeStopwatch FlushEventsWatch("FlushEvents (sequential)");
  SafeStopwatch::tag FlushEventsTag = FlushEventsWatch.Start();
  RefitSig.FlushEvents();
  FlushEventsWatch.Stop(FlushEventsTag);
  static SafeStopwatch WaitForFinisherWatch("Waiting for events to be finished at the end (sequential)");
  SafeStopwatch::tag WaitForFinisherTag = WaitForFinisherWatch.Start();
#ifdef USE_THREADS
  RefitSig.fProcessingIsDone.store(true, boost::memory_order_seq_cst);
  finishThread.join();
#else
  RefitSig.FinishEventThread();
#endif
  WaitForFinisherWatch.Stop(WaitForFinisherTag);
  OutputModule.ShutDown();
  WholeProgramWatch.Stop(WholeProgramTag);
}
