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
#include "EventFinisher.hh"
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
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
// To control order of destruction -- we'd like for all function-static objects to be destroyed
// before MPI_Finalize is called, so their destructors have a chance to execute side-effects.
// Also redirect output properly at this stage.
struct mpi_handler
{
  boost::mpi::environment env;
  boost::mpi::communicator comm;
  int rank;
  int numa_rank;
  std::string RankString;
  std::ofstream Output;
  mpi_handler(int argc, char** argv) :
  env(argc, argv)
  {
    assert(argc >= 2);
    rank = comm.rank();
#ifdef USE_PROCESSES
    numa_rank = rank/2;
#else
    numa_rank = rank;
#endif

    // Make sure every rank has its own output file.
    std::ostringstream OutFileName;
    OutFileName << argv[1] << "/outfile" << std::setw(4) << std::setfill('0') << rank << ".txt";
    Output.open(OutFileName.str().c_str());
    std::cout.rdbuf(Output.rdbuf()); // Redirect all output from this process to file.
    std::cerr.rdbuf(Output.rdbuf());

    // Make sure each numa_node reads in the proper input file.
    std::ostringstream ProcessRankString;
    ProcessRankString << std::setw(4) << std::setfill('0') << numa_rank << ".txt";
    RankString = ProcessRankString.str();
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
  std::ifstream OptionFile((std::string(argv[1]) + "/infile" + mpi.RankString).c_str());
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

#if defined(USE_MPI) && defined(USE_PROCESSES)
  if(mpi.rank % 2 == 1) {
    // This is an io process.
    EventFinisher& finisher = EventFinisher::Get(InputModule, RawFileName, OutFileName);
    finisher.Run();
    WholeProgramWatch.Stop(WholeProgramTag);
    return 0;
  }
  else {
    // This is a compute process.
#else
  // We need both io and compute objects.
  EventFinisher& finisher = EventFinisher::Get(InputModule, RawFileName, OutFileName);
#endif

  EXORefitSignals RefitSig
#ifndef USE_PROCESSES
    (finisher)
#endif
  ;
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

#if defined(USE_THREADS) && !defined(USE_PROCESSES)
  // Start a finish-up thread.
  boost::thread finishThread(&EventFinisher::Run, &finisher);
#endif

#ifdef USE_PROCESSES
  // Put out a continuing non-blocking request for messages from the io process.
  boost::mpi::request req = mpi.comm.irecv(mpi.rank+1, 1);
#endif

  for(Long64_t entryNum = StartEntry;
      NumEntries == -1 or entryNum < StartEntry + NumEntries;
      entryNum++) {
    if(entryNum % 10 == 0) std::cout << "Grabbing entry " << entryNum << std::endl;
#ifndef USE_PROCESSES
    if(entryNum % 100 == 0) {
#ifdef USE_THREADS
      while(finisher.GetFinishEventQueueLength() > 5000) {
        // If we're getting ahead of FinishEvent, wait a bit for it to catch up.
        // This is purely a memory concern of how many events we can save on the queue.
        std::cout<<"Stalling the computation threads for FinishEvent to catch up."<<std::endl;
        boost::this_thread::sleep_for(boost::chrono::seconds(5));
      }
#else
      finisher.Run();
#endif
    }
#endif
#ifdef USE_PROCESSES
    if(req.test()) { // We've been asked to pause.
      std::cout<<"Stalling the computation threads for FinishEvent to catch up."<<std::endl;
      mpi.comm.recv(mpi.rank+1, 0); // Wait until we get the go-ahead.
      req = mpi.comm.irecv(mpi.rank+1, 1); // Back to a passive wait for problems.
    }
#endif

#if defined(USE_THREADS) && !defined(USE_PROCESSES)
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
#if defined(USE_THREADS) && !defined(USE_PROCESSES)
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
#if defined(USE_THREADS) && !defined(USE_PROCESSES)
  finisher.SetProcessingIsFinished();
  finishThread.join();
#elif defined(USE_PROCESSES)
  mpi.comm.send(mpi.rank+1, 1, EventHandler());
  // Send a message with a non-zero tag -- the payload is unimportant.
#else
  finisher.Run();
#endif
  WaitForFinisherWatch.Stop(WaitForFinisherTag);
#ifdef USE_PROCESSES
  }
#endif
  WholeProgramWatch.Stop(WholeProgramTag);
}
