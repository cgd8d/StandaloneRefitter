#include "SafeStopwatch.hh"
#include <cassert>

SafeStopwatch::SafeStopwatch(std::string id)
{
  // Constructor -- *NOT* guaranteed thread-safe.
  // (You must finish constructing a stopwatch before anyone uses it.)

#ifdef USE_THREADS
  // Not strictly necessary, but certainly preferable to ensure we're not hurting performance elsewhere.
  assert(fNanoseconds.is_lock_free());
  fNanoseconds.store(0, boost::memory_order_relaxed);
#else
  fNanoseconds = 0;
#endif
  fID = id;
}

SafeStopwatch::~SafeStopwatch()
{
  // Destructor -- print the time accumulated in this stopwatch.
  Print();
}

void SafeStopwatch::Print() const
{
  // Print the time of the stopwatch, terminated by an endline.
  // Any timers which have not been stopped will have no effect on the time printed.
  // Printing from multiple threads may yield garbage,
  // but technically this function is thread-safe.
  std::cout<<fID<<" has accumulated ";
#ifdef USE_THREADS
  std::cout<<fNanoseconds.load(boost::memory_order_relaxed);
#else
  std::cout<<fNanoseconds;
#endif
  std::cout<<" ns."<<std::endl;
}

SafeStopwatch::tag SafeStopwatch::Start() const
{
  // Get a "tag", which identifies the start time.
  // This is what will let the class compute a duration when "stop" is later called.
  // This function is thread-safe.
#ifdef USE_THREADS
  return tag();
#else
  timespec ts;
  assert(clock_gettime(CLOCK_MONOTONIC, &ts) == 0);
  return ts;
#endif
}

void SafeStopwatch::Stop(SafeStopwatch::tag& t)
{
  // Increment the stopwatch time by the duration between when tag was created and now.
  // This function is thread-safe.
#ifdef USE_THREADS
  t.stop();
  underlying_rep duration = t.elapsed().wall;
  fNanoseconds.fetch_add(duration, boost::memory_order_relaxed);
#else
  timespec now;
  assert(clock_gettime(CLOCK_MONOTONIC, &now) == 0);
  fNanoseconds += 1000000000LL*(long long int)(now.tv_sec - t.tv_sec);
  fNanoseconds += (now.tv_nsec - t.tv_nsec);
#endif
}
