/*
Create a thread-safe stopwatch -- multiple threads can use this watch concurrently without any sort of locking.
The resulting time will be the sum of the times of all of the component stopwatches.
If USE_THREADS is undefined, this stopwatch will use POSIX calls instead of BOOST,
and no thread-safety properties are guaranteed.

Proper use:
void f() {
  static SafeStopwatch f_Watch("f()");
  SafeStopwatch::tag f_Tag = f_Watch.Start();
  ...
  f_Watch.Stop(f_Tag);
}

Upon program exit, the watch will print out its total accumulated time as it is destroyed.
Note that static local initialization is guaranteed thread-safe in C++11; it has also been
thread-safe in gcc for a long time.  Therefore, it is fine for our purposes.
*/
#ifndef SafeStopwatch_hh
#define SafeStopwatch_hh

#ifdef USE_THREADS
#include <boost/atomic.hpp>
#include <boost/timer/timer.hpp>
#else
#include <time.h>
#endif
#include <iostream>
#include <string>

class SafeStopwatch
{
 public:
#ifdef USE_THREADS
  typedef boost::timer::cpu_timer tag;
#else
  typedef timespec tag;
#endif
  SafeStopwatch(std::string id);
  ~SafeStopwatch();
  tag Start() const;
  void Stop(tag& t);
  void Print() const;
 private:
#ifdef USE_THREADS
  typedef boost::timer::nanosecond_type underlying_rep;
  typedef boost::atomic<underlying_rep> rep;
#else
  typedef long long int rep;
#endif
  rep fNanoseconds;
  std::string fID;
};
#endif
