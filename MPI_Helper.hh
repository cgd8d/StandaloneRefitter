#ifndef MPI_Helper_hh
#define MPI_Helper_hh

#include <cassert>

/**
 * Call the MPI routine MPIFunc with arguments Args (surrounded by
 * parentheses). If the result is not MPI_SUCCESS, use
 * boost::throw_exception to throw an exception or abort, depending on
 * BOOST_NO_EXCEPTIONS.
 */
#define MPI_CHECK_RESULT( Func, Args )                                  \
 {                                                                      \
   int _check_result = Func Args;                                       \
   assert(_check_result == MPI_SUCCESS);                                \
 }

#endif
