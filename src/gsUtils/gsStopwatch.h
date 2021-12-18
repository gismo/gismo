/** @file gsStopwatch.h

    @brief Timing functions.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): C. Hofreither, J. Vogl
*/

#pragma once

#include <ostream>

// includes for wall clocks
#if __cplusplus >= 201103L || _MSC_VER >= 1600
#  include <chrono>
#elif defined(__linux__)
#  include <sys/time.h>
#elif defined(_MSC_VER) || defined(__MINGW32__)
#  include <sys/timeb.h>
#else
#  include <ctime>
#endif

// includes for CPU clocks
#if defined(__linux__)
#  include <sys/resource.h>
#else
#  include <ctime>
#endif

namespace gismo
{

inline std::ostream& formatTime(std::ostream& os, double sec)
{
    int flo = (int)(sec);
    sec -= flo;
    int hh = flo / 3600;
    flo = flo  % 3600;
    int mm = flo / 60;
    double ss = (flo % 60) + sec;
    if (hh > 0) os << hh << "h ";
    if (mm > 0) os << mm << "m ";
    std::streamsize prec = os.precision();
    os.precision(2);
    os << ss << " s";
    os.precision(prec);
    return os;
}


/** @brief A Stopwatch object can be used to measure execution time of
    code, algorithms, etc.
    
    This GenericStopwatch can be initialized in two time-taking modes, CPU time
    and real time (wall time).
    
    - For wall time use the  \ref gsStopwatch typedef.

    - For CPU time use the \ref gsCPUStopwatch typedef.

    \tparam Clock clock (eg. system specific) to use for measuring the time.
    
    \ingroup Utils
*/
template <typename Clock>
class gsGenericStopwatch
{
public:

    /// Declares a stop-watch
    gsGenericStopwatch() : m_value(0) { restart(); }

    /// Start taking the time
    void restart() { m_start = Clock::getTime(); }

    /// Return elapsed time in seconds.
    double stop() { return m_value = Clock::getTime() - m_start; }

    /// Returns the last recorded elapsed time
    double elapsed() const { return m_value; };

    friend std::ostream& operator<< (std::ostream& os, const gsGenericStopwatch& sw)
    { return formatTime(os, sw.m_value); }

#if __cplusplus >= 201103L || _MSC_VER >= 1600   // C++11 //
    gsGenericStopwatch (const gsGenericStopwatch&) = delete;
    gsGenericStopwatch& operator=(const gsGenericStopwatch&) = delete;
#endif

private:
    double m_start;
    double m_value;
}; // class gsGenericStopwatch


/*
 * SYSTEM-SPECIFIC WALL CLOCKS
 */

#if __cplusplus >= 201103L || _MSC_VER >= 1600   // C++11 //
// highest-resolution wall clock time
struct CXX11WallClock
{
    static double getTime()
    {
        return ((std::chrono::duration<double>)std::chrono::high_resolution_clock::now().time_since_epoch()).count();
    }
};

/// @brief A stop-watch measuring real (wall) time
typedef gsGenericStopwatch<CXX11WallClock> gsStopwatch;

#elif defined(__linux__) // || defined(TARGET_OS_MAC)               // LINUX //

// higher resolution wall clock time
struct LinuxWallClock
{
    static double getTime()
    { timeval tv; gettimeofday(&tv, 0); return tv.tv_sec + 1e-6*tv.tv_usec; }
};

/// @brief A stop-watch measuring real (wall) time
typedef gsGenericStopwatch<LinuxWallClock> gsStopwatch;

#elif defined(_MSC_VER) || defined(__MINGW32__)                    // WINDOWS //

struct WindowsWallClock
{
    static double getTime()
    {
        _timeb tb;
       #ifdef __MINGW32__
        _ftime( &tb );
       #else
        _ftime_s( &tb );
       #endif
        return (double)tb.time + tb.millitm / 1000.0;
    }
};

/// @brief A stop-watch measuring real (wall) time
typedef gsGenericStopwatch<WindowsWallClock> gsStopwatch;

#else                                                             // PORTABLE //

// simple, low-resolution (only integer seconds) wall clock time
struct WallClock
{
    static double getTime()
    { time_t secs; time (&secs); return secs; }
};

/// @brief A stop-watch measuring real (wall) time
typedef gsGenericStopwatch<WallClock> gsStopwatch;

#endif


/*
 * SYSTEM-SPECIFIC CPU CLOCKS
 */

#if defined(__linux__) // || defined(TARGET_OS_MAC)                     // LINUX //

// A non-portable, more expensive, but higher resolution CPU clock
// which also does not suffer from the short wrap-around time of
// CPUClock.
struct LinuxCPUClock
{
    static double getTime()
    {
        rusage ru;
        getrusage(RUSAGE_SELF, &ru);
        return ru.ru_utime.tv_sec + 1.0e-6*ru.ru_utime.tv_usec;
    }
};

/// @brief A stop-watch measuring CPU time
typedef gsGenericStopwatch<LinuxCPUClock> gsCPUStopwatch;

#else                                                             // PORTABLE //

// A portable, but not very accurate CPU clock.  Note that on a
// typical Linux system, clock_t is just 4 bytes long and this will
// wrap around after about 36 min.
struct CPUClock
{
    static double getTime() { return (double) clock() / (double) CLOCKS_PER_SEC; }
};

/// @brief A stop-watch measuring CPU time
typedef gsGenericStopwatch<CPUClock> gsCPUStopwatch;

#endif


} // namespace gismo
