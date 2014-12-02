/** @file gsDebug.h
 * This file contains the debugging and messaging system of G+Smo. 
 *
 * @{
 */

#pragma once

#include <iomanip>
#include <stdexcept>
#include <typeinfo>

// See also about memory leak detection:
// http://msdn.microsoft.com/en-us/library/e5ewb1h3%28v=vs.90%29.aspx
#if _MSC_VER //>= 1400
#include <crtdbg.h>
//#include <errno.h>
#endif

#include <gsCore/gsLogging.h>

#ifdef GISMO_EXTRA_DEBUG
  #include <gsCore/gsStackWalker.h>
#endif


/** 
 *  Used for optional inclusion of .hpp header files in the .h files.
 *  Allows to drop dependencies on the .hpp files when using
 *  GISMO_HEADERS_ONLY for compiling a library instance.  
 *  When compiling as a pure template library CMake will detect
 *  dependency on .hpp files.
 */
#define GISMO_HPP_HEADER(x) #x

/**  
 *  Runtine assertions which display a message  
 *
 */

#ifndef NDEBUG

#   define GISMO_ASSERT(condition, message) \
    do {                                                                \
        if (! (condition) ) {                                           \
            gsDebug << "Assertion `" #condition "` failed in " << __FILE__ \
                    << " line " << __LINE__ << "\nMESSAGE :" << message << "\n"; \
                throw std::runtime_error("GISMO_ASSERT failure");       \
        }                                                               \
    } while (false)
#else
#   define GISMO_ASSERT(condition, message)
#endif

/**  
 *  Runtine check and display error message. This command is the same as
 *  GISMO_ASSERT but it is executed in release builds as well.
 *
 */

#   define GISMO_ENSURE(condition, message) \
        if (! (condition) ) {                                           \
            gsWarn  << "Condition `" #condition "` failed in " << __FILE__ \
                    << " line " << __LINE__ << ". MESSAGE:" << message << "\n"; \
            throw std::runtime_error("GISMO_ENSURE failure"); \
        }

/**  
 *  Denote a variable as unused, used to silence warnings in release
 *  mode builds.
 *
 */

#   define GISMO_UNUSED(x)  static_cast<void>(x)

/**  
 *  Runtine error message
 *
 */

#   define GISMO_ERROR(message)                 \
    {                                                                \
        gsInfo  << "Error in " << __FILE__                      \
                << " line " << __LINE__ << ". MESSAGE: " << message << "\n"; \
        throw std::runtime_error("GISMO_ERROR");	\
    }

/**  
 *  Runtine "no implementation" error happens when the user calls a
 *  virtual member function without a default implementetion.
 */
 
// TO DO: for GCC __PRETTY_FUNC__ is better
# define GISMO_NO_IMPLEMENTATION \
  gsInfo<<"Virtual member function \""<< __FUNCTION__ << "(..)\" declared in " \
  << __FILE__ <<" has not been implemented ("<<typeid(*this).name()<<").\n"; \
  throw std::runtime_error("GISMO_NO_IMPLEMENTATION");


/*
  Disable debug/abort popup windows on MS Windows
  
  See http://msdn.microsoft.com/en-us/library/1y71x448.aspx

  You might also need to disable "error reporting" on your windows
  system for popup-free runs.
*/
#if _MSC_VER //>= 1400 
static const int    gismo_CrtSetReportMode = _CrtSetReportMode( 
    _CRT_ASSERT, _CRTDBG_MODE_FILE   );
static const _HFILE gismo_CrtSetReportFile = _CrtSetReportFile( 
    _CRT_ASSERT, _CRTDBG_FILE_STDERR );
static const int  gismo_set_abort_behavior = _set_abort_behavior( 
    0x0, _WRITE_ABORT_MSG | _CALL_REPORTFAULT);
#endif

/*
  Disable some Warnings
*/

#ifdef _MSC_VER
// 4100 - unreferenced formal parameter
// 4101 - unreferenced local variable
// 4127 - conditional expression is constant (triggered by assertion macros)
// 4146 - unary minus operator applied to unsigned type, result still unsigned
// 4181 - qualifier applied to reference type ignored
// 4211 - nonstandard extension used : redefined extern to static
// 4244 - 'argument' : conversion from 'type1' to 'type2', possible loss of data
// 4251 - needs to have dll-interface to be used by clients of class
// 4273 - QtAlignedMalloc, inconsistent DLL linkage
// 4275 - non dll-interface base 
// 4324 - structure was padded due to declspec(align())
// 4428 - universal-character-name encountered in source
// 4503 - decorated name length exceeded
// 4505 - unreferenced local function has been removed
// 4512 - assignment operator could not be generated
// 4522 - 'class' : multiple assignment operators specified
// 4566 - character represented by universal-character-name cannot be represented in the current code page
// 4661 - no definition available
// 4700 - uninitialized local variable 'xyz' used
// 4702 - unreachable code
// 4714 - function marked as __forceinline not inlined
// 4717 - 'function' : recursive on all control paths, function will cause runtime stack overflow
// 4789 - destination of memory copy is too small (for Eigen)
// 4996 - 'sprintf': This function or variable may be unsafe. Consider using sprintf_s instead.
  #pragma warning( push )
  #pragma warning( disable : 4100 4127 4146 4251 4428 4275 4503 4505 4512 4566 4661 4714 4789 4996 )

#elif defined __INTEL_COMPILER
// 2196 - routine is both "inline" and "noinline" ("noinline" assumed)
//        ICC 12 generates this warning even without any inline keyword, when defining class methods 'inline' i.e. inside of class body
//        typedef that may be a reference type.
// 279  - controlling expression is constant
//        ICC 12 generates this warning on assert(constant_expression_depending_on_template_params) and frankly this is a legitimate use case.
// 161  - unrecognized pragma
// 175  - subscript out of range
//        to avoid warnings on #pragma GCC diagnostic          
  #pragma warning push
  #pragma warning disable 2196 279 161 175

#elif defined __clang__
// -Wconstant-logical-operand - warning: use of logical && with constant operand; switch to bitwise & or remove constant
// -Wbind-to-temporary-copy - warning: Warn about an unusable copy constructor when binding a reference to a temporary
    #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wconstant-logical-operand"
  #pragma clang diagnostic ignored "-Wbind-to-temporary-copy"

#elif defined __GNUC__ // major version
// typedef locally defined but not used [-Wunused-local-typedefs]
#if __GNUC_MINOR__ > 7
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#endif

#endif


/*
 *  Compile-time assertions: 
 * 
 * - in GISMO_STATIC_ASSERT(CONDITION,MSG) the parameter CONDITION must be a compile time boolean
 *    expression, and MSG an enum listed in struct internal::static_assertion<true>
 *
 *  - define GISMO_NO_STATIC_ASSERT to disable them (and save compilation time)
 *    in that case, the static assertion is converted to the following runtime assert:
 *      gismo_assert(CONDITION && "MSG")
 *
 *  - currently GISMO_STATIC_ASSERT can only be used in function scope
 *
 */

#ifndef GISMO_NO_STATIC_ASSERT

  #if defined(__GXX_EXPERIMENTAL_CXX0X__) || (defined(_MSC_VER) && (_MSC_VER >= 1600))

    // if native static_assert is enabled, let's use it
    #define GISMO_STATIC_ASSERT(X,MSG) static_assert(X,#MSG);

  #else // not CXX0X

    namespace gismo {

    namespace internal {

    template<bool condition>
    struct static_assertion {};

    template<>
    struct static_assertion<true>
    {
      enum {
          YOU_CALLED_AN_INVALID_CONSTRUCTOR,
          YOU_CALLED_A_FUNCTION_THAT_IS_NOT_IMPLEMENTED,
          OBJECT_ALLOCATED_ON_STACK_IS_TOO_BIG
      };
    };

// Check:
// #undef GISMO_NO_IMPLEMENTATION
//       template<class T>
//       struct gismo_false { };
// # define GISMO_NO_IMPLEMENTATION  
//       {gismo::internal::gismo_false<T>::YOU_CALLED_A_FUNCTION_THAT_IS_NOT_IMPLEMENTED;};


    } // end namespace internal

    } // end namespace gismo

    // Specialized implementation for MSVC to avoid "conditional
    // expression is constant" warnings.  This implementation doesn't
    // appear to work under GCC, hence the multiple implementations.
    #ifdef _MSC_VER

      #define GISMO_STATIC_ASSERT(CONDITION,MSG) \
        {gismo::internal::static_assertion<bool(CONDITION)>::MSG;}

    #else

      #define GISMO_STATIC_ASSERT(CONDITION,MSG) \
        if (gismo::internal::static_assertion<bool(CONDITION)>::MSG) {}

    #endif

  #endif // not CXX0X

#else // GISMO_NO_STATIC_ASSERT

#define GISMO_STATIC_ASSERT(CONDITION,MSG) gismo_assert((CONDITION) && #MSG);

#endif // GISMO_NO_STATIC_ASSERT

// static assertion failing if a function is that is not implemeted is called
#define GISMO_STATIC_NO_IMPLEMENTATION \
  GISMO_STATIC_ASSERT(false                      , \
                      YOU_CALLED_A_FUNCTION_THAT_IS_NOT_IMPLEMENTED)

// static assertion failing if the type \a TYPE is not fixed-size
#define GISMO_STATIC_ASSERT_FIXED_SIZE(TYPE) \
  GISMO_STATIC_ASSERT(TYPE::SizeAtCompileTime!=Gismo::Dynamic, \
                      YOU_CALLED_A_FIXED_SIZE_METHOD_ON_A_DYNAMIC_SIZE_MATRIX_OR_VECTOR)

// static assertion failing if the type \a TYPE is not dynamic-size
#define GISMO_STATIC_ASSERT_DYNAMIC_SIZE(TYPE) \
  GISMO_STATIC_ASSERT(TYPE::SizeAtCompileTime==Gismo::Dynamic, \
                      YOU_CALLED_A_DYNAMIC_SIZE_METHOD_ON_A_FIXED_SIZE_MATRIX_OR_VECTOR)

  #define GISMO_STATIC_ASSERT_NON_INTEGER(TYPE) \
    GISMO_STATIC_ASSERT(!NumTraits<TYPE>::IsInteger, THIS_FUNCTION_IS_NOT_FOR_INTEGER_NUMERIC_TYPES)

// static assertion failing if it is guaranteed at compile-time that
// the two matrix expression types have different sizes
#define GISMO_STATIC_ASSERT_SAME_MATRIX_SIZE(TYPE0,TYPE1) \
  GISMO_STATIC_ASSERT( \
     GISMO_PREDICATE_SAME_MATRIX_SIZE(TYPE0,TYPE1),\
    YOU_MIXED_MATRICES_OF_DIFFERENT_SIZES)

// Next line closes the DEBUG_GROUP of Doxygen
/** @} */



#ifdef __GNUC__
#define GS_DEPRECATED __attribute__((deprecated))
#elif defined(_MSC_VER)
#define GS_DEPRECATED __declspec(deprecated)
#else
#pragma message("WARNING: you will not be warned about deprecated functions with this compiler")
#define GS_DEPRECATED
#endif
