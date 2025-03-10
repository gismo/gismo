/** @file gsStackWalker.cpp

    @brief Provides implementation of gsStackwalker.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#include <gsCore/gsDebug.h>
#include "gsStackWalker.h"

#include <vector>
#include <stdio.h>
#include <stdlib.h>

#ifdef _WIN32
  #include <imagehlp.h>
#else
  #include <err.h>
  #include <execinfo.h>
  #include <dlfcn.h>
  #include <cxxabi.h>
#endif

#define MAX_STACK_FRAMES 64

namespace gismo {

namespace internal {


#if defined(__GNUC__) && !defined(_WIN32)

bool gsExceptionHandler()
{
    /* register our signal handlers */
    
    struct sigaction sig_action = {};
    sig_action.sa_sigaction = gsExceptionHook;
    sigemptyset(&sig_action.sa_mask);
    
#ifdef __APPLE__
    sig_action.sa_flags = SA_SIGINFO;
#else
    sig_action.sa_flags = SA_SIGINFO | SA_ONSTACK;
#endif
    
    if (sigaction(SIGSEGV, &sig_action, NULL) != 0) { err(1, "sigaction"); }
    if (sigaction(SIGFPE,  &sig_action, NULL) != 0) { err(1, "sigaction"); }
    if (sigaction(SIGINT,  &sig_action, NULL) != 0) { err(1, "sigaction"); }
    if (sigaction(SIGILL,  &sig_action, NULL) != 0) { err(1, "sigaction"); }
    if (sigaction(SIGTERM, &sig_action, NULL) != 0) { err(1, "sigaction"); }
    if (sigaction(SIGABRT, &sig_action, NULL) != 0) { err(1, "sigaction"); }

    return true;
}

void gsStackWalker(void *context)
{
    (void)context;
    void * stack_traces[25];
    int nSize = backtrace(stack_traces, MAX_STACK_FRAMES);
    char ** symbols = backtrace_symbols(stack_traces, nSize);
    
    gsInfo<< "CALL STACK:\n";
    for (int i = 2; i < nSize - 2; i++)
    {
        printGccDemangled(symbols[i]);
    }

    free(symbols);
    
    gsInfo.flush();
}

void gsExceptionHook(int sig, siginfo_t *siginfo, void *context)
{
    switch(sig)
    {
    case SIGSEGV:
        fputs("Caught SIGSEGV: Segmentation Fault\n", stderr);
        break;
    case SIGINT:
        fputs("Caught SIGINT: Interactive attention signal, (usually ctrl+c)\n", stderr);
        break;
    case SIGFPE:
        switch(siginfo->si_code)
        {
        case FPE_INTDIV:
            fputs("Caught SIGFPE: (integer divide by zero)\n", stderr);
            break;
        case FPE_INTOVF:
            fputs("Caught SIGFPE: (integer overflow)\n", stderr);
            break;
        case FPE_FLTDIV:
            fputs("Caught SIGFPE: (floating-point divide by zero)\n", stderr);
            break;
        case FPE_FLTOVF:
            fputs("Caught SIGFPE: (floating-point overflow)\n", stderr);
            break;
        case FPE_FLTUND:
            fputs("Caught SIGFPE: (floating-point underflow)\n", stderr);
            break;
        case FPE_FLTRES:
            fputs("Caught SIGFPE: (floating-point inexact result)\n", stderr);
            break;
        case FPE_FLTINV:
            fputs("Caught SIGFPE: (floating-point invalid operation)\n", stderr);
            break;
        case FPE_FLTSUB:
            fputs("Caught SIGFPE: (subscript out of range)\n", stderr);
            break;
        default:
            fputs("Caught SIGFPE: Arithmetic Exception\n", stderr);
            break;
        }
    case SIGILL:
        switch(siginfo->si_code)
        {
        case ILL_ILLOPC:
            fputs("Caught SIGILL: (illegal opcode)\n", stderr);
            break;
        case ILL_ILLOPN:
            fputs("Caught SIGILL: (illegal operand)\n", stderr);
            break;
        case ILL_ILLADR:
            fputs("Caught SIGILL: (illegal addressing mode)\n", stderr);
            break;
        case ILL_ILLTRP:
            fputs("Caught SIGILL: (illegal trap)\n", stderr);
            break;
        case ILL_PRVOPC:
            fputs("Caught SIGILL: (privileged opcode)\n", stderr);
            break;
        case ILL_PRVREG:
            fputs("Caught SIGILL: (privileged register)\n", stderr);
            break;
        case ILL_COPROC:
            fputs("Caught SIGILL: (coprocessor error)\n", stderr);
            break;
        case ILL_BADSTK:
            fputs("Caught SIGILL: (internal stack error)\n", stderr);
            break;
        default:
            fputs("Caught SIGILL: Illegal Instruction\n", stderr);
            break;
        }
        break;
    case SIGTERM:
        fputs("Caught SIGTERM: a termination request was sent to the program\n", stderr);
        break;
    case SIGABRT:
        fputs("Caught SIGABRT: usually caused by an abort() or assert()\n", stderr);
        break;
    default:
        break;
    }
    gsStackWalker(context);
    _Exit(1);
}


void printGccDemangled(char * symbol)
{
	char *begin_name = 0, *begin_offset = 0, *end_offset = 0;
    
    // allocate string which will be filled with the demangled
    // function name
    size_t funcnamesize = 256;
    char* funcname = (char*)malloc(funcnamesize);
    
	// find parentheses and +address offset surrounding the mangled
	// name:
	for (char *p = symbol; *p; ++p)
	{
	    if (*p == '(')
            begin_name = p;
	    else if (*p == '+')
            begin_offset = p;
	    else if (*p == ')' && begin_offset) {
            end_offset = p;
            break;
	    }
	}
    
	if (begin_name && begin_offset && end_offset
	    && begin_name < begin_offset)
	{
	    *begin_name++ = '\0';
	    *begin_offset++ = '\0';
	    *end_offset = '\0';
        
	    // mangled name is now in [begin_name, begin_offset) and caller
	    // offset in [begin_offset, end_offset). now apply
	    // __cxa_demangle():
	    int status;
	    char* ret = abi::__cxa_demangle(begin_name,
                                        funcname, &funcnamesize, &status);
	    if (status == 0) 
        {
            funcname = ret; // use possibly realloc()-ed string
            gsInfo<< 
                //symbol <<"(+"<<begin_offset<<"):\n"<< 
                funcname<<"\n";
	    }
	    else {
            // demangling failed. Output function name as a C function with
            // no arguments.
            gsInfo<< symbol <<"(+"<<begin_offset<<")\n";
	    }
	}
	else
	{
	    // couldn't parse the line? print the whole line.
        gsInfo<< symbol <<"\n";
	}
    
    free(funcname);
}
    
// bool printAddr(void *addr) 
// {
//     Dl_info info;
//     if (dladdr(addr, &info) && info.dli_sname)
//     {
//         gsInfo << info.dli_fname << ": " << 
//             abi::__cxa_demangle(info.dli_sname, 0, 0, 0) << "\n";
//             return true;
//     }
//     return false;
// }

#endif // __GNUC__


#if defined(_WIN32)

    
bool gsExceptionHandler()
{
    ::SetUnhandledExceptionFilter(gsExceptionHook);

    // Use gsExceptionHook for all exceptions
    BOOL bRet = PreventSetUnhandledExceptionFilter();

    return bRet != 0;
}

void gsStackWalker(CONTEXT* context)
{
    HANDLE processHandle = GetCurrentProcess();
    HANDLE threadHandle  = GetCurrentThread();
    DWORD64  dwDisplacement = 0;

    std::vector<char> buffer(sizeof(SYMBOL_INFO) + MAX_SYM_NAME);
    PSYMBOL_INFO pSymbol = (PSYMBOL_INFO)&buffer.front();
    pSymbol->SizeOfStruct = sizeof(SYMBOL_INFO);
    pSymbol->MaxNameLen = MAX_SYM_NAME;

    SymInitialize(processHandle, 0, true);
    
    STACKFRAME frame = { 0 };
    
    /* setup initial stack frame */
    frame.AddrPC.Mode           = AddrModeFlat;
    frame.AddrStack.Mode        = AddrModeFlat;
    frame.AddrFrame.Mode        = AddrModeFlat;
#ifdef _M_IX86
    frame.AddrPC.Offset         = context->Eip;
    frame.AddrStack.Offset      = context->Esp;
    frame.AddrFrame.Offset      = context->Ebp;
#endif
#ifdef _M_AMD64
    frame.AddrPC.Offset         = context->Rip;
    frame.AddrStack.Offset      = context->Rsp;
    frame.AddrFrame.Offset      = context->Rbp;
#endif

    gsInfo<< "CALL STACK:\n";
    // StackWalk64
    while (StackWalk(IMAGE_FILE_MACHINE_I386 ,
                     processHandle,
                     threadHandle,
                     &frame,
                     context,
                     0,
                     SymFunctionTableAccess,
                     SymGetModuleBase,
                     0 ) )
    {
        // SymGetSymFromAddr64
        if (SymFromAddr(processHandle, frame.AddrPC.Offset, &dwDisplacement, pSymbol) )
        {
            gsInfo << pSymbol->Name <<"\n";
            //<< std::hex << " +0x" << dwDisplacement <<"\n";
        }
        else
        {
            gsInfo << " ???";
        }
    }

    SymCleanup( processHandle );
}

LONG WINAPI gsExceptionHook(EXCEPTION_POINTERS * ExceptionInfo)
{
    switch(ExceptionInfo->ExceptionRecord->ExceptionCode)
    {
    case EXCEPTION_ACCESS_VIOLATION:
        fputs("Error: EXCEPTION_ACCESS_VIOLATION\n", stderr);
        break;
    case EXCEPTION_ARRAY_BOUNDS_EXCEEDED:
        fputs("Error: EXCEPTION_ARRAY_BOUNDS_EXCEEDED\n", stderr);
        break;
    case EXCEPTION_BREAKPOINT:
        fputs("Error: EXCEPTION_BREAKPOINT\n", stderr);
        break;
    case EXCEPTION_DATATYPE_MISALIGNMENT:
        fputs("Error: EXCEPTION_DATATYPE_MISALIGNMENT\n", stderr);
        break;
    case EXCEPTION_FLT_DENORMAL_OPERAND:
        fputs("Error: EXCEPTION_FLT_DENORMAL_OPERAND\n", stderr);
        break;
    case EXCEPTION_FLT_DIVIDE_BY_ZERO:
        fputs("Error: EXCEPTION_FLT_DIVIDE_BY_ZERO\n", stderr);
        break;
    case EXCEPTION_FLT_INEXACT_RESULT:
        fputs("Error: EXCEPTION_FLT_INEXACT_RESULT\n", stderr);
        break;
    case EXCEPTION_FLT_INVALID_OPERATION:
        fputs("Error: EXCEPTION_FLT_INVALID_OPERATION\n", stderr);
        break;
    case EXCEPTION_FLT_OVERFLOW:
        fputs("Error: EXCEPTION_FLT_OVERFLOW\n", stderr);
        break;
    case EXCEPTION_FLT_STACK_CHECK:
        fputs("Error: EXCEPTION_FLT_STACK_CHECK\n", stderr);
        break;
    case EXCEPTION_FLT_UNDERFLOW:
        fputs("Error: EXCEPTION_FLT_UNDERFLOW\n", stderr);
        break;
    case EXCEPTION_ILLEGAL_INSTRUCTION:
        fputs("Error: EXCEPTION_ILLEGAL_INSTRUCTION\n", stderr);
        break;
    case EXCEPTION_IN_PAGE_ERROR:
        fputs("Error: EXCEPTION_IN_PAGE_ERROR\n", stderr);
        break;
    case EXCEPTION_INT_DIVIDE_BY_ZERO:
        fputs("Error: EXCEPTION_INT_DIVIDE_BY_ZERO\n", stderr);
        break;
    case EXCEPTION_INT_OVERFLOW:
        fputs("Error: EXCEPTION_INT_OVERFLOW\n", stderr);
        break;
    case EXCEPTION_INVALID_DISPOSITION:
        fputs("Error: EXCEPTION_INVALID_DISPOSITION\n", stderr);
        break;
    case EXCEPTION_NONCONTINUABLE_EXCEPTION:
        fputs("Error: EXCEPTION_NONCONTINUABLE_EXCEPTION\n", stderr);
        break;
    case EXCEPTION_PRIV_INSTRUCTION:
        fputs("Error: EXCEPTION_PRIV_INSTRUCTION\n", stderr);
        break;
    case EXCEPTION_SINGLE_STEP:
        fputs("Error: EXCEPTION_SINGLE_STEP\n", stderr);
        break;
    case EXCEPTION_STACK_OVERFLOW:
        fputs("Error: EXCEPTION_STACK_OVERFLOW\n", stderr);
        break;
    default:
        // fputs("Error: Unrecognized Exception\n", stderr);
        break;
    }
    fflush(stderr);
    
    /* If this is a stack overflow then we can't walk the stack, so just show
       where the error happened */
    if (EXCEPTION_STACK_OVERFLOW != ExceptionInfo->ExceptionRecord->ExceptionCode)
    {
        gsStackWalker(ExceptionInfo->ContextRecord);
    }
    else
    {
        gsInfo<< "EXCEPTION_STACK_OVERFLOW\n";
    }
    
    return EXCEPTION_EXECUTE_HANDLER;
}

#if defined _M_X64 || defined _M_IX86
LPTOP_LEVEL_EXCEPTION_FILTER WINAPI 
MyDummySetUnhandledExceptionFilter(
    LPTOP_LEVEL_EXCEPTION_FILTER )
{ return NULL; }
#else
  #error "This code works only for x86 and x64!"
#endif

BOOL PreventSetUnhandledExceptionFilter()
{
    HMODULE hKernel32 = LoadLibrary("kernel32.dll");
    if (hKernel32 == NULL) 
        return FALSE;
    void *pOrgEntry = GetProcAddress(hKernel32, 
                      "SetUnhandledExceptionFilter");
    if(pOrgEntry == NULL) 
        return FALSE;
    
    DWORD dwOldProtect = 0;
    SIZE_T jmpSize = 5;
#ifdef _M_X64
    jmpSize = 13;
#endif
    BOOL bProt = VirtualProtect(pOrgEntry, jmpSize, 
                                PAGE_EXECUTE_READWRITE, &dwOldProtect);
    BYTE newJump[20];
    void *pNewFunc = &MyDummySetUnhandledExceptionFilter;
#ifdef _M_IX86
    DWORD dwOrgEntryAddr = (DWORD) pOrgEntry;
    dwOrgEntryAddr += jmpSize; // add 5 for 5 op-codes for jmp rel32
    DWORD dwNewEntryAddr = (DWORD) pNewFunc;
    DWORD dwRelativeAddr = dwNewEntryAddr - dwOrgEntryAddr;
    // JMP rel32: Jump near, relative, displacement relative to next instruction.
    newJump[0] = 0xE9;  // JMP rel32
    memcpy(&newJump[1], &dwRelativeAddr, sizeof(pNewFunc));
#elif _M_X64
    // We must use R10 or R11, because these are "scratch" registers 
    // which need not to be preserved accross function calls
    // For more info see: Register Usage for x64 64-Bit
    // http://msdn.microsoft.com/en-us/library/ms794547.aspx
    // Thanks to Matthew Smith!!!
    newJump[0] = 0x49;  // MOV R11, ...
    newJump[1] = 0xBB;  // ...
    memcpy(&newJump[2], &pNewFunc, sizeof (pNewFunc));
    //pCur += sizeof (ULONG_PTR);
    newJump[10] = 0x41;  // JMP R11, ...
    newJump[11] = 0xFF;  // ...
    newJump[12] = 0xE3;  // ...
#endif
    SIZE_T bytesWritten;
    BOOL bRet = WriteProcessMemory(GetCurrentProcess(),
                pOrgEntry, newJump, jmpSize, &bytesWritten);
    
    if (bProt != FALSE)
    {
        DWORD dwBuf;
        VirtualProtect(pOrgEntry, jmpSize, dwOldProtect, &dwBuf);
    }
    return bRet;
}

#endif //  WIN32


} // end namespace internal

} // end namespace gismo


#undef MAX_STACK_FRAMES
