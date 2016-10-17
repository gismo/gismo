/** @file gsJITCompiler.h

    @brief Provides declaration of JIT-compiler class.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): M. Moeller

    @note This class is based on the discussion on
    http://stackoverflow.com/questions/36040814/stdshared-ptr-and-dlopen-avoiding-undefined-behavior
*/
 
#pragma once

/*
#include <cstdlib>
#include <cstring>
#include <exception>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
*/

#ifdef __GNUC__ 
#include <cxxabi.h>
#endif

#include <gsIO/gsXml.h>

#if defined(_WIN32)
#include <windows.h>
#include <direct.h>
#define getcwd _getcwd 
#else //if defined(__APPLE__) || defined(__linux__) ||  defined(__unix)
#include <dlfcn.h>
#include <unistd.h>
#endif

namespace gismo {

/**
   @brief Struct definig a compiler configuration
   
   This class defines a compiler configuration that is used by the
   \ref gsJITCompiler class to perform just-in-time compilation.
*/
struct gsJITCompilerConfig
{
    /// Constructor (default)
    gsJITCompilerConfig()
    : cmd("missing"), flags("missing"), lang("missing")
    {
        char *env;
        env = getenv ("JIT_COMPILER_CMD");
        if(env!=NULL) cmd = env;
        
        env = getenv ("JIT_COMPILER_FLAGS");
        if(env!=NULL) flags = env;

        env = getenv ("JIT_COMPILER_LANG");
        if(env!=NULL) lang = env;

        out = "-o ";
    }

    /// Constructor (passing arguments as strings)
    gsJITCompilerConfig(const std::string& _cmd, const std::string& _flags,
                        const std::string& _lang, const std::string& _out)
    : cmd(_cmd), flags(_flags), lang(_lang), out(_out)
    {}
    
    /// Return compiler command
    virtual const std::string& getCmd() const { return cmd; }

    /// Return compiler flags
    virtual const std::string& getFlags() const { return flags; }

    /// Return compiler language
    virtual const std::string& getLang() const { return lang; }

    /// Return compiler output flag
    virtual const std::string& getOut() const { return out; }

    /// Set compiler command
    void setCmd(const std::string& _cmd)
    { this->cmd = _cmd; }

    /// Set compiler flags
    void setFlags(const std::string& _flags)
    { this->flags = _flags; }

    /// Set compiler language
    void setLang(const std::string& _lang)
    { this->lang = _lang; }
    
    /// Set compiler output flag
    void setOut(const std::string& _out)
    { this->out = _out; }

    /// Prints the object as a string
    std::ostream& print(std::ostream &os) const
    {
        os << "JITCompiler: " << cmd << ", flags: " << flags << ", language: " << lang <<"\n";
        return os;
    }

    /// Initialize to default GCC compiler
    static gsJITCompilerConfig gcc()
    {
        return gsJITCompilerConfig(
                                   "g++", // /usr/bin/g++
                                   "-fPIC -O3 -shared",
                                   "cxx", "-o ");
    }

    /// Initialize to default Clang compiler
    static gsJITCompilerConfig clang()
    {
        return gsJITCompilerConfig("clang++", // /usr/bin/clang++
                                   "-O3 -shared",
                                   "cxx", "-o ");
    }
    
    /// Initialize to default Intel compiler
    static gsJITCompilerConfig intel()
    {
        return gsJITCompilerConfig(
                                  "icpc",
                                  "-O3 -shared",
                                  "cxx", "-o ");
    }

    /// Initialize to default Clang compiler
    static gsJITCompilerConfig msvc()
    {
        return gsJITCompilerConfig("cl.exe",
                                   "/EHsc /Ox /LD",
                                   "cxx", "/Fe");//no space
    }

    /// Try to initialize compiler automatically based on the context
    static gsJITCompilerConfig guess()
    {        
#       if defined(_INTEL_COMPILER)
        return intel();
#       elif  defined(_MSC_VER)
        return msvc();
#       elif defined(__clang__)
        return clang();
#       elif defined(__GNUC__) || defined(__GNUG__)
        return gcc();
#       else
        GISMO_ERROR("Compiler not known");
#       endif
    }

protected:
    /// Members variables
    std::string cmd;
    std::string flags;
    std::string lang;
    std::string out;
};

/// Print (as string) operator to be used by all derived classes
std::ostream &operator<<(std::ostream &os, const gsJITCompilerConfig& c)
{return c.print(os); }

namespace internal
{

/** \brief Read a JITCompilerConfig from XML data
    \ingroup Core
*/
template<>
class gsXml< gsJITCompilerConfig >
{
private:
    gsXml() { }

public:
    GSXML_COMMON_FUNCTIONS(gsJITCompilerConfig)
    GSXML_GET_POINTER(gsJITCompilerConfig)
    static std::string tag () { return "JITCompilerConfig"; }
    static std::string type() { return ""; }

    static void get_into(gsXmlNode * node, gsJITCompilerConfig & result)
    {
        gsXmlAttribute * tmp = node->first_attribute("cmd");
        if (tmp!=NULL)
            result.setCmd(tmp->value());

        tmp = node->first_attribute("flags");
        if (tmp!=NULL)
            result.setFlags(tmp->value());

        tmp = node->first_attribute("lang");
        if (tmp!=NULL)
            result.setLang(tmp->value());

        tmp = node->first_attribute("out");
        if (tmp!=NULL)
            result.setOut(tmp->value());
    }

    static gsXmlNode * put (const gsJITCompilerConfig & obj, gsXmlTree & data)
    {
        // Make a new XML CompilerConfig node
        gsXmlNode * tmp = internal::makeNode("JITCompilerConfig", data);

        // Append the attributes
        tmp->append_attribute( makeAttribute("cmd"  , obj.getCmd()  , data) );
        tmp->append_attribute( makeAttribute("flags", obj.getFlags(), data) );
        tmp->append_attribute( makeAttribute("lang" , obj.getLang() , data) );
        tmp->append_attribute( makeAttribute("out"  , obj.getOut()  , data) );
        
        return tmp;
    }
};

} // namespace internal

namespace util
{

template<typename T> struct remove_pointer {typedef T type;};
template<typename T> struct remove_pointer<T*> {typedef typename remove_pointer<T>::type type;};


template<typename T>
struct type
{
public:
    static std::string name()
    {
#ifdef __GNUC__ 
        std::size_t len = 0;
        int status = 0;
        memory::unique<char>::ptr ptr(
            __cxxabiv1::__cxa_demangle( typeid(T).name(), NULL, &len, &status ) );
        return ptr.get();
#else
        return typeid(T).name();
#endif
    }
};

} // namespace util

/**
   @brief Class defining a dynamic library.

   This class stores a pointer to a dynamic library that can be
   compiled at runtime and provides extra functionality.
 */
struct gsDynamicLibrary
{
public:

    /// Constructor (using file name)
    gsDynamicLibrary(const char* filename, int flag)
    {
        gsDebug << "Loading dynamic library: " << filename << "\n";
        
#if defined(_WIN32)
        GISMO_UNUSED(flag);
        HMODULE dl = LoadLibrary(filename);
        if (!dl)
        {
            std::ostringstream err;
            err <<"LoadLibrary - error: " << GetLastError();
            throw std::runtime_error( err.str() );
        }
        handle.reset(dl, FreeLibrary);
#elif defined(__APPLE__) || defined(__linux__) || defined(__unix)        
        void * dl = ::dlopen(filename, flag);
        if (!dl)
            throw std::runtime_error( ::dlerror() );
        handle.reset(dl, ::dlclose);
#else
#error("Unsupported operating system")
#endif
    }

    /// Get symbol from dynamic library
    template<class T>
    T* getSymbol(const char* name) const
    {
        if (!handle)
            throw std::runtime_error("An error occured while accessing the dynamic library");
        
        T *symbol;
#if defined(_WIN32)
        *(void **)(&symbol) = GetProcAddress(handle.get(), name );
#elif defined(__APPLE__) || defined(__linux__) || defined(__unix)
        *(void **)(&symbol) = ::dlsym( handle.get(), name );
#endif
        if (!symbol)
            throw std::runtime_error("An error occured while getting symbolc from the dynamic library");
        
        return symbol;
    }
    
    /// Check if handle is assigned
    operator bool() const { return (bool)handle; }
    
private:

    /// Handle to dynamic library object
#if defined(_WIN32)
    memory::shared< util::remove_pointer<HMODULE>::type >::ptr handle;
#else //if defined(__APPLE__) || defined(__linux__) || defined(__unix)
    memory::shared<void>::ptr handle;
#endif
};

/**
   @brief Class defining a just-in-time compiler.

   This class compiles source code at runtime and links it to the
   running binary. This mechanism makes it possible to generate source
   code based on user-defined run-time parameters and still perform
   compile-time optimization.
*/
class gsJITCompiler
{    
public:
    /// Constructor (default)
    gsJITCompiler()
    : kernel(), config()
    { }

    /// Constructor (copy)
    gsJITCompiler(gsJITCompiler const& other)
    : config(other.config)
    {
        kernel << other.kernel.rdbuf();
    }

    /// Assignment operator (copy)
    gsJITCompiler& operator=(gsJITCompiler const& other)
    {
        kernel << other.kernel.rdbuf();
        config = other.config;
        return *this;
    }

    /*
    /// Constructor (move)
    gsJITCompiler(gsJITCompiler && other) :
    //kernel(std::move(other.kernel)),
    kernel(other.kernel.str()),
    config(std::move(other.config))
    {}


    /// Assignment operator (move)
    gsJITCompiler& operator=(gsJITCompiler && other)
    {
        kernel = std::move(other.kernel);
        config = std::move(other.config);
        return *this;
    }
    */
    
    /// Constructor (using compiler configuration)
    explicit gsJITCompiler(const gsJITCompilerConfig & _config)
    : kernel(), config(_config)
    {
        //gsInfo << config << "\n";
    }
    
    /// Input kernel source code from string
    gsJITCompiler & operator<<(const std::string & s)
    {
        kernel << s;
        return *this;
    }

    /// Input kernel source code from input stream
    gsJITCompiler & operator<<(std::istream & is)
    {
        kernel << is.rdbuf();
        return *this;
    }

    /// Compile kernel source code into dynamic library
    /// (determine filename from hash of kernel source code)
    gsDynamicLibrary build(bool force = false)
    {
        #if __cplusplus >= 201103L
        std::size_t h = std::hash<std::string>()(getKernel().str());
        return build(std::to_string(h), force);
        #else
        return build("JIT", true);
        #endif
    }
    
    /// Compile kernel source code into dynamic library
    /// (use given filename)
    gsDynamicLibrary build(const std::string &name, bool force = false)
    {
        // Prepare library name
        std::stringstream libName;
        memory::unique<char>::ptr path(getcwd(NULL,0));
#       if   defined(_WIN32)
        // DWORD psz = GetTempPath(MAX_PATH,// length of the buffer
        //                        path);    // buffer for path
        libName << path.get() << "\\." << name << ".dll";
#       elif defined(__APPLE__)
        libName << path.get() << "/.lib" << name << ".dylib";
#       elif defined(__unix)
        // char const * path = getenv("TMPDIR");
        // if (0 == folder) path = "/tmp";
        libName << path.get() << "/.lib" << name << ".so";
#       endif
        
        // Compile library (if required)
        std::ifstream libfile(libName.str().c_str());
        if(!libfile || force)
        {
            // Write kernel source code to file
            std::stringstream srcName;
#           ifdef _WIN32
            srcName<< path.get() << "\\." << name << "." << config.getLang();
            std::ofstream file(srcName.str().c_str());
            file   << "#define EXPORT extern \"C\" __declspec(dllexport)\n";
#           else
            srcName<< path.get() << "/." << name << "." << config.getLang();
            std::ofstream file(srcName.str().c_str());
            file   << "#define EXPORT extern \"C\"\n";
#           endif
            file << getKernel().str() <<"\n";
            file.close();
            
            // Compile kernel source code into library
            std::stringstream systemcall;
            /*
            systemcall << "del \""<< libName.str() <<"\"";
            (void)std::system(systemcall.str().c_str());
            systemcall.clear();systemcall.str("");
            //*/

            systemcall << config.getCmd()   << " "
                       << config.getFlags() << " \""
                       << srcName.str()     << "\" "
                       << config.getOut()   <<"\""<< libName.str()<<"\"";

            gsDebug << "Compiling dynamic library: " << systemcall.str() << "\n";
            if(std::system(systemcall.str().c_str()) != 0)
                throw std::runtime_error("An error occured while compiling the kernel source code");
        }

#ifdef _WIN32
        return gsDynamicLibrary( libName.str().c_str(), 0 );
#else
        return gsDynamicLibrary( libName.str().c_str(), RTLD_LAZY );
#endif
    }

    /// Clear kernel source code
    void clear()
    {
        kernel.clear();
        kernel.str(std::string());
    }

    /// Prints the object as a string
    std::ostream& print(std::ostream &os) const
    {
        os << kernel.str();
        return os;
    }

    /// Return kernel source code (as output stringstream)
    const std::ostringstream &getKernel() const
    {
        return kernel;
    }

    /// Return pointer to kernel source code
    std::ostringstream &getKernel()
    {
        return kernel;
    }
    
private:
    /// Kernel source code
    std::ostringstream kernel;
    
    /// Compiler configuration
    gsJITCompilerConfig config;
};

/// Print (as string) operator to be used by all derived classes
std::ostream &operator<<(std::ostream &os, const gsJITCompiler& c)
{ return c.print(os); }

namespace internal
{

/** \brief Read a JITCompiler from XML data
    \ingroup Core
*/
template<>
class gsXml< gsJITCompiler >
{
private:
    gsXml() { }

public:
    GSXML_COMMON_FUNCTIONS(gsJITCompiler)
    GSXML_GET_POINTER(gsJITCompiler)
    static std::string tag () { return "JITCompiler"; }
    static std::string type() { return ""; }

    static void get_into(gsXmlNode * node, gsJITCompiler & result)
    {
        result.getKernel() << node->value();
    }

    static gsXmlNode * put (const gsJITCompiler & obj, gsXmlTree & data)
    {
        // Make a new XML KnotVector node
        gsXmlNode * tmp = internal::makeNode("JITCompiler", obj.getKernel().str(), data);

        return tmp;
    }
};

} // namespace internal

} // namespace gismo
