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

#include <gsIO/gsXml.h>
#include <gsIO/gsFileManager.h>

#if defined(_WIN32)
#include <windows.h>
#else
#include <dlfcn.h>
#endif

#include <gsCore/gsMemory.h>

namespace gismo {

/**
   @brief Supported languages
*/
struct gsJITLang
{
    enum {
        C       = 0, ///< C
        CXX     = 1, ///< C++
        CUDA    = 2, ///< Cuda
        Fortran = 3  ///< Fortran
    };
};

/**
   @brief Struct definig a compiler configuration
   
   This class defines a compiler configuration that is used by the
   \ref gsJITCompiler class to perform just-in-time compilation.
*/
struct gsJITCompilerConfig
{
    /// Constructor (default)
    gsJITCompilerConfig()
    : cmd("missing"), flags("missing"), lang("missing"), out("-o "), temp(detectTemp())
    {
        char *env;
        env = getenv ("JIT_COMPILER_CMD");
        if(env!=NULL) cmd = env;
        
        env = getenv ("JIT_COMPILER_FLAGS");
        if(env!=NULL) flags = env;

        env = getenv ("JIT_COMPILER_LANG");
        if(env!=NULL) lang = env;

        env = getenv ("JIT_COMPILER_TEMP");
        if(env!=NULL) temp = env;
    }

    virtual ~gsJITCompilerConfig() { }
    
    /// Constructor (passing arguments as strings)
    gsJITCompilerConfig(const std::string& cmd,
                        const std::string& flags,
                        const std::string& lang,
                        const std::string& out,
                        const std::string& temp = detectTemp())
    : cmd(cmd), flags(flags), lang(lang), out(out), temp(temp)
    {}

    void swap(gsJITCompilerConfig & other)
    {
        std::swap(cmd  , other.cmd  );
        std::swap(flags, other.flags);
        std::swap(lang , other.lang );
        std::swap(out  , other.out  );
        std::swap(temp , other.temp );
    }
            
#   if __cplusplus >= 201103L || _MSC_VER >= 1600

    /// Constructor (copy)
    gsJITCompilerConfig(gsJITCompilerConfig const& other)
    { operator=(other);}

    /// Assignment operator
    gsJITCompilerConfig& operator=(const gsJITCompilerConfig & other)
    {
        cmd   = other.cmd;
        flags = other.flags;
        lang  = other.lang;
        out   = other.out;
        temp  = other.temp;
        return *this;
    }

    /// Constructor (move)
    gsJITCompilerConfig(gsJITCompilerConfig && other)
    : cmd(std::move(other.cmd)), flags(std::move(other.flags)),
      lang(std::move(other.lang)), out(std::move(other.out)),
      temp(std::move(other.temp))
    {}
    
    /// Assignment operator (move)
    gsJITCompilerConfig& operator=(gsJITCompilerConfig && other)
    {
        cmd   = std::move(other.cmd);
        flags = std::move(other.flags);
        lang  = std::move(other.lang);
        out   = std::move(other.out);
        temp  = std::move(other.temp);
        return *this;
    }
#else
    /// Assignment operator
    gsJITCompilerConfig& operator=(gsJITCompilerConfig other)
    {
        this->swap(other);
        return *this;
    }        
#   endif
    
    /// Return compiler command
    virtual const std::string& getCmd() const { return cmd; }

    /// Return compiler flags
    virtual const std::string& getFlags() const { return flags; }

    /// Return compiler language
    virtual const std::string& getLang() const { return lang; }

    /// Return compiler output flag
    virtual const std::string& getOut() const { return out; }
    
    /// Return compiler temporal directory
    virtual const std::string& getTemp() const { return temp; }

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
    
    /// Set compiler temporal directory
    void setTemp(const std::string& _temp)
    { this->temp = _temp; }

    /// Prints the object as a string
    std::ostream& print(std::ostream &os) const
    {
        os << "JIT Compiler.\n"
           << "  cmd:                " << cmd << "\n"
           << "  flags:              " << flags << "\n"
           << "  language:           " << lang << "\n"
           << "  output flag:        " << out << "\n"
           << "  temporal directory: " << temp << "\n";
        
        return os;
    }

    /// Reads compiler configuration from XML file by language
    void load(const std::string filename,
              const int _lang = gsJITLang::CXX)
    {
        GISMO_ENSURE(_lang >= gsJITLang::C && _lang <= gsJITLang::Fortran,
            "Error: Invalid compiler language.");

        gsFileData<real_t> f(filename);        
        gsJITCompilerConfig * cc = f.getId<gsJITCompilerConfig>(_lang).release();

        std::swap(*cc, *this);
        if (this->temp.empty()) this->temp=detectTemp();
        delete cc;
    }

    /// Reads compiler configuration from XML file by ID
    void load_id(const std::string filename,
                 const int id)
    {
        gsFileData<real_t> f(filename);        
        gsJITCompilerConfig * cc = f.getId<gsJITCompilerConfig>(id).release();

        std::swap(*cc, *this);
        if (this->temp.empty()) this->temp=detectTemp();
        delete cc;
    }
    
    /// Initialize to default Clang compiler
    static gsJITCompilerConfig clang(const int lang = gsJITLang::CXX)
    {
        switch(lang)
        {
        case (gsJITLang::C) :
            return gsJITCompilerConfig("clang",
                                       "-O3 -shared",
                                       "c",
                                       "-o ");
            break;
        case (gsJITLang::CXX) :
            return gsJITCompilerConfig("clang++",
                                       "-O3 -shared",
                                       "cxx",
                                       "-o ");
            break;
        case (gsJITLang::Fortran) :
            GISMO_ERROR("Error : Clang does not provide any Fortran compiler.");
            break;
        default :
            GISMO_ERROR("Error : Invalid compiler language.");
        }
    }
    
    /// Initialize to default GCC compiler
    static gsJITCompilerConfig gcc(const int lang = gsJITLang::CXX)
    {
        switch(lang)
        {
        case (gsJITLang::C) :
            return gsJITCompilerConfig("gcc",
                                       "-fPIC -O3 -shared",
                                       "c",
                                       "-o ");
            break;
        case (gsJITLang::CXX) :
            return gsJITCompilerConfig("g++",
                                       "-fPIC -O3 -shared",
                                       "cxx",
                                       "-o ");
            break;
        case (gsJITLang::Fortran) :
            return gsJITCompilerConfig("gfortran",
                                       "-fPIC -O3 -shared",
                                       "F90",
                                       "-o ");
            break;
        default :
            GISMO_ERROR("Error : Invalid compiler language.");
        }
    }
    
    /// Initialize to default Intel compiler
    static gsJITCompilerConfig intel(const int lang = gsJITLang::CXX)
    {
        switch(lang)
        {
        case (gsJITLang::C) :
#if         defined(_WIN32)
            return gsJITCompilerConfig("icl",
                                       "/O3 /dll",
                                       "c",
                                       "/Fo");//no space
#           else
            return gsJITCompilerConfig("icc",
                                       "-O3 -shared",
                                       "c",
                                       "-o ");
#           endif
            break;
        case (gsJITLang::CXX) :
#if         defined(_WIN32)
            return gsJITCompilerConfig("icl",
                                       "/O3 /dll",
                                       "c",
                                       "/Fo");//no space
#           else
            return gsJITCompilerConfig("icpc",
                                       "-O3 -shared",
                                       "cxx",
                                       "-o ");
#           endif
            break;
        case (gsJITLang::Fortran) :
#if         defined(_WIN32)
            return gsJITCompilerConfig("ifort",
                                       "/O3 /dll",
                                       "F90",
                                       "/Fo");//no space
#           else
            return gsJITCompilerConfig("ifort",
                                       "-O3 -shared",
                                       "F90",
                                       "-o ");
#           endif
            break;
        default :
            GISMO_ERROR("Error : Invalid compiler language.");
        }
    }

    /// Initialize to default Microsoft Visual Studio compiler
    static gsJITCompilerConfig msvc(const int lang = gsJITLang::CXX)
    {
        switch(lang)
        {
        case (gsJITLang::C) :
            return gsJITCompilerConfig("cl.exe",
                                       "/EHsc /Ox /LD",
                                       "c",
                                       "/Fe");//no space
            break;
        case (gsJITLang::CXX) :
            return gsJITCompilerConfig("cl.exe",
                                       "/EHsc /Ox /LD",
                                       "cxx",
                                       "/Fe");//no space
            break;
        default :
            GISMO_ERROR("Error : Invalid compiler language.");
        }   
    }

    /// Initialize to default NVIDIA nvcc compiler
    static gsJITCompilerConfig nvcc(const int lang = gsJITLang::CUDA)
    {
        switch(lang)
        {
        case (gsJITLang::CUDA) :
            return gsJITCompilerConfig("nvcc",
                                       "-O3 -shared",
                                       "cu",
                                       "-o ");
            break;
        default :
            GISMO_ERROR("Error : Invalid compiler language.");
        }
    }

    /// Initialize to default PGI compiler
    static gsJITCompilerConfig pgi(const int lang = gsJITLang::CXX)
    {
        switch(lang)
        {
        case (gsJITLang::C) :
            return gsJITCompilerConfig("pgcc",
                                       "-O3 -shared",
                                       "c",
                                       "-o ");
            break;
        case (gsJITLang::CXX) :
            return gsJITCompilerConfig("pgc++",
                                       "-O3 -shared",
                                       "cxx",
                                       "-o ");
            break;
        case (gsJITLang::Fortran) :
            return gsJITCompilerConfig("pgf90",
                                       "-O3 -shared",
                                       "F90",
                                       "-o ");
            break;
        default :
            GISMO_ERROR("Error : Invalid compiler language.");
        }
    }

    /// Initialize to default Oracle/SunStudio compiler
    static gsJITCompilerConfig sunstudio(const int lang = gsJITLang::CXX)
    {
        switch(lang)
        {
        case (gsJITLang::C) :
            return gsJITCompilerConfig("sunstudio",
                                       "-O3 -shared",
                                       "c",
                                       "-o ");
            break;
        case (gsJITLang::CXX) :
            return gsJITCompilerConfig("sunstudio",
                                       "-O3 -shared",
                                       "cxx",
                                       "-o ");
            break;
        case (gsJITLang::Fortran) :
            return gsJITCompilerConfig("sunstudio",
                                       "-O3 -shared",
                                       "F90",
                                       "-o ");
            break;
        default :
            GISMO_ERROR("Error : Invalid compiler language.");
        }
    }
    
    /// Try to initialize compiler automatically based on the context
    static gsJITCompilerConfig guess()
    {
#       if defined(__INTEL_COMPILER)
#       if defined(__ICC)
        return intel(gsJITLang::CXX);
#       else
        return intel(gsJITLang::Fortran);
#       endif
        
#       elif  defined(_MSC_VER)
        return msvc();
        
#       elif defined(__clang__)
        return clang();
        
#       elif defined(__GNUC__)
#       if defined(__cplusplus)
        return gcc(gsJITLang::CXX);
#       elif defined(__GFortran__)
        return gcc(gsJITLang::Fortran);
#       else
        return gcc(gsJITLang::C);
#       endif

#       elif defined(__PGIC__)
        return pgi();
        
#       elif defined(__SUNPRO_C)
        return sunstudio(gsJITLang::C);
#       elif defined(__SUNPRO_CC)
        return sunstudio(gsJITLang::CXX);
#       elif defined(__SUNPRO_F90) || defined(__SUNPRO_F95)
        return sunstudio(gsJITLang::Fortran);
            
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
    std::string temp;

private:

    /// Auto-detect temp directory
    static std::string detectTemp()
    {
        return gsFileManager::getTempPath();
    }
};

/// Print (as string) operator to be used by all derived classes
inline std::ostream &operator<<(std::ostream &os,
                                const gsJITCompilerConfig& c)
{ return c.print(os); }

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

        tmp = node->first_attribute("temp");
        if (tmp!=NULL)
            result.setTemp(tmp->value());
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
        tmp->append_attribute( makeAttribute("temp" , obj.getTemp() , data) );
        
        return tmp;
    }
};

} // namespace internal

/**
   @brief Class defining a dynamic library.

   This class stores a pointer to a dynamic library that can be
   compiled at runtime and provides extra functionality.
 */
struct gsDynamicLibrary
{
public:
    /// Default Constructor
    gsDynamicLibrary() {}

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
        *(void **)(&symbol) = (void*)GetProcAddress(handle.get(), name );
#elif defined(__APPLE__) || defined(__linux__) || defined(__unix)
        *(void **)(&symbol) = ::dlsym( handle.get(), name );
#endif
        if (!symbol)
            throw std::runtime_error("An error occured while getting symbol from the dynamic library");
        
        return symbol;
    }
    
    /// Check if handle is assigned
    operator bool() const { return (bool)handle; }
    
private:

    /// Handle to dynamic library object
#if defined(_WIN32)
    memory::shared_ptr< util::remove_pointer<HMODULE>::type > handle;
#else //if defined(__APPLE__) || defined(__linux__) || defined(__unix)
    memory::shared_ptr<void> handle;
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

    /// Constructor (using compiler configuration)
    explicit gsJITCompiler(const gsJITCompilerConfig & config)
    : kernel(), config(config)
    {}
    
    /// Assignment operator (copy)
    gsJITCompiler& operator=(gsJITCompiler const& other)
    {
        kernel << other.kernel.rdbuf();
        config = other.config;
        return *this;
    }

#   if __cplusplus >= 201103L || _MSC_VER >= 1600
    /// Constructor (move)
    gsJITCompiler(gsJITCompiler && other)
    : //kernel(std::move(other.kernel)),
      config(std::move(other.config))
    {
        kernel << other.kernel.rdbuf();
    }

    /// Assignment operator (move)
    gsJITCompiler& operator=(gsJITCompiler && other)
    {
        kernel << other.kernel.rdbuf();
        config = std::move(other.config);
        return *this;
    }
#   endif
    
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
#       if __cplusplus >= 201103L || _MSC_VER >= 1600
        size_t h = std::hash<std::string>()(getKernel().str() +
                                            config.getCmd()+config.getFlags() +
                                            config.getLang());
        return build(std::to_string(h), force);
#       else
        return build("JIT", true);
#       endif
    }
    
    /// Compile kernel source code into dynamic library
    /// (use given filename)
    gsDynamicLibrary build(const std::string &name, bool force = false)
    {
        // Prepare library name
        std::stringstream libName;

#       if   defined(_WIN32)
        libName << config.getTemp() << "." << name << ".dll";
        //(void)std::system("del /f " + libName);
        //force = true;
#       elif defined(__APPLE__)
        libName << config.getTemp() << ".lib" << name << ".dylib";
#       elif defined(unix) || defined(__unix__) || defined(__unix)
        libName << config.getTemp() << ".lib" << name << ".so";
#       else
#       error("Unsupported operating system")
#       endif
        
        // Compile library (if required)
        std::ifstream libfile(libName.str().c_str());
        if(!libfile || force)
        {
            // Write kernel source code to file
            std::stringstream srcName;
#           ifdef _WIN32
            srcName<< config.getTemp() << "\\." << name << "." << config.getLang();
            std::ofstream file(srcName.str().c_str());
            file << "#ifdef __cplusplus\n";
            file << "#define EXPORT extern \"C\" __declspec(dllexport)\n";
            file << "#endif\n";
#           else
            srcName<< config.getTemp() << "." << name << "." << config.getLang();
            std::ofstream file(srcName.str().c_str());
            file << "#ifdef __cplusplus\n";
            file << "#define EXPORT extern \"C\"\n";
            file << "#endif\n";
#           endif
            file << getKernel().str() <<"\n";
            file.close();

            // Compile kernel source code into library
            std::stringstream systemcall;

#           ifdef _WIN32
            // double quotes are better than single quotes..
            systemcall << "\"\"" << config.getCmd() << "\" "
                       << config.getFlags() << " \""
                       << srcName.str()     << "\" "
                       << config.getOut() << "\"" << libName.str() << "\"\"";
#           else
            systemcall << "\"" << config.getCmd() << "\" "
                       << config.getFlags() << " \""
                       << srcName.str()     << "\" "
                       << config.getOut() << "\"" << libName.str() << "\"";
#           endif
            
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
inline std::ostream &operator<<(std::ostream &os, const gsJITCompiler& c)
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
