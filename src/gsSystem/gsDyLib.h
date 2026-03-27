

namespace gismo
{


class gsDyLib
{
    void * handle = nullptr;

public:
    gsDyLib() : handle(nullptr) { }

    gsDyLib(const std::string& name)
    {
        std::string library_name = name;
        #ifdef _WIN32
		//no op
#elif __APPLE__
		library_name = std::string("./lib" + library_name + ".dylib");
#else
		//On unix the library name probably starts with `lib`
		library_name = std::string("./lib" + library_name + ".so");
#endif
        
        handle = 
#ifdef _WIN32
        reinterpret_cast<void*>(LoadLibraryA(library_name.c_str()));
#else
        dlopen(library_name.c_str(), RTLD_NOW);
#endif
        if(!handle) GISMO_ERROR("Module library did not load");
    }

    ~gsDyLib()
    {
        if(handle != nullptr)
        {
#ifdef _WIN32
            FreeLibrary(HMODULE(lib_handle));
#else
            dlclose(lib_handle);
#endif
            handle = nullptr;
        }
    }

    //deactivate copy
    gsDyLib& operator=(const gsDyLib&) = delete;
    gsDyLib(const gsDyLib&)			  = delete;

    //move assignment
    gsDyLib& operator=(gsDyLib&& other) noexcept
    {
        handle		 = other.handle;
        other.handle = nullptr;
        return *this;
    }

    gsDyLib(gsDyLib&& other) noexcept
    {
        *this = std::move(other);
    }

    plugin* getPluginInstance()
    {
        using boostrap_function = plugin* (*)();

        std::string& function_name =  ___YPLUGS_BOOTSTRAP_PROC_NAME_STR;
        auto fptr = reinterpret_cast<boostrap_function>(
#ifdef _WIN32
			GetProcAddress(HMODULE(handle), LPCSTR(function_name.c_str()));
#else
            dlsym(handle, function_name.c_str())
#endif
            );
    
        if(!fptr)
            throw std::runtime_error(
                "Did not find yplugs bootstrap function in plugin library");
        return fptr();
    }
};



} // namespace gismo

