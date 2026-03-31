
#pragma once

#include <gsCore/gsConfig.h>

class gsModule
{
public:
    virtual ~gsModule()   { }
    virtual bool load()   { GISMO_NO_IMPLEMENTATION }
    virtual bool unload() { GISMO_NO_IMPLEMENTATION }
};

/*
class gsModules
{
    template <typename T>
    typename gsSparseSolver<T>::uPtr SparseSolver(std::string name)
    {
        return std::make_unique<gsSparseSolver<T>>();
    }

    template <typename T>
    typename gsOptimizer<T>::uPtr Optimizer(std::string name)
    {
        return static_cast<gsOptimizer<T>*>(m_optimizers[name]);
    }

    template <typename T>
    void add(std::string name, gsOptimizer<T>* optimizer)
    {
        m_optimizers[name] = optimizer;
    }

protected:
    std::map<std::string, void*> m_optimizers; // registered optimizers
}
*/



/// in gsMumps
// How does gsMumps add itself to the list of optimizers in gsModules? Does it call gsModules::add() in its constructor or is there a factory pattern involved?


/* 

NOTES:

This is a classic "Circular Dependency" or "Inversion of Control" problem. If the Main Library knows about the Module, and the Module depends on the Main Library, you have a loop that CMake and the compiler will hate.

In a modular C++ framework, you handle this by **never letting the Main Library know the Module exists.** Instead, the Main Library defines a **contract (Interface)**, and the Module fulfills it.

---

## 1. The Design Pattern: The Registry
To allow the Main Library to use `OptionalClass` without knowing it exists at compile time, you use a **Registry Pattern**.

### Step A: The Abstract Interface (In Main Library)
The Main Library provides a header that defines what an "Optional Class" looks like, but doesn't implement it.

```cpp
// MainLib/OptionalInterface.h
class OptionalInterface {
public:
    virtual ~OptionalInterface() = default;
    virtual void doSomething() = 0;
};
```

### Step B: The Registry (In Main Library)
The Main Library holds a static pointer or a map to these interfaces.

```cpp
// MainLib/Registry.h
class Registry {
public:
    static void registerModule(std::unique_ptr<OptionalInterface> module) {
        instance().m_module = std::move(module);
    }
    static OptionalInterface* getModule() {
        return instance().m_module.get();
    }
private:
    static Registry& instance() { static Registry r; return r; }
    std::unique_ptr<OptionalInterface> m_module;
};
```

### Step C: The Module Implementation (In the Module)
The module includes the interface, implements it, and "hooks" itself into the library.

```cpp
// MyModule/OptionalClass.cpp
class OptionalClass : public OptionalInterface {
    void doSomething() override { 
    // Real logic here
    }
};

// When the module is loaded:
Registry::registerModule(std::make_unique<OptionalClass>());
```

---

## 2. How the Main Library knows it "Doesn't Have It"
In your code, the "existence" of the module is simply a null-check.

```cpp
void some_function_in_main_lib() {
    auto* plugin = Registry::getModule();
    if (plugin) {
        plugin->doSomething();
    } else {
        // Main library knows OptionalClass is missing!
        std::cout << "Optional feature not installed." << std::endl;
    }
}
```

---

## 3. How to handle this in CMake
To make this work seamlessly with your versioning and dependency requirements, you use **Compile Definitions**.

If you are building the Main Library and you want to know *at compile time* if the module is being included (for optimization or to enable certain code paths), you use `target_compile_definitions`.

```cmake
# In your "Master" CMakeLists.txt
if (ENABLE_OPTIONAL_MODULE)
    add_subdirectory(OptionalModule)
    target_compile_definitions(MainLibrary PUBLIC HAS_OPTIONAL_MODULE)
endif()
```

Then in your C++ code:
```cpp
#ifdef HAS_OPTIONAL_MODULE
    #include "OptionalClass.h" // Only safe if the module is built-in
#endif
```



---

## 4. Forward Declarations?
You generally **cannot** use a simple forward declaration (`class OptionalClass;`) if the Main Library needs to instantiate the object or call its methods. 
* **Forward declaration** only works if you are storing a *pointer* to the class and never touching its internals.
* **Interfaces (Abstract base classes)** are the correct way to let the Main Library interact with a class whose implementation lives in a separate, optional Git repo.

### Summary
1.  **Main Library** defines the `Interface`.
2.  **Module** implements the `Interface`.
3.  **Main Library** checks if a pointer to that interface is `nullptr`.
4.  **CMake** uses a `bool` option to decide whether to link the module or not.

**Would you like to see how to automate the "Registration" step so the Module registers itself automatically when the library is loaded?**

*/
