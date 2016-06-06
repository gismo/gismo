======================================================================
=====                G+Smo deployment instructions               ===== 
======================================================================


Assuming that you have the source tree at /path/to/gismo, here is how
to use G+Smo in your application:


* Choose a folder for building (eg. folder /path/to/build for demonstation):

$ mkdir /path/to/build
$ cd /path/to/build

Note: If you only need to experiment with a small main.cpp file, the
far easiest way to go is to copy the file "main.cpp" to the source
folder gismo/examples.  Then issue:

$ cmake /path/to/gismo
$ make -j2 main
$ ./bin/main
Hello G+Smo.

After executing cmake, the file is included directly to the "examples"
compilation list, and the binary file is created in the subfolder
./bin.  However, if you need to use the library as an external
dependency, then continue reading these instructions.


* Configure and install on a predefined location
  (eg. ./path/to/install for demonstration):

cmake /path/to/gismo -DCMAKE_INSTALL_PREFIX=/path/to/install -DGISMO_EXAMPLES=OFF
make -j2 install

The installation will create the folder /path/to/install, containing:

./include/gismo   Include path needed by the library (eg. gismo.h)
./lib             Shared library (eg. libgismo.so)
./bin             Compiled examples (only if GISMO_EXAMPLES=ON)

There are several solutions for using the library in your project:


* Solution 1: use the provided main.cpp and CMakeLists.txt (recommended)

Copy the provided "CMakeLists.txt" and main.cpp to a folder of your
choice. In that folder issue:

$ cd /path/to/myproject
$ cmake . -Dgismo_DIR=/path/to/install

This will create a "Makefile" (therefore will overwrite any existing makefiles).
For this solution the installation of G+Smo is not mandatory. We can
also point directly to the build folder (if install was not executed):

$ cmake . -Dgismo_DIR=/path/to/build

If the configuration is successfull, then typing "make" should create
the binary file "main", which can be executed as:

$ make
$ ./main
Hello G+Smo.

Note A: CMake can produce project files for Visual Studio, Ninja,
CodeBlocks, and many other build systems,
see https://cmake.org/cmake/help/v3.0/manual/cmake-generators.7.html.
Note B: On Windows, use the graphical tool cmake-gui to configure and
compile.


* Solution 2: use the provided main.cpp and Makefile (in unix-based systems)

Copy the provided "Makefile" and main.cpp to a folder of your
choice. You need to edit the Makefile and update the variable GISMODIR
to point to your install directory where G+Smo is found:

GISMODIR :=/path/to/install

Then typing "make" should create the binary file "main", which can be
executed as:

$ make
$ ./main
Hello G+Smo.


* Solution 3: Run your favorite compiler, eg. GCC, and link with libgismo.so:

$ g++ main.cpp -I/path/to/install/include/gismo -L/path/to/install/lib -Wl,-rpath=/path/to/install/lib -lgismo -o main
$ ./main
Hello G+Smo.

Or, on Windows, using the Microsoft compiler on a developer command prompt:

> cl /EHsc main.cpp /I C:\path\to\install\include\gismo /link C:\path\to\install\lib\gismo.lib
> main.exe
Hello G+Smo.

