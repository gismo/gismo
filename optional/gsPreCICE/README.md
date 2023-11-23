# PreCICE Adapter for G+Smo

... General introduction ...

## Installation instructions
**Prerequisite:** [Check where suits you best](https://precice.org/installation-overview.html) to install the preCICE dependencies.
<p><b>Dependencies:</b></p>
<html>
<table>
  <tr>
    <th>preCICE</th>
    <th>2.5.0</th>
  </tr>
  <tr>
    <th> C++ compiler</th>
    <th>C++14</th>
  </tr>
  <tr>
    <th>CMake</th>
    <th>>=3.16.1</th>
  </tr>
  <tr>
    <th>Eigen</th>
    <th>>=3.3.7</th>
  </tr>
  <tr>
    <th>PETSc</th>
    <th>>=3.12</th>
</table> 
</html>

1. Check [this page](https://precice.org/quickstart.html#installation) to install PreCICE and openFOAM. The dependencies for PreCICE can be found in the above table. For MacOS, you can install all dependencies using [Homebrew](https://brew.sh), as mentioned in PreCICE webpage.
```
brew install cmake eigen libxml2 boost petsc openmpi python3 numpy
```
2. Install few common dependencies
```
brew install pkg-config cmake git
```
3. Download and install OpenFOAM-preCICE adapter
```
git clone --branch=master --depth 1 https://github.com/precice/openfoam-adapter
cd openfoam-adapter
./Allwmake
cd ..
```
4. Enable the `gsPreCICE` submodule in G+Smo, e.g. by modifying the following line in `gismo/submodules.txt`:
```
set(SUBMODULES_TXT "<other submodules>;gsPreCICE")
```
5. Build gismo PreCICE
```
cd gismo
mkdir gismo_PreCICE
cd gismo_PreCICE
cmake -DCMAKE_BUILD_TYPE=Debug .. 
make
```
*NOTES:* macOS users might encounter `dyld[#]:library not loaded:@rpath/library_name.dylib` error. If the library exists but is not in a standard location, you may need to update your library paths. 
'''
find / -name library_name.11.dylib 2>/dev/null
'''
This gives you the directory of the library. 

Set the `DYLD_LIBRARY_PATH` environment variable to include the directory where `library_name.dylib` is located:
```
export DYLD_LIBRARY_PATH=/path/to/your/lib:$DYLD_LIBRARY_PATH
```
Replace `/path/to/your/lib` with the actual path where `library_name.dylib` is located.
6. Build a PreCICE example in G+Smo
```
make flow-over-heated-plate
```

## Running a tutorial

### 0. Vertical beam

(Requires `gsKLShell`)

0. Build the example file (`gismo/optional/gsPreCICE/examples/vertical-beam.cpp`)
```
cd <gismo/build/directory>
make vertical-beam
```
1. Go to the directory `gismo/optional/gsPreCICE/examples/vertical-beam`
```
cd gismo/optional/gsPreCICE/examples/vertical-beam
```
2. In one terminal window, run the generator
```
cd fluid-python
python generator.py
```
3. In another terminal window, run the G+Smo beam
```
cd solid-gismo
./run.sh
```

### 1. Flow over heated plate

0. Build the example file (`gismo/optional/gsPreCICE/examples/flow-over-heated-plate.cpp`)
```
cd <gismo/build/directory>
make flow-over-heated-plate
```
1. Go to the directory `gismo/optional/gsPreCICE/examples/flow-over-heated-plate`
```
cd gismo/optional/gsPreCICE/examples/flow-over-heated-plate
```
2. In one terminal window, run the fluid solver (requires openFOAM)
```
cd fluid-openfoam
./run.sh
```
3. In another terminal window, run the G+Smo heat model
```
cd solid-gismo
./run.sh
```

### 2. Partitioned heat conduction

### 3. Perpendicular flap

### 4. Turek-Hron Benchmark
WIP

### 5. Elastic tube 3D ([link](https://precice.org/tutorials-elastic-tube-3d.html))
WIP
