# PreCICE Adapter for G+Smo

... General introduction ...

## Installation instructions

1. Check [this page](https://precice.org/quickstart.html#installation) to install PreCICE and openFOAM   
```
some code
```
2. Enable the `gsPreCICE` submodule in G+Smo, e.g. by modifying the following line in `gismo/submodules.txt`:
```
set(SUBMODULES_TXT "<other submodules>;gsPreCICE")
```
3. Build a PreCICE example in G+Smo
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
