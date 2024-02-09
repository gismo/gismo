
# Start Matlab

Compile gismo (TARGET_ARCHITECTURE=none), then on linux start matlab as
```
LD_LIBRARY_PATH=/gismo/build/lib matlab
```

# Generating

Set correctly srcpath and binpath inside makegismo.m

Run:
```
cd /gismo//plugins/gsMex/lib
makegismo
build(definegismo)
```


# Using

```
addpath('/gismo//plugins/gsMex/lib')
```
