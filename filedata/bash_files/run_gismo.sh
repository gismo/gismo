# Run gismo files
echo First example 
make -j g1biharmonic_example && nohup /usr/bin/time -v ./bin/g1biharmonic_example -e 1e-03 -t 12 -g 1 -l 5 -k 4

echo Second example
make -j g1biharmonic_example && nohup /usr/bin/time -v ./bin/g1biharmonic_example -e 1e-03 -t 4 -g 1 -l 1 -k 9 --localEdge

echo Finished!
