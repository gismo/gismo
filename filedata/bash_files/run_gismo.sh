# Run gismo files
echo First example 4b
make -j g1biharmonic_example && nohup /usr/bin/time -v ./bin/g1biharmonic_example -e 1e-03 -t 12 -g 1 -l 5 -k 4 --latex --localEdge --localGd

echo Second example
make -j g1biharmonic_example && nohup /usr/bin/time -v ./bin/g1biharmonic_example -e 1e-03 -t 12 -g 1 -l 5 -k 4 --latex --localEdge --localGd -p 2

echo Third example 
make -j g1biharmonic_example && nohup /usr/bin/time -v ./bin/g1biharmonic_example -e 1e-03 -t 12 -g 1 -l 5 -k 4 --latex --localEdge --localGd -p 3

echo Forth example
make -j g1biharmonic_example && nohup /usr/bin/time -v ./bin/g1biharmonic_example -e 1e-03 -t 12 -g 1 -l 5 -k 4 --latex --localEdge --localGd -p 4

echo Finished!
