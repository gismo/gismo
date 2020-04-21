# Run gismo files
echo First example mixed
make -j g1biharmonic_example && nohup /usr/bin/time -v ./bin/g1biharmonic_example -e 1e-03 -t 12 -g 1 -l 5 -k 4 --latex --localGd --localEdge --localVertex 

echo Second example
make -j g1biharmonic_example && nohup /usr/bin/time -v ./bin/g1biharmonic_example -e 1e-03 -t 12 -g 1 -l 5 -k 4 --latex --localGd --localEdge --localVertex -p 2

echo Third example 
make -j g1biharmonic_example && nohup /usr/bin/time -v ./bin/g1biharmonic_example -e 1e-03 -t 12 -g 1 -l 5 -k 4 --latex --localGd --localEdge --localVertex -p 3

echo Forth example
make -j g1biharmonic_example && nohup /usr/bin/time -v ./bin/g1biharmonic_example -e 1e-03 -t 12 -g 1 -l 5 -k 4 --latex --localGd --localEdge --localVertex -p 4

echo Finished!
