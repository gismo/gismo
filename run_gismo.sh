# Run gismo files
echo "P = 3, R = 1"
#echo ./bin/biharmonic_argyris_example -g -1 -l 4
#./bin/biharmonic_argyris_example -g -1 -l 4

for i in 1 2 3 4 5 6 7 8
do
   echo ./bin/biharmonic_argyris_example -g -$i -l 3 -p 3
   ./bin/biharmonic_argyris_example -g -$i -l 3 -p 3
done

echo Finished!
