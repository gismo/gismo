# Run gismo files
echo "P = 3, R = 1"
#echo ./bin/biharmonic_argyris_example -g -1 -l 4
#./bin/biharmonic_argyris_example -g -1 -l 4

for i in 11 12 13 14 15 16 17 18
do
   echo ./bin/biharmonic_argyris_example -g -$i -l 3 -p 2 
   ./bin/biharmonic_argyris_example -g -$i -l 3 -p 2
done

echo Finished!
