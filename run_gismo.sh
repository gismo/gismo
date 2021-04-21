# Run gismo files
echo "Solve and save biharmonic example with Argyris space!"
#echo ./bin/biharmonic_argyris_example -g -1 -l 4
#./bin/biharmonic_argyris_example -g -1 -l 4

for rr in 1
do
    for pp in 4
    do
        for i in 1000 1012 1020 1100 1110 1120
        do
        echo ./bin/biharmonic_argyris_example -g $i -p $pp -r $rr -l 6 --simplified --csv
        ./bin/biharmonic_argyris_example -g $i -p $pp -r $rr -l 6 --simplified --csv

        echo ./bin/biharmonic_argyris_example -g $i -p $pp -r $rr -l 6 --simplified --interpolation --csv
        ./bin/biharmonic_argyris_example -g $i -p $pp -r $rr -l 6 --simplified --interpolation --csv
        done
    done
done

#for i in 1110 1120 1122
#do
#echo ./bin/biharmonic_argyris_example -g $i -p $pp -r 1 -l 0 --solution --mesh
#./bin/biharmonic_argyris_example -g $i -p $pp -r 1 -l 0 --solution --mesh
#done

echo Finished!
