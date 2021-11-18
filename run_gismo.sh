# Run gismo files
echo "Solve and save biharmonic example with Argyris space!"
#echo ./bin/biharmonic_argyris_example -g -1 -l 4
#./bin/biharmonic_argyris_example -g -1 -l 4


# for rr in 3 
# do
#     for pp in 4 5 
#     do
#         for i in 1123
#         do
#         echo ./bin/biharmonic_multiPatch_example -g $i -p 5 -r 4 -l 5 -P $pp -R $rr --latex
#         ./bin/biharmonic_multiPatch_example -g $i -p 5 -r 4 -l 5 -P $pp -R $rr --latex
#         done
#     done
# done

# for mu in 100 250 500 750 1000 2000 5000 10000 
# do
#     for i in 1100
#         do
#         echo ./bin/biharmonic_multiPatch_example -g $i -p 3 -r 2 -l 7 -m nitsche -y $mu --csv
#         ./bin/biharmonic_multiPatch_example -g $i -p 3 -r 2 -l 7 -m nitsche -y $mu --csv
#     done
# done
# 
# for i in 1100
#     do
#     echo ./bin/biharmonic_multiPatch_example -g $i -p 3 -r 2 -l 7 --csv
#     ./bin/biharmonic_multiPatch_example -g $i -p 3 -r 2 -l 7 --csv
# done
# 
# for mu in 0.001 0.01 0.1 0.5 1 2 5 10 20 50 100 250 500 750 1000 2000 5000 10000 
# do
#     for i in 1100
#         do
#         echo ./bin/biharmonic_multiPatch_example -g $i -p 4 -r 3 -l 7 -m nitsche -y $mu --csv
#         ./bin/biharmonic_multiPatch_example -g $i -p 4 -r 3 -l 7 -m nitsche -y $mu --csv
#     done
# done

# for i in 1311
#     do
#     echo ./bin/biharmonic_multiPatch_example -g $i -p 3 -r 2 -l 6 --csv
#     ./bin/biharmonic_multiPatch_example -g $i -p 3 -r 2 -l 6 --csv
# done

# for i in 1311
#     do
#     echo ./bin/biharmonic_multiPatch_example -g $i -p 4 -r 3 -l 6 --csv
#     ./bin/biharmonic_multiPatch_example -g $i -p 4 -r 3 -l 6 --csv
# done
# 
# for mu in 0.001 0.01 0.1 0.5 1 2 5 10 20 50 100 250 500 750 1000 2000 5000 10000 
# do
#     for i in 1311
#         do
#         echo ./bin/biharmonic_multiPatch_example -g $i -p 3 -r 2 -l 6 -m nitsche -y $mu --csv
#         ./bin/biharmonic_multiPatch_example -g $i -p 3 -r 2 -l 6 -m nitsche -y $mu --csv
#     done
# done

# for i in 1601
#     do
#     echo ./bin/biharmonic_multiPatch_example -g $i -p 3 -r 1 -l 5 --csv
#     ./bin/biharmonic_multiPatch_example -g $i -p 3 -r 1 -l 5 --csv
# done
# 
# for i in 1601
#     do
#     echo ./bin/biharmonic_multiPatch_example -g $i -p 4 -r 2 -l 5 --csv
#     ./bin/biharmonic_multiPatch_example -g $i -p 4 -r 2 -l 5 --csv
# done
# 
# for i in 1601
#     do
#     echo ./bin/biharmonic_multiPatch_example -g $i -p 5 -r 3 -l 5 --csv
#     ./bin/biharmonic_multiPatch_example -g $i -p 5 -r 3 -l 5 --csv
# done



# for mu in 0.001 0.01 0.1 0.5 1 2 5 10 20 50 100 250 500 750 1000 2000 5000 10000 
# do
#     for i in 1600
#         do
#         echo ./bin/biharmonic_multiPatch_example -g $i -p 3 -r 1 -l 4 -m nitsche -y $mu --csv
#         ./bin/biharmonic_multiPatch_example -g $i -p 3 -r 1 -l 4 -m nitsche -y $mu --csv
#     done
# done
# 
# for mu in 0.001 0.01 0.1 0.5 1 2 5 10 20 50 100 250 500 750 1000 2000 5000 10000 
# do
#     for i in 1600
#         do
#         echo ./bin/biharmonic_multiPatch_example -g $i -p 4 -r 2 -l 4 -m nitsche -y $mu --csv
#         ./bin/biharmonic_multiPatch_example -g $i -p 4 -r 2 -l 4 -m nitsche -y $mu --csv
#     done
# done
# 

for i in 1021 1121 1311 1601
    do
    echo ./bin/biharmonic_multiPatch_example -g $i -p 3 -r 2 -l 5 -m nitsche -y -2 -w 4 --csv
    ./bin/biharmonic_multiPatch_example -g $i -p 3 -r 2 -l 5 -m nitsche -y -2 -w 4 --csv

    echo ./bin/biharmonic_multiPatch_example -g $i -p 4 -r 3 -l 5 -m nitsche -y -2 -w 4 --csv
    ./bin/biharmonic_multiPatch_example -g $i -p 4 -r 3 -l 5 -m nitsche -y -2 -w 4 --csv

    echo ./bin/biharmonic_multiPatch_example -g $i -p 5 -r 4 -l 5 -m nitsche -y -2 -w 4 --csv
    ./bin/biharmonic_multiPatch_example -g $i -p 5 -r 4 -l 5 -m nitsche -y -2 -w 4 --csv
done


for mu in 0.001 0.01 0.1 0.5 1 2 5 10 20 50 100 250 500 750 1000 2000 5000 10000
do
    for i in 1021 1121 1311 1601
        do
        echo ./bin/biharmonic_multiPatch_example -g $i -p 3 -r 2 -l 5 -m nitsche -y $mu --csv
        ./bin/biharmonic_multiPatch_example -g $i -p 3 -r 2 -l 5 -m nitsche -y $mu --csv
    done
done

for mu in 0.001 0.01 0.1 0.5 1 2 5 10 20 50 100 250 500 750 1000 2000 5000 10000
do
    for i in 1021 1121 1311 1601
        do
        echo ./bin/biharmonic_multiPatch_example -g $i -p 4 -r 3 -l 5 -m nitsche -y $mu --csv
        ./bin/biharmonic_multiPatch_example -g $i -p 4 -r 3 -l 5 -m nitsche -y $mu --csv
    done
done

for mu in 0.001 0.01 0.1 0.5 1 2 5 10 20 50 100 250 500 750 1000 2000 5000 10000
do
    for i in 1021 1121 1311 1601
        do
        echo ./bin/biharmonic_multiPatch_example -g $i -p 5 -r 4 -l 5 -m nitsche -y $mu --csv
        ./bin/biharmonic_multiPatch_example -g $i -p 5 -r 4 -l 5 -m nitsche -y $mu --csv
    done
done

# for mu in 0.0000000001 0.000000001 0.00000001 0.0000001 0.000001 0.00001 0.0001
# do
#     for i in 1100
#         do
#         echo ./bin/biharmonic_multiPatch_example -g $i -p 3 -r 2 -l 5 -m nitsche -y $mu --csv
#         ./bin/biharmonic_multiPatch_example -g $i -p 3 -r 2 -l 5 -m nitsche -y $mu --csv
#     done
# done

# 

# for i in 1100
#     do
#     echo ./bin/biharmonic_multiPatch_example -g $i -p 3 -r 2 -l 5 --csv
#     ./bin/biharmonic_multiPatch_example -g $i -p 3 -r 2 -l 5 --csv
# done


# for i in 1021 1121 1311 1601
#     do
#     echo ./bin/biharmonic_multiPatch_example -g $i -p 4 -r 3 -l 5 --csv
#     ./bin/biharmonic_multiPatch_example -g $i -p 4 -r 3 -l 5 --csv
# done
#     
# for i in 1311 1601
#     do
#     echo ./bin/biharmonic_multiPatch_example -g $i -p 5 -r 4 -l 5 --csv
#     ./bin/biharmonic_multiPatch_example -g $i -p 5 -r 4 -l 5 --csv
# done

# for mu in -1
# do
#     for i in 1601
#         do
#         echo ./bin/biharmonic_multiPatch_example -g $i -p 3 -r 2 -l 5 -m nitsche -y $mu --csv
#         ./bin/biharmonic_multiPatch_example -g $i -p 3 -r 2 -l 5 -m nitsche -y $mu --csv
#     done
# done
# 
# 
# for mu in -1
# do
#     for i in 1021 1121 1311
#         do
#         echo ./bin/biharmonic_multiPatch_example -g $i -p 4 -r 3 -l 5 -m nitsche -y $mu --csv
#         ./bin/biharmonic_multiPatch_example -g $i -p 4 -r 3 -l 5 -m nitsche -y $mu --csv
#     done
# done
# 
# 
# for mu in -1
# do
#     for i in 1021 1121 1311
#         do
#         echo ./bin/biharmonic_multiPatch_example -g $i -p 5 -r 4 -l 5 -m nitsche -y $mu --csv
#         ./bin/biharmonic_multiPatch_example -g $i -p 5 -r 4 -l 5 -m nitsche -y $mu --csv
#     done
# done

# for ((mu=1; mu<5; mu+=1))
# do
#     for i in 1100
#     do
#     echo ./bin/biharmonic_multiPatch_example -g $i -p 3 -r 2 -l 6 -m nitsche -y $mu --csv
#     ./bin/biharmonic_multiPatch_example -g $i -p 3 -r 2 -l 6 -m nitsche -y $mu --csv
#     done
# done

#for i in 1110 1120 1122
#do
#echo ./bin/biharmonic_argyris_example -g $i -p $pp -r 1 -l 0 --solution --mesh
#./bin/biharmonic_argyris_example -g $i -p $pp -r 1 -l 0 --solution --mesh
#done

echo Finished!
