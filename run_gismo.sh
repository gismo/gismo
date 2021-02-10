# Run gismo files
echo "P = 4, R = 1"
echo -g 8 -k 4 -l 5 --latex_plot -P 1
./bin/g1biharmonicTwoPatch_example -g 8 -k 4 -l 5 --latex_plot -P 1

echo -g 8 -k 4 -l 5 --latex_plot -P 1 -p 2 -r 1
./bin/g1biharmonicTwoPatch_example -g 8 -k 4 -l 5 --latex_plot -P 1 -p 2 -r 1

echo -g 8 -k 4 -l 5 --latex_plot -P 1 -p 3 -r 2
./bin/g1biharmonicTwoPatch_example -g 8 -k 4 -l 5 --latex_plot -P 1 -p 3 -r 2

echo -g 8 -k 4 -l 5 --latex_plot -P 1 -p 4 -r 3
./bin/g1biharmonicTwoPatch_example -g 8 -k 4 -l 5 --latex_plot -P 1 -p 4 -r 3

echo Finished!
