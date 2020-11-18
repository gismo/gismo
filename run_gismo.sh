# Run gismo files
echo "P = 5, R = 1"
echo -g 1 -k 4 -l 6 --latex_plot -P 2 -p 4 -r 3
./bin/g1biharmonicTwoPatch_example -g 1 -k 4 -l 6 --latex_plot -P 2 -p 4 -r 3

echo -g 1 -k 4 -l 6 --latex_plot -P 2
./bin/g1biharmonicTwoPatch_example -g 1 -k 4 -l 6 --latex_plot -P 2

echo -g 1 -k 4 -l 6 --latex_plot -P 2 -p 2 -r 1
./bin/g1biharmonicTwoPatch_example -g 1 -k 4 -l 6 --latex_plot -P 2 -p 2 -r 1

echo -g 1 -k 4 -l 6 --latex_plot -P 2 -p 3 -r 2
./bin/g1biharmonicTwoPatch_example -g 1 -k 4 -l 6 --latex_plot -P 2 -p 3 -r 2

echo Finished!
