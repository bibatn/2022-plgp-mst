./gen_random -s $1 -k 2
./gen_valid_info -in random-$1
./mst_reference -in random-$1 -nIters 1
./validation -in_graph random-$1 -in_result random-$1.mst -in_valid random-$1.vinfo
