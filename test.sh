./gen_RMAT -s $1
./gen_valid_info -in rmat-$1
./mst_reference -in rmat-$1
./validation -in_graph rmat-$1 -in_result rmat-$1.mst -in_valid rmat-$1.vinfo
