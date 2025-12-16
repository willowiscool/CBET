#!/usr/bin/bash

echo strong scaling test!
rm ./*.o ./cbet
make
threads=$(($(lscpu -e | wc -l)-1))
for i in $(seq 1 $threads)
do
	echo "$i thread(s)"
	time OMP_NUM_THREADS="$i" ./cbet > /dev/null
	time OMP_NUM_THREADS="$i" ./cbet > /dev/null
	time OMP_NUM_THREADS="$i" ./cbet > /dev/null
done

