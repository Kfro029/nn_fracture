#!/bin/bash

#PBS -l walltime=1000:1:00,nodes=1:ppn=2
#PBS -N alatortsev_1st_sel
#PBS -q batch

function empty_folder {
	rm $1 -r
	mkdir -p $1
}

FRACTURES_DIR="fractures_1st_sel"

for txt_file in "$FRACTURES_DIR"/*.txt; do
	if [ -f "$txt_file" ]; then
		filename=$(basename "$txt_file")

		empty_folder "VTK"
		empty_folder "TXT"
		empty_folder "MESHES"

		rm "conf.conf"
		python3 Impulse_Generator.py
		python3 gen_conf.py "$filename" >> "conf.conf"
		#rm "receivers.csv"
		/home/faki/one_fracture_chimera/rect/build/rect "conf.conf"
	fi
done





