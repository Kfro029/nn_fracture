#!/bin/bash

#Author:  Alena V. Favorskaya, 2024 12 04

function empty_folder {
	rm $1 -r
	mkdir -p $1
}

FRACTURES_DIR="fractures"

for txt_file in "$FRACTURES_DIR"/*.txt; do
	if [ -f "$txt_file" ]; then
		filename=$(basename "$txt_file")

		empty_folder "VTK"
		empty_folder "TXT"
		empty_folder "MESHES"

		rm "conf.conf"
		python3 Impulse_Generator.py
		python3 gen_conf_Копия_2.py "$filename" >> "conf.conf"
		#rm "receivers.csv"
		/home/faki/one_fracture_chimera/rect/build/rect "conf.conf"
	fi
done





