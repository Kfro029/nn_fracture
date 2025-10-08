#!/bin/bash

#Author:  Alena V. Favorskaya, 2024 12 04

function empty_folder {
	rm $1 -r
	mkdir -p $1
}

empty_folder "VTK"
empty_folder "TXT"
empty_folder "MESHES"

rm "conf.conf"
python3 Impulse_Generator.py
python3 gen_conf_Копия_2.py >> "conf.conf"
rm "receivers.csv"
/home/faki/one_fracture_chimera/rect/build/rect "conf.conf"





