set -euo pipefail

BASE_DIR="/home/faki/one_fracture_chimera"
FRACTURES_DIR="$BASE_DIR/fractures_1st_sel" # _1st_sel

# Общие параметры генерации
N=3 #количество трещин в файле
VTK_TYPE="structured"   # или "polydata"
SPACING="1.0"
SAMPLES=1               # если >1 fracture_generator_vtk создаст суффиксы _0000 и т.д.

mkdir -p "$FRACTURES_DIR"

xs=(-150 -100 -50 0 50 100 150)
#xs=(50)

ys=(-100 -150 -200)
#ys=(-50)

dxs=(40 50 60)
#dxs=(40)

dys=(40 50 60)
#dys=(40)

MEAN_LENGTHs=(10)

MEAN_ANGLEs=($(seq -18 4 18))
#MEAN_ANGLEs=(10)

for x in "${xs[@]}"; do
    for y in "${ys[@]}"; do
        for dx in "${dxs[@]}"; do
            for dy in "${dys[@]}"; do
                for MEAN_LENGTH in "${MEAN_LENGTHs[@]}"; do
                    for MEAN_ANGLE in "${MEAN_ANGLEs[@]}"; do
                           
                        python3 fracture_gen_in_sq.py \
                        --xmin $((x-dx)) --xmax $((x+dx)) \
                        --ymin $((y-dy)) --ymax $((y+dy)) \
                        --N $N \
                        --out "$FRACTURES_DIR" \
                        --angle-mean $MEAN_ANGLE --angle-spread 0 \
                        --len-mean $MEAN_LENGTH --len-spread 0
                    
                    done
                done
            done
        done    
    done
done