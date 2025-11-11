set -euo pipefail

BASE_DIR="/home/faki/one_fracture_chimera"
FRACTURES_DIR="$BASE_DIR/fractures" # _1st_sel

# Общие параметры генерации
N=3 #количество трещин в файле
VTK_TYPE="structured"   # или "polydata"
SPACING="1.0"
SAMPLES=1               # если >1 fracture_generator_vtk создаст суффиксы _0000 и т.д.

mkdir -p "$FRACTURES_DIR"

#xs=(-150 -100 -50 0 50 100 150)
xs=(50)
#ys=(50 0 -50 -100)
ys=(-50)

#dxs=(30 40 50)
dxs=(50)

#dys=(30 40 50)
dys=(50)

MEAN_LENGTHs=(25)

#MEAN_ANGLEs=($(seq -18 2 18))
MEAN_ANGLEs=(6)

for x in "${xs[@]}"; do
    for y in "${ys[@]}"; do
        for dx in "${dxs[@]}"; do
            for dy in "${dys[@]}"; do
                for MEAN_LENGTH in "${MEAN_LENGTHs[@]}"; do
                    for MEAN_ANGLE in "${MEAN_ANGLEs[@]}"; do
                        BASE_NAME="fractures_${N}_${x}_${y}_${dx}_${dy}_${MEAN_LENGTH}_${MEAN_ANGLE}"

                        OUT_FILE="$FRACTURES_DIR/${BASE_NAME}.txt"
                        META_FILE="$FRACTURES_DIR/${BASE_NAME}.json"
    
                        python3 fracture_gen_simp.py \
                        --xmin $((x-dx)) --xmax $((x+dx)) \
                        --ymin $((y-dy)) --ymax $((y+dy)) \
                        --N $N \
                        --out "$OUT_FILE" \
                        --angle-mean $MEAN_ANGLE --angle-spread 0 \
                        --len-mean $MEAN_LENGTH --len-spread 0
                    
                    done
                done
            done
        done    
    done
done