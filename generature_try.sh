set -euo pipefail

BASE_DIR="/home/faki/one_fracture_chimera"
FRACTURES_DIR="$BASE_DIR/fractures_try"

# Общие параметры генерации
N=3 #количество трещин в файле
VTK_TYPE="structured"   # или "polydata"
SPACING="1.0"
SAMPLES=1               # если >1 fracture_generator_vtk создаст суффиксы _0000 и т.д.

mkdir -p "$FRACTURES_DIR"

xs=(-150 -100)
ys=(50)

dxs=(30)
dys=(30)

MEAN_LENGTHs=(25)

#MEAN_ANGLEs=($(seq -18 2 18))
MEAN_ANGLEs=(0)

for x in "${xs[@]}"; do
    for y in "${ys[@]}"; do
        for dx in "${dxs[@]}"; do
            for dy in "${dys[@]}"; do
                for MEAN_LENGTH in "${MEAN_LENGTHs[@]}"; do
                    for MEAN_ANGLE in "${MEAN_ANGLEs[@]}"; do
                        BASE_NAME="fractures_${N}_${x}_${y}_${dx}_${dy}_${MEAN_LENGTH}_${MEAN_ANGLE}"

                        OUT_FILE="$FRACTURES_DIR/${BASE_NAME}.txt"
                        META_FILE="$FRACTURES_DIR/${BASE_NAME}.json"
    
                        python3 fracture_generator_vtk.py \
                        --xmin $((x-dx)) --xmax $((x+dx)) \
                        --ymin $((y-dy)) --ymax $((y+dy)) \
                        --N $N \
                        --out "$OUT_FILE" \
                        --meta "$META_FILE" \
                        --vtk-type "$VTK_TYPE" \
                        --spacing "$SPACING" \
                        --angle-mean $MEAN_ANGLE --angle-spread 0
                    
                    done
                done
            done
        done    
    done
done