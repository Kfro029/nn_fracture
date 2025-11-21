set -euo pipefail

BASE_DIR="/home/faki/one_fracture_chimera"
FRACTURES_DIR="$BASE_DIR/fractures_try"

# Общие параметры генерации
N=3 #количество трещин в файле

mkdir -p "$FRACTURES_DIR"

xs=(-150)
ys=(-50)

dxs=(30)
dys=(30)

MEAN_LENGTHs=(10)

#MEAN_ANGLEs=($(seq -18 2 18))
MEAN_ANGLEs=(0)

for x in "${xs[@]}"; do
    for y in "${ys[@]}"; do
        for dx in "${dxs[@]}"; do
            for dy in "${dys[@]}"; do
                for MEAN_LENGTH in "${MEAN_LENGTHs[@]}"; do
                    for MEAN_ANGLE in "${MEAN_ANGLEs[@]}"; do
                        #BASE_NAME="fractures_${N}_${x}_${y}_${dx}_${dy}_${MEAN_LENGTH}_${MEAN_ANGLE}"
    
                        python3 fracture_gen_simp.py \
                        --xmin $((x-dx)) --xmax $((x+dx)) \
                        --ymin $((y-dy)) --ymax $((y+dy)) \
                        --N $N \
                        --out-folder "$FRACTURES_DIR" \
                        --angle-mean $MEAN_ANGLE --angle-spread 0
                    
                    done
                done
            done
        done    
    done
done
