set -euo pipefail

BASE_DIR="/home/faki/one_fracture_chimera"
FRACTURES_DIR="$BASE_DIR/fractures_2nd_sel"

# Общие параметры генерации
VTK_TYPE="structured"   # или "polydata"
SPACING="1.0"
SAMPLES=1               # если >1 fracture_generator_vtk создаст суффиксы _0000 и т.д.

mkdir -p "$FRACTURES_DIR"



Ns=(3 4 5) #количество трещин

MEAN_ANGLEs=($(seq -18 2 18)) #средний угол трещин
ANGLE_SPREADs=($(seq -1 0.5 1)) #разброс по углу

MEAN_LENGTHs=(25 30 35 40) #средняя длина трещины
LENGTH_SPREADs=(2 4 6) #разброс по длине

xs=(0)
ys=(-150)

dxs=(250)
dys=(150)


for x in "${xs[@]}"; do
    for y in "${ys[@]}"; do
        for dx in "${dxs[@]}"; do
            for dy in "${dys[@]}"; do
                for MEAN_LENGTH in "${MEAN_LENGTHs[@]}"; do
                    for MEAN_ANGLE in "${MEAN_ANGLEs[@]}"; do
                        for MEAN_ANGLE in "${MEAN_ANGLEs[@]}"; do
                            for N in "${Ns[@]}"; do
                                for LENGTH_SPREAD in "${LENGTH_SPREADs[@]}"; do
                           
                                    python3 fracture_gen_in_sq.py \
                                    --xmin -250 --xmax 250 \
                                    --ymin -300 --ymax 0 \
                                    --N $N \
                                    --out "$FRACTURES_DIR" \
                                    --angle-mean $MEAN_ANGLE --angle-spread $MEAN_ANGLE \
                                    --len-mean $MEAN_LENGTH --len-spread $LENGTH_SPREAD
                                done
                            done
                        done
                    done
                done
            done
        done    
    done
done