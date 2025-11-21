#!/usr/bin/env python3
"""
fracture_generator_vtk_simple.py

Простой скрипт для генерации трещин в 2D области.
Трещины — это просто отрезки, которые не пересекаются и не лежат слишком близко друг к другу.

Результат — файл .txt с координатами концов трещин:
x1 y1 x2 y2
x1 y1 x2 y2
...

Запуск:
python3 fracture_generator_vtk_simple.py \
    --xmin -200 --xmax 200 \
    --ymin -200 --ymax 200 \
    --N 10 \
    --out fractures.txt \
    --angle-mean 45 --angle-spread 10
"""

import argparse
import math
import random
import os
from params import N_chi  # параметр ширины буфера из params.py

# =======================================
# 1. Простейшая геометрия
# =======================================

def distance_point_to_segment(A, B, P):
    """Расстояние от точки P до отрезка AB"""
    # Вектор AB и AP
    ABx, ABy = B[0] - A[0], B[1] - A[1]
    APx, APy = P[0] - A[0], P[1] - A[1]
    AB_len2 = ABx**2 + ABy**2
    if AB_len2 == 0:
        # Отрезок вырожден в точку
        return math.hypot(APx, APy)
    # Находим проекцию точки P на отрезок
    t = (APx * ABx + APy * ABy) / AB_len2
    t = max(0.0, min(1.0, t))
    projx = A[0] + ABx * t
    projy = A[1] + ABy * t
    return math.hypot(P[0] - projx, P[1] - projy)


def distance_between_segments(A, B, C, D):
    """Минимальное расстояние между двумя отрезками AB и CD"""
    if segments_intersect(A, B, C, D):
        return 0.0
    return min(
        distance_point_to_segment(A, B, C),
        distance_point_to_segment(A, B, D),
        distance_point_to_segment(C, D, A),
        distance_point_to_segment(C, D, B)
    )


def ccw(P, Q, R):
    """Проверка ориентации трёх точек"""
    return (R[1] - P[1]) * (Q[0] - P[0]) > (Q[1] - P[1]) * (R[0] - P[0])


def segments_intersect(A, B, C, D):
    """Пересекаются ли два отрезка"""
    return (ccw(A, C, D) != ccw(B, C, D)) and (ccw(A, B, C) != ccw(A, B, D))


# =======================================
# 2. Размещение трещин
# =======================================

def generate_fractures(N, mean_angle, angle_spread, mean_length, length_spread, bbox, buffer):
    """
    Генерирует N трещин внутри области bbox = (xmin, xmax, ymin, ymax).
    Все трещины — отрезки, которые не пересекаются и находятся на расстоянии не меньше buffer.
    """
    xmin, xmax, ymin, ymax = bbox
    fractures = []
    
    #succes_fr = N

    for i in range(N):
        placed = False
        for attempt in range(100000):  # ограничим число попыток
            # случайный угол (в градусах)
            angle_deg = random.uniform(mean_angle - angle_spread, mean_angle + angle_spread)
            angle = math.radians(angle_deg)

            # случайная длина
            length = random.uniform(mean_length - length_spread, mean_length + length_spread)
            if length <= 0:
                continue

            dx = math.cos(angle) * (length / 2)
            dy = math.sin(angle) * (length / 2)

            # координаты центра трещины — чтобы она не выходила за границы
            cx = random.uniform(xmin + abs(dx) + buffer, xmax - abs(dx) - buffer)
            cy = random.uniform(ymin + abs(dy) + buffer, ymax - abs(dy) - buffer)

            A = (cx - dx, cy - dy)
            B = (cx + dx, cy + dy)

            # проверяем пересечения и близость
            conflict = False
            for (x1, y1, x2, y2) in fractures:
                C = (x1, y1)
                D = (x2, y2)
                if distance_between_segments(A, B, C, D) < buffer:
                    conflict = True
                    break

            if not conflict:
                fractures.append((A[0], A[1], B[0], B[1]))
                placed = True
                break

        if not placed:
            print(f"Не возможно разместить трещины!")
            #succes_fr -= 1
            #break

    return fractures


def main():
    parser = argparse.ArgumentParser(description="Простая генерация трещин без пересечений")
    parser.add_argument('--xmin', type=float, default=-100)
    parser.add_argument('--xmax', type=float, default=100)
    parser.add_argument('--ymin', type=float, default=-100)
    parser.add_argument('--ymax', type=float, default=100)
    parser.add_argument('--N', type=int, default=5, help='Количество трещин')
    parser.add_argument('--angle-mean', type=float, default=45)
    parser.add_argument('--angle-spread', type=float, default=10)
    parser.add_argument('--len-mean', type=float, default=20)
    parser.add_argument('--len-spread', type=float, default=5)
    parser.add_argument('--out-folder', type=str)
    args = parser.parse_args()

    
    buffer = 3 * math.sqrt(2) * N_chi

    bbox = (args.xmin, args.xmax, args.ymin, args.ymax)
    fractures = generate_fractures(
        args.N, args.angle_mean, args.angle_spread,
        args.len_mean, args.len_spread,
        bbox, buffer
    )
    out_file = args.out_folder + f"/fractures_{len(fractures)}_{(args.xmin + args.xmax)/2}_{(args.ymin + args.ymax)/2}_{(args.xmax - args.xmin)/2}_{(args.ymax - args.ymin)/2}_{args.len_mean}_{args.len_spread}_{args.angle_mean}_{args.angle_spread}.txt"

    # Сохраняем координаты трещин в файл
    with open(out_file, 'w', encoding='utf-8') as f:
        for (x1, y1, x2, y2) in fractures:
            f.write(f"{x1:.6f} {y1:.6f} {x2:.6f} {y2:.6f}\n")

    print(f"Создан файл {out_file} ({len(fractures)} трещин)")


if __name__ == '__main__':
    main()
