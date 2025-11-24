#!/usr/bin/env python3
"""
fracture_generator_vtk_simple.py

Скрипт для генерации трещин в 2D области.
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


def num_sqrt(N):
    """
    по чилу N возвращает ближайшее число больше N, и которое является квадратом
    """
    root = int(math.sqrt(N))
    if root * root == N:
        return root * root
    else:
        return (root + 1) * (root + 1)
    

def find_max_square(l, alpha_mean, d_alpha):
    """
    функция, которая находит максимальные пол-длины квадрата, в который влезает трещина со своей сеткой (размера N_chi)
    на вход l- длина трещины (в ячейках)
    aplha_mean и d_alpha- средний угол и разброс угла (все в градусах) нормали к трещине относительно вертикали 
    """

    alpha_mean_rad = alpha_mean / 180 * math.pi
    d_alpha_rad = d_alpha * 180 / math.pi

    phi = math.atan(N_chi / (N_chi + l)) #угол в радианах между диагональю и линией трещины

    l_diag = ((l + N_chi) ** 2 + N_chi ** 2) ** 0.5 #длина диагонали в прямоугольнике

    x_max = l_diag
    y_max = l_diag

    """
    пока сделано решение не оптимальное, и просто проверяем, что трещина под вообще любыми углами влезет в квадрат свой 
    в будущем надо будет дореализовать данный блок

    #проверяем может ли лежать диагональ на горизонтали
    if (abs(abs(alpha_mean_rad) - abs(d_alpha_rad))) <= phi:
        x_max = l_diag

        angle_max = abs(alpha_mean_rad) + abs(d_alpha_rad)
        y_max = l_diag * math.sin(angle_max + phi)

    #трещина наклонена против часовой стрелки
    elif (alpha_mean_rad - d_alpha_rad - phi) > 0:
        x_max = l_diag * math.cos(alpha_mean_rad - d_alpha_rad - phi)
        y_max = l_diag * math.sin(alpha_mean_rad + d_alpha_rad + phi)
    #трещина наклонена по часовой стрелке
    else: 
        x_max = l_diag * math.cos(abs(alpha_mean_rad) - d_alpha_rad - phi)
        y_max = l_diag * math.sin(abs(alpha_mean_rad) + d_alpha_rad + phi)
    """
    
    return x_max, y_max


# Размещение трещин

def generate_fractures(N, mean_angle, angle_spread, mean_length, length_spread, bbox, buffer):
    """
    Генерирует N трещин внутри области bbox = (xmin, xmax, ymin, ymax).
    Все трещины — отрезки, которые не пересекаются и находятся на расстоянии не меньше buffer.
    """
    xmin, xmax, ymin, ymax = bbox
    fractures = []

    placed = True
    
    #succes_fr = N

    
    N_len_sq = num_sqrt(N) # на сколько квадратиков делим квадрат
    nums_in_sq = range(N_len_sq) #номер квадратиков в квадрате (начинаются с 0 и заканчиваются N-1)
    sq_num_places = random.sample(nums_in_sq, N) #выбрали номера для квадратиков, в которых будут трещины

    cols_and_rows = int(N_len_sq ** 0.5) #сколько рядов и колонок будет в области
    len_x = (xmax - xmin) / cols_and_rows
    len_y = (ymax - ymin) / cols_and_rows

    max_len_x, max_len_y = find_max_square(mean_length + length_spread, mean_angle, angle_spread)
        

    if len_x / 2 < max_len_x or len_y / 2 < max_len_y:
            print(f"трещина с угломом: {mean_length}+-{length_spread} может не влезть. Необходимо изменять параметры")
            placed = False
            return fractures, placed
        
    for i in sq_num_places:
        num_x = i % cols_and_rows
        num_y = i // cols_and_rows 
            
        cx = xmin + (1/2 + num_x) * len_x #координата x центра трещины
        cy = ymin + (1/2 + num_y) * len_y #координата y центра трещины                               

        
        # случайный угол (в градусах)
        angle_deg = random.uniform(mean_angle - angle_spread, mean_angle + angle_spread)
        angle = math.radians(angle_deg)

        # случайная длина
        length = random.uniform(mean_length - length_spread, mean_length + length_spread)
        if length <= 0:
            continue

        dx = math.cos(angle) * (length / 2)
        dy = math.sin(angle) * (length / 2)

        A = (cx - dx, cy - dy)
        B = (cx + dx, cy + dy)

        fractures.append((A[0], A[1], B[0], B[1]))

        
 
        """"
        #placed = False
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
        """

        

    return fractures, placed


def main():
    parser = argparse.ArgumentParser(description="Генерация трещин без пересечений в квадратах")
    parser.add_argument('--xmin', type=float)
    parser.add_argument('--xmax', type=float)
    parser.add_argument('--ymin', type=float)
    parser.add_argument('--ymax', type=float)
    parser.add_argument('--N', type=int, help='Количество трещин')
    parser.add_argument('--angle-mean', type=float)
    parser.add_argument('--angle-spread', type=float)
    parser.add_argument('--len-mean', type=float)
    parser.add_argument('--len-spread', type=float)
    parser.add_argument('--out-folder', type=str)
    args = parser.parse_args()

    
    buffer = 3 * math.sqrt(2) * N_chi

    bbox = (args.xmin, args.xmax, args.ymin, args.ymax)
    fractures, placed = generate_fractures(
        args.N, args.angle_mean, args.angle_spread,
        args.len_mean, args.len_spread,
        bbox, buffer
    )
    out_file = args.out_folder + f"/fractures_{len(fractures)}_{(args.xmin + args.xmax)/2}_{(args.ymin + args.ymax)/2}_{(args.xmax - args.xmin)/2}_{(args.ymax - args.ymin)/2}_{args.len_mean}_{args.len_spread}_{args.angle_mean}_{args.angle_spread}.txt"

    # Сохраняем координаты трещин в файл
    if placed:
        with open(out_file, 'w', encoding='utf-8') as f:
            for (x1, y1, x2, y2) in fractures:
                f.write(f"{x1:.6f} {y1:.6f} {x2:.6f} {y2:.6f}\n")

        print(f"Создан файл {out_file} ({len(fractures)} трещин)")
    


if __name__ == '__main__':
    main()
