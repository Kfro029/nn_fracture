#!/usr/bin/env python3
"""
fracture_generator_vtk.py (обновлённая версия)

— По умолчанию генерирует VTK **STRUCTURED_GRID** с маской трещин (POINT_DATA scalar "fracture").
— Может писать POLYDATA (LINES) через флаг --vtk-type polydata.
— Поддерживает bbox из --xmin/--xmax/--ymin/--ymax или из --conf (origin/spacing + bottomLeftCorner/topRightCorner).

Использование:
  python3 fracture_generator_vtk.py \
      --xmin -200 --xmax 0 --ymin -200 --ymax 0 \
      --N 5 --out fractures_xxx.vtk


"""
import argparse
import math
import random
import json
import os
import sys
from typing import List, Tuple

try:
    import numpy as np
except Exception:
    print("Этот скрипт требует numpy. !pip install numpy")
    sys.exit(1)

# ---------------------------
# Простой парсер conf.conf
# ---------------------------
import re

def parse_conf(conf_path):
    conf = {}
    with open(conf_path, 'r', encoding='utf-8') as f:
        text = f.read()

    def find_tuple_f(key):
        m = re.search(rf"{key}\s*=\s*([\-0-9.,\s]+)", text)
        if not m: return None
        parts = [p.strip() for p in m.group(1).split(',') if p.strip()!='']
        return tuple(float(x) for x in parts)

    def find_tuple_i(key):
        m = re.search(rf"{key}\s*=\s*([0-9, \t]+)", text)
        if not m: return None
        parts = [p.strip() for p in m.group(1).split(',') if p.strip()!='']
        return tuple(int(x) for x in parts)

    conf['origin'] = find_tuple_f('origin')
    conf['spacing'] = find_tuple_f('spacing')
    conf['size'] = find_tuple_i('size')

    m_bl = re.search(r"bottomLeftCorner\s*=\s*([0-9]+)\s*,\s*([0-9]+)", text)
    m_tr = re.search(r"topRightCorner\s*=\s*([0-9]+)\s*,\s*([0-9]+)", text)
    if m_bl and m_tr:
        conf['bottomLeftCorner'] = (int(m_bl.group(1)), int(m_bl.group(2)))
        conf['topRightCorner'] = (int(m_tr.group(1)), int(m_tr.group(2)))

    return conf

# ---------------------------
# Геометрия
# ---------------------------

def seg_point_distance(A, B, P):
    ABx = B[0]-A[0]; ABy = B[1]-A[1]
    APx = P[0]-A[0]; APy = P[1]-A[1]
    ab2 = ABx*ABx + ABy*ABy
    if ab2 == 0:
        return math.hypot(APx, APy)
    t = (APx*ABx + APy*ABy) / ab2
    t = max(0.0, min(1.0, t))
    projx = A[0] + ABx * t
    projy = A[1] + ABy * t
    return math.hypot(P[0]-projx, P[1]-projy)


def segments_distance(A, B, C, D):
    if segments_intersect(A, B, C, D):
        return 0.0
    return min(seg_point_distance(A,B,C), seg_point_distance(A,B,D),
               seg_point_distance(C,D,A), seg_point_distance(C,D,B))


def ccw(P, Q, R):
    return (R[1]-P[1])*(Q[0]-P[0]) > (Q[1]-P[1])*(R[0]-P[0])


def segments_intersect(A,B,C,D):
    return (ccw(A,C,D) != ccw(B,C,D)) and (ccw(A,B,C) != ccw(A,B,D))

# ---------------------------
# Параметры
# ---------------------------

def sample_param(dist, mean, spread):
    if dist == 'uniform':
        low = mean - spread
        high = mean + spread
        return random.uniform(low, high)
    elif dist == 'lognormal':
        if mean <= 0:
            raise ValueError("mean must be > 0 for lognormal")
        mu = math.log(mean)
        sigma = max(1e-6, spread)
        return float(np.random.lognormal(mu, sigma))
    else:
        raise ValueError("Unknown dist")

# ---------------------------
# Spatial hash
# ---------------------------

class SpatialHash:
    def __init__(self, bbox, cell_size):
        self.xmin, self.xmax, self.ymin, self.ymax = bbox
        self.cs = cell_size
        self.cells = {}

    def _cell_coords(self, x, y):
        ix = int(math.floor((x - self.xmin) / self.cs))
        iy = int(math.floor((y - self.ymin) / self.cs))
        return ix, iy

    def _cells_for_seg(self, A, B):
        xmin = min(A[0], B[0]); xmax = max(A[0], B[0])
        ymin = min(A[1], B[1]); ymax = max(A[1], B[1])
        ix0, iy0 = self._cell_coords(xmin, ymin)
        ix1, iy1 = self._cell_coords(xmax, ymax)
        for ix in range(ix0-1, ix1+2):
            for iy in range(iy0-1, iy1+2):
                yield (ix, iy)

    def insert(self, seg_idx, A, B):
        for c in self._cells_for_seg(A,B):
            self.cells.setdefault(c, []).append(seg_idx)

    def query(self, A, B):
        seen = set()
        for c in self._cells_for_seg(A,B):
            for idx in self.cells.get(c, ( )):
                if idx not in seen:
                    seen.add(idx)
                    yield idx

# ---------------------------
# Раскладка трещин
# ---------------------------

def place_fractures(N, mean_angle, angle_spread, mean_len, len_spread,
                    angle_dist, len_dist,
                    cluster_bbox, buffer, max_attempts_per_seg=600):
    xmin,xmax,ymin,ymax = cluster_bbox
    placed = []
    avg_len = max(1.0, mean_len)
    cell_size = max(buffer * 2.0, avg_len)
    shash = SpatialHash(cluster_bbox, cell_size)

    for _ in range(N):
        success = False
        for _ in range(max_attempts_per_seg):
            ang = sample_param(angle_dist, mean_angle, angle_spread)
            length = sample_param(len_dist, mean_len, len_spread)
            if length <= 0.01:
                continue
            a = math.radians(ang)
            dx = math.cos(a) * (length / 2.0)
            dy = math.sin(a) * (length / 2.0)
            cx_min = xmin + abs(dx) + buffer
            cx_max = xmax - abs(dx) - buffer
            cy_min = ymin + abs(dy) + buffer
            cy_max = ymax - abs(dy) - buffer
            if cx_min > cx_max or cy_min > cy_max:
                continue
            cx = random.uniform(cx_min, cx_max)
            cy = random.uniform(cy_min, cy_max)
            A = (cx - dx, cy - dy)
            B = (cx + dx, cy + dy)
            conflict = False
            for idx in shash.query(A,B):
                seg = placed[idx]
                C = (seg[0], seg[1]); D = (seg[2], seg[3])
                d = segments_distance(A,B,C,D)
                if d < buffer - 1e-9:
                    conflict = True
                    break
            if conflict:
                continue
            placed.append( (A[0],A[1], B[0],B[1]) )
            shash.insert(len(placed)-1, A, B)
            success = True
            break
        if not success:
            return None
    return placed

# ---------------------------
# VTK writers
# ---------------------------

def write_vtk_polydata_lines(path, segments, header_comment="Created by fracture_generator_vtk.py"):
    import struct
    pts = []
    cells = []
    for seg in segments:
        A = (seg[0], seg[1], 0.0)
        B = (seg[2], seg[3], 0.0)
        i0 = len(pts)
        pts.append(A); pts.append(B)
        cells.append((i0, i0+1))

    with open(path, 'wb') as f:
        f.write(b"# vtk DataFile Version 3.0\n")
        f.write(header_comment.encode() + b"\n")
        f.write(b"BINARY\n")
        f.write(b"DATASET POLYDATA\n")
        f.write(f"POINTS {len(pts)} float\n".encode())
        for p in pts:
            f.write(struct.pack('>fff', p[0], p[1], p[2]))
        f.write(b"\n")
        total_ints = sum(1 + len(c) for c in cells)
        f.write(f"LINES {len(cells)} {total_ints}\n".encode())
        for c in cells:
            f.write(struct.pack('>i', len(c)))
            for index in c:
                f.write(struct.pack('>i', index))
        f.write(b"\n")
        f.write(f"POINT_DATA {len(pts)}\n".encode())
        f.write(b"SCALARS dummy float 1\n")
        f.write(b"LOOKUP_TABLE default\n")
        for _ in pts:
            f.write(struct.pack('>f', 0.0))


import struct
import numpy as np

def write_vtk_rectilinear_grid_mask(filename, segs, cluster_bbox, spacing, thickness, header_comment="Generated RECTILINEAR_GRID"):
    """
    Записывает бинарный VTK rectilinear grid с маской трещин.
    """
    xmin, xmax, ymin, ymax = cluster_bbox
    nx = int((xmax - xmin) / spacing) + 1
    ny = int((ymax - ymin) / spacing) + 1
    nz = 1

    # координаты по осям
    x_coords = np.linspace(xmin, xmax, nx, dtype=np.float32)
    y_coords = np.linspace(ymin, ymax, ny, dtype=np.float32)
    z_coords = np.array([0.0], dtype=np.float32)

    # маска трещин
    mask = np.zeros((ny, nx), dtype=np.uint8)
    for (x1, y1, x2, y2) in segs:
        ix1, iy1 = int((x1 - xmin) / spacing), int((y1 - ymin) / spacing)
        ix2, iy2 = int((x2 - xmin) / spacing), int((y2 - ymin) / spacing)
        length = max(abs(ix2 - ix1), abs(iy2 - iy1)) + 1
        for t in np.linspace(0, 1, length):
            ix = int(ix1 + t * (ix2 - ix1))
            iy = int(iy1 + t * (iy2 - iy1))
            if 0 <= ix < nx and 0 <= iy < ny:
                mask[iy, ix] = 1

    # запись в бинарный VTK
    with open(filename, "wb") as f:
        f.write(b"# vtk DataFile Version 3.0\n")
        f.write(f"{header_comment}\n".encode("utf-8"))
        f.write(b"BINARY\n")
        f.write(b"DATASET RECTILINEAR_GRID\n")
        f.write(f"DIMENSIONS {nx} {ny} {nz}\n".encode("utf-8"))

        f.write(f"X_COORDINATES {nx} float\n".encode("utf-8"))
        f.write(x_coords.byteswap().tobytes())
        f.write(b"\n")

        f.write(f"Y_COORDINATES {ny} float\n".encode("utf-8"))
        f.write(y_coords.byteswap().tobytes())
        f.write(b"\n")

        f.write(f"Z_COORDINATES {nz} float\n".encode("utf-8"))
        f.write(z_coords.byteswap().tobytes())
        f.write(b"\n")

        f.write(f"POINT_DATA {nx*ny*nz}\n".encode("utf-8"))
        f.write(b"SCALARS fracture_mask unsigned_char 1\n")
        f.write(b"LOOKUP_TABLE default\n")
        f.write(mask.byteswap().tobytes())


# ---------------------------
# CLI и main
# ---------------------------

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--out', default='fracture_mesh.vtk', help='output vtk file')
    p.add_argument('--meta', default='fracture_meta.json', help='output json metadata')
    p.add_argument('--conf', default=None, help='path to conf.conf (optional)')
    p.add_argument('--xmin', type=float, default=-200, help='cluster bbox xmin (m)')
    p.add_argument('--xmax', type=float, default=0, help='cluster bbox xmax (m)')
    p.add_argument('--ymin', type=float, default=-200, help='cluster bbox ymin (m)')
    p.add_argument('--ymax', type=float, default=0, help='cluster bbox ymax (m)')
    p.add_argument('--N', type=int, default=5, help='number of fractures to place')
    p.add_argument('--angle-mean', type=float, default=45.0, help='mean angle (deg)')
    p.add_argument('--angle-spread', type=float, default=10.0, help='angle spread parameter')
    p.add_argument('--len-mean', type=float, default=20.0, help='mean length (m)')
    p.add_argument('--len-spread', type=float, default=5.0, help='length spread parameter')
    p.add_argument('--dist', choices=['uniform','lognormal'], default='uniform',
                   help='distribution type for both angle and length')
    p.add_argument('--buffer-nodes', type=float, default=2.0,
                   help='buffer in grid nodes (multiplied by spacing)')
    p.add_argument('--spacing', type=float, default=1.0, help='grid spacing (m)')
    p.add_argument('--samples', type=int, default=1,
                   help='how many vtk samples to produce (appends _i)')
    p.add_argument('--seed', type=int, default=None)
    
    p.add_argument('--vtk-type', choices=['polydata','structured'], default='structured',
                   help='Тип выходного VTK: polydata (линии) или structured (сетка)')
                   
    args = p.parse_args()

    if args.seed is not None:
        random.seed(args.seed)
        np.random.seed(args.seed)

    # определяем bounding box
    if args.conf:
        conf = parse_conf(args.conf)
        origin = conf.get('origin', None)
        spacing = conf.get('spacing', (args.spacing, args.spacing))
        if 'bottomLeftCorner' in conf and 'topRightCorner' in conf:
            bl = conf['bottomLeftCorner']; tr = conf['topRightCorner']
            xmin = origin[0] + bl[0] * spacing[0]
            ymin = origin[1] + bl[1] * spacing[1]
            xmax = origin[0] + tr[0] * spacing[0]
            ymax = origin[1] + tr[1] * spacing[1]
            xmin, xmax = min(xmin,xmax), max(xmin,xmax)
            ymin, ymax = min(ymin,ymax), max(ymin,ymax)
            cluster_bbox = (xmin, xmax, ymin, ymax)
            spacing_val = spacing[0]
            print("Cluster bbox from conf:", cluster_bbox, "spacing", spacing_val)
        else:
            print("conf provided but corners not found; please provide xmin/xmax/ymin/ymax")
            sys.exit(1)
    else:
        cluster_bbox = (args.xmin, args.xmax, args.ymin, args.ymax)
        spacing_val = args.spacing
        print("Cluster bbox from args:", cluster_bbox, "spacing", spacing_val)

    from params import N_chi
    from math import sqrt
    buffer = N_chi * 2 * sqrt(2)

    # генерация
    base_out = args.out
    base_meta = args.meta
    for sample_i in range(args.samples):
    
        out_vtk = base_out if args.samples == 1 else f"{os.path.splitext(base_out)[0]}_{sample_i:04d}.vtk"
        out_meta = base_meta if args.samples == 1 else f"{os.path.splitext(base_meta)[0]}_{sample_i:04d}.json"

        segs = place_fractures(
            args.N,
            args.angle_mean, args.angle_spread,
            args.len_mean, args.len_spread,
            args.dist, args.dist,
            cluster_bbox,
            buffer
        )
        if segs is None:
            print("Failed to place segments with given parameters.")
            continue

        # ЗАКОММЕНТИРОВАННЫЙ БЛОК: создание VTK файлов
        # if args.vtk_type == 'polydata':
        #     write_vtk_polydata_lines(out_vtk, segs, header_comment=f"Generated sample {sample_i}")
        # else:
        #     # Пишем STRUCTURED_GRID — файл совместим с rect
        #     # Используем spacing_val, который определён выше (из args или conf)
        #     write_vtk_rectilinear_grid_mask(
        #         out_vtk, segs, cluster_bbox, spacing_val, thickness=0.6,
        #         header_comment=f"Generated RECTILINEAR_GRID sample {sample_i}"
        #     )
        #
        #
        # # metadata
        # meta = {
        #     'sample': sample_i,
        #     'N': args.N,
        #     'cluster_bbox': cluster_bbox,
        #     'segments': [{'Ax': float(s[0]), 'Ay': float(s[1]),
        #                   'Bx': float(s[2]), 'By': float(s[3])} for s in segs]
        # }
        # with open(out_meta, 'w', encoding='utf-8') as f:
        #     json.dump(meta, f, indent=2)
        # print(f"Wrote {out_vtk} and {out_meta}")
        
        # --- Дополнительно пишем простой TXT с координатами ---
        txt_file = os.path.splitext(out_meta)[0] + ".txt"
        with open(txt_file, 'w', encoding='utf-8') as ftxt:
            for s in segs:
                ftxt.write(f"{s[0]:.6f} {s[1]:.6f} {s[2]:.6f} {s[3]:.6f}\n")
        print(f"Wrote plain coords: {txt_file}")

if __name__ == '__main__':
    main()