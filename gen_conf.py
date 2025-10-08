#!/usr/bin/env python3
"""
gen_conf.py -- генерирует conf.conf для одного файла fractures_*.vtk

Вызов:
    python3 gen_conf.py --in fractures/fractures_1234_3tr.vtk > conf.conf

Требует: наличие рядом с .vtk мета-файла .json,
который создаёт fracture_generator_vtk.py (ключ "segments").
Использует параметры из params.py и функции из curved_grid_generator.
"""
import os
import sys
import argparse
import json
import glob
import numpy as np
from math import sqrt, atan
from params import *
from curved_grid_generator import (
    print_material, print_no_reflect_filler, print_free_boundary_2D,
    save_to_vtk, set_dummy_fracture,
    gen_interpolation_file_from_rotated_rect_to_rect_2D,
    gen_interpolation_file_from_rect_to_rotated_rect_2D,
    print_interpolation_corrector_with_given_file_name
)

def read_segments_from_meta(vtk_path):
    base = os.path.splitext(vtk_path)[0]
    meta_path = base + ".json"
    if not os.path.exists(meta_path):
        print(f"ERROR: meta .json not found for {vtk_path} (expected {meta_path})", file=sys.stderr)
        return []
    try:
        with open(meta_path, 'r', encoding='utf-8') as f:
            m = json.load(f)
    except Exception as e:
        print(f"ERROR: failed to read {meta_path}: {e}", file=sys.stderr)
        return []
    segs = []
    for s in m.get('segments', []):
        try:
            Ax = float(s['Ax']); Ay = float(s['Ay'])
            Bx = float(s['Bx']); By = float(s['By'])
            segs.append((Ax, Ay, Bx, By))
        except Exception:
            continue
    return segs

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--in', dest='in_vtk', required=True, help='input fractures_*.vtk (with .json meta)')
    p.add_argument('--out-dir', dest='out_dir', default='.', help='base dir for MESHES and TXT (default: project root)')
    args = p.parse_args()

    vtk_path = args.in_vtk
    if not os.path.exists(vtk_path):
        print(f"ERROR: input vtk not found: {vtk_path}", file=sys.stderr)
        sys.exit(2)

    out_base = args.out_dir
    MESHES_DIR = os.path.join(out_base, "MESHES")
    TXT_DIR = os.path.join(out_base, "TXT")
    os.makedirs(MESHES_DIR, exist_ok=True)
    os.makedirs(TXT_DIR, exist_ok=True)

    segments = read_segments_from_meta(vtk_path)
    if not segments:
        print(f"# Warning: no segments found in meta for {vtk_path}", file=sys.stderr)
        # proceed but no fracture grids will be printed

    # ---------------- print header and background grid ----------------
    size_x_m = (Nx - 1) * step
    size_y_m = (Ny - 1) * step
    tol = step / 4.0
    origin_background = np.zeros((2), dtype='float')
    origin_background[0] = 0.0 - size_x_m / 2.0
    origin_background[1] = 0.0 - size_y_m

    print("")
    print("verbose = true")
    print("")
    print(f"dt = {dt}")
    print("")
    print(f"steps = {Nt}")
    print("")
    print("[grids]")
    # background grid
    print("    [grid]")
    print("        id = background")
    print("        [node]")
    print("            name = ElasticMetaNode2D")
    print("        [/node]")
    print("        [material_node]")
    print("            name = ElasticMaterialMetaNode")
    print("        [/material_node]")
    print("        [material]")
    print_material(cP, cS, rho)
    print("        [/material]")
    print("        [factory]")
    print("            name = RectGridFactory")
    print(f"            size = {Nx}, {Ny}")
    print(f"            origin = {origin_background[0]}, {origin_background[1]}")
    print(f"            spacing = {step}, {step}")
    print("            curve = false")
    print("        [/factory]")
    print("        [schema]")
    print("            name = ElasticMatRectSchema2DRusanov3")
    print("        [/schema]")
    print("        [fillers]")
    print_no_reflect_filler(0,0)
    print_no_reflect_filler(0,1)
    print_no_reflect_filler(1,0)
    print_no_reflect_filler(1,1)
    print("        [/fillers]")
    print("        [correctors]")
    print_free_boundary_2D(1,0)
    source_gauss_w = 5
    print("            [corrector]")
    print("                name = PointSourceCorrector2D")
    print(f"                gauss_w = {source_gauss_w}")
    print("                compression = 1.0")
    print(f"                coords = {0.0}, {0.0 - (source_gauss_w + 1) * step}")
    print("                axis = 1")
    print(f"                eps = {tol}")
    print("                [impulse]")
    print("                    name = FileInterpolationImpulse")
    print("                    [interpolator]")
    print("                        name = PiceWiceInterpolator1D")
    print(f"                        file = TXT/impulse.txt")
    print("                    [/interpolator]")
    print("                [/impulse]")
    print("            [/corrector]")
    print("        [/correctors]")
    print("    [/grid]")

    # ---------------- build fracture grids for every segment ----------------
    contacts_to_print = []
    crack_infos = []
    fr_index = 0

    for seg in segments:
        fr_index += 1
        pos_A_x, pos_A_y, pos_B_x, pos_B_y = seg

        dx_full = pos_B_x - pos_A_x
        dy_full = pos_B_y - pos_A_y
        L = sqrt(dx_full*dx_full + dy_full*dy_full)
        N_fr = int(L/step + 0.5)
        if N_fr <= 0:
            # skip degenerate
            continue
        step_fr = L / float(N_fr)
        IF_FR_HOR = (abs(dy_full) < abs(dx_full))
        if IF_FR_HOR:
            N_fr_x = N_fr + 2 * N_chi + 1
            N_fr_y = 2 * N_chi + 1
            alpha = atan(dy_full / dx_full)
        else:
            N_fr_y = N_fr + 2 * N_chi + 1
            N_fr_x = 2 * N_chi + 1
            alpha = - atan(dx_full / dy_full)

        dx = dx_full / float(N_fr)
        dy = dy_full / float(N_fr)

        fracture_grid_x = np.zeros((N_fr_x, N_fr_y), dtype='float')
        fracture_grid_y = np.zeros((N_fr_x, N_fr_y), dtype='float')

        if IF_FR_HOR:
            for ii in range(0 - N_chi, N_fr + N_chi + 1):
                for jj in range(0 - N_chi, N_chi + 1):
                    fracture_grid_x[N_chi + ii, N_chi + jj] = pos_A_x + ii * dx - jj * dy
                    fracture_grid_y[N_chi + ii, N_chi + jj] = pos_A_y + ii * dy + jj * dx
        else:
            for ii in range(0 - N_chi, N_chi + 1):
                for jj in range(0 - N_chi, N_fr + N_chi + 1):
                    fracture_grid_x[N_chi + ii, N_chi + jj] = pos_A_x + ii * dy + jj * dx
                    fracture_grid_y[N_chi + ii, N_chi + jj] = pos_A_y - ii * dx + jj * dy

        # write local mesh and line (unique names)
        mesh_name = f"fracture_mesh_{fr_index:04d}.vtk"
        mesh_path = os.path.join(MESHES_DIR, mesh_name)
        line_name = f"fracture_line_{fr_index:04d}.vtk"
        line_path = os.path.join(MESHES_DIR, line_name)
        save_to_vtk(mesh_path, fracture_grid_x, fracture_grid_y, N_fr_x, N_fr_y)
        set_dummy_fracture(line_path, pos_A_x, pos_A_y, pos_B_x, pos_B_y)

        # unique ids and filenames
        fracture_grid_id_local = f"fracture_grid_{fr_index:04d}"
        interpol_file_fracture = os.path.join(TXT_DIR, f"interpol_file_fracture_{fr_index:04d}.txt")
        interpol_file_00 = os.path.join(TXT_DIR, f"interpol_file_{fr_index:04d}_00.txt")
        interpol_file_01 = os.path.join(TXT_DIR, f"interpol_file_{fr_index:04d}_01.txt")
        interpol_file_10 = os.path.join(TXT_DIR, f"interpol_file_{fr_index:04d}_10.txt")
        interpol_file_11 = os.path.join(TXT_DIR, f"interpol_file_{fr_index:04d}_11.txt")

        # prepare interpolation parameters and write interpolation txt files
        grid_limits = np.zeros((6), dtype='int')
        grid_limits[0] = 0
        grid_limits[1] = N_fr_x - 1
        grid_limits[2] = 0
        grid_limits[3] = N_fr_y - 1
        eps = step / 100.0
        spacing_background = np.array([step, step], dtype='float')
        spacing_fracture = np.array([step_fr, step_fr], dtype='float')
        size_background = np.array([Nx, Ny], dtype='int')
        point = np.zeros((2), dtype='float')
        point[0] = fracture_grid_x[0, 0]
        point[1] = fracture_grid_y[0, 0]

        gen_interpolation_file_from_rotated_rect_to_rect_2D(
            interpol_file_fracture,
            fracture_grid_id_local,
            grid_limits, eps, point, spacing_fracture, alpha,
            point, origin_background, size_background, spacing_background
        )

        # 00
        grid_limits[0] = -2; grid_limits[1] = -1
        grid_limits[2] = 0; grid_limits[3] = N_fr_y - 1
        gen_interpolation_file_from_rect_to_rotated_rect_2D(
            interpol_file_00, "background", grid_limits, eps, size_background,
            origin_background, spacing_background, alpha, point, point, spacing_fracture
        )
        # 01
        grid_limits[0] = N_fr_x; grid_limits[1] = N_fr_x + 1
        grid_limits[2] = 0; grid_limits[3] = N_fr_y - 1
        gen_interpolation_file_from_rect_to_rotated_rect_2D(
            interpol_file_01, "background", grid_limits, eps, size_background,
            origin_background, spacing_background, alpha, point, point, spacing_fracture
        )
        # 10
        grid_limits[0] = 0; grid_limits[1] = N_fr_x - 1
        grid_limits[2] = -2; grid_limits[3] = -1
        gen_interpolation_file_from_rect_to_rotated_rect_2D(
            interpol_file_10, "background", grid_limits, eps, size_background,
            origin_background, spacing_background, alpha, point, point, spacing_fracture
        )
        # 11
        grid_limits[0] = 0; grid_limits[1] = N_fr_x - 1
        grid_limits[2] = N_fr_y; grid_limits[3] = N_fr_y + 1
        gen_interpolation_file_from_rect_to_rotated_rect_2D(
            interpol_file_11, "background", grid_limits, eps, size_background,
            origin_background, spacing_background, alpha, point, point, spacing_fracture
        )

        # print grid block for this fracture grid
        print("    [grid]")
        print(f"        id = {fracture_grid_id_local}")
        print("        [node]")
        print("            name = ElasticMetaNode2D")
        print("        [/node]")
        print("        [material_node]")
        print("        [/material_node]")
        print("        [material]")
        print_material(cP, cS, rho)
        print("        [/material]")
        print("        [factory]")
        print("            name = VTKGridFactory")
        print(f"            file = {mesh_path}")
        print("            curve = true")
        print("        [/factory]")
        print("        [schema]")
        print("            name = ElasticCurveSchema2DRusanov3")
        print("        [/schema]")
        print("        [fillers]")
        print("        [/fillers]")
        print("        [correctors]")
        print("        [/correctors]")
        print("    [/grid]")

        # remember contacts to print later
        contacts_to_print.append(("background", fracture_grid_id_local, 1, 0, 0, interpol_file_00))
        contacts_to_print.append(("background", fracture_grid_id_local, 1, 0, 0, interpol_file_01))
        contacts_to_print.append(("background", fracture_grid_id_local, 1, 0, 1, interpol_file_10))
        contacts_to_print.append(("background", fracture_grid_id_local, 1, 0, 1, interpol_file_11))
        contacts_to_print.append((fracture_grid_id_local, "background", 0, 1, 1, interpol_file_fracture))
        contacts_to_print.append((fracture_grid_id_local, "background", 1, 0, 1, interpol_file_fracture))

        crack_infos.append({
            'grid_id': fracture_grid_id_local,
            'IF_FR_HOR': IF_FR_HOR,
            'N_fr': N_fr,
            'N_fr_x': N_fr_x,
            'N_fr_y': N_fr_y,
            'tol': tol
        })

    # end for segments

    # close grids and print contacts
    print("")
    print("[/grids]")
    print("")
    print("[contacts]")

    for it in contacts_to_print:
        bg_id, fr_id, pred_flag, corr_flag, axis_flag, fname = it
        # use helper to print the contact block
        print_interpolation_corrector_with_given_file_name(bg_id, fr_id, pred_flag, corr_flag, axis_flag, fname)

    # print crack-correctors
    for info in crack_infos:
        gid = info['grid_id']
        IFH = info['IF_FR_HOR']
        print("    [contact]")
        print("        name = CrackCorrector2D")
        if IFH:
            print("        axis = 1")
        else:
            print("        axis = 0")
        print(f"        grid1 = {gid}")
        print(f"        grid2 = {gid}")
        print(f"        tol = {info['tol']}")
        if IFH:
            print(f"        bottomLeftCorner = {N_chi}, {N_chi}")
            print(f"        topRightCorner = {N_chi + info['N_fr']}, {N_chi}")
        else:
            print(f"        bottomLeftCorner = {N_chi}, {N_chi}")
            print(f"        topRightCorner = {N_chi}, {N_chi + info['N_fr']}")
        print("        enableFiller = 1")
        print("        enableSaver = 0")
        print("        [schema]")
        print("            name = ElasticCurveSchema2DRusanov3")
        print("        [/schema]")
        print("        [contact]")
        print("            name = SlipRectElasticContact2D")
        print("        [/contact]")
        print("    [/contact]")

    # savers
    print("")
    print("[initials]")
    print("[/initials]")
    print("[savers]")
    print("    [saver]")
    print("        name = StructuredVTKSaver")
    print("        path = VTK/%g_%s.vtk")
    print("        order = 1")
    print(f"        save = {save}")
    print("        save_rank = true")
    print("        params = v")
    print("        norms = 1")
    print("    [/saver]")
    print("    [saver]")
    print("        name = SeismogramSaver")
    print("        path = receivers.csv")
    print(f"        r0 = {receivers_start_x}, {receivers_y}")
    print(f"        dr = {receivers_step}, 0.0")
    print(f"        number = {receivers_number}")
    print(f"        eps = {tol}")
    print("        order = 2")
    print(f"        save = {sei_save}")
    print(f"        update = {Nt+1}")
    print("    [/saver]")
    print("[/savers]")
    # done

if __name__ == '__main__':
    main()
дв
