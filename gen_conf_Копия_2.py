#By Alena Favorskaya, Jul 2025
import numpy as np
from params import *
import os
import glob
from curved_grid_generator import print_material, print_no_reflect_filler, print_free_boundary_2D, find_N, save_to_vtk, set_dummy_fracture, gen_interpolation_file_from_rotated_rect_to_rect_2D, gen_interpolation_file_from_rect_to_rotated_rect_2D, print_interpolation_corrector_with_given_file_name
from math import sqrt, atan
import sys


#то что к трещине уборать под цикл. e end s start. 1 - 1 раз, относится НЕ к трещине

background_id = "background" #один раз
size_x_m = (Nx - 1) * step #1
size_y_m = (Ny - 1) * step #1
source_gauss_w = 5 #1
tol = step / 4.0 #1

origin_background = np.zeros((2), dtype = 'float')
origin_background[0] = 0.0 - size_x_m / 2.0
origin_background[1] = 0.0 - size_y_m

print("")
print("verbose = true")
print("")
print("dt = %s" %(dt))
print("")
print("steps = %s" %(Nt))
print("")
print("[grids]")

print("	[grid]")
print("		id = " + background_id)
print("		[node]")
print("			name = ElasticMetaNode2D")
print("""		[/node]
		[material_node]		
		[/material_node]
		[material]""")
print_material(cP, cS, rho)
print("""		[/material]
		[factory]
			name = RectGridFactory""")
print("			size = %i, %i" %(Nx, Ny))
print("			origin = %s, %s" %(origin_background[0], origin_background[1]))	
print("			spacing = %s, %s" %(step, step))	
print("			curve = false")
print("""		[/factory]
		[schema]""")
print("			name = ElasticRectSchema2DRusanov3")  
print("""		[/schema]
		[fillers]""") 
print_no_reflect_filler(0,0)	
print_no_reflect_filler(0,1)
print_no_reflect_filler(1,0)
print_no_reflect_filler(1,1)
print("""		[/fillers]
		[correctors]""")
print_free_boundary_2D(1,0) #сверху свободная граница
# источник располагается тут максимально близко к поверхности, но с условием что он размазан в круге из source_gauss_w точек расчетной сетки.	
print("			[corrector]")
print("				name = PointSourceCorrector2D")
print("				gauss_w = %i" %source_gauss_w)
print("				compression = 1.0")
print("				coords = %s, %s" %(0.0, 0.0 - (source_gauss_w + 1) * step))
print("				axis = 1")
print("				eps = %s" %(tol))            
print("""				[impulse]
					name = FileInterpolationImpulse
					[interpolator]
						name = PiceWiceInterpolator1D
						file = SOURCES/impulse.txt
					[/interpolator]
				[/impulse]
			[/corrector]""")

print("""		[/correctors]
	[/grid]""")




# Папка с трещинами
FRACTURES_DIR = "/home/faki/one_fracture_chimera/fractures"

txt_files = glob.glob(os.path.join(FRACTURES_DIR, "*.txt"))
filename = sys.argv[1]
first_txt = os.path.join(FRACTURES_DIR, filename)
    
fractures = []
cnt = 0
with open(first_txt, 'r') as f:
    
    for line in f:
        fracture_path = f"MESHES/fracture_mesh_{cnt:05}.vtk" #к трещине , добавить индекс трещины в название файла
        fracture_line = f"MESHES/fracture_line_{cnt:05}.vtk" #к трещине , добавить индекс трещины в название файла
        
        fracture_grid_id = f"fracture_grid_{cnt:05}"
        pos_A_x, pos_A_y, pos_B_x, pos_B_y = (float(numer) for numer in line.split())
        #наша трещина -- это отрезок AB. Причем точка А -- левая и (или) нижняя. Если трещина близка к горизонатали - приоритетно левая. Если трещина близка к вертикали - приоритетно нижняя. Можно этот выбор точек автоматизировать.

        if pos_A_x >  pos_B_x:
            pos_A_x, pos_A_y, pos_B_x, pos_B_y = pos_B_x, pos_B_y, pos_A_x, pos_A_y

        fractures.append([pos_A_x, pos_A_y, pos_B_x, pos_B_y])
        
        #print(pos_A_x, pos_A_y, pos_B_x, pos_B_y)

        # генерация сетки с трещиной
        dx = pos_B_x - pos_A_x #координаты у каждой трещины свои
        dy = pos_B_y - pos_A_y
        L = sqrt(dx*dx + dy*dy)
        N_fr = (int)(L/step + 0.5)
        step_fr = L / ((float)(N_fr))

        IF_FR_HOR = (abs(dy) < abs(dx))

        if IF_FR_HOR:
            N_fr_x = N_fr + 2 * N_chi + 1
            N_fr_y = 2 * N_chi + 1
            alpha = atan(dy / dx)
        else:
            N_fr_y = N_fr + 2 * N_chi + 1
            N_fr_x = 2 * N_chi + 1
            alpha = 0.0 - atan(dx / dy)
	
        #sin_alpha = sin(alpha)
        #cos_alpha = cos(alpha)

        dx = dx / ((float)(N_fr))
        dy = dy / ((float)(N_fr))

        fracture_grid_x = np.zeros((N_fr_x,N_fr_y), dtype = 'float')
        fracture_grid_y = np.zeros((N_fr_x,N_fr_y), dtype = 'float')

        if IF_FR_HOR:
            for i in range(0-N_chi, N_fr + N_chi + 1):
                for j in range(0-N_chi, N_chi + 1):
                    fracture_grid_x[N_chi+i,N_chi+j] = pos_A_x + i * dx - j * dy
                    fracture_grid_y[N_chi+i,N_chi+j] = pos_A_y + i * dy + j * dx
        else:
            for i in range(0-N_chi, N_chi + 1):
                for j in range(0-N_chi, N_fr + N_chi + 1):
                    fracture_grid_x[N_chi+i,N_chi+j] = pos_A_x + i * dy + j * dx
                    fracture_grid_y[N_chi+i,N_chi+j] = pos_A_y - i * dx + j * dy

        save_to_vtk(fracture_path, fracture_grid_x, fracture_grid_y, N_fr_x, N_fr_y) #сохранение сетки вокруг трещины
        set_dummy_fracture(fracture_line, pos_A_x, pos_A_y, pos_B_x, pos_B_y) #сохранение отрезка на месте трещины для отрисовки

        #генерация интерполяционных файлов

        #, добавить индекс трещины в названия ВСЕХ 5 тхт файлов
        #имена интерполяционных файлов для интерполяции в разные стороны сетки вокруг трещины из бэкграунд сетки
        interpol_file_00 = f"TXT/interpol_file_00_{cnt:05}.txt"
        interpol_file_01 = f"TXT/interpol_file_01_{cnt:05}.txt"
        interpol_file_10 = f"TXT/interpol_file_10_{cnt:05}.txt"
        interpol_file_11 = f"TXT/interpol_file_11_{cnt:05}.txt"
        #интерполяционный файл для интерполяции из сетки вокруг трещины в бэкграунд сетку
        interpol_file_fracture = f"TXT/interpol_file_fracture_{cnt:05}.txt"


        #grid_limits = np.zeros((6), dtype = 'int')
        #grid_limits[0] = i_start, for example 0
        #grid_limits[1] = i_end, for example Nx-1
        #grid_limits[2] = j_start, for example 0
        #grid_limits[3] = j_end, for example Ny-1
        #grid_limits[4] = i_step, if 0, equals 1
        #grid_limits[5] = j_step, if 0, equals 1
        grid_limits = np.zeros((6), dtype = 'int')
        grid_limits[0] = 0
        grid_limits[1] = N_fr_x - 1
        grid_limits[2] = 0
        grid_limits[3] = N_fr_y - 1

        eps = step / 100.0

        spacing_background = np.zeros((2), dtype = 'float')
        spacing_fracture = np.zeros((2), dtype = 'float')
        size_background = np.zeros((2), dtype = 'int')
        size_fracture = np.zeros((2), dtype = 'int')
        point = np.zeros((2), dtype = 'float')

        #point is the point of rotation
        #angle is the angle of rotation of grid_to

        point[0] = fracture_grid_x[0,0]
        point[1] = fracture_grid_y[0,0]
        spacing_background[0] = step
        spacing_background[1] = step
        spacing_fracture[0] = step_fr
        spacing_fracture[1] = step_fr
        size_background[0] = Nx
        size_background[1] = Ny
        size_fracture[0] = N_fr_x
        size_fracture[1] = N_fr_y

        gen_interpolation_file_from_rotated_rect_to_rect_2D(interpol_file_fracture, fracture_grid_id, grid_limits, eps, point, spacing_fracture, alpha, point, origin_background, size_background, spacing_background)

        grid_limits[0] = -2
        grid_limits[1] = -1
        grid_limits[2] = 0
        grid_limits[3] = N_fr_y - 1

        gen_interpolation_file_from_rect_to_rotated_rect_2D(interpol_file_00, background_id, grid_limits, eps, size_background, origin_background, spacing_background, alpha, point, point, spacing_fracture)

        grid_limits[0] = N_fr_x
        grid_limits[1] = N_fr_x + 1
        grid_limits[2] = 0
        grid_limits[3] = N_fr_y - 1

        gen_interpolation_file_from_rect_to_rotated_rect_2D(interpol_file_01, background_id, grid_limits, eps, size_background, origin_background, spacing_background, alpha, point, point, spacing_fracture)

        grid_limits[0] = 0
        grid_limits[1] = N_fr_x - 1
        grid_limits[2] = -2
        grid_limits[3] = -1

        gen_interpolation_file_from_rect_to_rotated_rect_2D(interpol_file_10, background_id, grid_limits, eps, size_background, origin_background, spacing_background, alpha, point, point, spacing_fracture)

        grid_limits[0] = 0
        grid_limits[1] = N_fr_x - 1
        grid_limits[2] = N_fr_y
        grid_limits[3] = N_fr_y + 1

        gen_interpolation_file_from_rect_to_rotated_rect_2D(interpol_file_11, background_id, grid_limits, eps, size_background, origin_background, spacing_background, alpha, point, point, spacing_fracture)

        print("	[grid]")
        print("		id = " + fracture_grid_id) #
        print("		[node]")
        print("			name = ElasticMetaNode2D")
        print("""		[/node]
        		[material_node]
        		[/material_node]
        		[material]""")
        print_material(cP, cS, rho)
        print("""		[/material]
        		[factory]""")
        print("			name = VTKGridFactory")
        print("			file = " + fracture_path)
        print("			curve = true")
        print("""		[/factory]
        		[schema]""")
        print("			name = ElasticCurveSchema2DRusanov3")
        print("""		[/schema]
        		[fillers]""")
        print("""		[/fillers]
        		[correctors]""")
        print("""		[/correctors]
        	[/grid]""")  	
            
        cnt += 1





#1 s
print("""[/grids]

[contacts]""")
for i in range(cnt):

    pos_A_x, pos_A_y, pos_B_x, pos_B_y = fractures[i][0], fractures[i][1], fractures[i][2], fractures[i][3]
    dx = pos_B_x - pos_A_x #координаты у каждой трещины свои
    dy = pos_B_y - pos_A_y
    L = sqrt(dx*dx + dy*dy)
    N_fr = (int)(L/step + 0.5)
    step_fr = L / ((float)(N_fr))
    IF_FR_HOR = (abs(dy) < abs(dx))
    #это контакт собственно для расчета рассеяния на самой трещине
    fracture_grid_id = f"fracture_grid_{i:05}"
    interpol_file_fracture = f"TXT/interpol_file_fracture_{i:05}.txt"
    interpol_file_00 = f"TXT/interpol_file_00_{i:05}.txt"
    interpol_file_01 = f"TXT/interpol_file_01_{i:05}.txt"
    interpol_file_10 = f"TXT/interpol_file_10_{i:05}.txt"
    interpol_file_11 = f"TXT/interpol_file_11_{i:05}.txt"

    print_interpolation_corrector_with_given_file_name(background_id, fracture_grid_id, 1, 0, 0, interpol_file_00)
    print_interpolation_corrector_with_given_file_name(background_id, fracture_grid_id, 1, 0, 0, interpol_file_01)
    print_interpolation_corrector_with_given_file_name(background_id, fracture_grid_id, 1, 0, 1, interpol_file_10)
    print_interpolation_corrector_with_given_file_name(background_id, fracture_grid_id, 1, 0, 1, interpol_file_11)

    print_interpolation_corrector_with_given_file_name(fracture_grid_id, background_id, 0, 1, 1, interpol_file_fracture) #корректор оси 1 то же самое по этапу расчета, что и предиктор оси 0. При этом в данном варианте посчитанное в трещине поле неизвестных функций будет копироваться в бэкграунд сетку перед сохранением, а не после. В противном случае в бэкграунд сетке на месте трещины было бы отставание на полшага интегрирования по времени.
    print_interpolation_corrector_with_given_file_name(fracture_grid_id, background_id, 1, 0, 1, interpol_file_fracture) 
    print("""    [contact]
            name = CrackCorrector2D""")
    if IF_FR_HOR:
        print("        axis = 1")
         
    else:
        print("        axis = 0")
    print("        grid1 = " + fracture_grid_id)
    print("        grid2 = " + fracture_grid_id)
    print("        tol = %s" %tol)
    if IF_FR_HOR:
        print("        bottomLeftCorner = %s, %s" %(N_chi, N_chi))
        print("        topRightCorner = %s, %s" %(N_chi + N_fr, N_chi))
    else:
        print("        bottomLeftCorner = %s, %s" %(N_chi, N_chi))
        print("        topRightCorner = %s, %s" %(N_chi, N_chi + N_fr))
    print(f"        #L = {L}, dx = {dx}, dy = {dy}, N_fr = {N_fr}")
    print("""        enableFiller = 1
            enableSaver = 0
            [schema]
                name = ElasticCurveSchema2DRusanov3
            [/schema]

            [contact]
                name = SlipRectElasticContact2D
            [/contact]

        [/contact]""")


print("""[/contacts]

[initials]
[/initials]
[savers]""")

#сохранение волнового поля. чисто для теста. для генерации выборки оно не нужно
#сохранение волнового поля. чисто для теста. для генерации выборки оно не нужно
print("""    [saver]
        name = StructuredVTKSaver
        path = VTK/%g_%s.vtk
        order = 1""")
print("        save = %i" %save)
print("""        save_rank = true
        params = v
        norms = 1
    [/saver]""")

#сохранение сейсмограммы

print("""	[saver]					
		name = SeismogramSaver""")
print(f"		path = receivers_{filename}")					
print("		r0 = %s, %s" %(receivers_start_x, receivers_y))			
print("		dr = %s, 0.0" %(receivers_step))			
print("		number = %s" %(receivers_number))			
print("		eps = %s" %(tol))				
print("		order = 2")
print("		save = %s" %(sei_save))			
print("		update = %s" %(Nt+1))
print("	[/saver]")
print("""[/savers]
""")
