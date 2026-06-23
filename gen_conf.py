#By Alena Favorskaya, Jul 2025
import numpy as np
from params import *
import os
import glob
from curved_grid_generator import print_material, print_no_reflect_filler, print_free_boundary_2D, find_N, save_to_vtk, set_dummy_fracture, gen_interpolation_file_from_rotated_rect_to_rect_2D, gen_interpolation_file_from_rect_to_rotated_rect_2D, print_interpolation_corrector_with_given_file_name
from math import sqrt, atan
import sys


background_id = "background" #one time
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
print_free_boundary_2D(1,0) #on the top is free boundary
# the source is realy close to the surface, but with condition it is in circle from source_gauss_w points the calculating grid.	
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




# folder with fractures
FRACTURES_DIR = "/home/faki/one_fracture_chimera/fractures"

txt_files = glob.glob(os.path.join(FRACTURES_DIR, "*.txt"))
filename = sys.argv[1]
first_txt = os.path.join(FRACTURES_DIR, filename)
    
fractures = []
cnt = 0
with open(first_txt, 'r') as f:
    
    for line in f:
        fracture_path = "MESHES/fracture_mesh_{:05}.vtk".format(cnt)
        fracture_line = "MESHES/fracture_line_{:05}.vtk".format(cnt)
        
        fracture_grid_id = "fracture_grid_{:05}".format(cnt)
        pos_A_x, pos_A_y, pos_B_x, pos_B_y = (float(numer) for numer in line.split())
        #our fracture is AB. A is left and (or) bottom. If the fracture is near to horisintal, then more prioritize left. If near to vertical- bottom.

        if pos_A_x >  pos_B_x:
            pos_A_x, pos_A_y, pos_B_x, pos_B_y = pos_B_x, pos_B_y, pos_A_x, pos_A_y

        fractures.append([pos_A_x, pos_A_y, pos_B_x, pos_B_y])
        
        #print(pos_A_x, pos_A_y, pos_B_x, pos_B_y)

        # gen grid with fracture
        dx = pos_B_x - pos_A_x #coordinates the fracture
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

        save_to_vtk(fracture_path, fracture_grid_x, fracture_grid_y, N_fr_x, N_fr_y) # saving grid around the fracture
        set_dummy_fracture(fracture_line, pos_A_x, pos_A_y, pos_B_x, pos_B_y) # saving segment on the place of the fracture to rendering

        #gen interpolation files

        # add indexes in names ALL files
        #names if interp files to interpolate on all diractions grid around the background grid
        interpol_file_00 = "TXT/interpol_file_00_{:05}.txt".format(cnt)
        interpol_file_01 = "TXT/interpol_file_01_{:05}.txt".format(cnt)
        interpol_file_10 = "TXT/interpol_file_10_{:05}.txt".format(cnt)
        interpol_file_11 = "TXT/interpol_file_11_{:05}.txt".format(cnt)
        #interp file to interpolate from the grid around the fracture in background grid
        interpol_file_fracture = "TXT/interpol_file_fracture_{:05}.txt".format(cnt)


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
    dx = pos_B_x - pos_A_x #coord of each of fractures is its own
    dy = pos_B_y - pos_A_y
    L = sqrt(dx*dx + dy*dy)
    N_fr = (int)(L/step + 0.5)
    step_fr = L / ((float)(N_fr))
    IF_FR_HOR = (abs(dy) < abs(dx))
    #it is contact to calculate the scattering on the fracture
    fracture_grid_id = "fracture_grid_{:05}".format(i)
    interpol_file_fracture = "TXT/interpol_file_fracture_{:05}.txt".format(i)
    interpol_file_00 = "TXT/interpol_file_00_{:05}.txt".format(i)
    interpol_file_01 = "TXT/interpol_file_01_{:05}.txt".format(i)
    interpol_file_10 = "TXT/interpol_file_10_{:05}.txt".format(i)
    interpol_file_11 = "TXT/interpol_file_11_{:05}.txt".format(i)

    print_interpolation_corrector_with_given_file_name(background_id, fracture_grid_id, 1, 0, 0, interpol_file_00)
    print_interpolation_corrector_with_given_file_name(background_id, fracture_grid_id, 1, 0, 0, interpol_file_01)
    print_interpolation_corrector_with_given_file_name(background_id, fracture_grid_id, 1, 0, 1, interpol_file_10)
    print_interpolation_corrector_with_given_file_name(background_id, fracture_grid_id, 1, 0, 1, interpol_file_11)

    print_interpolation_corrector_with_given_file_name(fracture_grid_id, background_id, 0, 1, 1, interpol_file_fracture) #corrector axis 1 it is the same predictor axis 0. In this variant calculating unknown func will be copy in background grid before the saving, not after. On the other way in the background on the place of fractire will be backlog on dt/2
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
    print("        #L = {}, dx = {}, dy = {}, N_fr = {}".format(L, dx, dy, N_fr))
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

#saving wive field, for test. For sel it's not nessessary
# print("""    [saver]
#         name = StructuredVTKSaver
#         path = VTK/%g_%s.vtk
#         order = 1""")
# print("        save = %i" %save)
# print("""        save_rank = true
#         params = v
#         norms = 1
#     [/saver]""")

# saving the sesmogram

print("""	[saver]					
		name = SeismogramSaver""")
print("		path = receivers_{}".format(filename))					
print("		r0 = %s, %s" %(receivers_start_x, receivers_y))	#points to left receiver	
print("		dr = %s, 0.0" %(receivers_step)) #step for installing receivers
print("		number = %s" %(receivers_number))			
print("		eps = %s" %(tol))				
print("		order = 2")
print("		save = %s" %(sei_save))			
print("		update = %s" %(Nt+1))
print("	[/saver]")
print("""[/savers]
""")