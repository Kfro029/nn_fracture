#Author:  Alena V. Favorskaya, ver 30 06 2025 
#J >> K in 3D interpolation comments ??
#TODO разделить на части, чтобы не было в одном файле так много строк? 

import numpy as np 
from numpy import linalg as LA
from math import sqrt, asin, acos, sin, cos, pi, tan, exp, atan, fabs
from math import pow as power
import sys

#INNER PARAMS

def get_N_arc_length():
	return 10001 #1001 #10001 #for good arc lengt calcualtion, sholud be bigger than Nx and Ny 
def get_epsilon(): 
	return 1.0e-10 #less than small coordinate steps
def get_alpha_critical():
	return 10.0*pi/180.0 #in radians
def get_alpha_step():
	return 5.0*pi/180.0
def get_coef_arc_len_critical():
	return 0.1
def get_count_iteratons_critical():
	return 100
def get_if_save_iterations():
	return 0
def get_if_save_iterations_before_correction():
	return 0
def get_if_save_all_iterations():
	return 1
def get_need_test():
	return 1
def get_need_print_bad_cells_numbers():
	return 0
def get_interpolation_radius():
	return 3
def get_sin_critical():
	return sin(get_alpha_critical())
def get_cos_critical():
	return cos(get_alpha_critical())
def get_float_type_in_vtk_savers():
	return '>f4' #'float' #'>f4'
def get_print_grids():
	return 0

################
def save_grid_txt(path, points_x, points_y, Nx, Ny):
	f = open(path, "wt")
	f.write("%i %i \n" %(Nx, Ny))
	for j in range(Ny):
		for i in range(Nx):
			f.write("%s %s " %(points_x[i,j], points_y[i,j]))
		f.write("\n")
	f.close()
	return 1

def read_grid_txt(path):
	f = open(path, "rt")
	first = True
	for line in f:
		temp = line.split(" ")
		if first:
			first = False	
			Nx = (int)(temp[0])
			Ny = (int)(temp[1])
			points_x = np.zeros((Nx, Ny), dtype = 'float')
			points_y = np.zeros((Nx, Ny), dtype = 'float')
			j = 0
		else:
			for i in range(Nx):
				points_x[i,j] = (float)(temp[2*i])
				points_y[i,j] = (float)(temp[2*i+1])
			j = j + 1	
	f.close()
	return points_x, points_y, Nx, Ny

################

def convergence_order_ln_err(ln_err1, ln_err2, step1, step2):
    return (ln_err1 - ln_err2) / (log(step1) - log(step2))

def convergence_order(err1, err2, step1, step2):
#    return (log(err1/err2)) / (log(step1/step2)) #по варианту с делением вместо вычитания считает чуть дольше, зато может чуть точнее.
    return (log(err1) - log(err2)) / (log(step1) - log(step2))
	
##### save to vtk given array ######

#TODO еще функцию, которая рисует квадратики (треугольники, шестиугольники... ) в узлах сетки.

def save_to_vtk_line_width(width, grid_path, points_x, points_y, Nx, Ny):
	half_width = width / 2.0

	if ((Nx * Ny) == 1): #оба единица, одна точка
		coords = np.zeros((3 * 2 * 2 * 1), dtype=np.dtype('>f4'))

		coords[global_ind(0,0,2)] = points_x[0,0] - half_width
		coords[global_ind(0,0,2)+1] = points_y[0,0] - half_width

		coords[global_ind(1,0,2)] = points_x[0,0] + half_width
		coords[global_ind(1,0,2)+1] = points_y[0,0] - half_width

		coords[global_ind(0,1,2)] = points_x[0,0] - half_width
		coords[global_ind(0,1,2)+1] = points_y[0,0] + half_width

		coords[global_ind(1,1,2)] = points_x[0,0] + half_width
		coords[global_ind(1,1,2)+1] = points_y[0,0] + half_width
		
		f = write_header_vtk(grid_path + "_.vtk", 2, 2)
		coords.tofile(f)
		
		del(coords)
		f.close()
		return 1

	if (Nx > 1): #Nx == 1 -- одна вертикальная линия, и она сделается в другой части кода, а здесь ее делать не надо
	#горизонтальные линии
		coords = np.zeros((3 * Nx * 2 * 1), dtype=np.dtype('>f4'))
		for j in range(Ny):
			for i in range(Nx):
		
				if (i == 0):
					ind0 = i
				else:
					ind0 = i - 1

				if (i == (Nx-1)):
					ind1 = i
				else:
					ind1 = i + 1
			
				dx = points_x[ind1,j] - points_x[ind0,j]
				dy = points_y[ind1,j] - points_y[ind0,j]
				d = 1.0 / sqrt(dx*dx + dy*dy)
				dx = dx * d
				dy = dy * d

				coords[global_ind(i,0,Nx)] = points_x[i,j] + dy * half_width
				coords[global_ind(i,0,Nx)+1] = points_y[i,j] - dx * half_width
				coords[global_ind(i,1,Nx)] = points_x[i,j] - dy * half_width
				coords[global_ind(i,1,Nx)+1] = points_y[i,j] + dx * half_width			
			
			f = write_header_vtk(grid_path + "_ox_%i.vtk" %j, Nx, 2)
			coords.tofile(f)
			f.close()
		del(coords)

	if (Ny > 1): #Ny == 1 -- одна горизонтальная линия, и она сделана уже в другой части кода, а здесь ее делать не надо
	#вертикальные линии
		coords = np.zeros((3 * 2 * Ny * 1), dtype=np.dtype('>f4'))
		for i in range(Nx):
			for j in range(Ny):
		
				if (j == 0):
					ind0 = j
				else:
					ind0 = j - 1

				if (j == (Ny-1)):
					ind1 = j
				else:
					ind1 = j + 1		
		
				dx = points_x[i,ind1] - points_x[i,ind0]
				dy = points_y[i,ind1] - points_y[i,ind0]
				d = 1.0 / sqrt(dx*dx + dy*dy)
				dx = dx * d
				dy = dy * d

				coords[global_ind(0,j,2)] = points_x[i,j] - dy * half_width
				coords[global_ind(0,j,2)+1] = points_y[i,j] + dx * half_width
				coords[global_ind(1,j,2)] = points_x[i,j] + dy * half_width
				coords[global_ind(1,j,2)+1] = points_y[i,j] - dx * half_width			
			
			f = write_header_vtk(grid_path + "_oy_%i.vtk" %i, 2, Ny)
			coords.tofile(f)
			f.close()		
		del(coords)
	return 1
	
def global_ind(i,j,Nx):
	return 3 * (j * Nx + i)

def write_header_vtk(grid_path, Nx, Ny):
	# write header
	f = open(grid_path, "wt")
	f.write("# vtk DataFile Version 3.0\n")
	f.write("Created by Python generator written by Alena V. Favroskaya\n")
	f.write("BINARY\n")
	f.write("DATASET STRUCTURED_GRID\n")
	f.write("DIMENSIONS " + str(Nx) + " " + str(Ny) + " 1\n")
	f.write("POINTS " + str(Nx * Ny * 1) + " float\n")
	return f

#не работает так как при записи в файл пишется слишком длинный массив, ведь его размеры отдельно не подаются TODO
## coords = np.zeros(3 * Nx * Ny * 1, dtype=np.dtype('>f4'))
#def save_to_vtk_given_array(grid_path, points_x, points_y, Nx, Ny, coords):
#	# write header
#	f = open(grid_path, "wt")
#	f.write("# vtk DataFile Version 3.0\n")
#	f.write("Created by Python generator written by Alena V. Favroskaya\n")
#	f.write("BINARY\n")
#	f.write("DATASET STRUCTURED_GRID\n")
#	f.write("DIMENSIONS " + str(Nx) + " " + str(Ny) + " 1\n")
#	f.write("POINTS " + str(Nx * Ny * 1) + " float\n")

##	coords = np.zeros(3 * Nx * Ny * 1, dtype=np.dtype(get_float_type_in_vtk_savers()))

#	for i in range(Nx):
#		for j in range(Ny):
#			ind = 3 * (j * Nx + i)
#			coords[ind] = points_x[i,j]
#			coords[ind+1] = points_y[i,j]

#	coords.tofile(f)
#	del(coords)
#	f.close()
#	return 1

########## chimera_orthogonal ############

#points_x = np.zeros((Nx, Ny), dtype = 'float')
#points_y = np.zeros((Nx, Ny), dtype = 'float')
#boundary_points_x = np.zeros((N), dtype = 'float')
#boundary_points_y = np.zeros((N), dtype = 'float')

#axis - ось, условно перпендикулярная заданной границе
#axis = 0, N = Ny
#axis = 1, N = Nx
#side = 0, если заданная граница условно слева/снизу
#side = 1, если заданная граница условно справа/cверху

#N_arc_length = 100 * (N-1) + 1, рекомендовано, не менее чем, для сложной формы. чем сложнее форма, тем больше.
#points_x_to_arc_len_calc = np.zeros((N_arc_length), dtype = 'float')
#points_y_to_arc_len_calc = np.zeros((N_arc_length), dtype = 'float')
#other_points_x_to_arc_len_calc = np.zeros((N_arc_length), dtype = 'float')
#other_points_y_to_arc_len_calc = np.zeros((N_arc_length), dtype = 'float')
#arc_len = np.zeros((N_arc_length-1), dtype = 'float')

def chimera_orthogonal(step, Nx, Ny, points_x, points_y, side, axis, boundary_points_x, boundary_points_y):
	if (axis == 0):
		N = Ny
	else:
		N = Nx

	if ( (axis + side) == 1 ):
		s = 1
	else:
		s = -1
	step_x = 0.0 - s * step
	step_y = s * step	

	if (side == 0):
		ind = 0
		direction = 1
		if (axis == 0):
			final = Nx
		else:
			final = Ny
	else:
		if (axis == 0):
			ind = Nx - 1
		else:
			ind = Ny - 1
		direction = -1
		final = -1
	start = ind + direction

	for i in range(N):

		if (i == 0):
			ind0 = i
		else:
			ind0 = i - 1

		if (i == (N-1)):
			ind1 = i
		else:
			ind1 = i + 1

		dx = boundary_points_x[ind1] - boundary_points_x[ind0]
		dy = boundary_points_y[ind1] - boundary_points_y[ind0]
		d = 1.0 / sqrt(dx*dx + dy*dy)
		dx = dx * d
		dy = dy * d	

		if (axis == 0):
			points_x[ind,i] = boundary_points_x[i]
			points_y[ind,i] = boundary_points_y[i]

			k = 1
			for j in range(start, final, direction):
				points_x[j,i] = points_x[ind, i] + step_x * dy * k
				points_y[j,i] = points_y[ind, i] + step_y * dx * k
				k = k + 1 #TODO оптимизировать. можно ввести пересчет индексов для разных axis и side, и от k избавиться.

		else:
			points_x[i,ind] = boundary_points_x[i]
			points_y[i,ind] = boundary_points_y[i]

			k = 1
			for j in range(start, final, direction):
				points_x[i,j] = points_x[i, ind] + step_x * dy * k
				points_y[i,j] = points_y[i, ind] + step_y * dx * k
				k = k + 1 #TODO оптимизировать. можно ввести пересчет индексов для разных axis и side, и от k избавиться.
		
	return 1

def find_arc_len_chimera_ortho(width, axis, side, other_points_x_to_arc_len_calc, other_points_y_to_arc_len_calc, other_arc_len, diffs_x, diffs_y, points_x_to_arc_len_calc, points_y_to_arc_len_calc, N_arc_length):

	chimera_orthogonal_other_boundary(width, axis, side, N_arc_length, other_points_x_to_arc_len_calc, other_points_y_to_arc_len_calc, points_x_to_arc_len_calc, points_y_to_arc_len_calc)

	return get_total_arc_len_save_diffs(other_arc_len, diffs_x, diffs_y, other_points_x_to_arc_len_calc, other_points_y_to_arc_len_calc, N_arc_length)

def chimera_orthogonal_other_boundary(width, axis, side, N, other_boundary_points_x, other_boundary_points_y, boundary_points_x, boundary_points_y):

	if ( (axis + side) == 1 ):
		s = 1
	else:
		s = -1

	for i in range(N):

		if (i == 0):
			ind0 = i
		else:
			ind0 = i - 1

		if (i == (N-1)):
			ind1 = i
		else:
			ind1 = i + 1

		dx = boundary_points_x[ind1] - boundary_points_x[ind0]
		dy = boundary_points_y[ind1] - boundary_points_y[ind0]
		d = sqrt(dx*dx + dy*dy)
		dx = dx / d
		dy = dy / d
		
		other_boundary_points_x[i] = boundary_points_x[i] - s * dy * width
		other_boundary_points_y[i] = boundary_points_y[i] + s * dx * width

	return 1		

########## grid_size_based_Mach_number ############

#в момент времени time верхняя граница фронта плоской волны должна проходить через точку {poz_x, poz_y}

#в момент времени time конусы Маха от границ области не должны заходить в прямоугольник [grid_min_x, grid_max_x] x [grid_min_y, grid_max_y]

#wave_angle -- угол между фронтом волны и осью OX. Вектор движения волны с положительными компонентами {wave_sin, wave_cos}. Угол лежит в пределах от 0 до 45 градусов, иначе надо другую функцию использовать.

# wave_cos, wave_cos_y по идее можно было бы не давать на входе, но обычно они уже вычислены, и быстрее подать, чем вычислять заново

# c = cP * c_coef -- скорость рассматриваемой волны. единица для продольной либо cS / cP для поперечной

#int_coord_x, int_coord_y должно попадать в желаемую точку с целочисленными координатами будующей сетки, не обязательно {0,0}


def grid_size_based_Mach_number(wave_angle, wave_sin, wave_cos,  c, c_coef,  grid_min_x, grid_max_x,  grid_min_y, grid_max_y,  time, wave_length,  poz_x, poz_y,  step_x, step_y, int_coord_x, int_coord_y):

	mach_angle = atan(wave_cos/c_coef) #угол конуса Маха для волны вдоль вертикальной стенки
	mach_sin = sin(mach_angle)
	mach_cos = cos(mach_angle)

	X0, Y0, X1, Y1 = grid_size_based_Mach_number_wave_left_bottom_corner(mach_angle, mach_sin, mach_cos,  wave_angle, wave_sin, wave_cos,  c, c_coef,  grid_min_x, grid_max_x,  grid_min_y, grid_max_y,  time, wave_length,  poz_x, poz_y)
	
	#проверка, что угол волны не критический для заданного типа волны и c_coef.
	#угол поперечной волны критический, значит не надо выбирать, только "из угла"
	
	if ((c_coef - wave_sin) > 0.01):  #wave_sin <= c_coef. Угол не критический, выбирать надо.
	
		X0_along, Y1_along = cross_lines_coordinates(poz_x, poz_y,  wave_cos, 0.0 - wave_sin,   grid_min_x, grid_min_y,  mach_sin, 0.0 - mach_cos)

		if (((X1-X0) * (Y1-Y0)) > ((X1-X0_along) * (Y1_along-Y0))): #если "из угла" хуже, больше, чем "вдоль", то берем "вдоль".
			X0 = X0_along
			Y1 = Y1_along

	#проверяем, что наша область влезает тоже, а то мало ли где волна была в момент time
	X0 = min(grid_min_x, X0)
	Y0 = min(grid_min_y, Y0)
	X1 = max(grid_max_x, X1)
	Y1 = max(grid_max_y, Y1)	

	#делаем целочисленными с прибавлением шага и округлением в большую сторону. На всякий случай добавили еще по шагу сбоку.
	X0 = int_coord_x + signed_float_rounded_more(X0 - int_coord_x, step_x) - step_x
	Y0 = int_coord_y + signed_float_rounded_more(Y0 - int_coord_y, step_y) - step_y
	X1 = int_coord_x + signed_float_rounded_more(X1 - int_coord_x, step_x) + step_x
	Y1 = int_coord_y + signed_float_rounded_more(Y1 - int_coord_y, step_y) + step_y	
	
	return X0, Y0, X1, Y1

def grid_size_based_Mach_number_wave_left_bottom_corner(mach_angle, mach_sin, mach_cos,  wave_angle, wave_sin, wave_cos,  c, c_coef,  grid_min_x, grid_max_x,  grid_min_y, grid_max_y,  time, wave_length,  poz_x, poz_y):

	X1, Y0 = cross_lines_coordinates(poz_x, poz_y,  wave_cos, 0.0 - wave_sin,   grid_max_x, grid_min_y,  mach_sin, mach_cos)
	Y0 = Y0 - (c*time + wave_length)/wave_cos

	mach_angle_hor = atan(wave_sin/c_coef) #угол конуса Маха для волны вдоль горизонтальной стенки
	
	X0, Y1 = cross_lines_coordinates(poz_x, poz_y,  wave_cos, 0.0 - wave_sin,   grid_min_x, grid_max_y,  cos(mach_angle_hor), sin(mach_angle_hor))
	X0 = X0 - (c*time + wave_length)/wave_sin
	
	return X0, Y0, X1, Y1

def cross_lines_coordinates(x1, y1,  n1_x, n1_y,   x2, y2,  n2_x, n2_y):
	temp1 = n1_y * n2_x 
	temp2 = n1_x * n2_y
	temp = temp1 - temp2
	return ( x1 * temp1 - x2 * temp2 + (y2 - y1) * n1_x * n2_x ) / temp, ( y2 * temp1 - y1 * temp2 + (x1 - x2) * n1_y * n2_y ) / temp

def read_time_from_file(path):
	f = open(path, 'rt')
	for l in f:
		n = 0
		line_for_split = ""
		for c in l:
			if ((n==0)and(c=='m')):
				n = 1
			elif ((n==1)and(c=='a')):
				n = 2
			elif ((n==2)and(c=='i')):
				n = 3
			elif ((n==3)and(c=='n')):
				n = 4
			elif (n==4):
				line_for_split = line_for_split + c
		if (n==4):
			temp = line_for_split.split("	")
	f.close()
	return temp[2]

######################

def hat_boundary(t, center_start, amplitude, left, right, center_length):
	temp = t - center_start 
	if (temp < 0.0):
		if (left == 0.0):
			return 0.0
		temp = fabs(temp)/left #left
	else: 
		temp = temp - center_length
		if (temp <= 0.0):
			return amplitude #center
		if (right == 0.0):
			return 0.0
		temp = temp/right #right
	return amplitude * fi_hat_half_normalized(temp)

def float_data_from_file(file_name, row_index_I, column_index_J):
	f = open(file_name, "rt")
	i = 0
	for line in f:
		if (i == row_index_I):
			temp = line.split("	")
			result = (float)(temp[column_index_J])
			f.close()
			return result
		i = i + 1
	f.close()
	return "NOT FIND"

def param_to_file(file_id, name, value, if_int, ind):
	if (ind == -1):
		if if_int:
			file_id.write(name + "	%i\n" %value)
			return -1
		file_id.write(name + "	%s\n" %value)
		return -1
	if if_int:
		file_id.write("%i	" + name + "	%i\n" %(ind, value))
		return ind + 1
	file_id.write("%i	" + name + "	%s\n" %(ind, value))
	return ind + 1

def param_to_python_file(file_id, name, value, if_int):
	if if_int:
		file_id.write(name + " = %i\n" %value)
		return 1
	file_id.write(name + " = %s\n" %value)
	return 1

# (save_count < 0) => then save = 1
def find_save(Nt, save_count):
	if (save_count < 0):
		return 1
	return (int)( ((float)(Nt-1) / (float)(save_count)) + 0.5)

#plus 0.5 or 0.1 or 0.0 ...
def find_save_and_correct_Nt(Nt, save_count, plus):
	if (save_count < 0):
		return Nt, 1
	save = (int)( ((float)(Nt-1) / (float)(save_count)) + plus)
	return ((int)(save_count)) * Nt + 1, save

def print_free_boundary_2D(axis, side):
	print("""			[corrector]
				name = ForceRectElasticBoundary2D""")
	print("				axis = %i" %axis)
	print("				side = %i" %side)
	print("			[/corrector]")
	return 1

def print_free_boundary_3D(axis, side):
	print("""			[corrector]
				name = ForceRectElasticBoundary3D""")
	print("				axis = %i" %axis)
	print("				side = %i" %side)
	print("			[/corrector]")
	return 1
			
def cut_array_2D(Nx_to, Ny_to, array_to, array_from, origin_to_x, origin_to_y, origin_from_x, origin_from_y):
	for i in range(Nx_to):
		for j in range(Ny_to):
			array_to[i, j] = array_from[i + origin_to_x - origin_from_x, j + origin_to_y - origin_from_y]
	return 1

def cut_array_3D(Nx_to, Ny_to, Nz_to, array_to, array_from, origin_to_x, origin_to_y, origin_to_z, origin_from_x, origin_from_y, origin_from_z):
	for i in range(Nx_to):
		for j in range(Ny_to):
			for k in range(Nz_to):
				array_to[i, j, k] = array_from[i + origin_to_x - origin_from_x, j + origin_to_y - origin_from_y, k + origin_to_z - origin_from_z]
	return 1


def print_material(cP, cS, rho):
	print("			c1 = %s" %cP)	
	print("			c2 = %s" %cS)
	print("			rho = %s" %rho)
	return 1
	
def param_to_file(file_id, name, value, if_int, ind):
	if (ind == -1):
		if if_int:
			file_id.write("name	%i\n" %value)
			return -1
		file_id.write("name	%s\n" %value)
		return -1
	if if_int:
		file_id.write("%i	name	%i\n" %(ind, value))
		return ind + 1
	file_id.write("%i	name	%s\n" %(ind, value))
	return ind + 1

#IfPrint = False for debug
#ContactAxis == -1 if without axes and sides
def print_contact_2D(ContactType, ContactName1, ContactName2, ContactAxis, ContactSide1, ContactTol, IfRegular, IfPrint):
	if IfPrint:
		print("	[contact]")
		if (ContactType == "GLUE"):
			print("		name = GlueRectElasticContact2D")
		else: #SLIP
			print("		name = SlipRectElasticContact2D")
		print("		grid1 = "+ContactName1)
		print("		grid2 = "+ContactName2)
		print("		tol = %s" %(ContactTol))
		print("		currect_bounds = true")
		if IfRegular:
			print("		optimize_rect_grids = true")
		if (ContactAxis == -1):
			print("	[/contact]")
			return 1
		print("		axis1 = %i" %ContactAxis)
		print("		axis2 = %i" %ContactAxis)
		print("		side1 = %i" %ContactSide1)
		print("		side2 = %i" %(abs(ContactSide1-1)))
		print("	[/contact]")
	return 1
	
def print_contact_3D(ContactType, ContactName1, ContactName2, ContactAxis, ContactSide1, ContactTol, IfRegular, IfPrint):
	if IfPrint:
		print("	[contact]")
		if (ContactType == "GLUE"):
			print("		name = GlueRectElasticContact3D")
		else: #SLIP
			print("		name = SlipRectElasticContact3D")
		print("		grid1 = "+ContactName1)
		print("		grid2 = "+ContactName2)
		print("		tol = %s" %(ContactTol))
		print("		currect_bounds = true")
		if IfRegular:
			print("		optimize_rect_grids = true")
		if (ContactAxis == -1):
			print("	[/contact]")
			return 1
		print("		axis1 = %i" %ContactAxis)
		print("		axis2 = %i" %ContactAxis)
		print("		side1 = %i" %ContactSide1)
		print("		side2 = %i" %(abs(ContactSide1-1)))
		print("	[/contact]")
	return 1

def move_limits_2D(limits_new, limits_old, delta_i, delta_j, origin_i_new, origin_j_new, origin_i_old, origin_j_old):
	for i in range(4):
		if (i<2):
			limits_new[i] = limits_old[i] + origin_i_old - origin_i_new + delta_i
		else:
			limits_new[i] = limits_old[i] + origin_j_old - origin_j_new + delta_j
	return 1

def move_limits_3D(limits_new, limits_old, delta_i, delta_j, delta_k, origin_i_new, origin_j_new, origin_k_new, origin_i_old, origin_j_old, origin_k_old):
	for i in range(6):
		if (i<2):
			limits_new[i] = limits_old[i] + origin_i_old - origin_i_new + delta_i
		elif (i<4):
			limits_new[i] = limits_old[i] + origin_j_old - origin_j_new + delta_j
		else:
			limits_new[i] = limits_old[i] + origin_k_old - origin_k_new + delta_k
	return 1

def find_eff_step_min_max_average(points_x, points_y, Nx, Ny):
	first = True
	sum_steps = 0.0
	for i in range(Nx):
		for j in range(Ny):
			left_x = i - 1
			right_x = i + 1
			left_y = j - 1
			right_y = j + 1
			coef_x = 0.5
			coef_y = 0.5
			if (left_x < 0):
				left_x = 0
				coef_x = 1.0
			if (left_y < 0):
				left_y = 0
				coef_y = 1.0
			if (right_x > (Nx-1)):
				right_x = Nx-1
				coef_x = 1.0
			if (right_y > (Ny-1)):
				right_y = Ny-1
				coef_y = 1.0
			delta_X_x = coef_x * (points_x[right_x,j] - points_x[left_x,j])
			delta_X_y = coef_x * (points_y[right_x,j] - points_y[left_x,j])
			delta_Y_x = coef_y * (points_x[i,right_y] - points_x[i,left_y])
			delta_Y_y = coef_y * (points_y[i,right_y] - points_y[i,left_y])

			temp = abs(delta_X_x * delta_Y_y - delta_X_y * delta_Y_x)
#			eps = sqrt(delta_Y_x * delta_Y_x + delta_Y_y * delta_Y_y) / temp
#			eta = sqrt(delta_X_x * delta_X_x + delta_X_y * delta_X_y) / temp

#			if (eps > eta):
#				temp_min = 1.0 / eps
#				temp_max = 1.0 / eta
#			else:
#				temp_min = 1.0 / eta
#				temp_max = 1.0 / eps

			temp_X = delta_Y_x*delta_Y_x + delta_Y_y*delta_Y_y
			temp_Y = delta_X_x*delta_X_x + delta_X_y*delta_X_y

			if (temp_X < temp_Y):
				temp_max = temp_X
				temp_min = temp_Y
			else:
				temp_max = temp_Y
				temp_min = temp_X
			temp_min = temp / sqrt(temp_min)
			temp_max = temp / sqrt(temp_max)
			sum_steps = sum_steps + temp_min + temp_max

			if first:
				result_min = temp_min
				result_max = temp_max
				first = False
			else:
				result_min = min(result_min, temp_min)
				result_max = max(result_max, temp_max)

	return result_min, result_max, sum_steps / ((float)(2*Nx*Ny))

#rect/structured/src/rectgrid.cpp
#void RectGrid::calculateJacobian2D()
#{
#	vector2 *eps = getInvJ<vector2>(X);
#	vector2 *eta = getInvJ<vector2>(Y);
#	for (int_t j = 0; j < localSize(Y); j++)
#	for (int_t i = 0; i < localSize(X); i++) {
#		int_t ind = localIndex(i, j);
#		vector2 x = delta<vector2>(ind, X);
#		vector2 y = delta<vector2>(ind, Y);
#		real_t J_1 =  x.x * y.y - x.y * y.x;
#
#		eps[ind].x = y.y;
#		eps[ind].y = -y.x;
#		eps[ind] /= J_1;
#
#		eta[ind].x = -x.y;
#		eta[ind].y = x.x;
#		eta[ind] /= J_1;
#
#		if (i == 0 && j == 0) {
#			jmin[X] = jmax[X] = norm(eps[ind]);
#			jmin[Y] = jmax[Y] = norm(eta[ind]);
#		} else {
#			const real_t neps = norm(eps[ind]);
#			const real_t neta = norm(eta[ind]);
#			if (jmin[X] > neps) jmin[X] = neps;
#			if (jmax[X] < neps) jmax[X] = neps;
#			if (jmin[Y] > neta) jmin[Y] = neta;
#			if (jmax[Y] < neta) jmax[Y] = neta;
#		}
#	}
#	cout << "RectGrid::calculateJacobian2D: Min/max eps = " << jmin[X] << "/" << jmax[X] << endl;

def find_eff_step(points_x, points_y, Nx, Ny):
	first = True
	for i in range(Nx):
		for j in range(Ny):
			left_x = i - 1
			right_x = i + 1
			left_y = j - 1
			right_y = j + 1
			coef_x = 0.5
			coef_y = 0.5
			if (left_x < 0):
				left_x = 0
				coef_x = 1.0
			if (left_y < 0):
				left_y = 0
				coef_y = 1.0
			if (right_x > (Nx-1)):
				right_x = Nx-1
				coef_x = 1.0
			if (right_y > (Ny-1)):
				right_y = Ny-1
				coef_y = 1.0
			delta_X_x = coef_x * (points_x[right_x,j] - points_x[left_x,j])
			delta_X_y = coef_x * (points_y[right_x,j] - points_y[left_x,j])
			delta_Y_x = coef_y * (points_x[i,right_y] - points_x[i,left_y])
			delta_Y_y = coef_y * (points_y[i,right_y] - points_y[i,left_y])
			temp = delta_Y_x*delta_Y_x + delta_Y_y*delta_Y_y
			temp = max(temp, delta_X_x*delta_X_x + delta_X_y*delta_X_y)
			temp = abs(delta_X_x * delta_Y_y - delta_X_y * delta_Y_x) / sqrt(temp)
			if first:
				result = temp
				first = False
			else:
				result = min(result, temp)
	return result

def find_eff_step_3D(points_x, points_y, points_z, Nx, Ny, Nz):
	first = True
	for i in range(Nx):
		for j in range(Ny):
			for k in range(Nz):
				left_x = i - 1
				right_x = i + 1

				left_y = j - 1
				right_y = j + 1

				left_z = k - 1
				right_z = k + 1
				
				coef_x = 0.5
				coef_y = 0.5
				coef_z = 0.5

				if (left_x < 0):
					left_x = 0
					coef_x = 1.0
				if (left_y < 0):
					left_y = 0
					coef_y = 1.0
				if (left_z < 0):
					left_z = 0
					coef_z = 1.0
					
				if (right_x > (Nx-1)):
					right_x = Nx-1
					coef_x = 1.0
				if (right_y > (Ny-1)):
					right_y = Ny-1
					coef_y = 1.0
				if (right_z > (Nz-1)):
					right_z = Nz-1
					coef_z = 1.0
					
				ax = coef_x * (points_x[right_x,j,k] - points_x[left_x,j,k])
				ay = coef_x * (points_y[right_x,j,k] - points_y[left_x,j,k])
				az = coef_x * (points_z[right_x,j,k] - points_z[left_x,j,k])

				bx = coef_y * (points_x[i,right_y,k] - points_x[i,left_y,k])
				by = coef_y * (points_y[i,right_y,k] - points_y[i,left_y,k])
				bz = coef_y * (points_z[i,right_y,k] - points_z[i,left_y,k])

				cx = coef_z * (points_x[i,j,right_z] - points_x[i,j,left_z])
				cy = coef_z * (points_y[i,j,right_z] - points_y[i,j,left_z])
				cz = coef_z * (points_z[i,j,right_z] - points_z[i,j,left_z])

#				real_t J_1 = mixed_product(x, y, z);
#				eps[ind] = outer_product(y, z) / J_1;
#				eta[ind] = outer_product(z, x) / J_1;
#				ksi[ind] = outer_product(x, y) / J_1;
#					qmin = qmax = vector3(norm(eps[ind]), norm(eta[ind]), norm(ksi[ind]));

				temp = max(norm2_vector_product_3D(ax, ay, az, bx, by, bz), norm2_vector_product_3D(bx, by, bz, cx, cy, cz), norm2_vector_product_3D(cx, cy, cz, ax, ay, az)) 
				temp = abs(ax*by*cz + ay*bz*cx + az*bx*cy - cx*by*az - cy*bz*ax - cz*bx*ay) / sqrt(temp)
				
				if first:
					result = temp
					first = False
				else:
					result = min(result, temp)
	return result

def norm2_vector_product_3D(ax, ay, az, bx, by, bz):
	x = ay * bz - by * az
	y = az * bx - bz * ax
	z = ax * by - bx * ay
	return x*x + y*y + z*z

def find_min(array, N):
	res = array[0]
	for i in range(1,N):
		res = min(res, array[i])
	return res

def find_max(array, N):
	res = array[0]
	for i in range(1,N):
		res = max(res, array[i])
	return res

def IFN(folder_name, grid_from, grid_to, axis, side): #Interpolation File Name
	if (axis > -1):
		return folder_name + "/from_" + grid_from + "_to_" + grid_to + "_axis_%i_side_%i_.txt" %(axis, side)
	return folder_name + "/from_" + grid_from + "_to_" + grid_to + ".txt"

def print_interpolation_correctors_for_fracture_grid_with_dop_grids_2D(folder_name, gnM, gnF, fracture_axis, contact_type, tol):
	print_interpolation_correctors_for_fracture_grid_2D(folder_name, gnM, gnF)

	gn_mini_0 = gnF + "_A"
	gn_mini_1 = gnF + "_B"

	if (fracture_axis == 0):
		contact_axis = 1
	else:
		contact_axis = 0

	#from the main curved grid to the bottom/left grid
	print_interpolation_corrector_with_given_file_name(gnF, gn_mini_0, 1, 0, contact_axis, IFN(folder_name, gnF, gn_mini_0, -1, -1))
	#from the main curved grid to the top/right grid
	print_interpolation_corrector_with_given_file_name(gnF, gn_mini_1, 1, 0, contact_axis, IFN(folder_name, gnF, gn_mini_1, -1, -1))

	#from the the bottom/left grid to the main curved grid 
	print_interpolation_corrector_with_given_file_name(gn_mini_0, gnF, 1, 0, contact_axis, IFN(folder_name, gn_mini_0, gnF, -1, -1))
	#from the top/right grid to the main curved grid
	print_interpolation_corrector_with_given_file_name(gn_mini_1, gnF, 1, 0, contact_axis, IFN(folder_name, gn_mini_1, gnF, -1, -1))
	
	if (contact_type == "Slip"):
		print("	[contact]")
		print("		name = SlipRectElasticContact2D")
		print("		grid1 = " + gn_mini_1)
		print("		grid2 = " + gn_mini_0)
		print("		tol = %s" %tol)
		print("		axis1 = %i" %contact_axis)
		print("		axis2 = %i" %contact_axis)
		print("		side1 = 0")
		print("		side2 = 1")
		print("	[/contact]")
	elif (contact_type == "Gas"):
		sys.exit("Gas TODO as free boundaries")

	return 1

def print_interpolation_correctors_for_fracture_grid_2D(folder_name, gnM, gnF):
	print_interpolation_correctors_to_small_rect_grid_2D(folder_name, gnM, gnF, 1, 0)
	print_interpolation_corrector_with_given_file_name(gnF, gnM, 1, 0, 1, IFN(folder_name, gnF, gnM, -1, -1))
	print_interpolation_corrector_with_given_file_name(gnF, gnM, 0, 1, 1, IFN(folder_name, gnF, gnM, -1, -1))
	return 1

def print_interpolation_correctors_to_small_rect_grid_2D(folder_name, gnM, gnF, predictor_flag, corrector_flag):
	print_interpolation_corrector_with_given_file_name(gnM, gnF, predictor_flag, corrector_flag, 0, IFN(folder_name, gnM, gnF, 0, 0))
	print_interpolation_corrector_with_given_file_name(gnM, gnF, predictor_flag, corrector_flag, 0, IFN(folder_name, gnM, gnF, 0, 1))
	print_interpolation_corrector_with_given_file_name(gnM, gnF, predictor_flag, corrector_flag, 1, IFN(folder_name, gnM, gnF, 1, 0))
	print_interpolation_corrector_with_given_file_name(gnM, gnF, predictor_flag, corrector_flag, 1, IFN(folder_name, gnM, gnF, 1, 1))
	return 1
	
def print_interpolation_correctors_for_fracture_grid_3D(folder_name, gnM, gnF):
	print_interpolation_correctors_to_small_rect_grid_3D(folder_name, gnM, gnF, 1, 0)
	print_interpolation_corrector_with_given_file_name(gnF, gnM, 1, 0, 1, IFN(folder_name, gnF, gnM, -1, -1))
	print_interpolation_corrector_with_given_file_name(gnF, gnM, 1, 0, 2, IFN(folder_name, gnF, gnM, -1, -1))
	print_interpolation_corrector_with_given_file_name(gnF, gnM, 0, 1, 2, IFN(folder_name, gnF, gnM, -1, -1))
	return 1

def print_interpolation_correctors_to_small_rect_grid_3D(folder_name, gnM, gnF, predictor_flag, corrector_flag):
	print_interpolation_corrector_with_given_file_name(gnM, gnF, predictor_flag, corrector_flag, 0, IFN(folder_name, gnM, gnF, 0, 0))
	print_interpolation_corrector_with_given_file_name(gnM, gnF, predictor_flag, corrector_flag, 0, IFN(folder_name, gnM, gnF, 0, 1))
	print_interpolation_corrector_with_given_file_name(gnM, gnF, predictor_flag, corrector_flag, 1, IFN(folder_name, gnM, gnF, 1, 0))
	print_interpolation_corrector_with_given_file_name(gnM, gnF, predictor_flag, corrector_flag, 1, IFN(folder_name, gnM, gnF, 1, 1))
	print_interpolation_corrector_with_given_file_name(gnM, gnF, predictor_flag, corrector_flag, 2, IFN(folder_name, gnM, gnF, 2, 0))
	print_interpolation_corrector_with_given_file_name(gnM, gnF, predictor_flag, corrector_flag, 2, IFN(folder_name, gnM, gnF, 2, 1))
	return 1

def gen_interpolation_files_for_fracture_grid_with_dop_grids_2D(folder_name, gnM, gnF, array_corners):
	gen_interpolation_files_for_fracture_grid_2D(folder_name, gnM, gnF, array_corners)
	
	grid_to_limits = np.zeros((6), dtype = 'int')

	#from the main curved grid to the bottom/left grid
	gn_mini = gnF + "_A"

	if (array_corners[11] == 0):
		grid_to_limits[0] = 0
		grid_to_limits[1] = array_corners[2] - array_corners[0]
		grid_to_limits[2] = 0
		grid_to_limits[3] = 1

		def Func_I(i,j):
			return i + array_corners[0]
		
		def Func_J(i,j):
			return j - 4 + array_corners[1]
		
	else:
		grid_to_limits[0] = 0
		grid_to_limits[1] = 1
		grid_to_limits[2] = 0
		grid_to_limits[3] = array_corners[3] - array_corners[1]

		def Func_I(i,j):
			return i - 4 + array_corners[0]
		
		def Func_J(i,j):
			return j + array_corners[1]
		
	gen_interpolation_file_copy_rect_2D(IFN(folder_name, gnF, gn_mini, -1, -1), gnF, grid_to_limits, Func_I, Func_J)

	#from the main curved grid to the top/right grid
	gn_mini = gnF + "_B"

	def Func_I(i,j):
		return i + array_corners[0]
		
	def Func_J(i,j):
		return j + array_corners[1]

	if (array_corners[11] == 0):
		grid_to_limits[0] = 0
		grid_to_limits[1] = array_corners[2] - array_corners[0]
		grid_to_limits[2] = 3
		grid_to_limits[3] = 4

	else:
		grid_to_limits[0] = 3
		grid_to_limits[1] = 4
		grid_to_limits[2] = 0
		grid_to_limits[3] = array_corners[3] - array_corners[1]
	
	gen_interpolation_file_copy_rect_2D(IFN(folder_name, gnF, gn_mini, -1, -1), gnF, grid_to_limits, Func_I, Func_J)

	#from the bottom/left grid to the main curved grid
	gn_mini = gnF + "_A"

	if (array_corners[11] == 0):
		grid_to_limits[0] = array_corners[0]
		grid_to_limits[1] = array_corners[2]
		grid_to_limits[2] = array_corners[1] - 2
		grid_to_limits[3] = array_corners[3] - 1

		def Func_I(i,j):
			return i - array_corners[0]
		
		def Func_J(i,j):
			return j + 4 - array_corners[1]
		
	else:
		grid_to_limits[0] = array_corners[0] - 2
		grid_to_limits[1] = array_corners[2] - 1
		grid_to_limits[2] = array_corners[1]
		grid_to_limits[3] = array_corners[3]

		def Func_I(i,j):
			return i + 4 - array_corners[0]
		
		def Func_J(i,j):
			return j - array_corners[1]

	gen_interpolation_file_copy_rect_2D(IFN(folder_name, gn_mini, gnF, -1, -1), gn_mini, grid_to_limits, Func_I, Func_J)

	#from the top/right grid to the main curved grid
	gn_mini = gnF + "_B"

	if (array_corners[11] == 0):
		grid_to_limits[0] = array_corners[0]
		grid_to_limits[1] = array_corners[2]
		grid_to_limits[2] = array_corners[1]
		grid_to_limits[3] = array_corners[3] + 2

	else:
		grid_to_limits[0] = array_corners[0]
		grid_to_limits[1] = array_corners[2] + 2
		grid_to_limits[2] = array_corners[1]
		grid_to_limits[3] = array_corners[3]

	def Func_I(i,j):
		return i - array_corners[0]
		
	def Func_J(i,j):
		return j - array_corners[1]
	
	gen_interpolation_file_copy_rect_2D(IFN(folder_name, gn_mini, gnF, -1, -1), gn_mini, grid_to_limits, Func_I, Func_J)
	
	return 1

def gen_interpolation_files_for_fracture_grid_2D(folder_name, gnM, gnF, array_corners):
	gen_interpolation_files_to_small_rect_grid_2D(folder_name, gnM, gnF, array_corners[6], array_corners[7], array_corners[8], array_corners[9])

#plus 
#from fracture to medium
	grid_to_limits = np.zeros((6), dtype = 'int')

	grid_to_limits[0] = array_corners[8]
	grid_to_limits[1] = array_corners[8] + array_corners[6] - 1
	grid_to_limits[2] = array_corners[9]
	grid_to_limits[3] = array_corners[9] + array_corners[7] - 1

	def Func_I(i,j):
		return i - array_corners[8]
		
	def Func_J(i,j):
		return j - array_corners[9]

	gen_interpolation_file_copy_rect_2D(IFN(folder_name, gnF, gnM, -1, -1), gnF, grid_to_limits, Func_I, Func_J)

	return 1

#gnM - big grid (medium), grid name
#gnF - small grid (for example, curve fracture grid, origin - array_corners 8 and 9, size - 6 and 7.)

def gen_interpolation_files_to_small_rect_grid_2D(folder_name, gnM, gnF, size_small_x, size_small_y, int_origin_small_x, int_origin_small_y):
#from medium to fracture, axis 0, side 0
#from medium to fracture, axis 0, side 1
#from medium to fracture, axis 1, side 0
#from medium to fracture, axis 1, side 1
	
###

#grid_to_limits[0] = i_start, for example 0
#grid_to_limits[1] = i_end, for example Nx-1
#grid_to_limits[2] = j_start, for example 0
#grid_to_limits[3] = j_end, for example Ny-1
#grid_to_limits[4] = i_step, if 0, equals 1
#grid_to_limits[5] = j_step, if 0, equals 1

#Func_I(i,j) returns I in grid_from, if i and j in grid_to are known
#Func_J(i,j) returns J in grid_from, if i and j in grid_to are known

	grid_to_limits = np.zeros((6), dtype = 'int')

#from medium to fracture, axis 0, side 0

	grid_to_limits[0] = -2
	grid_to_limits[1] = -1
	grid_to_limits[3] = size_small_y - 1

	def Func_I(i,j):
		return i + int_origin_small_x
		
	def Func_J(i,j):
		return j + int_origin_small_y

	gen_interpolation_file_copy_rect_2D(IFN(folder_name, gnM, gnF, 0, 0), gnM, grid_to_limits, Func_I, Func_J)

#from medium to fracture, axis 0, side 1

	grid_to_limits[0] = size_small_x
	grid_to_limits[1] = size_small_x + 1

	gen_interpolation_file_copy_rect_2D(IFN(folder_name, gnM, gnF, 0, 1), gnM, grid_to_limits, Func_I, Func_J)

#from medium to fracture, axis 1, side 0

	grid_to_limits[0] = 0
	grid_to_limits[1] = size_small_x - 1
	grid_to_limits[2] = -2
	grid_to_limits[3] = -1

	gen_interpolation_file_copy_rect_2D(IFN(folder_name, gnM, gnF, 1, 0), gnM, grid_to_limits, Func_I, Func_J)

#from medium to fracture, axis 1, side 1

	grid_to_limits[2] = size_small_y
	grid_to_limits[3] = size_small_y + 1

	gen_interpolation_file_copy_rect_2D(IFN(folder_name, gnM, gnF, 1, 1), gnM, grid_to_limits, Func_I, Func_J)

	return 1

def print_no_reflect_filler(axis, side):
	print("""			[filler]
				name = RectNoReflectFiller""")
	print("				axis = %i" %axis)
	print("				side = %i" %side)
	print("			[/filler]")		
	return 1

def print_free_boundary(axis, side):
	print("""		[corrector]
				name = ForceRectElasticBoundary2D""")
	print("				axis = %i" %axis)
	print("				side = %i" %side)
	print("			[/corrector]")		
	return 1

def find_N(length, step):
	return (int)(length/step + 0.5) + 1

def find_int_rounded(length, step):
	return (int)(length/step + 0.5)

def find_N_and_correct_float_rounded(length, step):
	n = (int)(length/step + 0.5)
	return n + 1, n * step

def find_int_rounded_more_int_only(length, step):
	n = (int)(length/step + 0.5)
	if ((n * step) < length):
		n = n + 1
	return n

def find_int_rounded_more(length, step):
	n = (int)(length/step + 0.5)
	if ((n * step) < length):
		n = n + 1
	return n, n * step

def find_int_rounded_less(length, step):
	n = (int)(length/step + 0.5)
	if ((n * step) > length):
		n = n - 1
	return n, n * step

def float_rounded_more(length, step):
	n = (int)(length/step + 0.5)
	if ((n * step) < length):
		n = n + 1
	return n * step
	
def signed_float_rounded_more(distance, step):
	if (distance < 0.0):
		return 0.0 - float_rounded_more(abs(distance), step)
	return float_rounded_more(abs(distance), step)

def signed_float_rounded_int(distance, step):
	if (distance < 0.0):
		return 0.0 - find_int(abs(distance), step) * step
	return find_int(abs(distance), step) * step

def find_int(length, step):
	return (int)(length/step + 0.001)
	
def float_rounded(length, step):
	return ((int)(length/step + 0.5)) * step

def float_rounded_less(length, step):
	n = (int)(length/step + 0.5)
	if ((n>0) and ((n * step) > length)):
		n = n - 1
	return n * step

def find_int_N(length, step):
	return find_int(length, step) + 1
	
#Nt, save, real_time = find_Nt_and_save(time, dt, save_count)
def find_Nt_and_save(time, dt, save_count):
	save = (int)( time/(save_count*dt) + 0.5)
	Nt = save * (int)(save_count)
	return Nt + 1, save, Nt * dt
	
def find_step(length, N):
	return length / ((float)(N-1))

def novel_file(filename):
	try:
		f = open(filename)
		f.close()
	except FileNotFoundError:
		return True
	return False

def find_rounded_time_step_from_jacobian(jacobian, courant, dop_order):
	return find_rounded_value(courant/jacobian, dop_order)

def find_rounded_value(value, dop_order):
	dt = value
	ten = 1
	order = 0
	temp = dt
	
	def bad_temp(temp):
		temp_int = (int)(temp)
		if ((temp_int > 0) and (temp_int < 10)):
			return False
		return True
	
	def finish(dt, ten, key_div, order, dop_order):
		if key_div:
			curr_order = dop_order - order
			if (curr_order == 0):
				return (float)((int)(dt + 0.5))
			if (curr_order > 0):
				cur_ten = 10
				while i in range(curr_order - 1):
					cur_ten = cur_ten * 10
				curr_key_div = False
			else: #curr_order < 0
				cur_ten = 10
				abs_curr_order = 0 - curr_order
				for i in range(abs_curr_order - 1):
					cur_ten = cur_ten * 10
				curr_key_div = True
		else:
			curr_order = dop_order + order
			cur_ten = ten
			for i in range(dop_order):
				cur_ten = cur_ten * 10
			curr_key_div = False
			
		if curr_key_div:
			return (float)(((int)(dt / ((float)(cur_ten)) + 0.5)) * cur_ten)
		return ((float)((int)(dt * cur_ten + 0.5))) / ((float)(cur_ten))
			
	if not bad_temp(temp):
		dt = finish(dt, ten, False, order, dop_order)
		return dt
	
	if (((int)(temp)) == 0):
		key_div = False
	else:
		key_div = True
	
	while bad_temp(temp):
		ten = ten * 10
		order = order + 1
		if key_div:
			temp = dt/(float)(ten)
		else:
			temp = dt*ten

	dt = finish(dt, ten, key_div, order, dop_order)
	return dt

#TODO
#def find_rounded_more_value(value, dop_order):

def sign(x):
	if (x == 0.0):
		return 0
	if (x > 0.0):
		return 1
	return -1

def smooth_curve_grid_2D(path, Nx, Ny, ind_critical, ind_regular_start, ind_regular_end, hx, hy, func_y, y_start, y_end, eps_y):
	grid = np.zeros((Nx, Ny, 2), dtype = 'float')
	line = np.zeros((Ny), dtype = 'float')

#	N_step_bottom = ind_critical - ind_regular_start
#	N_step_top = N - ind_critical - ind_regular_end - 1
#	step_bottom = np.zeros((N_step_bottom), dtype = 'float')
#	step_top = np.zeros((N_step_top), dtype = 'float')

#here to save RAM, allocating np array one time per function
	N_step_bottom = ind_critical - ind_regular_start
	N_step_top = Ny - ind_critical - ind_regular_end - 1
	step_bottom = np.zeros((N_step_bottom), dtype = 'float')
	step_top = np.zeros((N_step_top), dtype = 'float')

	for i in range(Nx):
		grid[i,0,0] = i*hx
		for j in range(1,Ny):
			grid[i,j,0] = grid[i,0,0]
		get_smooth_line(line, Ny, ind_critical, ind_regular_start, ind_regular_end, hy, func_y(grid[i,0,0]), y_start, y_end, step_bottom, step_top, N_step_bottom, N_step_top, eps_y)
		for j in range(Ny):
			grid[i,j,1] = line[j]
		
	save_to_vtk_one_array(path, grid, Nx, Ny)
	del(grid)
	del(line)
	del(step_bottom)
	del(step_top)
	return 1	

def smooth_curve_grid_3D(path, Nx, Ny, Nz, ind_critical, ind_regular_start, ind_regular_end, hx, hy, hz, func_z, z_start, z_end, eps_z):
	grid = np.zeros((Nx, Ny, Nz, 3), dtype = 'float')
	line = np.zeros((Nz), dtype = 'float')

#	N_step_bottom = ind_critical - ind_regular_start
#	N_step_top = N - ind_critical - ind_regular_end - 1
#	step_bottom = np.zeros((N_step_bottom), dtype = 'float')
#	step_top = np.zeros((N_step_top), dtype = 'float')

#here to save RAM, allocating np array one time per function
	N_step_bottom = ind_critical - ind_regular_start
	N_step_top = Nz - ind_critical - ind_regular_end - 1
	step_bottom = np.zeros((N_step_bottom), dtype = 'float')
	step_top = np.zeros((N_step_top), dtype = 'float')

	for i in range(Nx):
		grid[i,0,0,0] = i*hx
		for j in range(Ny):
			for k in range(Nz):
				grid[i,j,k,0] = grid[i,0,0,0]

	for j in range(Ny):
		grid[0,j,0,1] = j*hy
		for i in range(Nx):
			for k in range(Nz):
				grid[i,j,k,1] = grid[0,j,0,1]

	for i in range(Nx):
		for j in range(Ny):
			get_smooth_line(line, Nz, ind_critical, ind_regular_start, ind_regular_end, hz, func_z(grid[i,0,0,0], grid[0,j,0,1]), z_start, z_end, step_bottom, step_top, N_step_bottom, N_step_top, eps_z)
			for k in range(Nz):
				grid[i,j,k,2] = line[k]
	
	save_to_vtk_3D(path, grid, Nx, Ny, Nz)
	del(grid)
	del(line)
	del(step_bottom)
	del(step_top)
	return 1		

def find_s_max(L, I, h, s_min, I1):
	return (2.0*L + (h-s_min)*I1)/((float)(I)) - h

def find_ind_cr_gr(L, I, h, s_min):
	debug = True #TODO delete after testing
	if debug:
		file_name = "FILES/debug.txt"

	a = 2.0*(h-s_min)
	b = 4.0*(L-h*I)
	mc = (2.0*L-(h+s_min)*I)*I

	if False: #(a < 0.2):
		I1 = (int)(mc/b + 0.5)
		if debug:
			f = open(file_name, "at")
			f.write("L = %s	I = %s	h = %s	s_min = %s\n" %(L, I, h, s_min))
			I1_float = mc/b
			f.write("a = %s	b = %s	mc = %s	error = %s\n" %(a, b, mc, a*I1_float*I1_float + b*I1_float - mc))
			s_max = find_s_max(L, I, h, s_min, I1_float)
			f.write("(Linear) I1 = %s	s_max = %s	k1 = %s	k2 = %s\n" %(I1_float, s_max, (s_max - s_min)/I1_float, (s_max - h)/(I - I1_float)))
			f.close()
	else:
		I1 = (int)((sqrt(b*b + 4.0*a*mc) - b)/(2.0*a) + 0.5)
		if debug:
			f = open(file_name, "at")
			f.write("L = %s	I = %s	h = %s	s_min = %s\n" %(L, I, h, s_min))
			I1_float = (sqrt(b*b + 4.0*a*mc) - b)/(2.0*a)
			f.write("a = %s	b = %s	mc = %s	error = %s\n" %(a, b, mc, a*I1_float*I1_float + b*I1_float - mc))
			s_max = find_s_max(L, I, h, s_min, I1_float)
			f.write("I1 = %s	s_max = %s	k1 = %s	k2 = %s\n" %(I1_float, s_max, (s_max - s_min)/I1_float, (s_max - h)/(I - I1_float)))
			f.close()
	if (I1 < 0):
		return 0, find_s_max(L, I, h, s_min, 0)
	if (I1 > I):
		return I, find_s_max(L, I, h, s_min, I)
	return I1, find_s_max(L, I, h, s_min, I1)

def get_smooth_line(line, N, ind_critical, ind_regular_start, ind_regular_end, h, critical, start, end, step_bottom, step_top, N_step_bottom, N_step_top, eps):

	regular_pos = ind_critical * h + start
	if (abs(regular_pos - critical) < eps):
		line[0] = start
		for i in range(1,ind_critical):
			line[i] = line[i-1] + h
		line[N-1] = end
		for i in range(N-2,ind_critical,-1):
			line[i] = line[i+1] - h
		line[ind_critical] = critical
		return 2

	line[0] = start
	for i in range(ind_regular_start):
		line[i+1] = line[i] + h 
	i_start = ind_regular_start + 1

	line[N-1] = end
	for i in range(ind_regular_end):
		line[N-2-i] = line[N-1-i] - h
	i_end = N-ind_regular_end - 1

	regular_pos = ind_critical * h + start
	if (regular_pos > critical):
		top_greater = True
	else:
		top_greater = False

	def bad(dif):
		if (abs(dif) < eps):
			return False
		return True

	if top_greater: 
		if (ind_regular_start == 0):
			step_min = (critical-start)/((float)(N_step_bottom))
			for i in range(i_start,ind_critical):
				line[i] = i*step_min + start
		else:
			step_min = 2.0 * (critical-start-ind_regular_start*h) / ((float)(N_step_bottom)) - h
			if (step_min < eps):
				sys.exit("bad ind_critical code 3")
			delta_1 = (h - step_min) / ((float)(N_step_bottom))
			for i in range(N_step_bottom):
				step_bottom[i] = h - (i+1) * delta_1
			for i in range(N_step_bottom):
				line[i_start + i] = line[i_start + i - 1] + step_bottom[i]
			dif = critical - line[ind_critical]
			if bad(dif):
				delta_1 = dif / ((float)(N_step_bottom))
				for i in range(N_step_bottom-1):
					line[i_start + i] = line[i_start + i] + delta_1*(i+1)
#			line[ind_critical] = critical
			step_min = critical - line[ind_critical-1]

		if (ind_regular_end == 0):
			step_max = 2*(end-critical) / ((float)(N_step_top)) - step_min
			if (step_max < eps):
				sys.exit("bad ind_critical code 4")
			ind_cr_gr = N_step_top
		else:
			ind_cr_gr, step_max = find_ind_cr_gr(end-critical-ind_regular_end*h, N_step_top, h, step_min)

		if (ind_cr_gr > 0):
			delta_1 = (step_max - step_min) / ((float)(ind_cr_gr))
			for i in range(ind_cr_gr):
				step_top[i] = step_min + delta_1 * (i+1)
		int_temp = N_step_top - ind_cr_gr
		if (int_temp > 0):
			delta_1 = (step_max - h) / ((float)(int_temp))
			for i in range(int_temp):
				step_top[N_step_top - i - 1] = h + delta_1 * (i + 1)
			
		for i in range(N_step_top):
			line[i_end - 1 - i] = line[i_end - i] - step_top[N_step_top - 1 - i]
		dif = critical - line[ind_critical]
		if bad(dif):
			delta_1 = dif / ((float)(N_step_top))
			for i in range(N_step_top-1):
				line[i_end - 1 - i] = line[i_end - 1 - i] + delta_1*(i+1)
		line[ind_critical] = critical
			
	else:
		if (ind_regular_end == 0):
			step_min = (end-critical) / ((float)(N_step_top))
			for i in range(i_end-1,ind_critical,-1):
				line[i] = end - (i_end-i)*step_min
		else:
			step_min = 2.0 * (end-critical-ind_regular_end*h) / ((float)(N_step_top)) - h
			if (step_min < eps):
				sys.exit("bad ind_critical code 1")
			delta_1 = (h - step_min) / ((float)(N_step_top))
			for i in range(N_step_top):
				step_top[i] = h - (i+1) * delta_1
			for i in range(N_step_top):
				line[i_end - 1 - i] = line[i_end - i] - step_top[i]
			dif = critical - line[ind_critical]
			if bad(dif):
				delta_1 = dif / ((float)(N_step_top))
				for i in range(N_step_top-1):
					line[i_end - 1 - i] = line[i_end - 1 - i] + delta_1*(i+1)
#			line[ind_critical] = critical
			step_min = line[ind_critical + 1] - critical
		
		if (ind_regular_start == 0):
			step_max = 2*(critical-start) / ((float)(N_step_bottom)) - step_min
			if (step_max < eps):
				sys.exit("bad ind_critical code 2")
			ind_cr_gr = N_step_bottom
		else:
			ind_cr_gr, step_max = find_ind_cr_gr(critical-start-ind_regular_start*h, N_step_bottom, h, step_min)

		if (ind_cr_gr > 0):
			delta_1 = (step_max - step_min) / ((float)(ind_cr_gr))
			for i in range(ind_cr_gr):
				step_bottom[i] = step_min + delta_1 * (i+1)
		int_temp = N_step_bottom - ind_cr_gr
		if (int_temp > 0):
			delta_1 = (step_max - h) / ((float)(int_temp))
			for i in range(int_temp):
				step_bottom[N_step_bottom - i - 1] = h + delta_1 * (i + 1)

		for i in range(N_step_bottom):
			line[i_start + i] = line[i_start - 1 + i] + step_bottom[N_step_bottom - 1 - i]
		dif = critical - line[ind_critical]
		if bad(dif):
			delta_1 = dif / ((float)(N_step_bottom))
			for i in range(N_step_bottom-1):
				line[i_start + i] = line[i_start + i] + delta_1*(i+1)
		line[ind_critical] = critical

	return 1

def list_to_file_float(file_name, data, count):
	f = open(file_name, "wt")
	for i in range(count):
		f.write("%s	" %(data[i]))
	f.write("\n")
	f.close()
	return 1

def list_from_file_float(file_name, data, count):
	f = open(file_name, "rt")
	line = f.readline()
	temp = line.split("	")
	for i in range(count):
		data[i] = (float)(temp[i])
	f.close()
	return 1

def list_to_file_int(file_name, data, count):
	f = open(file_name, "wt")
	for i in range(count):
		f.write("%i	" %(data[i]))
	f.write("\n")
	f.close()
	return 1

def list_from_file_int(file_name, data, count):
	f = open(file_name, "rt")
	line = f.readline()
	temp = line.split("	")
	for i in range(count):
		data[i] = (int)(temp[i])
	f.close()
	return 1

def list_to_file(file_name, data, count):
	f = open(file_name, "wt")
	for i in range(count):
		f.write("%s	" %(data[i]))
	f.write("\n")
	f.close()
	return 1

def list_from_file(file_name, data, count):
	f = open(file_name, "rt")
	line = f.readline()
	temp = line.split("	")
	for i in range(count):
		data[i] = temp[i]
	f.close()
	return 1

def get_txt_file_name(folder_name, grid_from, grid_to):
	return folder_name + "/from_" + grid_from + "_to_" + grid_to + ".txt"

def print_interpolation_corrector(grid_from, grid_to, predictor_flag, corrector_flag, axis, folder_name):
	print_interpolation_corrector_with_given_file_name(grid_from, grid_to, predictor_flag, corrector_flag, axis, get_txt_file_name(folder_name, grid_from, grid_to))
	return 1
	
def print_interpolation_corrector_with_given_file_name_IFN(grid_from, grid_to, axis, side, folder_name, predictor_flag, corrector_flag, axis_corrector):
	print_interpolation_corrector_with_given_file_name(grid_from, grid_to, predictor_flag, corrector_flag, axis_corrector, IFN(folder_name, grid_from, grid_to, axis, side))
	return 1

def print_interpolation_corrector_with_given_file_name(grid_from, grid_to, predictor_flag, corrector_flag, axis, file_name):
	print("""    [contact]
        name = RectGridInterpolationCorrector""")
	print("        interpolation_file = " + file_name)
	print("        grid1 = " + grid_from)
	print("        grid2 = " + grid_to)
	if (predictor_flag == 1):
		print("        predictor_flag = true")
	else:
		print("        predictor_flag = false")
	if (corrector_flag == 1):
		print("        corrector_flag = true")
	else:
		print("        corrector_flag = false")
	print("        axis = %s" %(axis))
	print("    [/contact]")
	return 1

#???
#in rect to rect linear only
#grid_limits = np.zeros((3), dtype = 'int') 
#grid_limits[0] = const_ind_start, for example 0
#grid_limits[1] = const_ind_end, for example Nx-1, if axis = 0
#grid_limits[2] = const_ind_step, if 0, equals 1

#size, origin and step are given for axis 1, if axis 0. And vice versa. size is int.
#END ???


#grid_to_limits = np.zeros((6), dtype = 'int')
#grid_to_limits[0] = i_start, for example 0
#grid_to_limits[1] = i_end, for example Nx-1
#grid_to_limits[2] = j_start, for example 0
#grid_to_limits[3] = j_end, for example Ny-1
#grid_to_limits[4] = i_step, if 0, equals 1
#grid_to_limits[5] = j_step, if 0, equals 1

#grid_float_params = np.zeros((4), dtype = 'float')
#grid_float_params[0] = origin_x
#grid_float_params[1] = origin_y
#grid_float_params[2] = step_x
#grid_float_params[3] = step_y

#grid_from_size = np.zeros((2), dtype = 'int')
#grid_from_size[0] = size_x
#grid_from_size[1] = size_y

def gen_interpolation_file_rect_to_rect_2D(file_name, grid_from_id, grid_from_size, grid_from_float_params, grid_to_limits, grid_to_float_params, eps):
	f = open(file_name, "wt")
	f.write(grid_from_id + "\n")
	f.write("2D\n")

	ind = np.zeros((2), dtype = 'int')
	int_step_to = np.zeros((2), dtype = 'int')

	for d in range(2):
		if (grid_to_limits[4+d]==0):
			int_step_to[d] = 1
		else:
			int_step_to[d] = grid_to_limits[4+d]

	m_eps = 0.0 - eps
		
	ind1_X = []	
	ind1_Y = []
	ind2_X = []	
	ind2_Y = []
	
	IND_X = []
	IND_Y = []
		
	coef_X = []
	IND_coef_X = []
	IND_1_min_coef_X = []
	
	coef_Y = []
	IND_coef_Y = []
	IND_1_min_coef_Y = []	

	for d in range(2):
		for ind[d] in range(grid_to_limits[2*d], grid_to_limits[2*d+1]+1, int_step_to[d]):
			not_case_of_copy = True

			delta = grid_to_float_params[d] + grid_to_float_params[2+d]*ind[d] - grid_from_float_params[d]
			if (delta < m_eps):
				continue
			
			if (delta < eps):
				not_case_of_copy = False
				IND = 0
				
			if not_case_of_copy:
				delta_right = delta - (grid_from_size[d]-1)*grid_from_float_params[2+d]
				if (delta_right > eps):
					continue
				if (delta_right > m_eps):
					not_case_of_copy = False
					IND = grid_from_size[d]-1

			if not_case_of_copy:					
				IND = (int)(delta/grid_from_float_params[2+d])
				if (IND<0):
					continue
				if (IND>(grid_from_size[d]-2)):
					continue			
				
			if not_case_of_copy:		
				delta_left = delta - grid_from_float_params[2+d]*IND
				delta_right = grid_from_float_params[2+d] - delta_left
			
				if (abs(delta_left) < eps):
					not_case_of_copy = False
				elif (abs(delta_right) < eps):
					IND = IND + 1
					not_case_of_copy = False
			
				elif ((delta_left < 0.0) or (delta_right < 0.0)):
					sys.exit("bad rect to rect interpolation, %s, %s, %s, %s, %s, file " %(d, ind[d], delta_left, delta_right, grid_from_float_params[2+d]) + file_name)

			if (not not_case_of_copy):
				if (d == 0):
					IND_X.append(IND)
					ind1_X.append(ind[d])
				else:	
					IND_Y.append(IND)
					ind1_Y.append(ind[d])
				continue
			
			if (delta_left < delta_right):
				coef = delta_left/grid_from_float_params[2+d]
				IND_coef = IND + 1
				IND_1_min_coef = IND
			else:
				coef = delta_right/grid_from_float_params[2+d]
				IND_coef = IND
				IND_1_min_coef = IND + 1
				
			if (d == 0):
				coef_X.append(coef)
				IND_coef_X.append(IND_coef)
				IND_1_min_coef_X.append(IND_1_min_coef)
				ind2_X.append(ind[d])
			else:
				coef_Y.append(coef)
				IND_coef_Y.append(IND_coef)
				IND_1_min_coef_Y.append(IND_1_min_coef)
				ind2_Y.append(ind[d])

#	f.write("%i\n" %(((int)((grid_to_limits[1]-grid_to_limits[0])//i_step)+1) * ((int)((grid_to_limits[3]-grid_to_limits[2])//j_step)+1) ))

#	for i in range(grid_to_limits[0], grid_to_limits[1]+1, int_step_to[0]):
#		for j in range(grid_to_limits[2], grid_to_limits[3]+1, int_step_to[1]):

	f.write("%i\n" %( (len(ind1_X) + len(ind2_X) ) * (len(ind1_Y) + len(ind2_Y) ) ))

	#case of copy			
	for i in range(len(ind1_X)):
		for j in range(len(ind1_Y)):
			f.write("%s %s 1 %s %s 1.0\n" %(ind1_X[i], ind1_Y[j],   IND_X[i], IND_Y[j] ))

	#case of bilinear			
	for i in range(len(ind2_X)):
		for j in range(len(ind2_Y)):
			f.write("%s %s 4 %s %s %s %s %s %s %s %s %s %s %s %s\n" %(ind2_X[i], ind2_Y[j],   IND_1_min_coef_X[i], IND_1_min_coef_Y[j], (1.0-coef_X[i])*(1.0-coef_Y[j]),   IND_1_min_coef_X[i], IND_coef_Y[j], (1.0-coef_X[i])*coef_Y[j],   IND_coef_X[i], IND_1_min_coef_Y[j], coef_X[i]*(1.0-coef_Y[j]),   IND_coef_X[i], IND_coef_Y[j], coef_X[i]*coef_Y[j] ))

	#case of linear at dir X
	for i in range(len(ind1_X)):
		for j in range(len(ind2_Y)):
			f.write("%s %s 2 %s %s %s %s %s %s\n" %(ind1_X[i], ind2_Y[j],   IND_X[i], IND_1_min_coef_Y[j], (1.0-coef_Y[j]),   IND_X[i], IND_coef_Y[j], coef_Y[j] ))

	#case of linear at dir Y
	for i in range(len(ind2_X)):
		for j in range(len(ind1_Y)):
			f.write("%s %s 2 %s %s %s %s %s %s\n" %(ind2_X[i], ind1_Y[j],   IND_1_min_coef_X[i], IND_Y[j], (1.0-coef_X[i]),   IND_coef_X[i], IND_Y[j], coef_X[i] ))
	
	del(int_step_to)
	del(ind)

	del(ind1_X)	
	del(ind1_Y)
	del(ind2_X)	
	del(ind2_Y)

	del(IND_X)
	del(IND_Y)
		
	del(coef_X)
	del(IND_coef_X)
	del(IND_1_min_coef_X)
	
	del(coef_Y)
	del(IND_coef_Y)
	del(IND_1_min_coef_Y)	
			
	f.close()
	return 1
	
def gen_interpolation_file_rect_to_rect_2D_grid_from_limits(file_name, grid_from_id, grid_from_size, grid_from_float_params, grid_from_limits, grid_to_float_params, eps):

	grid_to_limits = np.zeros((6), dtype = 'int')
	eps_half = eps/2.0
#grid_to_limits[0] = i_start, for example 0
#grid_to_limits[1] = i_end, for example Nx-1
#grid_to_limits[2] = j_start, for example 0
#grid_to_limits[3] = j_end, for example Ny-1
#grid_to_limits[4] = i_step, if 0, equals 1
#grid_to_limits[5] = j_step, if 0, equals 1

#grid_float_params = np.zeros((4), dtype = 'float')
#grid_float_params[0] = origin_x
#grid_float_params[1] = origin_y
#grid_float_params[2] = step_x
#grid_float_params[3] = step_y

	x_start_shift = grid_from_float_params[0] + grid_from_limits[0] * grid_from_float_params[2] - grid_to_float_params[0]
	x_end_shift = grid_from_float_params[0] + grid_from_limits[1] * grid_from_float_params[2] - grid_to_float_params[0]
	y_start_shift = grid_from_float_params[1] + grid_from_limits[2] * grid_from_float_params[3] - grid_to_float_params[1]
	y_end_shift = grid_from_float_params[1] + grid_from_limits[3] * grid_from_float_params[3] - grid_to_float_params[1]

	grid_to_limits[0] = (int)(x_start_shift / grid_from_to_params[2] + 0.5)
	if ((x_start_shift - grid_to_limits[0] * grid_from_to_params[2]) > eps_half):
		grid_to_limits[0] = grid_to_limits[0] + 1

	grid_to_limits[1] = (int)(x_end_shift / grid_from_to_params[2] + 0.5)
	if ((grid_to_limits[1] * grid_from_to_params[2] - x_end_shift) > eps_half):
		grid_to_limits[1] = grid_to_limits[1] - 1

	grid_to_limits[2] = (int)(y_start_shift / grid_from_to_params[3] + 0.5)
	if ((y_start_shift - grid_to_limits[2] * grid_from_to_params[3]) > eps_half):
		grid_to_limits[2] = grid_to_limits[2] + 1

	grid_to_limits[3] = (int)(y_end_shift / grid_from_to_params[3] + 0.5)
	if ((grid_to_limits[3] * grid_from_to_params[3] - y_end_shift) > eps_half):
		grid_to_limits[3] = grid_to_limits[3] - 1

	gen_interpolation_file_rect_to_rect_2D(file_name, grid_from_id, grid_from_size, grid_from_float_params, grid_to_limits, grid_to_float_params, eps)
	
	del(grid_to_limits)
	return 1

#grid_to_limits = np.zeros((9), dtype = 'int')
#grid_to_limits[0] = i_start, for example 0
#grid_to_limits[1] = i_end, for example Nx-1
#grid_to_limits[2] = j_start, for example 0
#grid_to_limits[3] = j_end, for example Ny-1
#grid_to_limits[4] = k_start, for example 0
#grid_to_limits[5] = k_end, for example Nz-1
#grid_to_limits[6] = i_step, if 0, equals 1
#grid_to_limits[7] = j_step, if 0, equals 1
#grid_to_limits[8] = k_step, if 0, equals 1

#grid_to_z_limits = np.zeros((3), dtype = 'int') 
#grid_to_z_limits[0] = z_ind_start
#grid_to_z_limits[1] = z_ind_end
#grid_to_z_limits[2] = z_step, if 0, equals 1

#delay and rect_step_z are defined as follows: z_index_in_from_grid = z_index_in_to_grid * rect_step_z + delay

def gen_interpolation_file_rect_to_rect_2_3D(file_name, grid_from_id, grid_from_size, grid_from_float_params, grid_to_limits, grid_to_float_params, eps, delay, rect_step_z, grid_to_z_limits):
	f = open(file_name, "wt")
	f.write(grid_from_id + "\n")
	f.write("3D\n")

	ind = np.zeros((2), dtype = 'int')
	int_step_to = np.zeros((2), dtype = 'int')

	for d in range(2):
		if (grid_to_limits[4+d]==0):
			int_step_to[d] = 1
		else:
			int_step_to[d] = grid_to_limits[4+d]

	if (grid_to_z_limits[2]==0):
		z_step = 1
	else:
		z_step = grid_to_z_limits[2]

	m_eps = 0.0 - eps
		
	ind1_X = []	
	ind1_Y = []
	ind2_X = []	
	ind2_Y = []
	
	IND_X = []
	IND_Y = []
		
	coef_X = []
	IND_coef_X = []
	IND_1_min_coef_X = []
	
	coef_Y = []
	IND_coef_Y = []
	IND_1_min_coef_Y = []	

	for d in range(2):
		for ind[d] in range(grid_to_limits[2*d], grid_to_limits[2*d+1]+1, int_step_to[d]):
			not_case_of_copy = True

			delta = grid_to_float_params[d] + grid_to_float_params[2+d]*ind[d] - grid_from_float_params[d]
			if (delta < m_eps):
				continue
			
			if (delta < eps):
				not_case_of_copy = False
				IND = 0
				
			if not_case_of_copy:
				delta_right = delta - (grid_from_size[d]-1)*grid_from_float_params[2+d]
				if (delta_right > eps):
					continue
				if (delta_right > m_eps):
					not_case_of_copy = False
					IND = grid_from_size[d]-1

			if not_case_of_copy:					
				IND = (int)(delta/grid_from_float_params[2+d])
				if (IND<0):
					continue
				if (IND>(grid_from_size[d]-2)):
					continue			
				
			if not_case_of_copy:		
				delta_left = delta - grid_from_float_params[2+d]*IND
				delta_right = grid_from_float_params[2+d] - delta_left
			
				if (abs(delta_left) < eps):
					not_case_of_copy = False
				elif (abs(delta_right) < eps):
					IND = IND + 1
					not_case_of_copy = False
			
				elif ((delta_left < 0.0) or (delta_right < 0.0)):
					sys.exit("bad rect to rect interpolation, %s, %s, %s, %s, %s, file " %(d, ind[d], delta_left, delta_right, grid_from_float_params[2+d]) + file_name)

			if (not not_case_of_copy):
				if (d == 0):
					IND_X.append(IND)
					ind1_X.append(ind[d])
				else:	
					IND_Y.append(IND)
					ind1_Y.append(ind[d])
				continue
			
			if (delta_left < delta_right):
				coef = delta_left/grid_from_float_params[2+d]
				IND_coef = IND + 1
				IND_1_min_coef = IND
			else:
				coef = delta_right/grid_from_float_params[2+d]
				IND_coef = IND
				IND_1_min_coef = IND + 1
				
			if (d == 0):
				coef_X.append(coef)
				IND_coef_X.append(IND_coef)
				IND_1_min_coef_X.append(IND_1_min_coef)
				ind2_X.append(ind[d])
			else:
				coef_Y.append(coef)
				IND_coef_Y.append(IND_coef)
				IND_1_min_coef_Y.append(IND_1_min_coef)
				ind2_Y.append(ind[d])

	f.write("%i\n" %( (len(ind1_X) + len(ind2_X) ) * (len(ind1_Y) + len(ind2_Y) ) * ((int)((grid_to_z_limits[1]-grid_to_z_limits[0])//z_step)+1) ))

	def Func_K(k):
		return k*rect_step_z + delay

	for k in range(grid_to_z_limits[0], grid_to_z_limits[1]+1, z_step):
		K = Func_K(k)

		#case of copy			
		for i in range(len(ind1_X)):
			for j in range(len(ind1_Y)):
				f.write("%s %s %s 1 %s %s %s 1.0\n" %(ind1_X[i], ind1_Y[j], k,   IND_X[i], IND_Y[j], K ))

		#case of bilinear			
		for i in range(len(ind2_X)):
			for j in range(len(ind2_Y)):
				f.write("%s %s %s 4 %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" %(ind2_X[i], ind2_Y[j], k,   IND_1_min_coef_X[i], IND_1_min_coef_Y[j], K, (1.0-coef_X[i])*(1.0-coef_Y[j]),   IND_1_min_coef_X[i], IND_coef_Y[j], K, (1.0-coef_X[i])*coef_Y[j],   IND_coef_X[i], IND_1_min_coef_Y[j], K, coef_X[i]*(1.0-coef_Y[j]),   IND_coef_X[i], IND_coef_Y[j], K, coef_X[i]*coef_Y[j] ))

		#case of linear at dir X
		for i in range(len(ind1_X)):
			for j in range(len(ind2_Y)):
				f.write("%s %s %s 2 %s %s %s %s %s %s %s %s\n" %(ind1_X[i], ind2_Y[j], k,   IND_X[i], IND_1_min_coef_Y[j], K, (1.0-coef_Y[j]),  IND_X[i], IND_coef_Y[j], K, coef_Y[j] ))

		#case of linear at dir Y
		for i in range(len(ind2_X)):
			for j in range(len(ind1_Y)):
				f.write("%s %s %s 2 %s %s %s %s %s %s %s %s\n" %(ind2_X[i], ind1_Y[j], k,   IND_1_min_coef_X[i], IND_Y[j], K, (1.0-coef_X[i]),   IND_coef_X[i], IND_Y[j], K, coef_X[i] ))
	
	del(int_step_to)
	del(ind)

	del(ind1_X)	
	del(ind1_Y)
	del(ind2_X)	
	del(ind2_Y)

	del(IND_X)
	del(IND_Y)
		
	del(coef_X)
	del(IND_coef_X)
	del(IND_1_min_coef_X)
	
	del(coef_Y)
	del(IND_coef_Y)
	del(IND_1_min_coef_Y)	
			
	f.close()
	return 1


#grid_limits = np.zeros((3), dtype = 'int') 
#grid_limits[0] = ind_start
#grid_limits[1] = ind_end
#grid_limits[2] = i_step, if 0, equals 1

def gen_interpolation_file_rect_to_rect_linear_2D(file_name, grid_from_id, axis, grid_from_limits, grid_from_origin, grid_from_step, grid_from_size, grid_to_limits, grid_to_origin, grid_to_step, grid_to_size, eps):
	f = open(file_name, "wt")
	f.write(grid_from_id + "\n")
	f.write("2D\n")

	if (grid_to_limits[2]==0):
		const_ind_step = 1
	else:
		const_ind_step = grid_to_limits[2]
		
	if (grid_from_limits[2]==0):
		const_IND_step = 1
	else:
		const_IND_step = grid_from_limits[2]

	m_eps = 0.0 - eps
	
	# TODO similarly everywhere add (int)
	f.write("%i\n" %(((int)((grid_to_limits[1]-grid_to_limits[0])/const_ind_step)+1)*grid_to_size ))

	for var_ind in range(grid_to_size):
		not_case_of_copy = True

		delta = grid_to_origin + grid_to_step*var_ind - grid_from_origin
		if (delta < m_eps):
			continue
			
		if (delta < eps):
			not_case_of_copy = False
			var_IND = 0
				
		if not_case_of_copy:
			delta_right = delta - (grid_from_size-1)*grid_from_step
			if (delta_right > eps):
				continue
			if (delta_right > m_eps):
				not_case_of_copy = False
				var_IND = grid_from_size-1

		if not_case_of_copy:					
			var_IND = (int)(delta/grid_from_step)
			if (var_IND<0):
				continue
			if (var_IND>(grid_from_size-2)):
				continue			
				
		if not_case_of_copy:		
			delta_left = delta - grid_from_step*var_IND
			delta_right = grid_from_step - delta_left
			
			if (abs(delta_left) < eps):
				not_case_of_copy = False
			elif (abs(delta_right) < eps):
				var_IND = var_IND + 1
				not_case_of_copy = False
			
			elif ((delta_left < 0.0) or (delta_right < 0.0)):
				sys.exit("bad rect to rect interpolation, %s, %s, %s, %s, file " %(var_ind, delta_left, delta_right, grid_from_step) + file_name)	
				
		if (not not_case_of_copy):
			if (axis == 0):
				for const_ind in range(grid_to_limits[0], grid_to_limits[1]+1, const_ind_step):
					const_IND = grid_from_limits[0] + const_IND_step*((int)((const_ind-grid_to_limits[0])/const_ind_step))
					f.write("%s %s 1 %s %s 1.0\n" %(const_ind, var_ind, const_IND, var_IND ))
				continue
			for const_ind in range(grid_to_limits[0], grid_to_limits[1]+1, const_ind_step):
				const_IND = grid_from_limits[0] + const_IND_step*((int)((const_ind-grid_to_limits[0])/const_ind_step))
				f.write("%s %s 1 %s %s 1.0\n" %(var_ind, const_ind, var_IND, const_IND )) #axis == 1
			continue
			
		if (delta_left < delta_right):
			coef = delta_left/grid_from_step
			var_IND_coef = var_IND + 1
			var_IND_1_min_coef = var_IND
		else:
			coef = delta_right/grid_from_step
			var_IND_coef = var_IND
			var_IND_1_min_coef = var_IND + 1
			
		if (axis == 0):	
			for const_ind in range(grid_to_limits[0], grid_to_limits[1]+1, const_ind_step):
				const_IND = grid_from_limits[0] + const_IND_step*((int)((const_ind-grid_to_limits[0])/const_ind_step))
				f.write("%s %s 2 %s %s %s %s %s %s\n" %(const_ind, var_ind, const_IND, var_IND_1_min_coef, (1.0 - coef), const_IND, var_IND_coef, coef )) 
			continue
		for const_ind in range(grid_to_limits[0], grid_to_limits[1]+1, const_ind_step):
			const_IND = grid_from_limits[0] + const_IND_step*((int)((const_ind-grid_to_limits[0])/const_ind_step))
			f.write("%s %s 2 %s %s %s %s %s %s\n" %(var_ind, const_ind, var_IND_1_min_coef, const_IND, (1.0 - coef), var_IND_coef, const_IND, coef )) # axis == 1			
			
	f.close()
	return 1

#Attention! axis 0 or 1, not 2

def gen_interpolation_file_rect_to_rect_linear_2_3D(file_name, grid_from_id, axis, grid_from_limits, grid_from_origin, grid_from_step, grid_from_size, grid_to_limits, grid_to_origin, grid_to_step, grid_to_size, eps, delay, rect_step_z, grid_to_z_limits):
	f = open(file_name, "wt")
	f.write(grid_from_id + "\n")
	f.write("2D\n")

	if (grid_to_limits[2]==0):
		const_ind_step = 1
	else:
		const_ind_step = grid_to_limits[2]
		
	if (grid_from_limits[2]==0):
		const_IND_step = 1
	else:
		const_IND_step = grid_from_limits[2]

	if (grid_to_z_limits[2]==0):
		z_step = 1
	else:
		z_step = grid_to_z_limits[2]

	m_eps = 0.0 - eps
	
	# TODO similarly everywhere add (int)
	f.write("%i\n" %(((int)((grid_to_limits[1]-grid_to_limits[0])/const_ind_step)+1)*grid_to_size*((int)((grid_to_z_limits[1]-grid_to_z_limits[0])//z_step)+1) ))

	def Func_K(k):
		return k*rect_step_z + delay

	for var_ind in range(grid_to_size):
		not_case_of_copy = True
			
		delta = grid_to_origin + grid_to_step*var_ind - grid_from_origin
		if (delta < m_eps):
			continue
			
		if (delta < eps):
			not_case_of_copy = False
			var_IND = 0
				
		if not_case_of_copy:
			delta_right = delta - (grid_from_size-1)*grid_from_step
			if (delta_right > eps):
				continue
			if (delta_right > m_eps):
				not_case_of_copy = False
				var_IND = grid_from_size-1

		if not_case_of_copy:					
			var_IND = (int)(delta/grid_from_step)
			if (var_IND<0):
				continue
			if (var_IND>(grid_from_size-2)):
				continue			
				
		if not_case_of_copy:		
			delta_left = delta - grid_from_step*var_IND
			delta_right = grid_from_step - delta_left
			
			if (abs(delta_left) < eps):
				not_case_of_copy = False
			elif (abs(delta_right) < eps):
				var_IND = var_IND + 1
				not_case_of_copy = False	
			
			elif ((delta_left < 0.0) or (delta_right < 0.0)):
				sys.exit("bad rect to rect interpolation, %s, %s, %s, %s, file " %(var_ind, delta_left, delta_right, grid_from_step) + file_name)
				
		if (not not_case_of_copy):
			if (axis == 0):
				for const_ind in range(grid_to_limits[0], grid_to_limits[1]+1, const_ind_step):
					const_IND = grid_from_limits[0] + const_IND_step*((int)((const_ind-grid_to_limits[0])/const_ind_step))
					for k in range(grid_to_z_limits[0], grid_to_z_limits[1]+1, z_step):
						f.write("%s %s %s 1 %s %s %s 1.0\n" %(const_ind, var_ind, k, const_IND, var_IND, Func_K(k) ))
				continue
			for const_ind in range(grid_to_limits[0], grid_to_limits[1]+1, const_ind_step):
				const_IND = grid_from_limits[0] + const_IND_step*((int)((const_ind-grid_to_limits[0])/const_ind_step))
				for k in range(grid_to_z_limits[0], grid_to_z_limits[1]+1, z_step):
					f.write("%s %s %s 1 %s %s %s 1.0\n" %(var_ind, const_ind, k, var_IND, const_IND, Func_K(k) )) #axis == 1
			continue
			
		if (delta_left < delta_right):
			coef = delta_left/grid_from_step
			var_IND_coef = var_IND + 1
			var_IND_1_min_coef = var_IND
		else:
			coef = delta_right/grid_from_step
			var_IND_coef = var_IND
			var_IND_1_min_coef = var_IND + 1
			
		if (axis == 0):	
			for const_ind in range(grid_to_limits[0], grid_to_limits[1]+1, const_ind_step):
				const_IND = grid_from_limits[0] + const_IND_step*((int)((const_ind-grid_to_limits[0])/const_ind_step))
				for k in range(grid_to_z_limits[0], grid_to_z_limits[1]+1, z_step):
					K = Func_K(k)
					f.write("%s %s %s 2 %s %s %s %s %s %s %s %s\n" %(const_ind, var_ind, k, const_IND, var_IND_1_min_coef, K, (1.0 - coef), const_IND, var_IND_coef, K, coef )) 
			continue
		for const_ind in range(grid_to_limits[0], grid_to_limits[1]+1, const_ind_step):
			const_IND = grid_from_limits[0] + const_IND_step*((int)((const_ind-grid_to_limits[0])/const_ind_step))
			for k in range(grid_to_z_limits[0], grid_to_z_limits[1]+1, z_step):
				K = Func_K(k)
				f.write("%s %s %s 2 %s %s %s %s %s %s %s %s\n" %(var_ind, const_ind, k, var_IND_1_min_coef, const_IND, K, (1.0 - coef), var_IND_coef, const_IND, K, coef )) # axis == 1			
			
	f.close()
	return 1



#2_3D and 3D COPYING INTERPOLATION TXT FILES GENERATION

#grid_to_limits = np.zeros((9), dtype = 'int')
#grid_to_limits[0] = i_start, for example 0
#grid_to_limits[1] = i_end, for example Nx-1
#grid_to_limits[2] = j_start, for example 0
#grid_to_limits[3] = j_end, for example Ny-1
#grid_to_limits[4] = k_start, for example 0
#grid_to_limits[5] = k_end, for example Nz-1
#grid_to_limits[6] = i_step, if 0, equals 1
#grid_to_limits[7] = j_step, if 0, equals 1
#grid_to_limits[8] = k_step, if 0, equals 1

#Func_I(i,j,k) returns I in grid_from, if i, j, and k in grid_to are known
#Func_J(i,j,k) returns J in grid_from, if i, j, and k in grid_to are known
#Func_K(i,j,k) returns K in grid_from, if i, j, and k in grid_to are known #TODO fixed J >> K in all places

def gen_interpolation_file_copy_rect_3D(file_name, grid_from_id, grid_to_limits, Func_I, Func_J, Func_K):
	f = open(file_name, "wt")
	f.write(grid_from_id + "\n")
	f.write("3D\n")
	
	if (grid_to_limits[6]==0):
		i_step = 1
	else:
		i_step = grid_to_limits[6]
		
	if (grid_to_limits[7]==0):
		j_step = 1
	else:
		j_step = grid_to_limits[7]

	if (grid_to_limits[8]==0):
		k_step = 1
	else:
		k_step = grid_to_limits[8]

	f.write("%i\n" %(((grid_to_limits[1]-grid_to_limits[0])//i_step+1)*((grid_to_limits[3]-grid_to_limits[2])//j_step+1)*((grid_to_limits[5]-grid_to_limits[4])//k_step+1)))
		
	for k in range(grid_to_limits[4], grid_to_limits[5]+1, k_step):
		for i in range(grid_to_limits[0], grid_to_limits[1]+1, i_step):
			for j in range(grid_to_limits[2], grid_to_limits[3]+1, j_step):
				f.write("%s %s %s 1 %s %s %s 1.0\n" %(i, j, k, Func_I(i,j,k), Func_J(i,j,k), Func_K(i,j,k) ))
	f.close()
	return 1

#grid_from_limits 1, 3, 5 does not matter
def gen_interpolation_file_copy_rect_3D_both_limits(file_name, grid_from_id, grid_from_limits, grid_to_limits):
	if (grid_from_limits[6] == 0):
		step_I = 1
	else:
		step_I = grid_from_limits[6]

	if (grid_to_limits[6] == 0):
		step_i = 1
	else:
		step_i = grid_to_limits[6]

	if (grid_from_limits[7] == 0):
		step_J = 1
	else:
		step_J = grid_from_limits[7]

	if (grid_to_limits[7] == 0):
		step_j = 1
	else:
		step_j = grid_to_limits[7]

	if (grid_from_limits[8] == 0):
		step_K = 1
	else:
		step_K = grid_from_limits[8]

	if (grid_to_limits[8] == 0):
		step_k = 1
	else:
		step_k = grid_to_limits[8]

	def Func_I(i,j,k):			
		return grid_from_limits[0] + step_I*((int)((i-grid_to_limits[0])/step_i))

	def Func_J(i,j,k):
		return grid_from_limits[2] + step_J*((int)((j-grid_to_limits[2])/step_j))
		
	def Func_K(i,j,k):
		return grid_from_limits[4] + step_K*((int)((k-grid_to_limits[4])/step_k))

	gen_interpolation_file_copy_rect_3D(file_name, grid_from_id, grid_to_limits, Func_I, Func_J, Func_K)
	return 1

#attention! Func_AXIS should return right index.
def gen_interpolation_file_copy_rect_3D_slices_odd_cell_number(file_name, grid_from_id, grid_to_limits, Func_I, Func_J, Func_K, axis):
	f = open(file_name, "wt")
	f.write(grid_from_id + "\n")
	f.write("3D\n")
	
	if (grid_to_limits[6]==0):
		i_step = 1
	else:
		i_step = grid_to_limits[6]
		
	if (grid_to_limits[7]==0):
		j_step = 1
	else:
		j_step = grid_to_limits[7]

	if (grid_to_limits[8]==0):
		k_step = 1
	else:
		k_step = grid_to_limits[8]

	f.write("%i\n" %(((grid_to_limits[1]-grid_to_limits[0])//i_step+1)*((grid_to_limits[3]-grid_to_limits[2])//j_step+1)*((grid_to_limits[5]-grid_to_limits[4])//k_step+1)))
		
	for k in range(grid_to_limits[4], grid_to_limits[5]+1, k_step):	
		for i in range(grid_to_limits[0], grid_to_limits[1]+1, i_step):
			for j in range(grid_to_limits[2], grid_to_limits[3]+1, j_step):
				I = Func_I(i,j,k)
				J = Func_J(i,j,k)
				K = Func_K(i,j,k)
				
				if (axis == 0):
					f.write("%s %s %s 2 %s %s %s 0.5 %s %s %s 0.5\n" %(i, j, k, I - 1, J, K, I, J, K ))
				elif (axis == 1):
					f.write("%s %s %s 2 %s %s %s 0.5 %s %s %s 0.5\n" %(i, j, k, I, J - 1, K, I, J, K ))
				else: #if (axis == 2):
					f.write("%s %s %s 2 %s %s %s 0.5 %s %s %s 0.5\n" %(i, j, k, I, J, K - 1, I, J, K ))

	f.close()
	return 1

#curved

def gen_interpolation_file_from_rect_to_curved_1_3D(file_name, grid_from_id, curved_limits_z, eps, points_z, size, origin_z, spacing_z):
	f = open(file_name, "wt")
	f.write(grid_from_id + "\n")
	f.write("3D\n")
	
	index = []
	i_one = []
	j_one = []
	k_curved_one = []
	i_two = []
	j_two = []
	k_curved_two = []
	index_1 = []
	index_2 = []
	coef_1 = []
	coef_2 = []

	for i in range(size[0]):
		for j in range(size[1]):
			for k_curved in range(curved_limits_z[0],curved_limits_z[1]+1):
				dz = points_z[i,j,k_curved] - origin_z
				k_rect_left = (int)(dz/spacing_z + 0.0001)
				k_rect_right = (int)(dz/spacing_z + 1.0001)
				if ((k_rect_right - k_rect_left) != 1):
					sys.exit("gen_interpolation_file_from_rect_to_curved_1_3D")
				if (k_rect_left<0):
					continue
				if (k_rect_left>(size[2]-2)):
					continue
					
				d_left = (dz - spacing_z * k_rect_left)
				if (abs(d_left) < eps):
					index.append(k_rect_left)
					i_one.append(i)
					j_one.append(j)
					k_curved_one.append(k_curved)
					continue
					
				d_right = (spacing_z * k_rect_right - dz)
				if (abs(d_right) < eps):
					index.append(k_rect_right)
					i_one.append(i)
					j_one.append(j)
					k_curved_one.append(k_curved)
					continue
					
				if (d_left < (spacing_z / 2.0)):
					t_coef_2 = d_left / spacing_z
					t_coef_1 = 1.0 - t_coef_2
				else:
					t_coef_1 = d_right / spacing_z
					t_coef_2 = 1.0 - t_coef_1
				coef_1.append(t_coef_1)
				coef_2.append(t_coef_2)
				i_two.append(i)
				j_two.append(j)
				index_1.append(k_rect_left)
				index_2.append(k_rect_right)
				k_curved_two.append(k_curved)					

	f.write("%i\n" %(len(i_one)+len(i_two)))					
					
	for ind in range(len(i_one)):
		f.write("%i %i %i 1 %i %i %i 1.0\n" %(i_one[ind], j_one[ind], k_curved_one[ind],    i_one[ind], j_one[ind], index[ind] ))

	for ind in range(len(i_two)):
		f.write("%i %i %i 2 %i %i %i %s %i %i %i %s\n" %(i_two[ind], j_two[ind], k_curved_two[ind],    i_two[ind], j_two[ind], index_1[ind], coef_1[ind],    i_two[ind], j_two[ind], index_2[ind], coef_2[ind] ))

	del(index)
	del(i_one)
	del(j_one)
	del(k_curved_one)
	del(i_two)
	del(j_two)
	del(k_curved_two)
	del(index_1)
	del(index_2)
	del(coef_1)
	del(coef_2)
	
	return 1

#будет работать только если координата по z в криволинейной сетке возрастает.
def gen_interpolation_file_from_curved_to_rect_1_3D(file_name, grid_from_id, curved_limits_z, eps, points_z, size, origin_z, spacing_z):
	interpolation_radius = get_interpolation_radius()

	f = open(file_name, "wt")
	f.write(grid_from_id + "\n")
	f.write("3D\n")

	index = []
	i_one = []
	j_one = []
	k_rect_one = []
	i_two = []
	j_two = []
	k_rect_two = []
	index_1 = []
	index_2 = []
	coef_1 = []
	coef_2 = []
	checked = np.zeros((size[2]+4), dtype = 'bool')

	for i in range(size[0]):
		for j in range(size[1]):
			for k_rect in range(size[2]+4):
				checked[k_rect] = 0
			
			for k_curved in range(curved_limits_z[0],curved_limits_z[1]):
				k_rect_start = max(-2, (int)((points_z[i,j,k_curved] - origin_z)/spacing_z + 0.5) - interpolation_radius)
				k_rect_end = min(size[2]+2, (int)((points_z[i,j,k_curved+1] - origin_z)/spacing_z + 0.5) + interpolation_radius + 1)
				
				for k_rect in range(k_rect_start, k_rect_end):
					if checked[k_rect+2]:
						continue

					rect_coord = k_rect * spacing_z + origin_z
					
					d_left = rect_coord - points_z[i,j,k_curved]
					if (abs(d_left) < eps):
						index.append(k_curved)
						i_one.append(i)
						j_one.append(j)
						k_rect_one.append(k_rect)
						checked[k_rect+2] = True
						continue

					if (d_left < 0.0):
						continue
						
					d_right = points_z[i,j,k_curved+1] - rect_coord
					if (abs(d_right) < eps):
						index.append(k_curved+1)
						i_one.append(i)
						j_one.append(j)
						k_rect_one.append(k_rect)
						checked[k_rect+2] = True
						continue
						
					if (d_right < 0.0):
						continue
						
					spacing_curved = points_z[i,j,k_curved+1] - points_z[i,j,k_curved]

					if (d_left < (spacing_curved / 2.0)):
						t_coef_2 = d_left / spacing_curved
						t_coef_1 = 1.0 - t_coef_2
					else:
						t_coef_1 = d_right / spacing_curved
						t_coef_2 = 1.0 - t_coef_1
					coef_1.append(t_coef_1)
					coef_2.append(t_coef_2)
					i_two.append(i)
					j_two.append(j)
					index_1.append(k_curved)
					index_2.append(k_curved + 1)
					k_rect_two.append(k_rect)
					checked[k_rect+2] = True

	f.write("%i\n" %(len(i_one)+len(i_two)))					
					
	for ind in range(len(i_one)):
		f.write("%i %i %i 1 %i %i %i 1.0\n" %(i_one[ind], j_one[ind], k_rect_one[ind],    i_one[ind], j_one[ind], index[ind] ))

	for ind in range(len(i_two)):
		f.write("%i %i %i 2 %i %i %i %s %i %i %i %s\n" %(i_two[ind], j_two[ind], k_rect_two[ind],    i_two[ind], j_two[ind], index_1[ind], coef_1[ind],    i_two[ind], j_two[ind], index_2[ind], coef_2[ind] ))
	
	del(checked)
	del(index)
	del(i_one)
	del(j_one)
	del(k_rect_one)
	del(i_two)
	del(j_two)
	del(k_rect_two)
	del(index_1)
	del(index_2)
	del(coef_1)
	del(coef_2)
	
	return 1


def gen_interpolation_file_from_rect_to_curved_1_2D(file_name, grid_from_id, curved_limits_y, eps, points_y, size, origin_y, spacing_y):
	f = open(file_name, "wt")
	f.write(grid_from_id + "\n")
	f.write("2D\n")
	
	index = []
	i_one = []
	j_curved_one = []
	i_two = []
	j_curved_two = []
	index_1 = []
	index_2 = []
	coef_1 = []
	coef_2 = []

	for i in range(size[0]):
		for j_curved in range(curved_limits_y[0],curved_limits_y[1]+1):
			dy = points_y[i,j_curved] - origin_y
			j_rect_left = (int)(dy/spacing_y + 0.0001)
			j_rect_right = (int)(dy/spacing_y + 1.0001)
			if ((j_rect_right - j_rect_left) != 1):
				sys.exit("gen_interpolation_file_from_rect_to_curved_1_2D")
			if (j_rect_left<0):
				continue
			if (j_rect_left>(size[1]-2)):
				continue
					
			d_left = (dy - spacing_y * j_rect_left)
			if (abs(d_left) < eps):
				index.append(j_rect_left)
				i_one.append(i)
				j_curved_one.append(j_curved)
				continue
					
			d_right = (spacing_y * j_rect_right - dy)
			if (abs(d_right) < eps):
				index.append(j_rect_right)
				i_one.append(i)
				j_curved_one.append(j_curved)
				continue
					
			if (d_left < (spacing_y / 2.0)):
				t_coef_2 = d_left / spacing_y
				t_coef_1 = 1.0 - t_coef_2
			else:
				t_coef_1 = d_right / spacing_y
				t_coef_2 = 1.0 - t_coef_1
			coef_1.append(t_coef_1)
			coef_2.append(t_coef_2)
			i_two.append(i)
			index_1.append(j_rect_left)
			index_2.append(j_rect_right)
			j_curved_two.append(j_curved)					

	f.write("%i\n" %(len(i_one)+len(i_two)))					
					
	for ind in range(len(i_one)):
		f.write("%i %i 1 %i %i 1.0\n" %(i_one[ind], j_curved_one[ind],    i_one[ind], index[ind] ))

	for ind in range(len(i_two)):
		f.write("%i %i 2 %i %i %s %i %i %s\n" %(i_two[ind], j_curved_two[ind],    i_two[ind], index_1[ind], coef_1[ind],    i_two[ind], index_2[ind], coef_2[ind] ))

	del(index)
	del(i_one)
	del(j_curved_one)
	del(i_two)
	del(j_curved_two)
	del(index_1)
	del(index_2)
	del(coef_1)
	del(coef_2)
	
	return 1

#будет работать только если координата по z в криволинейной сетке возрастает.
def gen_interpolation_file_from_curved_to_rect_1_2D(file_name, grid_from_id, curved_limits_y, eps, points_y, size, origin_y, spacing_y):
	interpolation_radius = get_interpolation_radius()

	f = open(file_name, "wt")
	f.write(grid_from_id + "\n")
	f.write("2D\n")

	index = []
	i_one = []
	j_rect_one = []
	i_two = []
	j_rect_two = []
	index_1 = []
	index_2 = []
	coef_1 = []
	coef_2 = []
	checked = np.zeros((size[1]+4), dtype = 'bool')

	for i in range(size[0]):
		for j_rect in range(size[1]+4):
			checked[j_rect] = 0
			
		for j_curved in range(curved_limits_y[0],curved_limits_y[1]):
			j_rect_start = max(-2, (int)((points_y[i,j_curved] - origin_y)/spacing_y + 0.5) - interpolation_radius)
			j_rect_end = min(size[1]+2, (int)((points_y[i,j_curved+1] - origin_y)/spacing_y + 0.5) + interpolation_radius + 1)
				
			for j_rect in range(j_rect_start, j_rect_end):
				if checked[j_rect+2]:
					continue

				rect_coord = j_rect * spacing_y + origin_y
					
				d_left = rect_coord - points_y[i,j_curved]
				if (abs(d_left) < eps):
					index.append(j_curved)
					i_one.append(i)
					j_rect_one.append(j_rect)
					checked[j_rect+2] = True
					continue

				if (d_left < 0.0):
					continue
						
				d_right = points_y[i,j_curved+1] - rect_coord
				if (abs(d_right) < eps):
					index.append(j_curved+1)
					i_one.append(i)
					j_rect_one.append(j_rect)
					checked[j_rect+2] = True
					continue
						
				if (d_right < 0.0):
					continue
						
				spacing_curved = points_y[i,j_curved+1] - points_y[i,j_curved]

				if (d_left < (spacing_curved / 2.0)):
					t_coef_2 = d_left / spacing_curved
					t_coef_1 = 1.0 - t_coef_2
				else:
					t_coef_1 = d_right / spacing_curved
					t_coef_2 = 1.0 - t_coef_1
				coef_1.append(t_coef_1)
				coef_2.append(t_coef_2)
				i_two.append(i)
				index_1.append(j_curved)
				index_2.append(j_curved + 1)
				j_rect_two.append(j_rect)
				checked[j_rect+2] = True

	f.write("%i\n" %(len(i_one)+len(i_two)))					
					
	for ind in range(len(i_one)):
		f.write("%i %i 1 %i %i 1.0\n" %(i_one[ind], j_rect_one[ind],    i_one[ind], index[ind] ))

	for ind in range(len(i_two)):
		f.write("%i %i 2 %i %i %s %i %i %s\n" %(i_two[ind], j_rect_two[ind],    i_two[ind], index_1[ind], coef_1[ind],    i_two[ind], index_2[ind], coef_2[ind] ))
	
	del(checked)
	del(index)
	del(i_one)
	del(j_rect_one)
	del(i_two)
	del(j_rect_two)
	del(index_1)
	del(index_2)
	del(coef_1)
	del(coef_2)
	
	return 1

#rect_contacts are deleted here in 2_3D case
#delay and rect_step_z are defined as follows: z_index_in_rect_grid = z_index_in_curved_grid * rect_step_z + delay
def gen_interpolation_file_from_rect_to_curved_2_3D(file_name, grid_from_id, curved_limits, eps, points_x, points_y, size, origin, spacing, delay, rect_step_z):
	f = open(file_name, "wt")
	f.write(grid_from_id + "\n")
	f.write("3D\n")
	
	eps_sq = eps*eps
	eps_99 = 0.0 - 0.99*eps
	
	if (curved_limits[6]==0):
		step_x = 1
	else:
		step_x = curved_limits[6]
		
	if (curved_limits[7]==0):
		step_y = 1
	else:
		step_y = curved_limits[7]
		
	if (curved_limits[8]==0):
		step_z = 1
	else:
		step_z = curved_limits[8]
		
	def Func_K(k):
		return k*rect_step_z + delay


#	list_to_file("curved_limits", curved_limits, 9)

	array = np.zeros((5,2), dtype = 'float')
	dis_sq = np.zeros((4), dtype = 'float')

	array1 = []
	array2_int = []
	array2_float = []
	array4_int = []
	array4_float = []
		
	for i in range(curved_limits[0], curved_limits[1]+1, step_x):
		for j in range(curved_limits[2], curved_limits[3]+1, step_y):
			dx = points_x[i,j] - origin[0]
			dy = points_y[i,j] - origin[1]

#			if ((dx < 0.0)or(dy < 0.0)):
#				continue

			#TODO without 0.99 and without the following double operations in these extreme cases
			if ((dx < eps_99)or(dy < eps_99)):
				continue
				
			I_flag = True
			J_flag = True
			
			if (dx < 0.0):
				I = 0
				I_flag = False
				
			if (dy < 0.0):
				J = 0
				J_flag = False
			
			if I_flag:
				dx_right = (size[0]-1)*spacing[0] - dx
				if (dx_right < eps_99):
					continue
				if (dx_right < eps):
					I = size[0] - 2
					I_flag = False
					
			if J_flag:
				dy_right = (size[1]-1)*spacing[1] - dy
				if (dy_right < eps_99):
					continue
				if (dy_right < eps):
					J = size[1] - 2
					J_flag = False
			
			if I_flag:
				I = (int)(dx/spacing[0])
				if (I<0):
					continue
				if (I>(size[0]-2)):
					continue
					
			if J_flag:
				J = (int)(dy/spacing[1])
				if (J<0):
					continue
				if (J>(size[1]-2)):
					continue		
		
			array[0,0] = points_x[i,j]
			array[0,1] = points_y[i,j]
			for num in range(1,5):
				get_rect_coords(num, ind_I(I,num), ind_J(J,num), array, origin, spacing)
			
			continue_flag = 0
			for num in range(1,5):
				dis_sq[num-1] = get_dis_sqrt(0,num,array)
				if (dis_sq[num-1] < eps_sq):
					array1.append(i)
					array1.append(j)
					array1.append(ind_I(I,num))
					array1.append(ind_J(J,num))
					continue_flag = 1
					break
			if (continue_flag == 1):
				continue
			
			continue_flag = 0
			for num in range(1,5):
				next_num = get_next_num(num)
				if (dis_side_rect(array, num, eps)):
					if ((num == 1) or (num == 3)):
						axis = 0
					else:
						axis = 1
					if (dis_sq[num-1]<dis_sq[next_num-1]):
						alpha = (array[0,axis]-array[num,axis])/spacing[axis]
						alpha = fabs(alpha)
						if (alpha>1.0):
							sys.exit("Bad interpolation coef " + file_name)
						array2_int.append(i)
						array2_int.append(j)
						array2_int.append(ind_I(I,num))
						array2_int.append(ind_J(J,num))
						array2_float.append(1.0-alpha)
						array2_int.append(ind_I(I,next_num))
						array2_int.append(ind_J(J,next_num))
						array2_float.append(alpha)
					else:
						alpha = (array[0,axis]-array[next_num,axis])/spacing[axis]
						alpha = fabs(alpha)
						if (alpha>1.0):
							sys.exit("Bad interpolation coef " + file_name)
						array2_int.append(i)
						array2_int.append(j)
						array2_int.append(ind_I(I,num))
						array2_int.append(ind_J(J,num))
						array2_float.append(alpha)
						array2_int.append(ind_I(I,next_num))
						array2_int.append(ind_J(J,next_num))
						array2_float.append(1.0-alpha)							
					continue_flag = 1
					break					
			if (continue_flag == 1):
				continue
				
			#if not inside_quadr(array):
			#	sys.exit("Not inside quadr " + file_name)
			
			if (dis_sq[0]<dis_sq[1]): # 1 and 2
				alpha_plus = (array[0,0] - array[1,0])/spacing[0]
				alpha = 1.0 - alpha_plus
			else:
				alpha = (array[2,0] - array[0,0])/spacing[0]
				alpha_plus = 1.0 - alpha
			if (bad_coef(alpha)):
				sys.exit("Bad interpolation coef " + file_name)
				
			if (dis_sq[0]<dis_sq[3]): # 1 and 4
				beta_plus = (array[0,1] - array[1,1])/spacing[1]
				beta = 1.0 - beta_plus
			else:
				beta = (array[4,1] - array[0,1])/spacing[1]
				beta_plus = 1.0 - beta
			if (bad_coef(beta)):
				sys.exit("Bad interpolation coef " + file_name)			
			
			array4_int.append(i)
			array4_int.append(j)
			array4_int.append(I)
			array4_int.append(J)
			array4_float.append(alpha*beta)
			array4_int.append(I+1)
			array4_int.append(J)
			array4_float.append(alpha_plus*beta)			
			array4_int.append(I)
			array4_int.append(J+1)
			array4_float.append(alpha*beta_plus)			
			array4_int.append(I+1)
			array4_int.append(J+1)
			array4_float.append(alpha_plus*beta_plus)

	count1 = int(len(array1)/4)
	count2 = int(len(array2_float)/2)
	count4 = int(len(array4_float)/4)
	
	f.write("%i\n" %(( count1+count2+count4 )*( (curved_limits[5]-curved_limits[4])//step_z+1 )) ) 

	for k in range(curved_limits[4], curved_limits[5]+1, step_z):
		K = Func_K(k)	
		for i in range(count1):
			f.write("%s %s %s 1 %s %s %s 1.0\n" %(array1[4*i], array1[4*i+1], k,   array1[4*i+2], array1[4*i+3], K ))
		for i in range(count2):
			f.write("%s %s %s 2 %s %s %s %s %s %s %s %s\n" %(array2_int[6*i], array2_int[6*i+1], k,   array2_int[6*i+2], array2_int[6*i+3], K, array2_float[2*i],   array2_int[6*i+4], array2_int[6*i+5], K, array2_float[2*i+1] ))
		for i in range(count4):
			f.write("%s %s %s 4 %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" %(array4_int[10*i], array4_int[10*i+1], k,    array4_int[10*i+2], array4_int[10*i+3], K, array4_float[4*i],    array4_int[10*i+4], array4_int[10*i+5], K, array4_float[4*i+1],    array4_int[10*i+6], array4_int[10*i+7], K, array4_float[4*i+2],    array4_int[10*i+8], array4_int[10*i+9], K, array4_float[4*i+3]))

	f.close()

	del(array1)
	del(array2_int)
	del(array2_float)
	del(array4_int)
	del(array4_float)

	#TODO the same arrays insert on lines in curved to rect 2D and 2_3D
	
	return 1

#curved_limits 6,7 does not matter
#eps_N - solve quadratic equation or use Newton's method. Value 0.3 is recommended, not more
#delay and rect_step_z are defined as follows: z_index_in_rect_grid = z_index_in_curved_grid * rect_step_z + delay
#if_test_axis_x, if_test_axis_y are 0 or 1
def gen_interpolation_file_from_curved_to_rect_2_3D(file_name, grid_from_id, curved_limits, eps, eps_N, points_x, points_y, size, origin, spacing, if_checked_to_draw, checked_path, delay, rect_step_z, if_test_axis_x, if_test_axis_y):

	interpolation_radius = get_interpolation_radius()

	f = open(file_name, "wt")
	f.write(grid_from_id + "\n")
	f.write("3D\n")

	if (curved_limits[8]==0):
		step_z = 1
	else:
		step_z = curved_limits[8]

	def Func_K(k):
		return k*rect_step_z + delay
	
	eps_sq = eps*eps
		
	array = np.zeros((5,2), dtype = 'float')
	dis_sq = np.zeros((4), dtype = 'float')
	height = np.zeros((4), dtype = 'float')
	length = np.zeros((4), dtype = 'float')
	checked = np.zeros((size[0]+4,size[1]+4), dtype = 'bool')

	line1_0 = [] #to write in file
	line1_1 = [] #to write in file

	line2_0 = [] #to write in file
	line2_1 = [] #to write in file
	line2_2 = [] #to write in file
	line2_3 = [] #to write in file

	line4_0 = [] #to write in file
	line4_1 = [] #to write in file
	line4_2 = [] #to write in file
	line4_3 = [] #to write in file
	line4_4 = [] #to write in file
	line4_5 = [] #to write in file

	res = np.zeros((2), dtype = 'float')
	info = np.zeros((5), dtype = 'int')
	XY = np.zeros((9), dtype = 'float')
	vectors = np.zeros((4,2), dtype = 'float')
	vectors0 = np.zeros((4,2), dtype = 'float')
	flags_vectors = np.ones((4), dtype = 'bool')
				
	#filling checked, koefs
	count1 = 0
	count2 = 0
	count4 = 0
	for I in range(curved_limits[0], curved_limits[1]):
		for J in range(curved_limits[2], curved_limits[3]):
			info[0] = 0 #new quadrangle
			for i in range(4):
				flags_vectors[i] = True
		
			for num in range(1,5):
				array[num,0] = points_x[ind_I(I,num),ind_J(J,num)]
				array[num,1] = points_y[ind_I(I,num),ind_J(J,num)]
		
			Xmin = min(array[1,0],array[2,0],array[3,0],array[4,0])
			Xmax = max(array[1,0],array[2,0],array[3,0],array[4,0])		
			Ymin = min(array[1,1],array[2,1],array[3,1],array[4,1])
			Ymax = max(array[1,1],array[2,1],array[3,1],array[4,1])	

			i_start = (int)((Xmin-origin[0])/spacing[0]+0.5) - interpolation_radius
			i_end = (int)((Xmax-origin[0])/spacing[0]+0.5) + interpolation_radius + 1
			j_start = (int)((Ymin-origin[1])/spacing[1]+0.5) - interpolation_radius
			j_end = (int)((Ymax-origin[1])/spacing[1]+0.5) + interpolation_radius + 1
			
			i_start = max(-2,i_start)
			i_end = min(size[0]+2,i_end)
			j_start = max(-2,j_start)
			j_end = min(size[1]+2,j_end)
			
			for i in range(i_start,i_end):
				for j in range(j_start,j_end):
					if (checked[i+2,j+2]):
						continue
				
					get_rect_coords(0, i, j, array, origin, spacing)
		
					continue_flag = 0
					for num in range(1,5):
						dis_sq[num-1] = get_dis_sqrt_vectors(0,num,array,vectors0[num-1])
						if (dis_sq[num-1] < eps_sq):
							line1_0.append("%s %s" %(i, j)) #"%s %s "" 1 %s %s " " 1.0\n"
							line1_1.append("%s %s" %(ind_I(I,num), ind_J(J,num) ))
							count1 = count1 + 1
							checked[i+2,j+2] = True
							continue_flag = 1
							break
					if (continue_flag == 1):
						continue
			
					continue_flag = 0
					for num in range(1,5):
						next_num = get_next_num(num)
						if (flags_vectors[num-1]):
							length[num-1] = sqrt(get_dis_sqrt_vectors(next_num, num, array, vectors[num-1]))
							flags_vectors[num-1] = False
						height[num-1] = fabs(area(array, num, next_num, 0)) / length[num-1]
						if (height[num-1] < eps):
							if (scalar_product(vectors0[num-1],vectors[num-1]) < 0.0):
								continue
							alpha = Linear_Interpolation_on_Line(dis_sq[num-1], height[num-1], length[num-1])
							if (alpha>1.0):
								continue_flag = 1
								break
							if (dis_sq[num-1]<dis_sq[next_num-1]):
								line2_0.append("%s %s"%(i, j))
								line2_1.append("%s %s" %(ind_I(I,num), ind_J(J,num) ))
								line2_2.append("%s %s %s" %(1.0-alpha, ind_I(I,next_num), ind_J(J,next_num) ))
								line2_3.append("%s" %(alpha ))
								count2 = count2 + 1
							else:
								alpha = Linear_Interpolation_on_Line(dis_sq[next_num-1], height[num-1], length[num-1])
								if (alpha>1.0):
									sys.exit("2 Bad interpolation coef " + file_name)
								line2_0.append("%s %s" %(i, j))
								line2_1.append("%s %s" %(ind_I(I,num), ind_J(J,num) ))
								line2_2.append("%s %s %s" %(alpha, ind_I(I,next_num), ind_J(J,next_num) ))
								line2_3.append("%s" %(1.0-alpha ))
								count2 = count2 + 1
							if ((sqrt(dis_sq[num-1])>length[num-1]) or (sqrt(dis_sq[next_num-1])>length[num-1])): #TODO investigate if it is needed to delete
								sys.exit("2 Bad interpolation coef " + file_name) #TODO investigate if it is needed to delete
							checked[i+2,j+2] = True
							continue_flag = 1
							break					
					if (continue_flag == 1):
						continue
				
					if not_inside_quadr_given_vectors(vectors, vectors0):
						continue
						
					if (Bilinear_Interpolation(res, array, vectors, vectors0, height, length, eps_sq, eps_N, info, XY) == -1):					
						sys.exit("Bad interpolation coef found using bilinear interpolation" + file_name)
					if (bad_coef(res[0]) or bad_coef(res[1])):
						sys.exit("Bad interpolation coef " + file_name)

					alpha_plus = res[0]
					alpha = 1.0 - alpha_plus
					beta_plus = res[1]
					beta = 1.0 - beta_plus

					line4_0.append("%s %s" %(i, j ))
					line4_1.append("%s %s " %(ind_I(I,info[1]), ind_J(J,info[1]) ))
					line4_2.append("%s %s %s" %(alpha*beta, ind_I(I,info[2]), ind_J(J,info[2]) ))
					line4_3.append("%s %s %s" %(alpha_plus*beta, ind_I(I,info[4]), ind_J(J,info[4]) ))
					line4_4.append("%s %s %s" %(alpha*beta_plus, ind_I(I,info[3]), ind_J(J,info[3]) ))
					line4_5.append("%s" %(alpha_plus*beta_plus))
					count4 = count4 + 1
					checked[i+2,j+2] = True

#TODO
#start testing
#	good_x = 1
#	good_y = 1
#
#	if (if_test_axis_x == 1):
#		break_flag = False
#		for j in range(size[1]+4):
#			total = 0
#			for i in range(size[0]+4):
#				if checked[i,j]:
#					total = total + 1
#				elif (total > 0):
#					if (total == 1):
#						good_x = 0
#						break_flag = True
#						break
#					else:
#						total = 0
			
#			if break_flag:
#				break
#			if (total == 1):
#				good_x = 0
#				break		

#	if (if_test_axis_y == 1):
#		break_flag = False
#		for i in range(size[0]+4):
#			total = 0
#			for j in range(size[1]+4):
#				if checked[i,j]:
#					total = total + 1
#				elif (total > 0):
#					if (total == 1):
#						good_y = 0
#						break_flag = True
#						break
#					else:
#						total = 0
#			
#			if break_flag:
#				break
#			if (total == 1):
#				good_y = 0
#				break	
	
#end testing

	f.write("%i\n" %((count1 + count2 + count4)*((curved_limits[5]-curved_limits[4])//step_z+1)))

	for k in range(curved_limits[4], curved_limits[5]+1, step_z):
		K = Func_K(k)
				
		for i in range(count1):
			f.write(line1_0[i] + " %s 1 " %(K) + line1_1[i] + " %s 1.0\n" %(k) )

		for i in range(count2):
			f.write(line2_0[i] + " %s 2 " %(K) + line2_1[i] + " %s " %(k) + line2_2[i] + " %s " %(k) + line2_3[i] + "\n" )			

		for i in range(count4):
			f.write(line4_0[i] + " %s 4 " %(K) + line4_1[i] + " %s " %(k) + line4_2[i] + " %s " %(k)+ line4_3[i] + " %s " %(k) + line4_4[i] + " %s " %(k) + line4_5[i] + "\n" )
			
	f.close()
	
	if (if_checked_to_draw == 1):
		checked_float = np.zeros((size[0]+4,size[1]+4), dtype = 'float')
		for i in range(size[0]+4):
			for j in range(size[1]+4):
				if (checked[i,j]):
					checked_float[i,j] = 1.0
		data_to_vtk(checked_path,checked_float,origin[0]-2.0*spacing[0],origin[1]-2.0*spacing[1],size[0]+4,size[1]+4,spacing[0],spacing[1])

	del(checked) #TODO add del in all old functions

	del(line1_0)
	del(line1_1)

	del(line2_0)
	del(line2_1)
	del(line2_2)
	del(line2_3)

	del(line4_0)
	del(line4_1)
	del(line4_2)
	del(line4_3)
	del(line4_4)
	del(line4_5)
	
	return 1 #good_x * good_y

#ARC LEN VARYING

#set_of_coefs = np.zeros((count, 5), dtype = 'float')
#set_of_coefs[i,0] = center_start
#set_of_coefs[i,1] = delta_height, should be in range (-1.0, 1.0)
#set_of_coefs[i,2] = delta_left, positive or zero
#set_of_coefs[i,3] = delta_right, positive or zero
#set_of_coefs[i,4] = delta_center, positive or zero

#temp = np.zeros((max_N), dtype = 'float') 
#max_N should be enough

#limits_to[0] = start, example 0
#limits_to[1] = end, example N-1
#limits_to[2] = step, if 0, equals 1
def gen_interpolation_file_2_3D_from_file(file_name_result, file_name_from, axis, limits_to, Func_IND):
	file_from = open(file_name_from, "rt")
	file_result = open(file_name_result, "wt")
	
	if (limits_to[2]==0):
		step = 1
	else:
		step = limits_to[2]
		
	if (axis == 0):
		line_ind = 0
		for line in file_from:
			if (line_ind == 0):
				file_result.write(line)
				line_ind = 1
				continue
			elif (line_ind == 1):
				file_result.write("3D\n")
				line_ind = 2
				continue
			elif (line_ind == 2):
				file_result.write("%i\n" %((int)(line)*((int)((limits_to[1] - limits_to[0])//step + 1)) ))
				line_ind = 3
				continue
			temp = line.split(" ")
			count = (int)(temp[2])

			if (count == 1):
				for ind in range(limits_to[0], limits_to[1] + 1, step):
					file_result.write("%i %s %s 1 %i %s %s 1.0\n" %(ind, temp[0], temp[1], Func_IND(ind), temp[3], temp[4]))

			elif (count == 2):
				for ind in range(limits_to[0], limits_to[1] + 1, step):
					IND = Func_IND(ind)
					file_result.write("%i %s %s 2 %i %s %s %s %i %s %s %s" %(ind, temp[0], temp[1], IND, temp[3], temp[4], temp[5], IND, temp[6], temp[7], temp[8]))

			else: #if (count == 4):
				for ind in range(limits_to[0], limits_to[1] + 1, step):
					IND = Func_IND(ind)
					file_result.write("%i %s %s 4 %i %s %s %s %i %s %s %s %i %s %s %s %i %s %s %s" %(ind, temp[0], temp[1], IND, temp[3], temp[4], temp[5], IND, temp[6], temp[7], temp[8], IND, temp[9], temp[10], temp[11], IND, temp[12], temp[13], temp[14]))

	elif (axis == 1):
		line_ind = 0
		for line in file_from:
			if (line_ind == 0):
				file_result.write(line)
				line_ind = 1
				continue
			elif (line_ind == 1):
				file_result.write("3D\n")
				line_ind = 2
				continue
			elif (line_ind == 2):
				file_result.write("%i\n" %((int)(line)*((int)((limits_to[1] - limits_to[0])//step + 1)) ))
				line_ind = 3
				continue
			temp = line.split(" ")
			count = (int)(temp[2])

			if (count == 1):
				for ind in range(limits_to[0], limits_to[1] + 1, step):
					file_result.write("%s %i %s 1 %s %i %s 1.0\n" %(temp[0], ind, temp[1], temp[3], Func_IND(ind), temp[4]))

			elif (count == 2):
				for ind in range(limits_to[0], limits_to[1] + 1, step):
					IND = Func_IND(ind)
					file_result.write("%s %i %s 2 %s %i %s %s %s %i %s %s" %(temp[0], ind, temp[1], temp[3], IND, temp[4], temp[5], temp[6], IND, temp[7], temp[8]))

			else: #if (count == 4):
				for ind in range(limits_to[0], limits_to[1] + 1, step):
					IND = Func_IND(ind)
					file_result.write("%s %i %s 4 %s %i %s %s %s %i %s %s %s %i %s %s %s %i %s %s" %(temp[0], ind, temp[1], temp[3], IND, temp[4], temp[5], temp[6], IND, temp[7], temp[8], temp[9], IND, temp[10], temp[11], temp[12], IND, temp[13], temp[14]))

	else: #if (axis == 2):
		line_ind = 0
		for line in file_from:
			if (line_ind == 0):
				file_result.write(line)
				line_ind = 1
				continue
			elif (line_ind == 1):
				file_result.write("3D\n")
				line_ind = 2
				continue
			elif (line_ind == 2):
				file_result.write("%i\n" %((int)(line)*((int)((limits_to[1] - limits_to[0])//step + 1)) ))
				line_ind = 3
				continue
			temp = line.split(" ")
			count = (int)(temp[2])

			if (count == 1):
				for ind in range(limits_to[0], limits_to[1] + 1, step):
					file_result.write("%s %s %i 1 %s %s %i 1.0\n" %(temp[0], temp[1], ind, temp[3], temp[4], Func_IND(ind)))

			elif (count == 2):
				for ind in range(limits_to[0], limits_to[1] + 1, step):
					IND = Func_IND(ind)
					file_result.write("%s %s %i 2 %s %s %i %s %s %s %i %s" %(temp[0], temp[1], ind, temp[3], temp[4], IND, temp[5], temp[6], temp[7], IND, temp[8]))

			else: #if (count == 4):
				for ind in range(limits_to[0], limits_to[1] + 1, step):
					IND = Func_IND(ind)
					file_result.write("%s %s %i 4 %s %s %i %s %s %s %i %s %s %s %i %s %s %s %i %s" %(temp[0], temp[1], ind, temp[3], temp[4], IND, temp[5], temp[6], temp[7], IND, temp[8], temp[9], temp[10], IND, temp[11], temp[12], temp[13], IND, temp[14]))
	
	file_from.close()
	file_result.close()
	return 1

def create_set_of_varying_steps_diff_number_of_points_linear(step, arc_len_start, arc_len_end, set_of_coefs, count, last_fixed_coord):
	return create_set_of_varying_steps_diff_number_of_points_func(step, arc_len_start, arc_len_end, set_of_coefs, count, varying_step_from_arc_len_linear, last_fixed_coord)

def create_set_of_varying_steps_diff_number_of_points_parabolic(step, arc_len_start, arc_len_end, set_of_coefs, count, last_fixed_coord):
	return create_set_of_varying_steps_diff_number_of_points_func(step, arc_len_start, arc_len_end, set_of_coefs, count, varying_step_from_arc_len, last_fixed_coord)

def create_set_of_varying_steps_diff_number_of_points_func(step, arc_len_start, arc_len_end, set_of_coefs, count, func, last_fixed_coord):
	sign_arc = sign(arc_len_end - arc_len_start) 
	steps = []
	points = []

	step_finish = step * func(arc_len_end, set_of_coefs, count)

	points.append(arc_len_start)
	while (((arc_len_end - points[len(points)-1])*sign_arc) > step_finish):
		ind = len(points)-1
		steps.append(step * func(points[ind], set_of_coefs, count))
		points.append(points[ind] + sign_arc*steps[ind])
	
	diff = arc_len_end - sign_arc*step_finish - points[len(points)-1]

	if (diff == 0.0):
		points.append(arc_len_end)
#		for i in range(len(points)):
#			temp[i] = points[i]
#		N = len(points)
#		del(points)
		del(steps)
		return points

	sign_diff = sign_arc*sign(diff) #### sign_arc*sign(diff)
	diff = abs(diff)

	if (((arc_len_start - last_fixed_coord)*sign_arc) >= 0.0):
		step_diff = sign_diff*diff/((float)(len(steps)))
		for i in range(len(steps)):
			steps[i] = steps[i] + step_diff		
			if (steps[i] < 0.0):
				sys.exit("Bad set_of_coef for function create_set_of_varying_steps_diff_number_of_points()")

	else:	
		numbers = []
		for i in range(len(points)):
			if (((points[i] - last_fixed_coord)*sign_arc) >= 0.0):
				if (i < len(steps)):
					numbers.append(i)

		if (len(numbers) > 0):	
			step_diff = sign_diff*diff/((float)(len(numbers)))	
			for i in range(len(numbers)):
				steps[numbers[i]] = steps[numbers[i]] + step_diff
				if (steps[numbers[i]] < 0.0):
					sys.exit("Bad set_of_coef for function create_set_of_varying_steps_diff_number_of_points()")
		else:
					sys.exit("Bad last_fixed_coord or set_of_coef for function create_set_of_varying_steps_diff_number_of_points()")
		del(numbers)
		
	for i in range(len(points)-2):
		points[i+1] = points[i] + sign_arc*steps[i]
	
	points[len(points)-1] = arc_len_end - sign_arc*step_finish
	points.append(arc_len_end)

#	for i in range(len(points)):
#		temp[i] = points[i]

#	N = len(points)
#	del(points)
	del(steps)
	return points

#points = np.zeros((number), dtype = 'float')
#total integral of product of hat funcions should be equal to (arc_len_end - arc_len_start)
def create_set_of_varying_steps(points, number, arc_len_start, arc_len_end, set_of_coefs, count, if_evenly):
	sign_arc = sign(arc_len_end - arc_len_start) 
	step = sign_arc*(arc_len_end - arc_len_start)/((float)(number-1))
	steps = np.zeros((number-1), dtype = 'float')

	points[0] = arc_len_start
	for i in range(number-1):
		steps[i] = step * varying_step_from_arc_len(points[i], set_of_coefs, count)
		points[i+1] = points[i] + sign_arc*steps[i]

	diff = arc_len_end - points[number-1]

	if (diff == 0.0):
		return 0.0

	sign_diff = sign_arc*sign(diff) #### sign_arc*sign(diff)
	diff = fabs(diff)

	if (if_evenly == 1):
		step_diff = sign_diff*diff/((float)(number-1))
		for i in range(number-1):
			steps[i] = steps[i] + step_diff		
			if (steps[i] < 0.0):
				sys.exit("Bad set_of_coef for function create_set_of_varying_steps()")

	else:	
		numbers = []
		for i in range(number-1):
			if ((sign_diff == -1) and (steps[i] < step)):
				numbers.append(i)
			elif ((sign_diff == 1) and (steps[i] > step)):
				numbers.append(i)

		if (len(numbers) > 0):	
			step_diff = sign_diff*diff/((float)(len(numbers)))	
			for i in range(len(numbers)):
				steps[numbers[i]] = steps[numbers[i]] + step_diff
				if (steps[numbers[i]] < 0.0):
					sys.exit("Bad set_of_coef for function create_set_of_varying_steps()")
		
	for i in range(number-1):
		points[i+1] = points[i] + sign_arc*steps[i]
	
	points[number-1] = arc_len_end
	return diff*sign_diff

#integral is equal to 0.5
def fi_hat_half_normalized(f_x): #will be in range [0.0, 1.0], non_zero in range = [0.0, 1.0], f_x SHOULD BE POZITIVE or equal to zero
	if (f_x > 1.0): 
		return 0.0
	elif (f_x < 0.5): 
		return 1.0 - 2.0*f_x*f_x 	
	temp = f_x - 1.0
	return 2.0*temp*temp 

def varying_step_from_arc_len(arc_len, set_of_coefs, count):
	res = 1.0
	for i in range(count):
		res = res * used_hat_func(arc_len, i, set_of_coefs)
	return res

def varying_step_from_arc_len_linear(arc_len, set_of_coefs, count):
	res = 1.0
	for i in range(count):
		res = res + used_hat_func_linear(arc_len, i, set_of_coefs)
	return res

#integral is equal to (delta_left + delta_right)*(1.0 + 0.5*delta_height) = 0.5*(delta_left + delta_right)*delta_height + 1.0*(delta_left + delta_right)
def used_hat_func(arc_len, i, set_of_coefs):
	f_x = arc_len - set_of_coefs[i,0] 
	if (f_x < 0.0):
		if (set_of_coefs[i,2] == 0.0):
			return 1.0
		f_x = fabs(f_x)/set_of_coefs[i,2] #left
	else: 
		f_x = f_x - set_of_coefs[i,4]
		if (f_x <= 0.0):
			return 1.0 + set_of_coefs[i,1]
		if (set_of_coefs[i,3] == 0.0):
			return 1.0
		f_x = f_x/set_of_coefs[i,3] #right
	return 1.0 + set_of_coefs[i,1] * fi_hat_half_normalized(f_x)

def used_hat_func_linear(arc_len, i, set_of_coefs):
	f_x = arc_len - set_of_coefs[i,0] 
	if (f_x < 0.0):
		if (set_of_coefs[i,2] == 0.0):
			return 0.0
		f_x = fabs(f_x)/set_of_coefs[i,2] #left
	else: 
		f_x = f_x - set_of_coefs[i,4]
		if (f_x <= 0.0):
			return set_of_coefs[i,1]
		if (set_of_coefs[i,3] == 0.0):
			return 0.0
		f_x = f_x/set_of_coefs[i,3] #right
	if (f_x >= 1.0):
		return 0.0
	return set_of_coefs[i,1] * (1.0 - f_x)

#TXT files generation for interpolation contacts

#grid_to_limits = np.zeros((6), dtype = 'int')
#grid_to_limits[0] = i_start, for example 0
#grid_to_limits[1] = i_end, for example Nx-1
#grid_to_limits[2] = j_start, for example 0
#grid_to_limits[3] = j_end, for example Ny-1
#grid_to_limits[4] = i_step, if 0, equals 1
#grid_to_limits[5] = j_step, if 0, equals 1

#Func_I(i,j) returns I in grid_from, if i and j in grid_to are known
#Func_J(i,j) returns J in grid_from, if i and j in grid_to are known

def gen_interpolation_file_copy_rect_2D(file_name, grid_from_id, grid_to_limits, Func_I, Func_J):
	f = open(file_name, "wt")
	f.write(grid_from_id + "\n")
	f.write("2D\n")
	
	if (grid_to_limits[4]==0):
		i_step = 1
	else:
		i_step = grid_to_limits[4]
		
	if (grid_to_limits[5]==0):
		j_step = 1
	else:
		j_step = grid_to_limits[5]

	f.write("%i\n" %(((grid_to_limits[1]-grid_to_limits[0])//i_step+1)*((grid_to_limits[3]-grid_to_limits[2])//j_step+1)))
		
	for i in range(grid_to_limits[0], grid_to_limits[1]+1, i_step):
		for j in range(grid_to_limits[2], grid_to_limits[3]+1, j_step):
			f.write("%s %s 1 %s %s 1.0\n" %(i, j, Func_I(i,j), Func_J(i,j) ))
	f.close()
	return 1
	
def gen_interpolation_file_copy_rect_2D_from_limits(file_name, grid_from_id, grid_from_limits, Func_i, Func_j):
	f = open(file_name, "wt")
	f.write(grid_from_id + "\n")
	f.write("2D\n")
	
	if (grid_from_limits[4]==0):
		I_step = 1
	else:
		I_step = grid_from_limits[4]
		
	if (grid_from_limits[5]==0):
		J_step = 1
	else:
		J_step = grid_from_limits[5]

	f.write("%i\n" %(((grid_from_limits[1]-grid_from_limits[0])//I_step+1)*((grid_from_limits[3]-grid_from_limits[2])//J_step+1)))
		
	for I in range(grid_from_limits[0], grid_from_limits[1]+1, I_step):
		for J in range(grid_from_limits[2], grid_from_limits[3]+1, J_step):
			f.write("%s %s 1 %s %s 1.0\n" %(Func_i(I,J), Func_j(I,J), I, J ))
	f.close()
	return 1

#grid_from_limits 1, 3 does not matter
def gen_interpolation_file_copy_rect_2D_both_limits(file_name, grid_from_id, grid_from_limits, grid_to_limits):
	if (grid_from_limits[4] == 0):
		step_I = 1
	else:
		step_I = grid_from_limits[4]

	if (grid_to_limits[4] == 0):
		step_i = 1
	else:
		step_i = grid_to_limits[4]

	if (grid_from_limits[5] == 0):
		step_J = 1
	else:
		step_J = grid_from_limits[5]

	if (grid_to_limits[5] == 0):
		step_j = 1
	else:
		step_j = grid_to_limits[5]

	def Func_I(i,j):			
		return grid_from_limits[0] + step_I*((int)((i-grid_to_limits[0])/step_i))

	def Func_J(i,j):
		return grid_from_limits[2] + step_J*((int)((j-grid_to_limits[2])/step_j))
		
	gen_interpolation_file_copy_rect_2D(file_name, grid_from_id, grid_to_limits, Func_I, Func_J)
	return 1

#point is the point of rotation
#angle is the angle of rotation of grid_to
def gen_interpolation_file_from_rect_to_rotated_rect_2D(file_name, grid_from_id, grid_to_limits, eps, size_from, origin_from, spacing_from, angle, point, origin_to, spacing_to):
	f = open(file_name, "wt")
	f.write(grid_from_id + "\n")
	f.write("2D\n")
	
	eps_sq = eps*eps
	eps_99 = 0.0 - 0.99*eps
	
	if (grid_to_limits[4]==0):
		step_x = 1
	else:
		step_x = grid_to_limits[4]
		
	if (grid_to_limits[5]==0):
		step_y = 1
	else:
		step_y = grid_to_limits[5]

	array = np.zeros((5,2), dtype = 'float')
	dis_sq = np.zeros((4), dtype = 'float')

	array1 = []
	array2_int = []
	array2_float = []
	array4_int = []
	array4_float = []
		
	for i in range(grid_to_limits[0], grid_to_limits[1]+1, step_x):
		for j in range(grid_to_limits[2], grid_to_limits[3]+1, step_y):
			rx, ry = rotated_2D(origin_to[0] + i*spacing_to[0], origin_to[1] + j*spacing_to[1], angle, point[0], point[1]) 
			
			dx = rx - origin_from[0]
			dy = ry - origin_from[1]

			#TODO without 0.99 and without the following double operations in these extreme cases
			if ((dx < eps_99)or(dy < eps_99)):
				continue
				
			I_flag = True
			J_flag = True
			
			if (dx < 0.0):
				I = 0
				I_flag = False
				
			if (dy < 0.0):
				J = 0
				J_flag = False
			
			if I_flag:
				dx_right = (size_from[0]-1)*spacing_from[0] - dx
				if (dx_right < eps_99):
					continue
				if (dx_right < eps):
					I = size_from[0] - 2
					I_flag = False
					
			if J_flag:
				dy_right = (size_from[1]-1)*spacing_from[1] - dy
				if (dy_right < eps_99):
					continue
				if (dy_right < eps):
					J = size_from[1] - 2
					J_flag = False
			
			if I_flag:
				I = (int)(dx/spacing_from[0])
				if (I<0):
					continue
				if (I>(size_from[0]-2)):
					continue
					
			if J_flag:
				J = (int)(dy/spacing_from[1])
				if (J<0):
					continue
				if (J>(size_from[1]-2)):
					continue						
		
			array[0,0] = rx
			array[0,1] = ry
			for num in range(1,5):
				get_rect_coords(num, ind_I(I,num), ind_J(J,num), array, origin_from, spacing_from)
			
			continue_flag = 0
			for num in range(1,5):
				dis_sq[num-1] = get_dis_sqrt(0,num,array)
				if (dis_sq[num-1] < eps_sq):
					array1.append(i)
					array1.append(j)
					array1.append(ind_I(I,num))
					array1.append(ind_J(J,num))
					continue_flag = 1
					break
			if (continue_flag == 1):
				continue
			
			continue_flag = 0
			for num in range(1,5):
				next_num = get_next_num(num)
				if (dis_side_rect(array, num, eps)):
					if ((num == 1) or (num == 3)):
						axis = 0
					else:
						axis = 1
					if (dis_sq[num-1]<dis_sq[next_num-1]):
						alpha = (array[0,axis]-array[num,axis])/spacing_from[axis]
						alpha = fabs(alpha)
						if (alpha>1.0):
							sys.exit("Bad interpolation coef " + file_name)
						array2_int.append(i)
						array2_int.append(j)
						array2_int.append(ind_I(I,num))
						array2_int.append(ind_J(J,num))
						array2_float.append(1.0-alpha)
						array2_int.append(ind_I(I,next_num))
						array2_int.append(ind_J(J,next_num))
						array2_float.append(alpha)
					else:
						alpha = (array[0,axis]-array[next_num,axis])/spacing_from[axis]
						alpha = fabs(alpha)
						if (alpha>1.0):
							sys.exit("Bad interpolation coef " + file_name)
						array2_int.append(i)
						array2_int.append(j)
						array2_int.append(ind_I(I,num))
						array2_int.append(ind_J(J,num))
						array2_float.append(alpha)
						array2_int.append(ind_I(I,next_num))
						array2_int.append(ind_J(J,next_num))
						array2_float.append(1.0-alpha)
					continue_flag = 1
					break					
			if (continue_flag == 1):
				continue
				
			#if not inside_quadr(array):
			#	sys.exit("Not inside quadr " + file_name)
			
			if (dis_sq[0]<dis_sq[1]): # 1 and 2
				alpha_plus = (array[0,0] - array[1,0])/spacing_from[0]
				alpha = 1.0 - alpha_plus
			else:
				alpha = (array[2,0] - array[0,0])/spacing_from[0]
				alpha_plus = 1.0 - alpha
			if (bad_coef(alpha)):
				sys.exit("Bad interpolation coef " + file_name)
				
			if (dis_sq[0]<dis_sq[3]): # 1 and 4
				beta_plus = (array[0,1] - array[1,1])/spacing_from[1]
				beta = 1.0 - beta_plus
			else:
				beta = (array[4,1] - array[0,1])/spacing_from[1]
				beta_plus = 1.0 - beta
			if (bad_coef(beta)):
				sys.exit("Bad interpolation coef " + file_name)			
			array4_int.append(i)
			array4_int.append(j)
			array4_int.append(I)
			array4_int.append(J)
			array4_float.append(alpha*beta)
			array4_int.append(I+1)
			array4_int.append(J)
			array4_float.append(alpha_plus*beta)			
			array4_int.append(I)
			array4_int.append(J+1)
			array4_float.append(alpha*beta_plus)			
			array4_int.append(I+1)
			array4_int.append(J+1)
			array4_float.append(alpha_plus*beta_plus)			

	count1 = int(len(array1)/4)
	count2 = int(len(array2_float)/2)
	count4 = int(len(array4_float)/4)
	
	f.write("%i\n" %( count1+count2+count4 ) ) 

	for i in range(count1):
		f.write("%s %s 1 %s %s 1.0\n" %(array1[4*i], array1[4*i+1], array1[4*i+2], array1[4*i+3] ))
	for i in range(count2):
		f.write("%s %s 2 %s %s %s %s %s %s\n" %(array2_int[6*i], array2_int[6*i+1], array2_int[6*i+2], array2_int[6*i+3], array2_float[2*i], array2_int[6*i+4], array2_int[6*i+5], array2_float[2*i+1] ))
	for i in range(count4):
		f.write("%s %s 4 %s %s %s %s %s %s %s %s %s %s %s %s\n" %(array4_int[10*i], array4_int[10*i+1], array4_int[10*i+2], array4_int[10*i+3], array4_float[4*i], array4_int[10*i+4], array4_int[10*i+5], array4_float[4*i+1], array4_int[10*i+6], array4_int[10*i+7], array4_float[4*i+2], array4_int[10*i+8], array4_int[10*i+9], array4_float[4*i+3]))

	f.close()

	del(array1)
	del(array2_int)
	del(array2_float)
	del(array4_int)
	del(array4_float)

	return 1

def print_array_test_06_06_22(array):
	for num in range(5):
		print("%i	%s	%s" %(num, array[num,0], array[num,1]))
	return 1

def to_file_array_test_06_06_22(file_name, array):
	f = open (file_name, 'at')
	for num in range(5):
		f.write("%i	%s	%s\n" %(num, array[num,0], array[num,1]))
	f.close()
	return 1

#origin у повернутых сеток дается в неповернутой системе координат, как в rect
def gen_interpolation_file_from_rotated_rect_to_rect_2D(file_name, grid_from_id, grid_from_limits, eps, origin_from, spacing_from, angle, point, origin_to, size_to, spacing_to):
	f = open(file_name, "wt")
	f.write(grid_from_id + "\n")
	f.write("2D\n")

	interpolation_radius = get_interpolation_radius()
	
	eps_sq = eps*eps
	eps_99 = 0.0 - 0.99*eps
	
	if (grid_from_limits[4]==0):
		step_x = 1
	else:
		step_x = grid_from_limits[4]
		
	if (grid_from_limits[5]==0):
		step_y = 1
	else:
		step_y = grid_from_limits[5]

	array = np.zeros((5,2), dtype = 'float')
	r_array = np.zeros((5,2), dtype = 'float')
	dis_sq = np.zeros((4), dtype = 'float')

	array1 = []
	array2_int = []
	array2_float = []
	array4_int = []
	array4_float = []
	checked = np.zeros((size_to[0]+4,size_to[1]+4), dtype = 'bool')
	
#TODO Можно оптимизировать. Так как много вычислений синусов и косинуов для поворотов, а угол один - можно один раз это делать.
#Для rect возможна та же самая оптимизация.
	for I in range(grid_from_limits[0], grid_from_limits[1], step_x):
		for J in range(grid_from_limits[2], grid_from_limits[3], step_y):
		
			for num in range(1,5):
				get_rect_coords(num, ind_I(I,num), ind_J(J,num), array, origin_from, spacing_from)
				r_array[num,0], r_array[num,1] = rotated_2D(array[num,0], array[num,1], angle, point[0], point[1])
		
			Xmin = min(r_array[1,0],r_array[2,0],r_array[3,0],r_array[4,0])
			Xmax = max(r_array[1,0],r_array[2,0],r_array[3,0],r_array[4,0])		
			Ymin = min(r_array[1,1],r_array[2,1],r_array[3,1],r_array[4,1])
			Ymax = max(r_array[1,1],r_array[2,1],r_array[3,1],r_array[4,1])	

			i_start = (int)((Xmin-origin_to[0])/spacing_to[0]+0.5) - interpolation_radius
			i_end = (int)((Xmax-origin_to[0])/spacing_to[0]+0.5) + interpolation_radius + 1
			j_start = (int)((Ymin-origin_to[1])/spacing_to[1]+0.5) - interpolation_radius
			j_end = (int)((Ymax-origin_to[1])/spacing_to[1]+0.5) + interpolation_radius + 1
			
			i_start = max(-2,i_start)
			i_end = min(size_to[0]+2,i_end)
			j_start = max(-2,j_start)
			j_end = min(size_to[1]+2,j_end)
			
			for i in range(i_start,i_end):
				for j in range(j_start,j_end):
					if (checked[i+2,j+2]):
						continue
				
					get_rect_coords(0, i, j, r_array, origin_to, spacing_to)		
					array[0,0], array[0,1] = rotated_2D(r_array[0,0], r_array[0,1], 0.0 - angle, point[0], point[1])
		
					if ((array[0,0] - array[1,0]) < eps_99):
						continue

					if ((array[2,0] - array[0,0]) < eps_99):
						continue

					if ((array[0,1] - array[1,1]) < eps_99):
						continue

					if ((array[4,1] - array[0,1]) < eps_99):
						continue

					continue_flag = 0
					for num in range(1,5):
						dis_sq[num-1] = get_dis_sqrt(0,num,array)
						if (dis_sq[num-1] < eps_sq):
							array1.append(i)
							array1.append(j)
							array1.append(ind_I(I,num))
							array1.append(ind_J(J,num))
							checked[i+2,j+2] = True
							continue_flag = 1
							break
					if (continue_flag == 1):
						continue
			
					continue_flag = 0
					for num in range(1,5):
						next_num = get_next_num(num)
						if (dis_side_rect(array, num, eps)):
							if ((num == 1) or (num == 3)):
								axis = 0
							else:
								axis = 1
							if (dis_sq[num-1]<dis_sq[next_num-1]):
								alpha = (array[0,axis]-array[num,axis])/spacing_from[axis]
								alpha = fabs(alpha)
								if (alpha>1.0):
									to_file_array_test_06_06_22("debug.txt", array) #TODO DEL
									to_file_array_test_06_06_22("debug.txt", r_array) #TODO DEL
									sys.exit("Bad interpolation coef 1 " + file_name)
								array2_int.append(i)
								array2_int.append(j)
								array2_int.append(ind_I(I,num))
								array2_int.append(ind_J(J,num))
								array2_float.append(1.0-alpha)
								array2_int.append(ind_I(I,next_num))
								array2_int.append(ind_J(J,next_num))
								array2_float.append(alpha)
							else:
								alpha = (array[0,axis]-array[next_num,axis])/spacing_from[axis]
								alpha = fabs(alpha)
								if (alpha>1.0):
									sys.exit("Bad interpolation coef 2 " + file_name)
								array2_int.append(i)
								array2_int.append(j)
								array2_int.append(ind_I(I,num))
								array2_int.append(ind_J(J,num))
								array2_float.append(alpha)
								array2_int.append(ind_I(I,next_num))
								array2_int.append(ind_J(J,next_num))
								array2_float.append(1.0-alpha)
							checked[i+2,j+2] = True
							continue_flag = 1
							break					
					if (continue_flag == 1):
						continue
				
			#if not inside_quadr(array):
			#	sys.exit("Not inside quadr " + file_name)
			
					if (dis_sq[0]<dis_sq[1]): # 1 and 2
						alpha_plus = (array[0,0] - array[1,0])/spacing_from[0]
						alpha = 1.0 - alpha_plus
					else:
						alpha = (array[2,0] - array[0,0])/spacing_from[0]
						alpha_plus = 1.0 - alpha
					if (bad_coef(alpha)):
						sys.exit("Bad interpolation coef 3 " + file_name)
				
					if (dis_sq[0]<dis_sq[3]): # 1 and 4
						beta_plus = (array[0,1] - array[1,1])/spacing_from[1]
						beta = 1.0 - beta_plus
					else:
						beta = (array[4,1] - array[0,1])/spacing_from[1]
						beta_plus = 1.0 - beta
					if (bad_coef(beta)):
						sys.exit("Bad interpolation coef 4 " + file_name)			
					array4_int.append(i)
					array4_int.append(j)
					array4_int.append(I)
					array4_int.append(J)
					array4_float.append(alpha*beta)
					array4_int.append(I+1)
					array4_int.append(J)
					array4_float.append(alpha_plus*beta)			
					array4_int.append(I)
					array4_int.append(J+1)
					array4_float.append(alpha*beta_plus)			
					array4_int.append(I+1)
					array4_int.append(J+1)
					array4_float.append(alpha_plus*beta_plus)
					checked[i+2,j+2] = True

	count1 = int(len(array1)/4)
	count2 = int(len(array2_float)/2)
	count4 = int(len(array4_float)/4)
	
	f.write("%i\n" %( count1+count2+count4 ) ) 

	for i in range(count1):
		f.write("%s %s 1 %s %s 1.0\n" %(array1[4*i], array1[4*i+1], array1[4*i+2], array1[4*i+3] ))
	for i in range(count2):
		f.write("%s %s 2 %s %s %s %s %s %s\n" %(array2_int[6*i], array2_int[6*i+1], array2_int[6*i+2], array2_int[6*i+3], array2_float[2*i], array2_int[6*i+4], array2_int[6*i+5], array2_float[2*i+1] ))
	for i in range(count4):
		f.write("%s %s 4 %s %s %s %s %s %s %s %s %s %s %s %s\n" %(array4_int[10*i], array4_int[10*i+1], array4_int[10*i+2], array4_int[10*i+3], array4_float[4*i], array4_int[10*i+4], array4_int[10*i+5], array4_float[4*i+1], array4_int[10*i+6], array4_int[10*i+7], array4_float[4*i+2], array4_int[10*i+8], array4_int[10*i+9], array4_float[4*i+3]))

	f.close()

	del(array1)
	del(array2_int)
	del(array2_float)
	del(array4_int)
	del(array4_float)
	del(checked)

	return 1


#c means the center of rotation
def rotated_2D(x, y, angle, cx, cy):
	dx = x - cx
	dy = y - cy
	return (cx + cos(angle)*dx - sin(angle)*dy ), (cy + sin(angle)*dx + cos(angle)*dy )

#curved

#size = np.zeros((2), dtype = 'int')
#origin = np.zeros((2), dtype = 'float')
#spacing = np.zeros((2), dtype = 'float')
def gen_interpolation_file_from_rect_to_curved_2D(file_name, grid_from_id, curved_limits, eps, points_x, points_y, size, origin, spacing):
	f = open(file_name, "wt")
	f.write(grid_from_id + "\n")
	f.write("2D\n")
	
	eps_sq = eps*eps
	eps_99 = 0.0 - 0.99*eps
	
	if (curved_limits[4]==0):
		step_x = 1
	else:
		step_x = curved_limits[4]
		
	if (curved_limits[5]==0):
		step_y = 1
	else:
		step_y = curved_limits[5]

	array = np.zeros((5,2), dtype = 'float')
	dis_sq = np.zeros((4), dtype = 'float')

	array1 = []
	array2_int = []
	array2_float = []
	array4_int = []
	array4_float = []
		
	for i in range(curved_limits[0], curved_limits[1]+1, step_x):
		for j in range(curved_limits[2], curved_limits[3]+1, step_y):
			dx = points_x[i,j] - origin[0]
			dy = points_y[i,j] - origin[1]

			#TODO without 0.99 and without the following double operations in these extreme cases
			if ((dx < eps_99)or(dy < eps_99)):
				continue
				
			I_flag = True
			J_flag = True
			
			if (dx < 0.0):
				I = 0
				I_flag = False
				
			if (dy < 0.0):
				J = 0
				J_flag = False
			
			if I_flag:
				dx_right = (size[0]-1)*spacing[0] - dx
				if (dx_right < eps_99):
					continue
				if (dx_right < eps):
					I = size[0] - 2
					I_flag = False
					
			if J_flag:
				dy_right = (size[1]-1)*spacing[1] - dy
				if (dy_right < eps_99):
					continue
				if (dy_right < eps):
					J = size[1] - 2
					J_flag = False
			
			if I_flag:
				I = (int)(dx/spacing[0])
				if (I<0):
					continue
				if (I>(size[0]-2)):
					continue
					
			if J_flag:
				J = (int)(dy/spacing[1])
				if (J<0):
					continue
				if (J>(size[1]-2)):
					continue						
		
			array[0,0] = points_x[i,j]
			array[0,1] = points_y[i,j]
			for num in range(1,5):
				get_rect_coords(num, ind_I(I,num), ind_J(J,num), array, origin, spacing)
			
			continue_flag = 0
			for num in range(1,5):
				dis_sq[num-1] = get_dis_sqrt(0,num,array)
				if (dis_sq[num-1] < eps_sq):
					array1.append(i)
					array1.append(j)
					array1.append(ind_I(I,num))
					array1.append(ind_J(J,num))
					continue_flag = 1
					break
			if (continue_flag == 1):
				continue
			
			continue_flag = 0
			for num in range(1,5):
				next_num = get_next_num(num)
				if (dis_side_rect(array, num, eps)):
					if ((num == 1) or (num == 3)):
						axis = 0
					else:
						axis = 1
					if (dis_sq[num-1]<dis_sq[next_num-1]):
						alpha = (array[0,axis]-array[num,axis])/spacing[axis]
						alpha = fabs(alpha)
						if (alpha>1.0):
							sys.exit("Bad interpolation coef " + file_name)
						array2_int.append(i)
						array2_int.append(j)
						array2_int.append(ind_I(I,num))
						array2_int.append(ind_J(J,num))
						array2_float.append(1.0-alpha)
						array2_int.append(ind_I(I,next_num))
						array2_int.append(ind_J(J,next_num))
						array2_float.append(alpha)
					else:
						alpha = (array[0,axis]-array[next_num,axis])/spacing[axis]
						alpha = fabs(alpha)
						if (alpha>1.0):
							sys.exit("Bad interpolation coef " + file_name)
						array2_int.append(i)
						array2_int.append(j)
						array2_int.append(ind_I(I,num))
						array2_int.append(ind_J(J,num))
						array2_float.append(alpha)
						array2_int.append(ind_I(I,next_num))
						array2_int.append(ind_J(J,next_num))
						array2_float.append(1.0-alpha)
					continue_flag = 1
					break					
			if (continue_flag == 1):
				continue
				
			#if not inside_quadr(array):
			#	sys.exit("Not inside quadr " + file_name)
			
			if (dis_sq[0]<dis_sq[1]): # 1 and 2
				alpha_plus = (array[0,0] - array[1,0])/spacing[0]
				alpha = 1.0 - alpha_plus
			else:
				alpha = (array[2,0] - array[0,0])/spacing[0]
				alpha_plus = 1.0 - alpha
			if (bad_coef(alpha)):
				sys.exit("Bad interpolation coef " + file_name)
				
			if (dis_sq[0]<dis_sq[3]): # 1 and 4
				beta_plus = (array[0,1] - array[1,1])/spacing[1]
				beta = 1.0 - beta_plus
			else:
				beta = (array[4,1] - array[0,1])/spacing[1]
				beta_plus = 1.0 - beta
			if (bad_coef(beta)):
				sys.exit("Bad interpolation coef " + file_name)			
			array4_int.append(i)
			array4_int.append(j)
			array4_int.append(I)
			array4_int.append(J)
			array4_float.append(alpha*beta)
			array4_int.append(I+1)
			array4_int.append(J)
			array4_float.append(alpha_plus*beta)			
			array4_int.append(I)
			array4_int.append(J+1)
			array4_float.append(alpha*beta_plus)			
			array4_int.append(I+1)
			array4_int.append(J+1)
			array4_float.append(alpha_plus*beta_plus)			

	count1 = int(len(array1)/4)
	count2 = int(len(array2_float)/2)
	count4 = int(len(array4_float)/4)
	
	f.write("%i\n" %( count1+count2+count4 ) ) 

	for i in range(count1):
		f.write("%s %s 1 %s %s 1.0\n" %(array1[4*i], array1[4*i+1], array1[4*i+2], array1[4*i+3] ))
	for i in range(count2):
		f.write("%s %s 2 %s %s %s %s %s %s\n" %(array2_int[6*i], array2_int[6*i+1], array2_int[6*i+2], array2_int[6*i+3], array2_float[2*i], array2_int[6*i+4], array2_int[6*i+5], array2_float[2*i+1] ))
	for i in range(count4):
		f.write("%s %s 4 %s %s %s %s %s %s %s %s %s %s %s %s\n" %(array4_int[10*i], array4_int[10*i+1], array4_int[10*i+2], array4_int[10*i+3], array4_float[4*i], array4_int[10*i+4], array4_int[10*i+5], array4_float[4*i+1], array4_int[10*i+6], array4_int[10*i+7], array4_float[4*i+2], array4_int[10*i+8], array4_int[10*i+9], array4_float[4*i+3]))

	f.close()

	del(array1)
	del(array2_int)
	del(array2_float)
	del(array4_int)
	del(array4_float)

	return 1

#curved_limits 4,5 does not matter
#eps_N - solve quadratic equation or use Newton's method. Value 0.3 is recommended, not more
#if_test_axis_x, if_test_axis_y are 0 or 1
def gen_interpolation_file_from_curved_to_rect_2D(file_name, grid_from_id, curved_limits, eps, eps_N, points_x, points_y, size, origin, spacing, if_checked_to_draw, checked_path):

	interpolation_radius = get_interpolation_radius()

	f = open(file_name, "wt")
	f.write(grid_from_id + "\n")
	f.write("2D\n")
	
	eps_sq = eps*eps
		
	array = np.zeros((5,2), dtype = 'float')
	dis_sq = np.zeros((4), dtype = 'float')
	height = np.zeros((4), dtype = 'float')
	length = np.zeros((4), dtype = 'float')
	checked = np.zeros((size[0]+4,size[1]+4), dtype = 'bool')
	line = [] #to write in file
	res = np.zeros((2), dtype = 'float')
	info = np.zeros((5), dtype = 'int')
	XY = np.zeros((9), dtype = 'float')
	vectors = np.zeros((4,2), dtype = 'float')
	vectors0 = np.zeros((4,2), dtype = 'float')
	flags_vectors = np.ones((4), dtype = 'bool')
				
	#filling checked, koefs
	count = 0
	for I in range(curved_limits[0], curved_limits[1]):
		for J in range(curved_limits[2], curved_limits[3]):
			info[0] = 0 #new quadrangle
			for i in range(4):
				flags_vectors[i] = True
		
			for num in range(1,5):
				array[num,0] = points_x[ind_I(I,num),ind_J(J,num)]
				array[num,1] = points_y[ind_I(I,num),ind_J(J,num)]
		
			Xmin = min(array[1,0],array[2,0],array[3,0],array[4,0])
			Xmax = max(array[1,0],array[2,0],array[3,0],array[4,0])		
			Ymin = min(array[1,1],array[2,1],array[3,1],array[4,1])
			Ymax = max(array[1,1],array[2,1],array[3,1],array[4,1])	

			i_start = (int)((Xmin-origin[0])/spacing[0]+0.5) - interpolation_radius
			i_end = (int)((Xmax-origin[0])/spacing[0]+0.5) + interpolation_radius + 1
			j_start = (int)((Ymin-origin[1])/spacing[1]+0.5) - interpolation_radius
			j_end = (int)((Ymax-origin[1])/spacing[1]+0.5) + interpolation_radius + 1
			
			i_start = max(-2,i_start)
			i_end = min(size[0]+2,i_end)
			j_start = max(-2,j_start)
			j_end = min(size[1]+2,j_end)
			
			for i in range(i_start,i_end):
				for j in range(j_start,j_end):
					if (checked[i+2,j+2]):
						continue
				
					get_rect_coords(0, i, j, array, origin, spacing)
		
					continue_flag = 0
					for num in range(1,5):
						dis_sq[num-1] = get_dis_sqrt_vectors(0,num,array,vectors0[num-1])
						if (dis_sq[num-1] < eps_sq):
							line.append("%s %s 1 %s %s 1.0\n" %(i, j, ind_I(I,num), ind_J(J,num) ))
							count = count + 1
							checked[i+2,j+2] = True
							continue_flag = 1
							break
					if (continue_flag == 1):
						continue
			
					continue_flag = 0
					for num in range(1,5):
						next_num = get_next_num(num)
						if (flags_vectors[num-1]):
							length[num-1] = sqrt(get_dis_sqrt_vectors(next_num, num, array, vectors[num-1]))
							flags_vectors[num-1] = False
						height[num-1] = fabs(area(array, num, next_num, 0)) / length[num-1]
						if (height[num-1] < eps):
							if (scalar_product(vectors0[num-1],vectors[num-1]) < 0.0):
								continue
							alpha = Linear_Interpolation_on_Line(dis_sq[num-1], height[num-1], length[num-1])
							if (alpha>1.0):
								continue_flag = 1
								break
							if (dis_sq[num-1]<dis_sq[next_num-1]):
								line.append("%s %s 2 %s %s %s %s %s %s\n" %(i, j, ind_I(I,num), ind_J(J,num), 1.0-alpha, ind_I(I,next_num), ind_J(J,next_num), alpha ))
							else:
								alpha = Linear_Interpolation_on_Line(dis_sq[next_num-1], height[num-1], length[num-1])
								if (alpha>1.0):
									sys.exit("2 Bad interpolation coef " + file_name)
								line.append("%s %s 2 %s %s %s %s %s %s\n" %(i, j, ind_I(I,num), ind_J(J,num), alpha, ind_I(I,next_num), ind_J(J,next_num), 1.0-alpha ))
							if ((sqrt(dis_sq[num-1])>length[num-1]) or (sqrt(dis_sq[next_num-1])>length[num-1])): #TODO investigate if it is needed to delete
								sys.exit("2 Bad interpolation coef " + file_name) #TODO investigate if it is needed to delete
							count = count + 1
							checked[i+2,j+2] = True
							continue_flag = 1
							break					
					if (continue_flag == 1):
						continue
				
					if not_inside_quadr_given_vectors(vectors, vectors0):
						continue
						
					if (Bilinear_Interpolation(res, array, vectors, vectors0, height, length, eps_sq, eps_N, info, XY) == -1):					
						sys.exit("Bad interpolation coef found using bilinear interpolation" + file_name)
					if (bad_coef(res[0]) or bad_coef(res[1])):
						sys.exit("Bad interpolation coef " + file_name)

					alpha_plus = res[0]
					alpha = 1.0 - alpha_plus
					beta_plus = res[1]
					beta = 1.0 - beta_plus

					line.append("%s %s 4 %s %s %s %s %s %s %s %s %s %s %s %s\n" %(i, j, ind_I(I,info[1]), ind_J(J,info[1]), alpha*beta, ind_I(I,info[2]), ind_J(J,info[2]), alpha_plus*beta, ind_I(I,info[4]), ind_J(J,info[4]), alpha*beta_plus, ind_I(I,info[3]), ind_J(J,info[3]), alpha_plus*beta_plus))
					count = count + 1
					checked[i+2,j+2] = True

#start testing
#	good_x = 1
#	good_y = 1

#	if (if_test_axis_x == 1):
#		break_flag = False
#		for j in range(size[1]+4):
#			total = 0
#			for i in range(size[0]+4):
#				if checked[i,j]:
#					total = total + 1
#				elif (total > 0):
#					if (total == 1):
#						good_x = 0
#						break_flag = True
#						break
#					else:
#						total = 0
#			
#			if break_flag:
#				break
#			if (total == 1):
#				good_x = 0
#				break		

#	if (if_test_axis_y == 1):
#		break_flag = False
#		for i in range(size[0]+4):
#			total = 0
#			for j in range(size[1]+4):
#				if checked[i,j]:
#					total = total + 1
#				elif (total > 0):
#					if (total == 1):
#						good_y = 0
#						break_flag = True
#						break
#					else:
#						total = 0
			
#			if break_flag:
#				break
#			if (total == 1):
#				good_y = 0
#				break	
	
#end testing	

	f.write("%i\n" %(count))
	for i in range(count):
		f.write(line[i])
	f.close()
	
	if (if_checked_to_draw == 1):
		checked_float = np.zeros((size[0]+4,size[1]+4), dtype = 'float')
		for i in range(size[0]+4):
			for j in range(size[1]+4):
				if (checked[i,j]):
					checked_float[i,j] = 1.0
		data_to_vtk(checked_path,checked_float,origin[0]-2.0*spacing[0],origin[1]-2.0*spacing[1],size[0]+4,size[1]+4,spacing[0],spacing[1])
		
	del(line)
	del(checked)
		
	return 1 #good_x * good_y

#TODO !!!!! не работает, расчет разваливается
#curved_limits 4,5 does not matter
#eps_N - solve quadratic equation or use Newton's method. Value 0.3 is recommended, not more
#checked = np.zeros((size[0]+4,size[1]+4), dtype = 'bool')
#size[0] = Nx_background
#size[1] = Nz_background
#origin[0] = origin_x_background
#origin[1] = origin_z_background
#spacing[0] = hx_background
#spacing[1] = hz_background
def gen_interpolation_file_from_fracture_grid_Y_to_rect_2_3D_given_checked(file_name, grid_from_id,   array_corners, dop_y, fracture_index_len_y, Ny, grid_y_start,   eps, eps_N,   points_x, points_z,   size, origin, spacing, checked):

	interpolation_radius = get_interpolation_radius()

	def Func_Y(y_index):
		return y_index + grid_y_start

	f = open(file_name, "wt")
	f.write(grid_from_id + "\n")
	f.write("2D\n")
	
	eps_sq = eps*eps
		
	array = np.zeros((5,2), dtype = 'float')
	dis_sq = np.zeros((4), dtype = 'float')
	height = np.zeros((4), dtype = 'float')
	length = np.zeros((4), dtype = 'float')
#	checked = np.zeros((size[0]+4,size[1]+4), dtype = 'bool')
	#to write in file
	line = []
	res = np.zeros((2), dtype = 'float')
	info = np.zeros((5), dtype = 'int')
	XY = np.zeros((9), dtype = 'float')
	vectors = np.zeros((4,2), dtype = 'float')
	vectors0 = np.zeros((4,2), dtype = 'float')
	flags_vectors = np.ones((4), dtype = 'bool')
	
	#filling checked, koefs
	count = 0
	for y_index in range(2, dop_y+1):
		if (y_index == dop_y):
			y_range = range(dop_y, dop_y + fracture_index_len_y + 1)
			delta_count = fracture_index_len_y + 1
		else:
			y_range = [y_index, Ny - y_index - 1]
			delta_count = 2

		for i in range(size[0]+4):
			for j in range(size[1]+4):
				checked[i,j] = False

		for i in range(5):
			info[i] = 0
		for i in range(4):
			flags_vectors[i] = True

		#J means K
		for I in range(0, array_corners[6]-1):
			for J in range(0, array_corners[7]-1):
		
				info[0] = 0 #new quadrangle
				for i in range(4):
					flags_vectors[i] = True
		
				for num in range(1,5):
					array[num,0] = points_x[ind_I(I,num),y_index,ind_J(J,num)]
					array[num,1] = points_z[ind_I(I,num),y_index,ind_J(J,num)]
		
				Xmin = min(array[1,0],array[2,0],array[3,0],array[4,0])
				Xmax = max(array[1,0],array[2,0],array[3,0],array[4,0])		
				Ymin = min(array[1,1],array[2,1],array[3,1],array[4,1])
				Ymax = max(array[1,1],array[2,1],array[3,1],array[4,1])	

				i_start = (int)((Xmin-origin[0])/spacing[0]+0.5) - interpolation_radius
				i_end = (int)((Xmax-origin[0])/spacing[0]+0.5) + interpolation_radius + 1
				j_start = (int)((Ymin-origin[1])/spacing[1]+0.5) - interpolation_radius
				j_end = (int)((Ymax-origin[1])/spacing[1]+0.5) + interpolation_radius + 1
			
				i_start = max(-2,i_start)
				i_end = min(size[0]+2,i_end)
				j_start = max(-2,j_start)
				j_end = min(size[1]+2,j_end)
			
				for i in range(i_start,i_end):
					for j in range(j_start,j_end):
						if (checked[i+2,j+2]):
							continue
				
						get_rect_coords(0, i, j, array, origin, spacing)
		
						continue_flag = 0
						for num in range(1,5):
							dis_sq[num-1] = get_dis_sqrt_vectors(0,num,array,vectors0[num-1])
							if (dis_sq[num-1] < eps_sq):
								for y_index_var in y_range:
									line.append("%i %i %i 1 %i %i %i 1.0\n" %(i, Func_Y(y_index_var), j, ind_I(I,num), y_index_var, ind_J(J,num) ))
								count = count + delta_count
								checked[i+2,j+2] = True
								continue_flag = 1
								break
						if (continue_flag == 1):
							continue
			
						continue_flag = 0
						for num in range(1,5):
							next_num = get_next_num(num)
							if (flags_vectors[num-1]):
								length[num-1] = sqrt(get_dis_sqrt_vectors(next_num, num, array, vectors[num-1]))
								flags_vectors[num-1] = False
							height[num-1] = fabs(area(array, num, next_num, 0)) / length[num-1]
							
							if (height[num-1] < eps):
								if (scalar_product(vectors0[num-1],vectors[num-1]) < 0.0):
									continue
								alpha = Linear_Interpolation_on_Line(dis_sq[num-1], height[num-1], length[num-1])
								if (alpha>1.0):
									continue_flag = 1
									break
								if (dis_sq[num-1]<dis_sq[next_num-1]):
									for y_index_var in y_range:
										line.append("%i %i %i 2 %i %i %i %s %i %i %i %s\n" %(i, Func_Y(y_index_var), j, ind_I(I,num), y_index_var, ind_J(J,num), 1.0-alpha, ind_I(I,next_num), y_index_var, ind_J(J,next_num), alpha ))
								else:
									alpha = Linear_Interpolation_on_Line(dis_sq[next_num-1], height[num-1], length[num-1])
									if (alpha>1.0):
										sys.exit("2 Bad interpolation coef " + file_name)
									for y_index_var in y_range:
										line.append("%i %i %i 2 %i %i %i %s %i %i %i %s\n" %(i, Func_Y(y_index_var), j, ind_I(I,num), y_index_var, ind_J(J,num), alpha, ind_I(I,next_num), y_index_var, ind_J(J,next_num), 1.0-alpha ))
								if ((sqrt(dis_sq[num-1])>length[num-1]) or (sqrt(dis_sq[next_num-1])>length[num-1])): #TODO investigate if it is needed to delete
									sys.exit("2 Bad interpolation coef " + file_name) #TODO investigate if it is needed to delete
								count = count + delta_count
								checked[i+2,j+2] = True
								continue_flag = 1
								break					
						if (continue_flag == 1):
							continue
				
						if not_inside_quadr_given_vectors(vectors, vectors0):
							continue
						
						if (Bilinear_Interpolation(res, array, vectors, vectors0, height, length, eps_sq, eps_N, info, XY) == -1):					
							sys.exit("Bad interpolation coef found using bilinear interpolation" + file_name)
						if (bad_coef(res[0]) or bad_coef(res[1])):
							sys.exit("Bad interpolation coef " + file_name)

						alpha_plus = res[0]
						alpha = 1.0 - alpha_plus
						beta_plus = res[1]
						beta = 1.0 - beta_plus
					
						for y_index_var in y_range:
							line.append("%i %i %i 4 %i %i %i %s %i %i %i %s %i %i %i %s %i %i %i %s\n" %(i, Func_Y(y_index_var), j, ind_I(I,info[1]), y_index_var, ind_J(J,info[1]), alpha*beta, ind_I(I,info[2]), y_index_var, ind_J(J,info[2]), alpha_plus*beta, ind_I(I,info[4]), y_index_var, ind_J(J,info[4]), alpha*beta_plus, ind_I(I,info[3]), y_index_var, ind_J(J,info[3]), alpha_plus*beta_plus))
						count = count + delta_count
						checked[i+2,j+2] = True	

	f.write("%i\n" %(count + 4 * array_corners[7] * array_corners[6]))

	for y_index in [0, 1, Ny - 2, Ny - 1]:
		for I in range(0, array_corners[6]):
			for J in range(0, array_corners[7]):	
				f.write("%i %i %i 1 %i %i %i 1.0\n" %(I + array_corners[8], Func_Y(y_index), J + array_corners[9], I, y_index, J))
		
	for i in range(count):
		f.write(line[i])

	f.close()
		
	del(line)
		
	return 1 


#TODO
def gen_interpolation_file_from_curved_to_rect_2D_rect_limits(file_name, grid_from_id, rect_limits, eps, eps_N, points_x, points_y, size, origin, spacing, if_checked_to_draw, checked_path, if_test_axis_x, if_test_axis_y):

	f = open(file_name, "wt")
	f.write(grid_from_id + "\n")
	f.write("2D\n")

	if (rect_limits[4]==0):
		i_step = 1
	else:
		i_step = rect_limits[4]
		
	if (rect_limits[5]==0):
		j_step = 1
	else:
		j_step = rect_limits[5]

	f.write("%i\n" %(((rect_limits[1]-rect_limits[0])//i_step+1)*((rect_limits[3]-rect_limits[2])//j_step+1)))

#HERE
	
	eps_sq = eps*eps
		
	array = np.zeros((5,2), dtype = 'float')
	dis_sq = np.zeros((4), dtype = 'float')
	height = np.zeros((4), dtype = 'float')
	length = np.zeros((4), dtype = 'float')
	checked = np.zeros((size[0]+4,size[1]+4), dtype = 'bool')
	line = [] #to write in file
	res = np.zeros((2), dtype = 'float')
	info = np.zeros((5), dtype = 'int')
	XY = np.zeros((9), dtype = 'float')
	vectors = np.zeros((4,2), dtype = 'float')
	vectors0 = np.zeros((4,2), dtype = 'float')
	flags_vectors = np.ones((4), dtype = 'bool')
				
	#filling checked, koefs
	count = 0
	for I in range(curved_limits[0], curved_limits[1]):
		for J in range(curved_limits[2], curved_limits[3]):
			info[0] = 0 #new quadrangle
			for i in range(4):
				flags_vectors[i] = True
		
			for num in range(1,5):
				array[num,0] = points_x[ind_I(I,num),ind_J(J,num)]
				array[num,1] = points_y[ind_I(I,num),ind_J(J,num)]
		
			Xmin = min(array[1,0],array[2,0],array[3,0],array[4,0])
			Xmax = max(array[1,0],array[2,0],array[3,0],array[4,0])		
			Ymin = min(array[1,1],array[2,1],array[3,1],array[4,1])
			Ymax = max(array[1,1],array[2,1],array[3,1],array[4,1])	

			i_start = (int)((Xmin-origin[0])/spacing[0]+0.5) - interpolation_radius
			i_end = (int)((Xmax-origin[0])/spacing[0]+0.5) + interpolation_radius + 1
			j_start = (int)((Ymin-origin[1])/spacing[1]+0.5) - interpolation_radius
			j_end = (int)((Ymax-origin[1])/spacing[1]+0.5) + interpolation_radius + 1
			
			i_start = max(-2,i_start)
			i_end = min(size[0]+2,i_end)
			j_start = max(-2,j_start)
			j_end = min(size[1]+2,j_end)
			
			for i in range(i_start,i_end):
				for j in range(j_start,j_end):
					if (checked[i+2,j+2]):
						continue
				
					get_rect_coords(0, i, j, array, origin, spacing)
		
					continue_flag = 0
					for num in range(1,5):
						dis_sq[num-1] = get_dis_sqrt_vectors(0,num,array,vectors0[num-1])
						if (dis_sq[num-1] < eps_sq):
							line.append("%s %s 1 %s %s 1.0\n" %(i, j, ind_I(I,num), ind_J(J,num) ))
							count = count + 1
							checked[i+2,j+2] = True
							continue_flag = 1
							break
					if (continue_flag == 1):
						continue
			
					continue_flag = 0
					for num in range(1,5):
						next_num = get_next_num(num)
						if (flags_vectors[num-1]):
							length[num-1] = sqrt(get_dis_sqrt_vectors(next_num, num, array, vectors[num-1]))
							flags_vectors[num-1] = False
						height[num-1] = fabs(area(array, num, next_num, 0)) / length[num-1]
						if (height[num-1] < eps):
							if (scalar_product(vectors0[num-1],vectors[num-1]) < 0.0):
								continue
							alpha = Linear_Interpolation_on_Line(dis_sq[num-1], height[num-1], length[num-1])
							if (alpha>1.0):
								continue_flag = 1
								break
							if (dis_sq[num-1]<dis_sq[next_num-1]):
								line.append("%s %s 2 %s %s %s %s %s %s\n" %(i, j, ind_I(I,num), ind_J(J,num), 1.0-alpha, ind_I(I,next_num), ind_J(J,next_num), alpha ))
							else:
								alpha = Linear_Interpolation_on_Line(dis_sq[next_num-1], height[num-1], length[num-1])
								if (alpha>1.0):
									sys.exit("2 Bad interpolation coef " + file_name)
								line.append("%s %s 2 %s %s %s %s %s %s\n" %(i, j, ind_I(I,num), ind_J(J,num), alpha, ind_I(I,next_num), ind_J(J,next_num), 1.0-alpha ))
							if ((sqrt(dis_sq[num-1])>length[num-1]) or (sqrt(dis_sq[next_num-1])>length[num-1])): #TODO investigate if it is needed to delete
								sys.exit("2 Bad interpolation coef " + file_name) #TODO investigate if it is needed to delete
							count = count + 1
							checked[i+2,j+2] = True
							continue_flag = 1
							break					
					if (continue_flag == 1):
						continue
				
					if not_inside_quadr_given_vectors(vectors, vectors0):
						continue
						
					if (Bilinear_Interpolation(res, array, vectors, vectors0, height, length, eps_sq, eps_N, info, XY) == -1):					
						sys.exit("Bad interpolation coef found using bilinear interpolation" + file_name)
					if (bad_coef(res[0]) or bad_coef(res[1])):
						sys.exit("Bad interpolation coef " + file_name)

					alpha_plus = res[0]
					alpha = 1.0 - alpha_plus
					beta_plus = res[1]
					beta = 1.0 - beta_plus

					line.append("%s %s 4 %s %s %s %s %s %s %s %s %s %s %s %s\n" %(i, j, ind_I(I,info[1]), ind_J(J,info[1]), alpha*beta, ind_I(I,info[2]), ind_J(J,info[2]), alpha_plus*beta, ind_I(I,info[4]), ind_J(J,info[4]), alpha*beta_plus, ind_I(I,info[3]), ind_J(J,info[3]), alpha_plus*beta_plus))
					count = count + 1
					checked[i+2,j+2] = True

#start testing
#	good_x = 1
#	good_y = 1

#	if (if_test_axis_x == 1):
#		break_flag = False
#		for j in range(size[1]+4):
#			total = 0
#			for i in range(size[0]+4):
#				if checked[i,j]:
#					total = total + 1
#				elif (total > 0):
#					if (total == 1):
#						good_x = 0
#						break_flag = True
#						break
#					else:
#						total = 0
#			
#			if break_flag:
#				break
#			if (total == 1):
#				good_x = 0
#				break		

#	if (if_test_axis_y == 1):
#		break_flag = False
#		for i in range(size[0]+4):
#			total = 0
#			for j in range(size[1]+4):
#				if checked[i,j]:
#					total = total + 1
#				elif (total > 0):
#					if (total == 1):
#						good_y = 0
#						break_flag = True
#						break
#					else:
#						total = 0
			
#			if break_flag:
#				break
#			if (total == 1):
#				good_y = 0
#				break	
	
#end testing	

	f.write("%i\n" %(count))
	for i in range(count):
		f.write(line[i])
	f.close()
	
	if (if_checked_to_draw == 1):
		checked_float = np.zeros((size[0]+4,size[1]+4), dtype = 'float')
		for i in range(size[0]+4):
			for j in range(size[1]+4):
				if (checked[i,j]):
					checked_float[i,j] = 1.0
		data_to_vtk(checked_path,checked_float,origin[0]-2.0*spacing[0],origin[1]-2.0*spacing[1],size[0]+4,size[1]+4,spacing[0],spacing[1])
		
	return 1 #good_x * good_y




def mesh_to_vtk_mesh_only(f,origin_x,origin_y,Nx,Ny,hx,hy):
	f.write("# vtk DataFile Version 2.0\n")
	f.write("Created by Python generator written by Alena V. Favroskaya\n")
	f.write("ASCII\n")
	f.write("DATASET STRUCTURED_POINTS\n")
	f.write("DIMENSIONS " + str(Nx) + " " + str(Ny) + " 1\n")
	f.write("ORIGIN " + str(origin_x) + " " + str(origin_y) + " 0.0\n")
	f.write("SPACING " + str(hx) + " " + str(hy) + " " + str(hx) + "\n")
	return 1

def data_to_vtk(path,array,origin_x,origin_y,Nx,Ny,hx,hy):
	# write header
	f = open(path, "wt")
	mesh_to_vtk_mesh_only(f,origin_x,origin_y,Nx,Ny,hx,hy)
	f.write("\n")
	f.write("POINT_DATA " + str(Nx * Ny * 1) + "\n")
	f.write("SCALARS vy float 1\n")
	f.write("LOOKUP_TABLE default\n")

	for j in range(Ny):
		for i in range(Nx):
			f.write("%s\n" %(array[i,j]))
	f.close()
	return 1

#array = np.zeros((5,2), dtype = 'float') 
#5 points, 2 coords, zero point is (x,y)
def inside_quadr(array):
	return  inside_triangle(array, 1, 2, 3) or  inside_triangle(array, 1, 4, 3)

def inside_triangle(array, num1, num2, num3):
	area12 = area(array, 0, num1, num2)
	area23 = area(array, 0, num2, num3)
	area31 = area(array, 0, num3, num1)
	return not (( (area12<0.0)or(area23<0.0)or(area31<0.0) ) and ( (area12>0.0)or(area23>0.0)or(area31>0.0) ))

def area(array, num1, num2, num3):
	return (array[num1,0] - array[num3,0]) * (array[num2,1] - array[num3,1]) - (array[num2,0] - array[num3,0]) * (array[num1,1] - array[num3,1])

def vector_product(vector1,vector2):
	return vector1[0]*vector2[1] - vector1[1]*vector2[0]

def scalar_product(vector1,vector2):
	return vector1[0]*vector2[0] + vector1[1]*vector2[1]

def not_inside_quadr_given_vectors(vectors, vectors0):
	sign0 = sign(vector_product(vectors[0], vectors0[0]))
	for i in [1,2,3]:
		if (sign0 != sign(vector_product(vectors[i], vectors0[i]))):
			return True
	return False

def get_rect_coords(num, I, J, array, origin, spacing):
	array[num,0] = origin[0] + I*spacing[0]
	array[num,1] = origin[1] + J*spacing[1]
	return 1

def ind_I(I,num):
	if ((num == 1) or (num == 4)):
		return I
	return I+1
		
def ind_J(J,num):
	if (num < 3):
		return J
	return J+1

def get_next_num(num):
	if (num == 4):
		return 1
	return num + 1

def get_dis_sqrt(num1,num2,array):
	dx = array[num1,0] - array[num2,0]
	dy = array[num1,1] - array[num2,1]
	return dx*dx + dy*dy

def get_dis_sqrt_vectors(num1,num2,array,vector):
	for i in [0,1]:
		vector[i] = array[num1,i] - array[num2,i]
	return vector[0]*vector[0] + vector[1]*vector[1]

def dis_side_rect(array, num, eps):
	if (num == 1):
		axis = 1
	elif (num == 2):
		axis = 0
	elif (num == 3):
		axis = 1
	else: #if (num == 4):
		axis = 0
	return fabs(array[0,axis] - array[num,axis]) < eps

def Linear_Interpolation_on_Line(dis_sq, height, length):
	return sqrt(dis_sq - height*height) / length

#info = np.zeros((5), dtype = 'int')
#XY = np.zeros((9), dtype = 'float')
#0 - Xa, 1 - Ya, 2 - Xb, 3 - Yb, 4 - Xab, 5 - Yab, 6 - a1, 7 - a2, 8 - S_ABD
def Bilinear_Interpolation(res, array, vectors, vectors0, height, length, eps_sq, eps_N, info, XY):
	if (info[0] == 0):
		len_sum_max = length[3] + length[0]
		info[1] = 1

		for i in range(3):
			len_sum_test = length[i] + length[i+1]
			if (len_sum_max < len_sum_test):
				info[1] = i+2
				len_sum_max = len_sum_test

		if (info[1] == 1):
			info[2] = 2
			info[3] = 3
			info[4] = 4

			XY[0] = vectors[0,0] #array[2,0] - array[1,0] #Xa
			XY[1] = vectors[0,1] #array[2,1] - array[1,1] #Ya
		
			XY[2] = 0.0 - vectors[3,0] #array[4,0] - array[1,0] #Xb
			XY[3] = 0.0 - vectors[3,1] #array[4,1] - array[1,1] #Yb

			XY[4] = vectors[1,0] + vectors[3,0]#array[1,0] - array[4,0] + array[3,0] - array[2,0] #Xab
			XY[5] = vectors[1,1] + vectors[3,1] #array[1,1] - array[4,1] + array[3,1] - array[2,1] #Yab

		elif (info[1] == 2):
			info[2] = 3
			info[3] = 4
			info[4] = 1

			XY[0] = vectors[1,0] #array[3,0] - array[2,0] #Xa
			XY[1] = vectors[1,1] #array[3,1] - array[2,1] #Ya
		
			XY[2] = 0.0 - vectors[0,0] #array[1,0] - array[2,0] #Xb
			XY[3] = 0.0 - vectors[0,1] #array[1,1] - array[2,1] #Yb

			XY[4] = vectors[0,0] + vectors[2,0] #array[2,0] - array[1,0] + array[4,0] - array[3,0] #Xab
			XY[5] = vectors[0,1] + vectors[2,1] #array[2,1] - array[1,1] + array[4,1] - array[3,1] #Yab

		elif (info[1] == 3):
			info[2] = 4
			info[3] = 1
			info[4] = 2

			XY[0] = vectors[2,0] #array[4,0] - array[3,0] #Xa
			XY[1] = vectors[2,1] #array[4,1] - array[3,1] #Ya
		
			XY[2] = 0.0 - vectors[1,0] #array[2,0] - array[3,0] #Xb
			XY[3] = 0.0 - vectors[1,1] #array[2,1] - array[3,1] #Yb

			XY[4] = vectors[1,0] + vectors[3,0] #array[3,0] - array[2,0] + array[1,0] - array[4,0] #Xab
			XY[5] = vectors[1,1] + vectors[3,1] #array[3,1] - array[2,1] + array[1,1] - array[4,1] #Yab

		else: #if (info[1] == 4):
			info[2] = 1
			info[3] = 2
			info[4] = 3

			XY[0] = vectors[3,0] #array[1,0] - array[4,0] #Xa
			XY[1] = vectors[3,1] #array[1,1] - array[4,1] #Ya
		
			XY[2] = 0.0 - vectors[2,0] #array[3,0] - array[4,0] #Xb
			XY[3] = 0.0 - vectors[2,1] #array[3,1] - array[4,1] #Yb

			XY[4] = vectors[0,0] + vectors[2,0] #array[2,0] - array[1,0] + array[4,0] - array[3,0] #Xab
			XY[5] = vectors[0,1] + vectors[2,1] #array[2,1] - array[1,1] + array[4,1] - array[3,1] #Yab

#			XY[0] = array[B,0] - array[A,0] #Xa
#			XY[1] = array[B,1] - array[A,1] #Ya
		
#			XY[2] = array[D,0] - array[A,0] #Xb
#			XY[3] = array[D,1] - array[A,1] #Yb

#			XY[4] = array[A,0] + array[C,0] - array[B,0] - array[D,0] #Xab
#			XY[5] = array[A,1] + array[C,1] - array[B,1] - array[D,1] #Yab

		info[0] = 1

	res[0] = height[info[4]-1]/(height[info[4]-1]+height[info[2]-1])
	res[1] = height[info[1]-1]/(height[info[1]-1]+height[info[3]-1])

	X = vectors0[info[1]-1,0]
	Y = vectors0[info[1]-1,1]
	
	if (test_point(res[0], res[1], XY, X, Y, eps_sq)): #close to parallelogram
		return 10
		
	alpha = find_alpha(res[1], XY, X, Y)
	if (test_point(alpha, res[1], XY, X, Y, eps_sq)): #trapezoid AB // DC or close to it
		res[0] = alpha
		return 12

	beta = find_beta(res[0], XY, X, Y)
	if (test_point(res[0], beta, XY, X, Y, eps_sq)): #trapezoid AD // BC or close to it
		res[1] = beta
		return 11

	if (info[0] == 1):
		XY[6] = XY[4] * XY[1] - XY[5] * XY[0] #a1 = Xab * Ya - Yab * Xa
		XY[7] = XY[4] * XY[3] - XY[5] * XY[2] #a2 = Xab * Yb - Yab * Xb
		XY[8] = XY[0] * XY[3] - XY[2] * XY[1] #Xa * Yb - Xb * Ya
		info[0] = 2

	#b1 = Xa * Yb - Xb * Ya + Xab * Y - Yab * X
	#b2 = Xb * Ya - Xa * Yb + Xab * Y - Yab * X

	b1 = XY[4] * Y - XY[5] * X
	b2 = b1 - XY[8]
	b1 = b1 + XY[8]
	
	k1 = XY[6]/b1
	k2 = XY[7]/b2
	
	if (fabs(k1)<fabs(k2)): 			
		if (fabs(k1)>eps_N): #quadratic equation using k2, eps < k1 < k2
			b2 = 0.5/k2 #b/2a		
			D = b2*b2 + (Y * XY[0] - X * XY[1]) / XY[7] #c2 = X * Ya - Y * Xa, c2/a2,
			#D = b2*b2 - c2
			
			if (D < 0.0):
				return -1
		
			D = sqrt(D)
			test1 = b2 + D	
			if bad_coef(test1): # beta1 is bad, beta2 should be good
				res[1] = b2 - D
				res[0] = find_alpha(res[1], XY, X, Y)
				return 211

			test2 = b2 - D
			if bad_coef(test2): # beta1 is good, and beta2 is bad
				res[1] = test1
				res[0] = find_alpha(res[1], XY, X, Y)
				return 212
			
			res[0] = find_alpha(test2, XY, X, Y) #beta1 is good and beta2 is good, may be beta1 = beta2, D is close to zero
			if (test_point(res[0], test2, XY, X, Y, eps_sq)):
				res[1] = test2
				return 213		
			
			res[0] = find_alpha(test1, XY, X, Y)
			if bad_coef(res[0]): #if alpha2, beta2 were good, previous case was. so alpha2 is bad. alpha1 should be good
				return -1
			
			res[1] = test1
			return 214			
		
		#Newton's method using k1, k1 < eps < k2
		c1 = (X * XY[3] - Y * XY[2]) / b1 #c1 = X * Yb - Y * Xb, c1/b1		
#		continue_flag = True
#		while (continue_flag):
		for i in range(100):
			k1_alpha0 = k1*res[0]
			res[0] = res[0] + (c1 - res[0] * (1.0 - k1_alpha0) ) / (1.0 - 2.0*k1_alpha0)
			res[1] = find_beta(res[0], XY, X, Y)
			if (test_point(res[0], res[1], XY, X, Y, eps_sq)):
				return 32
		return -1
		
	#The same for k2, k2 < k1
	if (fabs(k2)>eps_N): #quadratic equation using k1, eps < k2 < k1
		b1 = 0.5/k1 #b/2a		
		D = b1*b1 + (Y * XY[2] - X * XY[3]) / XY[6] #c1 = X * Yb - Y * Xb, c1/a1,
		#D = b1*b1 - c1
			
		if (D < 0.0):
			return -1
		
		D = sqrt(D)
		test1 = b1 + D	
		if bad_coef(test1): # alpha1 is bad, alpha2 should be good
			res[0] = b1 - D
			res[1] = find_beta(res[0], XY, X, Y)
			return 221

		test2 = b1 - D
		if bad_coef(test2): # alpha1 is good, and alpha2 is bad
			res[0] = test1
			res[1] = find_beta(res[0], XY, X, Y)
			return 222
			
		res[1] = find_beta(test2, XY, X, Y) #alpha1 is good and alpha2 is good, may be alpha1 = alpha2, D is close to zero
		if (test_point(test2, res[1], XY, X, Y, eps_sq)):
			res[0] = test2
			return 223		
			
		res[1] = find_beta(test1, XY, X, Y)
		if bad_coef(res[1]): #if alpha2, beta2 were good, previous case was. so beta2 is bad. beta1 should be good
			return -1
			
		res[0] = test1
		return 224	
		
	#Newton's method using k2, k2 < eps < k1
	c2 = (X * XY[1] - Y * XY[0]) / b2 #c2 = X * Ya - Y * Xa, c2/b2		
#	continue_flag = True
#	while (continue_flag):
	for i in range(100):
		k2_beta0 = k2*res[1]
		res[1] = res[1] + (c2 - res[1] * (1.0 - k2_beta0) ) / (1.0 - 2.0*k2_beta0)
		res[0] = find_alpha(res[1], XY, X, Y)
		if (test_point(res[0], res[1], XY, X, Y, eps_sq)):
			return 31
	return -1

def test_point(alpha, beta, XY, X, Y, eps_sq):
	alpha_beta = alpha*beta
	dX = get_XY_test(alpha, beta, alpha_beta, XY[0], XY[2], XY[4]) - X
	dY = get_XY_test(alpha, beta, alpha_beta, XY[1], XY[3], XY[5]) - Y	
	return ((dX*dX+dY*dY) < eps_sq)

def get_XY_test(alpha, beta, alpha_beta, XYa, XYb, XYab):
	return XYa*alpha + XYb*beta + XYab*alpha_beta

def find_alpha(beta, XY, X, Y):
	denominatorX = XY[0] + XY[4]*beta #Xa + Xab*beta
	denominatorY = XY[1] + XY[5]*beta #Ya + Yab*beta
	if (fabs(denominatorX) > fabs(denominatorY)):
		return (X - XY[2]*beta)/denominatorX #(X - Xb*beta)/denominatorX
	return (Y - XY[3]*beta)/denominatorY #(Y - Yb*beta)/denominatorY

def find_beta(alpha, XY, X, Y):
	denominatorX = XY[2] + XY[4]*alpha #Xb + Xab*alpha
	denominatorY = XY[3] + XY[5]*alpha #Yb + Yab*alpha
	if (fabs(denominatorX) > fabs(denominatorY)):
		return (X - XY[0]*alpha)/denominatorX #(X - Xa*alpha)/denominatorX
	return (Y - XY[1]*alpha)/denominatorY #(Y - Ya*alpha)/denominatorY

def bad_coef(coef):
	return (coef<0.0) or (coef>1.0)

#3D simple grid, TODO optimization

def generate_2_3D_grid(f_grid_path, f_points_x, f_points_y, f_Nx, f_Ny, f_Nz, origin_z, step_z):
	# write header
	f = open(f_grid_path, "wt")
	f.write("# vtk DataFile Version 3.0\n")
	f.write("Created by Python generator written by Alena V. Favroskaya\n")
	f.write("BINARY\n")
	f.write("DATASET STRUCTURED_GRID\n")
	f.write("DIMENSIONS " + str(f_Nx) + " " + str(f_Ny) + " " + str(f_Nz) + "\n")
	f.write("POINTS " + str(f_Nx * f_Ny * f_Nz) + " float\n")

	coords = np.zeros(3 * f_Nx * f_Ny * f_Nz, dtype=np.dtype(get_float_type_in_vtk_savers()))

	for k in range(f_Nz):
		coord_z = origin_z + step_z * k
		for i in range(f_Nx):
			for j in range(f_Ny):
				ind = 3 * (k * f_Nx * f_Ny + j * f_Nx + i)
				coords[ind] = f_points_x[i,j]
				coords[ind+1] = f_points_y[i,j]
				coords[ind+2] = coord_z
		
	coords.tofile(f)
	f.close()
	del(coords)
	return 1
	
def generate_2_3D_grid_Y(f_grid_path, f_points_x, f_points_z, f_Nx, f_Ny, f_Nz, origin_y, step_y):
	# write header
	f = open(f_grid_path, "wt")
	f.write("# vtk DataFile Version 3.0\n")
	f.write("Created by Python generator written by Alena V. Favroskaya\n")
	f.write("BINARY\n")
	f.write("DATASET STRUCTURED_GRID\n")
	f.write("DIMENSIONS " + str(f_Nx) + " " + str(f_Ny) + " " + str(f_Nz) + "\n")
	f.write("POINTS " + str(f_Nx * f_Ny * f_Nz) + " float\n")

	coords = np.zeros(3 * f_Nx * f_Ny * f_Nz, dtype=np.dtype(get_float_type_in_vtk_savers()))

	for j in range(f_Ny):
		coord_y = origin_y + step_y * j
		for i in range(f_Nx):
			for k in range(f_Nz):
				ind = 3 * (k * f_Nx * f_Ny + j * f_Nx + i)
				coords[ind] = f_points_x[i,k]
				coords[ind+1] = coord_y
				coords[ind+2] = f_points_z[i,k]
		
	coords.tofile(f)
	f.close()
	del(coords)
	return 1

def generate_2_3D_grid_Y_in_array(array_x, array_y, array_z, f_points_x, f_points_z, f_Nx, f_Ny, f_Nz, origin_y, step_y):
	for j in range(f_Ny):
		coord_y = origin_y + step_y * j
		for i in range(f_Nx):
			for k in range(f_Nz):
				array_x[i,j,k] = f_points_x[i,k]
				array_y[i,j,k] = coord_y
				array_z[i,j,k] = f_points_z[i,k]
	return 1


def generate_3D_curved_grid_from_21_array(grid_path, array_21, z20, x00, x01, y10, y11, Nx, Ny, Nz):
    points = np.zeros((Nx,Ny,Nz,3), dtype = 'float')
    
    hx = (x01 - x00)/((float)(Nx-1))
    hy = (y11 - y10)/((float)(Ny-1))
    Nzm1 = (float)(Nz-1)
    for i in range(Nx):
        x_l = x00 + i*hx
        for j in range(Ny):
            y_l = y10 + j*hy
            z21 = array_21[i,j]
            hz = (z21 - z20)/Nzm1
            for k in range(Nz):
                points[i,j,k,0] = x_l
                points[i,j,k,1] = y_l
                points[i,j,k,2] = z20 + k*hz
    
    save_to_vtk_3D(grid_path, points, Nx, Ny, Nz)
    
    del(points)
    return 1
    
def generate_3D_curved_grid_from_20_array(grid_path, array_20, z21, x00, x01, y10, y11, Nx, Ny, Nz):
    points = np.zeros((Nx,Ny,Nz,3), dtype = 'float')
    
    hx = (x01 - x00)/((float)(Nx-1))
    hy = (y11 - y10)/((float)(Ny-1))
    Nzm1 = (float)(Nz-1)
    for i in range(Nx):
        x_l = x00 + i*hx
        for j in range(Ny):
            y_l = y10 + j*hy
            z20 = array_20[i,j]
            hz = (z21 - z20)/Nzm1
            for k in range(Nz):
                points[i,j,k,0] = x_l
                points[i,j,k,1] = y_l
                points[i,j,k,2] = z20 + k*hz
    
    save_to_vtk_3D(grid_path, points, Nx, Ny, Nz)
    del(points)
    return 1
    
def generate_3D_curved_grid_from_20_and_21_arrays(grid_path, array_20, array_21, x00, x01, y10, y11, Nx, Ny, Nz):
    points = np.zeros((Nx,Ny,Nz,3), dtype = 'float')
    
    hx = (x01 - x00)/((float)(Nx-1))
    hy = (y11 - y10)/((float)(Ny-1))
    Nzm1 = (float)(Nz-1)
    for i in range(Nx):
        x_l = x00 + i*hx
        for j in range(Ny):
            y_l = y10 + j*hy
            z20 = array_20[i,j]
            z21 = array_21[i,j]
            hz = (z21 - z20)/Nzm1
            for k in range(Nz):
                points[i,j,k,0] = x_l
                points[i,j,k,1] = y_l
                points[i,j,k,2] = z20 + k*hz
    
    save_to_vtk_3D(grid_path, points, Nx, Ny, Nz)
    del(points)
    return 1

def generate_3D_curved_grid_from_21(grid_path, func_21, z20, x00, x01, y10, y11, Nx, Ny, Nz):
    points = np.zeros((Nx,Ny,Nz,3), dtype = 'float')
    
    hx = (x01 - x00)/((float)(Nx-1))
    hy = (y11 - y10)/((float)(Ny-1))
    Nzm1 = (float)(Nz-1)
    for i in range(Nx):
        x_l = x00 + i*hx
        for j in range(Ny):
            y_l = y10 + j*hy
            z21 = func_20(x_l, y_l)
            hz = (z21 - z20)/Nzm1
            for k in range(Nz):
                points[i,j,k,0] = x_l
                points[i,j,k,1] = y_l
                points[i,j,k,2] = z20 + k*hz
    
    save_to_vtk_3D(grid_path, points, Nx, Ny, Nz)
    del(points)
    return 1
    
def generate_3D_curved_grid_from_20(grid_path, func_20, z21, x00, x01, y10, y11, Nx, Ny, Nz):
    points = np.zeros((Nx,Ny,Nz,3), dtype = 'float')
    
    hx = (x01 - x00)/((float)(Nx-1))
    hy = (y11 - y10)/((float)(Ny-1))
    Nzm1 = (float)(Nz-1)
    for i in range(Nx):
        x_l = x00 + i*hx
        for j in range(Ny):
            y_l = y10 + j*hy
            z20 = func_20(x_l, y_l)
            hz = (z21 - z20)/Nzm1
            for k in range(Nz):
                points[i,j,k,0] = x_l
                points[i,j,k,1] = y_l
                points[i,j,k,2] = z20 + k*hz
    
    save_to_vtk_3D(grid_path, points, Nx, Ny, Nz)
    del(points)
    return 1
    
def generate_3D_curved_grid_from_20_and_21(grid_path, func_20, func_21, x00, x01, y10, y11, Nx, Ny, Nz):
    points = np.zeros((Nx,Ny,Nz,3), dtype = 'float')
    
    hx = (x01 - x00)/((float)(Nx-1))
    hy = (y11 - y10)/((float)(Ny-1))
    Nzm1 = (float)(Nz-1)
    for i in range(Nx):
        x_l = x00 + i*hx
        for j in range(Ny):
            y_l = y10 + j*hy
            z20 = func_20(x_l, y_l)
            z21 = func_21(x_l, y_l)
            hz = (z21 - z20)/Nzm1
            for k in range(Nz):
                points[i,j,k,0] = x_l
                points[i,j,k,1] = y_l
                points[i,j,k,2] = z20 + k*hz
    
    save_to_vtk_3D(grid_path, points, Nx, Ny, Nz)
    del(points)
    return 1

def save_to_vtk_3D_XYZ_arrays(f_grid_path, f_points_x, f_points_y, f_points_z, f_Nx, f_Ny, f_Nz):
	# write header
	f = open(f_grid_path, "wt")
	f.write("# vtk DataFile Version 3.0\n")
	f.write("Created by Python generator written by Alena V. Favroskaya\n")
	f.write("BINARY\n")
	f.write("DATASET STRUCTURED_GRID\n")
	f.write("DIMENSIONS " + str(f_Nx) + " " + str(f_Ny) + " " + str(f_Nz) + "\n")
	f.write("POINTS " + str(f_Nx * f_Ny * f_Nz) + " float\n")

	coords = np.zeros(3 * f_Nx * f_Ny * f_Nz, dtype=np.dtype(get_float_type_in_vtk_savers()))

	for i in range(f_Nx):
		for j in range(f_Ny):
			for k in range(f_Nz):
				ind = 3 * (k * f_Ny * f_Nx + j * f_Nx + i)
				coords[ind] = f_points_x[i,j,k]
				coords[ind+1] = f_points_y[i,j,k]
				coords[ind+2] = f_points_z[i,j,k]

	coords.tofile(f)
	f.close()
	del(coords)
	return 1
    
def save_to_vtk_3D(f_grid_path, f_points, f_Nx, f_Ny, f_Nz):
	# write header
	f = open(f_grid_path, "wt")
	f.write("# vtk DataFile Version 3.0\n")
	f.write("Created by Python generator written by Alena V. Favroskaya\n")
	f.write("BINARY\n")
	f.write("DATASET STRUCTURED_GRID\n")
	f.write("DIMENSIONS " + str(f_Nx) + " " + str(f_Ny) + " " + str(f_Nz) + "\n")
	f.write("POINTS " + str(f_Nx * f_Ny * f_Nz) + " float\n")

	coords = np.zeros(3 * f_Nx * f_Ny * f_Nz, dtype=np.dtype(get_float_type_in_vtk_savers()))

	for i in range(f_Nx):
		for j in range(f_Ny):
			for k in range(f_Nz):
				ind = 3 * (k * f_Ny * f_Nx + j * f_Nx + i)
				coords[ind] = f_points[i,j,k,0]
				coords[ind+1] = f_points[i,j,k,1]
				coords[ind+2] = f_points[i,j,k,2]

	coords.tofile(f)
	f.close()
    
	del(coords)
	return 1

def save_to_vtk_value_rect(path, value, origin_x, origin_y, Nx, Ny, hx, hy, label):
	header_to_vtk(path, origin_x, origin_y, Nx, Ny, hx, hy, label)
	f = open(path, "at")
	for j in range(Ny-1,-1,-1):
		for i in range(Nx):
			f.write("%s\n" %(value[i,j]))
	f.close()
	return 1
	
def mesh_to_vtk_mesh_only(f,origin_x, origin_y, Nx, Ny, hx, hy):
	f.write("# vtk DataFile Version 2.0\n")
	f.write("Created by Python generator written by Alena V. Favroskaya\n")
	f.write("ASCII\n")
	f.write("DATASET STRUCTURED_POINTS\n")
	f.write("DIMENSIONS " + str(Nx) + " " + str(Ny) + " 1\n")
	f.write("ORIGIN " + str(origin_x) + " " + str(origin_y) + " 0.0\n")
	f.write("SPACING " + str(hx) + " " + str(hy) + " " + str(hx) + "\n")
	return 1

def header_to_vtk(path, origin_x, origin_y, Nx, Ny, hx, hy, label):
	# write header
	f = open(path, "wt")
	mesh_to_vtk_mesh_only(f,origin_x,origin_y,Nx,Ny,hx,hy)
	f.write("\n")
	f.write("POINT_DATA " + str(Nx * Ny * 1) + "\n")
	f.write("SCALARS " + label + " float 1\n")
	f.write("LOOKUP_TABLE default\n")
	f.close()
	return 1

    
def save_to_vtk_3D_with_value(f_grid_path, f_points, f_Nx, f_Ny, f_Nz, f_value_name, f_values):
	# write header
	f = open(f_grid_path, "wt")
	f.write("# vtk DataFile Version 3.0\n")
	f.write("Created by Python generator written by Alena V. Favroskaya\n")
	f.write("BINARY\n")
	f.write("DATASET STRUCTURED_GRID\n")
	f.write("DIMENSIONS " + str(f_Nx) + " " + str(f_Ny) + " " + str(f_Nz) + "\n")
	f.write("POINTS " + str(f_Nx * f_Ny * f_Nz) + " float\n")

	coords = np.zeros(3 * f_Nx * f_Ny * f_Nz, dtype=np.dtype(get_float_type_in_vtk_savers()))

	for i in range(f_Nx):
		for j in range(f_Ny):
			for k in range(f_Nz):
				ind = 3 * (k * f_Ny * f_Nx + j * f_Nx + i)
				coords[ind] = f_points[i,j,k,0]
				coords[ind+1] = f_points[i,j,k,1]
				coords[ind+2] = f_points[i,j,k,2]

	coords.tofile(f)

	f.write("\n")
	f.write("POINT_DATA " + str(f_Nx * f_Ny * f_Nz) + "\n")
	f.write("SCALARS " + f_value_name + " float 1 \n")
	f.write("LOOKUP_TABLE "+ f_value_name +"_table \n")

	values_to_write = np.zeros(f_Nx * f_Ny * f_Nz, dtype=np.dtype(get_float_type_in_vtk_savers()))

	for i in range(f_Nx):
		for j in range(f_Ny):
			for k in range(f_Nz):
				ind = k * f_Ny * f_Nx + j * f_Nx + i
				values_to_write[ind] = f_values[i,j,k]

	values_to_write.tofile(f)
	f.close()
    
	del(coords)
	del(values_to_write)
	return 1
	
def save_to_vtk_2D_with_value(f_grid_path, f_points, f_Nx, f_Ny, f_value_name, f_values):
	# write header
	f = open(f_grid_path, "wt")
	f.write("# vtk DataFile Version 3.0\n")
	f.write("Created by Python generator written by Alena V. Favroskaya\n")
	f.write("BINARY\n")
	f.write("DATASET STRUCTURED_GRID\n")
	f.write("DIMENSIONS " + str(f_Nx) + " " + str(f_Ny) + " 1\n")
	f.write("POINTS " + str(f_Nx * f_Ny) + " float\n")

	coords = np.zeros(3 * f_Nx * f_Ny, dtype=np.dtype(get_float_type_in_vtk_savers()))

	for i in range(f_Nx):
		for j in range(f_Ny):
			ind = 3 * (j * f_Nx + i)
			coords[ind] = f_points[i,j,0]
			coords[ind+1] = f_points[i,j,1]

	coords.tofile(f)

	f.write("\n")
	f.write("POINT_DATA " + str(f_Nx * f_Ny) + "\n")
	f.write("SCALARS " + f_value_name + " float 1 \n")
	f.write("LOOKUP_TABLE "+ f_value_name +"_table \n")

	values_to_write = np.zeros(f_Nx * f_Ny, dtype=np.dtype(get_float_type_in_vtk_savers()))

	for i in range(f_Nx):
		for j in range(f_Ny):
			ind = j * f_Nx + i
			values_to_write[ind] = f_values[i,j]

	values_to_write.tofile(f)
	f.close()
    
	del(coords)
	del(values_to_write)
	return 1

# FRACTURE GENERATOR

def poly3(x_or_y, a, b, c, d):
	x_or_y_2 = x_or_y*x_or_y
	x_or_y_3 = x_or_y_2*x_or_y
	return a*x_or_y_3 + b*x_or_y_2 + c*x_or_y + d

def find_abcd(abcd, tan_a,x0,y0,x1,y1):
	matrix = np.zeros((4,4), dtype = 'float')
	matrix_inv = np.zeros((4,4), dtype = 'float')
	right_vector = np.zeros((4), dtype = 'float')
	
	x0_2 = x0*x0
	x0_3 = x0_2*x0
	x1_2 = x1*x1
	x1_3 = x1_2*x1
	
	matrix[0,0] = x0_3
	matrix[0,1] = x0_2
	matrix[0,2] = x0
	matrix[0,3] = 1.0
	right_vector[0] = y0

	matrix[1,0] = x1_3
	matrix[1,1] = x1_2
	matrix[1,2] = x1
	matrix[1,3] = 1.0
	right_vector[1] = y1
	
	matrix[2,0] = 3.0*x0_2
	matrix[2,1] = 2.0*x0
	matrix[2,2] = 1.0
	matrix[2,3] = 0.0
	right_vector[2] = tan_a

	matrix[3,0] = 3.0*x1_2
	matrix[3,1] = 2.0*x1
	matrix[3,2] = 1.0
	matrix[3,3] = 0.0
	right_vector[3] = 0.0
	
	matrix_inv = LA.inv(matrix)
	
	for i in range(4):
		abcd[i] = 0.0
		for j in range(4):
			abcd[i] = abcd[i] + matrix_inv[i,j]*right_vector[j]
	
	return 1

def func_delta_infinity(x):
	if (x<0.0):
		return 0.0
	if (x>1.0):
		return 1.0
	if (x<0.5):
		return 2.0*x*x
	d = 1.0 - x
	return 1.0 - 2.0*d*d

def func_delta(x):
	if (x<0.5):
		return 2.0*x*x
	d = 1.0 - x
	return 1.0 - 2.0*d*d

def func_delta_parabolic(x):
	return x * x

def get_coord(index,origin,h):
	return origin + index*h	

def test_hx_for_min(Nx, Ny, points_x, h_min):
	for j in range(Ny):
		for i in range(1,Nx):
			if ((points_x[i,j]-points_x[i-1,j])<h_min):
				return 1
	return 0
	
def test_hy_for_min(Nx, Ny, points_y, h_min):
	for i in range(Nx):
		for j in range(1,Ny):
			if ((points_y[i,j]-points_y[i,j-1])<h_min):
				return 1
	return 0

def test_line_for_min(N, line_x, line_y, alpha):
	for i in range((int)(N/2)):
		if (alpha_by_normalized_vector_product(line_x[i+1]-line_x[i], line_y[i+1]-line_y[i], line_x[i+2]-line_x[i+1], line_y[i+2]-line_y[i+1])>alpha):
			return 1
	return 0		

#проверка на alpha_min_degree не доабвлена, так как получается что для dop_y = 16 и адекватного угла трещин и так будет нормально.
def fracture_grid_2_3D_Y(points_x, points_y, points_z, fracture_index_len_y, fracture_index_start_y, array_corners, given_points_x, given_points_z, dop_y, hx, hy, hz, origin_x, origin_y, origin_z):
	Ny = fracture_index_len_y + 2 * dop_y + 1
	grid_y_start = fracture_index_start_y - dop_y
	for j in range(Ny):
		poz_y = origin_y + hy * (j + grid_y_start)
		for i in range(array_corners[6]):
			for k in range(array_corners[7]):
				points_y[i,j,k] = poz_y

	for j in range(fracture_index_len_y + 1):
		for i in range(array_corners[6]):
			for k in range(array_corners[7]):
				points_x[i,dop_y + j,k] = given_points_x[i,k]
				points_z[i,dop_y + j,k] = given_points_z[i,k]

	coef = 1.0 / ((float)(dop_y - 1))
	for i in range(array_corners[6]):
		poz_x = origin_x + (array_corners[8] + i) * hx
		for k in range(array_corners[7]):
			poz_z = origin_z + (array_corners[9] + k) * hz
			for j in [0, 1, Ny-2, Ny-1]:
				points_x[i,j,k] = poz_x
				points_z[i,j,k] = poz_z
			for j in range(2, dop_y):
				 delta_j = func_delta(coef * (j-1))
				 points_x[i,j,k] = poz_x + (given_points_x[i,k] - poz_x) * delta_j
				 points_z[i,j,k] = poz_z + (given_points_z[i,k] - poz_z) * delta_j
				 points_x[i,Ny - j - 1,k] = points_x[i,j,k] 
				 points_z[i,Ny - j - 1,k] = points_z[i,j,k] 

	return Ny, grid_y_start

#ATTENTION! If using dop grids, fracture_dop_x, fracture_dop_y should be >= 3

#fracture_topRight_y > fracture_bottomLeft_y
#fracture_topRight_x ? fracture_bottomLeft_x

#array_corners = np.zeros((12), dtype = 'int')
#if find eff step. array_corners[11] = fracture_axis

#array_corners = np.zeros((11), dtype = 'int')
#TODO comments about its components

def set_fracture_grid(filename, array_corners, fracture_bottomLeft_x, fracture_bottomLeft_y, fracture_topRight_x, fracture_topRight_y, fracture_dop_x, fracture_dop_y, origin_x, origin_y, hx, hy, if_smooth_step, if_resetting, h_min, alpha_min_degree, step, if_find_eff_step, if_dop_grids):
	epsilon = get_epsilon()

	if (fracture_topRight_y < fracture_bottomLeft_y):
		set_fracture_grid(filename, array_corners, fracture_topRight_x, fracture_topRight_y, fracture_bottomLeft_x, fracture_bottomLeft_y, fracture_dop_x, fracture_dop_y, origin_x, origin_y, hx, hy, if_smooth_step, if_resetting, h_min, alpha_min_degree, step, if_find_eff_step, if_dop_grids)
	
	dis_x_pos = fracture_topRight_x - fracture_bottomLeft_x
	if (dis_x_pos < 0.0):
		dis_x_pos = abs(dis_x_pos)
		inverse = 1
	else:
		inverse = 0 
	
	if ((fracture_topRight_y - fracture_bottomLeft_y) > dis_x_pos):
		fracture_axis = 1
	else:
		fracture_axis = 0	
	
	abcd_0 = np.zeros((4), dtype = 'float')
	abcd_1 = np.zeros((4), dtype = 'float')

	if (fracture_axis == 0):
		fracture_center_y_int = find_int_rounded((fracture_topRight_y + fracture_bottomLeft_y)/2.0-origin_y,hy)
		tan_a = (fracture_topRight_y - fracture_bottomLeft_y)/(fracture_topRight_x - fracture_bottomLeft_x)

		if (inverse == 0):
			fracture_Left_x_int = find_int_rounded(fracture_bottomLeft_x-origin_x,hx)
			fracture_Right_x_int = find_int_rounded(fracture_topRight_x-origin_x,hx)
			
			start_x = fracture_bottomLeft_x
			start_y = fracture_bottomLeft_y
			end_x = fracture_topRight_x
			end_y = fracture_topRight_y		

		else:
			fracture_Left_x_int = find_int_rounded(fracture_topRight_x-origin_x,hx)
			fracture_Right_x_int = find_int_rounded(fracture_bottomLeft_x-origin_x,hx)

			end_x = fracture_bottomLeft_x
			end_y = fracture_bottomLeft_y
			start_x = fracture_topRight_x
			start_y = fracture_topRight_y
			
		fracture_len_x_int = fracture_Right_x_int - fracture_Left_x_int

		delta_x_left = start_x - get_coord(fracture_Left_x_int,origin_x,hx)			
		delta_x_right = end_x - get_coord(fracture_Right_x_int,origin_x,hx)			

		fracture_dop_x_real = fracture_dop_x
		bad_grid = 1
		while (bad_grid>0):
			x0 = start_x
			y0 = start_y
			x1 = get_coord(fracture_Left_x_int-fracture_dop_x_real,origin_x,hx)
			y1 = get_coord(fracture_center_y_int,origin_y,hy)
			find_abcd(abcd_0,tan_a,x0,y0,x1,y1)
		
			x0 = end_x
			y0 = end_y
			x1 = get_coord(fracture_Right_x_int+fracture_dop_x_real,origin_x,hx)
			y1 = get_coord(fracture_center_y_int,origin_y,hy)
			find_abcd(abcd_1,tan_a,x0,y0,x1,y1)

			Nx = 2*fracture_dop_x_real + fracture_len_x_int + 3
		
			line_x = np.zeros((Nx), dtype = 'float')
			line_y = np.zeros((Nx), dtype = 'float')
		
			for i in range(2):
				line_x[i] = get_coord(fracture_Left_x_int - fracture_dop_x_real - 1 + i,origin_x,hx)
				line_y[i] = get_coord(fracture_center_y_int,origin_y,hy)
			current = 1
	
			for i in range(1,fracture_dop_x_real):
				delta_i = i/((float)(fracture_dop_x_real))
				line_x[current+i] = get_coord(fracture_Left_x_int - fracture_dop_x_real + i,origin_x,hx) + delta_x_left*delta_i
				line_y[current+i] = poly3(line_x[current+i], abcd_0[0], abcd_0[1], abcd_0[2], abcd_0[3])		
			current = current + fracture_dop_x_real
			dx = (end_x - start_x)/(float)(fracture_len_x_int)
			dy = (end_y - start_y)/(float)(fracture_len_x_int)
	
			line_x[current] = start_x
			line_y[current] = start_y
			for i in range(1,fracture_len_x_int):
				line_x[current+i] = line_x[current+i-1] + dx
				line_y[current+i] = line_y[current+i-1] + dy			
			current = current + fracture_len_x_int
			line_x[current] = end_x
			line_y[current] = end_y
	
			for i in range(1,fracture_dop_x_real):
				delta_i = (fracture_dop_x_real-i-1)/((float)(fracture_dop_x_real))
				line_x[current+i] = get_coord(fracture_Right_x_int + i,origin_x,hx) + delta_x_right*delta_i
				line_y[current+i] = poly3(line_x[current+i], abcd_1[0], abcd_1[1], abcd_1[2], abcd_1[3])		
			current = current + fracture_dop_x_real
		
			for i in range(2):
				line_x[current+i] = get_coord(fracture_Right_x_int + fracture_dop_x_real + i,origin_x,hx)
				line_y[current+i] = get_coord(fracture_center_y_int,origin_y,hy)

			if (if_resetting == 0):
				break
			bad_grid = test_line_for_min(Nx, line_x, line_y, alpha_min_degree*pi/180.0)
			if (bad_grid > 0):
				fracture_dop_x_real = fracture_dop_x_real + step	

		delta_x = np.zeros((Nx), dtype = 'float')
		delta_y = np.zeros((Nx), dtype = 'float')

		fracture_dop_y_real = fracture_dop_y + (int)(fracture_len_x_int*fabs(tan_a)/2.0) + 1
		bad_grid = 1
		while (bad_grid>0):
			Ny = 2*fracture_dop_y_real + 2 + 1
	
			points_x = np.zeros((Nx,Ny), dtype = 'float')
			points_y = np.zeros((Nx,Ny), dtype = 'float')
			for i in range(Nx):
				for j in range(Ny):
					points_x[i,j] = get_coord(fracture_Left_x_int - fracture_dop_x_real - 1 + i,origin_x,hx)
					points_y[i,j] = get_coord(fracture_center_y_int - fracture_dop_y_real - 1 + j,origin_y,hy)
		
			j = fracture_dop_y_real + 1
			for i in range(Nx):
				delta_x[i] = line_x[i] - points_x[i,j]
				delta_y[i] = line_y[i] - points_y[i,j]
				points_x[i,j] = line_x[i]
				points_y[i,j] = line_y[i]
		
			for j in range(fracture_dop_y_real):
				delta_j = j/((float)(fracture_dop_y_real))
				delta_jx = func_delta(delta_j)
				if (if_smooth_step == 1):
					delta_j = delta_jx
				for i in range(Nx):
					points_x[i,1+j] = points_x[i,1+j] + delta_jx * delta_x[i]
					points_y[i,1+j] = points_y[i,1+j] + delta_j * delta_y[i]
					points_x[i,Ny-2-j] = points_x[i,Ny-2-j] + delta_jx * delta_x[i]
					points_y[i,Ny-2-j] = points_y[i,Ny-2-j] + delta_j * delta_y[i]
			if (if_resetting == 0):
				array_corners[10] = test_hy_for_min(Nx, Ny, points_y, epsilon)
				break
			bad_grid = test_hy_for_min(Nx, Ny, points_y, h_min)
			if (bad_grid > 0):
				fracture_dop_y_real = fracture_dop_y_real + step
								
		array_corners[0] = fracture_dop_x_real + 1
		array_corners[1] = fracture_dop_y_real + 1
		array_corners[2] = fracture_dop_x_real + fracture_len_x_int + 1
		array_corners[3] = fracture_dop_y_real + 1
		array_corners[4] = fracture_dop_x_real
		array_corners[5] = fracture_dop_y_real - (int)(fracture_len_x_int*fabs(tan_a)/2.0) - 1	
	else:
		fracture_bottomLeft_y_int = find_int_rounded(fracture_bottomLeft_y-origin_y,hy)
		fracture_topRight_y_int = find_int_rounded(fracture_topRight_y-origin_y,hy)
		fracture_center_x_int = find_int_rounded((fracture_topRight_x + fracture_bottomLeft_x)/2.0-origin_x,hx)
		fracture_len_y_int = fracture_topRight_y_int - fracture_bottomLeft_y_int
		tan_a = (fracture_topRight_x - fracture_bottomLeft_x)/(fracture_topRight_y - fracture_bottomLeft_y)
	
		delta_y_bottom = fracture_bottomLeft_y - get_coord(fracture_bottomLeft_y_int,origin_y,hy)			
		delta_y_top = fracture_topRight_y - get_coord(fracture_topRight_y_int,origin_y,hy)			

		fracture_dop_y_real = fracture_dop_y
		bad_grid = 1
		while (bad_grid>0):
			x0 = fracture_bottomLeft_y
			y0 = fracture_bottomLeft_x
			x1 = get_coord(fracture_bottomLeft_y_int-fracture_dop_y_real,origin_y,hy)
			y1 = get_coord(fracture_center_x_int,origin_x,hx)
			find_abcd(abcd_0,tan_a,x0,y0,x1,y1)
		
			x0 = fracture_topRight_y
			y0 = fracture_topRight_x
			x1 = get_coord(fracture_topRight_y_int+fracture_dop_y_real,origin_y,hy)
			y1 = get_coord(fracture_center_x_int,origin_x,hx)
			find_abcd(abcd_1,tan_a,x0,y0,x1,y1)

			Ny = 2*fracture_dop_y_real + fracture_len_y_int + 3

			line_x = np.zeros((Ny), dtype = 'float')
			line_y = np.zeros((Ny), dtype = 'float')

			for i in range(2):
				line_y[i] = get_coord(fracture_bottomLeft_y_int - fracture_dop_y_real - 1 + i,origin_y,hy)
				line_x[i] = get_coord(fracture_center_x_int,origin_x,hx)
			current = 1

			for i in range(1,fracture_dop_y_real):
				delta_i = i/((float)(fracture_dop_y_real))
				line_y[current+i] = get_coord(fracture_bottomLeft_y_int - fracture_dop_y_real + i,origin_y,hy) + delta_y_bottom*delta_i
				line_x[current+i] = poly3(line_y[current+i], abcd_0[0], abcd_0[1], abcd_0[2], abcd_0[3])		
			current = current + fracture_dop_y_real
			dx = (fracture_topRight_x - fracture_bottomLeft_x)/(float)(fracture_len_y_int)
			dy = (fracture_topRight_y - fracture_bottomLeft_y)/(float)(fracture_len_y_int)
	
			line_x[current] = fracture_bottomLeft_x
			line_y[current] = fracture_bottomLeft_y
			for i in range(1,fracture_len_y_int):
				line_x[current+i] = line_x[current+i-1] + dx
				line_y[current+i] = line_y[current+i-1] + dy			
			current = current + fracture_len_y_int
			line_x[current] = fracture_topRight_x
			line_y[current] = fracture_topRight_y
		
			for i in range(1,fracture_dop_y_real):
				delta_i = (fracture_dop_y_real-i-1)/((float)(fracture_dop_y_real))
				line_y[current+i] = get_coord(fracture_topRight_y_int + i,origin_y,hy) + delta_y_top*delta_i
				line_x[current+i] = poly3(line_y[current+i], abcd_1[0], abcd_1[1], abcd_1[2], abcd_1[3])		
			current = current + fracture_dop_y_real
	
			for i in range(2):
				line_y[current+i] = get_coord(fracture_topRight_y_int + fracture_dop_y_real + i,origin_y,hy)
				line_x[current+i] = get_coord(fracture_center_x_int,origin_x,hx)

			if (if_resetting == 0):
				break
			bad_grid = test_line_for_min(Ny, line_x, line_y, alpha_min_degree*pi/180.0)
			if (bad_grid > 0):
				fracture_dop_y_real = fracture_dop_y_real + step	

		delta_x = np.zeros((Ny), dtype = 'float')
		delta_y = np.zeros((Ny), dtype = 'float')

		fracture_dop_x_real = fracture_dop_x + (int)(fracture_len_y_int*fabs(tan_a)/2.0) + 1
		bad_grid = 1
		while (bad_grid>0):
			Nx = 2*fracture_dop_x_real + 2 + 1

			points_x = np.zeros((Nx,Ny), dtype = 'float')
			points_y = np.zeros((Nx,Ny), dtype = 'float')
			for i in range(Nx):
				for j in range(Ny):
					points_x[i,j] = get_coord(fracture_center_x_int - fracture_dop_x_real - 1 + i,origin_x,hx)
					points_y[i,j] = get_coord(fracture_bottomLeft_y_int - fracture_dop_y_real - 1 + j,origin_y,hy)
				
			i = fracture_dop_x_real + 1
			for j in range(Ny):
				delta_x[j] = line_x[j] - points_x[i,j]
				delta_y[j] = line_y[j] - points_y[i,j]
				points_x[i,j] = line_x[j]
				points_y[i,j] = line_y[j]	

			for i in range(fracture_dop_x_real):
				delta_i = i/((float)(fracture_dop_x_real))
				delta_iy = func_delta(delta_i)
				if (if_smooth_step == 1):
					delta_i = delta_iy
				for j in range(Ny):
					points_x[1+i,j] = points_x[1+i,j] + delta_i * delta_x[j]
					points_y[1+i,j] = points_y[1+i,j] + delta_iy * delta_y[j]
					points_x[Nx-2-i,j] = points_x[Nx-2-i,j] + delta_i * delta_x[j]
					points_y[Nx-2-i,j] = points_y[Nx-2-i,j] + delta_iy * delta_y[j]

			if (if_resetting == 0):
				array_corners[10] = test_hx_for_min(Nx, Ny, points_x, epsilon)
				break
			bad_grid = test_hx_for_min(Nx, Ny, points_x, h_min)
			if (bad_grid > 0):
				fracture_dop_x_real = fracture_dop_x_real + step

		array_corners[0] = fracture_dop_x_real + 1
		array_corners[1] = fracture_dop_y_real + 1
		array_corners[2] = fracture_dop_x_real + 1
		array_corners[3] = fracture_dop_y_real + fracture_len_y_int + 1
		array_corners[4] = fracture_dop_x_real - (int)(fracture_len_y_int*fabs(tan_a)/2.0) - 1
		array_corners[5] = fracture_dop_y_real	

	array_corners[6]=Nx
	array_corners[7]=Ny
	array_corners[8]=find_int(points_x[0,0]-origin_x,hx)
	array_corners[9]=find_int(points_y[0,0]-origin_y,hy)
	
	if (if_dop_grids == 0):
		save_to_vtk(filename, points_x, points_y, Nx, Ny)
	else:
		save_to_vtk(filename + ".vtk", points_x, points_y, Nx, Ny)
		fracture_dop_x_real = array_corners[0] - 1
		fracture_dop_y_real = array_corners[1] - 1
		if (fracture_axis == 0):
			fracture_len_x_int = array_corners[2] - array_corners[0]
			save_part_to_vtk(filename + "_A.vtk", points_x, points_y, fracture_dop_x_real + 1, fracture_dop_x_real + fracture_len_x_int + 1, fracture_dop_y_real - 3, fracture_dop_y_real + 1)
			save_part_to_vtk(filename + "_B.vtk", points_x, points_y, fracture_dop_x_real + 1, fracture_dop_x_real + fracture_len_x_int + 1, fracture_dop_y_real + 1, fracture_dop_y_real + 5)
		else:
			fracture_len_y_int = array_corners[3] - array_corners[1]		
			save_part_to_vtk(filename + "_A.vtk", points_x, points_y, fracture_dop_x_real - 3, fracture_dop_x_real + 1, fracture_dop_y_real + 1, fracture_dop_y_real + fracture_len_y_int + 1)
			save_part_to_vtk(filename + "_B.vtk", points_x, points_y, fracture_dop_x_real + 1, fracture_dop_x_real + 5, fracture_dop_y_real + 1, fracture_dop_y_real + fracture_len_y_int + 1)		
	
	if (if_find_eff_step == 1):
		eff_step = find_eff_step(points_x, points_y, array_corners[6], array_corners[7])
	
	if (if_find_eff_step == 1):
		array_corners[11]=fracture_axis
		return eff_step	
	
	del(points_x)
	del(points_y)
	return fracture_axis

def set_fracture_grid_in_array(points_x, points_y, array_corners, fracture_bottomLeft_x, fracture_bottomLeft_y, fracture_topRight_x, fracture_topRight_y, fracture_dop_x, fracture_dop_y, origin_x, origin_y, hx, hy, if_smooth_step, if_resetting, h_min, alpha_min_degree, step):

	epsilon = get_epsilon()
	
	dis_x_pos = fracture_topRight_x - fracture_bottomLeft_x
	if (dis_x_pos < 0.0):
		dis_x_pos = abs(dis_x_pos)
		inverse = 1
	else:
		inverse = 0 
	
	if ((fracture_topRight_y - fracture_bottomLeft_y) > dis_x_pos):
		fracture_axis = 1
	else:
		fracture_axis = 0	
	
	abcd_0 = np.zeros((4), dtype = 'float')
	abcd_1 = np.zeros((4), dtype = 'float')

	if (fracture_axis == 0):
		fracture_center_y_int = find_int_rounded((fracture_topRight_y + fracture_bottomLeft_y)/2.0-origin_y,hy)
		tan_a = (fracture_topRight_y - fracture_bottomLeft_y)/(fracture_topRight_x - fracture_bottomLeft_x)

		if (inverse == 0):
			fracture_Left_x_int = find_int_rounded(fracture_bottomLeft_x-origin_x,hx)
			fracture_Right_x_int = find_int_rounded(fracture_topRight_x-origin_x,hx)
			
			start_x = fracture_bottomLeft_x
			start_y = fracture_bottomLeft_y
			end_x = fracture_topRight_x
			end_y = fracture_topRight_y		

		else:
			fracture_Left_x_int = find_int_rounded(fracture_topRight_x-origin_x,hx)
			fracture_Right_x_int = find_int_rounded(fracture_bottomLeft_x-origin_x,hx)

			end_x = fracture_bottomLeft_x
			end_y = fracture_bottomLeft_y
			start_x = fracture_topRight_x
			start_y = fracture_topRight_y
			
		fracture_len_x_int = fracture_Right_x_int - fracture_Left_x_int

		delta_x_left = start_x - get_coord(fracture_Left_x_int,origin_x,hx)			
		delta_x_right = end_x - get_coord(fracture_Right_x_int,origin_x,hx)			

		fracture_dop_x_real = fracture_dop_x
		bad_grid = 1
		while (bad_grid>0):
			x0 = start_x
			y0 = start_y
			x1 = get_coord(fracture_Left_x_int-fracture_dop_x_real,origin_x,hx)
			y1 = get_coord(fracture_center_y_int,origin_y,hy)
			find_abcd(abcd_0,tan_a,x0,y0,x1,y1)
		
			x0 = end_x
			y0 = end_y
			x1 = get_coord(fracture_Right_x_int+fracture_dop_x_real,origin_x,hx)
			y1 = get_coord(fracture_center_y_int,origin_y,hy)
			find_abcd(abcd_1,tan_a,x0,y0,x1,y1)

			Nx = 2*fracture_dop_x_real + fracture_len_x_int + 3
		
			line_x = np.zeros((Nx), dtype = 'float')
			line_y = np.zeros((Nx), dtype = 'float')
		
			for i in range(2):
				line_x[i] = get_coord(fracture_Left_x_int - fracture_dop_x_real - 1 + i,origin_x,hx)
				line_y[i] = get_coord(fracture_center_y_int,origin_y,hy)
			current = 1
	
			for i in range(1,fracture_dop_x_real):
				delta_i = i/((float)(fracture_dop_x_real))
				line_x[current+i] = get_coord(fracture_Left_x_int - fracture_dop_x_real + i,origin_x,hx) + delta_x_left*delta_i
				line_y[current+i] = poly3(line_x[current+i], abcd_0[0], abcd_0[1], abcd_0[2], abcd_0[3])		
			current = current + fracture_dop_x_real
			dx = (end_x - start_x)/(float)(fracture_len_x_int)
			dy = (end_y - start_y)/(float)(fracture_len_x_int)
	
			line_x[current] = start_x
			line_y[current] = start_y
			for i in range(1,fracture_len_x_int):
				line_x[current+i] = line_x[current+i-1] + dx
				line_y[current+i] = line_y[current+i-1] + dy			
			current = current + fracture_len_x_int
			line_x[current] = end_x
			line_y[current] = end_y
	
			for i in range(1,fracture_dop_x_real):
				delta_i = (fracture_dop_x_real-i-1)/((float)(fracture_dop_x_real))
				line_x[current+i] = get_coord(fracture_Right_x_int + i,origin_x,hx) + delta_x_right*delta_i
				line_y[current+i] = poly3(line_x[current+i], abcd_1[0], abcd_1[1], abcd_1[2], abcd_1[3])		
			current = current + fracture_dop_x_real
		
			for i in range(2):
				line_x[current+i] = get_coord(fracture_Right_x_int + fracture_dop_x_real + i,origin_x,hx)
				line_y[current+i] = get_coord(fracture_center_y_int,origin_y,hy)

			if (if_resetting == 0):
				break
			bad_grid = test_line_for_min(Nx, line_x, line_y, alpha_min_degree*pi/180.0)
			if (bad_grid > 0):
				fracture_dop_x_real = fracture_dop_x_real + step	

		delta_x = np.zeros((Nx), dtype = 'float')
		delta_y = np.zeros((Nx), dtype = 'float')

		fracture_dop_y_real = fracture_dop_y + (int)(fracture_len_x_int*fabs(tan_a)/2.0) + 1
		bad_grid = 1
		while (bad_grid>0):
			Ny = 2*fracture_dop_y_real + 2 + 1
	
#			points_x = np.zeros((Nx,Ny), dtype = 'float')
#			points_y = np.zeros((Nx,Ny), dtype = 'float')
			for i in range(Nx):
				for j in range(Ny):
					points_x[i,j] = get_coord(fracture_Left_x_int - fracture_dop_x_real - 1 + i,origin_x,hx)
					points_y[i,j] = get_coord(fracture_center_y_int - fracture_dop_y_real - 1 + j,origin_y,hy)
		
			j = fracture_dop_y_real + 1
			for i in range(Nx):
				delta_x[i] = line_x[i] - points_x[i,j]
				delta_y[i] = line_y[i] - points_y[i,j]
				points_x[i,j] = line_x[i]
				points_y[i,j] = line_y[i]
		
			for j in range(fracture_dop_y_real):
				delta_j = j/((float)(fracture_dop_y_real))
				delta_jx = func_delta(delta_j)
				if (if_smooth_step == 1):
					delta_j = delta_jx
				for i in range(Nx):
					points_x[i,1+j] = points_x[i,1+j] + delta_jx * delta_x[i]
					points_y[i,1+j] = points_y[i,1+j] + delta_j * delta_y[i]
					points_x[i,Ny-2-j] = points_x[i,Ny-2-j] + delta_jx * delta_x[i]
					points_y[i,Ny-2-j] = points_y[i,Ny-2-j] + delta_j * delta_y[i]
			if (if_resetting == 0):
				array_corners[10] = test_hy_for_min(Nx, Ny, points_y, epsilon)
				break
			bad_grid = test_hy_for_min(Nx, Ny, points_y, h_min)
			if (bad_grid > 0):
				fracture_dop_y_real = fracture_dop_y_real + step
								
		array_corners[0] = fracture_dop_x_real + 1
		array_corners[1] = fracture_dop_y_real + 1
		array_corners[2] = fracture_dop_x_real + fracture_len_x_int + 1
		array_corners[3] = fracture_dop_y_real + 1
		array_corners[4] = fracture_dop_x_real
		array_corners[5] = fracture_dop_y_real - (int)(fracture_len_x_int*fabs(tan_a)/2.0) - 1	
	else:
		fracture_bottomLeft_y_int = find_int_rounded(fracture_bottomLeft_y-origin_y,hy)
		fracture_topRight_y_int = find_int_rounded(fracture_topRight_y-origin_y,hy)
		fracture_center_x_int = find_int_rounded((fracture_topRight_x + fracture_bottomLeft_x)/2.0-origin_x,hx)
		fracture_len_y_int = fracture_topRight_y_int - fracture_bottomLeft_y_int
		tan_a = (fracture_topRight_x - fracture_bottomLeft_x)/(fracture_topRight_y - fracture_bottomLeft_y)
	
		delta_y_bottom = fracture_bottomLeft_y - get_coord(fracture_bottomLeft_y_int,origin_y,hy)			
		delta_y_top = fracture_topRight_y - get_coord(fracture_topRight_y_int,origin_y,hy)			

		fracture_dop_y_real = fracture_dop_y
		bad_grid = 1
		while (bad_grid>0):
			x0 = fracture_bottomLeft_y
			y0 = fracture_bottomLeft_x
			x1 = get_coord(fracture_bottomLeft_y_int-fracture_dop_y_real,origin_y,hy)
			y1 = get_coord(fracture_center_x_int,origin_x,hx)
			find_abcd(abcd_0,tan_a,x0,y0,x1,y1)
		
			x0 = fracture_topRight_y
			y0 = fracture_topRight_x
			x1 = get_coord(fracture_topRight_y_int+fracture_dop_y_real,origin_y,hy)
			y1 = get_coord(fracture_center_x_int,origin_x,hx)
			find_abcd(abcd_1,tan_a,x0,y0,x1,y1)

			Ny = 2*fracture_dop_y_real + fracture_len_y_int + 3

			line_x = np.zeros((Ny), dtype = 'float')
			line_y = np.zeros((Ny), dtype = 'float')

			for i in range(2):
				line_y[i] = get_coord(fracture_bottomLeft_y_int - fracture_dop_y_real - 1 + i,origin_y,hy)
				line_x[i] = get_coord(fracture_center_x_int,origin_x,hx)
			current = 1

			for i in range(1,fracture_dop_y_real):
				delta_i = i/((float)(fracture_dop_y_real))
				line_y[current+i] = get_coord(fracture_bottomLeft_y_int - fracture_dop_y_real + i,origin_y,hy) + delta_y_bottom*delta_i
				line_x[current+i] = poly3(line_y[current+i], abcd_0[0], abcd_0[1], abcd_0[2], abcd_0[3])		
			current = current + fracture_dop_y_real
			dx = (fracture_topRight_x - fracture_bottomLeft_x)/(float)(fracture_len_y_int)
			dy = (fracture_topRight_y - fracture_bottomLeft_y)/(float)(fracture_len_y_int)
	
			line_x[current] = fracture_bottomLeft_x
			line_y[current] = fracture_bottomLeft_y
			for i in range(1,fracture_len_y_int):
				line_x[current+i] = line_x[current+i-1] + dx
				line_y[current+i] = line_y[current+i-1] + dy			
			current = current + fracture_len_y_int
			line_x[current] = fracture_topRight_x
			line_y[current] = fracture_topRight_y
		
			for i in range(1,fracture_dop_y_real):
				delta_i = (fracture_dop_y_real-i-1)/((float)(fracture_dop_y_real))
				line_y[current+i] = get_coord(fracture_topRight_y_int + i,origin_y,hy) + delta_y_top*delta_i
				line_x[current+i] = poly3(line_y[current+i], abcd_1[0], abcd_1[1], abcd_1[2], abcd_1[3])		
			current = current + fracture_dop_y_real
	
			for i in range(2):
				line_y[current+i] = get_coord(fracture_topRight_y_int + fracture_dop_y_real + i,origin_y,hy)
				line_x[current+i] = get_coord(fracture_center_x_int,origin_x,hx)

			if (if_resetting == 0):
				break
			bad_grid = test_line_for_min(Ny, line_x, line_y, alpha_min_degree*pi/180.0)
			if (bad_grid > 0):
				fracture_dop_y_real = fracture_dop_y_real + step	

		delta_x = np.zeros((Ny), dtype = 'float')
		delta_y = np.zeros((Ny), dtype = 'float')

		fracture_dop_x_real = fracture_dop_x + (int)(fracture_len_y_int*fabs(tan_a)/2.0) + 1
		bad_grid = 1
		while (bad_grid>0):
			Nx = 2*fracture_dop_x_real + 2 + 1

#			points_x = np.zeros((Nx,Ny), dtype = 'float')
#			points_y = np.zeros((Nx,Ny), dtype = 'float')
			for i in range(Nx):
				for j in range(Ny):
					points_x[i,j] = get_coord(fracture_center_x_int - fracture_dop_x_real - 1 + i,origin_x,hx)
					points_y[i,j] = get_coord(fracture_bottomLeft_y_int - fracture_dop_y_real - 1 + j,origin_y,hy)
				
			i = fracture_dop_x_real + 1
			for j in range(Ny):
				delta_x[j] = line_x[j] - points_x[i,j]
				delta_y[j] = line_y[j] - points_y[i,j]
				points_x[i,j] = line_x[j]
				points_y[i,j] = line_y[j]	

			for i in range(fracture_dop_x_real):
				delta_i = i/((float)(fracture_dop_x_real))
				delta_iy = func_delta(delta_i)
				if (if_smooth_step == 1):
					delta_i = delta_iy
				for j in range(Ny):
					points_x[1+i,j] = points_x[1+i,j] + delta_i * delta_x[j]
					points_y[1+i,j] = points_y[1+i,j] + delta_iy * delta_y[j]
					points_x[Nx-2-i,j] = points_x[Nx-2-i,j] + delta_i * delta_x[j]
					points_y[Nx-2-i,j] = points_y[Nx-2-i,j] + delta_iy * delta_y[j]

			if (if_resetting == 0):
				array_corners[10] = test_hx_for_min(Nx, Ny, points_x, epsilon)
				break
			bad_grid = test_hx_for_min(Nx, Ny, points_x, h_min)
			if (bad_grid > 0):
				fracture_dop_x_real = fracture_dop_x_real + step

		array_corners[0] = fracture_dop_x_real + 1
		array_corners[1] = fracture_dop_y_real + 1
		array_corners[2] = fracture_dop_x_real + 1
		array_corners[3] = fracture_dop_y_real + fracture_len_y_int + 1
		array_corners[4] = fracture_dop_x_real - (int)(fracture_len_y_int*fabs(tan_a)/2.0) - 1
		array_corners[5] = fracture_dop_y_real	

	array_corners[6]=Nx
	array_corners[7]=Ny
	array_corners[8]=find_int(points_x[0,0]-origin_x,hx)
	array_corners[9]=find_int(points_y[0,0]-origin_y,hy)

	return fracture_axis	

def set_dummy_fracture(filename, fracture_bottomLeft_x, fracture_bottomLeft_y, fracture_topRight_x, fracture_topRight_y):
	d_points_x = np.zeros((2), dtype = 'float')
	d_points_y = np.zeros((2), dtype = 'float')
	d_points_x[0] = fracture_bottomLeft_x
	d_points_y[0] = fracture_bottomLeft_y
	d_points_x[1] = fracture_topRight_x
	d_points_y[1] = fracture_topRight_y
	generate_dummy_grid_from_boundary_points(filename, 2, d_points_x, d_points_y)
	return


# FUNCS FOR TESTING

def print_points(f_points_x,f_points_y,f_Nx,f_Ny):
	for j in range(f_Ny-1,-1,-1):
		string = "#!	"
		for i in range(f_Nx):
			string = string + "(%.1f,%.1f)	" %(f_points_x[i,j]*100.0,f_points_y[i,j]*100.0)
		print(string)
	print()

#HELP FUNCTIONS
# grid_path and curr_corr_count are not needed in the function, used only to scientific interest to print grid transformation in process. the algorythm might be non critical optmized

def get_good_angle(f_norm,f_dx,f_dy): # return angle between the vector and dir (1,0) TODO test this functoin
	if comparison(f_norm, 0.0): 
		return 0.0
	else:
		if (f_dx>0.0): #dx>0, dy>0 OR #dx>0, dy<0 
			return asin(f_dy/f_norm)
		else:
			if (f_dy>0.0): #dx<0, dy>0
				return acos(f_dx/f_norm)
			else: #dx<0, dy<0
				return asin(f_dy/f_norm) - pi

def if_bad_cell(curr_corr_count, f_index, f_side_axis, f_points_x, f_points_y, f_Nx, f_Ny, f_zero_opposite):

	epsilon = get_epsilon()
	count_iteratons_critical = get_count_iteratons_critical()
	if_save_all_iterations = get_if_save_all_iterations()
	sin_critical = get_sin_critical()
	cos_critical = get_cos_critical()

	if (f_side_axis == "00"):
		if (f_zero_opposite):
			f_x_const = f_points_x[0,f_index-1] - f_points_x[0,f_index]
			f_y_const = f_points_y[0,f_index-1] - f_points_y[0,f_index]
		else:
			f_x_const = f_points_x[0,f_index+1] - f_points_x[0,f_index]
			f_y_const = f_points_y[0,f_index+1] - f_points_y[0,f_index]
		f_x_to_rotate = f_points_x[1,f_index] - f_points_x[0,f_index]
		f_y_to_rotate = f_points_y[1,f_index] - f_points_y[0,f_index]
	elif (f_side_axis == "01"):
		if (f_zero_opposite):
			f_x_const = f_points_x[f_Nx-1,f_index-1] - f_points_x[f_Nx-1,f_index]
			f_y_const = f_points_y[f_Nx-1,f_index-1] - f_points_y[f_Nx-1,f_index]
		else:
			f_x_const = f_points_x[f_Nx-1,f_index+1] - f_points_x[f_Nx-1,f_index]
			f_y_const = f_points_y[f_Nx-1,f_index+1] - f_points_y[f_Nx-1,f_index]
		f_x_to_rotate = f_points_x[f_Nx-2,f_index] - f_points_x[f_Nx-1,f_index]
		f_y_to_rotate = f_points_y[f_Nx-2,f_index] - f_points_y[f_Nx-1,f_index]
	elif (f_side_axis == "10"):
		if (f_zero_opposite):
			f_x_const = f_points_x[f_index-1,0] - f_points_x[f_index,0]
			f_y_const = f_points_y[f_index-1,0] - f_points_y[f_index,0]
		else:
			f_x_const = f_points_x[f_index+1,0] - f_points_x[f_index,0]
			f_y_const = f_points_y[f_index+1,0] - f_points_y[f_index,0]
		f_x_to_rotate = f_points_x[f_index,1] - f_points_x[f_index,0]
		f_y_to_rotate = f_points_y[f_index,1] - f_points_y[f_index,0]
	else: #11
		if (f_zero_opposite):
			f_x_const = f_points_x[f_index-1,f_Ny-1] - f_points_x[f_index,f_Ny-1]
			f_y_const = f_points_y[f_index-1,f_Ny-1] - f_points_y[f_index,f_Ny-1]
		else:
			f_x_const = f_points_x[f_index+1,f_Ny-1] - f_points_x[f_index,f_Ny-1]
			f_y_const = f_points_y[f_index+1,f_Ny-1] - f_points_y[f_index,f_Ny-1]
		f_x_to_rotate = f_points_x[f_index,f_Ny-2] - f_points_x[f_index,f_Ny-1]
		f_y_to_rotate = f_points_y[f_index,f_Ny-2] - f_points_y[f_index,f_Ny-1]

	len_const = sqrt(f_x_const*f_x_const + f_y_const*f_y_const)
	len_to_rotate = sqrt(f_x_to_rotate*f_x_to_rotate + f_y_to_rotate*f_y_to_rotate)
	cos_const = f_x_const/len_const
	sin_const = f_y_const/len_const
	cos_to_rotate = f_x_to_rotate/len_to_rotate
	sin_to_rotate = f_y_to_rotate/len_to_rotate

	sign = 1.0
	if (f_zero_opposite):
		if ((f_side_axis == "01")or(f_side_axis == "10")):
			sign = -1.0
	else:		
		if ((f_side_axis == "00")or(f_side_axis == "11")):
			sign = -1.0
	
	vector_product_normalized = sign*(cos_const*sin_to_rotate - sin_const*cos_to_rotate)
	if (vector_product_normalized < (sin_critical-epsilon)):

		f_vector_x = 0.0 - sign*sin_const*sin_critical
		f_vector_y = sign*cos_const*sin_critical

		if (cos_const*cos_to_rotate + sin_const*sin_to_rotate>0.0):
			f_vector_x = f_vector_x + cos_const*cos_critical
			f_vector_y = f_vector_y + sin_const*cos_critical
		else:
			f_vector_x = f_vector_x - cos_const*cos_critical
			f_vector_y = f_vector_y - sin_const*cos_critical

		f_vector_x = f_vector_x - cos_to_rotate 
		f_vector_y = f_vector_y - sin_to_rotate 	

		f_vector_x = f_vector_x*len_to_rotate 
		f_vector_y = f_vector_y*len_to_rotate 

		axis = (f_side_axis == "10") or (f_side_axis == "11")

		if axis: #10 or 11
			delta_n_y = (float)(f_Ny-2)
		else: #00 or 01
			delta_n_x = (float)(f_Nx-2)

		for i in range(1,f_Nx-1):
			if axis: #10 or 11
				if (i == f_index):
					delta_n_x = 1.0 #real 0
					diff_x = 0
				elif (i < f_index):
					delta_n_x = (float)(f_index)
					diff_x = f_index - i
				else: #i > f_index
					delta_n_x = (float)(f_Nx-f_index-1)
					diff_x = i - f_index
			else: #00 or 01
				if (f_side_axis == "00"):
					diff_x = i - 1
				else:
					diff_x = f_Nx-2-i

			for j in range(1,f_Ny-1):
				if axis: #10 or 11
					if (f_side_axis == "10"):
						diff_y = j - 1
					else:
						diff_y = f_Ny-2-j
				else: #00 or 01
					if (j == f_index):
						delta_n_y = 1.0 #real 0
						diff_y = 0
					elif (j < f_index):
						delta_n_y = (float)(f_index)
						diff_y = f_index - j
					else: #j > f_index
						delta_n_y = (float)(f_Ny-f_index-1)
						diff_y = j - f_index

				coef = (1.0-diff_x/delta_n_x)*(1.0-diff_y/delta_n_y)
				f_points_x[i,j] = f_points_x[i,j] + f_vector_x*coef
				f_points_y[i,j] = f_points_y[i,j] + f_vector_y*coef

#		if (if_save_all_iterations and count_iteratons_critical):
#			save_to_vtk(grid_path+"all_correction_%s.vtk" %(curr_corr_count), f_points_x,f_points_y,f_Nx,f_Ny)
#			print("Iteration %s, Side and Axis %s, Index %s, if_oppozite %s, vector_product_normalized %.3f, sin_critical %.3f" %(curr_corr_count, f_side_axis, f_index, f_zero_opposite, vector_product_normalized, sin_critical))
		return 1
	else:
		return 0

def sin_alpha_by_normalized_vector_product(f_x1, f_y1, f_x2, f_y2):
	norm1 = sqrt(f_x1*f_x1+f_y1*f_y1)
	norm2 = sqrt(f_x2*f_x2+f_y2*f_y2)	
	return (f_x1/norm1) * (f_y2/norm2) - (f_y1/norm1) * (f_x2/norm2)

def alpha_by_normalized_vector_product(f_x1, f_y1, f_x2, f_y2):
	return asin(sin_alpha_by_normalized_vector_product(f_x1, f_y1, f_x2, f_y2))

def test_normalized_vector_product(f_x1, f_y1, f_x2, f_y2):
	return ( sin_alpha_by_normalized_vector_product(f_x1, f_y1, f_x2, f_y2) < sin_critical )

def test_normalized_vector_product_critical(critical, f_x1, f_y1, f_x2, f_y2):
	return ( sin_alpha_by_normalized_vector_product(f_x1, f_y1, f_x2, f_y2) < critical )

def comparison(f_a,f_b): #1 if equal, 0 if not equal 
	epsilon = get_epsilon()

	if (f_a>f_b):
		if ((f_a-f_b) > epsilon):
			return 0
		else:
			return 1
	else:
		if ((f_b-f_a) > epsilon):
			return 0
		else:
			return 1
	return -1
            
def comparison_with_zero(f_a,f_epsilon): #1 if equal, 0 if not equal 
	epsilon = get_epsilon()

	if (f_a>0.0):
		if (f_a > f_epsilon):
			return 0
		else:
			return 1
	else:
		if ((0.0-f_a) > f_epsilon):
			return 0
		else:
			return 1
	return -1

def get_total_arc_len(f_arc_len, f_points_x_to_arc_len_calc, f_points_y_to_arc_len_calc, N_arc_length):
	f_total_arc_len_curr = 0.0
#	N_arc_length = get_N_arc_length()
	for k in range(N_arc_length-1):
		help_curr_x = f_points_x_to_arc_len_calc[k+1]-f_points_x_to_arc_len_calc[k]
		help_curr_y = f_points_y_to_arc_len_calc[k+1]-f_points_y_to_arc_len_calc[k]
		help_curr_x = help_curr_x*help_curr_x
		help_curr_y = help_curr_y*help_curr_y
		f_arc_len[k] = sqrt(help_curr_x + help_curr_y)
		f_total_arc_len_curr = f_total_arc_len_curr + f_arc_len[k]
	return f_total_arc_len_curr 

def find_scalar_interpolated(f_coef, f_scalarA, f_scalarB): # A if coef = 0, B if coef = 1
	return f_scalarA + f_coef * (f_scalarB - f_scalarA)

def save_to_vtk_line(f_grid_path, f_points_x, f_points_y, f_N):
	# write header
	f = open(f_grid_path, "wt")
	f.write("# vtk DataFile Version 3.0\n")
	f.write("Created by Python generator written by Alena V. Favroskaya\n")
	f.write("BINARY\n")
	f.write("DATASET STRUCTURED_GRID\n")
	f.write("DIMENSIONS " + str(f_N) + " 1 1\n")
	f.write("POINTS " + str(f_N * 1) + " float\n")

	coords = np.zeros(3 * f_N, dtype=np.dtype(get_float_type_in_vtk_savers()))

	for i in range(f_N):
		ind = 3 * i
		coords[ind] = f_points_x[i]
		coords[ind+1] = f_points_y[i]

	coords.tofile(f)
	f.close()
	del(coords)
	return 1
	
def save_to_vtk_line_one_array(f_grid_path, f_points, f_N):
	# write header
	f = open(f_grid_path, "wt")
	f.write("# vtk DataFile Version 3.0\n")
	f.write("Created by Python generator written by Alena V. Favroskaya\n")
	f.write("BINARY\n")
	f.write("DATASET STRUCTURED_GRID\n")
	f.write("DIMENSIONS " + str(f_N) + " 1 1\n")
	f.write("POINTS " + str(f_N * 1) + " float\n")

	coords = np.zeros(3 * f_N, dtype=np.dtype(get_float_type_in_vtk_savers()))

	for i in range(f_N):
		ind = 3 * i
		coords[ind] = f_points[i,0]
		coords[ind+1] = f_points[i,1]

	coords.tofile(f)
	f.close()
	del(coords)
	return 1

def save_to_vtk_one_array(f_grid_path, f_points, f_Nx, f_Ny):
	# write header
	f = open(f_grid_path, "wt")
	f.write("# vtk DataFile Version 3.0\n")
	f.write("Created by Python generator written by Alena V. Favroskaya\n")
	f.write("BINARY\n")
	f.write("DATASET STRUCTURED_GRID\n")
	f.write("DIMENSIONS " + str(f_Nx) + " " + str(f_Ny) + " 1\n")
	f.write("POINTS " + str(f_Nx * f_Ny * 1) + " float\n")

	coords = np.zeros(3 * f_Nx * f_Ny * 1, dtype=np.dtype(get_float_type_in_vtk_savers()))

	for i in range(f_Nx):
		for j in range(f_Ny):
			ind = 3 * (j * f_Nx + i)
			coords[ind] = f_points[0,i,j]
			coords[ind+1] = f_points[1,i,j]

	coords.tofile(f)
	f.close()
	del(coords)
	return 1

def save_to_vtk(f_grid_path, f_points_x, f_points_y, f_Nx, f_Ny):
	# write header
	f = open(f_grid_path, "wt")
	f.write("# vtk DataFile Version 3.0\n")
	f.write("Created by Python generator written by Alena V. Favroskaya\n")
	f.write("BINARY\n")
	f.write("DATASET STRUCTURED_GRID\n")
	f.write("DIMENSIONS " + str(f_Nx) + " " + str(f_Ny) + " 1\n")
	f.write("POINTS " + str(f_Nx * f_Ny * 1) + " float\n")

	coords = np.zeros(3 * f_Nx * f_Ny * 1, dtype=np.dtype(get_float_type_in_vtk_savers()))

	for i in range(f_Nx):
		for j in range(f_Ny):
			ind = 3 * (j * f_Nx + i)
			coords[ind] = f_points_x[i,j]
			coords[ind+1] = f_points_y[i,j]

	coords.tofile(f)
	f.close()
	del(coords)
	return 1

#(f_Nx_end - f_Nx_start + 1) * (f_Ny_end - f_Ny_start + 1)
def save_part_to_vtk(f_grid_path, f_points_x, f_points_y, f_Nx_start, f_Nx_end, f_Ny_start, f_Ny_end):
	# write header
	f = open(f_grid_path, "wt")
	f.write("# vtk DataFile Version 3.0\n")
	f.write("Created by Python generator written by Alena V. Favroskaya\n")
	f.write("BINARY\n")
	f.write("DATASET STRUCTURED_GRID\n")
	Nx_count = f_Nx_end - f_Nx_start + 1
	Ny_count = f_Ny_end - f_Ny_start + 1
	f.write("DIMENSIONS " + str(Nx_count) + " " + str(Ny_count) + " 1\n")
	total = Nx_count * Ny_count * 1
	f.write("POINTS " + str(total) + " float\n")
	coords = np.zeros(3 * total, dtype=np.dtype(get_float_type_in_vtk_savers()))

	ind = 0
	for j in range(f_Ny_start, f_Ny_end + 1):
		for i in range(f_Nx_start, f_Nx_end + 1):
			coords[ind] = f_points_x[i,j]
			ind = ind + 1
			coords[ind] = f_points_y[i,j]
			ind = ind + 2

	coords.tofile(f)
	f.close()
	del(coords)
	return 1

def save_vector_field_to_vtk(f_grid_path,f_grid_x,f_grid_y,f_vector_x,f_vector_y,f_Nx,f_Ny):
	# write header
	f = open(f_grid_path, "wt")
	f.write("# vtk DataFile Version 3.0\n")
	f.write("Created by Python generator written by Alena V. Favroskaya\n")
	f.write("BINARY\n")
	f.write("DATASET STRUCTURED_GRID\n")
	f.write("DIMENSIONS " + str(f_Nx) + " " + str(f_Ny) + " 1\n")
	f.write("POINTS " + str(f_Nx * f_Ny * 1) + " float\n")

	coords = np.zeros(3 * f_Nx * f_Ny * 1, dtype=np.dtype(get_float_type_in_vtk_savers()))

	for i in range(f_Nx):
		for j in range(f_Ny):
			ind = 3 * (j * f_Nx + i)
			coords[ind] = f_grid_x[i,j]
			coords[ind+1] = f_grid_y[i,j]

	coords.tofile(f)

	f.write("\n")
	f.write("POINT_DATA " + str(f_Nx * f_Ny * 1) + "\n")
	f.write("VECTORS grad float\n")

	vectors = np.zeros(3 * f_Nx * f_Ny * 1, dtype=np.dtype('float'))

	for i in range(f_Nx):
		for j in range(f_Ny):
			ind = 3 * (j * f_Nx + i)
			vectors[ind] = f_vector_x[i,j]
			vectors[ind+1] = f_vector_y[i,j]

	vectors.tofile(f)
	f.close()
	del(coords)
	del(vectors)
	return 1

def save_vector_field_to_vtk_regular_grid(f_grid_path,f_h,f_vector_x,f_vector_y,f_Nx,f_Ny):
	# write header
	f = open(f_grid_path, "wt")
	f.write("# vtk DataFile Version 3.0\n")
	f.write("Created by Python generator written by Alena V. Favroskaya\n")
	f.write("BINARY\n")
	f.write("\n")
	f.write("DATASET STRUCTURED_POINTS\n")
	f.write("DIMENSIONS " + str(f_Nx) + " " + str(f_Ny) + " 1\n")
	f.write("ORIGIN 0.0 0.0 0.0\n")
	f.write("SPACING " + str(f_h) + " " + str(f_h) + " " + str(f_h) +"\n")
	f.write("\n")
	f.write("POINT_DATA " + str(f_Nx * f_Ny * 1) + "\n")
	f.write("VECTORS grad float\n")

	vectors = np.zeros(3 * f_Nx * f_Ny * 1, dtype=np.dtype(get_float_type_in_vtk_savers()))

	for i in range(f_Nx):
		for j in range(f_Ny):
			ind = 3 * (j * f_Nx + i)
			vectors[ind] = f_vector_x[i,j]
			vectors[ind+1] = f_vector_y[i,j]

	vectors.tofile(f)
	f.close()
	del(vectors)
	return 1

#use min and max
#def find_min(f_a, f_b):
#	if (f_a > f_b):
#		return f_b
#	return f_a

#def find_max(f_a, f_b):
#	if (f_a > f_b):
#		return f_a
#	return f_b

#CURVED GRID GENERATOR
def generate_curved_grid(grid_path, Nx, Ny, X00, Y00, X01, Y01, X10, Y10, X11, Y11, t_start00, t_end00, t_start01, t_end01, t_start10, t_end10, t_start11, t_end11, grid_correction):

	points_x_00 = np.zeros((Ny), dtype = 'float')
	points_y_00 = np.zeros((Ny), dtype = 'float')
	points_x_01 = np.zeros((Ny), dtype = 'float')
	points_y_01 = np.zeros((Ny), dtype = 'float')
	points_x_10 = np.zeros((Nx), dtype = 'float')
	points_y_10 = np.zeros((Nx), dtype = 'float')
	points_x_11 = np.zeros((Nx), dtype = 'float')
	points_y_11 = np.zeros((Nx), dtype = 'float')

	success = generate_boundary_points(grid_path, Nx, Ny, X00, Y00, X01, Y01, X10, Y10, X11, Y11, t_start00, t_end00, t_start01, t_end01, t_start10, t_end10, t_start11, t_end11, grid_correction, points_x_00, points_y_00, points_x_01, points_y_01, points_x_10, points_y_10, points_x_11, points_y_11)
	if not success:
		print("ERROR boundary points failed in grid " + grid_path)
		return 0
	success = generate_curved_grid_from_boundary_points_separate_setting(grid_path, Nx, Ny, points_x_00, points_y_00, points_x_01, points_y_01, points_x_10, points_y_10, points_x_11, points_y_11, grid_correction)
	if not success:
		print("ERROR failed in generation of grid " + grid_path)
		return 0
	return 1

###### BOUNDARY POINTS SETTING
def generate_boundary_points(grid_path, Nx, Ny, X00, Y00, X01, Y01, X10, Y10, X11, Y11, t_start00, t_end00, t_start01, t_end01, t_start10, t_end10, t_start11, t_end11, grid_correction, points_x_00, points_y_00, points_x_01, points_y_01, points_x_10, points_y_10, points_x_11, points_y_11):

	N_arc_length = get_N_arc_length()
	coef_arc_len_critical = get_coef_arc_len_critical()

	to_return = 0 
	if not comparison(X00(t_start00),X10(t_start10)):
		print("ERROR DATA point A X 00 10 %s %s IN GRID " %(X00(t_start00),X10(t_start10)) + grid_path)
		to_return = 1
	if not comparison(X01(t_start01),X10(t_end10)):
		print("ERROR DATA point B X 01 10 %s %s IN GRID " %(X01(t_start01),X10(t_end10)) + grid_path)
		to_return = 1
	if not comparison(X01(t_end01),X11(t_end11)):
		print("ERROR DATA point C X 01 11 %s %s IN GRID " %(X01(t_end01),X11(t_end11)) + grid_path)
		to_return = 1
	if not comparison(X00(t_end00),X11(t_start11)):
		print("ERROR DATA point D X 00 11 %s %s IN GRID " %(X00(t_end00),X11(t_start11)) + grid_path)
		to_return = 1
	if not comparison(Y00(t_start00),Y10(t_start10)):
		print("ERROR DATA point A Y 00 10 %s %s IN GRID " %(Y00(t_start00),Y10(t_start10)) + grid_path)
		to_return = 1
	if not comparison(Y01(t_start01),Y10(t_end10)):
		print("ERROR DATA point B Y 01 10 %s %s IN GRID " %(Y01(t_start01),Y10(t_end10)) + grid_path)
		to_return = 1
	if not comparison(Y01(t_end01),Y11(t_end11)):
		print("ERROR DATA point C Y 01 11 %s %s IN GRID " %(Y01(t_end01),Y11(t_end11)) + grid_path)
		to_return = 1
	if not comparison(Y00(t_end00),Y11(t_start11)):
		print("ERROR DATA point D Y 00 11 %s %s IN GRID " %(Y00(t_end00),Y11(t_start11)) + grid_path)
		to_return = 1
	if (to_return):
		return 0

	points_x_to_arc_len_calc_00 = np.zeros((N_arc_length), dtype = 'float')
	points_y_to_arc_len_calc_00 = np.zeros((N_arc_length), dtype = 'float')
	points_x_to_arc_len_calc_01 = np.zeros((N_arc_length), dtype = 'float')
	points_y_to_arc_len_calc_01 = np.zeros((N_arc_length), dtype = 'float')
	points_x_to_arc_len_calc_10 = np.zeros((N_arc_length), dtype = 'float')
	points_y_to_arc_len_calc_10 = np.zeros((N_arc_length), dtype = 'float')
	points_x_to_arc_len_calc_11 = np.zeros((N_arc_length), dtype = 'float')
	points_y_to_arc_len_calc_11 = np.zeros((N_arc_length), dtype = 'float')

	points_x_00[0] = X00(t_start00)
	points_y_00[0] = Y00(t_start00)
	points_x_00[Ny-1] = X00(t_end00)
	points_y_00[Ny-1] = Y00(t_end00)

	points_x_01[0] = X01(t_start01)
	points_y_01[0] = Y01(t_start01)
	points_x_01[Ny-1] = X01(t_end01)
	points_y_01[Ny-1] = Y01(t_end01)

	points_x_10[0] = X10(t_start10)
	points_y_10[0] = Y10(t_start10)
	points_x_10[Nx-1] = X10(t_end10)
	points_y_10[Nx-1] = Y10(t_end10)

	points_x_11[0] = X11(t_start11)
	points_y_11[0] = Y11(t_start11)
	points_x_11[Nx-1] = X11(t_end11)
	points_y_11[Nx-1] = Y11(t_end11)

	if (Ny>2):
		arc_len_00 = np.zeros((N_arc_length-1), dtype = 'float')

		points_x_to_arc_len_calc_00[0] = points_x_00[0]
		points_y_to_arc_len_calc_00[0] = points_y_00[0]
		points_x_to_arc_len_calc_00[N_arc_length-1] = points_x_00[Ny-1]
		points_y_to_arc_len_calc_00[N_arc_length-1] = points_y_00[Ny-1]

		plus_curr_00 = (t_end00 - t_start00)/((float)(N_arc_length-1))
		t_curr = t_start00 + plus_curr_00
		for k in range(1,N_arc_length-1):
			points_x_to_arc_len_calc_00[k] = X00(t_curr)
			points_y_to_arc_len_calc_00[k] = Y00(t_curr)
			t_curr = t_curr + plus_curr_00

		total_arc_len00 = get_total_arc_len(arc_len_00, points_x_to_arc_len_calc_00, points_y_to_arc_len_calc_00, N_arc_length)

		arc_len_01 = np.zeros((N_arc_length-1), dtype = 'float')

		points_x_to_arc_len_calc_01[0] = points_x_01[0]
		points_y_to_arc_len_calc_01[0] = points_y_01[0]
		points_x_to_arc_len_calc_01[N_arc_length-1] = points_x_01[Ny-1]
		points_y_to_arc_len_calc_01[N_arc_length-1] = points_y_01[Ny-1]

		plus_curr_01 = (t_end01 - t_start01)/((float)(N_arc_length-1))
		t_curr = t_start01 + plus_curr_01
		for k in range(1,N_arc_length-1):
			points_x_to_arc_len_calc_01[k] = X01(t_curr)
			points_y_to_arc_len_calc_01[k] = Y01(t_curr)
			t_curr = t_curr + plus_curr_01

		total_arc_len01 = get_total_arc_len(arc_len_01, points_x_to_arc_len_calc_01, points_y_to_arc_len_calc_01, N_arc_length)

	if (Nx>2):
		arc_len_10 = np.zeros((N_arc_length-1), dtype = 'float')

		points_x_to_arc_len_calc_10[0] = points_x_10[0]
		points_y_to_arc_len_calc_10[0] = points_y_10[0]
		points_x_to_arc_len_calc_10[N_arc_length-1] = points_x_10[Nx-1]
		points_y_to_arc_len_calc_10[N_arc_length-1] = points_y_10[Nx-1]

		plus_curr_10 = (t_end10 - t_start10)/((float)(N_arc_length-1))
		t_curr = t_start10 + plus_curr_10
		for k in range(1,N_arc_length-1):
			points_x_to_arc_len_calc_10[k] = X10(t_curr)
			points_y_to_arc_len_calc_10[k] = Y10(t_curr)
			t_curr = t_curr + plus_curr_10
	
		total_arc_len10 = get_total_arc_len(arc_len_10, points_x_to_arc_len_calc_10, points_y_to_arc_len_calc_10, N_arc_length)

		points_x_to_arc_len_calc_11[0] = points_x_11[0]
		points_y_to_arc_len_calc_11[0] = points_y_11[0]
		points_x_to_arc_len_calc_11[N_arc_length-1] = points_x_11[Nx-1]
		points_y_to_arc_len_calc_11[N_arc_length-1] = points_y_11[Nx-1]

		arc_len_11 = np.zeros((N_arc_length-1), dtype = 'float')

		plus_curr_11 = (t_end11 - t_start11)/((float)(N_arc_length-1))
		t_curr = t_start11 + plus_curr_11
		for k in range(1,N_arc_length-1):
			points_x_to_arc_len_calc_11[k] = X11(t_curr)
			points_y_to_arc_len_calc_11[k] = Y11(t_curr)
			t_curr = t_curr + plus_curr_11

		total_arc_len11 = get_total_arc_len(arc_len_11, points_x_to_arc_len_calc_11, points_y_to_arc_len_calc_11, N_arc_length)

	if ((grid_correction) and (Nx>2) and (Ny>2)):
		#Find good points near A, B, C, D in accordance with critical angle

		#Point A, 00 and 10
		arc_len_one_step_ver = total_arc_len00/((float)(Ny-1))
		arc_len_delta_ver = coef_arc_len_critical*arc_len_one_step_ver
		t_curr_ver = t_start00
		total_arc_len_curr_ver_start = 0.0
		k_start_00 = 0

		arc_len_one_step_hor = total_arc_len10/((float)(Nx-1))
		arc_len_delta_hor = coef_arc_len_critical*arc_len_one_step_hor
		t_curr_hor = t_start10
		total_arc_len_curr_hor_start = 0.0
		k_start_10 = 0

		First_try_A = 1
		bad_point = 1
		while bad_point:
			#find test point, axis 00
			total_arc_len_curr_ver = total_arc_len_curr_ver_start
			less = 1 
			while less:
				total_arc_len_curr_ver_start = total_arc_len_curr_ver
				total_arc_len_curr_ver = total_arc_len_curr_ver + arc_len_00[k_start_00]
				diff_curr = total_arc_len_curr_ver - arc_len_one_step_ver
				if comparison(diff_curr, 0.0):
					arc_len_start_00 = total_arc_len_curr_ver
					total_arc_len_curr_new_start_00 = total_arc_len_curr_ver
					points_x_00[1] = points_x_to_arc_len_calc_00[k_start_00+1] 
					points_y_00[1] = points_y_to_arc_len_calc_00[k_start_00+1]
					less = 0
				elif ( diff_curr > 0.0 ):
					t_curr_new = t_curr_ver + plus_curr_00*(1.0-diff_curr/arc_len_00[k_start_00])
					arc_len_start_00 = total_arc_len_curr_ver - diff_curr
					points_x_00[1] = X00(t_curr_new)
					points_y_00[1] = Y00(t_curr_new) 
					less = 0
				else:
					if (k_start_00>(N_arc_length-3)):
						print("ERROR bad angle at point A. Try to use bigger N_arc_length. Current is %s. Grid" %(N_arc_length) + grid_path)
						return 0
					t_curr_ver = t_curr_ver + plus_curr_00
					k_start_00 = k_start_00 + 1
 

			#find test point, axis 10
			total_arc_len_curr_hor = total_arc_len_curr_hor_start
			less = 1 
			while less:
				total_arc_len_curr_hor_start = total_arc_len_curr_hor
				total_arc_len_curr_hor = total_arc_len_curr_hor + arc_len_10[k_start_10]
				diff_curr = total_arc_len_curr_hor - arc_len_one_step_hor
				if comparison(diff_curr, 0.0):
					arc_len_start_10 = total_arc_len_curr_hor
					points_x_10[1] = points_x_to_arc_len_calc_10[k_start_10+1] 
					points_y_10[1] = points_y_to_arc_len_calc_10[k_start_10+1]
					less = 0
				elif ( diff_curr > 0.0 ):
					t_curr_new = t_curr_hor + plus_curr_10*(1.0-diff_curr/arc_len_10[k_start_10])
					arc_len_start_10 = total_arc_len_curr_hor - diff_curr
					points_x_10[1] = X10(t_curr_new)
					points_y_10[1] = Y10(t_curr_new) 
					less = 0
				else:
					if (k_start_10>(N_arc_length-3)):
						print("ERROR bad angle at point A. Try to use bigger N_arc_length. Current is %s. Grid" %(N_arc_length) + grid_path)
						return 0
					t_curr_hor = t_curr_hor + plus_curr_10
					k_start_10 = k_start_10 + 1

			bad_point = test_normalized_vector_product(points_x_10[1]-points_x_10[0], points_y_10[1]-points_y_10[0], points_x_00[1]-points_x_00[0], points_y_00[1]-points_y_00[0])
			arc_len_one_step_ver = arc_len_one_step_ver + arc_len_delta_ver
			arc_len_one_step_hor = arc_len_one_step_hor + arc_len_delta_hor
			if First_try_A and bad_point:
				First_try_A = 0
			if (arc_len_one_step_ver>total_arc_len00):
				print("ERROR bad angle at point A. Try to use bigger N_arc_length. Current is %s. Grid" %(N_arc_length) + grid_path)
				return 0	
			if (arc_len_one_step_hor>total_arc_len10):
				print("ERROR bad angle at point A. Try to use bigger N_arc_length. Current is %s. Grid" %(N_arc_length) + grid_path)
				return 0

		total_arc_len_curr_new_start_00 = total_arc_len_curr_ver_start
		total_arc_len_curr_new_start_10 = total_arc_len_curr_hor_start
		t_new_start_00 = t_curr_ver
		t_new_start_10 = t_curr_hor


		#Point B, 01 and 10
		arc_len_one_step_ver = total_arc_len01/((float)(Ny-1))
		arc_len_delta_ver = coef_arc_len_critical*arc_len_one_step_ver
		t_curr_ver = t_start01
		total_arc_len_curr_ver_start = 0.0
		k_start_01 = 0

		arc_len_one_step_hor = (total_arc_len10-arc_len_start_10)/((float)(Nx-2))
		arc_len_delta_hor = coef_arc_len_critical*arc_len_one_step_hor
		t_curr_hor = t_end10
		total_arc_len_curr_hor_start = 0.0
		k_end_10 = 0

		First_try_B = 1
		bad_point = 1
		while bad_point:
			#find test point, axis 01
			total_arc_len_curr_ver = total_arc_len_curr_ver_start
			less = 1 
			while less:
				total_arc_len_curr_ver_start = total_arc_len_curr_ver
				total_arc_len_curr_ver = total_arc_len_curr_ver + arc_len_01[k_start_01]
				diff_curr = total_arc_len_curr_ver - arc_len_one_step_ver
				if comparison(diff_curr, 0.0):
					arc_len_start_01 = total_arc_len_curr_ver
					points_x_01[1] = points_x_to_arc_len_calc_01[k_start_01+1] 
					points_y_01[1] = points_y_to_arc_len_calc_01[k_start_01+1]
					less = 0
				elif ( diff_curr > 0.0 ):
					t_curr_new = t_curr_ver + plus_curr_01*(1.0-diff_curr/arc_len_01[k_start_01])
					arc_len_start_01 = total_arc_len_curr_ver - diff_curr
					points_x_01[1] = X01(t_curr_new)
					points_y_01[1] = Y01(t_curr_new) 
					less = 0
				else:
					if (k_start_01>(N_arc_length-3)):
						print("ERROR bad angle at point B. Try to use bigger N_arc_length. Current is %s. Grid" %(N_arc_length) + grid_path)
						return 0
					t_curr_ver = t_curr_ver + plus_curr_01
					k_start_01 = k_start_01 + 1

			#find test point, axis 10
			total_arc_len_curr_hor = total_arc_len_curr_hor_start
			less = 1 
			while less:
				total_arc_len_curr_hor_start = total_arc_len_curr_hor
				total_arc_len_curr_hor = total_arc_len_curr_hor + arc_len_10[N_arc_length-k_end_10-2]
				diff_curr = total_arc_len_curr_hor - arc_len_one_step_hor
				if comparison(diff_curr, 0.0):
					arc_len_end_10 = total_arc_len_curr_hor
					points_x_10[Nx-2] = points_x_to_arc_len_calc_10[N_arc_length-k_end_10-2] 
					points_y_10[Nx-2] = points_y_to_arc_len_calc_10[N_arc_length-k_end_10-2]
					less = 0
				elif ( diff_curr > 0.0 ):
					t_curr_new = t_curr_hor - plus_curr_10*(1.0-diff_curr/arc_len_10[N_arc_length-k_end_10-2])
					arc_len_end_10 = total_arc_len_curr_hor - diff_curr
					points_x_10[Nx-2] = X10(t_curr_new)
					points_y_10[Nx-2] = Y10(t_curr_new) 
					less = 0
				else:
					if (k_end_10>(N_arc_length-3)):
						print("ERROR bad angle at point B. Try to use bigger N_arc_length. Current is %s. Grid" %(N_arc_length) + grid_path)
						return 0
					t_curr_hor = t_curr_hor - plus_curr_10
					k_end_10 = k_end_10 + 1

			bad_point = test_normalized_vector_product(points_x_01[1]-points_x_01[0], points_y_01[1]-points_y_01[0], points_x_10[Nx-2]-points_x_10[Nx-1], points_y_10[Nx-2]-points_y_10[Nx-1])
			arc_len_one_step_ver = arc_len_one_step_ver + arc_len_delta_ver
			arc_len_one_step_hor = arc_len_one_step_hor + arc_len_delta_hor
			if First_try_B and bad_point:
				First_try_B = 0
			if (arc_len_one_step_ver>total_arc_len01):
				print("ERROR bad angle at point B. Try to use bigger N_arc_length. Current is %s. Grid" %(N_arc_length) + grid_path)
				return 0	
			if (arc_len_one_step_hor>total_arc_len10):
				print("ERROR bad angle at point B. Try to use bigger N_arc_length. Current is %s. Grid" %(N_arc_length) + grid_path)
				return 0

		total_arc_len_curr_new_start_01 = total_arc_len_curr_ver_start
		t_new_start_01 = t_curr_ver

		#Point D, 00 and 11
		arc_len_one_step_ver = (total_arc_len00-arc_len_start_00)/((float)(Ny-2))
		arc_len_delta_ver = coef_arc_len_critical*arc_len_one_step_ver
		t_curr_ver = t_end00
		total_arc_len_curr_ver_start = 0.0
		k_end_00 = 0

		arc_len_one_step_hor = total_arc_len11/((float)(Nx-1))
		arc_len_delta_hor = coef_arc_len_critical*arc_len_one_step_hor
		t_curr_hor = t_start11
		total_arc_len_curr_hor_start = 0.0
		k_start_11 = 0

		First_try_D = 1
		bad_point = 1
		while bad_point:
			#find test point, axis 00
			total_arc_len_curr_ver = total_arc_len_curr_ver_start
			less = 1 
			while less:
				total_arc_len_curr_ver_start = total_arc_len_curr_ver
				total_arc_len_curr_ver = total_arc_len_curr_ver + arc_len_00[N_arc_length-k_end_00-2]
				diff_curr = total_arc_len_curr_ver - arc_len_one_step_ver
				if comparison(diff_curr, 0.0):
					arc_len_end_00 = total_arc_len_curr_ver
					points_x_00[Ny-2] = points_x_to_arc_len_calc_00[N_arc_length-k_end_00-2] 
					points_y_00[Ny-2] = points_y_to_arc_len_calc_00[N_arc_length-k_end_00-2]
					less = 0
				elif ( diff_curr > 0.0 ):
					t_curr_new = t_curr_ver - plus_curr_00*(1.0-diff_curr/arc_len_00[N_arc_length-k_end_00-2])
					arc_len_end_00 = total_arc_len_curr_ver - diff_curr
					points_x_00[Ny-2] = X00(t_curr_new)
					points_y_00[Ny-2] = Y00(t_curr_new) 
					less = 0
				else:
					if (k_end_00>(N_arc_length-3)):
						print("ERROR bad angle at point D. Try to use bigger N_arc_length. Current is %s. Grid" %(N_arc_length) + grid_path)
						return 0
					t_curr_ver = t_curr_ver - plus_curr_00
					k_end_00 = k_end_00 + 1

			#find test point, axis 11
			total_arc_len_curr_hor = total_arc_len_curr_hor_start
			less = 1 
			while less:
				total_arc_len_curr_hor_start = total_arc_len_curr_hor
				total_arc_len_curr_hor = total_arc_len_curr_hor + arc_len_11[k_start_11]
				diff_curr = total_arc_len_curr_hor - arc_len_one_step_hor
				if comparison(diff_curr, 0.0):
					arc_len_start_11 = total_arc_len_curr_hor
					points_x_11[1] = points_x_to_arc_len_calc_11[k_start_11+1] 
					points_y_11[1] = points_y_to_arc_len_calc_11[k_start_11+1]
					less = 0
				elif ( diff_curr > 0.0 ):
					t_curr_new = t_curr_hor + plus_curr_11*(1.0-diff_curr/arc_len_11[k_start_11])
					arc_len_start_11 = total_arc_len_curr_hor - diff_curr
					points_x_11[1] = X11(t_curr_new)
					points_y_11[1] = Y11(t_curr_new) 
					less = 0
				else:
					if (k_start_11>(N_arc_length-3)):
						print("ERROR bad angle at point D. Try to use bigger N_arc_length. Current is %s. Grid" %(N_arc_length) + grid_path)
						return 0
					t_curr_hor = t_curr_hor + plus_curr_11
					k_start_11 = k_start_11 + 1

			bad_point = test_normalized_vector_product(points_x_00[Ny-2]-points_x_00[Ny-1], points_y_00[Ny-2]-points_y_00[Ny-1], points_x_11[1]-points_x_11[0], points_y_11[1]-points_y_11[0])
			arc_len_one_step_ver = arc_len_one_step_ver + arc_len_delta_ver
			arc_len_one_step_hor = arc_len_one_step_hor + arc_len_delta_hor
			if First_try_D and bad_point:
				First_try_D = 0
			if (arc_len_one_step_ver>total_arc_len00):
				print("ERROR bad angle at point D. Try to use bigger N_arc_length. Current is %s. Grid" %(N_arc_length) + grid_path)
				return 0	
			if (arc_len_one_step_hor>total_arc_len11):
				print("ERROR bad angle at point D. Try to use bigger N_arc_length. Current is %s. Grid" %(N_arc_length) + grid_path)
				return 0

		total_arc_len_curr_new_start_11 = total_arc_len_curr_hor_start
		t_new_start_11 = t_curr_hor

		#Point C, 01 and 11
		arc_len_one_step_ver = (total_arc_len01-arc_len_start_01)/((float)(Ny-2))
		arc_len_delta_ver = coef_arc_len_critical*arc_len_one_step_ver
		t_curr_ver = t_end01
		total_arc_len_curr_ver_start = 0.0
		k_end_01 = 0

		arc_len_one_step_hor = (total_arc_len11-arc_len_start_11)/((float)(Nx-2))
		arc_len_delta_hor = coef_arc_len_critical*arc_len_one_step_hor
		t_curr_hor = t_end11
		total_arc_len_curr_hor_start = 0.0
		k_end_11 = 0

		First_try_C = 1
		bad_point = 1
		while bad_point:
			#find test point, axis 01
			total_arc_len_curr_ver = total_arc_len_curr_ver_start
			less = 1 
			while less:
				total_arc_len_curr_ver_start = total_arc_len_curr_ver
				total_arc_len_curr_ver = total_arc_len_curr_ver + arc_len_01[N_arc_length-k_end_01-2]
				diff_curr = total_arc_len_curr_ver - arc_len_one_step_ver
				if comparison(diff_curr, 0.0):
					arc_len_end_01 = total_arc_len_curr_ver
					points_x_01[Ny-2] = points_x_to_arc_len_calc_01[N_arc_length-k_end_01-2] 
					points_y_01[Ny-2] = points_y_to_arc_len_calc_01[N_arc_length-k_end_01-2]
					less = 0
				elif ( diff_curr > 0.0 ):
					t_curr_new = t_curr_ver - plus_curr_01*(1.0-diff_curr/arc_len_01[N_arc_length-k_end_01-2])
					arc_len_end_01 = total_arc_len_curr_ver - diff_curr
					points_x_01[Ny-2] = X01(t_curr_new)
					points_y_01[Ny-2] = Y01(t_curr_new) 
					less = 0
				else:
					if (k_end_01>(N_arc_length-3)):
						print("ERROR bad angle at point C. Try to use bigger N_arc_length. Current is %s. Grid" %(N_arc_length) + grid_path)
						return 0
					t_curr_ver = t_curr_ver - plus_curr_01
					k_end_01 = k_end_01 + 1

			#find test point, axis 11
			total_arc_len_curr_hor = total_arc_len_curr_hor_start
			less = 1 
			while less:
				total_arc_len_curr_hor_start = total_arc_len_curr_hor
				total_arc_len_curr_hor = total_arc_len_curr_hor + arc_len_11[N_arc_length-k_end_11-2]
				diff_curr = total_arc_len_curr_hor - arc_len_one_step_hor
				if comparison(diff_curr, 0.0):
					arc_len_end_11 = total_arc_len_curr_hor
					points_x_11[Nx-2] = points_x_to_arc_len_calc_11[N_arc_length-k_end_11-2] 
					points_y_11[Nx-2] = points_y_to_arc_len_calc_11[N_arc_length-k_end_11-2]
					less = 0
				elif ( diff_curr > 0.0 ):
					t_curr_new = t_curr_hor - plus_curr_11*(1.0-diff_curr/arc_len_11[N_arc_length-k_end_11-2])
					arc_len_end_11 = total_arc_len_curr_hor - diff_curr
					points_x_11[Nx-2] = X11(t_curr_new)
					points_y_11[Nx-2] = Y11(t_curr_new) 
					less = 0
				else:
					if (k_end_11>(N_arc_length-3)):
						print("ERROR bad angle at point C. Try to use bigger N_arc_length. Current is %s. Grid" %(N_arc_length) + grid_path)
						return 0
					t_curr_hor = t_curr_hor - plus_curr_11
					k_end_11 = k_end_11 + 1

			bad_point = test_normalized_vector_product(points_x_11[Nx-2]-points_x_11[Nx-1], points_y_11[Nx-2]-points_y_11[Nx-1], points_x_01[Ny-2]-points_x_01[Ny-1], points_y_01[Ny-2]-points_y_01[Ny-1])
			arc_len_one_step_ver = arc_len_one_step_ver + arc_len_delta_ver
			arc_len_one_step_hor = arc_len_one_step_hor + arc_len_delta_hor	
			if First_try_C and bad_point:
				First_try_C = 0
			if (arc_len_one_step_ver>total_arc_len01):
				print("ERROR bad angle at point C. Try to use bigger N_arc_length. Current is %s. Grid" %(N_arc_length) + grid_path)
				return 0	
			if (arc_len_one_step_hor>total_arc_len11):
				print("ERROR bad angle at point C. Try to use bigger N_arc_length. Current is %s. Grid" %(N_arc_length) + grid_path)
				return 0

		if (First_try_A and ((not First_try_B) or (not First_try_D))):
			#Point A, 00 and 10
			arc_len_one_step_ver = (total_arc_len00-arc_len_end_00)/((float)(Ny-2))
			arc_len_delta_ver = coef_arc_len_critical*arc_len_one_step_ver
			t_curr_ver = t_start00
			total_arc_len_curr_ver_start = 0.0
			k_start_00 = 0

			arc_len_one_step_hor = (total_arc_len10-arc_len_end_10)/((float)(Nx-2))
			arc_len_delta_hor = coef_arc_len_critical*arc_len_one_step_hor
			t_curr_hor = t_start10
			total_arc_len_curr_hor_start = 0.0
			k_start_10 = 0

			bad_point = 1
			while bad_point:
				#find test point, axis 00
				total_arc_len_curr_ver = total_arc_len_curr_ver_start
				less = 1 
				while less:
					total_arc_len_curr_ver_start = total_arc_len_curr_ver
					total_arc_len_curr_ver = total_arc_len_curr_ver + arc_len_00[k_start_00]
					diff_curr = total_arc_len_curr_ver - arc_len_one_step_ver
					if comparison(diff_curr, 0.0):
						arc_len_start_00 = total_arc_len_curr_ver
						total_arc_len_curr_new_start_00 = total_arc_len_curr_ver
						points_x_00[1] = points_x_to_arc_len_calc_00[k_start_00+1] 
						points_y_00[1] = points_y_to_arc_len_calc_00[k_start_00+1]
						less = 0
					elif ( diff_curr > 0.0 ):
						t_curr_new = t_curr_ver + plus_curr_00*(1.0-diff_curr/arc_len_00[k_start_00])
						arc_len_start_00 = total_arc_len_curr_ver - diff_curr
						points_x_00[1] = X00(t_curr_new)
						points_y_00[1] = Y00(t_curr_new) 
						less = 0
					else:
						t_curr_ver = t_curr_ver + plus_curr_00
						k_start_00 = k_start_00 + 1
 

				#find test point, axis 10
				total_arc_len_curr_hor = total_arc_len_curr_hor_start
				less = 1 
				while less:
					total_arc_len_curr_hor_start = total_arc_len_curr_hor
					total_arc_len_curr_hor = total_arc_len_curr_hor + arc_len_10[k_start_10]
					diff_curr = total_arc_len_curr_hor - arc_len_one_step_hor
					if comparison(diff_curr, 0.0):
						arc_len_start_10 = total_arc_len_curr_hor
						points_x_10[1] = points_x_to_arc_len_calc_10[k_start_10+1] 
						points_y_10[1] = points_y_to_arc_len_calc_10[k_start_10+1]
						less = 0
					elif ( diff_curr > 0.0 ):
						t_curr_new = t_curr_hor + plus_curr_10*(1.0-diff_curr/arc_len_10[k_start_10])
						arc_len_start_10 = total_arc_len_curr_hor - diff_curr
						points_x_10[1] = X10(t_curr_new)
						points_y_10[1] = Y10(t_curr_new) 
						less = 0
					else:
						t_curr_hor = t_curr_hor + plus_curr_10
						k_start_10 = k_start_10 + 1

				bad_point = test_normalized_vector_product(points_x_10[1]-points_x_10[0], points_y_10[1]-points_y_10[0], points_x_00[1]-points_x_00[0], points_y_00[1]-points_y_00[0])
				arc_len_one_step_ver = arc_len_one_step_ver + arc_len_delta_ver
				arc_len_one_step_hor = arc_len_one_step_hor + arc_len_delta_hor

			total_arc_len_curr_new_start_00 = total_arc_len_curr_ver_start
			total_arc_len_curr_new_start_10 = total_arc_len_curr_hor_start
			t_new_start_00 = t_curr_ver
			t_new_start_10 = t_curr_hor

		if (First_try_B and (not First_try_C)):
			#Point B, 01 and 10
			arc_len_one_step_ver = (total_arc_len01-arc_len_end_01)/((float)(Ny-1))
			arc_len_delta_ver = coef_arc_len_critical*arc_len_one_step_ver
			t_curr_ver = t_start01
			total_arc_len_curr_ver_start = 0.0
			k_start_01 = 0

			arc_len_one_step_hor = (total_arc_len10-arc_len_start_10)/((float)(Nx-2))
			arc_len_delta_hor = coef_arc_len_critical*arc_len_one_step_hor
			t_curr_hor = t_end10
			total_arc_len_curr_hor_start = 0.0
			k_end_10 = 0

			bad_point = 1
			while bad_point:
				#find test point, axis 01
				total_arc_len_curr_ver = total_arc_len_curr_ver_start
				less = 1 
				while less:
					total_arc_len_curr_ver_start = total_arc_len_curr_ver
					total_arc_len_curr_ver = total_arc_len_curr_ver + arc_len_01[k_start_01]
					diff_curr = total_arc_len_curr_ver - arc_len_one_step_ver
					if comparison(diff_curr, 0.0):
						arc_len_start_01 = total_arc_len_curr_ver
						points_x_01[1] = points_x_to_arc_len_calc_01[k_start_01+1] 
						points_y_01[1] = points_y_to_arc_len_calc_01[k_start_01+1]
						less = 0
					elif ( diff_curr > 0.0 ):
						t_curr_new = t_curr_ver + plus_curr_01*(1.0-diff_curr/arc_len_01[k_start_01])
						arc_len_start_01 = total_arc_len_curr_ver - diff_curr
						points_x_01[1] = X01(t_curr_new)
						points_y_01[1] = Y01(t_curr_new) 
						less = 0
					else:
						t_curr_ver = t_curr_ver + plus_curr_01
						k_start_01 = k_start_01 + 1

				#find test point, axis 10
				total_arc_len_curr_hor = total_arc_len_curr_hor_start
				less = 1 
				while less:
					total_arc_len_curr_hor_start = total_arc_len_curr_hor
					total_arc_len_curr_hor = total_arc_len_curr_hor + arc_len_10[N_arc_length-k_end_10-2]
					diff_curr = total_arc_len_curr_hor - arc_len_one_step_hor
					if comparison(diff_curr, 0.0):
						arc_len_end_10 = total_arc_len_curr_hor
						points_x_10[Nx-2] = points_x_to_arc_len_calc_10[N_arc_length-k_end_10-2] 
						points_y_10[Nx-2] = points_y_to_arc_len_calc_10[N_arc_length-k_end_10-2]
						less = 0
					elif ( diff_curr > 0.0 ):
						t_curr_new = t_curr_hor - plus_curr_10*(1.0-diff_curr/arc_len_10[N_arc_length-k_end_10-2])
						arc_len_end_10 = total_arc_len_curr_hor - diff_curr
						points_x_10[Nx-2] = X10(t_curr_new)
						points_y_10[Nx-2] = Y10(t_curr_new) 
						less = 0
					else:
						t_curr_hor = t_curr_hor - plus_curr_10
						k_end_10 = k_end_10 + 1

				bad_point = test_normalized_vector_product(points_x_01[1]-points_x_01[0], points_y_01[1]-points_y_01[0], points_x_10[Nx-2]-points_x_10[Nx-1], points_y_10[Nx-2]-points_y_10[Nx-1])
				arc_len_one_step_ver = arc_len_one_step_ver + arc_len_delta_ver
				arc_len_one_step_hor = arc_len_one_step_hor + arc_len_delta_hor

			total_arc_len_curr_new_start_01 = total_arc_len_curr_ver_start
			t_new_start_01 = t_curr_ver

		if (First_try_D and (not First_try_C)):		
			#Point D, 00 and 11
			arc_len_one_step_ver = (total_arc_len00-arc_len_start_00)/((float)(Ny-2))
			arc_len_delta_ver = coef_arc_len_critical*arc_len_one_step_ver
			t_curr_ver = t_end00
			total_arc_len_curr_ver_start = 0.0
			k_end_00 = 0

			arc_len_one_step_hor = (total_arc_len11-arc_len_end_11)/((float)(Nx-1))
			arc_len_delta_hor = coef_arc_len_critical*arc_len_one_step_hor
			t_curr_hor = t_start11
			total_arc_len_curr_hor_start = 0.0
			k_start_11 = 0

			bad_point = 1
			while bad_point:
				#find test point, axis 00
				total_arc_len_curr_ver = total_arc_len_curr_ver_start
				less = 1 
				while less:
					total_arc_len_curr_ver_start = total_arc_len_curr_ver
					total_arc_len_curr_ver = total_arc_len_curr_ver + arc_len_00[N_arc_length-k_end_00-2]
					diff_curr = total_arc_len_curr_ver - arc_len_one_step_ver
					if comparison(diff_curr, 0.0):
						arc_len_end_00 = total_arc_len_curr_ver
						points_x_00[Ny-2] = points_x_to_arc_len_calc_00[N_arc_length-k_end_00-2] 
						points_y_00[Ny-2] = points_y_to_arc_len_calc_00[N_arc_length-k_end_00-2]
						less = 0
					elif ( diff_curr > 0.0 ):
						t_curr_new = t_curr_ver - plus_curr_00*(1.0-diff_curr/arc_len_00[N_arc_length-k_end_00-2])
						arc_len_end_00 = total_arc_len_curr_ver - diff_curr
						points_x_00[Ny-2] = X00(t_curr_new)
						points_y_00[Ny-2] = Y00(t_curr_new) 
						less = 0
					else:
						t_curr_ver = t_curr_ver - plus_curr_00
						k_end_00 = k_end_00 + 1

				#find test point, axis 11
				total_arc_len_curr_hor = total_arc_len_curr_hor_start
				less = 1 
				while less:
					total_arc_len_curr_hor_start = total_arc_len_curr_hor
					total_arc_len_curr_hor = total_arc_len_curr_hor + arc_len_11[k_start_11]
					diff_curr = total_arc_len_curr_hor - arc_len_one_step_hor
					if comparison(diff_curr, 0.0):
						arc_len_start_11 = total_arc_len_curr_hor
						points_x_11[1] = points_x_to_arc_len_calc_11[k_start_11+1] 
						points_y_11[1] = points_y_to_arc_len_calc_11[k_start_11+1]
						less = 0
					elif ( diff_curr > 0.0 ):
						t_curr_new = t_curr_hor + plus_curr_11*(1.0-diff_curr/arc_len_11[k_start_11])
						arc_len_start_11 = total_arc_len_curr_hor - diff_curr
						points_x_11[1] = X11(t_curr_new)
						points_y_11[1] = Y11(t_curr_new) 
						less = 0
					else:
						t_curr_hor = t_curr_hor + plus_curr_11
						k_start_11 = k_start_11 + 1

				bad_point = test_normalized_vector_product(points_x_00[Ny-2]-points_x_00[Ny-1], points_y_00[Ny-2]-points_y_00[Ny-1], points_x_11[1]-points_x_11[0], points_y_11[1]-points_y_11[0])
				arc_len_one_step_ver = arc_len_one_step_ver + arc_len_delta_ver
				arc_len_one_step_hor = arc_len_one_step_hor + arc_len_delta_hor

			total_arc_len_curr_new_start_11 = total_arc_len_curr_hor_start
			t_new_start_11 = t_curr_hor

	if (Ny>2):
		if (grid_correction and (Nx>2) and (Ny>4)):
			arc_len_one_step = (total_arc_len00-arc_len_start_00-arc_len_end_00)/((float)(Ny-3))
			k = k_start_00
			k_max = N_arc_length - k_end_00 - 1
			t_curr = t_new_start_00
			total_arc_len_new_start = total_arc_len_curr_new_start_00
		else:
			arc_len_one_step = total_arc_len00/((float)(Ny-1))
			k = 0
			k_max = N_arc_length-1
			t_curr = t_start00
			total_arc_len_new_start = 0.0

		if ((not grid_correction) or (Ny>4)):
			index_start = 1 + grid_correction
			max_index = Ny-1 - grid_correction

			for index in range(index_start,max_index):
				not_find = 1 
				while not_find:
					k_start = k
					total_arc_len_curr = total_arc_len_new_start
					for k in range(k_start,k_max):
						total_arc_len_new_start = total_arc_len_curr
						total_arc_len_curr = total_arc_len_curr + arc_len_00[k]
						diff_curr = total_arc_len_curr-index*arc_len_one_step
						if comparison(diff_curr, 0.0):
							points_x_00[index] = points_x_to_arc_len_calc_00[k+1] 
							points_y_00[index] = points_y_to_arc_len_calc_00[k+1]
							not_find = 0 
							break
						elif ( diff_curr > 0.0 ):
							t_curr_new = t_curr + plus_curr_00*(1.0-diff_curr/arc_len_00[k])
							points_x_00[index] = X00(t_curr_new)
							points_y_00[index] = Y00(t_curr_new) 
							not_find = 0
							break
						t_curr = t_curr + plus_curr_00

		if (grid_correction and (Nx>2) and (Ny>4)):
			arc_len_one_step = (total_arc_len01-arc_len_start_01-arc_len_end_01)/((float)(Ny-3))
			k = k_start_01
			k_max = N_arc_length - k_end_01 - 1
			t_curr = t_new_start_01
			total_arc_len_new_start = total_arc_len_curr_new_start_01
		else:
			arc_len_one_step = total_arc_len01/((float)(Ny-1))
			k = 0
			k_max = N_arc_length-1
			t_curr = t_start01
			total_arc_len_new_start = 0.0

		if ((not grid_correction) or (Ny>4)):
			index_start = 1 + grid_correction
			max_index = Ny-1 - grid_correction

			for index in range(index_start,max_index):
				not_find = 1 
				while not_find:
					k_start = k
					total_arc_len_curr = total_arc_len_new_start
					for k in range(k_start,k_max):
						total_arc_len_new_start = total_arc_len_curr
						total_arc_len_curr = total_arc_len_curr + arc_len_01[k]
						diff_curr = total_arc_len_curr-index*arc_len_one_step
						if comparison(diff_curr, 0.0):
							points_x_01[index] = points_x_to_arc_len_calc_01[k+1] 
							points_y_01[index] = points_y_to_arc_len_calc_01[k+1]
							not_find = 0 
							break
						elif ( diff_curr > 0.0 ):
							t_curr_new = t_curr + plus_curr_01*(1.0-diff_curr/arc_len_01[k])
							points_x_01[index] = X01(t_curr_new)
							points_y_01[index] = Y01(t_curr_new) 
							not_find = 0
							break
						t_curr = t_curr + plus_curr_01

	if (Nx>2):
		if (grid_correction and (Ny>2) and (Nx>4)):
			arc_len_one_step = (total_arc_len10-arc_len_start_10-arc_len_end_10)/((float)(Nx-3))
			k = k_start_10
			k_max = N_arc_length - k_end_10 - 1
			t_curr = t_new_start_10
			total_arc_len_new_start = total_arc_len_curr_new_start_10
		else:
			arc_len_one_step = total_arc_len10/((float)(Nx-1))
			k = 0
			k_max = N_arc_length-1
			t_curr = t_start10
			total_arc_len_new_start = 0.0

		if ((not grid_correction) or (Nx>4)):
			index_start = 1 + grid_correction
			max_index = Nx-1 - grid_correction

			for index in range(index_start,max_index):
				not_find = 1 
				while not_find:
					k_start = k
					total_arc_len_curr = total_arc_len_new_start
					for k in range(k_start,k_max):
						total_arc_len_new_start = total_arc_len_curr
						total_arc_len_curr = total_arc_len_curr + arc_len_10[k]
						diff_curr = total_arc_len_curr-index*arc_len_one_step
						if comparison(diff_curr, 0.0):
							points_x_10[index] = points_x_to_arc_len_calc_10[k+1] 
							points_y_10[index] = points_y_to_arc_len_calc_10[k+1]
							not_find = 0 
							break
						elif ( diff_curr > 0.0 ):
							t_curr_new = t_curr + plus_curr_10*(1.0-diff_curr/arc_len_10[k])
							points_x_10[index] = X10(t_curr_new)
							points_y_10[index] = Y10(t_curr_new) 
							not_find = 0
							break
						t_curr = t_curr + plus_curr_10

		if (grid_correction and (Ny>2) and (Nx>4)):
			arc_len_one_step = (total_arc_len11-arc_len_start_11-arc_len_end_11)/((float)(Nx-3))
			k = k_start_11
			k_max = N_arc_length - k_end_11 - 1
			t_curr = t_new_start_11
			total_arc_len_new_start = total_arc_len_curr_new_start_11
		else:
			arc_len_one_step = total_arc_len11/((float)(Nx-1))
			k = 0
			k_max = N_arc_length-1
			t_curr = t_start11
			total_arc_len_new_start = 0.0

		if ((not grid_correction) or (Nx>4)):
			index_start = 1 + grid_correction
			max_index = Nx-1 - grid_correction
	
			for index in range(index_start,max_index):
				not_find = 1 
				while not_find:
					k_start = k
					total_arc_len_curr = total_arc_len_new_start
					for k in range(k_start,k_max):
						total_arc_len_new_start = total_arc_len_curr
						total_arc_len_curr = total_arc_len_curr + arc_len_11[k]
						diff_curr = total_arc_len_curr-index*arc_len_one_step
						if comparison(diff_curr, 0.0):
							points_x_11[index] = points_x_to_arc_len_calc_11[k+1] 
							points_y_11[index] = points_y_to_arc_len_calc_11[k+1]
							not_find = 0 
							break
						elif ( diff_curr > 0.0 ):
							t_curr_new = t_curr + plus_curr_11*(1.0-diff_curr/arc_len_11[k])
							points_x_11[index] = X11(t_curr_new)
							points_y_11[index] = Y11(t_curr_new) 
							not_find = 0
							break
						t_curr = t_curr + plus_curr_11
	return 1

### GRID GENERATION

def generate_quadrangle(grid_path, Nx, Ny, xA, yA, xB, yB, xC, yC, xD, yD):
	points_x = np.zeros((Nx,Ny), dtype = 'float')
	points_y = np.zeros((Nx,Ny), dtype = 'float')

	points_x[0,0] = xA
	points_y[0,0] = yA
	points_x[Nx-1,0] = xB
	points_y[Nx-1,0] = yB
	points_x[Nx-1,Ny-1] = xC
	points_y[Nx-1,Ny-1] = yC
	points_x[0,Ny-1] = xD
	points_y[0,Ny-1] = yD

	dxAB = (xB - xA) / (Nx - 1)
	dyAB = (yB - yA) / (Nx - 1)
	dxDC = (xC - xD) / (Nx - 1)
	dyDC = (yC - yD) / (Nx - 1)

	#trapezoid grid
	for i in range(Nx):
		if (i == 0):
			xAlocal = xA
			yAlocal = yA
			xDlocal = xD
			yDlocal = yD
		else:
			xAlocal = xA + dxAB * i
			yAlocal = yA + dyAB * i
			xDlocal = xD + dxDC * i
			yDlocal = yD + dyDC * i
		dxADlocal = (xDlocal - xAlocal) / (Ny - 1)
		dyADlocal = (yDlocal - yAlocal) / (Ny - 1)
		for j in range(Ny):
			if not ( (i==Nx-1)and( (j==0)or(j==Ny-1) ) ): #not B, not C. A and D is already without calculation
				if (j == 0):
					points_x[i,j] = xAlocal
					points_y[i,j] = yAlocal
				else:
					points_x[i,j] = xAlocal + dxADlocal * j
					points_y[i,j] = yAlocal + dyADlocal * j

	save_to_vtk(grid_path, points_x, points_y, Nx, Ny)
	return 1	

#	points_x = np.zeros((Nx,Ny), dtype = 'float')
#	points_y = np.zeros((Nx,Ny), dtype = 'float')
def trapezoid_grid(Nx, Ny, xA, yA, xB, yB, xC, yC, xD, yD, points_x, points_y, if_corners):
	if if_corners:
		points_x[0,0] = xA
		points_y[0,0] = yA
		points_x[Nx-1,0] = xB
		points_y[Nx-1,0] = yB
		points_x[Nx-1,Ny-1] = xC
		points_y[Nx-1,Ny-1] = yC
		points_x[0,Ny-1] = xD
		points_y[0,Ny-1] = yD
	if ((Nx == 2) and (Ny == 2)):
		return 1

#	if (Nx == 2):	
#		divide_line(Ny, xA, yA, xD, yD, points_x[0], points_y[0])
#		divide_line(Ny, xB, yB, xC, yC, points_x[1], points_y[1])
#		return 1
#	if (Ny == 2):
#		temp = np.zeros((2,Nx), dtype = 'float')
#		divide_line(Nx, xA, yA, xB, yB, temp[0], temp[1])
#		for i in range(Nx):
#			points_x[i,0] = temp[0,i]
#			points_y[i,0] = temp[1,i]
#		divide_line(Nx, xD, yD, xC, yC, temp[0], temp[1])
#		for i in range(Nx):
#			points_x[i,1] = temp[0,i]
#			points_y[i,1] = temp[1,i]			
#		del(temp)
#		return 1
		
	div_Nx = 1.0 / ((float)(Nx - 1))
	div_Ny = 1.0 / ((float)(Ny - 1))

	dxAB = (xB - xA) * div_Nx
	dyAB = (yB - yA) * div_Nx
	dxDC = (xC - xD) * div_Nx
	dyDC = (yC - yD) * div_Nx

	#trapezoid grid
	for i in range(Nx):
		if (i == 0):
			xAlocal = xA
			yAlocal = yA
			xDlocal = xD
			yDlocal = yD
		else:
			xAlocal = xA + dxAB * i
			yAlocal = yA + dyAB * i
			xDlocal = xD + dxDC * i
			yDlocal = yD + dyDC * i
		dxADlocal = (xDlocal - xAlocal) * div_Ny
		dyADlocal = (yDlocal - yAlocal) * div_Ny
		for j in range(Ny):
			if not ( (i==Nx-1)and( (j==0)or(j==Ny-1) ) ): #not B, not C. A and D is already without calculation
				if (j == 0):
					points_x[i,j] = xAlocal
					points_y[i,j] = yAlocal
				else:
					points_x[i,j] = xAlocal + dxADlocal * j
					points_y[i,j] = yAlocal + dyADlocal * j

	return 1	

####

def generate_curved_grid_from_boundary_points_in_array_from_axis_arcs(Nx, Ny, arcs_axis0_x, arcs_axis0_y, arcs_axis1_x, arcs_axis1_y, grid_correction, array_x, array_y, if_test_the_grid, if_test_corners):
	return generate_curved_grid_from_boundary_points_separate_setting_in_array(Nx, Ny, arcs_axis0_x[0], arcs_axis0_y[0], arcs_axis0_x[1], arcs_axis0_y[1], arcs_axis1_x[0], arcs_axis1_y[0], arcs_axis1_x[1], arcs_axis1_y[1], grid_correction, array_x, array_y, if_test_the_grid, if_test_corners)

def generate_curved_grid_from_boundary_points_in_array(Nx, Ny, points_x, points_y, grid_correction, array_x, array_y, if_test_the_grid, if_test_corners):
	return generate_curved_grid_from_boundary_points_separate_setting_in_array(Nx, Ny, points_x[0], points_y[0], points_x[1], points_y[1], points_x[2], points_y[2], points_x[3], points_y[3], grid_correction, array_x, array_y, if_test_the_grid, if_test_corners)

def generate_curved_grid_from_boundary_points(grid_path, Nx, Ny, points_x, points_y, grid_correction, if_test_the_grid, if_test_corners):
	return generate_curved_grid_from_boundary_points_separate_setting(grid_path, Nx, Ny, points_x[0], points_y[0], points_x[1], points_y[1], points_x[2], points_y[2], points_x[3], points_y[3], grid_correction, if_test_the_grid, if_test_corners)

def generate_curved_grid_from_boundary_points_separate_setting(grid_path, Nx, Ny, points_x_00, points_y_00, points_x_01, points_y_01, points_x_10, points_y_10, points_x_11, points_y_11, grid_correction, if_test_the_grid, if_test_corners):
	points_x = np.zeros((Nx,Ny), dtype = 'float')
	points_y = np.zeros((Nx,Ny), dtype = 'float')

	result = generate_curved_grid_from_boundary_points_separate_setting_in_array(Nx, Ny, points_x_00, points_y_00, points_x_01, points_y_01, points_x_10, points_y_10, points_x_11, points_y_11, grid_correction, points_x, points_y, if_test_the_grid, if_test_corners)

	save_to_vtk(grid_path, points_x, points_y, Nx, Ny)
	del(points_x)
	del(points_y)
	
	return result

def generate_curved_grid_from_boundary_points_separate_setting_in_array(Nx, Ny, points_x_00, points_y_00, points_x_01, points_y_01, points_x_10, points_y_10, points_x_11, points_y_11, grid_correction, points_x, points_y, if_test_the_grid, if_test_corners):

#	if_save_iterations = get_if_save_iterations()
#	if_save_iterations_before_correction = get_if_save_iterations_before_correction()
#	print_grids = get_print_grids()

	if (if_test_corners == 1):
		if not comparison(points_x_00[0],points_x_10[0]):
			sys.exit("ERROR DATA point A, lines 00 10, component X (%s,%s) IN GRID " %(points_x_00[0],points_x_10[0]))
		if not comparison(points_x_01[0],points_x_10[Nx-1]):
			sys.exit("ERROR DATA point B, lines 01 10, component X (%s,%s) IN GRID " %(points_x_01[0],points_x_10[Nx-1]))
		if not comparison(points_x_01[Ny-1],points_x_11[Nx-1]):
			sys.exit("ERROR DATA point C, lines 01 11, component X (%s,%s) IN GRID " %(points_x_01[Ny-1],points_x_11[Nx-1]))
		if not comparison(points_x_00[Ny-1],points_x_11[0]):
			sys.exit("ERROR DATA point D, lines 00 11, component X (%s,%s) IN GRID " %(points_x_00[Ny-1],points_x_11[0]))
		if not comparison(points_y_00[0],points_y_10[0]):
			sys.exit("ERROR DATA point A, lines 00 10, component Y (%s,%s) IN GRID " %(points_y_00[0],points_y_10[0]))
		if not comparison(points_y_01[0],points_y_10[Nx-1]):
			sys.exit("ERROR DATA point B, lines 01 10, component Y (%s,%s) IN GRID " %(points_y_01[0],points_y_10[Nx-1]))
		if not comparison(points_y_01[Ny-1],points_y_11[Nx-1]):
			sys.exit("ERROR DATA point C, lines 01 11, component Y (%s,%s) IN GRID " %(points_y_01[Ny-1],points_y_11[Nx-1]))
		if not comparison(points_y_00[Ny-1],points_y_11[0]):
			sys.exit("ERROR DATA point D, lines 00 11, component Y (%s,%s) IN GRID " %(points_y_00[Ny-1],points_y_11[0]))

	points_x[0,0] = points_x_00[0]
	points_y[0,0] = points_y_00[0]
	points_x[Nx-1,0] = points_x_01[0]
	points_y[Nx-1,0] = points_y_01[0]
	points_x[Nx-1,Ny-1] = points_x_11[Nx-1]
	points_y[Nx-1,Ny-1] = points_y_11[Nx-1]
	points_x[0,Ny-1] = points_x_11[0]
	points_y[0,Ny-1] = points_y_11[0]

	result = 1

	if ((Nx>2) and (Ny==2)):
		for index in range(1,Nx-1):
			points_x[index,0] = points_x_10[index]
			points_y[index,0] = points_y_10[index]
			points_x[index,1] = points_x_11[index]
			points_y[index,1] = points_y_11[index]

	elif ((Nx==2) and (Ny>2)):
		for index in range(1,Ny-1):
			points_x[0,index] = points_x_00[index]
			points_y[0,index] = points_y_00[index]
			points_x[1,index] = points_x_01[index]
			points_y[1,index] = points_y_01[index]

	else:
		conv00x = np.zeros((Ny-2), dtype = 'float')
		conv00y = np.zeros((Ny-2), dtype = 'float')
		conv01x = np.zeros((Ny-2), dtype = 'float')
		conv01y = np.zeros((Ny-2), dtype = 'float')
		conv10x = np.zeros((Nx-2), dtype = 'float')
		conv10y = np.zeros((Nx-2), dtype = 'float')
		conv11x = np.zeros((Nx-2), dtype = 'float')
		conv11y = np.zeros((Nx-2), dtype = 'float')

		xA = points_x[0,0]
		yA = points_y[0,0]
		xB = points_x[Nx-1,0]
		yB = points_y[Nx-1,0]
		xC = points_x[Nx-1,Ny-1]
		yC = points_y[Nx-1,Ny-1]
		xD = points_x[0,Ny-1]
		yD = points_y[0,Ny-1]

		trapezoid_grid(Nx, Ny, xA, yA, xB, yB, xC, yC, xD, yD, points_x, points_y, False)
#		if if_save_iterations:
#			save_to_vtk(grid_path + "_trapezoid.vtk", points_x,points_y,Nx,Ny)

		#convertation vectors caclulation
		for index in range(1,Nx-1):
			conv10x[index-1] = points_x_10[index] - points_x[index,0]
			conv10y[index-1] = points_y_10[index] - points_y[index,0]
			conv11x[index-1] = points_x_11[index] - points_x[index,Ny-1]
			conv11y[index-1] = points_y_11[index] - points_y[index,Ny-1]

		for index in range(1,Ny-1):
			conv00x[index-1] = points_x_00[index] - points_x[0,index]
			conv00y[index-1] = points_y_00[index] - points_y[0,index]
			conv01x[index-1] = points_x_01[index] - points_x[Nx-1,index]
			conv01y[index-1] = points_y_01[index] - points_y[Nx-1,index]

		#convertation trapezoid_grid to curved grid
		#copying in boundaries
		for index in range(1,Nx-1):
			points_x[index,0] = points_x_10[index]
			points_y[index,0] = points_y_10[index]
			points_x[index,Ny-1] = points_x_11[index]
			points_y[index,Ny-1] = points_y_11[index]
		for index in range(1,Ny-1):
			points_x[0,index] = points_x_00[index]
			points_y[0,index] = points_y_00[index]
			points_x[Nx-1,index] = points_x_01[index]
			points_y[Nx-1,index] = points_y_01[index]

		#using interpolation of convertation vectors for inner points
		delta_n = (float)(Nx-1)
		for i in range(1,Nx-1):
			coef = i/delta_n
			for j in range(1,Ny-1):
				points_x[i,j] = points_x[i,j] + find_scalar_interpolated(coef, conv00x[j-1], conv01x[j-1])
				points_y[i,j] = points_y[i,j] + find_scalar_interpolated(coef, conv00y[j-1], conv01y[j-1])

#		if if_save_iterations:
#			save_to_vtk(grid_path + "_hor_conv_only.vtk", points_x,points_y,Nx,Ny)

		delta_n = (float)(Ny-1)
		for j in range(1,Ny-1):
			coef = j/delta_n
			for i in range(1,Nx-1):
				points_x[i,j] = points_x[i,j] + find_scalar_interpolated(coef, conv10x[i-1], conv11x[i-1])
				points_y[i,j] = points_y[i,j] + find_scalar_interpolated(coef, conv10y[i-1], conv11y[i-1])
		
		del(conv00x)
		del(conv00y)
		del(conv01x)
		del(conv01y)
		del(conv10x)
		del(conv10y)
		del(conv11x)
		del(conv11y) 

#		if if_save_iterations_before_correction:
#			save_to_vtk(grid_path + "_before_correction.vtk", points_x,points_y,Nx,Ny)

		if (grid_correction):
			count_iteratons_critical = get_count_iteratons_critical()

			was_correction = 1
			cycle_count = 0
			total_iterations_count = 0 

			while was_correction:
				curr_corr_count = 0

				for i in range(1,Nx-1):
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "10", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "11", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "10", points_x, points_y, Nx, Ny, 0)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "11", points_x, points_y, Nx, Ny, 0)

				for j in range(1,Ny-1):
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "00", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "01", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "00", points_x, points_y, Nx, Ny, 0)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "01", points_x, points_y, Nx, Ny, 0)

				total_iterations_count = total_iterations_count + curr_corr_count
				if (not curr_corr_count):	
					was_correction = 0

#				if if_save_iterations:
#					save_to_vtk(grid_path + "_correction%s.vtk" %cycle_count, points_x,points_y,Nx,Ny)

				cycle_count = cycle_count + 1
				if (was_correction and (cycle_count>=count_iteratons_critical)):
					total_iterations_count = total_iterations_count + curr_corr_count
					sys.exit("ERROR, more iteratons than %s are needed to correct the grid" %count_iteratons_critical + "_%s_%s_" %(Nx, Ny))
					break

			result = "%s iterations including %s corrections was done to correct the grid " %(cycle_count - 1 + was_correction, total_iterations_count)

#	if print_grids:
#		print("$" + grid_path)
#		print_points(points_x,points_y,Nx,Ny)

	if (if_test_the_grid == 1):
		test_the_grid(Nx, Ny, points_x, points_y, "_%s_%s_" %(Nx, Ny))

	return result

###

def generate_curved_grid_from_boundary_points_to_array(the_array_points_x, the_array_points_y, Nx, Ny, points_x, points_y, grid_correction):
	return generate_curved_grid_from_boundary_points_separate_setting_to_array(the_array_points_x, the_array_points_y, Nx, Ny, points_x[0], points_y[0], points_x[1], points_y[1], points_x[2], points_y[2], points_x[3], points_y[3], grid_correction)

def generate_curved_grid_from_boundary_points_separate_setting_to_array(points_x, points_y, Nx, Ny, points_x_00, points_y_00, points_x_01, points_y_01, points_x_10, points_y_10, points_x_11, points_y_11, grid_correction):

	if_save_iterations = get_if_save_iterations()
	if_save_iterations_before_correction = get_if_save_iterations_before_correction()
	print_grids = get_print_grids()

	points_x[0,0] = points_x_00[0]
	points_y[0,0] = points_y_00[0]
	points_x[Nx-1,0] = points_x_01[0]
	points_y[Nx-1,0] = points_y_01[0]
	points_x[Nx-1,Ny-1] = points_x_11[Nx-1]
	points_y[Nx-1,Ny-1] = points_y_11[Nx-1]
	points_x[0,Ny-1] = points_x_11[0]
	points_y[0,Ny-1] = points_y_11[0]

	if ((Nx>2) and (Ny==2)):
		for index in range(1,Nx-1):
			points_x[index,0] = points_x_10[index]
			points_y[index,0] = points_y_10[index]
			points_x[index,1] = points_x_11[index]
			points_y[index,1] = points_y_11[index]

	elif ((Nx==2) and (Ny>2)):
		for index in range(1,Ny-1):
			points_x[0,index] = points_x_00[index]
			points_y[0,index] = points_y_00[index]
			points_x[1,index] = points_x_01[index]
			points_y[1,index] = points_y_01[index]

	else:
		conv00x = np.zeros((Ny-2), dtype = 'float')
		conv00y = np.zeros((Ny-2), dtype = 'float')
		conv01x = np.zeros((Ny-2), dtype = 'float')
		conv01y = np.zeros((Ny-2), dtype = 'float')
		conv10x = np.zeros((Nx-2), dtype = 'float')
		conv10y = np.zeros((Nx-2), dtype = 'float')
		conv11x = np.zeros((Nx-2), dtype = 'float')
		conv11y = np.zeros((Nx-2), dtype = 'float')

		xA = points_x[0,0]
		yA = points_y[0,0]
		xB = points_x[Nx-1,0]
		yB = points_y[Nx-1,0]
		xC = points_x[Nx-1,Ny-1]
		yC = points_y[Nx-1,Ny-1]
		xD = points_x[0,Ny-1]
		yD = points_y[0,Ny-1]

		trapezoid_grid(Nx, Ny, xA, yA, xB, yB, xC, yC, xD, yD, points_x, points_y, False)
		if if_save_iterations:
			save_to_vtk(grid_path + "_trapezoid.vtk", points_x,points_y,Nx,Ny)

		#convertation vectors caclulation
		for index in range(1,Nx-1):
			conv10x[index-1] = points_x_10[index] - points_x[index,0]
			conv10y[index-1] = points_y_10[index] - points_y[index,0]
			conv11x[index-1] = points_x_11[index] - points_x[index,Ny-1]
			conv11y[index-1] = points_y_11[index] - points_y[index,Ny-1]

		for index in range(1,Ny-1):
			conv00x[index-1] = points_x_00[index] - points_x[0,index]
			conv00y[index-1] = points_y_00[index] - points_y[0,index]
			conv01x[index-1] = points_x_01[index] - points_x[Nx-1,index]
			conv01y[index-1] = points_y_01[index] - points_y[Nx-1,index]

		#convertation trapezoid_grid to curved grid
		#copying in boundaries
		for index in range(1,Nx-1):
			points_x[index,0] = points_x_10[index]
			points_y[index,0] = points_y_10[index]
			points_x[index,Ny-1] = points_x_11[index]
			points_y[index,Ny-1] = points_y_11[index]
		for index in range(1,Ny-1):
			points_x[0,index] = points_x_00[index]
			points_y[0,index] = points_y_00[index]
			points_x[Nx-1,index] = points_x_01[index]
			points_y[Nx-1,index] = points_y_01[index]

		#using interpolation of convertation vectors for inner points
		delta_n = (float)(Nx-1)
		for i in range(1,Nx-1):
			coef = i/delta_n
			for j in range(1,Ny-1):
				points_x[i,j] = points_x[i,j] + find_scalar_interpolated(coef, conv00x[j-1], conv01x[j-1])
				points_y[i,j] = points_y[i,j] + find_scalar_interpolated(coef, conv00y[j-1], conv01y[j-1])

		if if_save_iterations:
			save_to_vtk(grid_path + "_hor_conv_only.vtk", points_x,points_y,Nx,Ny)

		delta_n = (float)(Ny-1)
		for j in range(1,Ny-1):
			coef = j/delta_n
			for i in range(1,Nx-1):
				points_x[i,j] = points_x[i,j] + find_scalar_interpolated(coef, conv10x[i-1], conv11x[i-1])
				points_y[i,j] = points_y[i,j] + find_scalar_interpolated(coef, conv10y[i-1], conv11y[i-1])

		if if_save_iterations_before_correction:
			save_to_vtk(grid_path + "_before_correction.vtk", points_x,points_y,Nx,Ny)

		if (grid_correction):
			count_iteratons_critical = get_count_iteratons_critical()
		
			was_correction = 1
			cycle_count = 0
			total_iterations_count = 0 

			while was_correction:
				curr_corr_count = 0

				for i in range(1,Nx-1):
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "10", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "11", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "10", points_x, points_y, Nx, Ny, 0)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "11", points_x, points_y, Nx, Ny, 0)

				for j in range(1,Ny-1):
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "00", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "01", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "00", points_x, points_y, Nx, Ny, 0)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "01", points_x, points_y, Nx, Ny, 0)

				total_iterations_count = total_iterations_count + curr_corr_count
				if (not curr_corr_count):	
					was_correction = 0

				if if_save_iterations:
					save_to_vtk(grid_path + "_correction%s.vtk" %cycle_count, points_x,points_y,Nx,Ny)

				cycle_count = cycle_count + 1
				if (was_correction and (cycle_count>=count_iteratons_critical)):
					total_iterations_count = total_iterations_count + curr_corr_count
					print("ERROR, more iteratons than %s are needed to correct the grid" %count_iteratons_critical + grid_path)
					break

			print("%s iterations including %s corrections was done to correct the grid " %(cycle_count - 1 + was_correction, total_iterations_count) + grid_path)

	if print_grids:
		print("$" + grid_path)
		print_points(points_x,points_y,Nx,Ny)

	return 1

###

def generate_curved_grid_from_boundary_points_00_01_in_array(Nx, Ny, points_x, points_y, grid_correction, array_x, array_y, if_test):
	return generate_curved_grid_from_boundary_points_00_01_separate_setting_in_array(Nx, Ny, points_x[0], points_y[0], points_x[1], points_y[1], grid_correction, array_x, array_y, if_test) #if_test was 1 here in prev vers

def generate_curved_grid_from_boundary_points_00_01(grid_path, Nx, Ny, points_x, points_y, grid_correction, if_test):
	return generate_curved_grid_from_boundary_points_00_01_separate_setting(grid_path, Nx, Ny, points_x[0], points_y[0], points_x[1], points_y[1], grid_correction, if_test)

def generate_curved_grid_from_boundary_points_00_01_separate_setting(grid_path, Nx, Ny, points_x_00, points_y_00, points_x_01, points_y_01, grid_correction, if_test):
	points_x = np.zeros((Nx,Ny), dtype = 'float')
	points_y = np.zeros((Nx,Ny), dtype = 'float')
	
	result = generate_curved_grid_from_boundary_points_00_01_separate_setting_in_array(Nx, Ny, points_x_00, points_y_00, points_x_01, points_y_01, grid_correction, points_x, points_y, if_test)
	
	save_to_vtk(grid_path, points_x, points_y, Nx, Ny)
	
	del(points_x)
	del(points_y)
	return 1

def generate_curved_grid_from_boundary_points_00_01_separate_setting_in_array(Nx, Ny, points_x_00, points_y_00, points_x_01, points_y_01, grid_correction, points_x, points_y, if_test):

#	if_save_iterations = get_if_save_iterations()
#	if_save_iterations_before_correction = get_if_save_iterations_before_correction()
#	print_grids = get_print_grids()

	points_x[0,0] = points_x_00[0]
	points_y[0,0] = points_y_00[0]
	points_x[Nx-1,0] = points_x_01[0]
	points_y[Nx-1,0] = points_y_01[0]
	points_x[Nx-1,Ny-1] = points_x_01[Ny-1]
	points_y[Nx-1,Ny-1] = points_y_01[Ny-1]
	points_x[0,Ny-1] = points_x_00[Ny-1]
	points_y[0,Ny-1] = points_y_00[Ny-1]

	result = 1

	if ((Nx==2) and (Ny>2)):
		for index in range(1,Ny-1):
			points_x[0,index] = points_x_00[index]
			points_y[0,index] = points_y_00[index]
			points_x[1,index] = points_x_01[index]
			points_y[1,index] = points_y_01[index]

	else:
		conv00x = np.zeros((Ny-2), dtype = 'float')
		conv00y = np.zeros((Ny-2), dtype = 'float')
		conv01x = np.zeros((Ny-2), dtype = 'float')
		conv01y = np.zeros((Ny-2), dtype = 'float')

		xA = points_x[0,0]
		yA = points_y[0,0]
		xB = points_x[Nx-1,0]
		yB = points_y[Nx-1,0]
		xC = points_x[Nx-1,Ny-1]
		yC = points_y[Nx-1,Ny-1]
		xD = points_x[0,Ny-1]
		yD = points_y[0,Ny-1]

		trapezoid_grid(Nx, Ny, xA, yA, xB, yB, xC, yC, xD, yD, points_x, points_y, False)
#		if if_save_iterations:
#			save_to_vtk(grid_path + "_trapezoid.vtk", points_x,points_y,Nx,Ny)

		#convertation vectors caclulation
		for index in range(1,Ny-1):
			conv00x[index-1] = points_x_00[index] - points_x[0,index]
			conv00y[index-1] = points_y_00[index] - points_y[0,index]
			conv01x[index-1] = points_x_01[index] - points_x[Nx-1,index]
			conv01y[index-1] = points_y_01[index] - points_y[Nx-1,index]

		#convertation trapezoid_grid to curved grid
		#copying in boundaries
		for index in range(1,Ny-1):
			points_x[0,index] = points_x_00[index]
			points_y[0,index] = points_y_00[index]
			points_x[Nx-1,index] = points_x_01[index]
			points_y[Nx-1,index] = points_y_01[index]

		#using interpolation of convertation vectors for inner points
		delta_n = (float)(Nx-1)
		for i in range(1,Nx-1):
			coef = i/delta_n
			for j in range(1,Ny-1):
				points_x[i,j] = points_x[i,j] + find_scalar_interpolated(coef, conv00x[j-1], conv01x[j-1])
				points_y[i,j] = points_y[i,j] + find_scalar_interpolated(coef, conv00y[j-1], conv01y[j-1])

		del(conv00x)
		del(conv00y)
		del(conv01x)
		del(conv01y)

		if (grid_correction):
			count_iteratons_critical = get_count_iteratons_critical()

			was_correction = 1
			cycle_count = 0
			total_iterations_count = 0 

			while was_correction:
				curr_corr_count = 0

				for i in range(1,Nx-1):
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "10", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "11", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "10", points_x, points_y, Nx, Ny, 0)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "11", points_x, points_y, Nx, Ny, 0)

				for j in range(1,Ny-1):
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "00", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "01", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "00", points_x, points_y, Nx, Ny, 0)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "01", points_x, points_y, Nx, Ny, 0)

				total_iterations_count = total_iterations_count + curr_corr_count
				if (not curr_corr_count):	
					was_correction = 0

#				if if_save_iterations:
#					save_to_vtk(grid_path + "_correction%s.vtk" %cycle_count, points_x,points_y,Nx,Ny)

				cycle_count = cycle_count + 1
				if (was_correction and (cycle_count>=count_iteratons_critical)):
					total_iterations_count = total_iterations_count + curr_corr_count
					sys.exit("ERROR, more iteratons than %s are needed to correct the grid" %count_iteratons_critical + "_00_01_%s_%s_" %(Nx, Ny))
					break

			result = "%s iterations including %s corrections was done to correct the grid " %(cycle_count - 1 + was_correction, total_iterations_count)

#	if print_grids:
#		print("$" + grid_path)
#		print_points(points_x,points_y,Nx,Ny)

	if (if_test == 1):
		test_the_grid(Nx, Ny, points_x, points_y, "_00_01_%s_%s_" %(Nx, Ny))
		
	return result

###

def generate_curved_grid_from_boundary_points_10_11_in_array(Nx, Ny, points_x, points_y, grid_correction, array_x, array_y, if_test):
	return generate_curved_grid_from_boundary_points_10_11_separate_setting_in_array(Nx, Ny, points_x[0], points_y[0], points_x[1], points_y[1], grid_correction, array_x, array_y, if_test)

def generate_curved_grid_from_boundary_points_10_11(grid_path, Nx, Ny, points_x, points_y, grid_correction, if_test):
	return generate_curved_grid_from_boundary_points_10_11_separate_setting(grid_path, Nx, Ny, points_x[0], points_y[0], points_x[1], points_y[1], grid_correction, if_test)

def generate_curved_grid_from_boundary_points_10_11_separate_setting(grid_path, Nx, Ny, points_x_10, points_y_10, points_x_11, points_y_11, grid_correction, if_test):
	points_x = np.zeros((Nx,Ny), dtype = 'float')
	points_y = np.zeros((Nx,Ny), dtype = 'float')
	
	result = generate_curved_grid_from_boundary_points_10_11_separate_setting_in_array(Nx, Ny, points_x_10, points_y_10, points_x_11, points_y_11, grid_correction, points_x, points_y, if_test)
	
	save_to_vtk(grid_path, points_x, points_y, Nx, Ny)
	
	del(points_x)
	del(points_y)
	return 1

def generate_curved_grid_from_boundary_points_10_11_separate_setting_in_array(Nx, Ny, points_x_10, points_y_10, points_x_11, points_y_11, grid_correction, points_x, points_y, if_test):

#	if_save_iterations = get_if_save_iterations()
#	if_save_iterations_before_correction = get_if_save_iterations_before_correction()
#	print_grids = get_print_grids()

	points_x[0,0] = points_x_10[0]
	points_y[0,0] = points_y_10[0]
	points_x[Nx-1,0] = points_x_10[Nx-1]
	points_y[Nx-1,0] = points_y_10[Nx-1]
	points_x[Nx-1,Ny-1] = points_x_11[Nx-1]
	points_y[Nx-1,Ny-1] = points_y_11[Nx-1]
	points_x[0,Ny-1] = points_x_11[0]
	points_y[0,Ny-1] = points_y_11[0]

	result = 1

	if ((Nx>2) and (Ny==2)):
		for index in range(1,Nx-1):
			points_x[index,0] = points_x_10[index]
			points_y[index,0] = points_y_10[index]
			points_x[index,1] = points_x_11[index]
			points_y[index,1] = points_y_11[index]

	else:
		conv10x = np.zeros((Nx-2), dtype = 'float')
		conv10y = np.zeros((Nx-2), dtype = 'float')
		conv11x = np.zeros((Nx-2), dtype = 'float')
		conv11y = np.zeros((Nx-2), dtype = 'float')

		xA = points_x[0,0]
		yA = points_y[0,0]
		xB = points_x[Nx-1,0]
		yB = points_y[Nx-1,0]
		xC = points_x[Nx-1,Ny-1]
		yC = points_y[Nx-1,Ny-1]
		xD = points_x[0,Ny-1]
		yD = points_y[0,Ny-1]

		trapezoid_grid(Nx, Ny, xA, yA, xB, yB, xC, yC, xD, yD, points_x, points_y, False)
#		if if_save_iterations:
#			save_to_vtk(grid_path + "_trapezoid.vtk", points_x,points_y,Nx,Ny)

		#convertation vectors caclulation
		for index in range(1,Nx-1):
			conv10x[index-1] = points_x_10[index] - points_x[index,0]
			conv10y[index-1] = points_y_10[index] - points_y[index,0]
			conv11x[index-1] = points_x_11[index] - points_x[index,Ny-1]
			conv11y[index-1] = points_y_11[index] - points_y[index,Ny-1]

		#convertation trapezoid_grid to curved grid
		#copying in boundaries
		for index in range(1,Nx-1):
			points_x[index,0] = points_x_10[index]
			points_y[index,0] = points_y_10[index]
			points_x[index,Ny-1] = points_x_11[index]
			points_y[index,Ny-1] = points_y_11[index]

		#using interpolation of convertation vectors for inner points
		delta_n = (float)(Ny-1)
		for j in range(1,Ny-1):
			coef = j/delta_n
			for i in range(1,Nx-1):
				points_x[i,j] = points_x[i,j] + find_scalar_interpolated(coef, conv10x[i-1], conv11x[i-1])
				points_y[i,j] = points_y[i,j] + find_scalar_interpolated(coef, conv10y[i-1], conv11y[i-1])

		del(conv10x)
		del(conv10y)
		del(conv11x)
		del(conv11y) 

		if (grid_correction):
			count_iteratons_critical = get_count_iteratons_critical()
		
			was_correction = 1
			cycle_count = 0
			total_iterations_count = 0 

			while was_correction:
				curr_corr_count = 0

				for i in range(1,Nx-1):
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "10", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "11", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "10", points_x, points_y, Nx, Ny, 0)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "11", points_x, points_y, Nx, Ny, 0)

				for j in range(1,Ny-1):
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "00", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "01", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "00", points_x, points_y, Nx, Ny, 0)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "01", points_x, points_y, Nx, Ny, 0)

				total_iterations_count = total_iterations_count + curr_corr_count
				if (not curr_corr_count):	
					was_correction = 0

#				if if_save_iterations:
#					save_to_vtk(grid_path + "_correction%s.vtk" %cycle_count, points_x,points_y,Nx,Ny)

				cycle_count = cycle_count + 1
				if (was_correction and (cycle_count>=count_iteratons_critical)):
					total_iterations_count = total_iterations_count + curr_corr_count
					sys.exit("ERROR, more iteratons than %s are needed to correct the grid" %count_iteratons_critical + "_10_11_%s_%s_" %(Nx, Ny))
					break

			result = "%s iterations including %s corrections was done to correct the grid " %(cycle_count - 1 + was_correction, total_iterations_count)

#	if print_grids:
#		print_points(points_x,points_y,Nx,Ny)
#		print("$" + grid_path)

	if (if_test == 1):
		test_the_grid(Nx, Ny, points_x, points_y, "_10_11_%s_%s_" %(Nx, Ny))

	return result

####

def generate_dummy_grid_from_boundary_points(grid_path, N, points_x, points_y):
	two_points_x = np.zeros((N,1), dtype = 'float')
	two_points_y = np.zeros((N,1), dtype = 'float')
	for i in range(N):
		two_points_x[i,0] = points_x[i]
		two_points_y[i,0] = points_y[i]
	save_to_vtk(grid_path, two_points_x, two_points_y, N, 1)
	return 1

def generate_boundary_points_in_curve_line(N, Func, t_start, t_end, points_x, points_y): #Func(f_vector_curr, f_t):
	vector_curr = np.zeros(2, dtype = 'float')

	N_arc_length = get_N_arc_length()

	points_x_to_arc_len_calc = np.zeros((N_arc_length), dtype = 'float')
	points_y_to_arc_len_calc = np.zeros((N_arc_length), dtype = 'float')

	Func(vector_curr,t_start)
	points_x[0] = vector_curr[0]
	points_y[0] = vector_curr[1]
	Func(vector_curr,t_end)
	points_x[N-1] = vector_curr[0]
	points_y[N-1] = vector_curr[1]

	if (N>2):
		arc_len = np.zeros((N_arc_length-1), dtype = 'float')

		points_x_to_arc_len_calc[0] = points_x[0]
		points_y_to_arc_len_calc[0] = points_y[0]
		points_x_to_arc_len_calc[N_arc_length-1] = points_x[N-1]
		points_y_to_arc_len_calc[N_arc_length-1] = points_y[N-1]

		plus_curr = (t_end - t_start)/((float)(N_arc_length-1))
		t_curr = t_start + plus_curr
		for k in range(1,N_arc_length-1):
			Func(vector_curr,t_curr)
			points_x_to_arc_len_calc[k] = vector_curr[0]
			points_y_to_arc_len_calc[k] = vector_curr[1]
			t_curr = t_curr + plus_curr

		total_arc_len = get_total_arc_len(arc_len, points_x_to_arc_len_calc, points_y_to_arc_len_calc, N_arc_length)

		arc_len_one_step = total_arc_len/((float)(N-1))
		k = 0
		k_max = N_arc_length-1
		t_curr = t_start
		total_arc_len_new_start = 0.0

		for index in range(1,N-1):
			not_find = 1 
			while not_find:
				k_start = k
				total_arc_len_curr = total_arc_len_new_start
				for k in range(k_start,k_max):
					total_arc_len_new_start = total_arc_len_curr
					total_arc_len_curr = total_arc_len_curr + arc_len[k]
					diff_curr = total_arc_len_curr-index*arc_len_one_step
					if comparison(diff_curr, 0.0):
						points_x[index] = points_x_to_arc_len_calc[k+1] 
						points_y[index] = points_y_to_arc_len_calc[k+1]
						not_find = 0 
						break
					elif ( diff_curr > 0.0 ):
						t_curr_new = t_curr + plus_curr*(1.0-diff_curr/arc_len[k])
						Func(vector_curr,t_curr_new)
						points_x[index] = vector_curr[0]
						points_y[index] = vector_curr[1] 
						not_find = 0
						break
					t_curr = t_curr + plus_curr
	return 1

#Рекомендуется создавать массивы для вычисления длин дуг единожды, чтобы не делать это много раз.
#N_arc_length = get_N_arc_length()
#points_x_to_arc_len_calc = np.zeros((N_arc_length), dtype = 'float')
#points_y_to_arc_len_calc = np.zeros((N_arc_length), dtype = 'float')
#arc_len = np.zeros((N_arc_length-1), dtype = 'float')
def find_arc_len(Func, t_start, t_end, points_x_to_arc_len_calc, points_y_to_arc_len_calc, arc_len, N_arc_length):
	f_vector_curr = np.zeros(2, dtype = 'float')
	Func(f_vector_curr, t_start)
	points_x_to_arc_len_calc[0] = f_vector_curr[0]
	points_y_to_arc_len_calc[0] = f_vector_curr[1]
	Func(f_vector_curr, t_end)
	points_x_to_arc_len_calc[N_arc_length-1] = f_vector_curr[0]
	points_y_to_arc_len_calc[N_arc_length-1] = f_vector_curr[1]
	plus_curr = (t_end - t_start)/((float)(N_arc_length-1))
	t_curr = t_start + plus_curr
	for k in range(1,N_arc_length-1):
		Func(f_vector_curr, t_curr)
		points_x_to_arc_len_calc[k] = f_vector_curr[0]
		points_y_to_arc_len_calc[k] = f_vector_curr[1]
		t_curr = t_curr + plus_curr
	return get_total_arc_len(arc_len, points_x_to_arc_len_calc, points_y_to_arc_len_calc, N_arc_length), plus_curr

def find_arc_len_save_diffs(Func, t_start, t_end, points_x_to_arc_len_calc, points_y_to_arc_len_calc, arc_len, diffs_x, diffs_y, N_arc_length):
	f_vector_curr = np.zeros(2, dtype = 'float')
	Func(f_vector_curr, t_start)
	points_x_to_arc_len_calc[0] = f_vector_curr[0]
	points_y_to_arc_len_calc[0] = f_vector_curr[1]
	Func(f_vector_curr, t_end)
	points_x_to_arc_len_calc[N_arc_length-1] = f_vector_curr[0]
	points_y_to_arc_len_calc[N_arc_length-1] = f_vector_curr[1]
	plus_curr = (t_end - t_start)/((float)(N_arc_length-1))
	t_curr = t_start + plus_curr
	for k in range(1,N_arc_length-1):
		Func(f_vector_curr, t_curr)
		points_x_to_arc_len_calc[k] = f_vector_curr[0]
		points_y_to_arc_len_calc[k] = f_vector_curr[1]
		t_curr = t_curr + plus_curr
	return get_total_arc_len_save_diffs(arc_len, diffs_x, diffs_y, points_x_to_arc_len_calc, points_y_to_arc_len_calc, N_arc_length), plus_curr

###

def get_total_arc_len_save_diffs(arc_len, diffs_x, diffs_y, points_x_to_arc_len_calc, points_y_to_arc_len_calc, N_arc_length):
	total_arc_len_curr = 0.0
	for k in range(N_arc_length-1):
		diffs_x[k] = points_x_to_arc_len_calc[k+1] - points_x_to_arc_len_calc[k]
		diffs_y[k] = points_y_to_arc_len_calc[k+1] - points_y_to_arc_len_calc[k]
		arc_len[k] = sqrt(diffs_x[k] * diffs_x[k] + diffs_y[k] * diffs_y[k])
		total_arc_len_curr = total_arc_len_curr + arc_len[k]
	return total_arc_len_curr

def generate_boundary_points_in_curve_line_with_given_arrays_piecewice_linear(N, points_x, points_y, points_x_to_arc_len_calc, points_y_to_arc_len_calc, arc_len, diffs_x, diffs_y, N_arc_length, total_arc_len):
	arc_len_one_step = total_arc_len/((float)(N-1))
	k = 0
	k_max = N_arc_length-1
	total_arc_len_new_start = 0.0
	eps = 0.01

	points_x[0] = points_x_to_arc_len_calc[0]
	points_y[0] = points_y_to_arc_len_calc[0]
	points_x[N-1] = points_x_to_arc_len_calc[N_arc_length-1]
	points_y[N-1] = points_y_to_arc_len_calc[N_arc_length-1]

	for index in range(1,N-1):
		not_find = True
		while not_find:
			k_start = k
			total_arc_len_curr = total_arc_len_new_start
			for k in range(k_start,k_max):
				total_arc_len_new_start = total_arc_len_curr
				total_arc_len_curr = total_arc_len_curr + arc_len[k]
				t = (total_arc_len_curr - index * arc_len_one_step) / arc_len[k]			
				if (abs(t) < eps): #численный ноль
					points_x[index] = points_x_to_arc_len_calc[k+1] 
					points_y[index] = points_y_to_arc_len_calc[k+1]
					not_find = False
					break
				if (t > 0.5): #near the second point k+1
					points_x[index] = points_x_to_arc_len_calc[k+1] - diffs_x[k] * t
					points_y[index] = points_y_to_arc_len_calc[k+1] - diffs_y[k] * t
					not_find = False
					break
				if (t > 0): #near the first point k
					t = 1.0 - t
					points_x[index] = points_x_to_arc_len_calc[k] + diffs_x[k] * t
					points_y[index] = points_y_to_arc_len_calc[k] + diffs_y[k] * t
					not_find = False
					break
				# а если меньше нуля, то надо ехать дальше по k, not between points k and k+1, more k+1
	return 1

###

def generate_boundary_points_in_curve_line_with_given_arrays(N, Func, t_start, t_end, points_x, points_y, points_x_to_arc_len_calc, points_y_to_arc_len_calc, arc_len, N_arc_length, total_arc_len, plus_curr):
	f_vector_curr = np.zeros(2, dtype = 'float')
	arc_len_one_step = total_arc_len/((float)(N-1))
	eps = arc_len_one_step * 0.001
	k = 0
	k_max = N_arc_length-1
	t_curr = t_start
	total_arc_len_new_start = 0.0

	points_x[0] = points_x_to_arc_len_calc[0]
	points_y[0] = points_y_to_arc_len_calc[0]
	points_x[N-1] = points_x_to_arc_len_calc[N_arc_length-1]	
	points_y[N-1] = points_y_to_arc_len_calc[N_arc_length-1]

	for index in range(1,N-1):
		not_find = True 
		while not_find:
			k_start = k
			total_arc_len_curr = total_arc_len_new_start
			for k in range(k_start,k_max):
				total_arc_len_new_start = total_arc_len_curr
				total_arc_len_curr = total_arc_len_curr + arc_len[k]
				diff_curr = total_arc_len_curr-index*arc_len_one_step
				if (abs(diff_curr) < eps):
					points_x[index] = points_x_to_arc_len_calc[k+1] 
					points_y[index] = points_y_to_arc_len_calc[k+1]
					not_find = False
					break
				if ( diff_curr > 0.0 ):
					t_curr_new = t_curr + plus_curr*(1.0-diff_curr/arc_len[k])
					Func(f_vector_curr, t_curr_new)					
					points_x[index] = f_vector_curr[0] 
					points_y[index] = f_vector_curr[1]
					not_find = False
					break
				t_curr = t_curr + plus_curr
	return 1

def test_the_grid_without_print_one_array(Nx, Ny, points):
	return test_the_grid_without_print(Nx, Ny, points[0], points[1])

def test_the_grid_without_print_one_array_sys_exit(Nx, Ny, points, message):
	test_the_grid_without_print_sys_exit(Nx, Ny, points[0], points[1], message)
	return 1

def test_the_grid_without_print_sys_exit(Nx, Ny, points_x, points_y, message):
	crit = get_epsilon() #sin_critical/3.0
	for i in range (Nx-1):
		for j in range (Ny-1):
			if (test_normalized_vector_product_critical(crit, points_x[i+1,j]-points_x[i,j], points_y[i+1,j]-points_y[i,j], points_x[i,j+1]-points_x[i,j], points_y[i,j+1]-points_y[i,j])):
				sys.exit("bad grid" + message)
	for i in range (Nx-1,0,-1):
		for j in range (Ny-1,0,-1):
			if (test_normalized_vector_product_critical(crit, points_x[i-1,j]-points_x[i,j], points_y[i-1,j]-points_y[i,j], points_x[i,j-1]-points_x[i,j], points_y[i,j-1]-points_y[i,j])):
				sys.exit("bad grid" + message)
	return 1

def test_the_grid_without_print(Nx, Ny, points_x, points_y):
	count_bad_cells = 0
	crit = get_epsilon() #sin_critical/3.0
	for i in range (Nx-1):
		for j in range (Ny-1):
			count_bad_cells = count_bad_cells + test_normalized_vector_product_critical(crit, points_x[i+1,j]-points_x[i,j], points_y[i+1,j]-points_y[i,j], points_x[i,j+1]-points_x[i,j], points_y[i,j+1]-points_y[i,j])
	for i in range (Nx-1,0,-1):
		for j in range (Ny-1,0,-1):
			count_bad_cells = count_bad_cells + test_normalized_vector_product_critical(crit, points_x[i-1,j]-points_x[i,j], points_y[i-1,j]-points_y[i,j], points_x[i,j-1]-points_x[i,j], points_y[i,j-1]-points_y[i,j])
	return count_bad_cells
	
def test_the_grid_without_print_given_critical(Nx, Ny, points_x, points_y, critical):
	count_bad_cells = 0
	for i in range (Nx-1):
		for j in range (Ny-1):
			count_bad_cells = count_bad_cells + test_normalized_vector_product_critical(critical, points_x[i+1,j]-points_x[i,j], points_y[i+1,j]-points_y[i,j], points_x[i,j+1]-points_x[i,j], points_y[i,j+1]-points_y[i,j])
	for i in range (Nx-1,0,-1):
		for j in range (Ny-1,0,-1):
			count_bad_cells = count_bad_cells + test_normalized_vector_product_critical(critical, points_x[i-1,j]-points_x[i,j], points_y[i-1,j]-points_y[i,j], points_x[i,j-1]-points_x[i,j], points_y[i,j-1]-points_y[i,j])
	return count_bad_cells
	#return (count_bad_cells > 0)

def test_the_grid(Nx, Ny, points_x, points_y, grid_path):
	need_test = get_need_test()
	need_print_bad_cells_numbers = get_need_print_bad_cells_numbers()

	if (need_test and (not need_print_bad_cells_numbers)):
		count_bad_cells = 0
		crit = get_epsilon() #sin_critical/3.0
		for i in range (Nx-1):
			for j in range (Ny-1):
				count_bad_cells = count_bad_cells + test_normalized_vector_product_critical(crit, points_x[i+1,j]-points_x[i,j], points_y[i+1,j]-points_y[i,j], points_x[i,j+1]-points_x[i,j], points_y[i,j+1]-points_y[i,j])
		for i in range (Nx-1,0,-1):
			for j in range (Ny-1,0,-1):
				count_bad_cells = count_bad_cells + test_normalized_vector_product_critical(crit, points_x[i-1,j]-points_x[i,j], points_y[i-1,j]-points_y[i,j], points_x[i,j-1]-points_x[i,j], points_y[i,j-1]-points_y[i,j])
	
		if count_bad_cells:
			print("TOTAL ERROR, %s total bad cells in the grid " %count_bad_cells + grid_path)

	if (need_test and need_print_bad_cells_numbers):
		count_bad_cells = 0
		crit = get_epsilon() #sin_critical/3.0
		for i in range (Nx-1):
			for j in range (Ny-1):
				test = test_normalized_vector_product_critical(crit, points_x[i+1,j]-points_x[i,j], points_y[i+1,j]-points_y[i,j], points_x[i,j+1]-points_x[i,j], points_y[i,j+1]-points_y[i,j])
				if test:
					print(("DIR+ %s %s") %(i,j))
				count_bad_cells = count_bad_cells + test
		for i in range (Nx-1,0,-1):
			for j in range (Ny-1,0,-1):
				test = test_normalized_vector_product_critical(crit, points_x[i-1,j]-points_x[i,j], points_y[i-1,j]-points_y[i,j], points_x[i,j-1]-points_x[i,j], points_y[i,j-1]-points_y[i,j])
				if test:
					print(("DIR- %s %s") %(i,j))
				count_bad_cells = count_bad_cells + test
		if count_bad_cells:
			print("TOTAL ERROR, %s total bad cells in the grid " %count_bad_cells + grid_path)
			print()

def divide_line(N, x_start, y_start, x_end, y_end, points_x, points_y):
	points_x[0] = x_start
	points_y[0] = y_start
	points_x[N-1] = x_end
	points_y[N-1] = y_end
	delta_x = (x_end - x_start)/( (float)(N-1) )
	delta_y = (y_end - y_start)/( (float)(N-1) )

	for i in range(1,N-1):
		points_x[i] = points_x[i-1] + delta_x
		points_y[i] = points_y[i-1] + delta_y
	
	return 1

### GRID GENERATION corners
#	0	1	corner	to set
#	00	10	A	C
#	01	10	B	D
#	01	11	C	A
#	00	11	D	B

def generate_curved_grid_from_boundary_points_00_10(grid_path, Nx, Ny, points_x, points_y, xC, yC, grid_correction):
	return generate_curved_grid_from_boundary_points_00_10_separate_setting(grid_path, Nx, Ny, points_x[0], points_y[0], points_x[1], points_y[1], xC, yC, grid_correction)

def generate_curved_grid_from_boundary_points_00_10_separate_setting(grid_path, Nx, Ny, points_x_00, points_y_00, points_x_10, points_y_10, xC, yC, grid_correction):

	if_save_iterations = get_if_save_iterations()
	if_save_iterations_before_correction = get_if_save_iterations_before_correction()
	print_grids = get_print_grids()

	to_return = 0 
	if not comparison(points_x_00[0],points_x_10[0]):
		print("ERROR DATA point A, lines 00 10, component X (%s,%s) IN GRID " %(points_x_00[0],points_x_10[0]) + grid_path)
		to_return = 1
	if not comparison(points_y_00[0],points_y_10[0]):
		print("ERROR DATA point A, lines 00 10, component Y (%s,%s) IN GRID " %(points_y_00[0],points_y_10[0]) + grid_path)
		to_return = 1
	if (to_return):
		return 0

	points_x = np.zeros((Nx,Ny), dtype = 'float')
	points_y = np.zeros((Nx,Ny), dtype = 'float')

	points_x[0,0] = points_x_00[0]
	points_y[0,0] = points_y_00[0]
	points_x[Nx-1,0] = points_x_10[Nx-1]
	points_y[Nx-1,0] = points_y_10[Nx-1]
	points_x[Nx-1,Ny-1] = xC
	points_y[Nx-1,Ny-1] = yC
	points_x[0,Ny-1] = points_x_00[Ny-1]
	points_y[0,Ny-1] = points_y_00[Ny-1]

	if ((Nx>2) or (Ny>2)):
		xB = points_x[Nx-1,0]
		yB = points_y[Nx-1,0]
		xD = points_x[0,Ny-1]
		yD = points_y[0,Ny-1]

	if ((Nx>2) and (Ny==2)):
		points_x_11 = np.zeros((Nx), dtype = 'float')
		points_y_11 = np.zeros((Nx), dtype = 'float')
		divide_line(Nx, xD, yD, xC, yC, points_x_11, points_y_11)
		for index in range(1,Nx-1):
			points_x[index,0] = points_x_10[index]
			points_y[index,0] = points_y_10[index]
			points_x[index,1] = points_x_11[index]
			points_y[index,1] = points_y_11[index]

	elif ((Nx==2) and (Ny>2)):
		points_x_01 = np.zeros((Ny), dtype = 'float')
		points_y_01 = np.zeros((Ny), dtype = 'float')
		divide_line(Ny, xB, yB, xC, yC, points_x_01, points_y_01)
		for index in range(1,Ny-1):
			points_x[0,index] = points_x_00[index]
			points_y[0,index] = points_y_00[index]
			points_x[1,index] = points_x_01[index]
			points_y[1,index] = points_y_01[index]

	else:
		xA = points_x[0,0]
		yA = points_y[0,0]

		conv00x = np.zeros((Ny-2), dtype = 'float')
		conv00y = np.zeros((Ny-2), dtype = 'float')
		conv10x = np.zeros((Nx-2), dtype = 'float')
		conv10y = np.zeros((Nx-2), dtype = 'float')

		trapezoid_grid(Nx, Ny, xA, yA, xB, yB, xC, yC, xD, yD, points_x, points_y, False)
		if if_save_iterations:
			save_to_vtk(grid_path + "_trapezoid.vtk", points_x,points_y,Nx,Ny)

		#convertation vectors caclulation
		for index in range(1,Nx-1):
			conv10x[index-1] = points_x_10[index] - points_x[index,0]
			conv10y[index-1] = points_y_10[index] - points_y[index,0]

		for index in range(1,Ny-1):
			conv00x[index-1] = points_x_00[index] - points_x[0,index]
			conv00y[index-1] = points_y_00[index] - points_y[0,index]

		#convertation trapezoid_grid to curved grid
		#copying in boundaries
		for index in range(1,Nx-1):
			points_x[index,0] = points_x_10[index]
			points_y[index,0] = points_y_10[index]
		for index in range(1,Ny-1):
			points_x[0,index] = points_x_00[index]
			points_y[0,index] = points_y_00[index]

		#using interpolation of convertation vectors for inner points
		delta_n = (float)(Nx-1)
		for i in range(1,Nx-1):
			coef = i/delta_n
			for j in range(1,Ny-1):
				points_x[i,j] = points_x[i,j] + find_scalar_interpolated(coef, conv00x[j-1], 0.0)
				points_y[i,j] = points_y[i,j] + find_scalar_interpolated(coef, conv00y[j-1], 0.0)

		if if_save_iterations:
			save_to_vtk(grid_path + "_hor_conv_only.vtk", points_x,points_y,Nx,Ny)

		delta_n = (float)(Ny-1)
		for j in range(1,Ny-1):
			coef = j/delta_n
			for i in range(1,Nx-1):
				points_x[i,j] = points_x[i,j] + find_scalar_interpolated(coef, conv10x[i-1], 0.0)
				points_y[i,j] = points_y[i,j] + find_scalar_interpolated(coef, conv10y[i-1], 0.0)

		if if_save_iterations_before_correction:
			save_to_vtk(grid_path + "_before_correction.vtk", points_x,points_y,Nx,Ny)

		if (grid_correction):
			count_iteratons_critical = get_count_iteratons_critical()
		
			was_correction = 1
			cycle_count = 0
			total_iterations_count = 0 

			while was_correction:
				curr_corr_count = 0

				for i in range(1,Nx-1):
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "10", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "11", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "10", points_x, points_y, Nx, Ny, 0)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "11", points_x, points_y, Nx, Ny, 0)

				for j in range(1,Ny-1):
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "00", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "01", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "00", points_x, points_y, Nx, Ny, 0)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "01", points_x, points_y, Nx, Ny, 0)

				total_iterations_count = total_iterations_count + curr_corr_count
				if (not curr_corr_count):	
					was_correction = 0

				if if_save_iterations:
					save_to_vtk(grid_path + "_correction%s.vtk" %cycle_count, points_x,points_y,Nx,Ny)

				cycle_count = cycle_count + 1
				if (was_correction and (cycle_count>=count_iteratons_critical)):
					total_iterations_count = total_iterations_count + curr_corr_count
					print("ERROR, more iteratons than %s are needed to correct the grid" %count_iteratons_critical + grid_path)
					break

			print("%s iterations including %s corrections was done to correct the grid " %(cycle_count - 1 + was_correction, total_iterations_count) + grid_path)

	if print_grids:
		print("$" + grid_path)
		print_points(points_x,points_y,Nx,Ny)

	test_the_grid(Nx, Ny, points_x, points_y, grid_path)
	save_to_vtk(grid_path, points_x, points_y, Nx, Ny)
	return 1

def generate_curved_grid_from_boundary_points_01_10(grid_path, Nx, Ny, points_x, points_y, xD, yD, grid_correction):
	return generate_curved_grid_from_boundary_points_01_10_separate_setting(grid_path, Nx, Ny, points_x[0], points_y[0], points_x[1], points_y[1], xD, yD, grid_correction)

def generate_curved_grid_from_boundary_points_01_10_separate_setting(grid_path, Nx, Ny, points_x_01, points_y_01, points_x_10, points_y_10, xD, yD, grid_correction):

	if_save_iterations = get_if_save_iterations()
	if_save_iterations_before_correction = get_if_save_iterations_before_correction()
	print_grids = get_print_grids()

	to_return = 0 
	if not comparison(points_x_01[0],points_x_10[Nx-1]):
		print("ERROR DATA point B, lines 01 10, component X (%s,%s) IN GRID " %(points_x_01[0],points_x_10[Nx-1]) + grid_path)
		to_return = 1
	if not comparison(points_y_01[0],points_y_10[Nx-1]):
		print("ERROR DATA point B, lines 01 10, component Y (%s,%s) IN GRID " %(points_y_01[0],points_y_10[Nx-1]) + grid_path)
		to_return = 1
	if (to_return):
		return 0

	points_x = np.zeros((Nx,Ny), dtype = 'float')
	points_y = np.zeros((Nx,Ny), dtype = 'float')

	points_x[0,0] = points_x_10[0]
	points_y[0,0] = points_y_10[0]
	points_x[Nx-1,0] = points_x_01[0]
	points_y[Nx-1,0] = points_y_01[0]
	points_x[Nx-1,Ny-1] = points_x_01[Ny-1]
	points_y[Nx-1,Ny-1] = points_y_01[Ny-1]
	points_x[0,Ny-1] = xD
	points_y[0,Ny-1] = yD

	if ((Nx>2) or (Ny>2)):
		xA = points_x[0,0]
		yA = points_y[0,0]
		xC = points_x[Nx-1,Ny-1]
		yC = points_y[Nx-1,Ny-1]

	if ((Nx>2) and (Ny==2)):
		points_x_11 = np.zeros((Nx), dtype = 'float')
		points_y_11 = np.zeros((Nx), dtype = 'float')
		divide_line(Nx, xD, yD, xC, yC, points_x_11, points_y_11)
		for index in range(1,Nx-1):
			points_x[index,0] = points_x_10[index]
			points_y[index,0] = points_y_10[index]
			points_x[index,1] = points_x_11[index]
			points_y[index,1] = points_y_11[index]

	elif ((Nx==2) and (Ny>2)):
		points_x_00 = np.zeros((Ny), dtype = 'float')
		points_y_00 = np.zeros((Ny), dtype = 'float')
		divide_line(Ny, xA, yA, xD, yD, points_x_00, points_y_00)
		for index in range(1,Ny-1):
			points_x[0,index] = points_x_00[index]
			points_y[0,index] = points_y_00[index]
			points_x[1,index] = points_x_01[index]
			points_y[1,index] = points_y_01[index]

	else:
		xB = points_x[Nx-1,0]
		yB = points_y[Nx-1,0]

		conv01x = np.zeros((Ny-2), dtype = 'float')
		conv01y = np.zeros((Ny-2), dtype = 'float')
		conv10x = np.zeros((Nx-2), dtype = 'float')
		conv10y = np.zeros((Nx-2), dtype = 'float')

		trapezoid_grid(Nx, Ny, xA, yA, xB, yB, xC, yC, xD, yD, points_x, points_y, False)
		if if_save_iterations:
			save_to_vtk(grid_path + "_trapezoid.vtk", points_x,points_y,Nx,Ny)

		#convertation vectors caclulation
		for index in range(1,Nx-1):
			conv10x[index-1] = points_x_10[index] - points_x[index,0]
			conv10y[index-1] = points_y_10[index] - points_y[index,0]

		for index in range(1,Ny-1):
			conv01x[index-1] = points_x_01[index] - points_x[Nx-1,index]
			conv01y[index-1] = points_y_01[index] - points_y[Nx-1,index]

		#convertation trapezoid_grid to curved grid
		#copying in boundaries
		for index in range(1,Nx-1):
			points_x[index,0] = points_x_10[index]
			points_y[index,0] = points_y_10[index]
		for index in range(1,Ny-1):
			points_x[Nx-1,index] = points_x_01[index]
			points_y[Nx-1,index] = points_y_01[index]

		#using interpolation of convertation vectors for inner points
		delta_n = (float)(Nx-1)
		for i in range(1,Nx-1):
			coef = i/delta_n
			for j in range(1,Ny-1):
				points_x[i,j] = points_x[i,j] + find_scalar_interpolated(coef, 0.0, conv01x[j-1])
				points_y[i,j] = points_y[i,j] + find_scalar_interpolated(coef, 0.0, conv01y[j-1])

		if if_save_iterations:
			save_to_vtk(grid_path + "_hor_conv_only.vtk", points_x,points_y,Nx,Ny)

		delta_n = (float)(Ny-1)
		for j in range(1,Ny-1):
			coef = j/delta_n
			for i in range(1,Nx-1):
				points_x[i,j] = points_x[i,j] + find_scalar_interpolated(coef, conv10x[i-1], 0.0)
				points_y[i,j] = points_y[i,j] + find_scalar_interpolated(coef, conv10y[i-1], 0.0)

		if if_save_iterations_before_correction:
			save_to_vtk(grid_path + "_before_correction.vtk", points_x,points_y,Nx,Ny)

		if (grid_correction):
			count_iteratons_critical = get_count_iteratons_critical()
		
			was_correction = 1
			cycle_count = 0
			total_iterations_count = 0 

			while was_correction:
				curr_corr_count = 0

				for i in range(1,Nx-1):
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "10", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "11", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "10", points_x, points_y, Nx, Ny, 0)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "11", points_x, points_y, Nx, Ny, 0)

				for j in range(1,Ny-1):
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "00", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "01", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "00", points_x, points_y, Nx, Ny, 0)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "01", points_x, points_y, Nx, Ny, 0)

				total_iterations_count = total_iterations_count + curr_corr_count
				if (not curr_corr_count):	
					was_correction = 0

				if if_save_iterations:
					save_to_vtk(grid_path + "_correction%s.vtk" %cycle_count, points_x,points_y,Nx,Ny)

				cycle_count = cycle_count + 1
				if (was_correction and (cycle_count>=count_iteratons_critical)):
					total_iterations_count = total_iterations_count + curr_corr_count
					print("ERROR, more iteratons than %s are needed to correct the grid" %count_iteratons_critical + grid_path)
					break

			print("%s iterations including %s corrections was done to correct the grid " %(cycle_count - 1 + was_correction, total_iterations_count) + grid_path)

	if print_grids:
		print("$" + grid_path)
		print_points(points_x,points_y,Nx,Ny)

	test_the_grid(Nx, Ny, points_x, points_y, grid_path)
	save_to_vtk(grid_path, points_x, points_y, Nx, Ny)
	return 1

def generate_curved_grid_from_boundary_points_01_11(grid_path, Nx, Ny, points_x, points_y, xA, yA, grid_correction):
	return generate_curved_grid_from_boundary_points_01_11_separate_setting(grid_path, Nx, Ny, points_x[0], points_y[0], points_x[1], points_y[1], xA, yA, grid_correction)

def generate_curved_grid_from_boundary_points_01_11_separate_setting(grid_path, Nx, Ny, points_x_01, points_y_01, points_x_11, points_y_11, xA, yA, grid_correction):

	if_save_iterations = get_if_save_iterations()
	if_save_iterations_before_correction = get_if_save_iterations_before_correction()
	print_grids = get_print_grids()

	to_return = 0 
	if not comparison(points_x_01[Ny-1],points_x_11[Nx-1]):
		print("ERROR DATA point C, lines 01 11, component X (%s,%s) IN GRID " %(points_x_01[Ny-1],points_x_11[Nx-1]) + grid_path)
		to_return = 1
	if not comparison(points_y_01[Ny-1],points_y_11[Nx-1]):
		print("ERROR DATA point C, lines 01 11, component Y (%s,%s) IN GRID " %(points_y_01[Ny-1],points_y_11[Nx-1]) + grid_path)
		to_return = 1
	if (to_return):
		return 0

	points_x = np.zeros((Nx,Ny), dtype = 'float')
	points_y = np.zeros((Nx,Ny), dtype = 'float')

	points_x[0,0] = xA
	points_y[0,0] = yA
	points_x[Nx-1,0] = points_x_01[0]
	points_y[Nx-1,0] = points_y_01[0]
	points_x[Nx-1,Ny-1] = points_x_11[Nx-1]
	points_y[Nx-1,Ny-1] = points_y_11[Nx-1]
	points_x[0,Ny-1] = points_x_11[0]
	points_y[0,Ny-1] = points_y_11[0]

	if ((Nx>2) or (Ny>2)):
		xB = points_x[Nx-1,0]
		yB = points_y[Nx-1,0]
		xD = points_x[0,Ny-1]
		yD = points_y[0,Ny-1]

	if ((Nx>2) and (Ny==2)):
		points_x_10 = np.zeros((Nx), dtype = 'float')
		points_y_10 = np.zeros((Nx), dtype = 'float')
		divide_line(Nx, xA, yA, xB, yB, points_x_11, points_y_11)
		for index in range(1,Nx-1):
			points_x[index,0] = points_x_10[index]
			points_y[index,0] = points_y_10[index]
			points_x[index,1] = points_x_11[index]
			points_y[index,1] = points_y_11[index]

	elif ((Nx==2) and (Ny>2)):
		points_x_00 = np.zeros((Ny), dtype = 'float')
		points_y_00 = np.zeros((Ny), dtype = 'float')
		divide_line(Ny, xA, yA, xD, yD, points_x_00, points_y_00)
		for index in range(1,Ny-1):
			points_x[0,index] = points_x_00[index]
			points_y[0,index] = points_y_00[index]
			points_x[1,index] = points_x_01[index]
			points_y[1,index] = points_y_01[index]

	else:
		xC = points_x[Nx-1,Ny-1]
		yC = points_y[Nx-1,Ny-1]

		conv01x = np.zeros((Ny-2), dtype = 'float')
		conv01y = np.zeros((Ny-2), dtype = 'float')
		conv11x = np.zeros((Nx-2), dtype = 'float')
		conv11y = np.zeros((Nx-2), dtype = 'float')

		trapezoid_grid(Nx, Ny, xA, yA, xB, yB, xC, yC, xD, yD, points_x, points_y, False)
		if if_save_iterations:
			save_to_vtk(grid_path + "_trapezoid.vtk", points_x,points_y,Nx,Ny)

		#convertation vectors caclulation
		for index in range(1,Nx-1):
			conv11x[index-1] = points_x_11[index] - points_x[index,Ny-1]
			conv11y[index-1] = points_y_11[index] - points_y[index,Ny-1]

		for index in range(1,Ny-1):
			conv01x[index-1] = points_x_01[index] - points_x[Nx-1,index]
			conv01y[index-1] = points_y_01[index] - points_y[Nx-1,index]

		#convertation trapezoid_grid to curved grid
		#copying in boundaries
		for index in range(1,Nx-1):
			points_x[index,Ny-1] = points_x_11[index]
			points_y[index,Ny-1] = points_y_11[index]
		for index in range(1,Ny-1):
			points_x[Nx-1,index] = points_x_01[index]
			points_y[Nx-1,index] = points_y_01[index]

		#using interpolation of convertation vectors for inner points
		delta_n = (float)(Nx-1)
		for i in range(1,Nx-1):
			coef = i/delta_n
			for j in range(1,Ny-1):
				points_x[i,j] = points_x[i,j] + find_scalar_interpolated(coef, 0.0, conv01x[j-1])
				points_y[i,j] = points_y[i,j] + find_scalar_interpolated(coef, 0.0, conv01y[j-1])

		if if_save_iterations:
			save_to_vtk(grid_path + "_hor_conv_only.vtk", points_x,points_y,Nx,Ny)

		delta_n = (float)(Ny-1)
		for j in range(1,Ny-1):
			coef = j/delta_n
			for i in range(1,Nx-1):
				points_x[i,j] = points_x[i,j] + find_scalar_interpolated(coef, 0.0, conv11x[i-1])
				points_y[i,j] = points_y[i,j] + find_scalar_interpolated(coef, 0.0, conv11y[i-1])

		if if_save_iterations_before_correction:
			save_to_vtk(grid_path + "_before_correction.vtk", points_x,points_y,Nx,Ny)

		if (grid_correction):
			count_iteratons_critical = get_count_iteratons_critical()
		
			was_correction = 1
			cycle_count = 0
			total_iterations_count = 0 

			while was_correction:
				curr_corr_count = 0

				for i in range(1,Nx-1):
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "10", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "11", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "10", points_x, points_y, Nx, Ny, 0)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "11", points_x, points_y, Nx, Ny, 0)

				for j in range(1,Ny-1):
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "00", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "01", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "00", points_x, points_y, Nx, Ny, 0)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "01", points_x, points_y, Nx, Ny, 0)

				total_iterations_count = total_iterations_count + curr_corr_count
				if (not curr_corr_count):	
					was_correction = 0

				if if_save_iterations:
					save_to_vtk(grid_path + "_correction%s.vtk" %cycle_count, points_x,points_y,Nx,Ny)

				cycle_count = cycle_count + 1
				if (was_correction and (cycle_count>=count_iteratons_critical)):
					total_iterations_count = total_iterations_count + curr_corr_count
					print("ERROR, more iteratons than %s are needed to correct the grid" %count_iteratons_critical + grid_path)
					break

			print("%s iterations including %s corrections was done to correct the grid " %(cycle_count - 1 + was_correction, total_iterations_count) + grid_path)

	if print_grids:
		print("$" + grid_path)
		print_points(points_x,points_y,Nx,Ny)

	test_the_grid(Nx, Ny, points_x, points_y, grid_path)
	save_to_vtk(grid_path, points_x, points_y, Nx, Ny)
	return 1

def generate_curved_grid_from_boundary_points_00_11(grid_path, Nx, Ny, points_x, points_y, xB, yB, grid_correction):
	return generate_curved_grid_from_boundary_points_00_11_separate_setting(grid_path, Nx, Ny, points_x[0], points_y[0], points_x[1], points_y[1], xB, yB, grid_correction)

def generate_curved_grid_from_boundary_points_00_11_separate_setting(grid_path, Nx, Ny, points_x_00, points_y_00, points_x_11, points_y_11, xB, yB, grid_correction):

	if_save_iterations = get_if_save_iterations()
	if_save_iterations_before_correction = get_if_save_iterations_before_correction()
	print_grids = get_print_grids()

	to_return = 0 
	if not comparison(points_x_00[Ny-1],points_x_11[0]):
		print("ERROR DATA point D, lines 00 11, component X (%s,%s) IN GRID " %(points_x_00[Ny-1],points_x_11[0]) + grid_path)
		to_return = 1
	if not comparison(points_y_00[Ny-1],points_y_11[0]):
		print("ERROR DATA point D, lines 00 11, component Y (%s,%s) IN GRID " %(points_y_00[Ny-1],points_y_11[0]) + grid_path)
		to_return = 1
	if (to_return):
		return 0

	points_x = np.zeros((Nx,Ny), dtype = 'float')
	points_y = np.zeros((Nx,Ny), dtype = 'float')

	points_x[0,0] = points_x_00[0]
	points_y[0,0] = points_y_00[0]
	points_x[Nx-1,0] = xB
	points_y[Nx-1,0] = yB
	points_x[Nx-1,Ny-1] = points_x_11[Nx-1]
	points_y[Nx-1,Ny-1] = points_y_11[Nx-1]
	points_x[0,Ny-1] = points_x_11[0]
	points_y[0,Ny-1] = points_y_11[0]

	if ((Nx>2) or (Ny>2)):
		xA = points_x[0,0]
		yA = points_y[0,0]
		xC = points_x[Nx-1,Ny-1]
		yC = points_y[Nx-1,Ny-1]

	if ((Nx>2) and (Ny==2)):
		points_x_10 = np.zeros((Nx), dtype = 'float')
		points_y_10 = np.zeros((Nx), dtype = 'float')
		divide_line(Nx, xA, yA, xB, yB, points_x_11, points_y_11)
		for index in range(1,Nx-1):
			points_x[index,0] = points_x_10[index]
			points_y[index,0] = points_y_10[index]
			points_x[index,1] = points_x_11[index]
			points_y[index,1] = points_y_11[index]

	elif ((Nx==2) and (Ny>2)):
		points_x_01 = np.zeros((Ny), dtype = 'float')
		points_y_01 = np.zeros((Ny), dtype = 'float')
		divide_line(Ny, xB, yB, xC, yC, points_x_01, points_y_01)
		for index in range(1,Ny-1):
			points_x[0,index] = points_x_00[index]
			points_y[0,index] = points_y_00[index]
			points_x[1,index] = points_x_01[index]
			points_y[1,index] = points_y_01[index]

	else:
		xD = points_x[0,Ny-1]
		yD = points_y[0,Ny-1]

		conv00x = np.zeros((Ny-2), dtype = 'float')
		conv00y = np.zeros((Ny-2), dtype = 'float')
		conv11x = np.zeros((Nx-2), dtype = 'float')
		conv11y = np.zeros((Nx-2), dtype = 'float')

		trapezoid_grid(Nx, Ny, xA, yA, xB, yB, xC, yC, xD, yD, points_x, points_y, False)
		if if_save_iterations:
			save_to_vtk(grid_path + "_trapezoid.vtk", points_x,points_y,Nx,Ny)

		#convertation vectors caclulation
		for index in range(1,Nx-1):
			conv11x[index-1] = points_x_11[index] - points_x[index,Ny-1]
			conv11y[index-1] = points_y_11[index] - points_y[index,Ny-1]

		for index in range(1,Ny-1):
			conv00x[index-1] = points_x_00[index] - points_x[0,index]
			conv00y[index-1] = points_y_00[index] - points_y[0,index]

		#convertation trapezoid_grid to curved grid
		#copying in boundaries
		for index in range(1,Nx-1):
			points_x[index,Ny-1] = points_x_11[index]
			points_y[index,Ny-1] = points_y_11[index]
		for index in range(1,Ny-1):
			points_x[0,index] = points_x_00[index]
			points_y[0,index] = points_y_00[index]

		#using interpolation of convertation vectors for inner points
		delta_n = (float)(Nx-1)
		for i in range(1,Nx-1):
			coef = i/delta_n
			for j in range(1,Ny-1):
				points_x[i,j] = points_x[i,j] + find_scalar_interpolated(coef, conv00x[j-1], 0.0)
				points_y[i,j] = points_y[i,j] + find_scalar_interpolated(coef, conv00y[j-1], 0.0)

		if if_save_iterations:
			save_to_vtk(grid_path + "_hor_conv_only.vtk", points_x,points_y,Nx,Ny)

		delta_n = (float)(Ny-1)
		for j in range(1,Ny-1):
			coef = j/delta_n
			for i in range(1,Nx-1):
				points_x[i,j] = points_x[i,j] + find_scalar_interpolated(coef, 0.0, conv11x[i-1])
				points_y[i,j] = points_y[i,j] + find_scalar_interpolated(coef, 0.0, conv11y[i-1])

		if if_save_iterations_before_correction:
			save_to_vtk(grid_path + "_before_correction.vtk", points_x,points_y,Nx,Ny)

		if (grid_correction):
			count_iteratons_critical = get_count_iteratons_critical()
		
			was_correction = 1
			cycle_count = 0
			total_iterations_count = 0 

			while was_correction:
				curr_corr_count = 0

				for i in range(1,Nx-1):
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "10", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "11", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "10", points_x, points_y, Nx, Ny, 0)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "11", points_x, points_y, Nx, Ny, 0)

				for j in range(1,Ny-1):
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "00", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "01", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "00", points_x, points_y, Nx, Ny, 0)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "01", points_x, points_y, Nx, Ny, 0)

				total_iterations_count = total_iterations_count + curr_corr_count
				if (not curr_corr_count):	
					was_correction = 0

				if if_save_iterations:
					save_to_vtk(grid_path + "_correction%s.vtk" %cycle_count, points_x,points_y,Nx,Ny)

				cycle_count = cycle_count + 1
				if (was_correction and (cycle_count>=count_iteratons_critical)):
					total_iterations_count = total_iterations_count + curr_corr_count
					print("ERROR, more iteratons than %s are needed to correct the grid" %count_iteratons_critical + grid_path)
					break

			print("%s iterations including %s corrections was done to correct the grid " %(cycle_count - 1 + was_correction, total_iterations_count) + grid_path)

	if print_grids:
		print("$" + grid_path)
		print_points(points_x,points_y,Nx,Ny)

	test_the_grid(Nx, Ny, points_x, points_y, grid_path)
	save_to_vtk(grid_path, points_x, points_y, Nx, Ny)
	return 1

#TODO 00_only and 01_only functions and testing, including Nx=2, Ny=2. If these functions are needed, ask Alena V. Favorskaya.
#TODO testing Nx=2, Ny=2 10_only, 11_only, total test of generate_curved_grid_from_boundary_points_10_only 

def generate_curved_grid_from_boundary_points_00_only(grid_path, Nx, Ny, points_x_00, points_y_00, xB, yB, xC, yC, grid_correction):

	if_save_iterations = get_if_save_iterations()
	if_save_iterations_before_correction = get_if_save_iterations_before_correction()
	print_grids = get_print_grids()

	points_x = np.zeros((Nx,Ny), dtype = 'float')
	points_y = np.zeros((Nx,Ny), dtype = 'float')

	points_x[0,0] = points_x_00[0]
	points_y[0,0] = points_y_00[0]
	points_x[Nx-1,0] = xB
	points_y[Nx-1,0] = yB
	points_x[Nx-1,Ny-1] = xC
	points_y[Nx-1,Ny-1] = yC
	points_x[0,Ny-1] = points_x_00[Ny-1]
	points_y[0,Ny-1] = points_y_00[Ny-1]

	if (Nx>2):
		xA = points_x[0,0]
		yA = points_y[0,0]
		xD = points_x[0,Ny-1]
		yD = points_y[0,Ny-1]

	if ((Ny>2) and (Nx==2)):
		points_x_01 = np.zeros((Ny), dtype = 'float')
		points_y_01 = np.zeros((Ny), dtype = 'float')
		divide_line(Ny, xB, yB, xC, yC, points_x_01, points_y_01)
		for index in range(1,Ny-1):
			points_x[0,index] = points_x_00[index]
			points_y[0,index] = points_y_00[index]
			points_x[1,index] = points_x_01[index]
			points_y[1,index] = points_y_01[index]

	elif ((Ny==2) and (Nx>2)):
		points_x_10 = np.zeros((Nx), dtype = 'float')
		points_y_10 = np.zeros((Nx), dtype = 'float')
		points_x_11 = np.zeros((Nx), dtype = 'float')
		points_y_11 = np.zeros((Nx), dtype = 'float')
		divide_line(Nx, xA, yA, xB, yB, points_x_10, points_y_10)
		divide_line(Nx, xD, yD, xC, yC, points_x_11, points_y_11)
		for index in range(1,Nx-1):
			points_x[index,0] = points_x_10[index]
			points_y[index,0] = points_y_10[index]
			points_x[index,1] = points_x_11[index]
			points_y[index,1] = points_y_11[index]

	else:
		conv00x = np.zeros((Ny-2), dtype = 'float')
		conv00y = np.zeros((Ny-2), dtype = 'float')

		trapezoid_grid(Nx, Ny, xA, yA, xB, yB, xC, yC, xD, yD, points_x, points_y, False)
		if if_save_iterations:
			save_to_vtk(grid_path + "_trapezoid.vtk", points_x,points_y,Nx,Ny)

		#convertation vectors caclulation
		for index in range(1,Ny-1):
			conv00x[index-1] = points_x_00[index] - points_x[0,index]
			conv00y[index-1] = points_y_00[index] - points_y[0,index]

		#convertation trapezoid_grid to curved grid
		#copying in boundaries
		for index in range(1,Ny-1):
			points_x[0,index] = points_x_00[index]
			points_y[0,index] = points_y_00[index]

		#using interpolation of convertation vectors for inner points
		delta_n = (float)(Nx-1)
		for i in range(1,Nx-1):
			coef = i/delta_n
			for j in range(1,Ny-1):
				points_x[i,j] = points_x[i,j] + find_scalar_interpolated(coef, conv00x[j-1], 0.0)
				points_y[i,j] = points_y[i,j] + find_scalar_interpolated(coef, conv00y[j-1], 0.0)

		if if_save_iterations_before_correction:
			save_to_vtk(grid_path + "_before_correction.vtk", points_x,points_y,Nx,Ny)

		if (grid_correction):
			count_iteratons_critical = get_count_iteratons_critical()
		
			was_correction = 1
			cycle_count = 0
			total_iterations_count = 0 

			while was_correction:
				curr_corr_count = 0

				for i in range(1,Nx-1):
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "10", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "11", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "10", points_x, points_y, Nx, Ny, 0)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "11", points_x, points_y, Nx, Ny, 0)

				for j in range(1,Ny-1):
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "00", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "01", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "00", points_x, points_y, Nx, Ny, 0)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "01", points_x, points_y, Nx, Ny, 0)

				total_iterations_count = total_iterations_count + curr_corr_count
				if (not curr_corr_count):	
					was_correction = 0

				if if_save_iterations:
					save_to_vtk(grid_path + "_correction%s.vtk" %cycle_count, points_x,points_y,Nx,Ny)

				cycle_count = cycle_count + 1
				if (was_correction and (cycle_count>=count_iteratons_critical)):
					total_iterations_count = total_iterations_count + curr_corr_count
					print("ERROR, more iteratons than %s are needed to correct the grid" %count_iteratons_critical + grid_path)
					break

			print("%s iterations including %s corrections was done to correct the grid " %(cycle_count - 1 + was_correction, total_iterations_count) + grid_path)

	if print_grids:
		print("$" + grid_path)
		print_points(points_x,points_y,Nx,Ny)

	test_the_grid(Nx, Ny, points_x, points_y, grid_path)
	save_to_vtk(grid_path, points_x, points_y, Nx, Ny)
	return 1

def generate_curved_grid_from_boundary_points_01_only(grid_path, Nx, Ny, points_x_01, points_y_01, xA, yA, xD, yD, grid_correction):

	if_save_iterations = get_if_save_iterations()
	if_save_iterations_before_correction = get_if_save_iterations_before_correction()
	print_grids = get_print_grids()

	points_x = np.zeros((Nx,Ny), dtype = 'float')
	points_y = np.zeros((Nx,Ny), dtype = 'float')

	points_x[0,0] = xA
	points_y[0,0] = yA
	points_x[Nx-1,0] = points_x_01[0]
	points_y[Nx-1,0] = points_y_01[0]
	points_x[Nx-1,Ny-1] = points_x_01[Ny-1]
	points_y[Nx-1,Ny-1] = points_y_01[Ny-1]
	points_x[0,Ny-1] = xD
	points_y[0,Ny-1] = yD

	if (Nx>2):
		xB = points_x[Nx-1,0]
		yB = points_y[Nx-1,0]
		xC = points_x[Nx-1,Ny-1]
		yC = points_y[Nx-1,Ny-1]

	if ((Ny>2) and (Nx==2)):
		points_x_00 = np.zeros((Ny), dtype = 'float')
		points_y_00 = np.zeros((Ny), dtype = 'float')
		divide_line(Nx, xA, yA, xD, yD, points_x_00, points_y_00)
		for index in range(1,Ny-1):
			points_x[0,index] = points_x_00[index]
			points_y[0,index] = points_y_00[index]
			points_x[1,index] = points_x_01[index]
			points_y[1,index] = points_y_01[index]

	elif ((Ny==2) and (Nx>2)):
		points_x_10 = np.zeros((Nx), dtype = 'float')
		points_y_10 = np.zeros((Nx), dtype = 'float')
		points_x_11 = np.zeros((Nx), dtype = 'float')
		points_y_11 = np.zeros((Nx), dtype = 'float')
		divide_line(Nx, xA, yA, xB, yB, points_x_10, points_y_10)
		divide_line(Nx, xD, yD, xC, yC, points_x_11, points_y_11)
		for index in range(1,Nx-1):
			points_x[index,0] = points_x_10[index]
			points_y[index,0] = points_y_10[index]
			points_x[index,1] = points_x_11[index]
			points_y[index,1] = points_y_11[index]

	else:
		conv01x = np.zeros((Ny-2), dtype = 'float')
		conv01y = np.zeros((Ny-2), dtype = 'float')

		trapezoid_grid(Nx, Ny, xA, yA, xB, yB, xC, yC, xD, yD, points_x, points_y, False)
		if if_save_iterations:
			save_to_vtk(grid_path + "_trapezoid.vtk", points_x,points_y,Nx,Ny)

		#convertation vectors caclulation
		for index in range(1,Ny-1):
			conv01x[index-1] = points_x_01[index] - points_x[Nx-1,index]
			conv01y[index-1] = points_y_01[index] - points_y[Nx-1,index]

		#convertation trapezoid_grid to curved grid
		#copying in boundaries
		for index in range(1,Ny-1):
			points_x[Nx-1,index] = points_x_11[index]
			points_y[Nx-1,index] = points_y_11[index]

		#using interpolation of convertation vectors for inner points
		delta_n = (float)(Nx-1)
		for i in range(1,Nx-1):
			coef = i/delta_n
			for j in range(1,Ny-1):
				points_x[i,j] = points_x[i,j] + find_scalar_interpolated(coef, 0.0, conv01x[j-1])
				points_y[i,j] = points_y[i,j] + find_scalar_interpolated(coef, 0.0, conv01y[j-1])

		if if_save_iterations_before_correction:
			save_to_vtk(grid_path + "_before_correction.vtk", points_x,points_y,Nx,Ny)

		if (grid_correction):
			count_iteratons_critical = get_count_iteratons_critical()
		
			was_correction = 1
			cycle_count = 0
			total_iterations_count = 0 

			while was_correction:
				curr_corr_count = 0

				for i in range(1,Nx-1):
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "10", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "11", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "10", points_x, points_y, Nx, Ny, 0)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "11", points_x, points_y, Nx, Ny, 0)

				for j in range(1,Ny-1):
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "00", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "01", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "00", points_x, points_y, Nx, Ny, 0)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "01", points_x, points_y, Nx, Ny, 0)

				total_iterations_count = total_iterations_count + curr_corr_count
				if (not curr_corr_count):	
					was_correction = 0

				if if_save_iterations:
					save_to_vtk(grid_path + "_correction%s.vtk" %cycle_count, points_x,points_y,Nx,Ny)

				cycle_count = cycle_count + 1
				if (was_correction and (cycle_count>=count_iteratons_critical)):
					total_iterations_count = total_iterations_count + curr_corr_count
					print("ERROR, more iteratons than %s are needed to correct the grid" %count_iteratons_critical + grid_path)
					break

			print("%s iterations including %s corrections was done to correct the grid " %(cycle_count - 1 + was_correction, total_iterations_count) + grid_path)

	if print_grids:
		print("$" + grid_path)
		print_points(points_x,points_y,Nx,Ny)

	test_the_grid(Nx, Ny, points_x, points_y, grid_path)
	save_to_vtk(grid_path, points_x, points_y, Nx, Ny)
	return 1

def generate_curved_grid_from_boundary_points_10_only(grid_path, Nx, Ny, points_x_10, points_y_10, xD, yD, xC, yC, grid_correction):

	if_save_iterations = get_if_save_iterations()
	if_save_iterations_before_correction = get_if_save_iterations_before_correction()
	print_grids = get_print_grids()

	points_x = np.zeros((Nx,Ny), dtype = 'float')
	points_y = np.zeros((Nx,Ny), dtype = 'float')

	points_x[0,0] = points_x_10[0]
	points_y[0,0] = points_y_10[0]
	points_x[Nx-1,0] = points_x_10[Nx-1]
	points_y[Nx-1,0] = points_y_10[Nx-1]
	points_x[Nx-1,Ny-1] = xC
	points_y[Nx-1,Ny-1] = yC
	points_x[0,Ny-1] = xD
	points_y[0,Ny-1] = yD

	if (Ny>2):
		xA = points_x[0,0]
		yA = points_y[0,0]
		xB = points_x[Nx-1,0]
		yB = points_y[Nx-1,0]

	if ((Nx>2) and (Ny==2)):
		points_x_11 = np.zeros((Nx), dtype = 'float')
		points_y_11 = np.zeros((Nx), dtype = 'float')
		divide_line(Nx, xD, yD, xC, yC, points_x_11, points_y_11)
		for index in range(1,Nx-1):
			points_x[index,0] = points_x_10[index]
			points_y[index,0] = points_y_10[index]
			points_x[index,1] = points_x_11[index]
			points_y[index,1] = points_y_11[index]

	elif ((Nx==2) and (Ny>2)):
		points_x_00 = np.zeros((Ny), dtype = 'float')
		points_y_00 = np.zeros((Ny), dtype = 'float')
		points_x_01 = np.zeros((Ny), dtype = 'float')
		points_y_01 = np.zeros((Ny), dtype = 'float')
		divide_line(Ny, xA, yA, xD, yD, points_x_00, points_y_00)
		divide_line(Ny, xB, yB, xC, yC, points_x_01, points_y_01)
		for index in range(1,Ny-1):
			points_x[0,index] = points_x_00[index]
			points_y[0,index] = points_y_00[index]
			points_x[1,index] = points_x_01[index]
			points_y[1,index] = points_y_01[index]

	else:
		conv10x = np.zeros((Nx-2), dtype = 'float')
		conv10y = np.zeros((Nx-2), dtype = 'float')

		trapezoid_grid(Nx, Ny, xA, yA, xB, yB, xC, yC, xD, yD, points_x, points_y, False)
		if if_save_iterations:
			save_to_vtk(grid_path + "_trapezoid.vtk", points_x,points_y,Nx,Ny)

		#convertation vectors caclulation
		for index in range(1,Nx-1):
			conv10x[index-1] = points_x_10[index] - points_x[index,0]
			conv10y[index-1] = points_y_10[index] - points_y[index,0]

		#convertation trapezoid_grid to curved grid
		#copying in boundaries
		for index in range(1,Nx-1):
			points_x[index,0] = points_x_10[index]
			points_y[index,0] = points_y_10[index]

		#using interpolation of convertation vectors for inner points
		delta_n = (float)(Ny-1)
		for j in range(1,Ny-1):
			coef = j/delta_n
			for i in range(1,Nx-1):
				points_x[i,j] = points_x[i,j] + find_scalar_interpolated(coef, conv10x[i-1], 0.0)
				points_y[i,j] = points_y[i,j] + find_scalar_interpolated(coef, conv10y[i-1], 0.0)

		if if_save_iterations_before_correction:
			save_to_vtk(grid_path + "_before_correction.vtk", points_x,points_y,Nx,Ny)

		if (grid_correction):
			count_iteratons_critical = get_count_iteratons_critical()
		
			was_correction = 1
			cycle_count = 0
			total_iterations_count = 0 

			while was_correction:
				curr_corr_count = 0

				for i in range(1,Nx-1):
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "10", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "11", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "10", points_x, points_y, Nx, Ny, 0)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "11", points_x, points_y, Nx, Ny, 0)

				for j in range(1,Ny-1):
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "00", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "01", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "00", points_x, points_y, Nx, Ny, 0)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "01", points_x, points_y, Nx, Ny, 0)

				total_iterations_count = total_iterations_count + curr_corr_count
				if (not curr_corr_count):	
					was_correction = 0

				if if_save_iterations:
					save_to_vtk(grid_path + "_correction%s.vtk" %cycle_count, points_x,points_y,Nx,Ny)

				cycle_count = cycle_count + 1
				if (was_correction and (cycle_count>=count_iteratons_critical)):
					total_iterations_count = total_iterations_count + curr_corr_count
					print("ERROR, more iteratons than %s are needed to correct the grid" %count_iteratons_critical + grid_path)
					break

			print("%s iterations including %s corrections was done to correct the grid " %(cycle_count - 1 + was_correction, total_iterations_count) + grid_path)

	if print_grids:
		print("$" + grid_path)
		print_points(points_x,points_y,Nx,Ny)

	test_the_grid(Nx, Ny, points_x, points_y, grid_path)
	save_to_vtk(grid_path, points_x, points_y, Nx, Ny)
	return 1

def find_nonlinear_norma_interpolation(f_coef, f_value, f_angle, f_h_regular):
	if (f_angle<0.0):
		abs_angle = 0.0 - f_angle
	else:
		abs_angle = f_angle
	curr_pow = 1.0 + atan(cos(abs_angle)*f_value/f_h_regular)/pi
#	abs_angle = cos(abs_angle)*pi/2.0
#	abs_angle = (pi/2.0 + abs_angle)/4.0
#	if comparison_with_zero(pi - abs_angle, 0.01):
#		abs_angle = pi - 0.01
#	curr_tan = tan(abs_angle)
#	curr_pow = exp(curr_tan)
#	curr_pow = 1.0 + tan(abs_angle)
	return f_value * power(f_coef, curr_pow)

def generate_curved_grid_from_boundary_points_10_only_regular_top(sides0_x, sides0_y, grid_path, Nx, Ny_total, Ny_regular, hy_regular, points_x_10, points_y_10, xD, xC, y_total, N_zero_dy_correction, grid_correction):

	if_save_iterations = get_if_save_iterations()
	if_save_iterations_before_correction = get_if_save_iterations_before_correction()
	print_grids = get_print_grids()

	to_return = 0 
	if not comparison(points_x_10[0],xD):
		print("ERROR DATA point D, component X (%s,%s) IN GRID " %(points_x_10[0],xD) + grid_path)
		to_return = 1
	if not comparison(points_x_10[Nx-1],xC):
		print("ERROR DATA point C, component X (%s,%s) IN GRID " %(points_x_10[Nx-1],xC) + grid_path)
		to_return = 1
	if (to_return):
		return 0

	Ny = Ny_total - Ny_regular
	yC = y_total - Ny_regular*hy_regular
	yD = y_total - Ny_regular*hy_regular

	points_x = np.zeros((Nx,Ny), dtype = 'float')
	points_y = np.zeros((Nx,Ny), dtype = 'float')

	points_x[0,0] = points_x_10[0]
	points_y[0,0] = points_y_10[0]
	points_x[Nx-1,0] = points_x_10[Nx-1]
	points_y[Nx-1,0] = points_y_10[Nx-1]
	points_x[Nx-1,Ny-1] = xC
	points_y[Nx-1,Ny-1] = yC
	points_x[0,Ny-1] = xD
	points_y[0,Ny-1] = yD

	if (Ny>2):
		xA = points_x[0,0]
		yA = points_y[0,0]
		xB = points_x[Nx-1,0]
		yB = points_y[Nx-1,0]

	if ((Nx>2) and (Ny==2)):
		points_x_11 = np.zeros((Nx), dtype = 'float')
		points_y_11 = np.zeros((Nx), dtype = 'float')
		divide_line(Nx, xD, yD, xC, yC, points_x_11, points_y_11)
		for index in range(1,Nx-1):
			points_x[index,0] = points_x_10[index]
			points_y[index,0] = points_y_10[index]
			points_x[index,1] = points_x_11[index]
			points_y[index,1] = points_y_11[index]

	elif ((Nx==2) and (Ny>2)):
		points_x_00 = np.zeros((Ny), dtype = 'float')
		points_y_00 = np.zeros((Ny), dtype = 'float')
		points_x_01 = np.zeros((Ny), dtype = 'float')
		points_y_01 = np.zeros((Ny), dtype = 'float')
		divide_line(Ny, xA, yA, xD, yD, points_x_00, points_y_00)
		divide_line(Ny, xB, yB, xC, yC, points_x_01, points_y_01)
		for index in range(1,Ny-1):
			points_x[0,index] = points_x_00[index]
			points_y[0,index] = points_y_00[index]
			points_x[1,index] = points_x_01[index]
			points_y[1,index] = points_y_01[index]

	else:
		alpha_conv10 = np.zeros((Nx-2), dtype = 'float')
		norm_conv10 = np.zeros((Nx-2), dtype = 'float')
		sign10x = np.zeros((Nx-2), dtype = 'int')
		sign10y = np.zeros((Nx-2), dtype = 'int')

		trapezoid_grid(Nx, Ny, xA, yA, xB, yB, xC, yC, xD, yD, points_x, points_y, False)
		if if_save_iterations:
			save_to_vtk(grid_path + "_trapezoid.vtk", points_x,points_y,Nx,Ny)

		#convertation vectors caclulation
		indices_zero = []
		N_indices_zero = 0
		indices_plus_minus = []
		N_indices_plus_minus = 0
		prev_sign = 0
		for index in range(1,Nx-1):
			dx = points_x_10[index] - points_x[index,0]
			dy = points_y_10[index] - points_y[index,0]
			norm_conv10[index-1] = sqrt(dx*dx + dy*dy)
			if ( comparison(norm_conv10[index-1], 0.0) or comparison_with_zero(norm_conv10[index-1]/hy_regular, 0.0001) ):
				alpha_conv10[index-1] = 0.0
				sign10x[index-1] = 0
				sign10y[index-1] = 0
				prev_sign = 0
				indices_zero.append(index)
				N_indices_zero = N_indices_zero + 1
			elif ( comparison(dy, 0.0) or comparison_with_zero(dy/norm_conv10[index-1], 0.0001) ):
				prev_sign = 0
				indices_zero.append(index)
				N_indices_zero = N_indices_zero + 1
				if (dx>0.0):
					alpha_conv10[index-1] = pi/2.0
					sign10x[index-1] = +1
					sign10y[index-1] = -1
				else:
					alpha_conv10[index-1] = pi/2.0
					sign10x[index-1] = -1
					sign10y[index-1] = -1

			elif (dy>0.0):
				alpha_conv10[index-1] = asin(dx/norm_conv10[index-1])
				sign10x[index-1] = +1
				sign10y[index-1] = +1
				if (prev_sign == -1):
					indices_plus_minus.append(index-1)
					N_indices_plus_minus = N_indices_plus_minus + 1
				prev_sign = +1
			else:
				alpha_conv10[index-1] = asin(dx/norm_conv10[index-1])
				sign10x[index-1] = +1
				sign10y[index-1] = -1
				if (prev_sign == +1):
					indices_plus_minus.append(index-1)
					N_indices_plus_minus = N_indices_plus_minus + 1
				prev_sign = -1
		
		#convertation trapezoid_grid to curved grid
		#copying in boundaries
		for index in range(1,Nx-1):
			points_x[index,0] = points_x_10[index]
			points_y[index,0] = points_y_10[index]

		#using interpolation of convertation vectors for inner points
		delta_n = (float)(Ny-1)
		for j in range(1,Ny-1):
			coef = j/delta_n
			for i in range(1,Nx-1):
				alpha_curr = find_scalar_interpolated(coef, alpha_conv10[i-1], 0.0)
				norm_curr = find_nonlinear_norma_interpolation(1.0 - coef, norm_conv10[i-1], alpha_conv10[i-1], hy_regular)
				points_x[i,j] = points_x[i,j] + norm_curr*sin(alpha_curr)*sign10x[i-1]
				points_y[i,j] = points_y[i,j] + norm_curr*cos(alpha_curr)*sign10y[i-1]

		#TODO good tests
		Nx_local = N_zero_dy_correction*2 + 1		
		N_max = find_max(Nx_local + 1, Ny)
		arcs_for_correction_curr_x = np.zeros((4, N_max), dtype = 'float')
		arcs_for_correction_curr_y = np.zeros((4, N_max), dtype = 'float')
		points_for_correction_curr_x = np.zeros((Nx_local + 1, Ny), dtype = 'float')
		points_for_correction_curr_y = np.zeros((Nx_local + 1, Ny), dtype = 'float')
		
		for index in range(N_indices_zero): #TODO make better
			for j in range(Ny):
				arcs_for_correction_curr_x[0,j] = points_x[indices_zero[index]-N_zero_dy_correction,j]
				arcs_for_correction_curr_y[0,j] = points_y[indices_zero[index]-N_zero_dy_correction,j]
				arcs_for_correction_curr_x[1,j] = points_x[indices_zero[index]+N_zero_dy_correction,j]
				arcs_for_correction_curr_y[1,j] = points_y[indices_zero[index]+N_zero_dy_correction,j]
			for i in range(Nx_local):
				arcs_for_correction_curr_x[2,i] = points_x[indices_zero[index]-N_zero_dy_correction+i,0]
				arcs_for_correction_curr_y[2,i] = points_y[indices_zero[index]-N_zero_dy_correction+i,0]
				arcs_for_correction_curr_x[3,i] = points_x[indices_zero[index]-N_zero_dy_correction+i,Ny-1]
				arcs_for_correction_curr_y[3,i] = points_y[indices_zero[index]-N_zero_dy_correction+i,Ny-1]
			generate_curved_grid_from_boundary_points_to_array(points_for_correction_curr_x, points_for_correction_curr_y, Nx_local, Ny, arcs_for_correction_curr_x, arcs_for_correction_curr_y, grid_correction)
			for j in range(1,Ny-1):
				for i in range(Nx_local):
					points_x[indices_zero[index]-N_zero_dy_correction+i,j] = points_for_correction_curr_x[i,j]
					points_y[indices_zero[index]-N_zero_dy_correction+i,j] = points_for_correction_curr_y[i,j]

		Nx_local = Nx_local + 1
		for index in range(N_indices_plus_minus): #TODO make better
			for j in range(Ny):
				arcs_for_correction_curr_x[0,j] = points_x[indices_plus_minus[index]-N_zero_dy_correction,j]
				arcs_for_correction_curr_y[0,j] = points_y[indices_plus_minus[index]-N_zero_dy_correction,j]
				arcs_for_correction_curr_x[1,j] = points_x[indices_plus_minus[index]+N_zero_dy_correction+1,j]
				arcs_for_correction_curr_y[1,j] = points_y[indices_plus_minus[index]+N_zero_dy_correction+1,j]
			for i in range(Nx_local):
				arcs_for_correction_curr_x[2,i] = points_x[indices_plus_minus[index]-N_zero_dy_correction+i,0]
				arcs_for_correction_curr_y[2,i] = points_y[indices_plus_minus[index]-N_zero_dy_correction+i,0]
				arcs_for_correction_curr_x[3,i] = points_x[indices_plus_minus[index]-N_zero_dy_correction+i,Ny-1]
				arcs_for_correction_curr_y[3,i] = points_y[indices_plus_minus[index]-N_zero_dy_correction+i,Ny-1]
			generate_curved_grid_from_boundary_points_to_array(points_for_correction_curr_x, points_for_correction_curr_y, Nx_local, Ny, arcs_for_correction_curr_x, arcs_for_correction_curr_y, grid_correction)
			for j in range(1,Ny-1):
				for i in range(Nx_local):
					points_x[indices_plus_minus[index]-N_zero_dy_correction+i,j] = points_for_correction_curr_x[i,j]
					points_y[indices_plus_minus[index]-N_zero_dy_correction+i,j] = points_for_correction_curr_y[i,j]

#old vers, bas for big meshes
#		for i in range(N_indices_zero): #TODO make better
#			for j in range(1,Ny-1):
#				points_x[indices_zero[i],j] = (points_x[indices_zero[i]-1,j] + points_x[indices_zero[i]+1,j])/2.0
#				points_y[indices_zero[i],j] = (points_y[indices_zero[i]-1,j] + points_y[indices_zero[i]+1,j])/2.0
				
#		for i in range(N_indices_plus_minus): #TODO make better
#			for j in range(1,Ny-1):
#				points_x[indices_plus_minus[i],j] = (2.0 * points_x[indices_plus_minus[i]-1,j] + points_x[indices_plus_minus[i]+2,j])/3.0
#				points_y[indices_plus_minus[i],j] = (2.0 * points_y[indices_plus_minus[i]-1,j] + points_y[indices_plus_minus[i]+2,j])/3.0
#				points_x[indices_plus_minus[i]+1,j] = (points_x[indices_plus_minus[i]-1,j] + 2.0 * points_x[indices_plus_minus[i]+2,j])/3.0
#				points_y[indices_plus_minus[i]+1,j] = (points_y[indices_plus_minus[i]-1,j] + 2.0 * points_y[indices_plus_minus[i]+2,j])/3.0

		if if_save_iterations_before_correction:
			save_to_vtk(grid_path + "_before_correction.vtk", points_x,points_y,Nx,Ny)

		if (grid_correction):
			count_iteratons_critical = get_count_iteratons_critical()
		
			was_correction = 1
			cycle_count = 0
			total_iterations_count = 0 

			while was_correction:
				curr_corr_count = 0

				for i in range(1,Nx-1):
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "10", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "11", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "10", points_x, points_y, Nx, Ny, 0)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "11", points_x, points_y, Nx, Ny, 0)

				for j in range(1,Ny-1):
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "00", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "01", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "00", points_x, points_y, Nx, Ny, 0)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "01", points_x, points_y, Nx, Ny, 0)

				total_iterations_count = total_iterations_count + curr_corr_count
				if (not curr_corr_count):	
					was_correction = 0

				if if_save_iterations:
					save_to_vtk(grid_path + "_correction%s.vtk" %cycle_count, points_x,points_y,Nx,Ny)

				cycle_count = cycle_count + 1
				if (was_correction and (cycle_count>=count_iteratons_critical)):
					total_iterations_count = total_iterations_count + curr_corr_count
					print("ERROR, more iteratons than %s are needed to correct the grid" %count_iteratons_critical + grid_path)
					break

			print("%s iterations including %s corrections was done to correct the grid " %(cycle_count - 1 + was_correction, total_iterations_count) + grid_path)

	if print_grids:
		print("$" + grid_path)
		print_points(points_x,points_y,Nx,Ny)

	test_the_grid(Nx, Ny, points_x, points_y, grid_path)

	points_x_final = np.zeros((Nx,Ny_total), dtype = 'float')
	points_y_final = np.zeros((Nx,Ny_total), dtype = 'float')

	for i in range(Nx):
		for j in range(Ny):
			points_x_final[i,j] = points_x[i,j]
			points_y_final[i,j] = points_y[i,j]

	for j in range(Ny_regular):
		for i in range(Nx):
			points_x_final[i,j+Ny] = points_x[i,Ny-1]
			points_y_final[i,j+Ny] = points_y[i,Ny-1] + (j+1)*hy_regular

	for j in range(Ny_total):
		sides0_x[0,j] = points_x_final[0,j]
		sides0_y[0,j] = points_y_final[0,j]
		sides0_x[1,j] = points_x_final[Nx-1,j]
		sides0_y[1,j] = points_y_final[Nx-1,j]

	test_the_grid(Nx, Ny, points_x, points_y, grid_path)
	save_to_vtk(grid_path, points_x_final, points_y_final, Nx, Ny_total)
	return 1

###

def generate_curved_grid_from_boundary_points_11_only(grid_path, Nx, Ny, points_x_11, points_y_11, xA, yA, xB, yB, grid_correction, if_test, if_save_iterations, path_save_iterations):
	points_x = np.zeros((Nx,Ny), dtype = 'float')
	points_y = np.zeros((Nx,Ny), dtype = 'float')
	
	result = generate_curved_grid_from_boundary_points_11_only_separate_setting_in_array(Nx, Ny, points_x_11, points_y_11, xA, yA, xB, yB, grid_correction, points_x, points_y, if_test, if_save_iterations, path_save_iterations)
	
	save_to_vtk(grid_path, points_x, points_y, Nx, Ny)
	
	del(points_x)
	del(points_y)
	return 1

#TODO когда понадобится переделать функцию, чтобы не только внизу сетка и  не только вертикально. И для остальных.
#conv_vertical = np.zeros((2, Ny - 2), dtype = 'float')
#N_modi = 1 -- подвинули только один ряд, только краешек сетки стал равномерным. N_modi <= Nx - 2
def move_points_11_vertical(N_modi, N_neighbor, side, points_x, points_y, Nx, Ny, if_save_iterations, path_save_iterations, ind_start, conv_vertical):
		
	if (side == 0):
		i0 = 1
		i_sign = 1
	else:
		i0 = Nx - 2
		i_sign = -1
	
	for j in range(N_neighbor):
		conv_vertical[0,j] = points_x[i0, 0] - points_x[i0, j+1] #по X должны точки сдвинуться на позицию соответствующей точки нижней стороны 10, где сетка равномерная.
		conv_vertical[1,j] = points_y[i0 - i_sign, j+1] - points_y[i0, j+1] #по Y должны точки сдвинуться на Y позицию граничной точки, где сетка тоже равномерная.
	
	div = 1.0 / ((float)(Ny - N_neighbor - 2))
	for j in range(N_neighbor, Ny-2):
		for k in range(2):
			conv_vertical[k,j] = conv_vertical[k,N_neighbor-1] * div * (Ny - 2 - j)

	if if_save_iterations:

		div = 1.0 / ((float)(N_modi))

		i = i0
		temp = div * (i0 + i_sign * N_modi - i_sign * i)
		for j in range(Ny-2):
			points_x[i,j+1] = points_x[i,j+1] + conv_vertical[0,j] * temp
			points_y[i,j+1] = points_y[i,j+1] + conv_vertical[1,j] * temp

		save_to_vtk(path_save_iterations + "_%i.vtk" %ind_start, points_x,points_y, Nx, Ny) 

		for i in range(i0 + i_sign, i0 + i_sign * N_modi, i_sign):
			temp = div * (i0 + i_sign * N_modi - i_sign * i)
			for j in range(Ny-2):
				points_x[i,j+1] = points_x[i,j+1] + conv_vertical[0,j] * temp
				points_y[i,j+1] = points_y[i,j+1] + conv_vertical[1,j] * temp

	else:
		div = 1.0 / ((float)(N_modi))
		for i in range(i0, i0 + i_sign * N_modi, i_sign):
			temp = div * (i0 + i_sign * N_modi - i_sign * i)
			for j in range(Ny-2):
				points_x[i,j+1] = points_x[i,j+1] + conv_vertical[0,j] * temp
				points_y[i,j+1] = points_y[i,j+1] + conv_vertical[1,j] * temp
			
	if if_save_iterations:
		save_to_vtk(path_save_iterations + "_%i.vtk" %(ind_start+1), points_x,points_y, Nx, Ny)	
	
	return ind_start + 2

#conv = np.zeros((2, Nx - 2), dtype = 'float')
def move_points_11_horisontal(step_y, points_x, points_y, Nx, Ny, if_save_iterations, path_save_iterations, ind_start, conv):
	temp = step_y + points_y[0,0]
	for i in range(Nx):
		conv[0,i-1] = points_x[i,0] - points_x[i,1] 
		conv[1,i-1] = temp - points_y[i,1]

	if if_save_iterations:

		div = 1.0 / ((float)(Ny - 2))

		j = 1
		temp = div * (Ny - 1 - j)
		for i in range(Nx):
			points_x[i,j] = points_x[i,j] + conv[0,i] * temp
			points_y[i,j] = points_y[i,j] + conv[1,i] * temp

		save_to_vtk(path_save_iterations + "_%i.vtk" %ind_start, points_x,points_y, Nx, Ny) 

		for j in range(2,Ny-1):
			temp = div * (Ny - 1 - j)
			for i in range(Nx):
				points_x[i,j] = points_x[i,j] + conv[0,i] * temp
				points_y[i,j] = points_y[i,j] + conv[1,i] * temp

	else:
		div = 1.0 / ((float)(Ny - 2))
		for j in range(1,Ny-1):
			temp = div * (Ny - 1 - j)
			for i in range(Nx):
				points_x[i,j] = points_x[i,j] + conv[0,i] * temp
				points_y[i,j] = points_y[i,j] + conv[1,i] * temp
			
	if if_save_iterations:
		save_to_vtk(path_save_iterations + "_%i.vtk" %(ind_start+1), points_x, points_y, Nx, Ny)	
	
	return ind_start + 2

#TODO find_scalar_interpolated убрать... сделать функции подвижек...
def generate_curved_grid_from_boundary_points_11_only_in_array(Nx, Ny, points_x_11, points_y_11, xA, yA, xB, yB, grid_correction, points_x, points_y, if_test, if_save_iterations, path_save_iterations):

#	if_save_iterations = get_if_save_iterations()
#	if_save_iterations_before_correction = get_if_save_iterations_before_correction()
#	print_grids = get_print_grids()

	points_x[0,0] = xA
	points_y[0,0] = yA
	points_x[Nx-1,0] = xB
	points_y[Nx-1,0] = yB
	points_x[Nx-1,Ny-1] = points_x_11[Nx-1]
	points_y[Nx-1,Ny-1] = points_y_11[Nx-1]
	points_x[0,Ny-1] = points_x_11[0]
	points_y[0,Ny-1] = points_y_11[0]

	if (Ny>2):
		xC = points_x[Nx-1,Ny-1]
		yC = points_y[Nx-1,Ny-1]
		xD = points_x[0,Ny-1]
		yD = points_y[0,Ny-1]

	if ((Nx>2) and (Ny==2)):
		points_x_10 = np.zeros((Nx), dtype = 'float')
		points_y_10 = np.zeros((Nx), dtype = 'float')
		divide_line(Nx, xA, yA, xB, yB, points_x_10, points_y_10)
		for index in range(1,Nx-1):
			points_x[index,0] = points_x_10[index]
			points_y[index,0] = points_y_10[index]
			points_x[index,1] = points_x_11[index]
			points_y[index,1] = points_y_11[index]

	elif ((Nx==2) and (Ny>2)):
		points_x_00 = np.zeros((Ny), dtype = 'float')
		points_y_00 = np.zeros((Ny), dtype = 'float')
		points_x_01 = np.zeros((Ny), dtype = 'float')
		points_y_01 = np.zeros((Ny), dtype = 'float')
		divide_line(Ny, xA, yA, xD, yD, points_x_00, points_y_00)
		divide_line(Ny, xB, yB, xC, yC, points_x_01, points_y_01)
		for index in range(1,Ny-1):
			points_x[0,index] = points_x_00[index]
			points_y[0,index] = points_y_00[index]
			points_x[1,index] = points_x_01[index]
			points_y[1,index] = points_y_01[index]

	else:
		conv11x = np.zeros((Nx-2), dtype = 'float')
		conv11y = np.zeros((Nx-2), dtype = 'float')

		trapezoid_grid(Nx, Ny, xA, yA, xB, yB, xC, yC, xD, yD, points_x, points_y, False)
		if if_save_iterations:
		#	save_to_vtk(path_save_iterations + "_trapezoid.vtk", points_x,points_y, Nx, Ny)
			save_to_vtk(path_save_iterations, points_x,points_y, Nx, Ny)

		#convertation vectors caclulation
		for index in range(1,Nx-1):
			conv11x[index-1] = points_x_11[index] - points_x[index,Ny-1]
			conv11y[index-1] = points_y_11[index] - points_y[index,Ny-1]

		#convertation trapezoid_grid to curved grid
		#copying in boundaries
		for index in range(1,Nx-1):
			points_x[index,Ny-1] = points_x_11[index]
			points_y[index,Ny-1] = points_y_11[index]

		#using interpolation of convertation vectors for inner points
		delta_n = (float)(Ny-1)
		for j in range(1,Ny-1):
			coef = j/delta_n
			for i in range(1,Nx-1):
				points_x[i,j] = points_x[i,j] + find_scalar_interpolated(coef, 0.0, conv11x[i-1])
				points_y[i,j] = points_y[i,j] + find_scalar_interpolated(coef, 0.0, conv11y[i-1])

#		if if_save_iterations_before_correction:
#			save_to_vtk(grid_path + "_before_correction.vtk", points_x,points_y,Nx,Ny)

		if (grid_correction):
			count_iteratons_critical = get_count_iteratons_critical()
		
			was_correction = 1
			cycle_count = 0
			total_iterations_count = 0 

			while was_correction:
				curr_corr_count = 0

				for i in range(1,Nx-1):
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "10", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "11", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "10", points_x, points_y, Nx, Ny, 0)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, i, "11", points_x, points_y, Nx, Ny, 0)

				for j in range(1,Ny-1):
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "00", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "01", points_x, points_y, Nx, Ny, 1)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "00", points_x, points_y, Nx, Ny, 0)
					curr_corr_count = curr_corr_count + if_bad_cell(total_iterations_count + curr_corr_count, j, "01", points_x, points_y, Nx, Ny, 0)

				total_iterations_count = total_iterations_count + curr_corr_count
				if (not curr_corr_count):	
					was_correction = 0

#				if if_save_iterations:
#					save_to_vtk(grid_path + "_correction%s.vtk" %cycle_count, points_x,points_y,Nx,Ny)

				cycle_count = cycle_count + 1
				if (was_correction and (cycle_count>=count_iteratons_critical)):
					total_iterations_count = total_iterations_count + curr_corr_count
					print("ERROR, more iteratons than %s are needed to correct the grid" %count_iteratons_critical + grid_path)
					break

			print("%s iterations including %s corrections was done to correct the grid " %(cycle_count - 1 + was_correction, total_iterations_count) + grid_path)

#	if print_grids:
#		print("$" + grid_path)
#		print_points(points_x,points_y,Nx,Ny)

	if (if_test == 1):
		test_the_grid(Nx, Ny, points_x, points_y, "_11_only_%s_%s_" %(Nx, Ny))
	return 1

###

def generate_dummy_grid_from_origin_and_size(grid_path, origin_x, origin_y, size_x, size_y): #size in meters
	all_dummy_x = np.zeros((5), dtype = 'float')
	all_dummy_y = np.zeros((5), dtype = 'float')

	#A
	all_dummy_x[0] = origin_x
	all_dummy_y[0] = origin_y

	#B
	all_dummy_x[1] = origin_x + size_x
	all_dummy_y[1] = origin_y

	#C
	all_dummy_x[2] = origin_x + size_x
	all_dummy_y[2] = origin_y + size_y

	#D
	all_dummy_x[3] = origin_x
	all_dummy_y[3] = origin_y + size_y

	#A'
	all_dummy_x[4] = origin_x
	all_dummy_y[4] = origin_y

	generate_dummy_grid_from_boundary_points(grid_path, 5, all_dummy_x, all_dummy_y)

def generate_dummy_grid_from_origin_and_length_and_axis(grid_path, origin_x, origin_y, length, axis): #length in meters, axis 0 if along OX, axis 1 if along OY
	all_dummy_x = np.zeros((2), dtype = 'float')
	all_dummy_y = np.zeros((2), dtype = 'float')

	#A
	all_dummy_x[0] = origin_x
	all_dummy_y[0] = origin_y

	#B
	all_dummy_x[1] = origin_x + (1-axis) * length
	all_dummy_y[1] = origin_y + axis * length

	generate_dummy_grid_from_boundary_points(grid_path, 2, all_dummy_x, all_dummy_y)
		
def generate_dop_points_line(points_x, points_y, N_total, N_regular, x0, y0, h0, q, dir_x, dir_y):
	h_x = np.zeros((N_total-1), dtype = 'float')
	h_y = np.zeros((N_total-1), dtype = 'float')
	h_x[0] = h0*dir_x
	h_y[0] = h0*dir_y

	for i in range(1, N_regular):
		h_x[i] = h_x[0]
		h_y[i] = h_y[0]

	for i in range(N_regular, N_total-1):
		if (i==0):
			h_x[i] = h_x[0]*q
			h_y[i] = h_y[0]*q
		else:
			h_x[i] = h_x[i-1]*q
			h_y[i] = h_y[i-1]*q

	points_x[0] = x0
	points_y[0] = y0
	for i in range(N_total-1):
		points_x[i+1] = points_x[i] + h_x[i]
		points_y[i+1] = points_y[i] + h_y[i]

	return 1

def generate_dop_hs(hs, N_total, N_regular, h0, q): #hs_x in range N_total-1
	hs[0] = h0

	for i in range(1, N_regular):
		hs[i] = hs[0]

	for i in range(N_regular, N_total-1):
		if (i==0):
			hs[i] = hs[0]*q
		else:
			hs[i] = hs[i-1]*q
	return 1

def generate_dop_points_grid_by_A(grid_path, N_total_x, N_total_y, N_regular_x, N_regular_y, xA, yA, hx0, hy0, qx, qy):
	points_x = np.zeros((N_total_x,N_total_y), dtype = 'float')
	points_y = np.zeros((N_total_x,N_total_y), dtype = 'float')

	hs_x = np.zeros((N_total_x-1), dtype = 'float')
	hs_y = np.zeros((N_total_y-1), dtype = 'float')

	generate_dop_hs(hs_x, N_total_x, N_regular_x, hx0, qx)
	generate_dop_hs(hs_y, N_total_y, N_regular_y, hy0, qy)

	for i in range(N_total_x):
		if (i==0):
			points_x[i,0] = xA
		else:
			points_x[i,0] = points_x[i-1,0] + hs_x[i-1]
		for j in range(1,N_total_y):
			points_x[i,j] = points_x[i,0]

	for j in range(N_total_y):
		if (j==0):
			points_y[0,j] = yA
		else:
			points_y[0,j] = points_y[0,j-1] + hs_y[j-1]
		for i in range(1,N_total_x):
			points_y[i,j] = points_y[0,j]
	
	save_to_vtk(grid_path, points_x, points_y, N_total_x, N_total_y)
	return 1

def generate_dop_points_grid_by_B(grid_path, N_total_x, N_total_y, N_regular_x, N_regular_y, xB, yB, hx0, hy0, qx, qy):
	points_x = np.zeros((N_total_x,N_total_y), dtype = 'float')
	points_y = np.zeros((N_total_x,N_total_y), dtype = 'float')

	hs_x = np.zeros((N_total_x-1), dtype = 'float')
	hs_y = np.zeros((N_total_y-1), dtype = 'float')

	generate_dop_hs(hs_x, N_total_x, N_regular_x, hx0, qx)
	generate_dop_hs(hs_y, N_total_y, N_regular_y, hy0, qy)

	for i in range(N_total_x):
		if (i==0):
			points_x[N_total_x-i-1,0] = xB
		else:
			points_x[N_total_x-i-1,0] = points_x[N_total_x-i,0] - hs_x[i-1]
		for j in range(1,N_total_y):
			points_x[N_total_x-i-1,j] = points_x[N_total_x-i-1,0]

	for j in range(N_total_y):
		if (j==0):
			points_y[0,j] = yB
		else:
			points_y[0,j] = points_y[0,j-1] + hs_y[j-1]
		for i in range(1,N_total_x):
			points_y[i,j] = points_y[0,j]
	
	save_to_vtk(grid_path, points_x, points_y, N_total_x, N_total_y)
	return 1

def generate_dop_points_grid_by_C(grid_path, N_total_x, N_total_y, N_regular_x, N_regular_y, xC, yC, hx0, hy0, qx, qy):
	points_x = np.zeros((N_total_x,N_total_y), dtype = 'float')
	points_y = np.zeros((N_total_x,N_total_y), dtype = 'float')

	hs_x = np.zeros((N_total_x-1), dtype = 'float')
	hs_y = np.zeros((N_total_y-1), dtype = 'float')

	generate_dop_hs(hs_x, N_total_x, N_regular_x, hx0, qx)
	generate_dop_hs(hs_y, N_total_y, N_regular_y, hy0, qy)

	for i in range(N_total_x):
		if (i==0):
			points_x[N_total_x-i-1,0] = xC
		else:
			points_x[N_total_x-i-1,0] = points_x[N_total_x-i,0] - hs_x[i-1]
		for j in range(1,N_total_y):
			points_x[N_total_x-i-1,j] = points_x[N_total_x-i-1,0]

	for j in range(N_total_y):
		if (j==0):
			points_y[0,N_total_y-j-1] = yC
		else:
			points_y[0,N_total_y-j-1] = points_y[0,N_total_y-j] - hs_y[j-1]
		for i in range(1,N_total_x):
			points_y[i,N_total_y-j-1] = points_y[0,N_total_y-j-1]
	
	save_to_vtk(grid_path, points_x, points_y, N_total_x, N_total_y)
	return 1

def generate_dop_points_grid_by_D(grid_path, N_total_x, N_total_y, N_regular_x, N_regular_y, xD, yD, hx0, hy0, qx, qy):
	points_x = np.zeros((N_total_x,N_total_y), dtype = 'float')
	points_y = np.zeros((N_total_x,N_total_y), dtype = 'float')

	hs_x = np.zeros((N_total_x-1), dtype = 'float')
	hs_y = np.zeros((N_total_y-1), dtype = 'float')

	generate_dop_hs(hs_x, N_total_x, N_regular_x, hx0, qx)
	generate_dop_hs(hs_y, N_total_y, N_regular_y, hy0, qy)

	for i in range(N_total_x):
		if (i==0):
			points_x[i,0] = xD
		else:
			points_x[i,0] = points_x[i-1,0] + hs_x[i-1]
		for j in range(1,N_total_y):
			points_x[i,j] = points_x[i,0]

	for j in range(N_total_y):
		if (j==0):
			points_y[0,N_total_y-j-1] = yD
		else:
			points_y[0,N_total_y-j-1] = points_y[0,N_total_y-j] - hs_y[j-1]
		for i in range(1,N_total_x):
			points_y[i,N_total_y-j-1] = points_y[0,N_total_y-j-1]
	
	save_to_vtk(grid_path, points_x, points_y, N_total_x, N_total_y)
	return 1

def generate_dop_points_grid_by_00(grid_path, N_total_x, N_regular_x, xAD, hx0, qx, Ny, yA, yD):
	points_x = np.zeros((N_total_x,Ny), dtype = 'float')
	points_y = np.zeros((N_total_x,Ny), dtype = 'float')

	hs_x = np.zeros((N_total_x-1), dtype = 'float')
	hy = (yD-yA)/((float)(Ny-1))

	generate_dop_hs(hs_x, N_total_x, N_regular_x, hx0, qx)

	for i in range(N_total_x):
		if (i==0):
			points_x[i,0] = xAD
		else:
			points_x[i,0] = points_x[i-1,0] + hs_x[i-1]
		for j in range(1,Ny):
			points_x[i,j] = points_x[i,0]

	for j in range(Ny):
		if (j==0):
			points_y[0,j] = yA
		elif (j==(Ny-1)):
			points_y[0,j] = yD
		else:
			points_y[0,j] = points_y[0,j-1] + hy
		for i in range(1,N_total_x):
			points_y[i,j] = points_y[0,j]
	
	save_to_vtk(grid_path, points_x, points_y, N_total_x, Ny)
	return 1

def generate_dop_points_grid_by_01(grid_path, N_total_x, N_regular_x, xBC, hx0, qx, Ny, yB, yC):
	points_x = np.zeros((N_total_x,Ny), dtype = 'float')
	points_y = np.zeros((N_total_x,Ny), dtype = 'float')

	hs_x = np.zeros((N_total_x-1), dtype = 'float')
	hy = (yC-yB)/((float)(Ny-1))

	generate_dop_hs(hs_x, N_total_x, N_regular_x, hx0, qx)

	for i in range(N_total_x):
		if (i==0):
			points_x[N_total_x-i-1,0] = xBC
		else:
			points_x[N_total_x-i-1,0] = points_x[N_total_x-i,0] - hs_x[i-1]
		for j in range(1,Ny):
			points_x[N_total_x-i-1,j] = points_x[N_total_x-i-1,0]

	for j in range(Ny):
		if (j==0):
			points_y[0,j] = yB
		elif (j==(Ny-1)):
			points_y[0,j] = yC
		else:
			points_y[0,j] = points_y[0,j-1] + hy
		for i in range(1,N_total_x):
			points_y[i,j] = points_y[0,j]
	
	save_to_vtk(grid_path, points_x, points_y, N_total_x, Ny)
	return 1

def generate_dop_points_grid_by_10(grid_path, N_total_y, N_regular_y, yAB, hy0, qy, Nx, xA, xB):
	points_x = np.zeros((Nx,N_total_y), dtype = 'float')
	points_y = np.zeros((Nx,N_total_y), dtype = 'float')

	hs_y = np.zeros((N_total_y-1), dtype = 'float')
	hx = (xB-xA)/((float)(Nx-1))

	generate_dop_hs(hs_y, N_total_y, N_regular_y, hy0, qy)

	for j in range(N_total_y):
		if (j==0):
			points_y[0,j] = yAB
		else:
			points_y[0,j] = points_y[0,j-1] + hs_y[j-1]
		for i in range(1,Nx):
			points_y[i,j] = points_y[0,j]

	for i in range(Nx):
		if (i==0):
			points_x[i,0] = xA
		elif (i==(Nx-1)):
			points_x[i,0] = xB
		else:
			points_x[i,0] = points_x[i-1,0] + hx
		for j in range(1,N_total_y):
			points_x[i,j] = points_x[i,0]
	
	save_to_vtk(grid_path, points_x, points_y, Nx, N_total_y)
	return 1

def generate_dop_points_grid_by_11(grid_path, N_total_y, N_regular_y, yDC, hy0, qy, Nx, xD, xC):
	points_x = np.zeros((Nx,N_total_y), dtype = 'float')
	points_y = np.zeros((Nx,N_total_y), dtype = 'float')

	hs_y = np.zeros((N_total_y-1), dtype = 'float')
	hx = (xC-xD)/((float)(Nx-1))

	generate_dop_hs(hs_y, N_total_y, N_regular_y, hy0, qy)

	for j in range(N_total_y):
		if (j==0):
			points_y[0,N_total_y-j-1] = yDC
		else:
			points_y[0,N_total_y-j-1] = points_y[0,N_total_y-j] - hs_y[j-1]
		for i in range(1,Nx):
			points_y[i,N_total_y-j-1] = points_y[0,N_total_y-j-1]

	for i in range(Nx):
		if (i==0):
			points_x[i,0] = xD
		elif (i==(Nx-1)):
			points_x[i,0] = xC
		else:
			points_x[i,0] = points_x[i-1,0] + hx
		for j in range(1,N_total_y):
			points_x[i,j] = points_x[i,0]
	
	save_to_vtk(grid_path, points_x, points_y, Nx, N_total_y)
	return 1

#TODO and test 10 and 11
def generate_dop_points_grid_by_points_00(grid_path, N_total_x, N_regular_x, hx0, qx, Ny, points_x_00, points_y_00):
	points_x = np.zeros((N_total_x,Ny), dtype = 'float')
	points_y = np.zeros((N_total_x,Ny), dtype = 'float')

	hs_x = np.zeros((N_total_x-1), dtype = 'float')
	generate_dop_hs(hs_x, N_total_x, N_regular_x, hx0, qx)

	for j in range(Ny):
		for i in range(N_total_x):
			if (i==0):
				points_x[i,j] = points_x_00[j]
			else:
				points_x[i,j] = points_x[i-1,j] + hs_x[i-1]

	for j in range(Ny):
		points_y[0,j] = points_y_00[j]
		for i in range(1,N_total_x):
			points_y[i,j] = points_y[0,j]
	
	save_to_vtk(grid_path, points_x, points_y, N_total_x, Ny)
	return 1

def generate_dop_points_grid_by_points_01(grid_path, N_total_x, N_regular_x, hx0, qx, Ny, points_x_01, points_y_01):
	points_x = np.zeros((N_total_x,Ny), dtype = 'float')
	points_y = np.zeros((N_total_x,Ny), dtype = 'float')

	hs_x = np.zeros((N_total_x-1), dtype = 'float')
	generate_dop_hs(hs_x, N_total_x, N_regular_x, hx0, qx)

	for j in range(Ny):
		for i in range(N_total_x):
			if (i==0):
				points_x[N_total_x-i-1,j] = points_x_01[j]
			else:
				points_x[N_total_x-i-1,j] = points_x[N_total_x-i,j] - hs_x[i-1]

	for j in range(Ny):
		points_y[0,j] = points_y_01[j]
		for i in range(1,N_total_x):
			points_y[i,j] = points_y[0,j]
	
	save_to_vtk(grid_path, points_x, points_y, N_total_x, Ny)
	return 1

