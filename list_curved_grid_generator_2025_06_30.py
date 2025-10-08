#Author:  Alena V. Favorskaya, ver 30 06 2025. 

#сделать еще сортированную версию этого списка?

#изменения в set_fracture_grid. Добавлен поиск эффективного минимального шага сетки.

def get_N_arc_length():
def get_epsilon(): 
def get_alpha_critical():
def get_alpha_step():
def get_coef_arc_len_critical():
def get_count_iteratons_critical():
def get_if_save_iterations():
def get_if_save_iterations_before_correction():
def get_if_save_all_iterations():
def get_need_test():
def get_need_print_bad_cells_numbers():
def get_interpolation_radius():
def get_sin_critical():
def get_cos_critical():
def get_float_type_in_vtk_savers():
def get_print_grids():

##NOT TESTED

#вроде не меняла

##NOT TESTED


##### save to vtk given array ######

def save_grid_txt(path, points_x, points_y, Nx, Ny):
def read_grid_txt(path):

def convergence_order_ln_err(ln_err1, ln_err2, step1, step2):
def convergence_order(e1, e2, h1, h2):

def save_to_vtk_line_width(width, grid_path, points_x, points_y, Nx, Ny):

#не работает так как при записи в файл пишется слишком длинный массив, ведь его размеры отдельно не подаются TODO

# coords = np.zeros(3 * Nx * Ny * 1, dtype=np.dtype('>f4'))
#def save_to_vtk_given_array(grid_path, points_x, points_y, Nx, Ny, coords):

# coords = np.zeros(max((3 * Nx * 2 * 1), (3 * 2 * Ny * 1)), dtype=np.dtype('>f4'))
# max((3 * Nx * 2 * 1), (3 * 2 * Ny * 1)) = 6 * max (Nx, Ny)
#def save_to_vtk_line_width_given_array(width, grid_path, points_x, points_y, Nx, Ny, coords):

def global_ind(i,j,Nx):

def write_header_vtk(grid_path, Nx, Ny):

#####

def get_total_arc_len_save_diffs(arc_len, diffs_x, diffs_y, points_x_to_arc_len_calc, points_y_to_arc_len_calc, N_arc_length):

def generate_boundary_points_in_curve_line_with_given_arrays_piecewice_linear(N, points_x, points_y, points_x_to_arc_len_calc, points_y_to_arc_len_calc, arc_len, diffs_x, diffs_y, N_arc_length, total_arc_len):

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

def find_arc_len_chimera_ortho(width, axis, side, other_points_x_to_arc_len_calc, other_points_y_to_arc_len_calc, other_arc_len, diffs_x, diffs_y, points_x_to_arc_len_calc, points_y_to_arc_len_calc, N_arc_length):

def chimera_orthogonal_other_boundary(width, axis, side, N, other_boundary_points_x, other_boundary_points_y, boundary_points_x, boundary_points_y):		

########## grid_size_based_Mach_number ############

#в момент времени time верхняя граница фронта плоской волны должна проходить через точку {pos_x, pos_y}

#в момент времени time конусы Маха от границ области не должны заходить в прямоугольник [grid_min_x, grid_max_x] x [grid_min_y, grid_max_y]

#wave_angle -- угол между фронтом волны и осью OX. Вектор движения волны с положительными компонентами {wave_sin, wave_cos}. Угол лежит в пределах от 0 до 45 градусов, иначе надо другую функцию использовать.

# wave_cos, wave_cos_y по идее можно было бы не давать на входе, но обычно они уже вычислены, и быстрее подать, чем вычислять заново

# c = cP * c_coef -- скорость рассматриваемой волны. единица для продольной либо cS / cP для поперечной

#int_coord_x, int_coord_y должно попадать в желаемую точку с целочисленными координатами будующей сетки, не обязательно {0,0}

def grid_size_based_Mach_number(wave_angle, wave_sin, wave_cos,  c, c_coef,  grid_min_x, grid_max_x,  grid_min_y, grid_max_y,  time, wave_length,  pos_x, pos_y,  step_x, step_y, int_coord_x, int_coord_y):

def grid_size_based_Mach_number_wave_left_bottom_corner(mach_angle, mach_sin, mach_cos,  wave_angle, wave_sin, wave_cos,  c, c_coef,  grid_min_x, grid_max_x,  grid_min_y, grid_max_y,  time, wave_length,  pos_x, pos_y):

def cross_lines_coordinates(x1, y1,  n1_x, n1_y,   x2, y2,  n2_x, n2_y):

def read_time_from_file(path):

############

#	points_x = np.zeros((Nx,Ny), dtype = 'float')
#	points_y = np.zeros((Nx,Ny), dtype = 'float')
def trapezoid_grid(Nx, Ny, xA, yA, xB, yB, xC, yC, xD, yD, points_x, points_y):

def hat_boundary(t, center_start, amplitude, left, right, center_length):

def float_data_from_file(file_name, row_index_I, column_index_J):

###

def param_to_file(file_id, name, value, if_int, ind):

def param_to_python_file(file_id, name, value, if_int):

# (save_count < 0) => then save = 1
def find_save(Nt, save_count):

#plus 0.5 or 0.1 or 0.0 ...
def find_save_and_correct_Nt(Nt, save_count, plus):

def print_free_boundary_2D(axis, side):
def print_free_boundary_3D(axis, side):

def print_interpolation_corrector_with_given_file_name_IFN(grid_from, grid_to, axis, side, folder_name, predictor_flag, corrector_flag, axis_corrector):

def cut_array_2D(Nx_to, Ny_to, array_to, array_from, origin_to_x, origin_to_y, origin_from_x, origin_from_y):

def cut_array_3D(Nx_to, Ny_to, Nz_to, array_to, array_from, origin_to_x, origin_to_y, origin_to_z, origin_from_x, origin_from_y, origin_from_z):

def print_material(cP, cS, rho):

def param_to_file(file_id, name, value, ind):

def generate_curved_grid_from_boundary_points_in_array_from_axis_arcs(Nx, Ny, arcs_axis0_x, arcs_axis0_y, arcs_axis1_x, arcs_axis1_y, grid_correction, array_x, array_y, if_test_the_grid, if_test_corners):

def move_limits_2D(limits_new, limits_old, delta_i, delta_j, origin_i_new, origin_j_new, origin_i_old, origin_j_old):

def move_limits_3D(limits_new, limits_old, delta_i, delta_j, delta_k, origin_i_new, origin_j_new, origin_k_new, origin_i_old, origin_j_old, origin_k_old):

#IfPrint = False for debug
#ContactAxis == -1 if without axes and sides
def print_contact_2D(ContactType, ContactName1, ContactName2, ContactAxis, ContactSide1, ContactTol, IfRegular, IfPrint):
def print_contact_3D(ContactType, ContactName1, ContactName2, ContactAxis, ContactSide1, ContactTol, IfRegular, IfPrint):

## old not tested

#Рекомендуется создавать массивы для вычисления длин дуг единожды, чтобы не делать это много раз.
#N_arc_length = get_N_arc_length()
#points_x_to_arc_len_calc = np.zeros((N_arc_length), dtype = 'float')
#points_y_to_arc_len_calc = np.zeros((N_arc_length), dtype = 'float')
#arc_len = np.zeros((N_arc_length-1), dtype = 'float')
def find_arc_len(Func, t_start, t_end, points_x_to_arc_len_calc, points_y_to_arc_len_calc): 
def find_arc_len_save_diffs(Func, t_start, t_end, points_x_to_arc_len_calc, points_y_to_arc_len_calc, arc_len, diffs_x, diffs_y, N_arc_length):
def get_total_arc_len_save_diffs(arc_len, diffs_x, diffs_y, points_x_to_arc_len_calc, points_y_to_arc_len_calc, N_arc_length):

def generate_boundary_points_in_curve_line_with_given_arrays(N, Func, t_start, t_end, points_x, points_y, points_x_to_arc_len_calc, points_y_to_arc_len_calc, arc_len, N_arc_length, total_arc_len, plus_curr):

## END of NOT TESTED

def find_eff_step_min_max_average(points_x, points_y, Nx, Ny):

def find_eff_step(points_x, points_y, Nx, Ny):

def find_eff_step_3D(points_x, points_y, points_z, Nx, Ny, Nz):

from curved_grid_generator import gen_interpolation_files_for_fracture_grid_with_dop_grids_2D, print_interpolation_correctors_for_fracture_grid_with_dop_grids_2D

#(f_Nx_end - f_Nx_start + 1) * (f_Ny_end - f_Ny_start + 1)
def save_part_to_vtk(f_grid_path, f_points_x,f_points_y,f_Nx_start,f_Nx_end,f_Ny_start,f_Ny_end):

def find_min(array, N):
def find_max(array, N):

def find_rounded_value(value, dop_order):

def IFN(folder_name, grid_from, grid_to, axis, side): #Interpolation File Name
def print_interpolation_correctors_for_fracture_grid_2D(folder_name, gnM, gnF):
def print_interpolation_correctors_to_small_rect_grid_2D(folder_name, gnM, gnF, predictor_flag, corrector_flag):
def gen_interpolation_files_for_fracture_grid_2D(folder_name, gnM, gnF, array_corners):
#gnM - big grid (medium), grid name
#gnF - small grid (for example, curve fracture grid, origin - array_corners 8 and 9, size - 6 and 7.)
def gen_interpolation_files_to_small_rect_grid_2D(folder_name, gnM, gnF, size_small_x, size_small_y, int_origin_small_x, int_origin_small_y):

def print_interpolation_correctors_for_fracture_grid_3D(folder_name, gnM, gnF):
def print_interpolation_correctors_to_small_rect_grid_3D(folder_name, gnM, gnF, predictor_flag, corrector_flag):

def print_no_reflect_filler(axis, side):

def find_N(length, step):
def find_int_rounded(length, step):
def find_int_rounded_more(length, step):
def find_int_rounded_more_int_only(length, step):
def find_int_rounded_less(length, step):
def float_rounded_more(length, step):
def signed_float_rounded_more(distance, step):
def signed_float_rounded_int(distance, step):
def find_int(length, step):
def float_rounded_less(length, step):
def float_rounded(length, step):
def find_int_N(length, step):
def find_N_and_correct_float_rounded(length, step):
#Nt, save, real_time = find_Nt_and_save(time, dt, save_count)
def find_Nt_and_save(time, dt, save_count):
def find_step(length, N):

def novel_file(filename):

def find_rounded_time_step_from_jacobian(jacobian, courant, dop_order):
def sign(x):

#нормально не доделана. ошибка вычматовая, не программистская. в текущем алгоритме шаг меняется не плавно
def smooth_curve_grid_2D(path, Nx, Ny, ind_critical, ind_regular_start, ind_regular_end, hx, hy, func_y, y_start, y_end, eps_y):
def smooth_curve_grid_3D(path, Nx, Ny, Nz, ind_critical, ind_regular_start, ind_regular_end, hx, hy, hz, func_z, z_start, z_end, eps_z):

def find_s_max(L, I, h, s_min, I1):
def find_ind_cr_gr(L, I, h, s_min):
def get_smooth_line(line, N, ind_critical, ind_regular_start, ind_regular_end, h, critical, start, end, step_bottom, step_top, N_step_bottom, N_step_top, eps):

def list_to_file_float(file_name, data, count):
def list_from_file_float(file_name, data, count):
def list_to_file_int(file_name, data, count):
def list_from_file_int(file_name, data, count):
def list_to_file(file_name, data, count):
def list_from_file(file_name, data, count):

def get_txt_file_name(folder_name, grid_from, grid_to):
def print_interpolation_corrector(grid_from, grid_to, predictor_flag, corrector_flag, axis, folder_name):
def ау:

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

#TODO
def gen_interpolation_file_rect_to_rect_2D_grid_from_limits(file_name, grid_from_id, grid_from_size, grid_from_float_params, grid_from_limits, grid_to_float_params, eps):


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

#grid_limits = np.zeros((3), dtype = 'int') 
#grid_limits[0] = ind_start
#grid_limits[1] = ind_end
#grid_limits[2] = i_step, if 0, equals 1


def gen_interpolation_file_rect_to_rect_linear_2D(file_name, grid_from_id, axis, grid_from_limits, grid_from_origin, grid_from_step, grid_from_size, grid_to_limits, grid_to_origin, grid_to_step, grid_to_size, eps):


#Attention! axis 0 or 1, not 2

def gen_interpolation_file_rect_to_rect_linear_2_3D(file_name, grid_from_id, axis, grid_from_limits, grid_from_origin, grid_from_step, grid_from_size, grid_to_limits, grid_to_origin, grid_to_step, grid_to_size, eps, delay, rect_step_z, grid_to_z_limits):


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

#grid_from_limits 1, 3, 5 does not matter
def gen_interpolation_file_copy_rect_3D_both_limits(file_name, grid_from_id, grid_from_limits, grid_to_limits):

#attention! Func_AXIS should return right index.
def gen_interpolation_file_copy_rect_3D_slices_odd_cell_number(file_name, grid_from_id, grid_to_limits, Func_I, Func_J, Func_K, axis):

#curved

def gen_interpolation_file_from_rect_to_curved_1_3D(file_name, grid_from_id, curved_limits_z, eps, points_z, size, origin_z, spacing_z):

#будет работать только если координата по z в криволинейной сетке возрастает.
def gen_interpolation_file_from_curved_to_rect_1_3D(file_name, grid_from_id, curved_limits_z, eps, points_z, size, origin_z, spacing_z):

def gen_interpolation_file_from_rect_to_curved_1_2D(file_name, grid_from_id, curved_limits_y, eps, points_y, size, origin_y, spacing_y):

def gen_interpolation_file_from_curved_to_rect_1_2D(file_name, grid_from_id, curved_limits_y, eps, points_y, size, origin_y, spacing_y):

#rect_contacts are deleted here in 2_3D case
#delay and rect_step_z are defined as follows: z_index_in_rect_grid = z_index_in_curved_grid * rect_step_z + delay
def gen_interpolation_file_from_rect_to_curved_2_3D(file_name, grid_from_id, curved_limits, eps, points_x, points_y, size, origin, spacing, delay, rect_step_z):

#curved_limits 6,7 does not matter
#eps_N - solve quadratic equation or use Newton's method. Value 0.3 is recommended, not more
#delay and rect_step_z are defined as follows: z_index_in_rect_grid = z_index_in_curved_grid * rect_step_z + delay
#if_test_axis_x, if_test_axis_y are 0 or 1
def gen_interpolation_file_from_curved_to_rect_2_3D(file_name, grid_from_id, curved_limits, eps, eps_N, points_x, points_y, size, origin, spacing, if_checked_to_draw, checked_path, delay, rect_step_z, if_test_axis_x, if_test_axis_y):

#limits_to[0] = start, example 0
#limits_to[1] = end, example N-1
#limits_to[2] = step, if 0, equals 1
def gen_interpolation_file_2_3D_from_file(file_name_result, file_name_from, axis, limits_to, Func_IND):

def create_set_of_varying_steps_diff_number_of_points_linear(step, arc_len_start, arc_len_end, set_of_coefs, count, last_fixed_coord):
def create_set_of_varying_steps_diff_number_of_points_parabolic(step, arc_len_start, arc_len_end, set_of_coefs, count, last_fixed_coord):
def create_set_of_varying_steps_diff_number_of_points_func(step, arc_len_start, arc_len_end, set_of_coefs, count, func, last_fixed_coord):
def create_set_of_varying_steps(points, number, arc_len_start, arc_len_end, set_of_coefs, count, if_evenly):

#integral is equal to 0.5
def fi_hat_half_normalized(f_x): #will be in range [0.0, 1.0], non_zero in range = [0.0, 1.0], f_x SHOULD BE POZITIVE or equal to zero

def varying_step_from_arc_len(arc_len, set_of_coefs, count):
def varying_step_from_arc_len_linear(arc_len, set_of_coefs, count):

#integral is equal to (delta_left + delta_right)*(1.0 + 0.5*delta_height) = 0.5*(delta_left + delta_right)*delta_height + 1.0*(delta_left + delta_right)
def used_hat_func(arc_len, i, set_of_coefs):

def used_hat_func_linear(arc_len, i, set_of_coefs):

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

def gen_interpolation_file_copy_rect_2D_from_limits(file_name, grid_from_id, grid_from_limits, Func_i, Func_j):

#grid_from_limits 1, 3 does not matter
def gen_interpolation_file_copy_rect_2D_both_limits(file_name, grid_from_id, grid_from_limits, grid_to_limits):

#point is the point of rotation
#angle is the angle of rotation of grid_to
def gen_interpolation_file_from_rect_to_rotated_rect_2D(file_name, grid_from_id, grid_to_limits, eps, size_from, origin_from, spacing_from, angle, point, origin_to, spacing_to):

def print_array_test_06_06_22(array):
def to_file_array_test_06_06_22(file_name, array):

#origin у повернутых сеток дается в неповернутой системе координат, как в rect
def gen_interpolation_file_from_rotated_rect_to_rect_2D(file_name, grid_from_id, grid_from_limits, eps, origin_from, spacing_from, angle, point, origin_to, size_to, spacing_to):

#c means the center of rotation
def rotated_2D(x, y, angle, cx, cy):

#curved

#size = np.zeros((2), dtype = 'int')
#origin = np.zeros((2), dtype = 'float')
#spacing = np.zeros((2), dtype = 'float')
def gen_interpolation_file_from_rect_to_curved_2D(file_name, grid_from_id, curved_limits, eps, points_x, points_y, size, origin, spacing):

#curved_limits 4,5 does not matter
#eps_N - solve quadratic equation or use Newton's method. Value 0.3 is recommended, not more
#if_test_axis_x, if_test_axis_y are 0 or 1
def gen_interpolation_file_from_curved_to_rect_2D(file_name, grid_from_id, curved_limits, eps, eps_N, points_x, points_y, size, origin, spacing, if_checked_to_draw, checked_path):

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

#TODO
def gen_interpolation_file_from_curved_to_rect_2D_rect_limits(file_name, grid_from_id, rect_limits, eps, eps_N, points_x, points_y, size, origin, spacing, if_checked_to_draw, checked_path, if_test_axis_x, if_test_axis_y):

def mesh_to_vtk_mesh_only(f,origin_x,origin_y,Nx,Ny,hx,hy):
def data_to_vtk(path,array,origin_x,origin_y,Nx,Ny,hx,hy):

#array = np.zeros((5,2), dtype = 'float') 
#5 points, 2 coords, zero point is (x,y)
def inside_quadr(array):

def inside_triangle(array, num1, num2, num3):
def area(array, num1, num2, num3):
def vector_product(vector1,vector2):
def scalar_product(vector1,vector2):
def not_inside_quadr_given_vectors(vectors, vectors0):
def get_rect_coords(num, I, J, array, origin, spacing):
def ind_I(I,num):
def ind_J(J,num):
def get_next_num(num):
def get_dis_sqrt(num1,num2,array):
def get_dis_sqrt_vectors(num1,num2,array,vector):
def dis_side_rect(array, num, eps):

def Linear_Interpolation_on_Line(dis_sq, height, length):

def Bilinear_Interpolation(res, array, vectors, vectors0, height, length, eps_sq, eps_N, info, XY):

def test_point(alpha, beta, XY, X, Y, eps_sq):
def get_XY_test(alpha, beta, alpha_beta, XYa, XYb, XYab):
def find_alpha(beta, XY, X, Y):
def find_beta(alpha, XY, X, Y):
def bad_coef(coef):

#3D simple grid, TODO optimization

def generate_2_3D_grid(f_grid_path, f_points_x, f_points_y, f_Nx, f_Ny, f_Nz, origin_z, step_z):
def generate_2_3D_grid_Y(f_grid_path, f_points_x, f_points_z, f_Nx, f_Ny, f_Nz, origin_y, step_y):
def generate_2_3D_grid_Y_in_array(array_x, array_y, array_z, f_points_x, f_points_z, f_Nx, f_Ny, f_Nz, origin_y, step_y):
def generate_3D_curved_grid_from_21_array(grid_path, array_21, z20, x00, x01, y10, y11, Nx, Ny, Nz):
def generate_3D_curved_grid_from_20_array(grid_path, array_20, z21, x00, x01, y10, y11, Nx, Ny, Nz):
def generate_3D_curved_grid_from_20_and_21_arrays(grid_path, array_20, array_21, x00, x01, y10, y11, Nx, Ny, Nz):
def generate_3D_curved_grid_from_21(grid_path, func_21, z20, x00, x01, y10, y11, Nx, Ny, Nz):
def generate_3D_curved_grid_from_20(grid_path, func_20, z21, x00, x01, y10, y11, Nx, Ny, Nz):
def generate_3D_curved_grid_from_20_and_21(grid_path, func_20, func_21, x00, x01, y10, y11, Nx, Ny, Nz):

def save_to_vtk_3D_XYZ_arrays(f_grid_path, f_points_x, f_points_y, f_points_z, f_Nx, f_Ny, f_Nz):
def save_to_vtk_3D(f_grid_path, f_points, f_Nx, f_Ny, f_Nz):
def save_to_vtk_3D_with_value(f_grid_path, f_points, f_Nx, f_Ny, f_Nz, f_value_name, f_values):
def save_to_vtk_2D_with_value(f_grid_path, f_points, f_Nx, f_Ny, f_value_name, f_values):

def save_to_vtk_value_rect(path, value, origin_x, origin_y, Nx, Ny, hx, hy, label):
def mesh_to_vtk_mesh_only(f,origin_x, origin_y, Nx, Ny, hx, hy):
def header_to_vtk(path, origin_x, origin_y, Nx, Ny, hx, hy, label):

# FRACTURE GENERATOR

def poly3(x_or_y, a, b, c, d):
def find_abcd(abcd, tan_a,x0,y0,x1,y1):
def func_delta_infinity(x):
def func_delta(x):
def func_delta_parabolic(x):
def get_coord(index,origin,h):
def test_hx_for_min(Nx, Ny, points_x, h_min):
def test_hy_for_min(Nx, Ny, points_y, h_min):
def test_line_for_min(N, line_x, line_y, alpha):

#проверка на alpha_min_degree не доабвлена, так как получается что для dop_y = 16 и адекватного угла трещин и так будет нормально.
def fracture_grid_2_3D_Y(points_x, points_y, points_z, fracture_index_len_y, fracture_index_start_y, array_corners, given_points_x, given_points_z, dop_y, hx, hy, hz, origin_x, origin_y, origin_z):

#ATTENTION! If using dop grids, fracture_dop_x, fracture_dop_y should be >= 3

#fracture_topRight_y > fracture_bottomLeft_y
#fracture_topRight_x ? fracture_bottomLeft_x

#array_corners = np.zeros((12), dtype = 'int')
#if find eff step. array_corners[11] = fracture_axis

#array_corners = np.zeros((11), dtype = 'int')
#TODO comments about its components
def set_fracture_grid(filename, array_corners, fracture_bottomLeft_x, fracture_bottomLeft_y, fracture_topRight_x, fracture_topRight_y, fracture_dop_x, fracture_dop_y, origin_x, origin_y, hx, hy, if_smooth_step, if_resetting, h_min, alpha_min_degree, step, if_find_eff_step, if_dop_grids):

def set_fracture_grid_in_array(points_x, points_y, array_corners, fracture_bottomLeft_x, fracture_bottomLeft_y, fracture_topRight_x, fracture_topRight_y, fracture_dop_x, fracture_dop_y, origin_x, origin_y, hx, hy, if_smooth_step, if_resetting, h_min, alpha_min_degree, step):

#array_corners = np.zeros((11), dtype = 'int')
#TODO comments about its components
#OLD before 20 Nov 2022
#def set_fracture_grid(filename, array_corners, fracture_bottomLeft_x, fracture_bottomLeft_y, fracture_topRight_x, fracture_topRight_y, fracture_dop_x, fracture_dop_y, origin_x, origin_y, hx, hy, if_smooth_step, if_resetting, h_min, alpha_min_degree, step):

def set_dummy_fracture(filename, fracture_bottomLeft_x, fracture_bottomLeft_y, fracture_topRight_x, fracture_topRight_y):

# FUNCS FOR TESTING

def print_points(f_points_x,f_points_y,f_Nx,f_Ny):

#HELP FUNCTIONS
# grid_path and curr_corr_count are not needed in the function, used only to scientific interest to print grid transformation in process. the algorythm might be non critical optmized

def get_good_angle(f_norm,f_dx,f_dy): # return angle between the vector and dir (1,0) TODO test this functoin
def if_bad_cell(curr_corr_count, f_index, f_side_axis, f_points_x, f_points_y, f_Nx, f_Ny, f_zero_opposite):
def sin_alpha_by_normalized_vector_product(f_x1, f_y1, f_x2, f_y2):
def alpha_by_normalized_vector_product(f_x1, f_y1, f_x2, f_y2):
def test_normalized_vector_product(f_x1, f_y1, f_x2, f_y2):
def test_normalized_vector_product_critical(critical, f_x1, f_y1, f_x2, f_y2):
def comparison(f_a,f_b): #1 if equal, 0 if not equal 
def comparison_with_zero(f_a,f_epsilon): #1 if equal, 0 if not equal 
def get_total_arc_len(f_arc_len, f_points_x_to_arc_len_calc, f_points_y_to_arc_len_calc,N_arc_length):
def find_scalar_interpolated(f_coef, f_scalarA, f_scalarB): # A if coef = 0, B if coef = 1

def save_to_vtk_line(f_grid_path, f_points_x, f_points_y, f_N):
def save_to_vtk_line_one_array(f_grid_path, f_points, f_N):
def save_to_vtk_one_array(f_grid_path, f_points, f_Nx, f_Ny):
def save_to_vtk(f_grid_path, f_points_x,f_points_y,f_Nx,f_Ny):
def save_vector_field_to_vtk(f_grid_path,f_grid_x,f_grid_y,f_vector_x,f_vector_y,f_Nx,f_Ny):
def save_vector_field_to_vtk_regular_grid(f_grid_path,f_h,f_vector_x,f_vector_y,f_Nx,f_Ny):

#TODO --- use min and max
#def find_min(f_a, f_b):
#def find_max(f_a, f_b):

#CURVED GRID GENERATOR
def generate_curved_grid(grid_path, Nx, Ny, X00, Y00, X01, Y01, X10, Y10, X11, Y11, t_start00, t_end00, t_start01, t_end01, t_start10, t_end10, t_start11, t_end11, grid_correction):
def generate_boundary_points(grid_path, Nx, Ny, X00, Y00, X01, Y01, X10, Y10, X11, Y11, t_start00, t_end00, t_start01, t_end01, t_start10, t_end10, t_start11, t_end11, grid_correction, points_x_00, points_y_00, points_x_01, points_y_01, points_x_10, points_y_10, points_x_11, points_y_11):

### GRID GENERATION

def generate_quadrangle(grid_path, Nx, Ny, xA, yA, xB, yB, xC, yC, xD, yD):



def generate_curved_grid_from_boundary_points_in_array(Nx, Ny, points_x, points_y, grid_correction, array_x, array_y, if_test_the_grid, if_test_corners):
def generate_curved_grid_from_boundary_points(grid_path, Nx, Ny, points_x, points_y, grid_correction, if_test_the_grid, if_test_corners):
def generate_curved_grid_from_boundary_points_separate_setting(grid_path, Nx, Ny, points_x_00, points_y_00, points_x_01, points_y_01, points_x_10, points_y_10, points_x_11, points_y_11, grid_correction, if_test_the_grid, if_test_corners):
def generate_curved_grid_from_boundary_points_separate_setting_in_array(Nx, Ny, points_x_00, points_y_00, points_x_01, points_y_01, points_x_10, points_y_10, points_x_11, points_y_11, grid_correction, points_x, points_y, if_test_the_grid, if_test_corners):

#to delete?
def generate_curved_grid_from_boundary_points_to_array(the_array_points_x, the_array_points_y, Nx, Ny, points_x, points_y, grid_correction):
def generate_curved_grid_from_boundary_points_separate_setting_to_array(points_x, points_y, Nx, Ny, points_x_00, points_y_00, points_x_01, points_y_01, points_x_10, points_y_10, points_x_11, points_y_11, grid_correction):
#end to delete?

def generate_curved_grid_from_boundary_points_00_01_in_array(Nx, Ny, points_x, points_y, grid_correction, array_x, array_y, if_test):
def generate_curved_grid_from_boundary_points_00_01(grid_path, Nx, Ny, points_x, points_y, grid_correction, if_test):
def generate_curved_grid_from_boundary_points_00_01_separate_setting(grid_path, Nx, Ny, points_x_00, points_y_00, points_x_01, points_y_01, grid_correction, if_test):
def generate_curved_grid_from_boundary_points_00_01_separate_setting_in_array(Nx, Ny, points_x_00, points_y_00, points_x_01, points_y_01, grid_correction, points_x, points_y, if_test):

def generate_curved_grid_from_boundary_points_10_11_in_array(Nx, Ny, points_x, points_y, grid_correction, array_x, array_y, if_test):
def generate_curved_grid_from_boundary_points_10_11(grid_path, Nx, Ny, points_x, points_y, grid_correction, if_test):
def generate_curved_grid_from_boundary_points_10_11_separate_setting(grid_path, Nx, Ny, points_x_10, points_y_10, points_x_11, points_y_11, grid_correction, if_test):
def generate_curved_grid_from_boundary_points_10_11_separate_setting_in_array(Nx, Ny, points_x_10, points_y_10, points_x_11, points_y_11, grid_correction, points_x, points_y, if_test):



def generate_dummy_grid_from_boundary_points(grid_path, N, points_x, points_y):

def generate_boundary_points_in_curve_line(N, Func, t_start, t_end, points_x, points_y): #Func(f_vector_curr, f_t):

def test_the_grid_without_print_one_array_sys_exit(Nx, Ny, points, message):
def test_the_grid_without_print_sys_exit(Nx, Ny, points_x, points_y, message):
def test_the_grid_without_print_one_array(Nx, Ny, points):
def test_the_grid_without_print(Nx, Ny, points_x, points_y):
def test_the_grid(Nx, Ny, points_x, points_y, grid_path):
def test_the_grid_without_print_given_critical(Nx, Ny, points_x, points_y, critical):

def divide_line(N, x_start, y_start, x_end, y_end, points_x, points_y):

### GRID GENERATION corners
#	0	1	corner	to set
#	00	10	A	C
#	01	10	B	D
#	01	11	C	A
#	00	11	D	B

def generate_curved_grid_from_boundary_points_00_10(grid_path, Nx, Ny, points_x, points_y, xC, yC, grid_correction):
def generate_curved_grid_from_boundary_points_00_10_separate_setting(grid_path, Nx, Ny, points_x_00, points_y_00, points_x_10, points_y_10, xC, yC, grid_correction):
def generate_curved_grid_from_boundary_points_01_10(grid_path, Nx, Ny, points_x, points_y, xD, yD, grid_correction):
def generate_curved_grid_from_boundary_points_01_10_separate_setting(grid_path, Nx, Ny, points_x_01, points_y_01, points_x_10, points_y_10, xD, yD, grid_correction):
def generate_curved_grid_from_boundary_points_01_11(grid_path, Nx, Ny, points_x, points_y, xA, yA, grid_correction):
def generate_curved_grid_from_boundary_points_01_11_separate_setting(grid_path, Nx, Ny, points_x_01, points_y_01, points_x_11, points_y_11, xA, yA, grid_correction):
def generate_curved_grid_from_boundary_points_00_11(grid_path, Nx, Ny, points_x, points_y, xB, yB, grid_correction):
def generate_curved_grid_from_boundary_points_00_11_separate_setting(grid_path, Nx, Ny, points_x_00, points_y_00, points_x_11, points_y_11, xB, yB, grid_correction):

#TODO 00_only and 01_only functions and testing, including Nx=2, Ny=2. If these functions are needed, ask Alena V. Favorskaya.
#TODO testing Nx=2, Ny=2 10_only, 11_only, total test of generate_curved_grid_from_boundary_points_10_only 

def generate_curved_grid_from_boundary_points_00_only(grid_path, Nx, Ny, points_x_00, points_y_00, xB, yB, xC, yC, grid_correction):
def generate_curved_grid_from_boundary_points_01_only(grid_path, Nx, Ny, points_x_01, points_y_01, xA, yA, xD, yD, grid_correction):
def generate_curved_grid_from_boundary_points_10_only(grid_path, Nx, Ny, points_x_10, points_y_10, xD, yD, xC, yC, grid_correction):

def find_nonlinear_norma_interpolation(f_coef, f_value, f_angle, f_h_regular):

def generate_curved_grid_from_boundary_points_10_only_regular_top(sides0_x, sides0_y, grid_path, Nx, Ny_total, Ny_regular, hy_regular, points_x_10, points_y_10, xD, xC, y_total, N_zero_dy_correction, grid_correction):


###

#TODO когда понадобится переделать функцию, чтобы не только внизу сетка и  не только вертикально. И для остальных.
#conv_vertical = np.zeros((2, Ny - 2), dtype = 'float')
#N_modi = 1 -- подвинули только один ряд, только краешек сетки стал равномерным. N_modi <= Nx - 2
def move_points_11_vertical(N_modi, N_neighbor, side, points_x, points_y, Nx, Ny, if_save_iterations, path_save_iterations, ind_start, conv_vertical):

#conv = np.zeros((2, Nx - 2), dtype = 'float')
def move_points_11_horisontal(step_y, points_x, points_y, Nx, Ny, if_save_iterations, path_save_iterations, ind_start, conv):


def generate_curved_grid_from_boundary_points_11_only(grid_path, Nx, Ny, points_x_11, points_y_11, xA, yA, xB, yB, grid_correction, if_test, if_save_iterations, path_save_iterations):

def generate_curved_grid_from_boundary_points_11_only_in_array(Nx, Ny, points_x_11, points_y_11, xA, yA, xB, yB, grid_correction, points_x, points_y, if_test, if_save_iterations, path_save_iterations):

### TODO other

def generate_dummy_grid_from_origin_and_size(grid_path, origin_x, origin_y, size_x, size_y): #size in meters

def generate_dummy_grid_from_origin_and_length_and_axis(grid_path, origin_x, origin_y, length, axis): #length in meters, axis 0 if along OX, axis 1 if along OY

def generate_dop_points_line(points_x, points_y, N_total, N_regular, x0, y0, h0, q, dir_x, dir_y):
def generate_dop_hs(hs, N_total, N_regular, h0, q): #hs_x in range N_total-1
def generate_dop_points_grid_by_A(grid_path, N_total_x, N_total_y, N_regular_x, N_regular_y, xA, yA, hx0, hy0, qx, qy):
def generate_dop_points_grid_by_B(grid_path, N_total_x, N_total_y, N_regular_x, N_regular_y, xB, yB, hx0, hy0, qx, qy):
def generate_dop_points_grid_by_C(grid_path, N_total_x, N_total_y, N_regular_x, N_regular_y, xC, yC, hx0, hy0, qx, qy):
def generate_dop_points_grid_by_D(grid_path, N_total_x, N_total_y, N_regular_x, N_regular_y, xD, yD, hx0, hy0, qx, qy):
def generate_dop_points_grid_by_00(grid_path, N_total_x, N_regular_x, xAD, hx0, qx, Ny, yA, yD):
def generate_dop_points_grid_by_01(grid_path, N_total_x, N_regular_x, xBC, hx0, qx, Ny, yB, yC):
def generate_dop_points_grid_by_10(grid_path, N_total_y, N_regular_y, yAB, hy0, qy, Nx, xA, xB):
def generate_dop_points_grid_by_11(grid_path, N_total_y, N_regular_y, yDC, hy0, qy, Nx, xD, xC):
def generate_dop_points_grid_by_points_00(grid_path, N_total_x, N_regular_x, hx0, qx, Ny, points_x_00, points_y_00):
def generate_dop_points_grid_by_points_01(grid_path, N_total_x, N_regular_x, hx0, qx, Ny, points_x_01, points_y_01):
