#By Alena Favorskaya, Jul 2025

Nt = 601
Nx = 401
Ny = 201
step = 1.0
cP = 2.0 #это км в сек. Или м в мс. то есть раз у нас все в м, то время в мс измерятеся. и тогда скорости звука такие
cS = 1.0
rho = 2.0
dt = 0.4 #это мс
N_chi = 10
save = 6

receivers_number_half = 200
receivers_start_x = (0-receivers_number_half) * step
receivers_end_x = receivers_number_half * step
receivers_step = step
receivers_y = 0.0 - step
sei_save = 1
receivers_number = receivers_number_half * 2 + 1

#наша трещина -- это отрезок AB. Причем точка А -- левая и (или) нижняя. Если трещина близка к горизонатали - приоритетно левая. Если трещина близка к вертикали - приоритетно нижняя. Можно этот выбор точек автоматизировать.
"""
var = 2

if (var == 0):
	pos_A_x = -25.0
	pos_A_y = -71.0
	pos_B_x = 25.0
	pos_B_y = -40.0
elif (var == 1):
	pos_A_x = -25.0
	pos_A_y = -50.0
	pos_B_x = 20.0
	pos_B_y = -70.0
elif (var == 2):
	pos_A_x = -10.0
	pos_A_y = -100.0
	pos_B_x = 10.0
	pos_B_y = -40.0
elif (var == 3):
	pos_A_x = 10.0
	pos_A_y = -100.0
	pos_B_x = -10.0
	pos_B_y = -40.0
"""
freq = 0.03
