#Author:  Alena V. Favorskaya 07 Jun 22

#запуск

#python3 job_maker.py run.sh job_name folder_name job_file.job

#run.sh - имя скрипта, который будем переделывать в job_file.job

#folder_name - имя папки на сервере, в которую все нужно будет скопировать
#scp -r folder_name_from aleanera@geophys.mipt.ru:~/folder_name

import sys

fw = open(sys.argv[4], 'wt')

fw.write("#!/bin/bash\n")
fw.write("\n")
fw.write("#PBS -l walltime=1000:01:00,nodes=1:ppn=1\n")
fw.write("#PBS -N " + sys.argv[2] + "\n")
fw.write("#PBS -q batch\n")
fw.write("\n")
fw.write("cd $PBS_O_WORKDIR\n")
fw.write("cd ~/" + sys.argv[3] + "\n")
fw.write("\n")

f = open(sys.argv[1], 'rt')
for l in f:

	change_flag = False
	for c in l:
		if (c == '	'):
			continue
		if (c == '/'): # /home/aleanera/work/rect/build/rect >> ~/rect/build/rect
			change_flag = True
			break
		break

	if change_flag:
		flag = 0
		new_line = ""
		for c in l:
			if (flag == 0): #before
				if (c == '/'):
					flag = 1
					continue
				new_line = new_line + c
				continue
			if (flag == 1): #skipping
				if (c == 'k'):
					flag = 2
				continue
			if (flag == 2): #changing
				flag = 3
				new_line = new_line + '~'
				continue
			#if (flag == 3): #after
			new_line = new_line + c
		fw.write(new_line)
	else:
		fw.write(l)
		
fw.close()
f.close()



