import os
import matplotlib.pyplot as plt
import numpy as np

def parse_filename(filename):
    """
    Извлекает количество трещин и угол из имени файла через split('_').
    Ожидаемый формат: rec_frac_<кол-во>_<угол>.txt
    """
    parts = filename.split('_')
    # parts = ['rec', 'frac', '<кол-во>', '<угол>.txt']
    count = int(parts[-9])
    angle = float(parts[-2].replace('.npy', ''))
    return count, angle


directory = os.path.join('', '2nd_selection')
counts = []
angles = []

for filename in os.listdir(directory):
    if filename.endswith('.npy'):
        # Предполагаем, что все .txt файлы имеют нужный формат
        cnt, ang = parse_filename(filename)
        counts.append(cnt)
        angles.append(ang)


print(angles)

# Построение двух бар-чартов
plt.figure(figsize=(12, 5))

# 1. Распределение по количеству трещин
plt.subplot(1, 2, 1)
unique_counts, freq_counts = np.unique(counts, return_counts=True)
bars = plt.bar(unique_counts, freq_counts)
plt.bar_label(bars)  # автоматически подписывает значения высоты
plt.bar(unique_counts, freq_counts, color='skyblue', edgecolor='black')
plt.xlabel('Количество трещин')
plt.ylabel('Частота')
plt.title('Распределение количества трещин')
plt.grid(axis='y', linestyle='--', alpha=0.7)



# 2. Распределение по углам наклона
plt.subplot(1, 2, 2)
unique_angles, freq_angles = np.unique(angles, return_counts=True)
bars = plt.bar(unique_angles, freq_angles)
plt.bar_label(bars)  # автоматически подписывает значения высоты
plt.bar(unique_angles, freq_angles, color='salmon', edgecolor='black')
plt.xlabel('Угол наклона (градусы)')
plt.ylabel('Частота')
plt.title('Распределение углов наклона')
plt.grid(axis='y', linestyle='--', alpha=0.7)

plt.tight_layout()
plt.savefig('cnt_ang.png')

plt.show()


