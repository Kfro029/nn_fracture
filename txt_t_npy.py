import numpy as np
import os
import io

input_dir = "1st_sel_comp"
output_dir = "1st_sel_comp_npy"

os.makedirs(output_dir, exist_ok=True)

files = [f for f in os.listdir(input_dir) if f.endswith(".txt")]

for fname in files:
    path = os.path.join(input_dir, fname)
    with open(path) as f:
        lines = [line.rstrip('\n').rstrip(';') for line in f][1:]
    data = np.loadtxt(io.StringIO('\n'.join(lines)), delimiter=';')
    values = data[:, 1:]
    vx = values[:, 0::2]
    vy = values[:, 1::2]
    x = np.stack([vx, vy], axis=0).astype(np.float32)
    out_path = os.path.join(output_dir, fname.replace(".txt", ".npy"))
    np.save(out_path, x)
