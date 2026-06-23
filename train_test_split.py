from pathlib import Path
import random
import json

p = Path.cwd()

num = '2nd'

path2fold = p / f'{num}_sel_comp_npy'

train_files = []
test_files = []

random.seed(42)
train_ratio = 0.8

train_files = []
test_files = []

all_files = list([str(txt) for txt in path2fold.iterdir() if txt.is_file()])
random.shuffle(all_files)

split_idx = int(len(all_files) * train_ratio)
train_files = all_files[:split_idx]
test_files = all_files[split_idx:]

to_json = {'train': train_files,
           'test': test_files}


with open(f'{num}_sel.json', 'w', encoding='utf-8') as f:
    json.dump(to_json, f, ensure_ascii=False, indent=4)
