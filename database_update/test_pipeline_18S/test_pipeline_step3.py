import os
import pandas as pd

input_dir = "blast_results"
output_dir = "best_hits"
os.makedirs(output_dir, exist_ok=True)

for file in os.listdir(input_dir):
    if file.endswith(".txt"):
        input_path = os.path.join(input_dir, file)
    output_path = os.path.join(output_dir, f"{file.replace('.txt', '_top5.txt')}")
    with open(input_path, "r") as f:
        lines = [line.strip().split("\t") for line in f]
    # сортировка по e-value и взятие 5 первых
    lines_sorted = sorted(lines, key=lambda x: float(x[10]))
    top5 = lines_sorted[:5]
    with open(output_path, "w") as f:
        for line in top5:
            f.write("\t".join(line) + "\n")
