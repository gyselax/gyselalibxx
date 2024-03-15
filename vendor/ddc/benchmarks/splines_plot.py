# First execute :
# ./benchmarks/ddc_benchmark_splines --benchmark_format=json --benchmark_out=splines_bench.json
# then execute this code will be able to plot results:
# python3 splines_plot.py /path/to/splines_bench.json

import argparse
import matplotlib.pyplot as plt
import json
import numpy as np

parser = argparse.ArgumentParser(description="Plot bytes_per_second from a JSON file.")
parser.add_argument("json_file", help="Path to the JSON file")
args = parser.parse_args()

with open(args.json_file, 'r') as file:
        data = json.load(file);

# Extract the values at the end of "name" and corresponding "bytes_per_second"
nx_values = sorted(set(int(benchmark["name"].split("/")[1]) for benchmark in data["benchmarks"]))
data_groups = {nx: {"ny": [], "cols_per_chunk": [], "preconditionner_max_block_size": [], "bytes_per_second": [], "gpu_mem_occupancy": []} for nx in nx_values}

for benchmark in data["benchmarks"]:
    nx = int(benchmark["name"].split("/")[1])
    data_groups[nx]["ny"].append(int(benchmark["name"].split("/")[2]))
    data_groups[nx]["cols_per_chunk"].append(int(benchmark["name"].split("/")[3]))
    data_groups[nx]["preconditionner_max_block_size"].append(int(benchmark["name"].split("/")[4]))
    data_groups[nx]["bytes_per_second"].append(benchmark["bytes_per_second"])
    data_groups[nx]["gpu_mem_occupancy"].append(benchmark["gpu_mem_occupancy"])

########
## ny ##
########

# Plotting the data for each group
plt.figure(figsize=(8, 6))
for nx, group_data in data_groups.items():
    ny = group_data["ny"]
    throughput = [group_data["bytes_per_second"][i] for i in range(len(ny))]
    plt.plot(ny, throughput, marker='o', markersize=5, label=f'nx={nx}')

x = np.linspace(min(ny), 20*min(ny))
plt.plot(x, np.mean([data_groups[nx]["bytes_per_second"][0] for nx in nx_values])/min(ny)*x, linestyle='--', color='black', label='perfect scaling')

# Plotting the data
plt.grid()
plt.xscale("log")
plt.xlabel("ny")
plt.ylabel("Throughput [B/s]")
plt.title("Throughput on "+str.upper(data["context"]["chip"]));
plt.legend()
plt.savefig("throughput_ny.png")

#gpu_mem
plt.figure(figsize=(8, 6))
for nx, group_data in data_groups.items():
    ny = [group_data["ny"][i] for i in range(len(group_data["ny"])) if group_data["ny"][i]>=8e3]
    gpu_mem_overhead = [(group_data["gpu_mem_occupancy"][i]-nx*group_data["ny"][i]*8)/(nx*group_data["ny"][i]*8)*100 for i in range(len(group_data["ny"])) if group_data["ny"][i]>=8e3]
    plt.plot(ny, gpu_mem_overhead, marker='o', markersize=5, label=f'nx={nx}')

# Plotting the data
plt.grid()
plt.xscale("log")
plt.xlabel("ny")
plt.ylabel("Relative memory overhead [%]")
plt.title("Relative memory occupancy overhead")
plt.legend()
plt.savefig("gpu_mem_occupancy.png")

########################
## cols_per_chunk ##
########################

# Plotting the data for each group
plt.figure(figsize=(8, 6))
for nx, group_data in data_groups.items():
    cols_per_chunk = group_data["cols_per_chunk"]
    throughput = [group_data["bytes_per_second"][i] for i in range(len(cols_per_chunk))]
    plt.plot(cols_per_chunk, throughput, marker='o', markersize=5, label=f'nx={nx}')

x = [(int)(data["context"]["cols_per_chunk_ref"]), (int)(data["context"]["cols_per_chunk_ref"])*1.001];
plt.plot(x, [0.99*min([min(group_data["bytes_per_second"]) for nx, group_data in data_groups.items()]), 1.01*max([max(group_data["bytes_per_second"]) for nx, group_data in data_groups.items()])], linestyle='dotted', color='black', label='reference config')

# Plotting the data
plt.grid()
plt.xscale("log")
plt.xlabel("cols_per_chunk")
plt.ylabel("Throughput [B/s]")
plt.title("Throughput on "+str.upper(data["context"]["chip"])+" (with ny=100000)");
plt.legend()
plt.savefig("throughput_cols.png")

#####################
## preconditionner ##
#####################

# Plotting the data for each group
plt.figure(figsize=(8, 6))
for nx, group_data in data_groups.items():
    preconditionner_max_block_size = group_data["preconditionner_max_block_size"]
    throughput = [group_data["bytes_per_second"][i] for i in range(len(preconditionner_max_block_size))]
    plt.plot(preconditionner_max_block_size, throughput, marker='o', markersize=5, label=f'nx={nx}')

x = [(int)(data["context"]["preconditionner_max_block_size_ref"]), (int)(data["context"]["preconditionner_max_block_size_ref"])*1.001];
plt.plot(x, [0.99*min([min(group_data["bytes_per_second"]) for nx, group_data in data_groups.items()]), 1.01*max([max(group_data["bytes_per_second"]) for nx, group_data in data_groups.items()])], linestyle='dotted', color='black', label='reference config')

# Plotting the data
plt.grid()
plt.xscale("log")
plt.xlabel("preconditionner_max_block_size")
plt.ylabel("Throughput [B/s]")
plt.title("Throughput on "+str.upper(data["context"]["chip"])+" (with ny=100000)");
plt.legend()
plt.savefig("throughput_precond.png")

plt.close();
