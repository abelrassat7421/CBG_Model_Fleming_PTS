"""Generates the seeds for the simulations and the position of 
the STN neurons as defined in generate_stn_xy_pos_new.py"""

import random
import numpy as np

random.seed(42)

num_seeds = 19
min_seed = 1
max_seed = 10000

seeds = [random.randint(min_seed, max_seed) for _ in range(num_seeds)]

# Convert the lists to comma-separated strings
seeds = ",".join(str(i) for i in seeds) 

print("Hello")
with open("Cortex_BasalGanglia_DBS_model/random_seeds.txt", "w+") as f:
    f.write(seeds)

seeds_read = np.loadtxt("Cortex_BasalGanglia_DBS_model/random_seeds.txt", delimiter=",")

print(seeds_read)
print("Hello")