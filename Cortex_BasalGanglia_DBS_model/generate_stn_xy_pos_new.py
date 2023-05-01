import random
import numpy as np
import matplotlib.pyplot as plt

# Define the range of x and y values
y_range = [-2000, 2000]
x_range_left = [-2000, -500]
x_range_right = [500, 2000]
x_range_full = [-2000, 2000]
num_positions = 100
electrode_depth = -1550

# randomly sampling y coordinates 
y_points = np.random.uniform(y_range[0], y_range[1], size=100)
x_points = []

# decides whether neurons that are not "under" the electrode are on the left or the right (true and false respectively) 
values = [True, False]
electrode_side = [random.choice(values) for _ in range(num_positions)]

for i in range(num_positions):
   if y_points[i] > electrode_depth: 
      if electrode_side[i]:
        x_points.append(random.uniform(x_range_left[0], x_range_left[1]))
      else:
        x_points.append(random.uniform(x_range_right[0], x_range_right[1]))
   else: 
      x_points.append(random.uniform(x_range_full[0], x_range_full[1]))

y_points = y_points.tolist()

# Convert the lists to comma-separated strings
x_points = ",".join(str(i) for i in x_points) 
x_points += "\n"
y_points = ",".join(str(i) for i in y_points)

with open("Cortex_BasalGanglia_DBS_model/STN_xy_pos_new.txt", "w+") as f:
    f.write(x_points)
    f.write(y_points)

#print(x_points)




      


