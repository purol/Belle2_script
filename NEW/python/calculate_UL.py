import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
from scipy.interpolate import interp1d

# Step 1: Load the CSV file
df = pd.read_csv("total_scan_destroyer.csv", header=None)
x = df.iloc[:, 0]
y = np.exp(- df.iloc[:, 2])

# Step 2: Sort x and y based on x values
sorted_indices = np.argsort(x)
x_sorted = x.iloc[sorted_indices].values
y_sorted = y.iloc[sorted_indices].values

# Step 3: Filter x and y to only include points where x >= 0
x_non_negative = x_sorted[x_sorted >= 0]
y_non_negative = y_sorted[x_sorted >= 0]

# Step 4: Plot the data with the horizontal line y = exp(-0.5)
target_y = np.exp(-0.5)
plt.plot(x_sorted, y_sorted)
plt.axhline(target_y, color='red', linestyle='--', label=f"y = exp(-0.5)")
plt.xlabel(r"$\mu$")
plt.ylabel("profile likelihood")
plt.title("Profile likelihood scan")
plt.legend()
plt.show()

# Step 5: Calculate the total area under the curve from x = 0 to the end
total_area = simps(y_non_negative, x_non_negative)

# Step 6: Calculate the cumulative area at each point, starting from x = 0
cumulative_area = np.array([simps(y_non_negative[:i+1], x_non_negative[:i+1]) for i in range(len(x_non_negative))])

# Step 7: Interpolate to find the precise `c` where cumulative area reaches 90% of the total area
target_area = 0.9 * total_area
interp_func = interp1d(cumulative_area, x_non_negative)
c_value = float(interp_func(target_area))

# Step 8: Find the x value where y is highest within the non-negative x range
max_y_index = np.argmax(y_sorted)
x_at_max_y = x_sorted[max_y_index]
y_max_value = y_sorted[max_y_index]

# Step 9: Find all x values where y intersects with y = exp(-0.5) within non-negative x values
x_intersections = []
if y_sorted.min() <= target_y <= y_sorted.max():
    intersect_indices = np.where(np.diff(np.sign(y_sorted - target_y)))[0]
    for i in intersect_indices:
        # Interpolate to find exact x where the intersection occurs
        interp = interp1d([y_sorted[i], y_sorted[i + 1]], [x_sorted[i], x_sorted[i + 1]])
        x_intersections.append(float(interp(target_y)))
else:
    print(f"No intersections found for y = {target_y} within the data range.")

print(f"The total area under the curve from x=0: {total_area}")
print(f"The value of c such that the area from 0 to c is 90% of the total area: {c_value}")
print(f"The x value where y is highest: {x_at_max_y} (y = {y_max_value})")
print(f"The x values where y intersects with y = exp(-0.5): {x_intersections}")
