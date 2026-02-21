import matplotlib.pyplot as plt

# Data
label = [r"$M_{bc}^{tag}$", "electronID", "muonID"]
result = [1.136825, 1.225300, 1.227647]
uncertainty = [0.028931, 0.035814, 0.032830]

weights = []
weight_times_values = []
for value, sigma in zip(result, uncertainty):
    weight = 1/(sigma * sigma)
    weights.append(weight)
    weight_times_values.append(weight * value)
    
weighted_average = sum(weight_times_values) / sum(weights)
width = 0.06

# Plot
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot()

# Plotting results with uncertainties
plt.errorbar(label, result, yerr=uncertainty, fmt='o', color='black', label='')

# Adding labels and title
plt.xlabel('type of sideband', fontsize=14)
plt.ylabel('data/MC', fontsize=14)
plt.title(("%.2f" % weighted_average) + r"$\pm$" + ("%.2f" % width), fontsize=16)

# draw correction factor
plt.axhline(y=weighted_average, color='r', linestyle='-', label='')
#plt.fill_between(label, 1.1 - 0.1, 1.1 + 0.1, color='green', alpha=0.3, label='Band around y=1.1')
ax.axhspan(weighted_average - width, weighted_average + width, facecolor=(0.0294117647, 1, 0), alpha=0.5)
ax.axhspan(weighted_average + width * 1, weighted_average + width * 2, facecolor=(1, 0.93333333333, 0), alpha=0.5)
ax.axhspan(weighted_average - width * 1, weighted_average - width * 2, facecolor=(1, 0.93333333333, 0), alpha=0.5)

# Adding legend
#plt.legend()

# Show plot
plt.grid(True, linestyle='--')
plt.tight_layout()
plt.show()
