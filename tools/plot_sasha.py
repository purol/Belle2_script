import matplotlib.pyplot as plt

# Data
FOM = [0.812066844, 0.821596005, 0.794924485, 0.765978086, 0.67447189, 0.498989698]
Cut = [r"$M_{bc}^{tag}>5.27 \, \mathrm{GeV}$", r"$M_{bc}^{tag}>5.272 \, \mathrm{GeV}$", 
       r"$M_{bc}^{tag}>5.274 \, \mathrm{GeV}$", r"$M_{bc}^{tag}>5.276 \, \mathrm{GeV}$", 
       r"$M_{bc}^{tag}>5.278 \, \mathrm{GeV}$", r"$M_{bc}^{tag}>5.28 \, \mathrm{GeV}$"]

# Plot
plt.figure(figsize=(10, 6))
plt.plot(Cut, FOM, marker='o', color='b', linestyle='-')

# Adding labels and title
plt.xlabel('Cut', fontsize=14)
plt.ylabel('FOM', fontsize=14)
plt.title('Figure of Merit vs Cut', fontsize=16)

# Adding grid
plt.grid(True)

# Rotating x-axis labels for better readability
plt.xticks(rotation=45, ha='right')

# Show plot
plt.tight_layout()
plt.show()
