import numpy as np
import matplotlib.pyplot as plt

# Enable LaTeX rendering
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Times"],
    "axes.labelsize": 14,
    "font.size": 12,
    "legend.fontsize": 14,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
})

# Load data (assuming whitespace-delimited file)
data = np.loadtxt("LB_Results/Data_particle/data_particle_1.dat")

# Extract required columns
# Column indices in Python are 0-based: col1 = [:,0], col4 = [:,3]
nu_LB =0.05
Ellipse_Dia = 50
Ny = 200
x = data[:, 0] * nu_LB / (Ellipse_Dia**2)
y = data[:, 4] / Ny

xRange = 18

# Create plot
fig, ax = plt.subplots(figsize=(7, 5))
ax.set_xlim(0, xRange) 
ax.set_ylim(0.4, 0.7) 

ax.set_xticks(range(0, xRange + 1, 2))   # ticks at 0, 2, 4, ..., 24


ax.plot(x, y, label=r"Y vs T",
        color='blue', marker='', markersize=3, linestyle='-')

#### Xu (2018) data
# ----- Load second dataset -----
data2 = np.loadtxt("Xu2018_Y_vs_t.txt")

# Columns: col1 (index 0) vs col2 (index 1)
x2 = data2[:, 0]
y2 = data2[:, 1]

# ----- Plot on same axes using symbols -----
ax.plot(
    x2, y2,
    linestyle='',        # no line
    marker='o',          # choose symbol
    markersize=5,
    color='red',
    label=r"Xu data"
)


#############################

ax.set_xlabel(r"Time [s]")
ax.set_ylabel(r"Y pos [-]")
ax.set_title(r"\textbf{Data Plot}")
# Enable major & minor ticks and grids
ax.minorticks_on()
ax.grid(which='major', linestyle='-', linewidth=0.8, color='gray')
ax.grid(which='minor', linestyle='--', linewidth=0.5, color='lightgray')
ax.legend()

# Save the figure
plt.tight_layout()
fig.savefig("data_plot.pdf", format="pdf")

plt.show()
