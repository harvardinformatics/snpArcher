import csv
import gzip
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

x_data = []
y_data = []

with gzip.open('.test/ccgp/results/genome1/CCGP/test_qc.stat.gz', 'rt') as file:
    # Skip the first line
    next(file)
    reader = csv.reader(file, delimiter='\t')
    for row in reader:
        x_data.append(float(row[0]))
        y_data.append(float(row[1]))

#print(x_data)
#print(y_data)

# Define the exponential decay function
def decay_func(x, a, b, c):
    return a * np.exp(-b * x) + c

# Fit the decay function to the data
popt, pcov = curve_fit(decay_func, x_data, y_data)

# Extract the decay constant (b) from the fitted parameters
decay_constant = popt[1]

# Calculate the half-life
half_life = np.log(2) / decay_constant

# Generate x values for the decay curve
x_fit = np.linspace(0, 6, 100)

# Compute the decay curve
y_fit = decay_func(x_fit, *popt)

# Plot the data and the decay curve
plt.scatter(x_data, y_data, label='Data')
plt.plot(x_fit, y_fit, 'r-', label='Decay curve')

# Mark the half-life point on the graph
plt.axhline(y=0.5, color='g', linestyle='--', label='Half-life')
plt.axvline(x=half_life, color='g', linestyle='--')

# Add labels and a legend
plt.xlabel('Time')
plt.ylabel('Intensity')
plt.legend()

# Show the plot
#plt.show()

# Print the half-life
print("Half-life:", half_life)
