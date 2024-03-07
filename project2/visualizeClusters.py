import json
import matplotlib.pyplot as plt

# Load clusters from JSON file
with open('clusters.json', 'r') as f:
    clusters = json.load(f)

# Create a new plot
fig, ax = plt.subplots()

# Plot each cluster
for i, cluster in enumerate(clusters):
    # Convert cluster to list of points
    points = [(patient['x_coord'], patient['y_coord']) for patient in cluster]
    xs, ys = zip(*points)

    # Plot points
    ax.scatter(xs, ys, label=f'Cluster {i+1}')

# Show the legend and display the plot
ax.legend()
plt.show()
