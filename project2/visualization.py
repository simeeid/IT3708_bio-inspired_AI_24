import json
import matplotlib.pyplot as plt
import numpy as np

# Will visualize the depot, patients and a solution

# Load data from json file
with open('file2.json', 'r') as f:
    data = json.load(f)

# Extract the data and the hub
people = json.loads(data['data'])
hub = json.loads(data['hub'])
travel_time = json.loads(data['travel_time'])

# Count the total number of travelers and the number of travelers not visiting any cities
total_travelers = len(people)
travelers_not_visiting = sum(1 for person in people if not person)

# Create a color map
colors = plt.cm.tab20(np.linspace(0, 1, 20))

# Plot each person's path
# for i, person in enumerate(people):
#     if person:  # Check if the person array is not empty
#         person = np.array(person)
#         x = np.concatenate(([hub[1]], person[:, 1], [hub[1]]))
#         y = np.concatenate(([hub[2]], person[:, 2], [hub[2]]))
#         color = colors[i % 20]  # Cycle through the colors
#         plt.scatter(x, y, color=color)  # plot points
#         plt.plot(x, y, color=color)  # plot lines
# Plot each person's path
for i, person in enumerate(people):
    if person:  # Check if the person array is not empty
        person = np.array(person)
        x = np.concatenate(([hub[1]], person[:, 1], [hub[1]]))
        y = np.concatenate(([hub[2]], person[:, 2], [hub[2]]))
        color = colors[i % 20]  # Cycle through the colors
        linewidth = 1 + (i % 5) * 0.2  # Cycle through the line widths with a step of 0.2
        plt.scatter(x, y, color=color)  # plot points
        plt.plot(x, y, color=color, linewidth=linewidth)  # plot lines

# Plot the hub on top
plt.scatter(hub[1], hub[2], color='black', zorder=5)

# Add a limited legend
plt.plot([], [], color='black', label='Hub')
for i in range(min(3, total_travelers)):  # Adjust this value to include more or fewer travelers in the legend
    plt.plot([], [], color=colors[i % 20], label=f'Person {i+1}')
if total_travelers > 3:
    plt.plot([], [], color='gray', label='...')

plt.title(f'Total travelers: {total_travelers}\nTravelers not visiting any cities: {travelers_not_visiting}\nTotal travel time: {travel_time}')
plt.legend()
plt.show()
