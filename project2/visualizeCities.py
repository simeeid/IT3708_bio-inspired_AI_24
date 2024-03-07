import matplotlib.pyplot as plt
import json

# will only visualize depot and patients, no solution

# Load your data
with open('src/main/resources/train/train_5.json', 'r') as f:
    data = json.load(f)

# Extract depot coordinates
depot_x = data['depot']['x_coord']
depot_y = data['depot']['y_coord']

# Extract city coordinates
city_x = [data['patients'][str(i)]['x_coord'] for i in data['patients'].keys()]
city_y = [data['patients'][str(i)]['y_coord'] for i in data['patients'].keys()]

# Create plot
plt.figure(figsize=(10,10))
plt.scatter(city_x, city_y, c='b', label='Cities')
plt.scatter(depot_x, depot_y, c='r', label='Depot')
plt.title('Depot and Cities')
plt.xlabel('X')
plt.ylabel('Y')
plt.legend()
plt.show()
