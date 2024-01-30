
import numpy as np
import matplotlib.pyplot as plt

bitstring_size = 15
scaling_factor = 2 ** (7 - bitstring_size)

num_parents = 2
constraint = True

def generate_initial_population(population_size):
    population = []
    for _ in range(population_size):
        bitstring = np.random.randint(2, size=bitstring_size)
        population.append(bitstring)
    return population

def bitstring_to_float(bitstring):
    return bitstring.dot(2 ** np.arange(bitstring_size)[::-1]) * scaling_factor

# def bitstring_to_float(bitstring):
#     return 2 ** sum(bitstring) * scaling_factor

def fitness_function(bitstring):
    x = bitstring_to_float(bitstring)
    if constraint and (x < 5 or x > 10):
        return np.clip(np.sin(x) - 0.5, -1, np.inf)
    return np.sin(x)

def find_parents(population, num_parents): 
    fitnesses = [fitness_function(bitstring) for bitstring in population]
    parents = []
    for _ in range(num_parents):
        max_fitness_idx = np.argmax(fitnesses)
        parents.append(population[max_fitness_idx])
        fitnesses[max_fitness_idx] = -np.inf
    return parents

def mutate(bitstring):
    for i in range(bitstring_size):
        if np.random.rand() < 1/bitstring_size:
            bitstring[i] = 1 - bitstring[i]
    return bitstring

def crossover(parents):
    parent1, parent2 = parents
    crossover_point = np.random.randint(1, bitstring_size)
    child1 = np.concatenate((parent1[:crossover_point], parent2[crossover_point:]))
    child2 = np.concatenate((parent2[:crossover_point], parent1[crossover_point:]))
    return child1, child2

def survivor_selection(population, population_size):
    fitnesses = [fitness_function(bitstring) for bitstring in population]
    fitnesses = np.array(fitnesses)
    min_fitness = []
    for _ in range(len(population) - population_size):
        min_fitness.append(np.argmin(fitnesses))
        fitnesses[min_fitness[-1]] = np.inf
    population = np.delete(population, min_fitness, axis=0)
    return population

global num
num = 1
plt.ion()
def plot_population(population):
    global num
    plt.clf()
    x = np.linspace(0, 128, 1000)
    y = np.sin(x)
    plt.plot(x, y)
    x = [bitstring_to_float(bitstring) for bitstring in population]
    y = [fitness_function(bitstring) for bitstring in population]
    plt.scatter(x, y, c='r')
    # plt.show()
    plt.title("Generation " + str(num))
    plt.draw()
    plt.pause(0.1)
    num += 1
plt.show(block=True)

def genetic_algorithm(population_size, num_generations):
    population = generate_initial_population(population_size)
    global initial_population
    initial_population = population.copy()
    for _ in range(num_generations):
        plot_population(population)
        parents = find_parents(population, num_parents)
        children = crossover(parents)
        children = [mutate(child) for child in children]
        population = np.concatenate((population, children))
        population = survivor_selection(population, population_size)
    return population

def best_solution(population):
    fitnesses = [fitness_function(bitstring) for bitstring in population]
    if constraint:
        feasible = False
        for i in range(len(population)):
            if 5 <= bitstring_to_float(population[i]) <= 10:
                feasible = True
            else:
                fitnesses[i] = -np.inf
        if feasible:
            max_fitness_idx = np.argmax(fitnesses)
            print('Max fitness:', fitnesses[max_fitness_idx])
            print('Bitstring:', population[max_fitness_idx])
            print('Float:', bitstring_to_float(population[max_fitness_idx]))
        else:
            print('No feasible solutions')
    else:
        max_fitness_idx = np.argmax(fitnesses)
        print('Max fitness:', fitnesses[max_fitness_idx])
        print('Bitstring:', population[max_fitness_idx])
        print('Float:', bitstring_to_float(population[max_fitness_idx]))

if __name__ == '__main__':
    population_size = 100
    num_generations = 100
    population = genetic_algorithm(population_size, num_generations)
    best_solution(population)
    
