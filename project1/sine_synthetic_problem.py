
import numpy as np
import matplotlib.pyplot as plt

def generate_initial_population(population_size):
    population = []
    for _ in range(population_size):
        bitstring = np.random.randint(2, size=bitstring_size)
        population.append(bitstring)
    return population

def bitstring_to_float(bitstring):
    return bitstring.dot(2 ** np.arange(bitstring_size)[::-1]) * scaling_factor

def fitness_function(bitstring):
    x = bitstring_to_float(bitstring)
    if constraint and (x < 5 or x > 10):
        if x < 5:
            return np.clip(np.sin(x) - (5 - x), -1, 1)
        return np.clip(np.sin(x) - (x - 10), -1, 1)
    return np.sin(x)

def find_parents(population, num_parents=2): 
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

def survivor_selection_with_crowding(population, population_size):
    fitnesses = [fitness_function(bitstring) for bitstring in population]
    fitnesses = np.array(fitnesses)

    sigma = np.pi / 2
    a = 1
    def sh(d): return 1 - (d / sigma) ** a if d < sigma else 0
    share_fitnesses = np.zeros(len(population))
    for i in range(len(population)):
        x = bitstring_to_float(population[i])
        for j in range(len(population)):
            y = bitstring_to_float(population[j])
            share_fitnesses[i] += sh(np.abs(x - y))

    fitnesses = fitnesses / share_fitnesses

    min_fitness = []
    for _ in range(len(population) - population_size):
        min_fitness.append(np.argmin(fitnesses))
        fitnesses[min_fitness[-1]] = np.inf
    population = np.delete(population, min_fitness, axis=0)
    return population

def calculate_entropy(population):
    bit_string_probability = np.sum(population, axis=0) / len(population)
    bit_string_probability = bit_string_probability[bit_string_probability != 0]
    return -np.sum(bit_string_probability * np.log2(bit_string_probability))    

global entropy_list
entropy_list = []

fig, axs = plt.subplots(1, 2, figsize=(10, 5), gridspec_kw={'width_ratios': [2, 1]})

def plot_population(population):
    # First subplot
    axs[0].cla()  # clear first subplot
    axs[0].plot(np.linspace(0, 128, 1000), np.sin(np.linspace(0, 128, 1000)))
    x = [bitstring_to_float(bitstring) for bitstring in population]
    y = [fitness_function(bitstring) for bitstring in population]
    axs[0].scatter(x, y, c='r')
    axs[0].set_title("Generation " + str(len(entropy_list)+1))
    axs[0].set_xlabel("Float value of individuals")
    axs[0].set_ylabel("Fitness")

    # Second subplot
    axs[1].cla()  # clear second subplot
    entropy_list.append(calculate_entropy(population))
    axs[1].plot(entropy_list, 'b-')  # plot num_list as blue line
    axs[1].set_title("Entropy")
    axs[1].set_ylim([0, 10])
    axs[1].set_xlim([0, num_generations])
    axs[1].set_xlabel("Generation")
    axs[1].set_ylabel("Entropy")

    plt.draw()
    plt.pause(0.001)

def genetic_algorithm(population_size, num_generations):
    population = generate_initial_population(population_size)
    global initial_population
    initial_population = population.copy()
    for _ in range(num_generations):
        plot_population(population)
        parents = find_parents(population)
        children = crossover(parents)
        children = [mutate(child) for child in children]
        population = np.concatenate((population, children))
        if crowding:
            population = survivor_selection_with_crowding(population, population_size)
        else:
            population = survivor_selection(population, population_size)
        # plt.show(block=False)
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

    bitstring_size = 15
    scaling_factor = 2 ** (7 - bitstring_size)

    constraint = False
    crowding = False

    population = genetic_algorithm(population_size, num_generations)
    best_solution(population)
    
