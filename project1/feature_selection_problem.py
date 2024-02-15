
import numpy as np
import matplotlib.pyplot as plt
import LinReg
import pandas as pd

data = pd.read_csv('dataset.txt', header=None)
regressor = LinReg.LinReg()

def generate_initial_population(population_size):
    population = []
    for _ in range(population_size):
        bitstring = np.random.randint(2, size=bitstring_size)
        population.append(bitstring)
    return population

def fitness_function(bitstring):
    # change the last bit of bitstring to 1
    bitstring[-1] = 1
    X = regressor.get_columns(data.values, bitstring)
    return regressor.get_fitness(X[:,:-1], X[:,-1])

def find_parents(population, fitnesses, num_parents=2): 
    parents = []
    for _ in range(num_parents):
        max_fitness_idx = np.argmin(fitnesses)
        parents.append(population[max_fitness_idx])
        fitnesses[max_fitness_idx] = np.inf
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
    original_fitnesses = fitnesses
    fitnesses = np.array(fitnesses)
    min_fitness = []
    for _ in range(len(population) - population_size):
        min_fitness.append(np.argmax(fitnesses))
        fitnesses[min_fitness[-1]] = -np.inf
    population = np.delete(population, min_fitness, axis=0)
    original_fitnesses = np.delete(original_fitnesses, min_fitness, axis=0)
    return population, original_fitnesses

def normalize_data(data):
    return (data - np.min(data)) / (np.max(data) - np.min(data))

def survivor_selection_with_crowding(population, population_size):
    fitnesses = [fitness_function(bitstring) for bitstring in population]
    fitnesses = np.array(fitnesses)

    # similarity is the number of 1's in the same position
    def similarity(x, y): return np.sum(np.where(x == y, x, 0)) # np.sum(x == y)
    share_fitnesses = np.zeros(len(population))
    i = -1
    for x in population:
        i += 1
        for y in population:
            share_fitnesses[i] += similarity(x, y)
    
    normalized_share_fitnesses = normalize_data(share_fitnesses)
    scaled_fitnesses = fitnesses * normalized_share_fitnesses

    min_fitness = []
    for _ in range(len(population) - population_size):
        min_fitness.append(np.argmax(scaled_fitnesses))
        scaled_fitnesses[min_fitness[-1]] = -np.inf
    population = np.delete(population, min_fitness, axis=0)
    fitnesses = np.delete(fitnesses, min_fitness, axis=0)
    return population, fitnesses

def calculate_entropy(population):
    bit_string_probability = np.sum(population, axis=0) / len(population)
    bit_string_probability = bit_string_probability[bit_string_probability != 0]
    return -np.sum(bit_string_probability * np.log2(bit_string_probability))    

global entropy_list
entropy_list = []  # list to store all entropy values
fig, axs = plt.subplots(1, 2, figsize=(10, 5), gridspec_kw={'width_ratios': [2, 1]})  # specify width ratios

def plot_population(population, fitnesses):
    # First subplot
    axs[0].cla()  # clear first subplot
    axs[0].set_ylim([0, 0.5]) #10])  # set y-axis limits
    axs[0].set_xlim([0, 100])  # set x-axis limits
    x = [np.sum(population, axis=1)]
    y = fitnesses
    axs[0].scatter(x, y, c='r')
    axs[0].set_title("Generation " + str(len(entropy_list)+1))
    axs[0].set_xlabel("Number of columns used in individuals")
    axs[0].set_ylabel("Fitness")

    # Second subplot
    axs[1].cla()  # clear second subplot
    entropy_list.append(calculate_entropy(population))  # append entropy value to list
    axs[1].plot(entropy_list, 'b-')  # plot entropy
    axs[1].set_title("Entropy")
    axs[1].set_ylim([0, 60]) #10])  # set y-axis limits
    axs[1].set_xlim([0, num_generations])  # set x-axis limits
    axs[1].set_xlabel("Generation")
    axs[1].set_ylabel("Entropy")

    plt.draw()
    plt.pause(0.001)

def genetic_algorithm():
    global_best_individual = [np.inf, [], 0]
    generation_num = 1
    population = generate_initial_population(population_size)
    fitnesses = [fitness_function(bitstring) for bitstring in population]
    for _ in range(num_generations):
        plot_population(population, fitnesses)
        min_fitness = np.min(fitnesses)
        if global_best and min_fitness < global_best_individual[0]:
            global_best_individual = [min_fitness, population[np.argmin(fitnesses)], generation_num]
        if early_stopping and min_fitness < threshold:
            break
        parents = find_parents(population, fitnesses)
        children = crossover(parents)
        children = [mutate(child) for child in children]
        population = np.concatenate((population, children))
        if crowding:
            population, fitnesses = survivor_selection_with_crowding(population, population_size)
        else:
            population, fitnesses = survivor_selection(population, population_size)
        generation_num += 1
    return population, fitnesses, generation_num, global_best_individual

def best_solution(population, fitnesses, generation_num):
    max_fitness_idx = np.argmin(fitnesses)
    print('Max fitness:', fitnesses[max_fitness_idx])
    print('Bitstring:', population[max_fitness_idx])
    print('Generation:', generation_num)

if __name__ == '__main__':
    population_size = 20
    num_generations = 50

    bitstring_size = data.shape[1]

    early_stopping = False
    threshold = 0.124
    global_best = True
    crowding = False

    population, fitnesses, generation_num, global_best_individual = genetic_algorithm()
    if global_best:
        best_solution([global_best_individual[1]], [global_best_individual[0]], global_best_individual[2])
    else:
        best_solution(population, fitnesses, generation_num)
    
