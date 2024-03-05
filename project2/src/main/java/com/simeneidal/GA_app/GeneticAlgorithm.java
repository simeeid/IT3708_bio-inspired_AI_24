package com.simeneidal.GA_app;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.ArrayList;
import java.util.Arrays;

import com.google.gson.Gson;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.reflect.Array;

public class GeneticAlgorithm {

    private int populationSize;
    private String filePath;

    private int nbrNurses;
    private int capacityNurse;
    private JsonData.Depot depot;
    private int nbrPatients;
    private Map<String, JsonData.Patient> patients;
    private double[][] travelTimes;
    private Individual[] population;

    public GeneticAlgorithm(int populationSize, String filePath, Individual[] localPopulation) {
        this.populationSize = populationSize;
        this.filePath = filePath;

        String instanceName = "";
        nbrNurses = 0;
        capacityNurse = 0;
        double benchmark = 0;
        depot = null;
        patients = null;
        travelTimes = null;
        nbrPatients = 0;

        try (FileReader reader = new FileReader(filePath)) {
            // Create Gson instance
            Gson gson = new Gson();

            // Deserialize JSON to Java object
            JsonData jsonData = gson.fromJson(reader, JsonData.class);

            // Access data from the Java object
            instanceName = jsonData.getInstanceName();
            nbrNurses = jsonData.getNbrNurses();
            capacityNurse = jsonData.getCapacityNurse();
            benchmark = jsonData.getBenchmark();
            depot = jsonData.getDepot();
            patients = jsonData.getPatients();
            travelTimes = jsonData.getTravelTimes(); 

            // find the number of patients
            nbrPatients = patients.size();
        } catch (Exception e) {
            e.printStackTrace();
        }
        
    }

    // method that will actually run the GA
    public void runGA(Individual[] localPopulation, int nbrGenerations) {
        // System.out.println("Running the GA");
        
        if (Objects.isNull(localPopulation)) {
            // generatePopulation();
            generateGreedyPopulation();
        } else {
            population = localPopulation;
        }

        Individual nextBestIndividual = population[0];
        Individual bestIndividual = population[0];
        boolean flag = false;
        boolean secondFlag = true;
        
        for (int i = 0; i < nbrGenerations; i++) {
            if (false && (i + 1) % 100000 == 0) {
                // find the next best individual
                // Individual bestIndividual = population[0];
                double bestFitness = bestIndividual.getFitness();
                double nextBestFitness = nextBestIndividual.getFitness();
                for (Individual individual : population) {
                    if (individual.getFitness() < bestFitness) {
                        nextBestIndividual = bestIndividual;
                        bestIndividual = individual;
                    } else if (individual.getFitness() < nextBestIndividual.getFitness()) {
                        nextBestIndividual = individual;
                    }
                }
                if (nextBestFitness == nextBestIndividual.getFitness()) {// && bestFitness == bestIndividual.getFitness()) {
                    if (flag) {
                        if (true) {//secondFlag) {
                            // replace the entire population except for the most fit individual, with random (greedy) individuals
                            //generatePopulation();
                            generateGreedyPopulation();
                            population[0] = bestIndividual;
                            System.out.println("Population reset");
                            flag = false;
                            secondFlag = false;
                        } else {
                            // make copies of the best individual and use destroy operators to create new individuals
                            cutLongestTravel(bestIndividual);
                            // population[populationSize - 1] = bestIndividual;
                            System.out.println("Population made from best individual");
                            flag = false;
                            secondFlag = true;
                        }
                    } else {
                        flag = true;
                        System.out.println("Flag set");
                    }
                }
            }

            Individual[] parents = parentSelection();
            if (Math.random() < 0.1) {
                Individual parent1 = new Individual(parents[0].getChromosome().clone(), 0);
                Individual parent2 = new Individual(parents[1].getChromosome().clone(), 0);

                mutation(parent1); mutation(parent2);
                calculateFitness(parent1); calculateFitness(parent2);
                survivorSelection(parent1, parent2);
            } else {
                Individual[] children = crossover(parents[0].getChromosome(), parents[1].getChromosome());
                mutation(children[0]); mutation(children[1]);
                calculateFitness(children[0]); calculateFitness(children[1]);
                survivorSelection(children[0], children[1]);
            }
        }
    }
    

    public void cutLongestTravel(Individual bestIndividual) {
        Individual[] parentCopies = new Individual[populationSize];
        int[] chromosome = bestIndividual.getChromosome().clone();
        parentCopies[0] = bestIndividual;
        for (int j = 1; j < populationSize; j++) {
            parentCopies[j] = new Individual(chromosome.clone(), 0);
        }

        int[] chromosomeCopy = chromosome.clone();
        for (int k = 0; k < chromosomeCopy.length; k++) {
            if (chromosomeCopy[k] < 0) {
                chromosomeCopy[k] = 0;
            }
        }

        int[] indexesBroken = new int[populationSize];
        for (int j = 1; j < populationSize; j++) {
            indexesBroken[j] = 0;
            double longestTravelTime = 0;
            for (int k = 0; k < chromosome.length - 1; k++) {
                if (chromosomeCopy[k] > 0 && chromosomeCopy[k+1] > 0 &&travelTimes[chromosomeCopy[k]][chromosomeCopy[k+1]] > longestTravelTime) {
                    // check that k is not in indexesBroken
                    boolean flag = false;
                    for (int l = 0; l < indexesBroken.length; l++) {
                        if (k + 1 == indexesBroken[l]) {
                            flag = true;
                            break;
                        }
                    }
                    if (!flag) {
                        longestTravelTime = travelTimes[chromosomeCopy[k]][chromosomeCopy[k+1]];
                        indexesBroken[j] = k+1;
                    }
                }
            }

            int[] tempArray = new int[chromosome.length];
            int pointer = indexesBroken[j];
            int count = 0;
            while (chromosome[pointer] > 0) {
                tempArray[count] = chromosome[pointer];
                count++;
                pointer++;
            }

            boolean spaceFlag = false;
            int previous = -1;
            int emptyIndex = 0;
            for (int k = 0; k < chromosome.length; k++) {
                int current = chromosome[k];
                if (current < 0 && previous < 0) {
                    spaceFlag = true;
                    break;
                }
                previous = current;
                emptyIndex = k;
            }
            if (spaceFlag) {
                // System.out.println("Empty index: " + emptyIndex + ", indexes broken: " + indexesBroken[j] + ", count: " + count);
                if (emptyIndex < indexesBroken[j]) {
                    for (int k = emptyIndex; k < emptyIndex + count; k++) {
                        parentCopies[j].getChromosome()[k] = tempArray[k - emptyIndex];
                    }
                    for (int k = emptyIndex + count; k < indexesBroken[j] + count; k++) {
                        parentCopies[j].getChromosome()[k] = chromosome[k - count];
                    }
                } else {
                    for (int k = indexesBroken[j]; k < indexesBroken[j] - count; k++) {
                        parentCopies[j].getChromosome()[k] = chromosome[k + count];
                    }
                    for (int k = emptyIndex + count; k < indexesBroken[j] - count; k++) {
                        parentCopies[j].getChromosome()[k] = tempArray[k + emptyIndex];
                    }
                } 
            }
        }

        for (int j = 1; j < populationSize; j++) {
            calculateFitness(parentCopies[j]);
        }

        if (true) {
            population = parentCopies; // replace the entire population with the copies

            // generateGreedyPopulation();
            // for (int j = 0; j < populationSize; j++) {
            //     if (j % 2 == 0) {
            //         population[j] = parentCopies[j/2];
            //     }
            // }
        } else {
            // find the worst individual in the population and replace this with the best individual in parentCopies
            Individual bestCopyIndividual = parentCopies[1];
            for (int i = 1; i < populationSize; i++) {
                if (parentCopies[i].getFitness() < bestCopyIndividual.getFitness()) {
                    bestCopyIndividual = parentCopies[i];
                }
            }

            int worstIndex = 0;
            double worstFitness = population[0].getFitness();
            for (int i = 0; i < populationSize; i++) {
                if (population[i].getFitness() > worstFitness) {
                    worstFitness = population[i].getFitness();
                    worstIndex = i;
                }
            }

            population[worstIndex] = bestCopyIndividual; // replace only the worst individual with the best copy
        }
    }

    public void generatePopulation() {
        population = new Individual[populationSize];
        for (int i = 0; i < populationSize; i++) {
            population[i] = new Individual(new int[nbrNurses + nbrPatients], 0);
        }

        int[] patientList = new int[nbrPatients];
        for (int i = 0; i < nbrPatients; i++) {
            patientList[i] = i + 1;
        }

        for (int i = 0; i < populationSize; i++) {
            // randomly shuffle the patient list
            List<Integer> list = new ArrayList<>();
            for (int patient : patientList) {
                list.add(patient);
            }
            // add nbr_nurses - 1 0s to the list
            for (int j = 0; j < nbrNurses - 1; j++) {
                list.add(-j - 1);
            }
            Collections.shuffle(list);
            list.add(-nbrNurses); // add a 0 to the end of the list
            population[i].setChromosome(list.stream().mapToInt(Integer::intValue).toArray());
        }

        for (Individual individual : population) {
            calculateFitness(individual);
        }
    }

    public void generateGreedyPopulation() {
        population = new Individual[populationSize];
        for (int i = 0; i < populationSize; i++) {
            population[i] = new Individual(new int[nbrNurses + nbrPatients], 0);
        }

        int[] patientList = new int[nbrPatients];
        for (int i = 0; i < nbrPatients; i++) {
            patientList[i] = i + 1;
        }

        for (int i = 0; i < populationSize; i++) {
            // randomly shuffle the patient list
            List<Integer> list = new ArrayList<>();
            for (int patient : patientList) {
                list.add(patient);
            }
            // add nbr_nurses - 1 0s to the list
            // for (int j = 0; j < nbrNurses - 1; j++) {
            //     // list.add(-j - 1);
            //     list.add(0);
            // }
            Collections.shuffle(list);
            // list.add(-nbrNurses); // add a 0 to the end of the list
            // population[i].setChromosome(list.stream().mapToInt(Integer::intValue).toArray());
            
            int[] temp = list.stream().mapToInt(Integer::intValue).toArray();
            int[] chromosome = new int[temp.length + 25];

            // add nbrNurses - 1 0s to the list
            for (int j = 0; j < nbrNurses - 1; j++) {
                chromosome[j] = 0;
            }

            // System.out.println(Arrays.toString(temp));

            //System.out.println("Temp length: " + temp.length);
            int count = nbrNurses - 1;
            int start = 0;
            // iterate though temp and place each int where they cause the least increase in fitness
            for (int j = 0; j < temp.length; j++) {
                count++;
                double minIncrease = 10000;
                int minIncreaseIndex = 0;
                for (int k = start; k < count; k++) {
                    //System.out.println("k: " + k + " j: " + j);
                    if (travelTimes[chromosome[k]][temp[j]] + travelTimes[chromosome[k+1]][temp[j]] < minIncrease) {
                        minIncrease = travelTimes[chromosome[k]][temp[j]] + travelTimes[chromosome[k+1]][temp[j]];
                        minIncreaseIndex = k;
                    }
                }

                for (int k = count; k > minIncreaseIndex; k--) {
                    chromosome[k] = chromosome[k-1];
                }

                chromosome[minIncreaseIndex] = temp[j];
                if (j == 0) {
                    start = 1;
                } else if ( 0 < j && j < 26) {
                    start = j * 2;
                } else {
                    start = 0;
                }
                // if (j % 2 == 0) {
                //     if (j == 0) {
                //         start = 1;
                //     } else if ( 0 < j && j < 26) {
                //         start = j * 2;
                //     } else {
                //         start = 0;
                //     }
                // }
            }
            // set the last element to 0
            // chromosome[chromosome.length - 1] = 0;

            int num = -1;
            for (int j = 0; j < chromosome.length; j++) {
                if (chromosome[j] == 0) {
                    chromosome[j] = num;
                    num--;
                }
            }

            population[i].setChromosome(chromosome);

            /* System.out.println(Arrays.toString(chromosome));
            System.out.println(chromosome.length); */

            // check that chromosome contains all numbers 1 - 100
            // for (int j = -25; j < 101; j++) {
            //     boolean flag = false;
            //     for (int k = 0; k < chromosome.length; k++) {
            //         if (chromosome[k] == j) {
            //             flag = true;
            //             break;
            //         }
            //     }
            //     if (!flag) {
            //         System.out.println("Number " + j + " is missing from the chromosome");
            //     }
            // }
            // System.out.println("Done");
        }

        for (Individual individual : population) {
            calculateFitness(individual);
        }
    }

    public void calculateFitness(Individual individual, boolean... flags) {
        // calculate the fitness of the individual
        // each individual consists of a list of patients and nurses
        // each nurse travels to each patient in their group
        // each nurse has a capacity and the patients have a demand, total demand cannot exceed capacity for a nurse
        // each patient has a time window, and the nurse must arrive within the time window, else nurse has to wait
        // each patient has a service time, and the nurse must spend that much time with the patient
        // the fitness is the total travel time (only) taken by all nurses to visit all patients

        Object[] result = calculateTimeAndTimeRequirement(individual);
        double totalTravelTime = (double) result[0];
        int numberOfPatientsMissed = (int) result[1];

        int[] groupwiseDemand = totalGroupwiseDemand(individual);

        int numberOfNursesOverCapacity = 0;
        for (int demand : groupwiseDemand) {
            if (demand > capacityNurse) {
                numberOfNursesOverCapacity++;
            }
        }

        // for (int demand : groupwiseDemand) {
        //     System.out.println(demand);
        // }
        // System.out.println(totalTravelTime);
        // System.out.println(numberOfPatientsMissed);
        // System.out.println(numberOfNursesOverCapacity);
        // number++;
        // boolean flag = flags.length > 0 ? true : false;

        individual.setFitness(totalTravelTime + numberOfPatientsMissed * 500 + numberOfNursesOverCapacity * 500);
        // System.out.println(individual.getFitness());
    }

    public int[] totalGroupwiseDemand(Individual individual) {
        // calculate the total demand of each group of patients
        // each nurse has a capacity and the patients have a demand, 
        // total demand cannot exceed capacity for a nurse (checked later)

        int[] individualChromosome = individual.getChromosome();
        int[] groupwiseDemand = new int[nbrNurses];

        int demand = 0;
        int group = 0;
        for (int patient : individualChromosome) {
            if (patient < 0) {
                groupwiseDemand[group] = demand;
                group++;
                demand = 0;
            } else {
                demand += patients.get(String.valueOf(patient)).getDemand();
            }
        }

        return groupwiseDemand;
    }

    public Object[] calculateTimeAndTimeRequirement(Individual individual) {
        // check if the individual meets the time requirement of each patient
        // each patient has a time window, and the nurse must arrive within the time window, else nurse has to wait
        // each patient has a service time, and the nurse must spend that much time with the patient

        int[] individualChromosome = individual.getChromosome().clone();
        int numberOfPatientsMissed = 0;

        // replace all <0 with 0
        for (int i = 0; i < individualChromosome.length; i++) {
            if (individualChromosome[i] < 0) {
                individualChromosome[i] = 0;
            }
        }

        int group = 0;
        double time = 0;

        int previous = individualChromosome[0];
        int current;
        double totalTravelTime = 0;
        
        // if the first element is a patient, then the nurse has to travel from the depot to the patient
        if (previous != 0) {
            totalTravelTime += travelTimes[0][previous];
            time += travelTimes[0][previous];

            int careTime = patients.get(String.valueOf(previous)).getCareTime();
            int endTime = patients.get(String.valueOf(previous)).getEndTime();
            if (time > endTime || time + careTime > endTime) {
                numberOfPatientsMissed++;
            } else if (time < patients.get(String.valueOf(previous)).getStartTime()) {
                time = patients.get(String.valueOf(previous)).getStartTime();
            }
            time += careTime;
        }

        // for each pair of patients, the nurse has to travel the distance between them
        // a nurse is denoted by a 0 in the individual, so when we encounter a 0 this will
        // indicate that the nurse has to travel to depot, and a new nurse will start from the depot
        for (int i = 1; i < individualChromosome.length; i++) {
            current = individualChromosome[i];

            totalTravelTime += travelTimes[previous][current];
            time += travelTimes[previous][current];

            if (current != 0) {
                int careTime = patients.get(String.valueOf(current)).getCareTime();
                int endTime = patients.get(String.valueOf(current)).getEndTime();
                if (time > endTime || time + careTime > endTime) {
                    numberOfPatientsMissed++;
                } else if (time < patients.get(String.valueOf(current)).getStartTime()) {
                    time = patients.get(String.valueOf(current)).getStartTime();
                }
                time += careTime;
            } else {
                group++;
                time = 0;
            }

            previous = current;
        }
        
        Object[] result = new Object[2];
        result[0] = totalTravelTime;
        result[1] = numberOfPatientsMissed; 
        return result;
    }

    public Individual[] parentSelection() {
        if (true) { // tournament selection
            // randomly select 5 individuals from the population
            List<Integer> list = new ArrayList<>();
            for (int i = 0; i < populationSize; i++) {
                list.add(i);
            }
            Collections.shuffle(list);
            int[] indices = list.stream().mapToInt(Integer::intValue).toArray();

            Individual[] tournament = new Individual[5];
            for (int i = 0; i < 5; i++) {
                tournament[i] = population[indices[i]];
            }

            return findBestIndividuals(tournament);
        } else {
            return findBestIndividuals(population);
        }
    }

    public Individual[] generationalParentSelection() {
        int lambda = 3 * populationSize;
        Individual[] parents = new Individual[lambda];
        int current_member = 1; int i = 1;

        

        // set double r to a value [0, 1/lambda]
        double r = Math.random() * (1/lambda);

        while (current_member <= lambda) {
            while (r <= 1) {

            }
        }


        return findBestIndividuals(population);
    }

    public Individual findOneBestIndividual(Individual[] localPopulation) {
        Individual bestIndividual = localPopulation[0];
        for (Individual individual : localPopulation) {
            if (individual.getFitness() < bestIndividual.getFitness()) {
                bestIndividual = individual;
            }
        }
        return bestIndividual;
    }

    public Individual[] findBestIndividuals(Individual[] localPopulation) {
        // find the two individuals in the population with the best fitness
        Individual parent1 = localPopulation[0];
        Individual parent2 = localPopulation[1];
        for (Individual individual : localPopulation) {
            if (individual.getFitness() < parent1.getFitness()) {
                parent2 = parent1;
                parent1 = individual;
            } else if (individual.getFitness() < parent2.getFitness()) {
                parent2 = individual;
            }
        }
        return new Individual[] {parent1, parent2};
    }

    public Individual[] crossover(int[] parent1, int[] parent2) {
        // remove the last nurse from both parents
        parent1 = Arrays.copyOf(parent1, parent1.length - 1);
        parent2 = Arrays.copyOf(parent2, parent2.length - 1);

        int chromosomeLength = parent1.length;

        int cuttingPoint1 = (int) (Math.random() * (chromosomeLength) - 1);
        int cuttingPoint2 = (int) (Math.random() * (chromosomeLength));

        if (cuttingPoint1 == cuttingPoint2) {
            cuttingPoint2++;
        }
        if (cuttingPoint1 > cuttingPoint2) {
            cuttingPoint1 = cuttingPoint1 + cuttingPoint2;
            cuttingPoint2 = cuttingPoint1 - cuttingPoint2;
            cuttingPoint1 = cuttingPoint1 - cuttingPoint2;
        }
        if (cuttingPoint1 == 0 && cuttingPoint2 == chromosomeLength - 1) {
            if (Math.random() < 0.5) {
                cuttingPoint1 += 2;
            } else {
                cuttingPoint2 -= 2;
            }
        }

        int[] child1 = new int[chromosomeLength];
        int[] child2 = new int[chromosomeLength];

        for (int i = 0; i < chromosomeLength; i++) {
            child1[i] = 0;
            child2[i] = 0;
        }

        for (int i = cuttingPoint1; i < cuttingPoint2 + 1; i++) {
            child1[i] = parent1[i];
            child2[i] = parent2[i];
        }

        // System.out.println("Cutting point 1: " + cuttingPoint1 + " Cutting point 2: " + cuttingPoint2);
        // System.out.println("Parent 1: " + Arrays.toString(parent1) + " Parent 2: " + Arrays.toString(parent2));
        // System.out.println("Child 1:  " + Arrays.toString(child1) + " Child 2:  " + Arrays.toString(child2));

        int parentIndex = cuttingPoint2 + 1;
        int childIndex = cuttingPoint2 + 1;
        if (cuttingPoint2 == chromosomeLength - 1) {
            parentIndex = 0;
            childIndex = 0;
        }

        if(false) {
            child1 = orderOneCrossover(child1, parent2, cuttingPoint1, cuttingPoint2, parentIndex, childIndex, chromosomeLength);
            child2 = orderOneCrossover(child2, parent1, cuttingPoint1, cuttingPoint2, parentIndex, childIndex, chromosomeLength);
        } else {
            child1 = partiallyMappedCrossover(child1, parent1, parent2, cuttingPoint1, cuttingPoint2, chromosomeLength);
            child2 = partiallyMappedCrossover(child2, parent2, parent1, cuttingPoint1, cuttingPoint2, chromosomeLength);
        }
        
        // add the last nurse back to the children
        child1 = Arrays.copyOf(child1, child1.length + 1);
        child2 = Arrays.copyOf(child2, child2.length + 1);
        child1[child1.length - 1] = -nbrNurses;
        child2[child2.length - 1] = -nbrNurses;

        // check that both children have all numbers -nbrNurses -> nbrPatients, print error message if not
        // for (int i = -nbrNurses; i < 101; i++) {
        //     boolean flag1 = false;
        //     boolean flag2 = false;
        //     for (int j = 0; j < chromosomeLength+1; j++) {
        //         if (child1[j] == i) {
        //             flag1 = true;
        //         }
        //         if (child2[j] == i) {
        //             flag2 = true;
        //         }
        //     }
        //     if (!flag1 && i!=0) {
        //         System.out.println("Number " + i + " is missing from child 1");
        //     }
        //     if (!flag2 && i!=0) {
        //         System.out.println("Number " + i + " is missing from child 2");
        //     }
        // }

        return new Individual[] {new Individual(child1, 0), new Individual(child2, 0)};
    }

    public int[] partiallyMappedCrossover(int[] child1, int[] parent1, int[] parent2, int cuttingPoint1, int cuttingPoint2, int chromosomeLength) {
        int parentIndex = cuttingPoint1;
        int parentIndex2;
        while (parentIndex < cuttingPoint2 + 1) {
            // check if parent2[parentIndex] is in child1
            boolean flag = false;
            for (int i = 0; i < chromosomeLength; i++) {
                if (child1[i] == parent2[parentIndex]) {
                    flag = true;
                    break;
                }
            }
            if (flag) {
                parentIndex++;
            } else {
                parentIndex2 = parentIndex;
                while (true) {
                    // find the index of int num = parent1[parentIndex2] in parent2
                    int num = parent1[parentIndex2];
                    int index = -1;
                    for (int j = 0; j < chromosomeLength; j++) {
                        if (parent2[j] == num) {
                            index = j;
                            break;
                        }
                    }
                    
                    if (index < cuttingPoint1 || index > cuttingPoint2) {
                        child1[index] = parent2[parentIndex];
                        parentIndex++;
                        break;
                    } else {
                        parentIndex2 = index;
                    }
                }
            }
        }

        for (int i = 0; i < chromosomeLength; i++) {
            if (child1[i] == 0) {
                child1[i] = parent2[i];
            }
        }

        return child1;
    }

    public int[] orderOneCrossover(int[] child1, int[] parent2, int cuttingPoint1, int cuttingPoint2, int parentIndex, int childIndex, int chromosomeLength) {
        while (true) {
            boolean flag = false;

            for (int i = 0; i < chromosomeLength; i++) {
                if (child1[i] == parent2[parentIndex]) {
                    flag = true;
                    break;
                }
            }

            if (flag) {
                if (parentIndex == chromosomeLength - 1) {
                    parentIndex = 0;
                } else {
                    parentIndex++;
                }
            } else {
                child1[childIndex] = parent2[parentIndex];
                childIndex++;

                if (childIndex == chromosomeLength) {
                    childIndex = 0;
                } 
                if (childIndex == cuttingPoint1) {
                    break;
                }
            }
        }

        // // check that all numbers 1 to n are present in the child
        // boolean noNumMissing = true;
        // for (int i = -24; i < 101; i++) {
        //     for (int j = 0; j < chromosomeLength; j++) {
        //         if (child1[j] == i) {
        //             break;
        //         }
        //         if (j == chromosomeLength - 1) {
        //             System.out.println("Number " + i + " is missing from the child");
        //             noNumMissing = false;
        //         }
        //     }
        // }
        // if(noNumMissing) {
        //     System.out.println("All numbers are present in the child");
        // }
        // System.out.println("Length: " + child1.length);

        return child1;
    }
    
    public void mutation(Individual child) {
        // 50% chance of mutation
        // if (Math.random() < 1) {
        int[] chromosome = child.getChromosome();

        // int temp = chromosome[index1];
        // chromosome[index1] = chromosome[index2];
        // chromosome[index2] = temp;
        // child.setChromosome(chromosome);
        // if (Math.random() < 0.5) { 
        // while (true) {

        //     double prob = Math.random();

        //     if (prob < 0.1) { // randomly swap two elements in the same group
        //         intraSwapMutation(child, chromosome);
        //     } else if (prob < 0.2) { // randomly move an element in the same group
        //         intraMoveMutation(child, chromosome);
        //     } else if (prob < 0.3) { // randomly swap two elements in the chromosome
        //         interSwapMutation(child, chromosome);
        //     } else if (prob < 0.4) { // randomly move an element from one group to another
        //         interMoveMutation(child, chromosome);
        //     } else if (prob < 0.45) { // pick two random indices and reverse the elements between them
        //         reverseMutation(child, chromosome);
        //     }

            if (Math.random() < 0.6) { // randomly swap two elements in the same group
                intraSwapMutation(child, chromosome);
            }
            if (Math.random() < 0.6) { // randomly move an element in the same group
                intraMoveMutation(child, chromosome);
            }
            if (Math.random() < 0.6) { // randomly swap two elements in the chromosome
                interSwapMutation(child, chromosome);
            }
            if (Math.random() < 0.6) { // randomly move an element from one group to another
                interMoveMutation(child, chromosome);
            }
            if (Math.random() < 0.1) { // pick two random indices and reverse the elements between them
                reverseMutation(child, chromosome);
            }
        //     if (Math.random() < 0.75){
        //         break;
        //     }
        // }
    }
    
    public void intraSwapMutation(Individual child, int[] chromosome) {
        int group = (int) (Math.random() * nbrNurses);
        int count = 0;
        int start = 0;
        int end = 0;
        if (group == 0) {
            start = 0;
            for (int i = 0; i < chromosome.length; i++) {
                if (chromosome[i] < 0) {
                    end = i;
                    break;
                }
            }
        } else {
            for (int i = 0; i < chromosome.length; i++) {
                if (chromosome[i] < 0) {
                    count++;
                    if (count == group - 1) {
                        start = i + 1;
                    } else if (count == group) {
                        end = i;
                        break;
                    }
                }
                
            }
        }
        if (end - start > 2) {
            // find two random indices in the group and swap them
            int index1 = (int) (Math.random() * (end - start - 1) + start);
            int index2 = (int) (Math.random() * (end - start - 1) + start);
            int temp = chromosome[index1];
            chromosome[index1] = chromosome[index2];
            chromosome[index2] = temp;
            child.setChromosome(chromosome);
        }
    }

    public void intraMoveMutation(Individual child, int[] chromosome) {
        int group = (int) (Math.random() * nbrNurses);
        int count = 0;
        int start = 0;
        int end = 0;
        if (group == 0) {
            start = 0;
            for (int i = 0; i < chromosome.length; i++) {
                if (chromosome[i] < 0) {
                    end = i;
                    break;
                }
            }
        } else {
            for (int i = 0; i < chromosome.length; i++) {
                if (chromosome[i] < 0) {
                    count++;
                    if (count == group - 1) {
                        start = i + 1;
                    } else if (count == group) {
                        end = i;
                        break;
                    }
                }
                
            }
        }
        if (end - start > 2) {
            // find two random indices in the group and swap them
            int index1 = (int) (Math.random() * (end - start - 1) + start);
            int index2 = (int) (Math.random() * (end - start - 1) + start);

            int temp = chromosome[index1];

            if (index1 > index2) {
                // move all elements between index1 and index2 one step to the right
                for (int i = index1; i > index2; i--) {
                    chromosome[i] = chromosome[i - 1];
                }
            } else if (index1 < index2) {
            // move all elements between index1 and index2 one step to the left
                for (int i = index1; i < index2; i++) {
                    chromosome[i] = chromosome[i + 1];
                }
            }

            chromosome[index2] = temp;
            child.setChromosome(chromosome);
        }
    }

    public void interSwapMutation(Individual child, int[] chromosome) {
        int index1 = (int) (Math.random() * (chromosome.length - 1));
        int index2 = (int) (Math.random() * (chromosome.length - 1));

        int temp = chromosome[index1];
        chromosome[index1] = chromosome[index2];
        chromosome[index2] = temp;
        child.setChromosome(chromosome);
    }

    public void interMoveMutation(Individual child, int[] chromosome) {
        int index1 = (int) (Math.random() * (chromosome.length - 1));
        int index2 = (int) (Math.random() * (chromosome.length - 1));

        int temp = chromosome[index1];

        if (index1 > index2) {
            // move all elements between index1 and index2 one step to the right
            for (int i = index1; i > index2; i--) {
                chromosome[i] = chromosome[i - 1];
            }
        } else if (index1 < index2) {
        // move all elements between index1 and index2 one step to the left
            for (int i = index1; i < index2; i++) {
                chromosome[i] = chromosome[i + 1];
            }
        }

        chromosome[index2] = temp;
        child.setChromosome(chromosome);
    }

    public void reverseMutation(Individual child, int[] chromosome) {
        int index1 = (int) (Math.random() * (chromosome.length - 1));
        int index2 = (int) (Math.random() * (chromosome.length - 1));

        if (index1 > index2) {
            index1 = index1 + index2;
            index2 = index1 - index2;
            index1 = index1 - index2;
        }
        int[] temp = Arrays.copyOfRange(chromosome, index1, index2);
        for (int i = 0; i < temp.length; i++) {
            chromosome[index1 + i] = temp[temp.length - 1 - i];
        }
        child.setChromosome(chromosome);
    }

    public void survivorSelection(Individual child1, Individual child2) {
        // find the two individuals in the population with the worst fitness and replace them with the children
        double worstFitness1 = population[0].getFitness();
        double worstFitness2 = population[1].getFitness();
        int worstFitnessIndex1 = 0;
        int worstFitnessIndex2 = 1;

        for (int i = 2; i < population.length; i++) {
            if (population[i].getFitness() > worstFitness1) {
                worstFitnessIndex2 = worstFitnessIndex1;
                worstFitness2 = worstFitness1;
                worstFitnessIndex1 = i;
                worstFitness1 = population[i].getFitness();
            } else if (population[i].getFitness() > worstFitness2) {
                worstFitnessIndex2 = i;
                worstFitness2 = population[i].getFitness();
            }
        }

        if (child2.getFitness() < worstFitness1) {
            population[worstFitnessIndex1] = child2;
        } else if (child2.getFitness() < worstFitness2) {
            population[worstFitnessIndex2] = child2;
        }

        if (child1.getFitness() < worstFitness2) {
            population[worstFitnessIndex2] = child1;
        } 
    }

    public void printSolution(Individual bestIndividual) {
        System.out.println();
        System.out.print("[[");
        if (bestIndividual.getChromosome()[0] < 0) {
            System.out.print("], [");
        } else {
            System.out.print(bestIndividual.getChromosome()[0]);
        }
        int previous = bestIndividual.getChromosome()[0];
        for (int j = 1; j < bestIndividual.getChromosome().length - 1; j++) {
            if (bestIndividual.getChromosome()[j] < 0) {
                System.out.print("], [");
            } else if (previous > 0) {
                System.out.print(", " + bestIndividual.getChromosome()[j]);
            } else {
                System.out.print(bestIndividual.getChromosome()[j]);
            }
            previous = bestIndividual.getChromosome()[j];
        }
        System.out.println("]]");
        System.out.println();
    }

    public void saveSolution(Individual bestIndividual, double travelTime) {
        // save the best individual to a file
        // start by building the String that represents the solution
        StringBuilder solution = new StringBuilder();

        solution.append("[[");
        if (bestIndividual.getChromosome()[0] < 0) {
            solution.append("], [");
        } else {
            solution.append("[" + bestIndividual.getChromosome()[0] + "," + patients.get(String.valueOf(bestIndividual.getChromosome()[0])).getXCoord() + "," + patients.get(String.valueOf(bestIndividual.getChromosome()[0])).getYCoord() + "]");
        }
        int previous = bestIndividual.getChromosome()[0];
        for (int j = 1; j < bestIndividual.getChromosome().length - 1; j++) {
            if (bestIndividual.getChromosome()[j] < 0) {
                solution.append("], [");
            } else if (previous > 0) {
                solution.append(", " + "[" + bestIndividual.getChromosome()[j] + "," + patients.get(String.valueOf(bestIndividual.getChromosome()[j])).getXCoord() + "," + patients.get(String.valueOf(bestIndividual.getChromosome()[j])).getYCoord() + "]");
            } else {
                solution.append("[" + bestIndividual.getChromosome()[j] + "," + patients.get(String.valueOf(bestIndividual.getChromosome()[j])).getXCoord() + "," + patients.get(String.valueOf(bestIndividual.getChromosome()[j])).getYCoord() + "]");
            }
            previous = bestIndividual.getChromosome()[j];
        }
        solution.append("]]");

        StringBuilder hub = new StringBuilder();
        hub.append("[0, " + depot.getXCoord() + ", " + depot.getYCoord() + "]");

        Map<String, Object> solutionMap = new HashMap<>();
        solutionMap.put("data", solution.toString());
        solutionMap.put("hub", hub.toString());
        solutionMap.put("travel_time", String.valueOf(travelTime));

        Gson gson = new Gson();
        String json = gson.toJson(solutionMap);

        try (FileWriter file = new FileWriter("file2.json")) {
            file.write(json);
            System.out.println("Successfully Copied JSON Object to File...");
            // System.out.println("\nJSON Object: " + json);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public Individual[] getPopulation() {
        return population;
    }

    public void setPopulation(Individual[] population) {
        this.population = population;
    }
    public static void main(String[] args) {
        GeneticAlgorithm GA = new GeneticAlgorithm(15, "src/main/resources/train/train_9.json", null);

        // Individual individ = new Individual(new int[]{10, 11, -16, 55, 54, 53, -19, 56, 58, -20, 32, 33, 37, 34, 22, -18, 57, 31, 35, 38, 39, 36, 52, -22, 48, 51, -5, 43, 42, 40, 44, 46, 47, -14, 21, -24, 81, 78, 61, 64, 68, 66, -8, 18, 19, 16, 14, 2, 1, 75, -7, 27, 29, -6, 98, 96, 97, 100, 99, -4, 3, 7, 8, 9, 6, 4, 90, 87, 95, 94, 92, 93, -11, 74, 72, 60, 59, -10, 24, 23, -2, 5, 86, 73, 77, 80, -13, 25, 30, 28, 26, 50, 49, -1, 62, -23, 76, 71, 70, 79, 89, 91, -12, 67, 41, 45, 69, -3, 20, -21, 65, 63, -15, 83, 82, 84, -9, 85, 88, -17, 13, 17, 15, 12, -25}, 0);
        // System.out.println(Arrays.toString(GA.cutLongestTravel(individ).getChromosome()));


        // GA.generateGreedyPopulation();
        // System.out.println(Arrays.toString(GA.getPopulation()[0].getChromosome()));
        // System.out.println(GA.getPopulation()[0].getFitness());

        GA.generatePopulation();
        Individual[] OGpopulation = GA.getPopulation();

        Individual bestTotalIndividual = OGpopulation[0];
        for (int i = 0; i < 10; i++) {
            // double bestFitnessBefore = population[0].getFitness();
            // for (Individual individual : population) {
            //     if (individual.getFitness() < bestFitnessBefore) {
            //         bestFitnessBefore = individual.getFitness();
            //     }
            // }
            // System.out.println("Best fitness before: " + bestFitnessBefore);
            
            GA.runGA(null, 1000000);

            // Individual individ = new Individual(new int[]{10, 11, -16, 55, 54, 53, 56, 58, -20, 32, 33, 37, 34, 22, -18, 57, 31, 35, 38, 39, 36, 52, -22, 48, 51, -5, 43, 42, 40, 44, 46, 47, -14, 21, -24, 81, 78, 61, 64, 68, 66, -8, 18, 19, 16, 14, 2, 1, 75, -7, 27, 29, -6, 98, 96, 97, 100, 99, -4, 3, 7, 8, 9, 6, 4, -15, -9, 90, 87, 95, 94, 92, 93, -11, 74, 72, 60, 59, -10, 24, 23, -2, -19, 5, 86, 73, 77, 80, -13, 25, 30, 28, 26, 50, 49, -1, 62, -23, 76, 71, 70, 79, 89, 91, -12, 67, 41, 45, 69, -3, 20, -21, 65, 63, 83, 82, 84, 85, 88, -17, 13, 17, 15, 12, -25}, 0);
            // GA.saveSolution(individ);
            // GA.calculateFitness(individ);
            // System.out.println(individ.getFitness());

            // Individual individ = new Individual(new int[]{1, 2, 3, 4, -2, 5, 6, 7, 8, -3}, 0);
            // GA.interMoveMutation(individ, individ.getChromosome());
            // GA.intraMoveMutation(individ, individ.getChromosome());
            // System.out.println(Arrays.toString(individ.getChromosome()));

            Individual[] population = GA.getPopulation();
            Individual bestIndividual = population[0];
            double bestFitnessAfter = population[0].getFitness();
            for (Individual individual : population) {
                if (individual.getFitness() < bestFitnessAfter) {
                    bestIndividual = individual;
                    bestFitnessAfter = individual.getFitness();
                }
            }
            System.out.println("Best fitness after: " + bestFitnessAfter);
            if (bestFitnessAfter < bestTotalIndividual.getFitness()) {
                bestTotalIndividual = bestIndividual;
            }
        }
        GA.saveSolution(bestTotalIndividual, bestTotalIndividual.getFitness());
        GA.printSolution(bestTotalIndividual); 
    }
}

