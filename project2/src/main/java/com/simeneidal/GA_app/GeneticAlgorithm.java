package com.simeneidal.GA_app;
import org.apache.commons.math3.ml.clustering.*;
import org.apache.commons.math3.ml.distance.EuclideanDistance;

import java.util.*;
import java.util.stream.Collectors;

import com.google.gson.Gson;

import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;

import java.util.concurrent.*;

public class GeneticAlgorithm implements Callable<Individual> {

    private int counting = 1;

    private int populationSize;

    private int nbrNurses;
    private int capacityNurse;
    private JsonData.Depot depot;
    private int nbrPatients;
    private Map<String, JsonData.Patient> patients;
    private double[][] travelTimes;
    private Individual[] population;

    private boolean tournamentParentSelection;
    private boolean orderOneCrossover;

    Individual[] localPopulation;
    int nbrGenerations;
    boolean replacePopulation;
    boolean replaceGreedy;
    boolean parentClone;
    int populationGeneration;
    int numberOfClusters;
    int numberOfIndividualsInClusters;
    boolean diversity;
    boolean shuffleNurses;

    /**
     * Constructs a GeneticAlgorithm object with the specified parameters.
     *
     * @param populationSize The size of the population for the genetic algorithm.
     * @param filePath The file path to the JSON data used for initializing the genetic algorithm.
     * @param tournamentParentSelection A boolean value indicating whether tournament parent selection is used.
     * @param orderOneCrossover A boolean value indicating whether order one crossover is used.
     */
    public GeneticAlgorithm(int populationSize, String filePath, boolean tournamentParentSelection, boolean orderOneCrossover, Individual[] localPopulation, int nbrGenerations, boolean replacePopulation, boolean replaceGreedy, boolean parentClone, int populationGeneration, int numberOfClusters, int numberOfIndividualsInClusters, boolean diversity, boolean shuffleNurses) {
        this.populationSize = populationSize;
        this.tournamentParentSelection = tournamentParentSelection;
        this.orderOneCrossover = orderOneCrossover;
        this.localPopulation = localPopulation;
        this.nbrGenerations = nbrGenerations;
        this.replacePopulation = replacePopulation;
        this.replaceGreedy = replaceGreedy;
        this.parentClone = parentClone;
        this.populationGeneration = populationGeneration;
        this.numberOfClusters = numberOfClusters;
        this.numberOfIndividualsInClusters = numberOfIndividualsInClusters;
        this.diversity = diversity;
        this.shuffleNurses = shuffleNurses;

        nbrNurses = 0;
        capacityNurse = 0;
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
            nbrNurses = jsonData.getNbrNurses();
            capacityNurse = jsonData.getCapacityNurse();
            depot = jsonData.getDepot();
            patients = jsonData.getPatients();
            travelTimes = jsonData.getTravelTimes(); 

            // find the number of patients
            nbrPatients = patients.size();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    @Override
    public Individual call() throws Exception {
        // Used for multithreading
        // Run the GA and return the best individual
        runGA();
        return findOneBestIndividual(population);
    }

    /**
     * Runs the Genetic Algorithm with the given parameters.
     *
     * @param localPopulation   The initial population of individuals. If null, a new population will be generated.
     * @param nbrGenerations    The number of generations to run the algorithm.
     * @param replacePopulation Specifies whether to replace the entire population at regular intervals (if stale).
     * @param replaceGreedy     Specifies whether to replace the population with random (greedy) individuals when replacing.
     * @param parentClone       Specifies whether to possibly clone parents before applying mutation.
     * @param populationGeneration The method used to generate the initial population.
     * @param numberOfClusters The number of clusters to generate in the clustered population.
     * @param numberOfIndividualsInClusters The number of individuals to generate for each cluster in the clustered population.
     */
    public void runGA() {        
        if (Objects.isNull(localPopulation)) {
            if (populationGeneration == 0) {
                generatePopulation();
            } else if (populationGeneration == 1) {
                generateGreedyPopulation();
            } else if (populationGeneration == 2) {
                generateClusteredPopulation(numberOfClusters, numberOfIndividualsInClusters, diversity);
            } else {
                throw new IllegalArgumentException("Invalid population generation method");
            }
        } else {
            population = localPopulation;
        }

        Individual nextBestIndividual = population[0];
        Individual bestIndividual = population[0];
        boolean flag = false;
        
        for (int i = 0; i < nbrGenerations; i++) {
            
            if(true && (i+1) % 100000 == 0) {
                // find the best individual and print its fitness
                bestIndividual = findOneBestIndividual(population);
                System.out.println("Fitness: " + bestIndividual.getFitness());
            }

            if ((i+1) % 1000000 == 0) {
                System.out.println("Start shuffling");
                for (Individual individual : population) {
                    randomizeOrderOfNurses(individual);
                }
            }
            if ((i+1) % 500000 == 0) {
                System.out.println("Start combining");
                bestIndividual = findOneBestIndividual(population);
                combineRoutes(bestIndividual);
            }
            if ((i+1) % 750000 == 0) {
                System.out.println("Start splitting");
                bestIndividual = findOneBestIndividual(population);
                counting = 1;
                divideRoutes(bestIndividual.getChromosome().clone());
            }

            if (replacePopulation && i % 100000 == 0) {
                // find the next best individual
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
                if (nextBestFitness == nextBestIndividual.getFitness()) {
                    if (flag) {
                        if (replaceGreedy) {
                            // replace the entire population except for the most fit individual, with random (greedy) individuals
                            //generatePopulation();
                            generateGreedyPopulation();
                            population[0] = bestIndividual;
                            System.out.println("Population reset");
                            flag = false;
                        } else {
                            // make copies of the best individual and use destroy operators to create new individuals
                            cutLongestTravel(bestIndividual);
                            System.out.println("Population made from best individual");
                            flag = false;
                        }
                    } else {
                        flag = true;
                        System.out.println("Flag set");
                    }
                }
            }

            Individual[] parents = parentSelection();
            if (parentClone && Math.random() < 0.1) {
                Individual parent1 = new Individual(parents[0].getChromosome().clone(), 0);
                Individual parent2 = new Individual(parents[1].getChromosome().clone(), 0);

                mutation(parent1); mutation(parent2);

                if (shuffleNurses) {
                    randomizeOrderOfNurses(parent1); randomizeOrderOfNurses(parent2);
                }
                
                calculateFitness(parent1); calculateFitness(parent2);
                survivorSelection(parent1, parent2);
            } else {
                Individual[] children = crossover(parents[0].getChromosome(), parents[1].getChromosome());
                mutation(children[0]); mutation(children[1]);
                if (shuffleNurses) {
                    randomizeOrderOfNurses(children[0]); randomizeOrderOfNurses(children[1]);
                }
                
                calculateFitness(children[0]); calculateFitness(children[1]);
                survivorSelection(children[0], children[1]);
            }
        }
    }
    
    /**
     * Cuts the longest travel in the best individual's chromosome and updates the population accordingly.
     * 
     * @param bestIndividual The best individual in the population.
     */
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

    /**
     * Generates the initial population of individuals for the genetic algorithm.
     * Each individual is created with a chromosome representing the assignment of nurses to patients.
     * The chromosome is randomly generated by shuffling a list of patients and adding zeros to represent nurses.
     * The fitness of each individual is calculated after generation.
     */
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

    /**
     * Generates a population of individuals using a greedy approach.
     * Each individual's chromosome is created by iteratively placing patients in a way that minimizes the increase in fitness.
     * The population size is determined by the 'populationSize' variable.
     * 
     * @throws IllegalArgumentException if 'populationSize', 'nbrNurses', or 'nbrPatients' is less than or equal to zero.
     */
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
            Collections.shuffle(list);

            int[] temp = list.stream().mapToInt(Integer::intValue).toArray();
            int[] chromosome = new int[temp.length + 25];

            // add nbrNurses - 1 0s to the list
            for (int j = 0; j < nbrNurses - 1; j++) {
                chromosome[j] = 0;
            }

            int count = nbrNurses - 1;
            int start = 0;

            // iterate though temp and place each int where they cause the least increase in fitness
            for (int j = 0; j < temp.length; j++) {
                count++;
                double minIncrease = 10000;
                int minIncreaseIndex = 0;
                for (int k = start; k < count; k++) {
                    if (travelTimes[chromosome[k]][temp[j]] + travelTimes[chromosome[k+1]][temp[j]] < minIncrease) {
                        minIncrease = travelTimes[chromosome[k]][temp[j]] + travelTimes[chromosome[k+1]][temp[j]];
                        minIncreaseIndex = k;
                    }
                }

                for (int k = count; k > minIncreaseIndex; k--) {
                    chromosome[k] = chromosome[k-1];
                }

                chromosome[minIncreaseIndex] = temp[j];

                // Spread individuals out in the chromosome
                if (j == 0) {
                    start = 1;
                } else if ( 0 < j && j < 26) {
                    start = j * 2;
                } else {
                    start = 0;
                }
            }

            int num = -1;
            for (int j = 0; j < chromosome.length; j++) {
                if (chromosome[j] == 0) {
                    chromosome[j] = num;
                    num--;
                }
            }

            population[i].setChromosome(chromosome);
        }

        for (Individual individual : population) {
            calculateFitness(individual);
        }
    }

    /**
     * Generates a clustered population of individuals based on patient data.
     * Each individual represents a clustering solution obtained through the K-means clustering algorithm.
     * The method converts patient data into points for clustering, performs clustering, and converts the resulting clusters back to patients.
     * It then creates four individuals for each cluster, each with a different sorting criterion.
     * The generated population is stored in the 'population' array.
     * 
     * @param numberOfClusters The number of clusters to generate in the population.
     * @param numberOfIndividualsInClusters The number of individuals to generate for each cluster.
     */
    public void generateClusteredPopulation(int numberOfClusters, int numberOfIndividualsInClusters, boolean diversity) {
        Individual[] tempPopulation = new Individual[populationSize];

        // for all patients in the map, id is set to the key value
        for (Map.Entry<String, JsonData.Patient> entry : patients.entrySet()) {
            String id = entry.getKey();
            JsonData.Patient patient = entry.getValue();
            // id must be an integer
            patient.setId(Integer.parseInt(id));
        }
        // Convert patients to points for clustering
        List<Clusterable> patientPoints = new ArrayList<>();
        for (JsonData.Patient patient : patients.values()) {
            patientPoints.add(new DoublePoint(new double[]{patient.getXCoord(), patient.getYCoord()}));
        }

        for (int i = 0; i < numberOfClusters; i++) {
            if (diversity) {
                for (int j = 0; j < numberOfIndividualsInClusters * 4; j+=4) {
                    // Perform clustering
                    KMeansPlusPlusClusterer<Clusterable> clusterer = new KMeansPlusPlusClusterer<>(3 + i, 5000, new EuclideanDistance());
                    List<CentroidCluster<Clusterable>> clusters = clusterer.cluster(patientPoints);
    
                    // Convert clusters back to patients
                    List<List<JsonData.Patient>> patientClusters = new ArrayList<>();
                    for (CentroidCluster<Clusterable> cluster : clusters) {
                        List<JsonData.Patient> patientCluster = new ArrayList<>();
                        for (Clusterable point : cluster.getPoints()) {
                            patientCluster.add(findPatient(patients, point));
                        }
                        patientClusters.add(patientCluster);
                    }
    
                    // add 4 individuals to the population for diversity
                    Individual individual1 = makeIndividualOfPatientClusters(patientClusters);
                    // randomizeOrderOfNurses(individual1);
                    tempPopulation[i*40 + j] = individual1;
    
                    sortOnEndTime(patientClusters);
                    Individual individual2 = makeIndividualOfPatientClusters(patientClusters);
                    // randomizeOrderOfNurses(individual2);
                    // tempPopulation[i*40 + j + 1] = individual2;
                    tempPopulation[i*40 + j + 1] = individual2;
    
    
                    sortOnEndTimeAndDistance(patientClusters);
                    Individual individual3 = makeIndividualOfPatientClusters(patientClusters);
                    // randomizeOrderOfNurses(individual3);
                    tempPopulation[i*40 + j + 2] = individual3;
    
                    sortOnNearestNeighbor(patientClusters);
                    Individual individual4 = makeIndividualOfPatientClusters(patientClusters);
                    // randomizeOrderOfNurses(individual4);
                    tempPopulation[i*40 + j + 3] = individual4;
                }
            } else {
                for (int j = 0; j < numberOfIndividualsInClusters; j++) {
                    // Perform clustering
                    KMeansPlusPlusClusterer<Clusterable> clusterer = new KMeansPlusPlusClusterer<>(3 + i, 5000, new EuclideanDistance());
                    List<CentroidCluster<Clusterable>> clusters = clusterer.cluster(patientPoints);

                    // Convert clusters back to patients
                    List<List<JsonData.Patient>> patientClusters = new ArrayList<>();
                    for (CentroidCluster<Clusterable> cluster : clusters) {
                        List<JsonData.Patient> patientCluster = new ArrayList<>();
                        for (Clusterable point : cluster.getPoints()) {
                            patientCluster.add(findPatient(patients, point));
                        }
                        patientClusters.add(patientCluster);
                    }

                    sortOnEndTime(patientClusters);
                    Individual individual2 = makeIndividualOfPatientClusters(patientClusters);
                    tempPopulation[i*10 + j] = individual2;
                }
            }
        }

        population = tempPopulation;
    }

    /**
     * Helper method for generating a clustered population of individuals.
     * Represents an individual in the genetic algorithm.
     */
    private Individual makeIndividualOfPatientClusters(List<List<JsonData.Patient>> patientClusters) {
        int[] tempArray = new int[nbrPatients + nbrNurses];
        int index = 0;
        int number = 1;
        // for each cluster, add the patients to the chromosome, divide each cluster with a negative number
        for (List<JsonData.Patient> cluster : patientClusters) {
            for (JsonData.Patient patient : cluster) {
                tempArray[index] = patient.getId();
                index++;
            }
            tempArray[index] = -number;
            number++;
            index++;
        }

        for (int i = index; i < tempArray.length; i++) {
            tempArray[i] = -number;
            number++;
        }

        Individual individual = new Individual(tempArray, 0);
        calculateFitness(individual);

        return individual;
    }

    /**
     * Helper method for generating a clustered population of individuals.
     * Represents a patient in the JSON data.
     */
    private JsonData.Patient findPatient(Map<String, JsonData.Patient> patients, Clusterable point) {
        for (JsonData.Patient patient : patients.values()) {
            if (patient.getXCoord() == point.getPoint()[0] && patient.getYCoord() == point.getPoint()[1]) {
                return patient;
            }
        }
        return null;
    }

    /**
     * Sorts the patient clusters based on the end time of each patient.
     *
     * @param patientClusters the list of patient clusters to be sorted
     */
    public void sortOnEndTime(List<List<JsonData.Patient>> patientClusters) {
        for (List<JsonData.Patient> cluster : patientClusters) {
            // Sort the patients based on end_time
            cluster.sort(Comparator.comparing(JsonData.Patient::getEndTime));
        }
    }

    /**
     * Sorts the patient clusters based on the end time and distance.
     * 
     * @param patientClusters the list of patient clusters to be sorted
     */
    public void sortOnEndTimeAndDistance(List<List<JsonData.Patient>> patientClusters) {
        for (List<JsonData.Patient> cluster : patientClusters) {
            if (cluster.size() < 3) {
                continue;
            }
            // Sort patients based on distance to depot
            cluster.sort(Comparator.comparingDouble(patient -> calculateDistance(depot, patient)));

            // Get two closest patients
            JsonData.Patient closestPatient1 = cluster.get(0);
            JsonData.Patient closestPatient2 = cluster.get(1);

            // Remove two closest patients from cluster
            cluster.remove(closestPatient1);
            cluster.remove(closestPatient2);

            // Sort the rest of the patients based on end_time
            cluster.sort(Comparator.comparing(JsonData.Patient::getEndTime));

            // Add the closest patients back to the cluster
            if (closestPatient1.getEndTime() < closestPatient2.getEndTime()) {
                cluster.add(0, closestPatient1);
                cluster.add(closestPatient2);
            } else {
                cluster.add(0, closestPatient2);
                cluster.add(closestPatient1);
            }
        }
    }

    /**
     * Sorts the patient clusters based on the nearest neighbor algorithm.
     * Each cluster is sorted in ascending order of the distance between the depot and each patient.
     * The patients within each cluster are then sorted in ascending order of the distance to their nearest neighbor.
     * The sorted paths are then added back to their respective clusters.
     *
     * @param patientClusters the list of patient clusters to be sorted
     */
    public void sortOnNearestNeighbor(List<List<JsonData.Patient>> patientClusters) {
        for (List<JsonData.Patient> cluster : patientClusters) {
            List<JsonData.Patient> path = new ArrayList<>();
            
            cluster.sort(Comparator.comparingDouble(patient -> calculateDistance(depot, patient)));
            JsonData.Patient closestPatient = cluster.get(0);
            path.add(closestPatient);
            cluster.remove(closestPatient);

            while (!cluster.isEmpty()) {
                final JsonData.Patient current = closestPatient;

                // Sort patients based on distance to current patient
                cluster.sort(Comparator.comparingDouble(patient -> calculateDistance(current, patient)));

                // Get the closest patient and add it to the path
                closestPatient = cluster.get(0);
                path.add(closestPatient);

                // Remove the closest patient from the cluster
                cluster.remove(closestPatient);
            }

            // Now 'path' contains the patients sorted based on distance to each other
            // Replace the original cluster with the sorted path
            cluster.addAll(path);
        }
    }

    /**
     * Calculates the Euclidean distance between two patients.
     *
     * @param patient1 the first patient
     * @param patient2 the second patient
     * @return the Euclidean distance between the two patients
     */
    private double calculateDistance(JsonData.Patient patient1, JsonData.Patient patient2) {
        double dx = (patient1.getXCoord() - patient2.getXCoord());
        double dy = (patient1.getYCoord() - patient2.getYCoord());
        return Math.sqrt(dx * dx + dy * dy);
    }

    /**
     * Calculates the Euclidean distance between the depot and a patient.
     *
     * @param depot    the depot object representing the starting point
     * @param patient2 the patient object representing the destination point
     * @return the Euclidean distance between the depot and the patient
     */
    private double calculateDistance(JsonData.Depot depot, JsonData.Patient patient2) {
        double dx = (depot.getXCoord() - patient2.getXCoord());
        double dy = (depot.getYCoord() - patient2.getYCoord());
        return Math.sqrt(dx * dx + dy * dy);
    }

    /**
     * Calculates the fitness of an individual in the genetic algorithm.
     * The fitness is determined by the total travel time taken by all nurses to visit all patients.
     * Each individual consists of a list of patients and nurses.
     * Each nurse travels to each patient in their group.
     * Each nurse has a capacity and the patients have a demand.
     * The total demand cannot exceed the capacity for a nurse.
     * Each patient has a time window, and the nurse must arrive within the time window, else the nurse has to wait.
     * Each patient has a service time, and the nurse must spend that much time with the patient.
     * The fitness is calculated based on the total travel time, the number of patients missed, the number of nurses over the capacity,
     * and the number of nurses over the return time.
     *
     * @param individual The individual for which to calculate the fitness.
     */
    public void calculateFitness(Individual individual) {
        Object[] result = calculateTimeAndTimeRequirement(individual);
        double totalTravelTime = (double) result[0];
        int numberOfPatientsMissed = (int) result[1];
        int numberOfNursesOverReturnTime = (int) result[2];

        int[] groupwiseDemand = totalGroupwiseDemand(individual);

        int numberOfNursesOverCapacity = 0;
        for (int demand : groupwiseDemand) {
            if (demand > capacityNurse) {
                numberOfNursesOverCapacity++;
            }
        }

        individual.setFitness(totalTravelTime + numberOfPatientsMissed * 500 + numberOfNursesOverCapacity * 500 + numberOfNursesOverReturnTime * 500);
    }

    /**
     * Calculates the total demand of each group of patients in the given individual.
     * Each nurse has a capacity and the patients have a demand. The total demand cannot exceed the capacity for a nurse (checked later).
     *
     * @param individual The individual for which to calculate the total groupwise demand.
     * @return An array of integers representing the total demand of each group of patients.
     */
    public int[] totalGroupwiseDemand(Individual individual) {
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

    /**
     * Calculates the total travel time, number of patients missed, and number of nurses over return time
     * based on the given individual's chromosome.
     *
     * @param individual The individual for which to calculate the time and time requirement.
     * @return An array of objects containing the total travel time, number of patients missed,
     *         and number of nurses over return time.
     */
    public Object[] calculateTimeAndTimeRequirement(Individual individual) {
        int[] individualChromosome = individual.getChromosome().clone();
        int numberOfPatientsMissed = 0;
        int numberOfNursesOverReturnTime = 0;

        // replace all <0 with 0
        for (int i = 0; i < individualChromosome.length; i++) {
            if (individualChromosome[i] < 0) {
                individualChromosome[i] = 0;
            }
        }

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
                if (time > depot.getReturnTime()) {
                    numberOfNursesOverReturnTime++;
                }
                time = 0;
            }

            previous = current;
        }
        
        Object[] result = new Object[3];
        result[0] = totalTravelTime;
        result[1] = numberOfPatientsMissed;
        result[2] = numberOfNursesOverReturnTime;
        return result;
    }

    /**
     * Performs parent selection to choose individuals for reproduction.
     * If tournament parent selection is enabled, randomly selects 30 individuals from the population
     * and returns the best individuals from the tournament.
     * If tournament parent selection is disabled, returns the best individuals from the entire population.
     *
     * @return An array of Individual objects representing the selected parents for reproduction.
     */
    public Individual[] parentSelection() {
        if (tournamentParentSelection) { // tournament selection
            // randomly select 30 individuals from the population
            List<Integer> list = new ArrayList<>();
            for (int i = 0; i < populationSize; i++) {
                list.add(i);
            }
            Collections.shuffle(list);
            int[] indices = list.stream().mapToInt(Integer::intValue).toArray();

            Individual[] tournament = new Individual[30];
            for (int i = 0; i < 30; i++) {
                tournament[i] = population[indices[i]];
            }

            return findBestIndividuals(tournament);
        } else {
            return findBestIndividuals(population);
        }
    }

    /**
     * Finds the individual in localPopulation with the best fitness, and returns it
     *
     * @param localPopulation list of Individual objects to search for the best individual
     * @return the Individual with the best fitness
     */
    public Individual findOneBestIndividual(Individual[] localPopulation) {
        Individual bestIndividual = localPopulation[0];
        for (Individual individual : localPopulation) {
            if (individual.getFitness() < bestIndividual.getFitness()) {
                bestIndividual = individual;
            }
        }
        return bestIndividual;
    }

    /**
     * Finds the two individuals in the population with the best fitness.
     *
     * @param localPopulation the population of individuals
     * @return an array containing the two individuals with the best fitness
     */
    public Individual[] findBestIndividuals(Individual[] localPopulation) {
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

    /**
     * Performs crossover operation on two parent individuals to create two child individuals.
     * 
     * @param parent1 The first parent individual represented as an array of integers.
     * @param parent2 The second parent individual represented as an array of integers.
     * @return An array of two child individuals created through crossover.
     */
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

        int parentIndex = cuttingPoint2 + 1;
        int childIndex = cuttingPoint2 + 1;
        if (cuttingPoint2 == chromosomeLength - 1) {
            parentIndex = 0;
            childIndex = 0;
        }

        if(orderOneCrossover) {
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

        return new Individual[] {new Individual(child1, 0), new Individual(child2, 0)};
    }

    /**
     * Performs partially mapped crossover between two parent arrays to generate a child array.
     *
     * @param child1          The child array to be generated.
     * @param parent1         The first parent array.
     * @param parent2         The second parent array.
     * @param cuttingPoint1   The starting index of the crossover segment.
     * @param cuttingPoint2   The ending index of the crossover segment.
     * @param chromosomeLength The length of the chromosome arrays.
     * @return The child array generated through partially mapped crossover.
     */
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

    /**
     * Performs order one crossover between two parent arrays to generate a child array.
     * 
     * @param child1 The first child array.
     * @param parent2 The second parent array.
     * @param cuttingPoint1 The first cutting point for crossover.
     * @param cuttingPoint2 The second cutting point for crossover.
     * @param parentIndex The index of the parent array to start from.
     * @param childIndex The index of the child array to start from.
     * @param chromosomeLength The length of the chromosome.
     * @return The child array after performing order one crossover.
     */
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

        return child1;
    }
    
    /**
     * Randomizes the order of nurses in the given individual's chromosome.
     * Each nurse is marked as a negative number, and all patients to the left of this number belong to the nurse.
     * This method randomizes the groups of nurses, so that the content within each of them remains the same.
     *
     * @param child The individual whose chromosome needs to be modified.
     */
    public void randomizeOrderOfNurses(Individual child) {
        int[] chromosome = child.getChromosome();
        int[] chromosomeCopy = new int[chromosome.length];

        Map<String, int[]> nurseMap = new HashMap<>();
        
        int count = 0;
        int start = 0;
        int end = 0;
        for (int i = 0; i < chromosome.length; i++) {
            if (chromosome[i] < 0) {
                end = i;
                int[] tempArray = new int[end - start];
                for (int j = 0; j < end - start; j++) {
                    tempArray[j] = chromosome[start + j];
                }
                nurseMap.put(String.valueOf(count), tempArray);

                start = end + 1;
                count++;
            }
        }

        List<Integer> list = new ArrayList<>();
        for (int i = 0; i < nbrNurses; i++) {
            list.add(i);
        }
        Collections.shuffle(list);
        int[] indices = list.stream().mapToInt(Integer::intValue).toArray();

        int index = 0;
        for (int i = 0; i < nbrNurses; i++) {
            int[] tempArray = nurseMap.get(String.valueOf(indices[i]));
            for (int j = 0; j < tempArray.length; j++) {
                chromosomeCopy[index] = tempArray[j];
                index++;
            }
            chromosomeCopy[index] = -i - 1;
            index++;
        }

        child.setChromosome(chromosomeCopy);
    }

    /**
     * Performs mutation on the given individual by randomly modifying its chromosome.
     * Mutation involves swapping and moving elements within and between groups, as well as reversing a portion of the chromosome.
     * The probability of each mutation operation is controlled by the specified probabilities.
     *
     * @param child the individual to be mutated
     */
    public void mutation(Individual child) {
        int[] chromosome = child.getChromosome();

        if (Math.random() < 0.65) { // randomly swap two elements in the same group
            intraSwapMutation(child, chromosome);
        }
        if (Math.random() < 0.65) { // randomly move an element in the same group
            intraMoveMutation(child, chromosome);
        }
        if (Math.random() < 0.65) { // randomly swap two elements in the chromosome
            interSwapMutation(child, chromosome);
        }
        if (Math.random() < 0.65) { // randomly move an element from one group to another
            interMoveMutation(child, chromosome);
        }
        if (Math.random() < 0.1) { // pick two random indices and reverse the elements between them
            reverseMutation(child, chromosome);
        }
    }
    
    /**
     * Performs an intra-group swap mutation on the given child individual's chromosome.
     * An intra-group swap mutation randomly selects a group of patients in the chromosome and swaps two random patients within the group.
     * The group is determined by selecting a random index in the chromosome where a negative value is encountered.
     * If the selected group has more than two patients, two random patients within the group are swapped.
     * 
     * @param child The individual for which the mutation is performed.
     * @param chromosome The chromosome of the child individual.
     */
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

    /**
     * Performs an intra-move mutation on the given child individual's chromosome.
     * An intra-move mutation randomly selects a group of elements in the chromosome and
     * moves one element within the group to a different position.
     *
     * @param child      the individual to mutate
     * @param chromosome the chromosome of the individual
     */
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
            // find two random indices in the group and move the element at index1 to index2
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

    /**
     * Performs an inter-swap mutation on the given child individual.
     * Randomly selects two indices in the chromosome and swaps the values at those indices.
     * Updates the chromosome of the child individual with the mutated chromosome.
     *
     * @param child      The individual to be mutated.
     * @param chromosome The chromosome of the individual to be mutated.
     */
    public void interSwapMutation(Individual child, int[] chromosome) {
        int index1 = (int) (Math.random() * (chromosome.length - 1));
        int index2 = (int) (Math.random() * (chromosome.length - 1));

        int temp = chromosome[index1];
        chromosome[index1] = chromosome[index2];
        chromosome[index2] = temp;
        child.setChromosome(chromosome);
    }

    /**
     * Performs an inter-move mutation on the given child individual's chromosome.
     * Inter-move mutation picks two indices in the chromosome and moves the element at index1 to index2.
     *
     * @param child      the individual to perform the mutation on
     * @param chromosome the chromosome of the individual
     */
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

    /**
     * Performs a reverse mutation on the given child individual's chromosome.
     * The reverse mutation swaps a random subsequence of the chromosome in-place.
     *
     * @param child      the individual to perform the reverse mutation on
     * @param chromosome the chromosome of the child individual
     */
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

    /**
     * Performs survivor selection by replacing the two individuals in the population with the worst fitness
     * with the given child individuals, if they have better fitness.
     *
     * @param child1 The first child individual.
     * @param child2 The second child individual.
     */
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
            // randomizeOrderOfNurses(child2);
            population[worstFitnessIndex1] = child2;
        } else if (child2.getFitness() < worstFitness2) {
            // randomizeOrderOfNurses(child2);
            population[worstFitnessIndex2] = child2;
        }

        if (child1.getFitness() < worstFitness2) {
            // randomizeOrderOfNurses(child1);
            population[worstFitnessIndex2] = child1;
        } 
    }

    /**
     * Prints the solution represented by the given best individual.
     * The solution is a double array, where each subarray represents the route of a nurse.
     * This is printed to console in the format [[nurse1_route], [nurse2_route], ...].
     *
     * @param bestIndividual the best individual representing the solution
     */
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

    /**
     * Saves the best individual and travel time to a JSON file.
     * 
     * @param bestIndividual The best individual found by the genetic algorithm.
     * @param travelTime The travel time of the best individual's solution.
     */
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
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Prints the solution in detail for a given individual.
     * 
     * @param individual The individual for which the solution details are printed.
     */
    public void printSolutionInDetail(Individual individual) {
        DecimalFormat df = new DecimalFormat("#.##");
        
        System.out.println("Nurse capacity: " + capacityNurse);
        System.out.println("Depot return time: " + depot.getReturnTime());
        System.out.println("-----------------------------------");
        
        Object[] sequenceAndDuration = registerPatientSequence(individual);
        String[] patientSequence = (String[]) sequenceAndDuration[0];
        String[] routeDuration = (String[]) sequenceAndDuration[1];
        double totalDuration = (double) sequenceAndDuration[2];

        int[] totalGroupwiseDemand = totalGroupwiseDemand(individual);

        for (int i = 0; i < patientSequence.length; i++) {
            System.out.println("Nurse " + (i + 1) + " " + routeDuration[i] + " " + totalGroupwiseDemand[i] + " " + patientSequence[i]);
        }
        System.out.println("-----------------------------------");
        System.out.println("Objective value (total duration): " + df.format(totalDuration));
        System.out.println("Total travel time: " + df.format(individual.getFitness()));
    }

    /**
     * Helper method to printSolutionInDetail that registers the patient sequence for a given individual.
     * 
     * @param individual The individual for which to register the patient sequence.
     * @return An array containing the patient sequence, route durations, and total duration.
     */
    public Object[] registerPatientSequence(Individual individual) {
        int[] individualChromosome = individual.getChromosome().clone();

        String[] patientSequence = new String[nbrNurses];
        String[] routeDuration = new String[nbrNurses];
        double totalDuration = 0;

        DecimalFormat df = new DecimalFormat("#.##");

        // define a string builder
        StringBuilder solution = new StringBuilder();
        solution.append("D (0)");

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

        int start;
        
        // if the first element is a patient, then the nurse has to travel from the depot to the patient
        if (previous != 0) {
            time += travelTimes[0][previous];

            int careTime = patients.get(String.valueOf(previous)).getCareTime();
            if (time < patients.get(String.valueOf(previous)).getStartTime()) {
                time = patients.get(String.valueOf(previous)).getStartTime();
            }
            solution.append(" -> " + previous + " (" + df.format(time) + "-");
            time += careTime;
            solution.append(df.format(time) + ") [" + patients.get(String.valueOf(previous)).getStartTime() + "-" + patients.get(String.valueOf(previous)).getEndTime() + "]");
        
            start = 1;
        } else {
            start = 0;
        }

        // for each pair of patients, the nurse has to travel the distance between them
        // a nurse is denoted by a 0 in the individual, so when we encounter a 0 this will
        // indicate that the nurse has to travel to depot, and a new nurse will start from the depot
        for (int i = start; i < individualChromosome.length; i++) {
            current = individualChromosome[i];

            time += travelTimes[previous][current];

            if (current != 0) {
                int careTime = patients.get(String.valueOf(current)).getCareTime();
                if (time < patients.get(String.valueOf(current)).getStartTime()) {
                    time = patients.get(String.valueOf(current)).getStartTime();
                }
                solution.append(" -> " + current + " (" + df.format(time) + "-");
                time += careTime;
                solution.append(df.format(time) + ") [" + patients.get(String.valueOf(current)).getStartTime() + "-" + patients.get(String.valueOf(current)).getEndTime() + "]");    
            } else {
                solution.append(" -> D (" + df.format(time) + ")");
                patientSequence[group] = solution.toString();
                routeDuration[group] = df.format(time);

                solution = new StringBuilder();
                solution.append("D (0)");
                totalDuration += time;
                group++;
                time = 0;
            }

            previous = current;
        }
        
        Object[] result = new Object[3];
        result[0] = patientSequence;
        result[1] = routeDuration;
        result[2] = totalDuration;
        return result;
    }

    /**
     * Tries every combination of putting two nurse rouse together as one.
     * The new route is ordered on end-time, and undergoes survivor selection.
     * 
     * @param individual The individual to combine routes for.
     */
    public void combineRoutes(Individual individual) {
        // try every combination of two nurses
        Map<String, int[]> nurseRoutes = new HashMap<>();

        int index = 0;
        int start = 0;
        int end = 0;
        for (int i = 0; i < individual.getChromosome().length; i++) {
            if (individual.getChromosome()[i] < 0) {
                end = i;
                int[] tempArray = new int[end - start];
                for (int j = 0; j < end - start; j++) {
                    tempArray[j] = individual.getChromosome()[start + j];
                }
                nurseRoutes.put(String.valueOf(index), tempArray);

                start = end + 1;
                index++;
            }
        }

        for (int i = 0; i < nbrNurses; i++) {
            for (int j = i + 1; j < nbrNurses; j++) {
                // get int[] at map index i and j, and combine them into one new int
                int[] tempArray = new int[nurseRoutes.get(String.valueOf(i)).length + nurseRoutes.get(String.valueOf(j)).length];
                for (int k = 0; k < nurseRoutes.get(String.valueOf(i)).length; k++) {
                    tempArray[k] = nurseRoutes.get(String.valueOf(i))[k];
                }
                for (int k = 0; k < nurseRoutes.get(String.valueOf(j)).length; k++) {
                    tempArray[k + nurseRoutes.get(String.valueOf(i)).length] = nurseRoutes.get(String.valueOf(j))[k];
                }

                // sort the content of int[] tempArray on patients.get(String.valueOf(int)).getEndTime()
                Integer[] tempArray2 = Arrays.stream(tempArray).boxed().toArray(Integer[]::new);
                Arrays.sort(tempArray2, Comparator.comparing(o -> patients.get(String.valueOf(o)).getEndTime()));
                tempArray = Arrays.stream(tempArray2).mapToInt(Integer::intValue).toArray();

                // create a new chromosome with the combined int[] tempArray
                int[] chromosome = new int[individual.getChromosome().length];
                for (int k = 0; k < tempArray.length; k++) {
                    chromosome[k] = tempArray[k];
                }

                int group = 0;
                for (int k = tempArray.length; k < chromosome.length; k++) {
                    if (group != i && group != j) {
                        for (int l = 0; l < nurseRoutes.get(String.valueOf(group)).length; l++) {
                            chromosome[k] = nurseRoutes.get(String.valueOf(group))[l];
                            k++;
                        }
                    }
                    chromosome[k] = -group - 1;
                    group++;
                }

                Individual newIndividual = new Individual(chromosome, 0);
                calculateFitness(newIndividual);
                survivorSelection(newIndividual, newIndividual);
            }
        }
    }

    /**
     * Helper method to divideRoutes.
     * Generates all combinations of dividing the patients of a nurse between two nurses.
     * Each new combination results in a individual that may undergo survivor selection.
     * Maximum 1 million different combinations, for time concerns.
     * 
     * @param arr The array of patients to divide between two nurses.
     * @param list1 The list of patients for the first nurse.
     * @param list2 The list of patients for the second nurse.
     * @param index The current index in the array.
     * @param nurse1 The first nurse.
     * @param nurse2 The second nurse.
     * @param otherNurses The other nurses in the route.
     * @param otherRoutes The routes of the other nurses.
     */
    public void generateCombinations(int[] arr, List<Integer> list1, List<Integer> list2, int index, int nurse1, int nurse2, List<Integer> otherNurses, List<List<Integer>> otherRoutes) {
        if (counting > 1000000) {
            return;
        }
        counting++;

        if (index == arr.length) {
            List<Integer> newRoute = new ArrayList<>();
            for (int i = 0; i < otherNurses.size(); i++) {
                newRoute.addAll(otherRoutes.get(i));
                newRoute.add(otherNurses.get(i));
            }
            newRoute.addAll(list1);
            newRoute.add(nurse1);
            newRoute.addAll(list2);
            newRoute.add(nurse2);
            for (int i = 1; i <= 25; i++) {
                if (!newRoute.contains(-i)) {
                    newRoute.add(-i);
                }
            }

            int[] newRouteInt = newRoute.stream().mapToInt(i -> i).toArray();
            Individual ind = new Individual(newRouteInt, 0);
            calculateFitness(ind);
            survivorSelection(ind, ind);

            return;
        }

        list1.add(arr[index]);
        generateCombinations(arr, list1, list2, index + 1, nurse1, nurse2, otherNurses, otherRoutes);
        list1.remove(list1.size() - 1);

        if (!list1.isEmpty()) {
            list2.add(arr[index]);
            generateCombinations(arr, list1, list2, index + 1, nurse1, nurse2, otherNurses, otherRoutes);
            list2.remove(list2.size() - 1);
        }
    }

    /**
     * Iterates through the nurses of route, trying all combinations of dividing the patients of each nurse
     * between two nurses. 
     * 
     * @param route The route (chromosome) to make new individuals of.
     */
    public void divideRoutes(int[] route) {
        List<Integer> emptyNurses = new ArrayList<>();
        Map<Integer, List<Integer>> nurseRoutes = new HashMap<>();

        List<Integer> currentRoute = new ArrayList<>();
        int currentNurse = 0;
        for (int i = 0; i < route.length; i++) {
            if (route[i] < 0) {
                if (currentNurse != 0) {
                    if (!currentRoute.isEmpty()) {
                        nurseRoutes.put(currentNurse, new ArrayList<>(currentRoute));
                        currentRoute.clear();
                    } else {
                        emptyNurses.add(currentNurse);
                    }
                }
                currentNurse = route[i];
            } else {
                currentRoute.add(route[i]);
            }
        }
        if (!currentRoute.isEmpty()) {
            nurseRoutes.put(currentNurse, new ArrayList<>(currentRoute));
        } else {
            emptyNurses.add(currentNurse);
        }

        for (Map.Entry<Integer, List<Integer>> entry : nurseRoutes.entrySet()) {
            if (!emptyNurses.isEmpty()) {
                Integer emptyNurse = emptyNurses.remove(0); // Use and remove the first empty nurse
                List<Integer> otherNurses = new ArrayList<>(nurseRoutes.keySet());
                otherNurses.remove(entry.getKey());
                List<List<Integer>> otherRoutes = otherNurses.stream().map(nurse -> nurseRoutes.get(nurse)).collect(Collectors.toList());
                counting++;
                generateCombinations(entry.getValue().stream().mapToInt(i -> i).toArray(), new ArrayList<>(), new ArrayList<>(), 0, entry.getKey(), emptyNurse, otherNurses, otherRoutes);
                emptyNurses.add(emptyNurse); // Put the empty nurse back
            }
        }
    }

    public Individual[] getPopulation() {
        return population;
    }

    public void setPopulation(Individual[] population) {
        this.population = population;
    }
    public static void main(String[] args) {
        boolean replacePopulation = false; // true=replace population if stale, false=no action
        boolean replaceGreedy = false; // true=replace greedily, false=replace with longest travel cut
        boolean parentClone = true; // true=possibly clone, false=only crossover
        boolean tournamentParentSelection = true; // true=tournament, false=best
        boolean orderOneCrossover = false; // true=order one, false=partially mapped
        int populationGeneration = 2; // 0=random, 1=greedy, 2=heuristic
        boolean diversity = false; // true=diversity in clutsering, false=no diversity // another for multi
        boolean shuffleNurses = false; // true=order of nurses are shuffled, false=no shuffle // another for multi

        int populationSize = 150; //600; // numberOfClusters * numberOfIndividualsInClusters * 4 // another for multi
        String filePath = "src/main/resources/train/train_9.json";
        int numberOfRuns = 1;
        int numberOfGenerations = 1500000;
        int numberOfClusters = 15;
        int numberOfIndividualsInClusters = 10; // will be timed by 4 for diversity

        GeneticAlgorithm GA = new GeneticAlgorithm(populationSize, filePath, tournamentParentSelection, orderOneCrossover, null, numberOfGenerations, replacePopulation, replaceGreedy, parentClone, populationGeneration, numberOfClusters, numberOfIndividualsInClusters, diversity, shuffleNurses);
        


        // GA.generatePopulation();
        // Individual bestTotalIndividual = GA.getPopulation()[0];

        // for (int i = 0; i < numberOfRuns; i++) {
        //     GA.runGA();

        //     Individual[] population = GA.getPopulation();
        //     Individual bestIndividual = population[0];
        //     double bestFitnessAfter = population[0].getFitness();

        //     for (Individual individual : population) {
        //         if (individual.getFitness() < bestFitnessAfter) {
        //             bestIndividual = individual;
        //             bestFitnessAfter = individual.getFitness();
        //         }
        //     }
        //     System.out.println("Best fitness of generation " + i + ": " + bestFitnessAfter);
        //     if (bestFitnessAfter < bestTotalIndividual.getFitness()) {
        //         bestTotalIndividual = bestIndividual;
        //     }
        // }

        // GA.saveSolution(bestTotalIndividual, bestTotalIndividual.getFitness());
        // GA.printSolution(bestTotalIndividual);
        // GA.printSolutionInDetail(bestTotalIndividual);



        ExecutorService executor = Executors.newFixedThreadPool(3); // Create a thread pool with 3 threads
        List<Future<Individual>> futures = new ArrayList<>();

        int[] popSizes = new int[] {150, 150, 600};
        boolean[] diversities = new boolean[] {false, false, true};
        boolean[] shuffleNurseses = new boolean[] {false, true, false};

        for (int i = 0; i < 3; i++) {
            GeneticAlgorithm ga = new GeneticAlgorithm(popSizes[i], filePath, tournamentParentSelection, orderOneCrossover, null, numberOfGenerations, replacePopulation, replaceGreedy, parentClone, populationGeneration, numberOfClusters, numberOfIndividualsInClusters, diversities[i], shuffleNurseses[i]); // Create a new GA instance
            Future<Individual> future = executor.submit(ga); // Submit it to be run by the thread pool
            futures.add(future);
        }

        Individual bestIndividual = null;
        for (Future<Individual> future : futures) {
            try {
                Individual ind = future.get(); // Get the result of the GA run
                if (bestIndividual == null || ind.getFitness() < bestIndividual.getFitness()) {
                    bestIndividual = ind; // Update the best individual if necessary
                }
            } catch (InterruptedException | ExecutionException e) {
                e.printStackTrace();
            }
        }

        executor.shutdown(); 

        GA.saveSolution(bestIndividual, bestIndividual.getFitness());
        GA.printSolution(bestIndividual);
        GA.printSolutionInDetail(bestIndividual);
    }
}
