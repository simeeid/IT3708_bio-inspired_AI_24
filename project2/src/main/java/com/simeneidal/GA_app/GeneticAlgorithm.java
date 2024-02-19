package com.simeneidal.GA_app;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.ArrayList;
import java.util.Arrays;

import com.google.gson.Gson;
import java.io.FileReader;

public class GeneticAlgorithm {

    private int populationSize;
    private String filePath;

    private int nbrNurses;
    private Map<String, JsonData.Patient> patients;
    private double[][] travelTimes;
    private Individual[] population;

    public GeneticAlgorithm(int populationSize, String filePath) {
        this.populationSize = populationSize;
        this.filePath = filePath;
        this.population = null;

        String instanceName = "";
        nbrNurses = 0;
        int capacityNurse = 0;
        double benchmark = 0;
        JsonData.Depot depot = null;
        patients = null;
        travelTimes = null;
        int nbrPatients = 0;

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

        generatePopulation(nbrNurses, nbrPatients, populationSize);
        Individual[] parents = parentSelection();
        crossover(parents[0].getChromosome(), parents[1].getChromosome());
    }
    
    public void generatePopulation(int nbrNurses, int nbrPatients, int populationSize) {
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
                list.add(0);
            }
            Collections.shuffle(list);
            list.add(0); // add a 0 to the end of the list
            population[i].setChromosome(list.stream().mapToInt(Integer::intValue).toArray());
        }

        for (Individual individual : population) {
            calculateFitness(individual);
        }
    }
    
    public void calculateFitness(Individual individual) {
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
            if (demand > nbrNurses) {
                numberOfNursesOverCapacity++;
            }
        }

        // for (int demand : groupwiseDemand) {
        //     System.out.println(demand);
        // }
        // System.out.println(totalTravelTime);
        // System.out.println(numberOfPatientsMissed);
        // System.out.println(numberOfNursesOverCapacity);

        individual.setFitness(totalTravelTime + numberOfPatientsMissed * 50 + numberOfNursesOverCapacity * 50);
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
            if (patient == 0) {
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

        int[] individualChromosome = individual.getChromosome();
        int numberOfPatientsMissed = 0;

        int group = 0;
        double time = 0;

        int previous = individualChromosome[0];
        int current;
        double totalTravelTime = 0;
        
        // if the first element is a patient, then the nurse has to travel from the depot to the patient
        if (previous != 0) {
            totalTravelTime += travelTimes[0][previous];
            time += travelTimes[0][previous];
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
        // find the two individuals in the population with the best fitness
        // these two individuals will be the parents for the next generation
        Individual parent1 = population[0];
        Individual parent2 = population[1];
        for (Individual individual : population) {
            if (individual.getFitness() < parent1.getFitness()) {
                parent2 = parent1;
                parent1 = individual;
            } else if (individual.getFitness() < parent2.getFitness()) {
                parent2 = individual;
            }
        }
        return new Individual[] {parent1, parent2};
    }

    public void crossover(int[] parent1, int[] parent2) {
        if (true) { // order 1 crossover
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
                child1[i] = -1;
                child2[i] = -1;
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

            child1 = orderOneCrossover(child1, parent2, cuttingPoint1, cuttingPoint2, parentIndex, childIndex, chromosomeLength);
            child2 = orderOneCrossover(child2, parent1, cuttingPoint1, cuttingPoint2, parentIndex, childIndex, chromosomeLength);

            // add the last nurse back to the children
            child1 = Arrays.copyOf(child1, child1.length + 1);
            child2 = Arrays.copyOf(child2, child2.length + 1);
            child1[child1.length - 1] = 0;
            child2[child2.length - 1] = 0;
        }
    }

    public int[] orderOneCrossover(int[] child1, int[] parent2, int cuttingPoint1, int cuttingPoint2, int parentIndex, int childIndex, int chromosomeLength) {
        while (true) {
            boolean flag = false;
            int nbrZeroes = 0;

            if (parent2[parentIndex] == 0) {
                for (int i = 0; i < chromosomeLength; i++) {
                    if (child1[i] == 0) {
                        nbrZeroes++;
                        if (nbrZeroes == nbrNurses - 1) {
                            flag = true;
                            break;
                        }
                    }
                }
            } else {
                for (int i = 0; i < chromosomeLength; i++) {
                    if (child1[i] == parent2[parentIndex]) {
                        flag = true;
                        break;
                    }
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
        // // print the number of 0s in the child
        // int nbrZeroes = 0;
        // for (int i = 0; i < chromosomeLength; i++) {
        //     if (child1[i] == 0) {
        //         nbrZeroes++;
        //     }
        // }
        // System.out.println("Number of 0s in the child: " + nbrZeroes);

        // // check that all numbers 1 to n are present in the child
        // boolean noNumMissing = true;
        // for (int i = 1; i < 101; i++) {
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
    
    public static void mutation() {
        System.out.println("Mutation");
    }
    
    public void survivorSelection(Individual child1, Individual child2) {
        System.out.println("Survivor selection");
    }

    public static void main(String[] args) {
        GeneticAlgorithm GA = new GeneticAlgorithm(10, "src/main/resources/train/train_0.json");
        GA.parentSelection();
    }
}

