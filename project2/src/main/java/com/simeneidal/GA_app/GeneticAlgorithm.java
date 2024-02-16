package com.simeneidal.GA_app;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.ArrayList;

import com.google.gson.Gson;
import java.io.FileReader;

public class GeneticAlgorithm {

    private static int populationSize = 10;
    private static String filePath = "src/main/resources/train/train_0.json";

    private static int nbrNurses;
    private static Map<String, JsonData.Patient> patients;
    private static double[][] travelTimes;

    public static void main(String[] args) {
        
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

        int[][] population = generatePopulation(nbrNurses, nbrPatients, populationSize);

        Object[] result = calculateTimeAndTimeRequirement(population[0]);
        double totalTravelTime = (double) result[0];
        boolean[] meetTimeRequirement = (boolean[]) result[1];

        int[] groupwiseDemand = totalGroupwiseDemand(population[0]);
        
        for (int demand : groupwiseDemand) {
            System.out.println(demand);
        }
        // System.out.println(totalTravelTime);
        // for (boolean b : meetTimeRequirement) {
        //     System.out.println(b);
        // }
    }
    
    public static int[][] generatePopulation(int nbrNurses, int nbrPatients, int populationSize) {
        int[][] population = new int[populationSize][nbrNurses + nbrPatients];
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
            population[i] = list.stream().mapToInt(Integer::intValue).toArray();
        }
        return population;
    }
    
    public static int fitness(int[] individual) {
        // calculate the fitness of the individual
        // each individual consists of a list of patients and nurses
        // each nurse travels to each patient in their group
        // each nurse has a capacity and the patients have a demand, total demand cannot exceed capacity for a nurse
        // each patient has a time window, and the nurse must arrive within the time window, else nurse has to wait
        // each patient has a service time, and the nurse must spend that much time with the patient
        // the fitness is the total travel time (only) taken by all nurses to visit all patients

        return 0;
    }

    public static double totalTravelTime(int[] individual) {
        // calculate the total travel time of the individual
        int previous = individual[0];
        int current;
        double totalTravelTime = 0;
        
        // if the first element is a patient, then the nurse has to travel from the depot to the patient
        if (previous != 0) {
            totalTravelTime += travelTimes[0][previous];
            // System.out.println("Depot to " + previous + " " + travelTimes[0][previous] + " " + totalTravelTime);
        }

        // for each pair of patients, the nurse has to travel the distance between them
        // a nurse is denoted by a 0 in the individual, so when we encounter a 0 this will
        // indicate that the nurse has to travel to depot, and a new nurse will start from the depot
        for (int i = 1; i < individual.length; i++) {
            current = individual[i];
            totalTravelTime += travelTimes[previous][current];
            // System.out.println(previous + " " + current + " " + travelTimes[previous][current] + " " + totalTravelTime);
            previous = current;
        }
        return totalTravelTime;
    }

    public static int[] totalGroupwiseDemand(int[] individual) {
        // calculate the total demand of each group of patients
        // each nurse has a capacity and the patients have a demand, 
        // total demand cannot exceed capacity for a nurse (checked later)

        int[] groupwiseDemand = new int[nbrNurses];

        int total = 0;

        int demand = 0;
        int group = 0;
        for (int patient : individual) {
            if (patient == 0) {
                groupwiseDemand[group] = demand;
                group++;
                demand = 0;
            } else {
                demand += patients.get(String.valueOf(patient)).getDemand();
                total += patients.get(String.valueOf(patient)).getDemand();
            }
        }

        return groupwiseDemand;
    }

    public static Object[] calculateTimeAndTimeRequirement(int[] individual) {
        // check if the individual meets the time requirement of each patient
        // each patient has a time window, and the nurse must arrive within the time window, else nurse has to wait
        // each patient has a service time, and the nurse must spend that much time with the patient

        boolean[] timeDemand = new boolean[nbrNurses];
        for (int i = 0; i < nbrNurses; i++) {
            timeDemand[i] = true;
        }

        int group = 0;
        double time = 0;

        int previous = individual[0];
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
        for (int i = 1; i < individual.length; i++) {
            current = individual[i];

            totalTravelTime += travelTimes[previous][current];
            time += travelTimes[previous][current];

            if (current != 0) {
                int careTime = patients.get(String.valueOf(current)).getCareTime();
                int endTime = patients.get(String.valueOf(current)).getEndTime();
                if (time > endTime || time + careTime > endTime) {
                    timeDemand[group] = false;
                } else if (time < patients.get(String.valueOf(current)).getStartTime()) {
                    time = patients.get(String.valueOf(current)).getStartTime();
                }
                time += careTime;
            } else {
                group++;
            }

            previous = current;
        }
        
        Object[] result = new Object[2];
        result[0] = totalTravelTime;
        result[1] = timeDemand; 
        return result;
    }

    public static void crossover() {
        System.out.println("Crossover");
    }
    
    public static void mutation() {
        System.out.println("Mutation");
    }
    
    public static void survivorSelection() {
        System.out.println("Selection");
    }
    
    
    public static void termination() {
        System.out.println("Termination");
    }
    
    public static void initialization() {
        System.out.println("Initialization");
    }
    
    public static void evaluation() {
        System.out.println("Evaluation");
    }
    
    public static void replacement() {
        System.out.println("Replacement");
    }
    
    public static void reproduction() {
        System.out.println("Reproduction");
    }
    
}
