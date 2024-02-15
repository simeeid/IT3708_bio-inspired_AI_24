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

    private static double[][] travelTimes;

    public static void main(String[] args) {
        
        String instanceName = "";
        int nbrNurses = 0;
        int capacityNurse = 0;
        double benchmark = 0;
        JsonData.Depot depot = null;
        Map<String, JsonData.Patient> patients = null;
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
        // for (int i = 0; i < 1; i++) {
        //     for (int j = 0; j < nbrNurses + nbrPatients; j++) {
        //         System.out.print(population[i][j] + " ");
        //     }
        //     System.out.println();
        // }
        double travelTime = totalTravelTime(population[0]);
        System.out.println();
        System.out.println("Total time: " + travelTime);
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
            System.out.println("Depot to " + previous + " " + travelTimes[0][previous] + " " + totalTravelTime);
        }

        // for each pair of patients, the nurse has to travel the distance between them
        // a nurse is denoted by a 0 in the individual, so when we encounter a 0 this will
        // indicate that the nurse has to travel to depot, and a new nurse will start from the depot
        for (int i = 1; i < individual.length; i++) {
            current = individual[i];
            totalTravelTime += travelTimes[previous][current];
            System.out.println(previous + " " + current + " " + travelTimes[previous][current] + " " + totalTravelTime);
            previous = current;
        }
        return totalTravelTime;
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
