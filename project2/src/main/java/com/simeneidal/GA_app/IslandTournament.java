package com.simeneidal.GA_app;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class IslandTournament {
    private Individual[][] islandPopulations;
    GeneticAlgorithm ga;

    // constructor
    public IslandTournament(int populationSize, String filepath, int nbrIslands) {
        islandPopulations = new Individual[nbrIslands][populationSize];
        ga = new GeneticAlgorithm(populationSize, filepath, null);

        for (int i = 0; i < nbrIslands; i++) {
            ga.generatePopulation();
            islandPopulations[i] = ga.getPopulation();
        }
    }

    // run island simulation
    public void runIslandSimulation(int nbrGenerations, int nbrIterations) {
        for (int i = 0; i < nbrIterations; i++) {
            for (int j = 0; j < islandPopulations.length; j++) {
                ga.runGA(islandPopulations[j], nbrGenerations);
                islandPopulations[j] = ga.getPopulation();
            }
            // migrate individuals between islands
            List<Integer> list = new ArrayList<>();
            for (int j = 0; j < islandPopulations.length; j++) {
                list.add(j);
            }
            Collections.shuffle(list);
            int[] indices = list.stream().mapToInt(Integer::intValue).toArray();

            for (int j = 0; j < islandPopulations.length; j++) {
                Individual bestIndividual;
                bestIndividual = ga.findOneBestIndividual(islandPopulations[indices[j]]);
                // find the worst individual
                Individual worstIndividual = islandPopulations[j][0];
                int worstIndex = 0;
                for (int k = 0; k < islandPopulations[j].length; k++) {
                    if (islandPopulations[j][k].getFitness() > worstIndividual.getFitness()) {
                        worstIndividual = islandPopulations[j][k];
                        worstIndex = k;
                    }
                }
                // replace the worst individual with the best individual
                islandPopulations[j][worstIndex] = bestIndividual;
            }
        }

        Individual[] bestIndividuals = new Individual[islandPopulations.length];
        for (int j = 0; j < islandPopulations.length; j++) {
            // find the best individual from each island
            bestIndividuals[j] = islandPopulations[j][0];
            for (int k = 0; k < islandPopulations[j].length; k++) {
                if (islandPopulations[j][k].getFitness() < bestIndividuals[j].getFitness()) {
                    bestIndividuals[j] = islandPopulations[j][k];
                }
            }
        }
        // print the best fitnesses
        for (int j = 0; j < bestIndividuals.length; j++) {
            // ga.calculateFitness(bestIndividuals[j], true);
            System.out.println("Island " + j + " best fitness: " + bestIndividuals[j].getFitness());
        }
        // print the best individual
        Individual bestIndividual = bestIndividuals[0];
        for (int j = 0; j < bestIndividuals.length; j++) {
            if (bestIndividuals[j].getFitness() < bestIndividual.getFitness()) {
                bestIndividual = bestIndividuals[j];
            }
        }
        System.out.println("Best individual fitness: " + bestIndividual.getFitness());
        // System.out.println("Best individual chromosome: " + Arrays.toString(bestIndividual.getChromosome()));
        
        ga.printSolution(bestIndividual);

        // System.out.println();
        // System.out.print("[[");
        // if (bestIndividual.getChromosome()[0] < 0) {
        //     System.out.print("], [");
        // } else {
        //     System.out.print(bestIndividual.getChromosome()[0]);
        // }
        // int previous = bestIndividual.getChromosome()[0];
        // for (int j = 1; j < bestIndividual.getChromosome().length - 1; j++) {
        //     if (bestIndividual.getChromosome()[j] < 0) {
        //         System.out.print("], [");
        //     } else if (previous > 0) {
        //         System.out.print(", " + bestIndividual.getChromosome()[j]);
        //     } else {
        //         System.out.print(bestIndividual.getChromosome()[j]);
        //     }
        //     previous = bestIndividual.getChromosome()[j];
        // }
        // System.out.println("]]");
        // System.out.println();
    }

    public static void main(String[] args) {
        IslandTournament it = new IslandTournament(10, "src/main/resources/train/train_9.json", 8);
        it.runIslandSimulation(5000, 10);
    }
}
