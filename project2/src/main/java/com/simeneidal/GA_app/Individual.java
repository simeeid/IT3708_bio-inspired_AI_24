package com.simeneidal.GA_app;

public class Individual {
    private int[] chromosome;
    private double fitness;

    public Individual(int[] chromosome, double fitness) {
        this.chromosome = chromosome;
        this.fitness = fitness;
    }

    public int[] getChromosome() {
        return chromosome;
    }

    public void setChromosome(int[] chromosome) {
        this.chromosome = chromosome;
    }

    public double getFitness() {
        return fitness;
    }

    public void setFitness(double fitness) {
        this.fitness = fitness;
    }
}

