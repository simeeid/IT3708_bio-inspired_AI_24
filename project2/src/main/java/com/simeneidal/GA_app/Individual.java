package com.simeneidal.GA_app;

/**
 * Represents an individual in a genetic algorithm population.
 */
public class Individual {
    private int[] chromosome;
    private double fitness;

    /**
     * Constructs a new individual with the given chromosome and fitness.
     *
     * @param chromosome the chromosome of the individual
     * @param fitness the fitness value of the individual
     */
    public Individual(int[] chromosome, double fitness) {
        this.chromosome = chromosome;
        this.fitness = fitness;
    }

    /**
     * Returns the chromosome of the individual.
     *
     * @return the chromosome of the individual
     */
    public int[] getChromosome() {
        return chromosome;
    }

    /**
     * Sets the chromosome of the individual.
     *
     * @param chromosome the new chromosome for the individual
     */
    public void setChromosome(int[] chromosome) {
        this.chromosome = chromosome;
    }

    /**
     * Returns the fitness value of the individual.
     *
     * @return the fitness value of the individual
     */
    public double getFitness() {
        return fitness;
    }

    /**
     * Sets the fitness value of the individual.
     *
     * @param fitness the new fitness value for the individual
     */
    public void setFitness(double fitness) {
        this.fitness = fitness;
    }
}

