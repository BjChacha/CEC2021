package etmo.metaheuristics.matmy3.models;

public abstract class AbstractDistribution {
    public abstract double[] sample();

    public abstract double[][] sample(int count); 
}
