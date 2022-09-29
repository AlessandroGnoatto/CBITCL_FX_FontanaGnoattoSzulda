package org.calibrationframework.montecarlo.process;

import org.calibrationframework.stochastic.*;

import net.finmath.time.TimeDiscretization;
import net.finmath.stochastic.RandomVariable;

/**
 * This interface has to be implemented by every class that is supposed to stand for a Monte Carlo simulation of a CBI process 
 * (or a family of CBI processes like the flow of tempered processes) along some time discretization by means of some pseudo random number generator of some seed.
 * The number of simulation is stocked inside the variable numberOfPaths,
 * and the sample paths are stocked in every increment of process value via the RandomVariable class.
 * @author Szulda Guillaume
 *
 */
public interface MonteCarloCBIProcessInterface extends ProcessInterface {
	
	public double getTimeHorizon();
	
	public int getNumberOfTimeSteps();
	
	public CBIProcessInterface getCBIProcess();
	
	public long getSeed();
	
	public RandomVariable getCBIProcessValue(double time, int factorIndex);
	
	public RandomVariable getCBIProcessValue(int timeIndex, int factorIndex);

	public int getNumberOfPaths();
	
	public MonteCarloCBIProcessInterface getCloneWithModifiedTimeDiscretization(TimeDiscretization newTimeDiscretization);
	
	public RandomVariable getRandomVariableForConstant(double value);
	
	public MonteCarloCBIProcessInterface getCloneWithModifiedSeed(int newSeed);
	
	public RandomVariable getIncrement(int timeIndex, int factor);
	
	public RandomVariable getMonteCarloWeights(double time);
	
}
