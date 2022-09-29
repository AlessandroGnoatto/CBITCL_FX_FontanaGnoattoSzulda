/*
 * (c) Copyright Christian P. Fries, Germany. Contact: email@christian-fries.de.
 *
 * Created on 20.05.2005
 */
package org.calibrationframework.montecarlo.models;

import java.util.Map;

import org.calibrationframework.modelling.ModelInterface;
import net.finmath.stochastic.RandomVariable;
import net.finmath.time.TimeDiscretization;

/**
 * The interface implemented by a simulation of an SDE.
 * Provides the dimension of the SDE and the the time discretization of the
 * simulation.
 *
 * @author Christian Fries
 */
public interface MonteCarloSimulationInterface extends ModelInterface {

	/**
	 * Returns the numberOfPaths.
	 *
	 * @return Returns the numberOfPaths.
	 */
	int getNumberOfPaths();

	/**
	 * Returns the timeDiscretization.
	 *
	 * @return Returns the timeDiscretization.
	 */
	TimeDiscretization getTimeDiscretization();

	/**
	 * Returns the time for a given time index.
	 *
	 * @param timeIndex Time index
	 * @return Returns the time for a given time index.
	 */
	double getTime(int timeIndex);

	/**
	 * Returns the time index for a given time.
	 *
	 * @param time The time.
	 * @return Returns the time index for a given time.
	 */
	int getTimeIndex(double time);

	/**
	 * Returns a random variable which is initialized to a constant,
	 * but has exactly the same number of paths or discretization points as the ones used by this <code>MonteCarloSimulationInterface</code>.
	 *
	 * @param value The constant value to be used for initialized the random variable.
	 * @return A new random variable.
	 */
	RandomVariable getRandomVariableForConstant(double value);

	/**
	 * This method returns the weights of a weighted Monte Carlo method (the probability density).
	 *
	 * @param timeIndex Time index at which the process should be observed
	 * @return A vector of positive weights which sums up to one
	 */
	RandomVariable getMonteCarloWeights(int timeIndex);

	/**
	 * This method returns the weights of a weighted Monte Carlo method (the probability density).
	 *
	 * @param time Time at which the process should be observed
	 * @return A vector of positive weights which sums up to one
	 */
	RandomVariable getMonteCarloWeights(double time);

	/**
	 * Create a clone of this simulation modifying some of its properties (if any).
	 *
	 * @param dataModified The data which should be changed in the new model
	 * @return Returns a clone of this model, with some data modified (then it is no longer a clone :-)
	 */
	MonteCarloSimulationInterface getCloneWithModifiedData(Map<String, Object> dataModified);
}
