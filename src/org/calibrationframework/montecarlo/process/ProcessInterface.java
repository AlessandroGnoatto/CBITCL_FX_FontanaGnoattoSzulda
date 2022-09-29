/*
 * (c) Copyright Christian P. Fries, Germany. Contact: email@christian-fries.de.
 *
 * Created on 29.02.2008
 */
package org.calibrationframework.montecarlo.process;

import net.finmath.stochastic.RandomVariable;
import net.finmath.time.TimeDiscretization;

/**
 * The interface for a stochastic process <i>X</i>.
 *
 * @author Christian Fries
 */
public interface ProcessInterface {

	/**
	 * This method returns the realization of a component of the process for a given time index.
	 *
	 * @param timeIndex Time index at which the process should be observed
	 * @param component Component index of the process
	 * @return The process component realizations (given as <code>RandomVariable</code>)
	 */
	RandomVariable getProcessValue(int timeIndex, int component);

	/**
	 * This method returns the weights of a weighted Monte Carlo method (the probability density).
	 *
	 * @param timeIndex Time index at which the process should be observed
	 * @return A vector of positive weights which sums up to one
	 */
	RandomVariable getMonteCarloWeights(int timeIndex);

	/**
	 * @return Returns the numberOfComponents.
	 */
	int getNumberOfComponents();

	/**
	 * @return Returns the timeDiscretization.
	 */
	TimeDiscretization getTimeDiscretization();

	/**
	 * @param timeIndex Time index.
	 * @return Returns the time for a given time index.
	 */
	double getTime(int timeIndex);

	/**
	 * Returns the time index for a given simulation time.
	 * @param time The given simulation time.
	 * @return Returns the time index for a given time
	 */
	int getTimeIndex(double time);

	/**
	 * Create and return a clone of this process. The clone is not tied to any model, but has the same
	 * process specification, that is, if the model is the same, it would generate the same paths.
	 *
	 * @return Clone of the process
	 */
	ProcessInterface clone();

}
