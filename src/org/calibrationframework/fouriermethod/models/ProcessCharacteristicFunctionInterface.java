/*
 * (c) Copyright Christian P. Fries, Germany. Contact: email@christian-fries.de.
 *
 * Created on 24.03.2014
 */

package org.calibrationframework.fouriermethod.models;

import org.calibrationframework.fouriermethod.CharacteristicFunctionInterface;
import org.calibrationframework.modelling.ModelInterface;

/**
 * Interface which has to be implemented by models providing the
 * characteristic function of the underlying stochastic process.
 * 
 * @author Christian Fries
 */
@FunctionalInterface
public interface ProcessCharacteristicFunctionInterface extends ModelInterface {

	/**
	 * Returns the characteristic function of X(t), where X is <code>this</code> stochastic process.
	 * 
	 * @param time The time at which the stochastic process is observed.
	 * @return The characteristic function of X(t).
	 */
	CharacteristicFunctionInterface apply(double time);
	
}
