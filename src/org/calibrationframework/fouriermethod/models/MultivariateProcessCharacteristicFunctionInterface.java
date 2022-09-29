package org.calibrationframework.fouriermethod.models;

import org.calibrationframework.fouriermethod.CharacteristicFunctionInterface;
import org.calibrationframework.modelling.ModelInterface;

/**
 * Interface which has to be implemented by models providing the
 * characteristic function of some stochastic process among others,
 * associated to a name of the user's choice.
 * 
 */
@FunctionalInterface
public interface MultivariateProcessCharacteristicFunctionInterface extends ModelInterface {
	
	/**
	 * Returns the characteristic function of X(t), where X is <code>this</code> stochastic process.
	 * X is multivariate process and each component of X has a name that uniquely identifies it.
	 * 
	 * @param time The time at which the stochastic process is observed.
	 * @return The characteristic function of X(t).
	 */
	CharacteristicFunctionInterface apply(double time, String underlyingName);

}
