package org.calibrationframework.fouriermethod.calibration.constraints;

/**
 * Absence of constraints.
 * 
 * @author Alessandro Gnoatto
 *
 */
public class Unconstrained extends BoundConstraint{

	public Unconstrained() {
		super(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
	}

	@Override
	public double applyConstraint(double parameterToTest) {
		return parameterToTest;
	}	

}
