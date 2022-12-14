package org.calibrationframework.fouriermethod.calibration.constraints;

/**
 * Positivity constraint for calibration parameters
 * @author Alessandro Gnoatto
 *
 */
public class PositivityConstraint extends BoundConstraint {

	public PositivityConstraint() {
		super(0.0, Double.POSITIVE_INFINITY);
	}
	
	@Override
	public double applyConstraint(double parameterToTest) {
		return Math.abs(parameterToTest);
	}

}
