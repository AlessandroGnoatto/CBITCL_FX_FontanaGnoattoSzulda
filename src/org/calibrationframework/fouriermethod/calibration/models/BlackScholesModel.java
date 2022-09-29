/*
 * (c) Copyright Christian P. Fries, Germany. Contact: email@christian-fries.de.
 *
 * Created on 23.03.2014
 */

package org.calibrationframework.fouriermethod.calibration.models;

import org.apache.commons.math3.complex.Complex;

import org.calibrationframework.fouriermethod.CharacteristicFunctionInterface;
import org.calibrationframework.fouriermethod.calibration.constraints.*;
import org.calibrationframework.fouriermethod.models.ProcessCharacteristicFunctionInterface;

/**
 * Implements the characteristic function of a Black Scholes model.
 * 
 * @author Christian Fries
 * @author Alessandro Gnoatto
 */
public class BlackScholesModel implements CalibrableProcessInterface {

	private double initialValue;
	private double riskFreeRate;		// Actually the same as the drift (which is not stochastic)
	private double volatility;
	 
	private ScalarParameterInformationInterface volatilityInfo = new ScalarParameterInformation(true, new PositivityConstraint()); 
	
	private final double[] parameterLowerBounds = {0.0};
	private final double[] parameterUpperBounds = {1E6};

	public BlackScholesModel(double initialValue, double riskFreeRate, double volatility) {
		super();
		this.initialValue = initialValue;
		this.riskFreeRate = riskFreeRate;
		this.volatility = volatility;
	}
	
	public BlackScholesModel getCloneForModifiedParameters(double[] params) {	
		double newVol = this.volatilityInfo.getConstraint().applyConstraint(params[0]);
		return new BlackScholesModel(this.initialValue,this.riskFreeRate,newVol);
	}
	
	public double getInitialValue() {
		return this.initialValue;
	}
	
	public void setInitialValue(double newInitialValue) {
		this.initialValue = newInitialValue;
	}
	
	public double getVolatility() {
		return this.volatility;
	}
	
	public double getRiskFreeRate() {
		return this.riskFreeRate;
	}

	@Override
	public CharacteristicFunctionInterface apply(double time) {
		return argument -> {
			Complex iargument = argument.multiply(Complex.I);
			return	iargument
					.multiply(
							iargument
							.multiply(0.5*volatility*volatility*time)
							.add(Math.log(initialValue)-0.5*volatility*volatility*time+riskFreeRate*time))
					.add(-riskFreeRate*time)
					.exp();
		};
	}

	@Override
	public ProcessCharacteristicFunctionInterface getCharacteristiFunction() {
		return this;
	}

	@Override
	public double[] getParameterUpperBounds() {
		return this.parameterUpperBounds;
	}

	@Override
	public double[] getParameterLowerBounds() {
		return this.parameterLowerBounds;
	}

}
