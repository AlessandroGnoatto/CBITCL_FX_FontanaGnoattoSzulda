package org.calibrationframework.fouriermethod.products;

import java.util.Map;

import org.apache.commons.math3.complex.Complex;

import org.calibrationframework.exception.CalculationException;
import org.calibrationframework.fouriermethod.models.MultivariateProcessCharacteristicFunctionInterface;

public abstract class EuropeanOptionSmileMultiAsset implements SmileByIntegralTransformMultiAsset {
	
	private final String underlyingName;
	private final double maturity;
	private final double[] strikes;
	
	public EuropeanOptionSmileMultiAsset(String underlyingName, double maturity, double[] strikes) {
		super();
		this.underlyingName = underlyingName;
		this.maturity = maturity;
		this.strikes = strikes;
	}
	
	public EuropeanOptionSmileMultiAsset(double maturity, double[] strikes) {
		this(null,maturity,strikes);      
	}
	
	@Override
	public double getMaturity() {
		return this.maturity;
	}
	
	public double[] getStrikes() {
		return this.strikes;
	}
	
	public String getUnderlyingName() {
		return this.underlyingName;
	}

	@Override
	public double getIntegrationDomainImagLowerBound() {
		return -1;
	}

	@Override
	public double getIntegrationDomainImagUpperBound() {
		return 0;
	}

	public abstract Map<Double, Double> getValue(MultivariateProcessCharacteristicFunctionInterface model) throws CalculationException;
	
	/**
	 * This method allows us to reuse the same pricer (same pricing algorithm) over different option smiles.
	 * @param underlyingName
	 * @param maturity
	 * @param strikes
	 * @return the same pricer now referring to a different smile.
	 */
	public abstract EuropeanOptionSmileMultiAsset getCloneWithModifiedParameters(String underlyingName, double maturity, double[] strikes);
	

	@Override
	public Complex apply(Complex z) {
		return ((z.subtract(Complex.I)).multiply(z)).negate();
	}

}
