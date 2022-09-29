package org.calibrationframework.fouriermethod.products;

import java.util.Map;

import org.apache.commons.math3.complex.Complex;

import org.calibrationframework.exception.CalculationException;
import org.calibrationframework.fouriermethod.models.ProcessCharacteristicFunctionInterface;

/**
 * This is an abstract base class for all Fourier-based methodologies for the pricing of a smile of options.
 * 
 * Concrete different Fourier methodologies should provide different implementations of the getValue method, which is left
 * here as abstract.
 * 
 * @author Alessandro Gnoatto
 *
 */
public abstract class EuropeanOptionSmile implements SmileByIntegralTransform {
	
	private final double maturity;
	private final double[] strikes;
	
	public EuropeanOptionSmile(double maturity, double[] strikes) {
		super();
		this.maturity = maturity;
		this.strikes = strikes;
	}
	
	@Override
	public double getMaturity() {
		return this.maturity;
	}
	
	public double[] getStrikes() {
		return this.strikes;
	}

	@Override
	public double getIntegrationDomainImagLowerBound() {
		return -1;
	}

	@Override
	public double getIntegrationDomainImagUpperBound() {
		return 0;
	}

	public abstract Map<Double, Double> getValue(ProcessCharacteristicFunctionInterface model) throws CalculationException;
	
	/**
	 * This method allows us to reuse the same pricer (same pricing algorithm) over different option smiles.
	 * @param maturity
	 * @param strikes
	 * @return the same pricer now referring to a different smile.
	 */
	public abstract EuropeanOptionSmile getCloneWithModifiedParameters(double maturity, double[] strikes);
	

	@Override
	public Complex apply(Complex z) {
		return ((z.subtract(Complex.I)).multiply(z)).negate();
	}

}
