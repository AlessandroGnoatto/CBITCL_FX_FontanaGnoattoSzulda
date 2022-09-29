package org.calibrationframework.fouriermethod.quantization;

import java.util.HashMap;
import java.util.Map;

import net.finmath.stochastic.RandomVariable;

/**
 * Multi curve caplet pricer based on quantization:
 * Given a QuantizableMultiDimProcessInterface model and one of its components denoted by String underlyingName,
 * this pricer provides a whole smile of caplet prices for a given maturity and a given array of strikes, 
 * via the method "Map<Double, Double> getArrayOfValue(QuantizableMultiDimProcessInterface model)",
 * through the quantization of the "underlyingName" component of QuantizableMultiDimProcessInterface model at maturity.
 * 
 * @author Szulda Guillaume
 */

public class QuantizationMultiCurveSmileOfCaplets extends PricingByQuadraticQuantization {
	
	private final String underlyingName;  //standing for the underlying tenor/curve of the caplet.
	private final double maturity;
	private final double[] strikes; 
	
	public QuantizationMultiCurveSmileOfCaplets(String underlyingName, double maturity, double[] strikes, String currency) {
		super(currency);
		this.strikes = strikes;
		this.maturity = maturity;
		this.underlyingName = underlyingName;
	}
	
	public QuantizationMultiCurveSmileOfCaplets(String underlyingName, double maturity, double[] strikes) {
		super();
		this.strikes = strikes;
		this.maturity = maturity;
		this.underlyingName = underlyingName;
	}
	
	@Override
	public Double getValue(Quantizable1DProcessInterface model) throws IllegalArgumentException {
		return null;
	}
	
	@Override
	public Double getValue(QuantizableMultiDimProcessInterface model) throws IllegalArgumentException {
		if(model instanceof QuantizableCBIDrivenMultiCurveModel) {
			try {
				return getValue((QuantizableCBIDrivenMultiCurveModel)model);
			} catch (IllegalArgumentException e) {
				return null;
			}
		}
		else {
			throw new IllegalArgumentException("The product " + this.getClass()
			+ " cannot be valued against a model " + model.getClass() + "."
			+ "It requires a model of type " + QuantizableCBIDrivenMultiCurveModel.class + ".");
		}
	}
	
	public Double getValue(QuantizableCBIDrivenMultiCurveModel model) throws IllegalArgumentException {
		if(this.maturity != model.getMaturity() || underlyingToTenor(this.underlyingName) != model.getTenorLength()) {
			throw new IllegalArgumentException("The maturity along with the tenor of the payoff must be equal to the evaluation time at which the underlying process has been quantized and its associated tenor (respectively).");
		} else {
			/* Declaration of the payoff of the product to price by means of the Voronoi quantization (QuantizableCBIDrivenMultiCurveModel) of its underlying process (CBIDrivenMultiCurveModel).*/ 
			RandomVariable payoff = ((model.getVoronoiQuantization()).sub(1+this.strikes[0]*model.getTenorLength())).floor(0);
			/* Use of the cubature formula to approximate the time-0 price of the caplet by means of the expectation of the payoff function of the quantization,
			 * whose distribution is known as a discrete law over the corresponding quantization grid with the associated companion weights.
			 * */
			return payoff.getAverage(model.getCompanionWeights())*model.getLevel()*(model.getDiscountCurve()).getDiscountFactor(this.maturity+model.getTenorLength());
		}
	}
	
	public QuantizationMultiCurveSmileOfCaplets getCloneWithModifiedParameters(String underlyingName, double maturity, double[] strikes) {
		return new QuantizationMultiCurveSmileOfCaplets(underlyingName, maturity, strikes);
	}
	
	public Map<Double, Double> getArrayOfValue(QuantizableMultiDimProcessInterface model) throws IllegalArgumentException {
		if(model instanceof QuantizableCBIDrivenMultiCurveModel) {
			try {
				return getArrayOfValue((QuantizableCBIDrivenMultiCurveModel)model);
			} catch (IllegalArgumentException e) {
				return null;
			} 
		}
		else {
			throw new IllegalArgumentException("The product " + this.getClass()
			+ " cannot be valued against a model " + model.getClass() + "."
			+ "It requires a model of type " + QuantizableCBIDrivenMultiCurveModel.class + ".");
		}
	}
	
	public Map<Double, Double> getArrayOfValue(QuantizableCBIDrivenMultiCurveModel model) throws IllegalArgumentException {
		if(this.maturity != model.getMaturity() || underlyingToTenor(this.underlyingName) != model.getTenorLength()) {
			throw new IllegalArgumentException("The maturity along with the tenor of the payoff must be equal to the evaluation time at which the underlying process has been quantized and its associated tenor (respectively).");
		} else {
			HashMap<Double, Double> results = new HashMap<Double, Double>();
			int numberOfStrikes = (this.strikes).length;
			for(int k = 0; k<numberOfStrikes; k++) {
				RandomVariable kthPayoff = ((model.getVoronoiQuantization()).sub(1+this.strikes[k]*model.getTenorLength())).floor(0);
				double kthPrice = kthPayoff.getAverage(model.getCompanionWeights())*model.getLevel()*(model.getDiscountCurve()).getDiscountFactor(this.maturity+model.getTenorLength());
				results.put(this.strikes[k], kthPrice);
			}
			return results;
		}
	}
	
	private double underlyingToTenor(String underlyingName) {
		if(underlyingName.contains("3M")) {
			return 0.25;
		}else if(underlyingName.contains("6M")) {
			return 0.5;
		}else if(underlyingName.contains("12M")) {
			return 1.0;
			
		}else if(underlyingName.contains("1Y")) {
			return 1.0;
			
		}else {
			throw new IllegalArgumentException("Tenor not recognized.");
		}
	}

}
