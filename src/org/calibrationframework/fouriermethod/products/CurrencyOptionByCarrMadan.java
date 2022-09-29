package org.calibrationframework.fouriermethod.products;

import java.util.*;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.*;

import org.calibrationframework.exception.CalculationException;
import org.calibrationframework.fouriermethod.CharacteristicFunctionInterface;
import org.calibrationframework.fouriermethod.models.MultivariateProcessCharacteristicFunctionInterface;
import net.finmath.interpolation.RationalFunctionInterpolation;
import net.finmath.interpolation.RationalFunctionInterpolation.ExtrapolationMethod;
import net.finmath.interpolation.RationalFunctionInterpolation.InterpolationMethod;

public class CurrencyOptionByCarrMadan extends EuropeanOptionSmileMultiAsset {
	
	//Fields
	private final double lineOfIntegration;
	private final int numberOfPoints;
	private final double gridSpacing;
	private final InterpolationMethod intMethod;
	private final ExtrapolationMethod extMethod;
	
	//Constructors
	public CurrencyOptionByCarrMadan(String underlyingName, double maturity, double[] strikes, double lineOfIntegration) throws IllegalArgumentException {
		super(underlyingName, maturity, strikes);
		if(lineOfIntegration > 1) {
			this.lineOfIntegration = lineOfIntegration;
			this.numberOfPoints = 1024;
			this.gridSpacing = 0.05;
			this.intMethod =InterpolationMethod.LINEAR /*HARMONIC_SPLINE*/;
			this.extMethod = ExtrapolationMethod.CONSTANT;
		} else {
			throw new IllegalArgumentException("The line of integration must be strictly greater than one");
		}
	}
		
	public CurrencyOptionByCarrMadan(String underlyingName, double maturity, double[] strikes, double lineOfIntegration, int numberOfPoints,
			double gridSpacing, InterpolationMethod intMethod, ExtrapolationMethod extMethod) throws IllegalArgumentException {
		super(underlyingName, maturity, strikes);
		if(lineOfIntegration > 1) {
			this.lineOfIntegration = lineOfIntegration;
			this.numberOfPoints = numberOfPoints;
			this.gridSpacing = gridSpacing;
			this.intMethod = intMethod;
			this.extMethod = extMethod;
		} else {
			throw new IllegalArgumentException("The line of integration must be strictly greater than one");
		}	
	}

	public Map<Double, Double> getValue(MultivariateProcessCharacteristicFunctionInterface model) throws CalculationException {

		CharacteristicFunctionInterface modelCF = model.apply(getMaturity(), this.getUnderlyingName());
		
		String underlyingNameDomestic = this.getUnderlyingName().subSequence(3,6).toString();

		double lambda = 2*Math.PI/(numberOfPoints*gridSpacing); //Equation 23 Carr and Madan
		double upperBound = (numberOfPoints * lambda)/2.0; //Equation 20 Carr and Madan

		Complex[] integrandEvaluations = new Complex[numberOfPoints];

		for(int i = 0; i<numberOfPoints; i++) {

			double u = gridSpacing * i;

			Complex z = new Complex(lineOfIntegration, u);
			
			Complex numerator = modelCF.apply(z);
			
			Complex denominator = (z.subtract(1)).multiply(z);
			Complex ratio = numerator.divide(denominator);
			
			if( underlyingNameDomestic.contains("JPY") ) {
				ratio = (ratio.multiply(((Complex.I).multiply( u*( upperBound - Math.log(110.0) ) )).exp())).multiply(gridSpacing);
			} else {
				ratio = (ratio.multiply(((Complex.I).multiply(upperBound*u)).exp())).multiply(gridSpacing);
			}

			double delta;
			if (i==0){
				delta=1.0;
			}else{
				delta = 0.0;
			}	
			double simpsonWeight = (3+Math.pow(-1,i+1)-delta)/3;

			integrandEvaluations[i] = ratio.multiply(simpsonWeight);
			
		}

		//Compute the FFT
		Complex[] transformedVector = new Complex[numberOfPoints];
		FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);
		transformedVector = fft.transform(integrandEvaluations, TransformType.FORWARD);

		//Find relevant prices via interpolation
		double[] logStrikeVector = new double[numberOfPoints];
		double[] strikeVector = new double[numberOfPoints];
		double[] optionPriceVector = new double[numberOfPoints];

		for(int j = 0; j<numberOfPoints; j++) {
			
			if( underlyingNameDomestic.contains("JPY") ) {
				logStrikeVector[j] = Math.log(110.0)-upperBound+lambda*j;
			} else {
				logStrikeVector[j] = -upperBound+lambda*j;
			}
			
			strikeVector[j] = Math.exp(logStrikeVector[j]);
			optionPriceVector[j] = ( transformedVector[j].multiply( Math.exp( (1-lineOfIntegration)*logStrikeVector[j] ) ) ).getReal()/Math.PI;
			
		}

		RationalFunctionInterpolation interpolation = new RationalFunctionInterpolation(strikeVector, optionPriceVector, intMethod, extMethod);

		double[] strikes = getStrikes();

		int numberOfStrikes = strikes.length;
		HashMap<Double, Double> results = new HashMap<Double, Double>();

		for(int k = 0; k < numberOfStrikes; k++) {
			results.put( strikes[k], Math.abs( interpolation.getValue(strikes[k]) ) );
		}

		return results;
		
	}
	
	@Override
	public EuropeanOptionSmileMultiAsset getCloneWithModifiedParameters(String underlyingName, double maturity, double[] strikes) {
		
		return new CurrencyOptionByCarrMadan(underlyingName, maturity, strikes, this.lineOfIntegration, this.numberOfPoints,
				this.gridSpacing, this.intMethod, this.extMethod);
		
	}

}
