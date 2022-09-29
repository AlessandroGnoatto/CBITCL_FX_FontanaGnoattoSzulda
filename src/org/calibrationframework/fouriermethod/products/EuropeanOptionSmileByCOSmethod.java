package org.calibrationframework.fouriermethod.products;

import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DiagonalMatrix;
import org.calibrationframework.exception.CalculationException;
import org.calibrationframework.fouriermethod.CharacteristicFunctionInterface;
import org.calibrationframework.fouriermethod.calibration.models.BlackScholesModel;
import org.calibrationframework.fouriermethod.models.ProcessCharacteristicFunctionInterface;

public class EuropeanOptionSmileByCOSmethod extends EuropeanOptionSmile {
	
	//Fields
	private final int L; 
	private final int numberOfPoints;
			
	//Constructors
	public EuropeanOptionSmileByCOSmethod(double maturity, double[] strikes) throws IllegalArgumentException {
			
		super(maturity, strikes);
			
		this.L = 10;
		this.numberOfPoints = 256;
			
	}
				
	public EuropeanOptionSmileByCOSmethod(double maturity, double[] strikes, int L, int numberOfPoints) throws IllegalArgumentException {
			
		super(maturity, strikes);
			
		this.L = L;
		this.numberOfPoints = numberOfPoints;
			
	}

	public Map<Double, Double> getValue(ProcessCharacteristicFunctionInterface model) throws CalculationException {

				CharacteristicFunctionInterface modelCF = model.apply(this.getMaturity());
				
				double a = 0;
				
				double b = 0;
				
				if (model instanceof BlackScholesModel) {
					
					double mu = Math.log( ((BlackScholesModel)model).getInitialValue() );
					
					double sigma = Math.sqrt( this.getMaturity()*((BlackScholesModel)model).getVolatility()*((BlackScholesModel)model).getVolatility() );
					
					a = mu - 3*sigma;
					
					b = mu + 3*sigma;
					
				}
				
				double upperBound = b-a;

				double[] diagOfSummand = new double[this.getStrikes().length];
				
				double u = 0;
				
				for(int k = 0; k < this.getStrikes().length; k++) {
					
					u = (1.0/upperBound)*( (Math.exp(b)/this.getStrikes()[k]) - 1.0 - b + Math.log(this.getStrikes()[k]) );
					
					diagOfSummand[k] = 0.5*u;
							
				}
				
				DiagonalMatrix summand = new DiagonalMatrix(diagOfSummand);
				
				double realPart = 0;

				for(int i = 1; i<numberOfPoints; i++) {
					
					realPart = ( ( modelCF.apply( Complex.I.multiply(i*Math.PI/upperBound) ) ).
							multiply( ( Complex.I.multiply( -i*a*Math.PI/upperBound ) ).exp() ) ).getReal();
					
					for(int k = 0; k < this.getStrikes().length; k++) {
						
						u = (2.0/upperBound)*( ( 1.0 /( 1.0+(i*Math.PI/upperBound)*(i*Math.PI/upperBound) ) )*( 
								Math.pow(-1,i)*(Math.exp(b)/this.getStrikes()[k]) - Math.cos( i*Math.PI*( a-Math.log(this.getStrikes()[k]) )/upperBound ) + ( i*Math.PI/upperBound )*Math.sin( i*Math.PI*( a-Math.log(this.getStrikes()[k]) )/upperBound ) ) - 
								( upperBound/(i*Math.PI) )*Math.sin( i*Math.PI*( a-Math.log(this.getStrikes()[k]) )/upperBound ) );
						
						diagOfSummand[k] = realPart*u;
								
					}
					
					summand = summand.add( new DiagonalMatrix(diagOfSummand) );
					
				}
				
				ArrayRealVector callPrices = new ArrayRealVector(summand.operate(this.getStrikes()));
				
				double[] prices = callPrices.toArray();

				HashMap<Double, Double> results = new HashMap<Double, Double>();

				for(int k = 0; k < this.getStrikes().length; k++) {
					results.put( this.getStrikes()[k], Math.abs( prices[k] ) );
				}

				return results;
				
			}
			
			@Override
			public EuropeanOptionSmile getCloneWithModifiedParameters(double maturity, double[] strikes) {
				
				return new EuropeanOptionSmileByCOSmethod(maturity, strikes, this.L, this.numberOfPoints);
				
			}

}
