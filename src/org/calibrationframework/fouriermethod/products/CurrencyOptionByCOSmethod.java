package org.calibrationframework.fouriermethod.products;

import java.util.*;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.DiagonalMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;

import org.calibrationframework.exception.CalculationException;
import org.calibrationframework.fouriermethod.CharacteristicFunctionInterface;
import org.calibrationframework.fouriermethod.models.*;

import net.finmath.marketdata.model.curves.*;

/**
 * This class implements the COS method (see Fang and Oosterlee (2008)).
 * 
 * It provides an approximation of the price of a currency option.
 * 
 * This can be seen as an alternative to CurrencyOptionByCarrMadan,
 * where no lineOfIntegration is required to proceed.
 *  
 * @author Szulda Guillaume
 *
 */
public class CurrencyOptionByCOSmethod extends EuropeanOptionSmileMultiAsset {
	
	//Fields
		private final int L;
		private final int numberOfPoints;
		
		private final double[] maturities;
		
		private final DiscountCurve discountCurveUSD;
		private final DiscountCurve forwardCurveEURUSD;
		private final DiscountCurve discountCurveJPY;
		private final DiscountCurve forwardCurveEURJPY;
		private final DiscountCurve forwardCurveUSDJPY;
			
		//Constructors
		public CurrencyOptionByCOSmethod(String underlyingName, double maturity, double[] strikes, 
				DiscountCurve discountCurveUSD, DiscountCurve forwardCurveEURUSD, 
				DiscountCurve discountCurveJPY, DiscountCurve forwardCurveEURJPY, DiscountCurve forwardCurveUSDJPY,
				double[] maturities) throws IllegalArgumentException {
			
			super(underlyingName, maturity, strikes);
			
			this.maturities = maturities;
			
			String underlyingNameDomestic = underlyingName.subSequence(3,6).toString();
			
			if(underlyingNameDomestic.contains("JPY")) {
				
				if(maturity == maturities[0]) {
					
					this.L = 250;
					this.numberOfPoints = 40;
					
				} else {
					
					this.L = 250;
					this.numberOfPoints = 10;
				}
				
			} else {
				
				this.L = 100;
				this.numberOfPoints = 100;
				
			}
			
			this.discountCurveUSD = discountCurveUSD;
			this.forwardCurveEURUSD = forwardCurveEURUSD;
			this.discountCurveJPY = discountCurveJPY;
			this.forwardCurveEURJPY = forwardCurveEURJPY;
			this.forwardCurveUSDJPY = forwardCurveUSDJPY;
			
		}
				
		public CurrencyOptionByCOSmethod(String underlyingName, double maturity, double[] strikes, int L, int numberOfPoints,
				DiscountCurve discountCurveUSD, DiscountCurve forwardCurveEURUSD, 
				DiscountCurve discountCurveJPY, DiscountCurve forwardCurveEURJPY, DiscountCurve forwardCurveUSDJPY) throws IllegalArgumentException {
			
			super(underlyingName, maturity, strikes);
			
			this.maturities = null;
			
			this.L = L;
			this.numberOfPoints = numberOfPoints;
			
			this.discountCurveUSD = discountCurveUSD;
			this.forwardCurveEURUSD = forwardCurveEURUSD;
			this.discountCurveJPY = discountCurveJPY;
			this.forwardCurveEURJPY = forwardCurveEURJPY;
			this.forwardCurveUSDJPY = forwardCurveUSDJPY;
			
		}

		public Map<Double, Double> getValue(MultivariateProcessCharacteristicFunctionInterface model) throws CalculationException {

				CharacteristicFunctionInterface modelCF = model.apply(this.getMaturity(), this.getUnderlyingName());
				
				String underlyingNameForeign = this.getUnderlyingName().subSequence(0,3).toString();
				String underlyingNameDomestic = this.getUnderlyingName().subSequence(3,6).toString();
				
				double discountFactor;
				double forward;
				
				if( underlyingNameDomestic.contains("USD") ) {
					discountFactor = this.discountCurveUSD.getDiscountFactor(this.getMaturity());
					forward = this.forwardCurveEURUSD.getDiscountFactor(this.getMaturity());
				} else {
					discountFactor = this.discountCurveJPY.getDiscountFactor(this.getMaturity());
					if( underlyingNameForeign.contains("USD") ) {
						forward = this.forwardCurveUSDJPY.getDiscountFactor(this.getMaturity());
					} else {
						forward = this.forwardCurveEURJPY.getDiscountFactor(this.getMaturity());
					}
				}
				
				/*
				 * Determination of the truncation range through the computation of the cumulants via forward finite differences,
				 * see Financial modelling Theory, Implementation and Practice with Matlab Source, book by Jörg Kienitz and Daniel Wetterau.
				 */
				double h = 1E-6;
				
				//First cumulant
				double c1 = (1.0/h)*( Math.log( modelCF.apply( new Complex(h) ).getReal() ) - Math.log(discountFactor) );
						
				//Second cumulant
				double c2 = (1.0/h*h)*( Math.log( modelCF.apply( new Complex(2*h) ).getReal() ) - 2*Math.log( modelCF.apply( new Complex(h) ).getReal() ) + Math.log(discountFactor) );
				
				//Fourth cumulant 
				double c4 = (1.0/h*h*h*h)*( Math.log( modelCF.apply( new Complex(4*h) ).getReal() ) - 
						4*Math.log( modelCF.apply( new Complex(3*h) ).getReal() ) + 6*Math.log( modelCF.apply( new Complex(2*h) ).getReal() ) - 
						4*Math.log( modelCF.apply( new Complex(h) ).getReal() ) + Math.log(discountFactor) );
				
				//Sixth cumulant 
				double c6 = (1.0/h*h*h*h*h*h)*( Math.log( modelCF.apply( new Complex(6*h) ).getReal() ) - 
						6*Math.log( modelCF.apply( new Complex(5*h) ).getReal() ) + 15*Math.log( modelCF.apply( new Complex(4*h) ).getReal() ) - 
						20*Math.log( modelCF.apply( new Complex(3*h) ).getReal() ) + 15*Math.log( modelCF.apply( new Complex(2*h) ).getReal() ) -
						6*Math.log( modelCF.apply( new Complex(h) ).getReal() ) + Math.log(discountFactor) );
				
				//Eighth cumulant
				double c8 = (1.0/h*h*h*h*h*h*h*h)*( Math.log( modelCF.apply( new Complex(8*h) ).getReal() ) - 
						8*Math.log( modelCF.apply( new Complex(7*h) ).getReal() ) + 28*Math.log( modelCF.apply( new Complex(6*h) ).getReal() ) - 
						56*Math.log( modelCF.apply( new Complex(5*h) ).getReal() ) + 70*Math.log( modelCF.apply( new Complex(4*h) ).getReal() ) -
						56*Math.log( modelCF.apply( new Complex(3*h) ).getReal() ) + 28*Math.log( modelCF.apply( new Complex(2*h) ).getReal() ) - 
						8*Math.log( modelCF.apply( new Complex(h) ).getReal() ) + Math.log(discountFactor) );
				
				//Tenth cumulant
				double c10 = (1.0/h*h*h*h*h*h*h*h*h*h)*( Math.log( modelCF.apply( new Complex(10*h) ).getReal() ) - 
						10*Math.log( modelCF.apply( new Complex(9*h) ).getReal() ) + 45*Math.log( modelCF.apply( new Complex(8*h) ).getReal() ) - 
						120*Math.log( modelCF.apply( new Complex(7*h) ).getReal() ) + 210*Math.log( modelCF.apply( new Complex(6*h) ).getReal() ) -
						252*Math.log( modelCF.apply( new Complex(5*h) ).getReal() ) + 210*Math.log( modelCF.apply( new Complex(4*h) ).getReal() ) - 
						120*Math.log( modelCF.apply( new Complex(3*h) ).getReal() ) + 45*Math.log( modelCF.apply( new Complex(2*h) ).getReal() ) - 
						10*Math.log( modelCF.apply( new Complex(h) ).getReal() ) + Math.log(discountFactor) );
				
				double a = c1 - this.L*Math.sqrt( Math.abs(c2) + Math.sqrt( Math.abs(c4) + Math.sqrt( 
						Math.abs(c6) + Math.sqrt( Math.abs(c8) + Math.sqrt( Math.abs(c10) ) ) ) ) );
				
				double b = c1 + this.L*Math.sqrt( Math.abs(c2) + Math.sqrt( Math.abs(c4) + Math.sqrt( 
						Math.abs(c6) + Math.sqrt( Math.abs(c8) + Math.sqrt( Math.abs(c10) ) ) ) ) );
				
				double upperBound = b-a;

				double[] diagOfSummand = new double[this.getStrikes().length];
				
				Complex cf = modelCF.apply( new Complex(0) );
				
				double u = 0;
				
				for(int k = 0; k < this.getStrikes().length; k++) {
					
					u = (1.0/upperBound)*( (Math.exp(a)/this.getStrikes()[k]) - 1.0 - a + Math.log(this.getStrikes()[k]) );
					
					diagOfSummand[k] = cf.getReal()*u;
							
				}
				
				DiagonalMatrix summand = new DiagonalMatrix(diagOfSummand);

				for(int i = 1; i<numberOfPoints; i++) {
					
					cf = modelCF.apply( Complex.I.multiply(i*Math.PI/upperBound) );
					
					for(int k = 0; k < this.getStrikes().length; k++) {
						
						u = (2.0/upperBound)*( ( 1.0 /( 1.0+(i*Math.PI/upperBound)*(i*Math.PI/upperBound) ) )*( 
								(Math.exp(a)/this.getStrikes()[k]) - Math.cos( i*Math.PI*( a-Math.log(this.getStrikes()[k]) )/upperBound ) + ( i*Math.PI/upperBound )*Math.sin( i*Math.PI*( a-Math.log(this.getStrikes()[k]) )/upperBound ) ) - 
								( upperBound/(i*Math.PI) )*Math.sin( i*Math.PI*( a-Math.log(this.getStrikes()[k]) )/upperBound ) );
						
						diagOfSummand[k] = ( cf.multiply( 
								( Complex.I.multiply( -i*a*Math.PI/upperBound ) ).exp() ) ).getReal()*u;
								
					}
					
					summand = summand.add( new DiagonalMatrix(diagOfSummand) );
					
				}
				
				ArrayRealVector putPrices = new ArrayRealVector(summand.operate(this.getStrikes()));
				
				//Call-Put parity (see Eq. (50), Fang and Oosterlee (2008)):
				ArrayRealVector StockMinusStrike = ( new ArrayRealVector( this.getStrikes().length, discountFactor*forward ) ).subtract(
						new ArrayRealVector( this.getStrikes() ).mapMultiply(discountFactor) );
				
				ArrayRealVector callPrices = putPrices.add(StockMinusStrike);
				
				double[] prices = callPrices.toArray();

				HashMap<Double, Double> results = new HashMap<Double, Double>();

				for(int k = 0; k < this.getStrikes().length; k++) {
					results.put( this.getStrikes()[k], Math.abs( prices[k] ) );
				}

				return results;
				
			}
			
			@Override
			public EuropeanOptionSmileMultiAsset getCloneWithModifiedParameters(String underlyingName, double maturity, double[] strikes) {
				
				if(this.maturities == null) {
					
					return new CurrencyOptionByCOSmethod(underlyingName, maturity, strikes, this.L, this.numberOfPoints, 
							this.discountCurveUSD, this.forwardCurveEURUSD, this.discountCurveJPY, this.forwardCurveEURJPY, this.forwardCurveUSDJPY);
					
				} else {
					
					return new CurrencyOptionByCOSmethod(underlyingName, maturity, strikes, 
							this.discountCurveUSD, this.forwardCurveEURUSD, this.discountCurveJPY, this.forwardCurveEURJPY, this.forwardCurveUSDJPY, 
							this.maturities);
					
				}
				
			}
	/*
	//Fields
	private final int L; 
	private final int numberOfPoints;
		
	//Constructors
	public CurrencyOptionByCOSmethod(String underlyingName, double maturity, double[] strikes) throws IllegalArgumentException {
		
		super(underlyingName, maturity, strikes);
		
		this.L = 130;
		this.numberOfPoints = 1024;
		
	}
			
	public CurrencyOptionByCOSmethod(String underlyingName, double maturity, double[] strikes, int L, int numberOfPoints) throws IllegalArgumentException {
		
		super(underlyingName, maturity, strikes);
		
		this.L = L;
		this.numberOfPoints = numberOfPoints;
		
	}

	public Map<Double, Double> getValue(MultivariateProcessCharacteristicFunctionInterface model) throws CalculationException {

			CharacteristicFunctionInterface modelCF = model.apply(this.getMaturity(), this.getUnderlyingName());
			
			/* Determination of the truncation range through the computation of the cumulants via forward finite differences,
			 * see Financial modelling Theory, Implementation and Practice with Matlab Source, book by Jörg Kienitz and Daniel Wetterau.
			 *//*
			double h = 1E-9;
			
			//First cumulant
			double c1 = (1.0/h)*Math.log( modelCF.apply( new Complex(h) ).getReal() );
					
			//Second cumulant
			double c2 = (1.0/h*h)*( Math.log( modelCF.apply( new Complex(2*h) ).getReal() ) - 2*Math.log( modelCF.apply( new Complex(h) ).getReal() ) );
			
			//Fourth cumulant 
			double c4 = (1.0/h*h*h*h)*( Math.log( modelCF.apply( new Complex(4*h) ).getReal() ) - 
					4*Math.log( modelCF.apply( new Complex(3*h) ).getReal() ) + 6*Math.log( modelCF.apply( new Complex(2*h) ).getReal() ) - 
					4*Math.log( modelCF.apply( new Complex(h) ).getReal() ) );
			
			//Sixth cumulant 
			double c6 = (1.0/h*h*h*h*h*h)*( Math.log( modelCF.apply( new Complex(6*h) ).getReal() ) - 
					6*Math.log( modelCF.apply( new Complex(5*h) ).getReal() ) + 15*Math.log( modelCF.apply( new Complex(4*h) ).getReal() ) - 
					20*Math.log( modelCF.apply( new Complex(3*h) ).getReal() ) + 15*Math.log( modelCF.apply( new Complex(2*h) ).getReal() ) -
					6*Math.log( modelCF.apply( new Complex(h) ).getReal() ) );
			
			//Eighth cumulant
			double c8 = (1.0/h*h*h*h*h*h*h*h)*( Math.log( modelCF.apply( new Complex(8*h) ).getReal() ) - 
					8*Math.log( modelCF.apply( new Complex(7*h) ).getReal() ) + 28*Math.log( modelCF.apply( new Complex(6*h) ).getReal() ) - 
					56*Math.log( modelCF.apply( new Complex(5*h) ).getReal() ) + 70*Math.log( modelCF.apply( new Complex(4*h) ).getReal() ) -
					56*Math.log( modelCF.apply( new Complex(3*h) ).getReal() ) + 28*Math.log( modelCF.apply( new Complex(2*h) ).getReal() ) - 
					8*Math.log( modelCF.apply( new Complex(h) ).getReal() ) );
			
			//Tenth cumulant
			double c10 = (1.0/h*h*h*h*h*h*h*h*h*h)*( Math.log( modelCF.apply( new Complex(10*h) ).getReal() ) - 
					10*Math.log( modelCF.apply( new Complex(9*h) ).getReal() ) + 45*Math.log( modelCF.apply( new Complex(8*h) ).getReal() ) - 
					120*Math.log( modelCF.apply( new Complex(7*h) ).getReal() ) + 210*Math.log( modelCF.apply( new Complex(6*h) ).getReal() ) -
					252*Math.log( modelCF.apply( new Complex(5*h) ).getReal() ) + 210*Math.log( modelCF.apply( new Complex(4*h) ).getReal() ) - 
					120*Math.log( modelCF.apply( new Complex(3*h) ).getReal() ) + 45*Math.log( modelCF.apply( new Complex(2*h) ).getReal() ) - 
					10*Math.log( modelCF.apply( new Complex(h) ).getReal() ) );
			
			double a = c1 - this.L*Math.sqrt( Math.abs(c2) + Math.sqrt( Math.abs(c4) + Math.sqrt( 
					Math.abs(c6) + Math.sqrt( Math.abs(c8) + Math.sqrt( Math.abs(c10) ) ) ) ) );
			
			double b = c1 + this.L*Math.sqrt( Math.abs(c2) + Math.sqrt( Math.abs(c4) + Math.sqrt( 
					Math.abs(c6) + Math.sqrt( Math.abs(c8) + Math.sqrt( Math.abs(c10) ) ) ) ) );
			
			double upperBound = b-a;
			
			System.out.println("Upper bound = " + upperBound);

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
		public EuropeanOptionSmileMultiAsset getCloneWithModifiedParameters(String underlyingName, double maturity, double[] strikes) {
			
			return new CurrencyOptionByCOSmethod(underlyingName, maturity, strikes, this.L, this.numberOfPoints);
			
		}*/

}
