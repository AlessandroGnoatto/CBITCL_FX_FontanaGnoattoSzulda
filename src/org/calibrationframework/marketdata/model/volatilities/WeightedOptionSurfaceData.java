package org.calibrationframework.marketdata.model.volatilities;

import java.time.LocalDate;
import java.util.HashMap;

import net.finmath.marketdata.model.AnalyticModel;
import net.finmath.marketdata.model.curves.DiscountCurve;

/**
 * An option quote surface taking into account, for each option quote of the surface, the associated weight.
 * 
 * Most often, this weight is computed as the inverse bid-ask spread, measuring the liquidity of the option.
 * 
 * Make sure that the weights are all normalized before inputting them, they are then comprised between zero and one.
 * 
 * Normalization has to be performed all over the option quote surface.
 * 
 * @author Szulda Guillaume.
 *
 */
public class WeightedOptionSurfaceData extends OptionSurfaceData {
	
	private final HashMap<Double, WeightedOptionSmileData> weightedSurface;
	
	/**
	 * This is a very restrictive constructor that assumes that for each maturity we have the same number of option quotes.
	 * @param underlying
	 * @param referenceDate
	 * @param strikes
	 * @param maturities
	 * @param values
	 * @param weights
	 * @param convention
	 * @param discountCurve
	 * @param equityForwardCurve
	 * @throws IllegalArgumentException
	 */
	public WeightedOptionSurfaceData(String underlying, LocalDate referenceDate, double[] strikes,
			double[] maturities, double[][] values, double[][] weights,
			QuotingConvention convention,DiscountCurve discountCurve, DiscountCurve equityForwardCurve) 
	throws IllegalArgumentException {
		
		super(underlying, referenceDate, strikes, maturities, values, convention, discountCurve, equityForwardCurve); 
		
		if(weights.length != values.length || weights[0].length != values[0].length ) {
			throw new IllegalArgumentException("Inconsistent number of weights and values");
		}else {
			
			this.weightedSurface = new HashMap<Double, WeightedOptionSmileData>();
			
			for(int j = 0; j< maturities.length; j++) {
				
				double[] valuesOfInterest = new double[strikes.length];
				double[] weightsOfInterest = new double[strikes.length];
				
				for(int i= 0; i< strikes.length; i++) {
					valuesOfInterest[i] = values[i][j];
					weightsOfInterest[i] = weights[i][j];
				}
				
				WeightedOptionSmileData jthSmile = new WeightedOptionSmileData(underlying, referenceDate, strikes, maturities[j], valuesOfInterest, weightsOfInterest, convention);
				this.weightedSurface.put(maturities[j],jthSmile);
				
			}
								
		}
				
	}
	
	
	/**
	 * Creates a weighted equity option surface from an array of smiles consisting of weighted option quotes.
	 * @param smiles
	 * @param discountCurve
	 * @param equityForwardCurve
	 */
	public WeightedOptionSurfaceData(WeightedOptionSmileData[] smiles, DiscountCurve discountCurve,
			DiscountCurve equityForwardCurve) throws IllegalArgumentException {
		
		super(smiles, discountCurve, equityForwardCurve);
		
		this.weightedSurface = new HashMap<Double, WeightedOptionSmileData>();
		
		for(int t = 0; t<smiles.length;t++) {
			
			double maturity = smiles[t].getMaturity();
			
			this.weightedSurface.put(maturity, smiles[t]);
		}
		
	}
	
	public HashMap<Double, WeightedOptionSmileData> getWeightedSurface(){
		return this.weightedSurface;
	}
	
	public double getWeight(double maturity, double strike){
		return getWeight(null, maturity, strike);
	}
	
	public double getWeight(AnalyticModel model, double maturity, double strike) {
		
			WeightedOptionSmileData relevantSmile = this.weightedSurface.get(maturity);
			
			return relevantSmile.getWeightedSmile().get(strike).getWeight();
	
	}
	
	public WeightedOptionSmileData getWeightedSmile(double maturity) {
		return this.weightedSurface.get(maturity);
	}
	
	public double getSmileWeight(double maturity) {
		return this.weightedSurface.get(maturity).getSmileWeight();
	}

}
