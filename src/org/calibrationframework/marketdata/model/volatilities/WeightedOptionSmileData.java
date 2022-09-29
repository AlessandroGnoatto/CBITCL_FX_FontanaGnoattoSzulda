package org.calibrationframework.marketdata.model.volatilities;

import java.time.LocalDate;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import org.calibrationframework.marketdata.model.volatilities.VolatilitySurfaceInterface.QuotingConvention;

/**
 * This class implements, for a given marturity, a collection of couples of the form:
 * option prices or implied volatilities coupled with its liquidity weight.
 * 
 * The normalization of the weights has to be performed all over the surface, not over one smile.
 * 
 * @author Szulda Guillaume
 *
 */
public class WeightedOptionSmileData extends OptionSmileData {
	
	private final HashMap<Double, WeightedOptionData> weightedSmile;
	
	public WeightedOptionSmileData(String underlying, LocalDate referenceDate, double[] strikes, double maturity, 
			double[] values, double[] weights, QuotingConvention convention) throws IllegalArgumentException {
		
		super(underlying, referenceDate, strikes, maturity, values, convention);
		
		if(weights.length != values.length) {
			throw new IllegalArgumentException("Weights and market quotes do not coincide");
		}else {
			int numberOfQuotes = strikes.length;
			this.weightedSmile = new HashMap<Double, WeightedOptionData>();
			for(int i = 0; i< numberOfQuotes; i++) {
				this.weightedSmile.put(strikes[i], new WeightedOptionData(underlying, referenceDate, strikes[i], maturity, values[i], weights[i], convention)) ;
			}
		}
		
	}

	public HashMap<Double, WeightedOptionData> getWeightedSmile() {
		return this.weightedSmile;
	}
	
	public WeightedOptionData getWeightedOption(double strike) {
		return this.weightedSmile.get(strike);
	}
	
	public double getSmileWeight() {
		
		double weight = 0;
	
		Iterator<Map.Entry<Double, WeightedOptionData>> iterator = this.weightedSmile.entrySet().iterator();
	
		while(iterator.hasNext()) {
					
			Map.Entry<Double, WeightedOptionData> pair = 
					(Map.Entry<Double, WeightedOptionData>) iterator.next();
			
			WeightedOptionData option = pair.getValue();
			
			weight = weight + option.getWeight();
			
		}
		
		return weight;
		
	}

}
