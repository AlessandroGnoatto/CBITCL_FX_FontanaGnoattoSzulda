package org.calibrationframework.marketdata.model.volatilities;

import java.time.LocalDate;

import org.calibrationframework.marketdata.model.volatilities.VolatilitySurfaceInterface.QuotingConvention;

/**
 * This class represents an equity option quote associated to its liquidity weight.
 * 
 * The weight comprises between zero and one.
 * 
 * Its normalization has been performed all over the corresponding surface, not only over one smile.
 * 
 * @author Szulda Guillaume
 *
 */
public class WeightedOptionData extends OptionData {
	
	private final double weight;
	
	public WeightedOptionData(String underlying, LocalDate referenceDate, double strike, double maturity, 
			double value, double weight, QuotingConvention convention) {
		super(underlying, referenceDate, strike, maturity, value, convention);
		this.weight = weight;
	}

	public double getWeight() {
		return this.weight;
	}

}
