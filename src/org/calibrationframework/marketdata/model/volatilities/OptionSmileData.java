package org.calibrationframework.marketdata.model.volatilities;

import java.time.LocalDate;
import java.util.Arrays;
import java.util.HashMap;

import org.calibrationframework.marketdata.model.volatilities.OptionData;
import org.calibrationframework.marketdata.model.volatilities.VolatilitySurfaceInterface.QuotingConvention;

/**
 * A collection of option prices or implied volatilities for a given maturity.
 * 
 * This object is natural in conjunction with FFT methods that run over a whole smile for a fixed maturity.
 * 
 * @author Alessandro Gnoatto
 *
 */
public class OptionSmileData implements SmileInterface {
	
	private final String underlying;	
	private final LocalDate referenceDate;
	private final double[] strikes;
	private final double maturity;
	private final QuotingConvention convention;
	private final HashMap<Double,OptionData> smile;
	
	public OptionSmileData(String underlying, LocalDate referenceDate, double[] strikes, double maturity, 
			double[] values, QuotingConvention convention) throws IllegalArgumentException {
		
		if(strikes.length != values.length) {
			throw new IllegalArgumentException("Number of strikes and market quotes does not coincide");
		}else {
			int numberOfQuotes = strikes.length;
			smile = new HashMap<Double,OptionData>();
			for(int i = 0; i< numberOfQuotes; i++) {
				smile.put(strikes[i],new OptionData(underlying, referenceDate, strikes[i],maturity,values[i], convention)) ;
			}
			this.underlying = underlying;
			this.referenceDate = referenceDate;
			this.strikes = strikes;
			this.maturity = maturity;
			this.convention = convention;
		}
		
	}
	
	public HashMap<Double, OptionData> getSmile() {
		return smile;
	}
	
	public OptionData getOption(double strike) {
		return smile.get(strike);
	}

	@Override
	public String getUnderlying() {
		return underlying;
	}
	
	@Override
	public LocalDate getReferenceDate() {
		return referenceDate;
	}

	@Override
	public double[] getStrikes() {
		return strikes;
	}

	@Override
	public double getMaturity() {
		return maturity;
	}
	
	@Override
	public QuotingConvention getQuotingConvention() {
		return convention;
	}
	
	@Override
	public int getSize() {
		return smile.size();
	}
	
	@Override
	public String toString() {
		return "EquityOptionSmile [underlying=" + underlying + ", strikes=" + Arrays.toString(strikes)
				+ ", maturity=" + maturity + ", smile=" + smile + "]";
	}
	
}
