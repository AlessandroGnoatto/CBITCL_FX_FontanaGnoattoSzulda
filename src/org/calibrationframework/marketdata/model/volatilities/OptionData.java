package org.calibrationframework.marketdata.model.volatilities;

import java.time.LocalDate;

import org.calibrationframework.marketdata.model.volatilities.VolatilitySurfaceInterface.QuotingConvention;

/**
 * An Equity option quote is a function of strike and maturity. The quote can be represented in terms of prices or volatilities.
 * Concerning the strike: being a double, one might decide to store there a moneyness instead of a price, i.e. a relative strike where ATM = 0.
 * 
 * @author Alessandro Gnoatto
 *
 */
public class OptionData implements OptionInterface {
	
	private final String underlying;
	private final LocalDate referenceDate;
	private final double strike;
	private final double maturity;
	private final double value;
	private final QuotingConvention convention;
	
	public OptionData(String underlying, LocalDate referenceDate, double strike, double maturity, double value, QuotingConvention convention) {
		super();
		this.underlying = underlying;
		this.referenceDate = referenceDate;
		this.strike = strike;
		this.maturity = maturity;
		this.value = value;
		this.convention = convention;
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
	public double getStrike() {
		return strike;
	}

	@Override
	public double getMaturity() {
		return maturity;
	}

	@Override
	public double getValue() {
		return value;
	}

	@Override
	public QuotingConvention getConvention() {
		return convention;
	}

	@Override
	public String toString() {
		return "EquityOptionQuote [underlying=" + underlying + ", referenceDate=" + referenceDate + ", strike="
				+ strike + ", maturity=" + maturity + ", value=" + value + ", convention=" + convention + "]";
	}
	
}