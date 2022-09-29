package org.calibrationframework.marketdata.model.volatilities;

import java.time.LocalDate;

import org.calibrationframework.marketdata.model.volatilities.VolatilitySurfaceInterface.QuotingConvention;

/**
 * Representation of a collection of option prices or implied volatilities for a given maturity.
 * 
 * @author Szulda Guillaume
 *
 */
public interface SmileInterface {
	
	public String getUnderlying();
	
	public LocalDate getReferenceDate();

	public double[] getStrikes();

	public double getMaturity();
	
	public QuotingConvention getQuotingConvention();
	
	public int getSize();

}
