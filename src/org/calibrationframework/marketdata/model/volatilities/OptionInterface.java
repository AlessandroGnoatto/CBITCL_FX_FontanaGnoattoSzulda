package org.calibrationframework.marketdata.model.volatilities;

import java.time.LocalDate;

import org.calibrationframework.marketdata.model.volatilities.VolatilitySurfaceInterface.QuotingConvention;

/**
 * Representation of an equity quote option, depending on the maturity and strike, 
 * 
 * the date on which the product is issued, the quoting convention along with,
 * 
 * the name of the underlying over which the option is written.
 * 
 * @author Szulda Guillaume
 *
 */
public interface OptionInterface {
	
	public String getUnderlying();

	public LocalDate getReferenceDate();

	public double getStrike();

	public double getMaturity();

	public double getValue();

	public QuotingConvention getConvention();

}
