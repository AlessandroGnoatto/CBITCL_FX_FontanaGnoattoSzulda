package org.calibrationframework.marketdata.model.volatilities;

import java.time.LocalDate;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import net.finmath.marketdata.model.AnalyticModel;
import net.finmath.marketdata.model.curves.DiscountCurve;


/**
 * An option quote surface with the ability to query option quotes for different strikes and maturities.
 * 
 * The surface is constructed as a collection of smiles. The choice of this dimension is convenient in view of calibration via FFT methods.
 * 
 * This class does not perform any interpolation of market quotes. It merely represents a container of information.
 * 
 * The class provides also the ability to perform the conversion among different quoting conventions and hence can be used both for a calibration on prices or implied volatilities.
 * 
 * The class currently does not cover normal volatilities. Lognormal volatilities are more common in the equity space. The extension is not problematic.
 * 
 * @author Alessandro Gnoatto
 *
 */
public class OptionSurfaceData implements VolatilitySurfaceInterface {
	
	private final String underlying; 
	private final LocalDate referenceDate;
	private final DiscountCurve discountCurve; //\exp(-r*T) needed e.g. for application of the B&S formula
	private final DiscountCurve equityForwardCurve; //S0*\exp((r-d)*T)
	private final QuotingConvention convention; //either price or volatility (lognormal/normal)
	private final HashMap<Double, OptionSmileData> surface;
	private final double[] maturities;
		
	/**
	 * This is a very restrictive constructor that assumes that for each maturity we have the same number of option quotes.
	 * @param underlying
	 * @param referenceDate
	 * @param strikes
	 * @param maturities
	 * @param values
	 * @param convention
	 * @param discountCurve
	 * @param equityForwardCurve
	 * @throws IllegalArgumentException
	 */
	public OptionSurfaceData(String underlying, LocalDate referenceDate, double[] strikes,
			double[] maturities, double[][] values, QuotingConvention convention,DiscountCurve discountCurve,
			DiscountCurve equityForwardCurve) throws IllegalArgumentException {
		
		if(strikes.length != values.length || maturities.length != values[0].length ) {
			throw new IllegalArgumentException("Inconsistent number of strikes, maturities or values");
		}else {
			
			this.surface = new HashMap<Double,OptionSmileData>();
			
			for(int j = 0; j< maturities.length; j++) {
				
				double[] valuesOfInterest = new double[strikes.length];
				
				for(int i= 0; i< strikes.length; i++) {
					valuesOfInterest[i] = values[i][j];
				}
				
				OptionSmileData jthSmile = new OptionSmileData(underlying, referenceDate, strikes, maturities[j], valuesOfInterest, convention);
				this.surface.put(maturities[j],jthSmile);
				
			}
			
			this.underlying = underlying;
			this.referenceDate = referenceDate;
			this.discountCurve = discountCurve;
			this.equityForwardCurve = equityForwardCurve;
			this.convention = convention;
			this.maturities = maturities;
								
		}
				
	}
	
	
	/**
	 * Creates an equity option surface from an array of smiles.
	 * @param smiles
	 * @param discountCurve
	 * @param equityForwardCurve
	 */
	public OptionSurfaceData(OptionSmileData[] smiles, DiscountCurve discountCurve,
			DiscountCurve equityForwardCurve) throws IllegalArgumentException  {
		
		OptionSmileData firstSmile = smiles[0];
		String myUnderlying = firstSmile.getUnderlying();
		LocalDate myReferenceDate = firstSmile.getReferenceDate();
		QuotingConvention myConvention = firstSmile.getSmile().get(firstSmile.getStrikes()[0]).getConvention();
		
		HashMap<Double, OptionSmileData> mySurface = new HashMap<Double, OptionSmileData>();
		double[] mats = new double[smiles.length];
		
		for(int t = 0; t<smiles.length;t++) {
			double maturity = smiles[t].getMaturity();
			mats[t] = maturity;
			
			if(!(smiles[t].getReferenceDate().equals(myReferenceDate)))
				throw new IllegalArgumentException("All reference dates must be equal");
			
			if(!(smiles[t].getUnderlying().equals(myUnderlying)))
				throw new IllegalArgumentException("Option must be written on the same underlying");
			
			QuotingConvention testConvention = smiles[t].getSmile().get(smiles[t].getStrikes()[0]).getConvention();
			
			if(!(testConvention.equals(myConvention)))
				throw new IllegalArgumentException("Convention must be the same for all points in the surface");
			
			mySurface.put(maturity, smiles[t]);
		}
		this.underlying = myUnderlying;
		this.referenceDate = myReferenceDate;
		this.discountCurve = discountCurve;
		this.equityForwardCurve = equityForwardCurve;
		this.surface = mySurface;
		this.convention = myConvention;
		this.maturities = mats;
		
	}
	
	public HashMap<Double, OptionSmileData> getSurface() {
		return this.surface;
	}

	public OptionSmileData getSmile(double maturity) {
		return this.surface.get(maturity);
	}
	
	public DiscountCurve getDiscountCurve() {
		return this.discountCurve;
	}
	
	public DiscountCurve getEquityForwardCurve() {
		return this.equityForwardCurve;
	}
	
	@Override
	public String getName() {
		return this.underlying;
	}

	@Override
	public LocalDate getReferenceDate() {
		return this.referenceDate;
	}
	
	@Override
	public QuotingConvention getQuotingConvention() {
		return this.convention;
	}

	@Override
	public double[] getMaturities() {
		return this.maturities;
	}
	
	@Override
	public double getValue(double maturity, double strike){
		return getValue(maturity, strike, this.convention);
	}
	
	@Override
	public double getValue(double maturity, double strike, QuotingConvention quotingConvention)  {

		return getValue(null,maturity,strike,quotingConvention);
	}
	
	@Override
	public double getValue(AnalyticModel model, double maturity, double strike,
			QuotingConvention quotingConvention) {
		if(quotingConvention.equals(this.convention)) {
			OptionSmileData relevantSmile = this.surface.get(maturity);
			
			return relevantSmile.getSmile().get(strike).getValue();
		}else {
			if(quotingConvention == QuotingConvention.PRICE && this.convention == QuotingConvention.VOLATILITYLOGNORMAL) {
				
				double forwardPrice = equityForwardCurve.getValue(maturity);
				double discountBond = discountCurve.getValue(maturity);
				OptionSmileData relevantSmile = this.surface.get(maturity);
				double volatility = relevantSmile.getSmile().get(strike).getValue();
				return net.finmath.functions.AnalyticFormulas.blackScholesGeneralizedOptionValue(forwardPrice, volatility, maturity, strike, discountBond);
								
			}else if(quotingConvention == QuotingConvention.VOLATILITYLOGNORMAL && this.convention == QuotingConvention.PRICE) {
				double forwardPrice = equityForwardCurve.getValue(maturity);
				double discountBond = discountCurve.getValue(maturity);
				OptionSmileData relevantSmile = this.surface.get(maturity);
				double price = relevantSmile.getSmile().get(strike).getValue();				
				return net.finmath.functions.AnalyticFormulas.blackScholesOptionImpliedVolatility(forwardPrice,maturity,strike,discountBond,price);
				
			}
			return 0.0;
		}
	}

	@Override
	public int getSurfaceSize() {
		
		int size = 0;
		
		Iterator< Map.Entry< Double, OptionSmileData > > iterator = this.surface.entrySet().iterator();
		
		while (iterator.hasNext()) {
			
			Map.Entry< Double, OptionSmileData > pair = 
					(Map.Entry< Double, OptionSmileData >) iterator.next();
			
			OptionSmileData smile = pair.getValue();
			
			size = size + smile.getSize();
			
		}
		
		return size;
		
	}

	@Override
	public int getSmileSize(double maturity) {
		return this.surface.get(maturity).getSize();
	}


	@Override
	public int getNumberOfSmiles() {
		return this.surface.size();
	}
	
}