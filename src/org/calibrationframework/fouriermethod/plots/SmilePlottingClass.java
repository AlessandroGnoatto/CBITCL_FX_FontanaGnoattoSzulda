package org.calibrationframework.fouriermethod.plots;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.*;
import java.util.function.*;

import org.apache.commons.lang3.ArrayUtils;
import org.calibrationframework.exception.CalculationException;
import org.calibrationframework.fouriermethod.models.MultivariateProcessCharacteristicFunctionInterface;
import org.calibrationframework.fouriermethod.products.EuropeanOptionSmileMultiAsset;
import org.calibrationframework.marketdata.model.volatilities.OptionSurfaceData;
import org.calibrationframework.marketdata.model.volatilities.VolatilitySurfaceInterface.QuotingConvention;

import net.finmath.plots.*;

/**
 * This class provides several methods for drawing purposes.
 * 
 * For a chosen model coupled with a pricer along the input surfaces (providing the axes for the plots), it allows to do the following:
 * 
 * 	- generate a scatter graph of the chosen maturity and surface via scatterSmile;
 * 
 * 	- draw the entire curve representing the smile of the chosen maturity for the chosen surface via plotSmile;
 * 
 * 	- draw all the smiles of the chosen surface in the same graph via plotSmiles;
 * 
 * 	- plot directly the entire surface of the input name, for a chosen number of strike points and a chosen number of maturity points via plotSurface.
 * 
 * @author Szulda Guillaume
 *
 */
public class SmilePlottingClass {
	
	private final LinkedHashMap<String, OptionSurfaceData> surfaces;
	private final MultivariateProcessCharacteristicFunctionInterface model;
	private final EuropeanOptionSmileMultiAsset pricer;
	
	public SmilePlottingClass(LinkedHashMap<String, OptionSurfaceData> surfaces, 
			MultivariateProcessCharacteristicFunctionInterface model, EuropeanOptionSmileMultiAsset pricer) {
		
		this.surfaces = surfaces;
		this.model = model;
		this.pricer = pricer;
		
	}
	
	/**
	 * This method generates a scatter graph of the smile of the input maturity for the chosen surface.
	 * 
	 * @param surfaceName
	 * @param maturity
	 * @throws IllegalArgumentException
	 */
	public void scatterSmile(String surfaceName, double maturity) throws IllegalArgumentException {
		
		OptionSurfaceData surface = this.surfaces.get(surfaceName);
		QuotingConvention convention = surface.getQuotingConvention();
		double[] maturities = surface.getMaturities();
		boolean contains = DoubleStream.of(maturities).anyMatch(x -> x == maturity);
		
		if(!contains) {
			throw new IllegalArgumentException("maturity must be part of the maturities of the selected surface");
		}
		
		double[] strikes = surface.getSmile(maturity).getStrikes();
		EuropeanOptionSmileMultiAsset newPricer = this.pricer.getCloneWithModifiedParameters(surfaceName, maturity, strikes);
		
		List<Double> vals = new ArrayList<Double>(strikes.length);
		
		try {
			
			Map<Double, Double> modelPrices = newPricer.getValue(this.model);
					
			for(int i = 0; i < strikes.length; i++) {
			
				if(convention.equals(QuotingConvention.VOLATILITYLOGNORMAL)) {
							
					double forward = surface.getEquityForwardCurve().getDiscountFactor(maturity);
					double optionMaturity = maturity;
					double optionStrike = strikes[i];
					double payoffUnit = surface.getDiscountCurve().getDiscountFactor(maturity);
					double optionValue = modelPrices.get(strikes[i]);
					double implVol = Math.abs(net.finmath.functions.AnalyticFormulas.blackScholesOptionImpliedVolatility(forward, optionMaturity, optionStrike, payoffUnit, optionValue));
					if(implVol < 2.0) {
						vals.add(implVol);
					} else {
						vals.add( surface.getSurface().get(optionMaturity).getSmile().get(optionStrike).getValue() - 1E-6 );
					}
					
				} else if(convention.equals(QuotingConvention.VOLATILITYNORMAL)) {
							
					double forward = surface.getEquityForwardCurve().getDiscountFactor(maturity);
					double optionMaturity = maturity;
					double optionStrike = strikes[i];
					double payoffUnit = surface.getDiscountCurve().getDiscountFactor(maturity);
					double optionValue = modelPrices.get(strikes[i]);
					double implVol = Math.abs(net.finmath.functions.AnalyticFormulas.bachelierOptionImpliedVolatility(forward, optionMaturity, optionStrike, payoffUnit, optionValue));
					if(implVol < 2.0) {
						vals.add(implVol);
					} else {
						vals.add( surface.getSurface().get(optionMaturity).getSmile().get(optionStrike).getValue() - 1E-6 );
					}
							
				} else {
					
					vals.add(modelPrices.get(strikes[i]));
							
				}						
						
			}	
					
		} catch (CalculationException e) {
			e.printStackTrace();
		}
		
		double[] values = ArrayUtils.toPrimitive( vals.toArray(new Double[vals.size()]) );
		
		Plot smileScatter = Plots.createScatter(strikes, values, strikes[0], strikes[strikes.length-1], 10)
				.setTitle("Smile of maturity T = " + maturity + " Surface " + surfaceName)
				.setXAxisLabel("Strikes")
				.setYAxisLabel("Implied volatilities")
				.setIsLegendVisible(true);
		
		try {
			smileScatter.show();
		} catch (Exception e) {
			e.printStackTrace();
		}
			
	}
	
	/**
	 * This method draws the curve representing the smile of the chosen maturity for the chosen surface,
	 * with an input number of points that will be used for displaying the curve.
	 * 
	 * @param surfaceName
	 * @param maturity
	 * @param numberOfPoints
	 * @throws IllegalArgumentException
	 */
	public void plotSmile(String surfaceName, double maturity, int numberOfPoints) throws IllegalArgumentException {
		
		OptionSurfaceData surface = this.surfaces.get(surfaceName);
		QuotingConvention convention = surface.getQuotingConvention();
		double[] maturities = surface.getMaturities();
		boolean contains = DoubleStream.of(maturities).anyMatch(x -> x == maturity);
		
		if(!contains) {
			throw new IllegalArgumentException("maturity must be part of the maturities of the selected surface");
		}
		
		double[] strikes = surface.getSmile(maturity).getStrikes();
		
		DoubleUnaryOperator pricingFunction = strike -> {
			
			double value = 0;
			
			EuropeanOptionSmileMultiAsset newPricer = this.pricer.getCloneWithModifiedParameters(surfaceName,
					maturity, new double[] {strike});
			
			try {
				
				Map<Double, Double> modelPrice = newPricer.getValue(this.model);
				
				if(convention.equals(QuotingConvention.VOLATILITYLOGNORMAL)) {
								
					double forward = surface.getEquityForwardCurve().getDiscountFactor(maturity);
					double optionMaturity = maturity;
					double optionStrike = strike;
					double payoffUnit = surface.getDiscountCurve().getDiscountFactor(maturity);
					double optionValue = modelPrice.get(strike);
					double implVol = Math.abs(net.finmath.functions.AnalyticFormulas.blackScholesOptionImpliedVolatility(forward, optionMaturity, optionStrike, payoffUnit, optionValue));
					if(implVol < 2.0) {
						value = implVol;
					} else {
						double distance = Math.abs(strikes[0] - optionStrike);
						int idx = 0;
						for(int c = 1; c < strikes.length; c++) {
							double cdistance = Math.abs(strikes[c] - optionStrike);
							if(cdistance < distance){
								idx = c;
								distance = cdistance;
							}
						}
						value = surface.getSmile(optionMaturity).getOption(strikes[idx]).getValue() - 1E-6;	
					}
						
				} else if(convention.equals(QuotingConvention.VOLATILITYNORMAL)) {
								
					double forward = surface.getEquityForwardCurve().getDiscountFactor(maturity);
					double optionMaturity = maturity;
					double optionStrike = strike;
					double payoffUnit = surface.getDiscountCurve().getDiscountFactor(maturity);
					double optionValue = modelPrice.get(strike);
					double implVol = Math.abs(net.finmath.functions.AnalyticFormulas.bachelierOptionImpliedVolatility(forward, optionMaturity, optionStrike, payoffUnit, optionValue));
					if(implVol < 2.0) {
						value = implVol;
					} else {
						double distance = Math.abs(strikes[0] - optionStrike);
						int idx = 0;
						for(int c = 1; c < strikes.length; c++) {
							double cdistance = Math.abs(strikes[c] - optionStrike);
							if(cdistance < distance){
								idx = c;
								distance = cdistance;
							}
						}
						value = surface.getSmile(optionMaturity).getOption(strikes[idx]).getValue() - 1E-6;	
					}
								
				} else {
						
					value = modelPrice.get(strike);
								
				}							
						
			} catch (CalculationException e) {
				e.printStackTrace();
			}
			
			return value;
			
		};
		
		Plot smilePlot = new Plot2D(strikes[0], strikes[strikes.length-1], numberOfPoints, pricingFunction)
				.setTitle("Smile of maturity T = " + maturity + " Surface " + surfaceName)
				.setXAxisLabel("Strikes")
				.setYAxisLabel("Implied volatilities")
				.setIsLegendVisible(true);
		
		try {
			smilePlot.show();
		} catch (Exception e) {
			e.printStackTrace();
		}
			
	}
	
	/**
	 * This method draws all the smiles of the chosen surface in the same graph, 
	 * by means of the input number of points for representing the curves.
	 * 
	 * @param surfaceName
	 * @param numberOfPoints
	 * @throws IllegalArgumentException
	 */
	public void plotSmiles(String surfaceName, int numberOfPoints) throws IllegalArgumentException {
		
		OptionSurfaceData surface = this.surfaces.get(surfaceName);
		QuotingConvention convention = surface.getQuotingConvention();
		
		double[] maturities = surface.getMaturities();
		List<Named<DoubleUnaryOperator>> pricingFunctions = new ArrayList<Named<DoubleUnaryOperator>>(maturities.length);
		
		for(int i = 0; i < maturities.length; i++) {
			
			double mat = maturities[i];
			
			double[] strikes = surface.getSmile(mat).getStrikes();
			
			DoubleUnaryOperator pricingFunction = strike -> {
				
				double value = 0;
				
				EuropeanOptionSmileMultiAsset newPricer = this.pricer.getCloneWithModifiedParameters(surfaceName,
						mat, new double[] {strike});
				
				try {
					
					Map<Double, Double> modelPrice = newPricer.getValue(this.model);
					
					if(convention.equals(QuotingConvention.VOLATILITYLOGNORMAL)) {
									
						double forward = surface.getEquityForwardCurve().getDiscountFactor(mat);
						double optionMaturity = mat;
						double optionStrike = strike;
						double payoffUnit = surface.getDiscountCurve().getDiscountFactor(mat);
						double optionValue = modelPrice.get(strike);
						double implVol = Math.abs(net.finmath.functions.AnalyticFormulas.blackScholesOptionImpliedVolatility(forward, optionMaturity, optionStrike, payoffUnit, optionValue));
						if(implVol < 2.0) {
							value = implVol;
						} else {
							double distance = Math.abs(strikes[0] - optionStrike);
							int idx = 0;
							for(int c = 1; c < strikes.length; c++) {
								double cdistance = Math.abs(strikes[c] - optionStrike);
								if(cdistance < distance){
									idx = c;
									distance = cdistance;
								}
							}
							value = surface.getSmile(optionMaturity).getOption(strikes[idx]).getValue() - 1E-6;	
						}
							
					} else if(convention.equals(QuotingConvention.VOLATILITYNORMAL)) {
									
						double forward = surface.getEquityForwardCurve().getDiscountFactor(mat);
						double optionMaturity = mat;
						double optionStrike = strike;
						double payoffUnit = surface.getDiscountCurve().getDiscountFactor(mat);
						double optionValue = modelPrice.get(strike);
						double implVol = Math.abs(net.finmath.functions.AnalyticFormulas.bachelierOptionImpliedVolatility(forward, optionMaturity, optionStrike, payoffUnit, optionValue));
						if(implVol < 2.0) {
							value = implVol;
						} else {
							double distance = Math.abs(strikes[0] - optionStrike);
							int idx = 0;
							for(int c = 1; c < strikes.length; c++) {
								double cdistance = Math.abs(strikes[c] - optionStrike);
								if(cdistance < distance){
									idx = c;
									distance = cdistance;
								}
							}
							value = surface.getSmile(optionMaturity).getOption(strikes[idx]).getValue() - 1E-6;	
						}
									
					} else {
							
						value = modelPrice.get(strike);
									
					}							
							
				} catch (CalculationException e) {
					e.printStackTrace();
				}
				
				return value;
				
			};
			
			pricingFunctions.add( new Named<DoubleUnaryOperator>("T = " + mat, pricingFunction));
			
		}
		
		double[] strikes = surface.getSmile(maturities[maturities.length-1]).getStrikes();
		
		Plot smilePlot = new Plot2D(strikes[0], strikes[strikes.length-1], numberOfPoints, pricingFunctions)
				.setTitle("Surface " + surfaceName + ": Smiles")
				.setXAxisLabel("Strikes")
				.setYAxisLabel("Implied volatilities")
				.setIsLegendVisible(true);
		
		try {
			smilePlot.show();
		} catch (Exception e) {
			e.printStackTrace();
		}
			
	}
	
	/**
	 * This method directly plots the entire surface of the input name, for a chosen number of strike points and a chosen number of maturity points.
	 * 
	 * @param surfaceName
	 * @param numStrikePoints
	 * @param numMaturityPoints
	 * @throws IllegalArgumentException
	 */
	public void plotSurface(String surfaceName, int numStrikePoints, int numMaturityPoints) throws IllegalArgumentException {
		
		OptionSurfaceData surface = this.surfaces.get(surfaceName);
		QuotingConvention convention = surface.getQuotingConvention();
		double[] maturities = surface.getMaturities();
		
		DoubleBinaryOperator pricingFunction = (strike, maturity) -> {
			
			double value = 0;
			
			EuropeanOptionSmileMultiAsset newPricer = this.pricer.getCloneWithModifiedParameters(surfaceName,
					maturity, new double[] {strike});
			
			try {
				
				Map<Double, Double> modelPrice = newPricer.getValue(this.model);
				
				if(convention.equals(QuotingConvention.VOLATILITYLOGNORMAL)) {
								
					double forward = surface.getEquityForwardCurve().getDiscountFactor(maturity);
					double optionMaturity = maturity;
					double optionStrike = strike;
					double payoffUnit = surface.getDiscountCurve().getDiscountFactor(maturity);
					double optionValue = modelPrice.get(strike);
					double implVol = Math.abs(net.finmath.functions.AnalyticFormulas.blackScholesOptionImpliedVolatility(forward, optionMaturity, optionStrike, payoffUnit, optionValue));
					if(implVol < 2.0) {
						value = implVol;
					} else {
						double distance1 = Math.abs(maturities[0] - optionMaturity);
						int idx1 = 0;
						for(int c = 1; c < maturities.length; c++) {
							double cdistance = Math.abs(maturities[c] - optionMaturity);
							if(cdistance < distance1){
								idx1 = c;
								distance1 = cdistance;
							}
						}
						double[] strikes = surface.getSmile(maturities[idx1]).getStrikes();
						double distance2 = Math.abs(strikes[0] - optionStrike);
						int idx2 = 0;
						for(int c = 1; c < strikes.length; c++) {
							double cdistance = Math.abs(strikes[c] - optionStrike);
							if(cdistance < distance2){
								idx2 = c;
								distance2 = cdistance;
							}
						}
						value = surface.getSmile(maturities[idx1]).getOption(strikes[idx2]).getValue() - 1E-6;	
					}
						
				} else if(convention.equals(QuotingConvention.VOLATILITYNORMAL)) {
								
					double forward = surface.getEquityForwardCurve().getDiscountFactor(maturity);
					double optionMaturity = maturity;
					double optionStrike = strike;
					double payoffUnit = surface.getDiscountCurve().getDiscountFactor(maturity);
					double optionValue = modelPrice.get(strike);
					double implVol = Math.abs(net.finmath.functions.AnalyticFormulas.bachelierOptionImpliedVolatility(forward, optionMaturity, optionStrike, payoffUnit, optionValue));
					if(implVol < 2.0) {
						value = implVol;
					} else {
						double distance1 = Math.abs(maturities[0] - optionMaturity);
						int idx1 = 0;
						for(int c = 1; c < maturities.length; c++) {
							double cdistance = Math.abs(maturities[c] - optionMaturity);
							if(cdistance < distance1){
								idx1 = c;
								distance1 = cdistance;
							}
						}
						double[] strikes = surface.getSmile(maturities[idx1]).getStrikes();
						double distance2 = Math.abs(strikes[0] - optionStrike);
						int idx2 = 0;
						for(int c = 1; c < strikes.length; c++) {
							double cdistance = Math.abs(strikes[c] - optionStrike);
							if(cdistance < distance2){
								idx2 = c;
								distance2 = cdistance;
							}
						}
						value = surface.getSmile(maturities[idx1]).getOption(strikes[idx2]).getValue() - 1E-6;	
					}
								
				} else {
						
					value = modelPrice.get(strike);
								
				}						
						
			} catch (CalculationException e) {
				e.printStackTrace();
			}
			
			return value;
			
		};
		
		double[] strikes = surface.getSmile(maturities[maturities.length-1]).getStrikes();
		
		Plot smilePlot = new Plot3D(strikes[0], strikes[strikes.length-1],
				maturities[0], maturities[maturities.length-1],
				numStrikePoints, numMaturityPoints,
				pricingFunction)
				.setTitle("Surface " + surfaceName)
				.setXAxisLabel("Strikes")
				.setYAxisLabel("Maturities")
				.setZAxisLabel("Implied volatilities")
				.setIsLegendVisible(true);
		
		try {
			smilePlot.show();
		} catch (Exception e) {
			e.printStackTrace();
		}
			
	}
	
}
