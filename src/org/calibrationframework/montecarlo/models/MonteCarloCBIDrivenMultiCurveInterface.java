package org.calibrationframework.montecarlo.models;

import java.util.Map;

import net.finmath.marketdata.model.AnalyticModel;
import net.finmath.marketdata.model.curves.*;
import net.finmath.stochastic.RandomVariable;

import org.calibrationframework.montecarlo.process.*;

/**
 * This interface has to be implemented by a class that is intended to represent a Monte Carlo simulation of a multiple yield curve model,
 * then it will have as an attribute the Monte Carlo simulation of its driving process (here a single CBI process or several like a flow),
 * which will be used to compute some data that usually characterize a multiple yield curve model (ZC bonds, multiplicative spreads),
 * as well as the usual features that an interest rate model must have at disposal (the initial term structure for instance).
 * @author Szulda Guillaume
 *
 */
public interface MonteCarloCBIDrivenMultiCurveInterface extends MonteCarloSimulationInterface {
	
	public MonteCarloCBIProcessInterface getMonteCarloCBIProcess();
	
	public double getTimeHorizon();
	
	public int getNumberOfTimeSteps();
	
	public int getDimensionOfCBI();
	
	public double getTenorLength(String tenorName);
	
	public AnalyticModel getAnalyticModel();
	
	public Map<String,Curve> getInitialCurves();
	
	public DiscountCurve getDiscountCurve();
	
	public ForwardCurve getForwardCurve(String tenorName);
	
	public MonteCarloCBIDrivenMultiCurveInterface getCloneWithModifiedSeed(int seed);
	
	public RandomVariable getSpreadValue(double time, String tenorName);
	
	public RandomVariable getZCBond(double time, double maturity);
	
    public RandomVariable getNumeraire(double time);
	
}
