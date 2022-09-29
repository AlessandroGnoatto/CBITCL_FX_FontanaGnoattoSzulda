package org.calibrationframework.montecarlo.models;

import java.util.Map;
import java.util.function.*;
import java.lang.Math;

import org.calibrationframework.montecarlo.process.*;
import org.calibrationframework.stochastic.*;

import net.finmath.montecarlo.*;
import net.finmath.stochastic.RandomVariable;
import net.finmath.marketdata.model.*;
import net.finmath.marketdata.model.curves.*;
import net.finmath.time.TimeDiscretization;

/**
 * This class stands for the Monte Carlo simulation of a multiple yield curve model whose driving process is a (or a family of) CBI process.
 * @author Szulda Guillaume
 */
public class MonteCarloCBIDrivenMultiCurveModel implements MonteCarloCBIDrivenMultiCurveInterface {
	
	private AnalyticModel curves;
	private MonteCarloCBIProcessInterface mcCBIProcess;
	private MultiCurveTenor[] tenors;
	private DiscountCurve initialDC;
	private ForwardCurve[] initialFC;
	private DoubleBinaryOperator integralOfFunctionL;
	private DoubleUnaryOperator[] functionsC;
	
	/**
	 * First constructor, creates a Monte Carlo simulation of a multiple yield curve model using the simulation of the driving process @param mcCBIProcess,
	 * and also based on the other parameters of the model:
	 * @param curves
	 * @param mcCBIProcess
	 * @param tenors
	 * @throws IllegalArgumentException
	 */
	public MonteCarloCBIDrivenMultiCurveModel(AnalyticModel curves, MonteCarloCBIProcessInterface mcCBIProcess, MultiCurveTenor[] tenors) throws IllegalArgumentException {
		this.mcCBIProcess = mcCBIProcess;
		if((this.mcCBIProcess).getNumberOfComponents() != tenors.length) {
			throw new IllegalArgumentException("The dimension of the underlying CBI process and the number of tenors must be euqal.");
		} else {
			this.curves = curves;
			this.initialDC = (this.curves).getDiscountCurve("discount-EUR-OIS");
			this.tenors = tenors;
			this.initialFC = new ForwardCurve[(this.tenors).length];
			for(int i = 0; i < mcCBIProcess.getNumberOfComponents(); i++) {
				if((this.tenors[i]).getTenorName() == "3M") {
					this.initialFC[i] = (this.curves).getForwardCurve("forward-EUR-3M");
				} else if((this.tenors[i]).getTenorName() == "6M") {
					this.initialFC[i] = (this.curves).getForwardCurve("forward-EUR-6M");
				} else {
					throw new IllegalArgumentException("The curve is not available");
				}
			}
			this.functionsC = new DoubleUnaryOperator[mcCBIProcess.getNumberOfComponents()];
			achieveFitToInitialCurves();
		}	
	}
	
	/**
	 * Second constructor, creates a Monte Carlo simulation of a multiple yield curve model using the simulation of the driving process @param mcCBIProcess,
	 * and also based on the other parameters of the model:
	 * @param curves
	 * @param mcCBIProcess
	 * @param tenorLengths
	 * @param tenorNames
	 * @throws IllegalArgumentException
	 */
	public MonteCarloCBIDrivenMultiCurveModel(AnalyticModel curves, MonteCarloCBIProcessInterface mcCBIProcess, double[] tenorLengths, String[] tenorNames) throws IllegalArgumentException {
		this.mcCBIProcess = mcCBIProcess;
		if((this.mcCBIProcess).getNumberOfComponents() != tenorNames.length || (this.mcCBIProcess).getNumberOfComponents() != tenorLengths.length) {
			throw new IllegalArgumentException("The dimension of the underlying CBI process and the number of tenors must be euqal.");
		} else {
			this.curves = curves;
			this.initialDC = (this.curves).getDiscountCurve("discount-EUR-OIS");
			this.tenors = new MultiCurveTenor[mcCBIProcess.getNumberOfComponents()];
			this.initialFC = new ForwardCurve[this.tenors.length];
			for(int i = 0; i < mcCBIProcess.getNumberOfComponents(); i++) {
				this.tenors[i] = new MultiCurveTenor(tenorLengths[i],tenorNames[i]);
				if((this.tenors[i]).getTenorName() == "3M") {
					this.initialFC[i] = (this.curves).getForwardCurve("forward-EUR-3M");
				} else if((this.tenors[i]).getTenorName() == "6M") {
					this.initialFC[i] = (this.curves).getForwardCurve("forward-EUR-6M");
				} else {
					throw new IllegalArgumentException("The curve is not available");
				}
			}
			this.functionsC = new DoubleUnaryOperator[mcCBIProcess.getNumberOfComponents()];
			achieveFitToInitialCurves();
		}	
	}
	
	/**
	 * Third constructor, creates a Monte Carlo simulation of a multiple yield curve model while simulating the driving process @param mcCBIProcess at once,
	 * doing so via the constructor of the MonteCarloFlowOfTemperedCBIProcess class, the simulation of a flow of tempered processes then.
	 * Here are the rest of the input parameters:
	 * @param seed
	 * @param numberOfPaths
	 * @param timeDiscretization
	 * @param cbiProcess
	 * @param curves
	 * @param mcCBIProcess
	 * @param tenorLengths
	 * @param tenorNames
	 * @throws IllegalArgumentException
	 */
	public MonteCarloCBIDrivenMultiCurveModel(int seed, int numberOfPaths, TimeDiscretization timeDiscretization, FlowOfTemperedAlphaStableCBIprocess cbiProcess, AnalyticModel curves, MonteCarloCBIProcessInterface mcCBIProcess, double[] tenorLengths, String[] tenorNames) throws IllegalArgumentException {
		this.mcCBIProcess = new MonteCarloFlowOfTemperedCBIProcess(seed, numberOfPaths, timeDiscretization, cbiProcess);
		if((this.mcCBIProcess).getNumberOfComponents() != tenorNames.length || (this.mcCBIProcess).getNumberOfComponents() != tenorLengths.length) {
			throw new IllegalArgumentException("The dimension of the underlying CBI process and the number of tenors must be euqal.");
		} else {
			this.curves = curves;
			this.initialDC = (this.curves).getDiscountCurve("discount-EUR-OIS");
			this.tenors = new MultiCurveTenor[mcCBIProcess.getNumberOfComponents()];
			this.initialFC = new ForwardCurve[this.tenors.length];
			for(int i = 0; i < mcCBIProcess.getNumberOfComponents(); i++) {
				this.tenors[i] = new MultiCurveTenor(tenorLengths[i],tenorNames[i]);
				if((this.tenors[i]).getTenorName() == "3M") {
					this.initialFC[i] = (this.curves).getForwardCurve("forward-EUR-3M");
				} else if((this.tenors[i]).getTenorName() == "6M") {
					this.initialFC[i] = (this.curves).getForwardCurve("forward-EUR-6M");
				} else {
					throw new IllegalArgumentException("The curve is not available");
				}
			}
			this.functionsC = new DoubleUnaryOperator[mcCBIProcess.getNumberOfComponents()];
			achieveFitToInitialCurves();
		}	
	}
	
	/**
	 * Fourth and ultimate constructor, creates a Monte Carlo simulation of a multiple yield curve model while simulating the driving process @param mcCBIProcess at once,
	 * doing so via the other constructor of the MonteCarloFlowOfTemperedCBIProcess class, which is by using all the parameters required to set up the flow of tempered CBI processes.
	 * Here are the rest of the input parameters:
	 * @param timeHorizon
	 * @param numberOfTimeSteps
	 * @param seed
	 * @param numberOfPaths
	 * @param timeDiscretization
	 * @param initialValues
	 * @param immigrationRates
	 * @param b
	 * @param sigma
	 * @param eta
	 * @param zeta
	 * @param alpha
	 * @param lambda
	 * @param curves
	 * @param tenorLengths
	 * @param tenorNames
	 * @throws IllegalArgumentException
	 */
	public MonteCarloCBIDrivenMultiCurveModel(double timeHorizon, int numberOfTimeSteps, int seed, int numberOfPaths, TimeDiscretization timeDiscretization, double[] initialValues, double[] immigrationRates, double b, double sigma, double eta, double zeta, double alpha, double[] lambda, AnalyticModel curves, double[] tenorLengths, String[] tenorNames) throws IllegalArgumentException {
		this.mcCBIProcess = new MonteCarloFlowOfTemperedCBIProcess(timeHorizon, numberOfTimeSteps, seed, numberOfPaths, timeDiscretization, initialValues, immigrationRates, b, sigma, eta, zeta, alpha, lambda);
		if((this.mcCBIProcess).getNumberOfComponents() != tenorNames.length || (this.mcCBIProcess).getNumberOfComponents() != tenorLengths.length) {
			throw new IllegalArgumentException("The dimension of the underlying CBI process and the number of tenors must be euqal.");
		} else {
			this.curves = curves;
			this.initialDC = (this.curves).getDiscountCurve("discount-EUR-OIS");
			this.tenors = new MultiCurveTenor[mcCBIProcess.getNumberOfComponents()];
			this.initialFC = new ForwardCurve[this.tenors.length];
			for(int i = 0; i < mcCBIProcess.getNumberOfComponents(); i++) {
				this.tenors[i] = new MultiCurveTenor(tenorLengths[i],tenorNames[i]);
				if((this.tenors[i]).getTenorName() == "3M") {
					this.initialFC[i] = (this.curves).getForwardCurve("forward-EUR-3M");
				} else if((this.tenors[i]).getTenorName() == "6M") {
					this.initialFC[i] = (this.curves).getForwardCurve("forward-EUR-6M");
				} else {
					throw new IllegalArgumentException("The curve is not available");
				}
			} 
			this.functionsC = new DoubleUnaryOperator[mcCBIProcess.getNumberOfComponents()];
			achieveFitToInitialCurves();
		}
	}
	
	@Override 
	public double getTimeHorizon() {
		return this.mcCBIProcess.getTimeHorizon();
	}
	
	@Override
	public int getNumberOfTimeSteps() {
		return this.mcCBIProcess.getNumberOfTimeSteps();
	}
	
	@Override
	public int getDimensionOfCBI() {
		return (this.mcCBIProcess).getNumberOfComponents();
	}
	
	@Override
	public int getNumberOfPaths() {
		return mcCBIProcess.getNumberOfPaths();
	}
	
	@Override
	public TimeDiscretization getTimeDiscretization() {
		return mcCBIProcess.getTimeDiscretization();
	}
	
	@Override
	public double getTime(int timeIndex) {
		return mcCBIProcess.getTime(timeIndex);
	}
	
	@Override
	public int getTimeIndex(double time) {
		return mcCBIProcess.getTimeIndex(time);
	}
	
	@Override
	public MonteCarloCBIProcessInterface getMonteCarloCBIProcess() {
		return this.mcCBIProcess;
	}
	
	@Override
	public double getTenorLength(String tenorName) {
		
		int tenor;
		
		if(tenorName.equals("forward-EUR-3M")) {
			tenor = 0;
		}else if(tenorName.equals("forward-EUR-6M")){
			tenor = 1;
		}else {
			throw new IllegalArgumentException("The curve is not available");
		}
		
		return getTenorLength(tenor);
		
	}
	
	public double getTenorLength(int i) {
	    return (this.tenors[i]).getTenorLength();
	}
	
	@Override
	public AnalyticModel getAnalyticModel() {
		return this.curves;
	}
	
	@Override
	public Map<String,Curve> getInitialCurves() {
		return (this.curves).getCurves();
	}
	
	@Override
	public DiscountCurve getDiscountCurve() {
		return this.initialDC;
	}
	
	@Override
	public ForwardCurve getForwardCurve(String tenorName) {

        int tenor;
		
		if(tenorName.equals("forward-EUR-3M")) {
			tenor = 0;
		}else if(tenorName.equals("forward-EUR-6M")){
			tenor = 1;
		}else {
			throw new IllegalArgumentException("The curve is not available");
		}
		
		return getForwardCurve(tenor);
		
	}
	
	public ForwardCurve getForwardCurve(int i) {
		return this.initialFC[i];
	}
	
	@Override
	public MonteCarloSimulationInterface getCloneWithModifiedData(Map<String, Object> dataModified) {
		return new MonteCarloCBIDrivenMultiCurveModel((AnalyticModel)(dataModified.get("curves")), (MonteCarloCBIProcessInterface)(dataModified.get("mcCBIProcess")), (double[])(dataModified.get("tenorLengths")), (String[])(dataModified.get("tenorNames")));
	}
	
	@Override
	public RandomVariable getRandomVariableForConstant(double value) {
		return new RandomVariableFromDoubleArray(0.0, value);
	}

	@Override
	public RandomVariable getMonteCarloWeights(int timeIndex) {
		return mcCBIProcess.getMonteCarloWeights(timeIndex);
	}

	@Override
	public RandomVariable getMonteCarloWeights(double time) {
		return mcCBIProcess.getMonteCarloWeights(time);
	}

	@Override
	public MonteCarloCBIDrivenMultiCurveInterface getCloneWithModifiedSeed(int seed) {
	    return new MonteCarloCBIDrivenMultiCurveModel(this.curves, (MonteCarloCBIProcessInterface)((this.mcCBIProcess).getCloneWithModifiedSeed(seed)), this.tenors);
	}
	
	@Override
	public RandomVariable getSpreadValue(double time, String tenorName) {
		
        int tenor;
		
		if(tenorName.equals("forward-EUR-3M")) {
			tenor = 0;
		}else if(tenorName.equals("forward-EUR-6M")){
			tenor = 1;
		}else {
			throw new IllegalArgumentException("The curve is not available");
		}
		
		return getSpreadValue(time, tenor);
	}

	public RandomVariable getSpreadValue(double time, int j) {
		RandomVariable s = new RandomVariableFromDoubleArray(time, 0);
		for(int i = 0; i <= j; i++) {
			s = s.add(this.mcCBIProcess.getCBIProcessValue(time, i));
		}
		return (s.add(this.functionsC[j].applyAsDouble(time))).exp();
	}
	
	@Override
	public RandomVariable getZCBond(double time, double maturity) {
		RandomVariable s = new RandomVariableFromDoubleArray(time, 0);
		for(int i = 0; i < this.getDimensionOfCBI(); i++) {
			s = s.add(this.mcCBIProcess.getCBIProcessValue(time, i).mult(this.mcCBIProcess.getCBIProcess().getFunctionsVZero()[i].getValue(maturity - time)).add(this.mcCBIProcess.getCBIProcess().getImmigrationRates()[i]*this.mcCBIProcess.getCBIProcess().getFunctionsVZero()[i].getIntegral(0, maturity - time)));
		}
		return (s.mult(-1).sub(this.integralOfFunctionL.applyAsDouble(time, maturity))).exp();
	}

	@Override
	public RandomVariable getNumeraire(double time) {
		RandomVariable s = new RandomVariableFromDoubleArray(time, 0);
		for(int factor = 0; factor < getDimensionOfCBI(); factor++) {
			RandomVariable sum = new RandomVariableFromDoubleArray(time, 0);
			for(int index = 0; index < getTimeIndex(time); index++) {
				sum = sum.add(this.mcCBIProcess.getCBIProcessValue(index, factor).mult(this.mcCBIProcess.getTimeDiscretization().getTimeStep(index)));
			}
			s = s.add(sum.mult(this.mcCBIProcess.getCBIProcess().getLambda()[factor]));
		}
		return (s.add(this.integralOfFunctionL.applyAsDouble(0, time))).exp();
	}
	
	private void achieveFitToInitialCurves() {
		for(int k = 0; k < this.getDimensionOfCBI(); k++) {
			int i = k;
			this.functionsC[i] = t -> {
				double a = Math.log(1 + this.tenors[i].getTenorLength()*this.initialFC[i].getForward(this.getAnalyticModel(), t)) + 
						Math.log(this.initialDC.getDiscountFactor(t + this.tenors[i].getTenorLength())) - 
						Math.log(this.initialDC.getDiscountFactor(t));
				double sum = 0;
				for(int j = 0; j <= i; j++) {
					sum = sum + 
							this.mcCBIProcess.getCBIProcess().getInitialValues()[j]*(this.mcCBIProcess.getCBIProcess().getFunctionsVMinusOne()[j].getValue(t) - this.mcCBIProcess.getCBIProcess().getFunctionsVZero()[j].getValue(t)) +
							this.mcCBIProcess.getCBIProcess().getImmigrationRates()[j]*(this.mcCBIProcess.getCBIProcess().getFunctionsVMinusOne()[j].getIntegral(0, t) - this.mcCBIProcess.getCBIProcess().getFunctionsVZero()[j].getIntegral(0, t));
				}
				return a + sum;
			};
		}
		this.integralOfFunctionL = (a,b) -> {
			double x = Math.log(this.initialDC.getDiscountFactor(a)) - Math.log(this.initialDC.getDiscountFactor(b));
			double sum = 0;
			for(int i = 0; i < this.getDimensionOfCBI(); i ++) {
				sum = sum + this.mcCBIProcess.getCBIProcess().getImmigrationRates()[i]*this.mcCBIProcess.getCBIProcess().getFunctionsVZero()[i].getIntegral(a, b) +
						this.mcCBIProcess.getCBIProcess().getInitialValues()[i]*(this.mcCBIProcess.getCBIProcess().getFunctionsVZero()[i].getValue(b) - this.mcCBIProcess.getCBIProcess().getFunctionsVZero()[i].getValue(a)); 
			}
			return x - sum;
		};
	}
	
}

