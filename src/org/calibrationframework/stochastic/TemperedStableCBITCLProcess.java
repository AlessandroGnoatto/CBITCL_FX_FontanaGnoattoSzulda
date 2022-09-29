package org.calibrationframework.stochastic;

import java.util.*;
import java.util.function.UnaryOperator;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.complex.Complex;

import org.calibrationframework.fouriermethod.calibration.constraints.*;
import org.calibrationframework.timeseries.FunctionV;

import org.nd4j.common.primitives.Pair;

/**
 * This class represents any tempered stable CBITCL process $(v_t, V_t, L_{V_t})_{t\geq0}$ where $v$ is a tempered stable CBI process.
 * 
 * The whole process is defined up to timeHorizon along a time discretization with a certain numberOfTimeSteps.
 * 
 * The class will be specified by the type of the Lévy process $L$, defined by its Lévy-Khintchine triplet Lévy exponent.
 * 
 * The main feature will be given by the computation of the joint Laplace-Fourier transform of the process, 
 * which is by the method getLaplaceFourierTransform(double time, Complex u1, Complex u2, Complex u3).
 * 
 * @author Szulda Guillaume
 */
public abstract class TemperedStableCBITCLProcess implements CBITCLProcess {

	protected final double timeHorizon;
	protected final int numberOfTimeSteps;
	
	protected double v0;
	
	protected double beta;
	
	protected double b;
	protected double sigma;
	protected double eta;
	protected double theta;
	protected double alpha;
	
	protected UnaryOperator<Complex> branchingMechanism;
	protected UnaryOperator<Complex> immigrationRate;
	
	protected ScalarParameterInformationInterface v0Info;
	
	protected ScalarParameterInformationInterface betaInfo;
	
	protected ScalarParameterInformationInterface bInfo;
	protected ScalarParameterInformationInterface sigmaInfo;
	protected ScalarParameterInformationInterface etaInfo;
	protected ScalarParameterInformationInterface thetaInfo;
	protected ScalarParameterInformationInterface alphaInfo;
	
	/*
	 * Upper and lower bounds are collected here for convenience:
	 * such vectors are then passed to the factory of the optimization algorithm.
	 * In this way we guarantee consistency between the constraints in the model
	 * and the constraints in the optimizer factory.
	 */
	private final double[] parameterLowerBounds;
	private final double[] parameterUpperBounds;
	
	/**
	 * First constructor creates a textbook specification of tempered stable CBITCL process $(v_t, V_t, L_{V_t})_{t\geq0}$,
	 * where $v$ defines a tempered-stable CBI process. In particular, we have $\alpha \in (1, 2)$.
	 * 
	 * @param timeHorizon
	 * @param numberOfTimeSteps
	 * @param v0
	 * @param beta
	 * @param b
	 * @param sigma
	 * @param eta
	 * @param theta
	 * @param alpha
	 */
	public TemperedStableCBITCLProcess(double timeHorizon, int numberOfTimeSteps, double v0, double beta, 
			double b, double sigma, double eta, double theta, double alpha) throws IllegalArgumentException {
		
		if(theta > 0.0 && eta > 0.0 && alpha > 1.0 && alpha < 2.0) {
			
			this.timeHorizon = timeHorizon;
			this.numberOfTimeSteps = numberOfTimeSteps;
			
			this.v0 = v0;
			
			this.beta = beta;
			
			this.b = b;
			this.sigma = sigma;
			this.eta = eta;
			this.theta = theta;
			this.alpha = alpha;

			/* The next one will determine the information on the parameters that was predetermined for this specification of the model and that will be used for calibration,
			 * which consists of the constraints that make the model admissible and whether one parameter has to be calibrated or not. */
			generateInformationForParameters();
			this.parameterLowerBounds = extractLowerBounds();
			this.parameterUpperBounds = extractUpperBounds();
			
			//Computation of the branching mechanism and immigration rate of the model:
			this.branchingMechanism = x -> {
				Complex a = x.multiply(-this.b);
				Complex r = (x.multiply(x)).multiply(0.5*Math.pow(this.sigma, 2));
				Complex l = ( x.multiply( this.alpha*Math.pow(this.theta, this.alpha-1) ) ).subtract( Math.pow(this.theta, this.alpha) );
				Complex c = ( new Complex(this.theta).subtract(x) ).pow(this.alpha);
				return (a.add(r)).add( ( l.add(c) ).multiply( Math.pow(this.eta, this.alpha) ) );
			};
			
			this.immigrationRate = x -> x.multiply(this.beta);	
			
		} else {
			
			throw new IllegalArgumentException("We must have alpha included in (1, 2) as well as theta, eta > 0");
			
		}
				
	}
	
	/**
	 * This second constructor creates a tempered stable CBITCL process $(v_t, V_t, L_{V_t})_{t\geq0}$ with more general information on parameters.
	 * 
	 * @param timeHorizon
	 * @param numberOfTimeSteps
	 * @param v0
	 * @param beta
	 * @param b
	 * @param sigma
	 * @param eta
	 * @param theta
	 * @param alpha
	 * @param v0Info
	 * @param betaInfo
	 * @param bInfo
	 * @param sigmaInfo
	 * @param etaInfo
	 * @param thetaInfo
	 * @param alphaInfo
	 */
	public TemperedStableCBITCLProcess(double timeHorizon, int numberOfTimeSteps, double v0, double beta, double b, double sigma, double eta, double theta, double alpha,
			ScalarParameterInformationInterface v0Info, ScalarParameterInformationInterface betaInfo, 
			ScalarParameterInformationInterface bInfo, ScalarParameterInformationInterface sigmaInfo, ScalarParameterInformationInterface etaInfo, 
			ScalarParameterInformationInterface thetaInfo, ScalarParameterInformationInterface alphaInfo) throws IllegalArgumentException {
		
		if(theta > 0.0 && eta > 0.0 && alpha > 1.0 && alpha < 2.0) {
			
			this.timeHorizon = timeHorizon;
			this.numberOfTimeSteps = numberOfTimeSteps;
			
			this.v0 = v0;
			
			this.beta = beta;
			
			this.b = b;
			this.sigma = sigma;
			this.eta = eta;
			this.theta = theta;
			this.alpha = alpha;
				
			this.v0Info = v0Info;
			
			this.betaInfo = betaInfo;
			
			this.bInfo = bInfo;
			this.sigmaInfo = sigmaInfo;
			this.etaInfo = etaInfo;
			this.thetaInfo = thetaInfo;
			this.alphaInfo = alphaInfo;
			
			this.parameterLowerBounds = extractLowerBounds();
			this.parameterUpperBounds = extractUpperBounds();
	
			//Computation of the branching mechanism and immigration rate of the model:
			this.branchingMechanism = x -> {
				Complex a = x.multiply(-this.b);
				Complex r = (x.multiply(x)).multiply(0.5*Math.pow(this.sigma, 2));
				Complex l = ( x.multiply( this.alpha*Math.pow(this.theta, this.alpha-1) ) ).subtract( Math.pow(this.theta, this.alpha) );
				Complex c = ( new Complex(this.theta).subtract(x) ).pow(this.alpha);
				return (a.add(r)).add( ( l.add(c) ).multiply( Math.pow(this.eta, this.alpha) ) );
			};
			
			this.immigrationRate = x -> x.multiply(this.beta);
			
		} else {
			
			throw new IllegalArgumentException("We must have alpha included in (1, 2) as well as theta, eta > 0");
			
		}
		
	}
	
	@Override
	public double getTimeHorizon() {
		return this.timeHorizon;
	}
	
	@Override 
	public int getNumberOfTimeSteps() {
		return this.numberOfTimeSteps;
	}
	
	@Override
	public int getNumberOfParameters() {
		return 7;
	}
	
	@Override
	public int getDimension() {
		return 3;
	}
	
	@Override
    public double getInitialValueOfCBI() {
    	return this.v0;
    }
	
	public double getBeta() {
		return this.beta;
	}
	
	public double getb() {
		return this.b;
	}
	
	public double getSigma() {
		return this.sigma;
	}
	
	public double getEta() {
		return this.eta;
	}
	
	public double getTheta() {
		return this.theta;
	}
	
	public double getAlpha() {
		return this.alpha;
	}
	
	@Override
	public UnaryOperator<Complex> getBranchingMechanism() {
		return this.branchingMechanism;
	}
	
	@Override
	public UnaryOperator<Complex> getImmigrationRate() {
		return this.immigrationRate;
	}
	
	@Override
	public abstract UnaryOperator<Complex> getLevyExponent();
	
	@Override
    public abstract FunctionV getFunctionV(Complex u1, Complex u2, Complex u3) throws IllegalArgumentException;
	
	@Override
	public abstract Complex getLaplaceFourierTransform(double time, Complex u1, Complex u2, Complex u3) throws IllegalArgumentException;
	
	@Override
	public abstract TemperedStableCBITCLProcess getCloneForModifiedParameters(double[] parameters);
		
	@Override
	public double[] getParameterLowerBounds() {
		return parameterLowerBounds;
	}
	
	@Override
	public double[] getParameterUpperBounds() {
		return parameterUpperBounds;
	}
	
	/**
	 * This method allows to impose the usual constraints on the parameters of the model.
	 * It also informs whether one parameter or another has to be calibrated or not.
	 **/
	private void generateInformationForParameters() {
		
		this.v0Info = new ScalarParameterInformation(true, new PositivityConstraint());
		
		this.betaInfo = new ScalarParameterInformation(true, new PositivityConstraint());
		
		this.bInfo = new ScalarParameterInformation(true);
		this.sigmaInfo = new ScalarParameterInformation(true, new PositivityConstraint());
		this.etaInfo = new ScalarParameterInformation(true, new PositivityConstraint());
		this.thetaInfo = new ScalarParameterInformation(true, new PositivityConstraint());
		this.alphaInfo = new ScalarParameterInformation(true, new BoundConstraint(1, 2));
	
	}
	
	private double[] extractUpperBounds() {
		
		double[] upperBounds = new double[this.getNumberOfParameters()];
		double threshold = 1E6;
		
		upperBounds[0] = this.v0Info.getConstraint().getUpperBound() > threshold ? threshold : this.v0Info.getConstraint().getUpperBound();
		
		upperBounds[1] = this.betaInfo.getConstraint().getUpperBound() > threshold ? threshold : this.betaInfo.getConstraint().getUpperBound();
		
		upperBounds[2] = this.bInfo.getConstraint().getUpperBound() > threshold ? threshold : this.bInfo.getConstraint().getUpperBound();
		upperBounds[3] = this.sigmaInfo.getConstraint().getUpperBound() > threshold ? threshold : this.sigmaInfo.getConstraint().getUpperBound();
		upperBounds[4] = this.etaInfo.getConstraint().getUpperBound() > threshold ? threshold : this.etaInfo.getConstraint().getUpperBound();
		upperBounds[5] = this.thetaInfo.getConstraint().getUpperBound() > threshold ? threshold : this.thetaInfo.getConstraint().getUpperBound();
		upperBounds[6] = this.alphaInfo.getConstraint().getUpperBound() > threshold ? threshold : this.alphaInfo.getConstraint().getUpperBound();
		
		return upperBounds;
		
	}

	private double[] extractLowerBounds() {

		double[] lowerBounds = new double[this.getNumberOfParameters()];
		double threshold = -1E6;
		
		lowerBounds[0] = this.v0Info.getConstraint().getLowerBound() < threshold ? threshold : this.v0Info.getConstraint().getLowerBound();
		
		lowerBounds[1] = this.betaInfo.getConstraint().getLowerBound() < threshold ? threshold : this.betaInfo.getConstraint().getLowerBound();
		
		lowerBounds[2] = this.bInfo.getConstraint().getLowerBound() < threshold ? threshold : this.bInfo.getConstraint().getLowerBound();
		lowerBounds[3] = this.sigmaInfo.getConstraint().getLowerBound() < threshold ? threshold : this.sigmaInfo.getConstraint().getLowerBound();
		lowerBounds[4] = this.etaInfo.getConstraint().getLowerBound() < threshold ? threshold : this.etaInfo.getConstraint().getLowerBound();
		lowerBounds[5] = this.thetaInfo.getConstraint().getLowerBound() < threshold ? threshold : this.thetaInfo.getConstraint().getLowerBound();
		lowerBounds[6] = this.alphaInfo.getConstraint().getLowerBound() < threshold ? threshold : this.alphaInfo.getConstraint().getLowerBound();
		
		return lowerBounds;
		
	}
	
	@Override
	public double[] getParameters() {
		  
		  List<Double> params = new ArrayList<Double>();
		  
		  params.add(this.getInitialValueOfCBI());
		  
		  params.add(this.getBeta());
		  
		  params.add(this.getb());
		  params.add(this.getSigma());
		  params.add(this.getEta());
		  params.add(this.getTheta());
		  params.add(this.getAlpha());
		  
		  Double[] array = params.toArray(new Double[params.size()]);
		  
		  return ArrayUtils.toPrimitive(array);
		  
	}
	
	/**
	 * This method will be implemented by each specification. 
	 * It provides a couple consisting of a sample copy of the textbook specification of the tempered stable CBITCL process and its associated set of parameters.
	 * The parameters will be generated randomly by means of the Mersenne Twister pseudo-random number generator.
	 * 
	 * @param sizeDataSet
	 */
	public abstract Pair< List<Double>, TemperedStableCBITCLProcess > generateSamplePair();
	
}
