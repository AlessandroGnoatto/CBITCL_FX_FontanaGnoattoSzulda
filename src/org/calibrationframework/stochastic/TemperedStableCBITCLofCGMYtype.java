package org.calibrationframework.stochastic;

import java.util.ArrayList;
import java.util.List;
import java.util.function.UnaryOperator;

import org.apache.commons.lang3.ArrayUtils;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.random.*;

import org.calibrationframework.fouriermethod.calibration.constraints.*;
import org.calibrationframework.timeseries.FunctionV;

import org.nd4j.common.primitives.Pair;

/**
 * This class specifies TemperedStableCBITCLProcess by choosing a CGMY model for the Lévy process $L$.
 * 
 * It then corresponds to the class of tempered stable CBITCL processes of CGMY type.
 * 
 * We also fix $C = \frac{1}{\Gamma(-Y)}$ for normalization purposes.
 * 
 * @author Szulda Guillaume
 */
public class TemperedStableCBITCLofCGMYtype extends TemperedStableCBITCLProcess {
	
	private final double P; //P > 1 needed for pricing by Carr & Madan.

	private double betaL;
	private double G;
	private double M;
	private double Y;
	
	private UnaryOperator<Complex> levyExponent; 
	
	private ScalarParameterInformationInterface betaLInfo;
	
	private ScalarParameterInformationInterface GInfo;
	private ScalarParameterInformationInterface MInfo;
	private ScalarParameterInformationInterface YInfo;
	
	/*
	 * Upper and lower bounds are collected here for convenience:
	 * such vectors are then passed to the factory of the optimization algorithm.
	 * In this way we guarantee consistency between the constraints in the model
	 * and the constraints in the optimizer factory.
	 */
	private final double[] parameterLowerBounds;
	private final double[] parameterUpperBounds;

	/**
	 * First constructor creates a textbook specification of the tempered stable CBITCL process of CGMY type $(v_t, V_t, L_{V_t})_{t\geq0}$.
	 * 
	 * @param P
	 * @param timeHorizon
	 * @param numberOfTimeSteps
	 * @param v0
	 * @param beta
	 * @param b
	 * @param sigma
	 * @param eta
	 * @param theta
	 * @param alpha
	 * @param betaL
	 * @param G
	 * @param M
	 * @param Y
	 */
	public TemperedStableCBITCLofCGMYtype(double P, double timeHorizon, int numberOfTimeSteps, double v0, double beta, 
			double b, double sigma, double eta, double theta, double alpha,
			double betaL, double G, double M, double Y) throws IllegalArgumentException {
	
		super(timeHorizon, numberOfTimeSteps, v0, beta, b, sigma, eta, theta, alpha);	
		
		if(P > 1.0 && Y > 1.0 && Y < 2.0) {

			this.P = P;
			
			this.betaL = betaL;
					
			this.G = G;
			this.M = M;
			this.Y = Y;
					
			/* The next one will determine the information on the parameters that was predetermined for this specification of the model that will be used for calibration,
			 * which consists of the constraints that make the model admissible and whether one parameter has to be calibrated or not. 
			 */
			generateInformationForParameters();
			this.parameterLowerBounds = extractLowerBounds();
			this.parameterUpperBounds = extractUpperBounds();	
			
			if(this.G <= this.P) //Corresponds to having G > P.
				this.G = this.P + 1E-9;
			
			this.levyExponent = x -> { 
					Complex f = x.multiply( this.Y*( Math.pow(this.M, this.Y-1) - Math.pow(this.G, this.Y-1) ) );
					Complex k = f.add( ( ( x.add(this.G) ).pow(this.Y) ).subtract( Math.pow(this.G, this.Y) ) );         
					Complex l = k.add( ( ( x.negate().add(this.M) ).pow(this.Y) ).subtract( Math.pow(this.M, this.Y) ) ); 
					return x.multiply(this.betaL).add(l);
			};
			
		} else {
			
			throw new IllegalArgumentException("P must be greater than one and Y strictly included between one and two");
			
		}	
	
	}
	
	/**
	 * This second constructor creates a tempered stable CBITCL process of CGMY type with more general information on parameters.
	 * 
	 * @param P
	 * @param timeHorizon
	 * @param numberOfTimeSteps
	 * @param v0
	 * @param beta
	 * @param b
	 * @param sigma
	 * @param eta
	 * @param theta
	 * @param alpha
	 * @param betaL
	 * @param G
	 * @param M
	 * @param Y
	 * @param v0Info
	 * @param betaInfo
	 * @param bInfo
	 * @param sigmaInfo
	 * @param etaInfo
	 * @param thetaInfo
	 * @param alphaInfo
	 * @param betaLInfo
	 * @param GInfo
	 * @param MInfo
	 * @param YInfo
	 */
	public TemperedStableCBITCLofCGMYtype(double P, double timeHorizon, int numberOfTimeSteps, double v0, double beta, 
			double b, double sigma, double eta, double theta, double alpha, double betaL, double G, double M, double Y,
			ScalarParameterInformationInterface v0Info, ScalarParameterInformationInterface betaInfo, 
			ScalarParameterInformationInterface bInfo, ScalarParameterInformationInterface sigmaInfo, ScalarParameterInformationInterface etaInfo, 
			ScalarParameterInformationInterface thetaInfo, ScalarParameterInformationInterface alphaInfo, ScalarParameterInformationInterface betaLInfo, 
			ScalarParameterInformationInterface GInfo, ScalarParameterInformationInterface MInfo, ScalarParameterInformationInterface YInfo) throws IllegalArgumentException {
		
		super(timeHorizon, numberOfTimeSteps, v0, beta, b, sigma, eta, theta, alpha, v0Info, betaInfo, bInfo, sigmaInfo, etaInfo, thetaInfo, alphaInfo);
		
		if(P > 1.0 && Y > 1.0 && Y < 2.0) {
			
			this.P = P;
	
			this.betaL = betaL;
			
			this.G = G;
			this.M = M;
			this.Y = Y;
			
			this.betaLInfo = betaLInfo;
			
			this.GInfo = GInfo;
			this.MInfo = MInfo;
			this.YInfo = YInfo;
			
			this.parameterLowerBounds = extractLowerBounds();
			this.parameterUpperBounds = extractUpperBounds();
			
			if(this.G <= this.P) //Corresponds to having G > P.
				this.G = this.P + 1E-9;
			
			this.levyExponent = x -> { 
				Complex f = x.multiply( this.Y*( Math.pow(this.M, this.Y-1) - Math.pow(this.G, this.Y-1) ) );
				Complex k = f.add( ( ( x.add(this.G) ).pow(this.Y) ).subtract( Math.pow(this.G, this.Y) ) );         
				Complex l = k.add( ( ( x.negate().add(this.M) ).pow(this.Y) ).subtract( Math.pow(this.M, this.Y) ) ); 
				return x.multiply(this.betaL).add(l);
			};
			
			} else {
			
			throw new IllegalArgumentException("P must be greater than one and Y strictly included between one and two");
			
		}	
			
	}
	
	public double getP() {
		return this.P;
	}
	
	@Override
	public int getNumberOfParameters() {
		return super.getNumberOfParameters()+4;
	}
	
	public double getBetaL() {
		return this.betaL;
	}
	
	public double getG() {
		return this.G;
	}
	
	public double getM() {
		return this.M;
	}
	
	public double getY() {
		return this.Y;
	}
	
	@Override
	public UnaryOperator<Complex> getLevyExponent() {
		return this.levyExponent;
	}
	
	@Override
    public FunctionV getFunctionV(Complex u1, Complex u2, Complex u3) throws IllegalArgumentException {
		if( u1.getReal() >= this.theta || u3.getReal() >= this.M || u3.getReal() <= -this.G ) {
			throw new IllegalArgumentException("The real part of the complex number u1 must be in the interior of D_1"
					+ " and that of u3 must be in the interior of D_2");
		} else {
			return new FunctionV(this.getTimeHorizon(), this.getNumberOfTimeSteps(), this.getBranchingMechanism(), this.getImmigrationRate(), this.getLevyExponent(), u1, u2, u3, this.theta);
		}
    }
	
	@Override
	public Complex getLaplaceFourierTransform(double time, Complex u1, Complex u2, Complex u3) throws IllegalArgumentException {
		if(time <= this.getTimeHorizon()) {
			if( u1.getReal() >= this.theta || u3.getReal() >= this.M || u3.getReal() <= -this.G ) {
				throw new IllegalArgumentException("The real part of the complex number u1 must be in the interior of D_1"
						+ " and that of u3 must be in the interior of D_2");
			} else {
				FunctionV vFunction = this.getFunctionV(u1, u2, u3);
				return ( vFunction.getIntegralOfImmigrationRate(0.001, time).add( vFunction.getValue(time).multiply( this.getInitialValueOfCBI() ) ) ).exp();
			}
		} else {	
			throw new IllegalArgumentException("The time at which the Laplace-Fourier transform is considered must be lower than the time horizon");	
		}
	}
	
	@Override
	public TemperedStableCBITCLofCGMYtype getCloneForModifiedParameters(double[] parameters) {
		
		/* For each parameter, we check whether it has to be calibrated or not. If so, we replace it for the new one to which we apply the corresponding constraint.
		 * If not, this parameter is not modified.
		 */
		double newV0 = this.v0Info.getIsParameterToCalibrate() == true ? this.v0Info.getConstraint().applyConstraint(parameters[0]) : this.v0;
		
		double newBeta = this.betaInfo.getIsParameterToCalibrate() == true ? this.betaInfo.getConstraint().applyConstraint(parameters[1]) : this.beta;
		
		double newb = this.bInfo.getIsParameterToCalibrate() == true ? this.bInfo.getConstraint().applyConstraint(parameters[2]) : this.b;
		double newSigma = this.sigmaInfo.getIsParameterToCalibrate() == true ? this.sigmaInfo.getConstraint().applyConstraint(parameters[3]) : this.sigma;
		double newEta = this.etaInfo.getIsParameterToCalibrate() == true ? this.etaInfo.getConstraint().applyConstraint(parameters[4]) : this.eta;
		double newTheta = this.thetaInfo.getIsParameterToCalibrate() == true ? this.thetaInfo.getConstraint().applyConstraint(parameters[5]) : this.theta;
		double newAlpha = this.alphaInfo.getIsParameterToCalibrate() == true ? this.alphaInfo.getConstraint().applyConstraint(parameters[6]) : this.alpha;
		
		double newBetaL = this.betaLInfo.getIsParameterToCalibrate() == true ? this.betaLInfo.getConstraint().applyConstraint(parameters[7]) : this.betaL;
		double newG = this.GInfo.getIsParameterToCalibrate() == true ? this.GInfo.getConstraint().applyConstraint(parameters[8]) : this.G;
		double newM = this.MInfo.getIsParameterToCalibrate() == true ? this.MInfo.getConstraint().applyConstraint(parameters[9]) : this.M;
		double newY = this.YInfo.getIsParameterToCalibrate() == true ? this.YInfo.getConstraint().applyConstraint(parameters[10]) : this.Y;
	
		return new TemperedStableCBITCLofCGMYtype(this.P, this.timeHorizon, this.numberOfTimeSteps, newV0, newBeta, newb, newSigma, newEta, newTheta, newAlpha, 
				newBetaL, newG, newM, newY, this.v0Info, this.betaInfo, this.bInfo, this.sigmaInfo, this.etaInfo, this.thetaInfo, this.alphaInfo, 
				this.betaLInfo, this.GInfo, this.MInfo, this.YInfo);
		
	}
		
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
		
		this.betaLInfo = new ScalarParameterInformation(true);
		this.GInfo = new ScalarParameterInformation(true, new PositivityConstraint());
		this.MInfo = new ScalarParameterInformation(true, new PositivityConstraint());
		this.YInfo = new ScalarParameterInformation(true, new BoundConstraint(1, 2));
	
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
		
		upperBounds[7] = this.betaLInfo.getConstraint().getUpperBound() > threshold ? threshold : this.betaLInfo.getConstraint().getUpperBound();
		upperBounds[8] = this.GInfo.getConstraint().getUpperBound() > threshold ? threshold : this.GInfo.getConstraint().getUpperBound();
		upperBounds[9] = this.MInfo.getConstraint().getUpperBound() > threshold ? threshold : this.MInfo.getConstraint().getUpperBound();
		upperBounds[10] = this.YInfo.getConstraint().getUpperBound() > threshold ? threshold : this.YInfo.getConstraint().getUpperBound();
		
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
		
		lowerBounds[7] = this.betaLInfo.getConstraint().getLowerBound() < threshold ? threshold : this.betaLInfo.getConstraint().getLowerBound();
		lowerBounds[8] = this.GInfo.getConstraint().getLowerBound() < threshold ? threshold : this.GInfo.getConstraint().getLowerBound();
		lowerBounds[9] = this.MInfo.getConstraint().getLowerBound() < threshold ? threshold : this.MInfo.getConstraint().getLowerBound();
		lowerBounds[10] = this.YInfo.getConstraint().getLowerBound() < threshold ? threshold : this.YInfo.getConstraint().getLowerBound();
		
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
		  
		  params.add(this.getBetaL());
		  params.add(this.getG());
		  params.add(this.getM());
		  params.add(this.getY());
		  
		  Double[] array = params.toArray(new Double[params.size()]);
		  
		  return ArrayUtils.toPrimitive(array);
		  
	}

	@Override
	public Pair< List<Double>, TemperedStableCBITCLProcess > generateSamplePair() {

		//Creation of the Mersenne Twister pseudo-random number generator
		RandomGenerator rng = new MersenneTwister();
		RandomDataGenerator rndNumGen = new RandomDataGenerator(rng);

		//Creation of the parameter list 
		List<Double> params = new ArrayList<Double>(this.getNumberOfParameters());
			
		//Random generation of the parameters of TimeChangedCGMYmodel
		double newV0 = rndNumGen.nextExponential(1.0);
		double newBeta = rndNumGen.nextExponential(1.0);
		double newb = rndNumGen.nextGaussian(0.0, 1.0);
		double newSigma = rndNumGen.nextExponential(1.0);
		double newEta = rndNumGen.nextExponential(1.0);
		double newTheta = rndNumGen.nextExponential(1.0);
		double newAlpha = rndNumGen.nextUniform(1.0, 2.0);
		double newBetaL = rndNumGen.nextGaussian(0.0, 1.0);
		double newG = rndNumGen.nextExponential(1.0);
		double newM = rndNumGen.nextExponential(1.0);
		double newY = rndNumGen.nextUniform(1.0, 2.0);
			
		//Creation of the CBI-time-changed CGMY model
		TemperedStableCBITCLofCGMYtype process = new TemperedStableCBITCLofCGMYtype(this.P, this.timeHorizon, this.numberOfTimeSteps, 
					newV0, newBeta, newb, newSigma, newEta, newTheta, newAlpha, newBetaL, newG, newM, newY);
			
		params.add(newV0);
		
		params.add(newBeta);
		
		params.add(newb);
		params.add(newSigma);
		params.add(newEta);
		params.add(newTheta);
		params.add(newAlpha);
		
		params.add(newBetaL);
		params.add(newG);
		params.add(newM);
		params.add(newY);
			
		return Pair.< List<Double>, TemperedStableCBITCLProcess >create(params, process);
		
	}

}
