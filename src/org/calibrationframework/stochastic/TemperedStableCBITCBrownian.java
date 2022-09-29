package org.calibrationframework.stochastic;

import java.util.*;
import java.util.function.UnaryOperator;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.random.*;
import org.calibrationframework.fouriermethod.calibration.constraints.ScalarParameterInformationInterface;
import org.calibrationframework.timeseries.FunctionV;

import org.nd4j.common.primitives.Pair;

/**
 * This class specifies TemperedStableCBITCLProcess by choosing a Brownian motion $W$ for the Lévy process $L$.
 * 
 * @author Szulda Guillaume
 */
public class TemperedStableCBITCBrownian extends TemperedStableCBITCLProcess {

	private UnaryOperator<Complex> levyExponent; 

	public TemperedStableCBITCBrownian(double timeHorizon, int numberOfTimeSteps, double v0, double beta, double b, double sigma, double eta, double theta, double alpha) {
		
			super(timeHorizon, numberOfTimeSteps, v0, beta, b, sigma, eta, theta, alpha);	
			
			this.levyExponent = x -> {
				Complex h = x.multiply(x).multiply( 0.5 );
				return h;
			};
	}
	
	public TemperedStableCBITCBrownian(double timeHorizon, int numberOfTimeSteps,
			double v0, double beta, double b, double sigma, double eta, double theta, double alpha, 
			ScalarParameterInformationInterface v0Info,
			ScalarParameterInformationInterface betaInfo, 
			ScalarParameterInformationInterface bInfo, 
			ScalarParameterInformationInterface sigmaInfo, 
			ScalarParameterInformationInterface etaInfo, 
			ScalarParameterInformationInterface thetaInfo, 
			ScalarParameterInformationInterface alphaInfo) {
		
		super(timeHorizon, numberOfTimeSteps, v0, beta, b, sigma, eta, theta, alpha, v0Info, betaInfo, bInfo, sigmaInfo, etaInfo, thetaInfo, alphaInfo);
		
		this.levyExponent = x -> {
			Complex h = x.multiply(x).multiply( 0.5 );
			return h;
		};
		
	}
	
	public void setNewValueOfDrift(double newDrift) {
		
		this.levyExponent = x -> {
			Complex h = x.multiply(x).multiply( 0.5 );
			return h.add(x.multiply(newDrift));
		};
		
	}
	
	@Override
	public UnaryOperator<Complex> getLevyExponent() {
		return this.levyExponent;
	}
	
	@Override
	public FunctionV getFunctionV(Complex u1, Complex u2, Complex u3) throws IllegalArgumentException {
		if(u1.getReal() >= this.theta) {
			throw new IllegalArgumentException("The real part of the complex number u1 must be in the interior of D_1");
		} else {
			return new FunctionV(this.getTimeHorizon(), this.getNumberOfTimeSteps(), this.getBranchingMechanism(), this.getImmigrationRate(), this.getLevyExponent(), u1, u2, u3, this.theta);
		}
	}

	@Override
	public Complex getLaplaceFourierTransform(double time, Complex u1, Complex u2, Complex u3)
			throws IllegalArgumentException {
		if(time <= this.getTimeHorizon()) {
			if(u1.getReal() >= this.theta) {
				throw new IllegalArgumentException("The real part of the complex number u1 must be in the interior of D_1");
			} else {
				FunctionV vFunction = this.getFunctionV(u1, u2, u3);
				return ( vFunction.getIntegralOfImmigrationRate(0.001, time).add( vFunction.getValue(time).multiply( this.getInitialValueOfCBI() ) ) ).exp();
			}
		} else {	
			throw new IllegalArgumentException("The time at which the Laplace-Fourier transform is considered must be lower than the time horizon");	
		}
	}

	@Override
	public TemperedStableCBITCBrownian getCloneForModifiedParameters(double[] parameters) {
		
		/* For each parameter, we check whether it has to be calibrated or not. If so, we replace it for the new one to which we apply the corresponding constraint.
		 * If not, this parameter is not modified.
		 */
		double newV0 = this.v0Info.getIsParameterToCalibrate() == true ? this.v0Info.getConstraint().applyConstraint(parameters[0]) : this.v0;
		
		double newBeta = this.betaInfo.getIsParameterToCalibrate() == true ? this.betaInfo.getConstraint().applyConstraint(parameters[1]) : this.beta;
		
		double newb = this.bInfo.getIsParameterToCalibrate() == true ? this.bInfo.getConstraint().applyConstraint(parameters[2]) : this.b;
		double newSigma = this.sigmaInfo.getIsParameterToCalibrate() == true ? this.sigmaInfo.getConstraint().applyConstraint(parameters[3]) : this.sigma;
		double newEta = this.etaInfo.getIsParameterToCalibrate() == true ? this.etaInfo.getConstraint().applyConstraint(parameters[4]) : this.eta;
		double newTheta = this.thetaInfo.getIsParameterToCalibrate() == true ? this.thetaInfo.getConstraint().applyConstraint(parameters[5]) : this.theta;
		double newAlphaY = this.alphaInfo.getIsParameterToCalibrate() == true ? this.alphaInfo.getConstraint().applyConstraint(parameters[6]) : this.alpha;
			
		return new TemperedStableCBITCBrownian(this.timeHorizon, this.numberOfTimeSteps, newV0, newBeta, newb, newSigma, newEta, newTheta, newAlphaY,
					this.v0Info, this.betaInfo, this.bInfo, this.sigmaInfo, this.etaInfo, this.thetaInfo, this.alphaInfo);
		
	}
	
	@Override
	public Pair< List<Double>, TemperedStableCBITCLProcess > generateSamplePair() {
		
		//Creation of the Mersenne Twister pseudo-random number generator
		RandomGenerator rng = new MersenneTwister();
		RandomDataGenerator rndNumGen = new RandomDataGenerator(rng);
			
		//Creation of the parameter list
		List<Double> params = new ArrayList<Double>(this.getNumberOfParameters());
			
		//Random generation of the parameters of TimeChangedBrownian
		double newV0 = rndNumGen.nextExponential(1.0);
		double newBeta = rndNumGen.nextExponential(1.0);
		double newb = rndNumGen.nextGaussian(0.0, 1.0);
		double newSigma = rndNumGen.nextExponential(1.0);
		double newEta = rndNumGen.nextExponential(1.0);
		double newTheta = rndNumGen.nextExponential(1.0);
		double newAlpha = rndNumGen.nextUniform(1.0, 2.0);
			
		//Creation of the CBI-time-changed Brownian
		TemperedStableCBITCBrownian process = new TemperedStableCBITCBrownian(this.timeHorizon, this.numberOfTimeSteps, 
					newV0, newBeta, newb, newSigma, newEta, newTheta, newAlpha);
			
		params.add(newV0);
		
		params.add(newBeta);
		
		params.add(newb);
		params.add(newSigma);
		params.add(newEta);
		params.add(newTheta);
		params.add(newAlpha);
			
		return Pair.< List<Double>, TemperedStableCBITCLProcess >create(params, process);
		
	}

}
