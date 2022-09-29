package org.calibrationframework.fouriermethod.calibration.models;

import java.util.List;

import org.calibrationframework.fouriermethod.models.MultivariateProcessCharacteristicFunctionInterface;

import org.nd4j.common.primitives.Pair;

public interface MultivariateCalibrableProcessInterface extends MultivariateProcessCharacteristicFunctionInterface{
	
	/**
	 * The following two methods return the features of the time discretization over which the model is defined.
	 * The time horizon along with the number of time steps.
	 * 
	 */
	double getTimeHorizon();
	
	int getNumberOfTimeSteps();
	
	/**
	 * This method is needed for the corresponding model to be passed to the FFT pricing routine.
	 * @return
	 */
	MultivariateProcessCharacteristicFunctionInterface getCharacteristiFunction();
	
	/**
	 * Calibration substitutes in the model the parameters of the process with calibrated ones.
	 * Market observables such as the initial stock value should not be changed.
	 * @param parameters
	 * @return a clone of the original model with modified parameters.
	 */
	MultivariateCalibrableProcessInterface getCloneForModifiedParameters(double[] parameters);

	/*
	 * Upper and lower bounds have to be collected for them to be passed to the factory of the optimization algorithm.
	 * In this way we guarantee consistency between the constraints in the model
	 * and the constraints in the optimizer factory.
	 */
	double[] getParameterUpperBounds();
	
	double[] getParameterLowerBounds();
	
	double[] getParameters();
	
	/**
	 * The following method generates, for each specification of MultivariateCalibrableProcessInterface,
	 * a sample pair consisting of a randomly-generated parameter set and its associated MultivariateCalibrableProcessInterface model.
	 * 
	 * @param sizeDataSet
	 * @return
	 */
	Pair< List<Double>, MultivariateCalibrableProcessInterface > generateSamplePair();

}
