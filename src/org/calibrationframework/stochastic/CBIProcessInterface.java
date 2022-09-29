package org.calibrationframework.stochastic;

import java.util.function.DoubleUnaryOperator;
import java.util.function.UnaryOperator;

import org.apache.commons.math3.complex.Complex;

import org.calibrationframework.timeseries.*;

/**
 * This interface has to be implemented by each (multi-dimensional) CBI process, 
 * which we are going to use within this project. 
 *
 * @author Szulda Guillaume
 */
public interface CBIProcessInterface extends AffineProcessInterface {
	
	public double[] getInitialValues();
	
	public double[] getImmigrationRates();
	
	public double[] getLambda();

	public DoubleUnaryOperator getBranchingMechanism();
	
	public UnaryOperator<Complex> getComplexBranchingMechanism();
	
	public FunctionVZero[] getFunctionsVZero();
	
	public FunctionVMinusOne[] getFunctionsVMinusOne();
	
	public FunctionW[] getFunctionsW(Complex[] u);
	
	@Override
	public CBIProcessInterface getCloneForModifiedParameters(double[] parameters);

}
