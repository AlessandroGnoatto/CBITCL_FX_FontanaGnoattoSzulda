package org.calibrationframework.stochastic;

import java.util.function.UnaryOperator;

import org.apache.commons.math3.complex.Complex;

import org.calibrationframework.timeseries.*;

/**
 * This interface has to be implemented by each CBITCL process $(v_t, V_t, X_t)_{t\geq0}$,
 * which is driven by a CBI process $v=(v_t)_{t \geq 0}$,
 * its time integral $V_t := \int_0^t {v_s ds}$ and a time-changed Levy process $X_t := L_{V_t}$.
 *
 * @author Szulda Guillaume
 */
public interface CBITCLProcess extends AffineProcessInterface {
	
	public double getInitialValueOfCBI();
	
	public UnaryOperator<Complex> getBranchingMechanism();
	
	public UnaryOperator<Complex> getImmigrationRate();
	
	public UnaryOperator<Complex> getLevyExponent();
	
	public FunctionV getFunctionV(Complex u1, Complex u2, Complex u3);
	
	public Complex getLaplaceFourierTransform(double T, Complex u1, Complex u2, Complex u3);
	
	@Override
	public CBITCLProcess getCloneForModifiedParameters(double[] parameters);
	
}
