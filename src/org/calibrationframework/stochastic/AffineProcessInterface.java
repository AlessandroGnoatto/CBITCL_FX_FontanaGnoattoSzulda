package org.calibrationframework.stochastic;

/**
 * This interface has to be implemented by each Affine process, 
 * which we are going to use within this project. 
 *
 * @author Szulda Guillaume
 */
public interface AffineProcessInterface {
	
	public double getTimeHorizon();
	
	public int getNumberOfTimeSteps();
	
	public int getNumberOfParameters();
	
	public double[] getParameterLowerBounds();
	
	public double[] getParameterUpperBounds();
	
	public double[] getParameters();
	
	public int getDimension();
	
    public AffineProcessInterface getCloneForModifiedParameters(double[] parameters);
    
}
