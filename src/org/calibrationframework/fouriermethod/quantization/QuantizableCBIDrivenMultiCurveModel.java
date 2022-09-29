package org.calibrationframework.fouriermethod.quantization;

import java.util.Arrays;
import java.util.Map;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

import net.finmath.marketdata.model.*;
import net.finmath.marketdata.model.curves.*;
import net.finmath.montecarlo.*;
import net.finmath.stochastic.RandomVariable;

import org.calibrationframework.stochastic.*;
import net.finmath.integration.TrapezoidalRealIntegrator;
import org.calibrationframework.fouriermethod.calibration.models.*;
import org.calibrationframework.functions.ComplexSpecialFunctions;

/**
 * This class, implementing the QuantizableMultiDimProcessInterface interface, deals with the quantization of a multiple yield curve model. 
 * Its underlying class, termed "CBIDrivenMultiCurveModel" or "CBIMultiCurveModelForNegativeRates", implements the MultivariateCalibrableProcessInterface interface,
 * which provides a characteristic function known in closed form of the underlying model,
 * making pricing by quantization towards calibration possible.
 * 
 * @author Szulda Guillaume
 */
public class QuantizableCBIDrivenMultiCurveModel implements QuantizableMultiDimProcessInterface {
	
	private MultivariateCalibrableProcessInterface model;
	private AnalyticModel curves;
	private DiscountCurve initialDC;
	private ForwardCurve initialFC;
	private double maturity;
	private MultiCurveTenor tenor;
	private int level;
	private double upperThreshold;
	private double lowerThreshold;
	private double[] quantizer;
	private double[] companionWeights;

	/**
	 * First constructor, creates an instance of the QuantizableCBIDrivenMultiCurveModel class, 
	 * for the quantization at some level of the CBIDrivenMultiCurveModel model at maturity for some tenor.
	 * The method call "generateQuantizer()" allows for the generation of the quantization grid.
	 * "generateCompanionWeights()" generates all the companions weights for the quantizer,
	 * thus enabling to create the Voronoi quantization of the CBIDrivenMultiCurveModel model at maturity.
	 * This discrete variate can be obtained through the method "getVoronoiQuantization", the corresponding quantization error via "getQuantizationError", 
	 * its discrete state space with "getQuantizationGrid" and its law (probabilities) by "getCompanionWeights".
	 * 
	 * @param maturity
	 * @param level
	 * @param model
	 * @param tenor
	 * @param initialGuess
	 */
	public QuantizableCBIDrivenMultiCurveModel(double maturity, int level, CBIDrivenMultiCurveModel model, MultiCurveTenor tenor, double[] initialGuess) throws IllegalArgumentException {
		if(maturity + tenor.getTenorLength() <= model.getTimeHorizon()) {
				this.curves = model.getAnalyticModel();
				this.initialDC = model.getDiscountCurve();
				this.initialFC = model.getForwardCurve(tenor.getTenorName());
				this.maturity = maturity;
				this.level = level;
				this.model = model;
				this.tenor = tenor;
				this.quantizer = initialGuess;
				this.companionWeights = new double[level];
				generateQuantizer();
				generateCompanionWeights();
		} else {
			throw new IllegalArgumentException("The time at which the variate to quantize is considered must be inside the validity domain of the CBI process.");
		}	
	}
	
	/**
	 * Second constructor, creates an instance of the QuantizableCBIDrivenMultiCurveModel class, 
	 * for the quantization at some level of the CBIDrivenMultiCurveModel model at maturity for some tenor of name tenorName.
	 * The method call "generateQuantizer()" allows for the generation of the quantization grid.
	 * "generateCompanionWeights()" generates all the companions weights for the quantizer,
	 * thus enabling to create the Voronoi quantization of the CBIDrivenMultiCurveModel model at maturity.
	 * This discrete variate can be obtained through the method "getVoronoiQuantization", the corresponding quantization error via "getQuantizationError", 
	 * its discrete state space with "getQuantizationGrid" and its law (probabilities) by "getCompanionWeights".
	 * 
	 * @param maturity
	 * @param level
	 * @param model
	 * @param tenorName
	 * @param initialGuess
	 */
	public QuantizableCBIDrivenMultiCurveModel(double maturity, int level, CBIDrivenMultiCurveModel model, String tenorName, double[] initialGuess) throws IllegalArgumentException {
		this.tenor = new MultiCurveTenor(underlyingToTenor(tenorName), tenorName);
		if(maturity + tenor.getTenorLength() <= model.getTimeHorizon()) {
				this.curves = model.getAnalyticModel();
				this.initialDC = model.getDiscountCurve();
				this.initialFC = model.getForwardCurve(tenor.getTenorName());
				this.maturity = maturity;
				this.level = level;
				this.model = model;
				this.quantizer = initialGuess;
				this.companionWeights = new double[level];
				generateQuantizer();
				generateCompanionWeights();
		} else {
			throw new IllegalArgumentException("The time at which the variate to quantize is considered must be inside the validity domain of the CBI process.");
		}	
	}
	
	/**
	 * Third constructor, creates an instance of the QuantizableCBIDrivenMultiCurveModel class, 
	 * for the quantization at some level of the CBIMultiCurveModelForNegativeRates model at maturity for some tenor.
	 * The method call "generateQuantizer()" allows for the generation of the quantization grid.
	 * "generateCompanionWeights()" generates all the companions weights for the quantizer,
	 * thus enabling to create the Voronoi quantization of the CBIMultiCurveModelForNegativeRates model at maturity.
	 * This discrete variate can be obtained through the method "getVoronoiQuantization", the corresponding quantization error via "getQuantizationError", 
	 * its discrete state space with "getQuantizationGrid" and its law (probabilities) by "getCompanionWeights".
	 * 
	 * @param maturity
	 * @param level
	 * @param model
	 * @param tenor
	 * @param initialGuess
	 */
	public QuantizableCBIDrivenMultiCurveModel(double maturity, int level, CBIMultiCurveModelForNegativeRates model, MultiCurveTenor tenor, double[] initialGuess) throws IllegalArgumentException {
		if(maturity + tenor.getTenorLength() <= model.getTimeHorizon()) {
				this.curves = model.getAnalyticModel();
				this.initialDC = model.getDiscountCurve();
				this.initialFC = model.getForwardCurve(tenor.getTenorName());
				this.maturity = maturity;
				this.level = level;
				this.model = model;
				this.tenor = tenor;
				this.quantizer = initialGuess;
				this.companionWeights = new double[level];
				generateQuantizer();
				generateCompanionWeights();
		} else {
			throw new IllegalArgumentException("The time at which the variate to quantize is considered must be inside the validity domain of the CBI process.");
		}	
	}
	
	/**
	 * Fourth constructor, creates an instance of the QuantizableCBIDrivenMultiCurveModel class, 
	 * for the quantization at some level of the CBIMultiCurveModelForNegativeRates model at maturity for some tenor of name tenorName.
	 * The method call "generateQuantizer()" allows for the generation of the quantization grid.
	 * "generateCompanionWeights()" generates all the companions weights for the quantizer,
	 * thus enabling to create the Voronoi quantization of the CBIMultiCurveModelForNegativeRates model at maturity.
	 * This discrete variate can be obtained through the method "getVoronoiQuantization", the corresponding quantization error via "getQuantizationError", 
	 * its discrete state space with "getQuantizationGrid" and its law (probabilities) by "getCompanionWeights".
	 * 
	 * @param maturity
	 * @param level
	 * @param model
	 * @param tenorName
	 * @param initialGuess
	 */
	public QuantizableCBIDrivenMultiCurveModel(double maturity, int level, CBIMultiCurveModelForNegativeRates model, String tenorName, double[] initialGuess) throws IllegalArgumentException {
		this.tenor = new MultiCurveTenor(underlyingToTenor(tenorName), tenorName);
		if(maturity + tenor.getTenorLength() <= model.getTimeHorizon()) {
				this.curves = model.getAnalyticModel();
				this.initialDC = model.getDiscountCurve();
				this.initialFC = model.getForwardCurve(tenor.getTenorName());
				this.maturity = maturity;
				this.level = level;
				this.model = model;
				this.quantizer = initialGuess;
				this.companionWeights = new double[level];
				generateQuantizer();
				generateCompanionWeights();
		} else {
			throw new IllegalArgumentException("The time at which the variate to quantize is considered must be inside the validity domain of the CBI process.");
		}	
	}
	
	public QuantizableCBIDrivenMultiCurveModel(double maturity, int level, MultivariateCalibrableProcessInterface model, MultiCurveTenor tenor, double[] initialGuess) throws IllegalArgumentException {
		if(model instanceof CBIDrivenMultiCurveModel) {
			if(maturity + tenor.getTenorLength() <= ((CBIDrivenMultiCurveModel)model).getTimeHorizon()) {
				this.curves = ((CBIDrivenMultiCurveModel)model).getAnalyticModel();
				this.initialDC = ((CBIDrivenMultiCurveModel)model).getDiscountCurve();
				this.initialFC = ((CBIDrivenMultiCurveModel)model).getForwardCurve(tenor.getTenorName());
				this.maturity = maturity;
				this.level = level;
				this.model = (CBIDrivenMultiCurveModel)model;
				this.tenor = tenor;
				this.quantizer = initialGuess;
				this.companionWeights = new double[level];
				generateQuantizer();
				generateCompanionWeights();
			} else {
				throw new IllegalArgumentException("The time at which the variate to quantize is considered must be inside the validity domain of the CBI process.");
			}	
		} else if(model instanceof CBIMultiCurveModelForNegativeRates) {
			if(maturity + tenor.getTenorLength() <= ((CBIMultiCurveModelForNegativeRates)model).getTimeHorizon()) {
				this.curves = ((CBIMultiCurveModelForNegativeRates)model).getAnalyticModel();
				this.initialDC = ((CBIMultiCurveModelForNegativeRates)model).getDiscountCurve();
				this.initialFC = ((CBIMultiCurveModelForNegativeRates)model).getForwardCurve(tenor.getTenorName());
				this.maturity = maturity;
				this.level = level;
				this.model = (CBIMultiCurveModelForNegativeRates)model;
				this.tenor = tenor;
				this.quantizer = initialGuess;
				this.companionWeights = new double[level];
				generateQuantizer();
				generateCompanionWeights();
			} else {
				throw new IllegalArgumentException("The time at which the variate to quantize is considered must be inside the validity domain of the CBI process.");
			}		
		} else {
			throw new IllegalArgumentException("The product " + this.getClass()
			+ " cannot be valued against a model " + model.getClass() + "."
			+ "It requires a model of type " + CBIDrivenMultiCurveModel.class + " or a model of type " + CBIMultiCurveModelForNegativeRates.class + ".");
		}
	}
	
	public QuantizableCBIDrivenMultiCurveModel(double maturity, int level, MultivariateCalibrableProcessInterface model, String tenorName, double[] initialGuess) throws IllegalArgumentException {
		if(model instanceof CBIDrivenMultiCurveModel) {
			this.tenor = new MultiCurveTenor(underlyingToTenor(tenorName), tenorName);
			if(maturity + tenor.getTenorLength() <= ((CBIDrivenMultiCurveModel)model).getTimeHorizon()) {
					this.curves = ((CBIDrivenMultiCurveModel)model).getAnalyticModel();
					this.initialDC = ((CBIDrivenMultiCurveModel)model).getDiscountCurve();
					this.initialFC = ((CBIDrivenMultiCurveModel)model).getForwardCurve(tenor.getTenorName());
					this.maturity = maturity;
					this.level = level;
					this.model = (CBIDrivenMultiCurveModel)model;
					this.quantizer = initialGuess;
					this.companionWeights = new double[level];
					generateQuantizer();
					generateCompanionWeights();
			} else {
				throw new IllegalArgumentException("The time at which the variate to quantize is considered must be inside the validity domain of the CBI process.");
			}	
		} else if(model instanceof CBIMultiCurveModelForNegativeRates) {
			this.tenor = new MultiCurveTenor(underlyingToTenor(tenorName), tenorName);
			if(maturity + tenor.getTenorLength() <= ((CBIMultiCurveModelForNegativeRates)model).getTimeHorizon()) {
					this.curves = ((CBIMultiCurveModelForNegativeRates)model).getAnalyticModel();
					this.initialDC = ((CBIMultiCurveModelForNegativeRates)model).getDiscountCurve();
					this.initialFC = ((CBIMultiCurveModelForNegativeRates)model).getForwardCurve(tenor.getTenorName());
					this.maturity = maturity;
					this.level = level;
					this.model = (CBIMultiCurveModelForNegativeRates)model;
					this.quantizer = initialGuess;
					this.companionWeights = new double[level];
					generateQuantizer();
					generateCompanionWeights();
			} else {
				throw new IllegalArgumentException("The time at which the variate to quantize is considered must be inside the validity domain of the CBI process.");
			}		
		} else {
			throw new IllegalArgumentException("The product " + this.getClass()
			+ " cannot be valued against a model " + model.getClass() + "."
			+ "It requires a model of type " + CBIDrivenMultiCurveModel.class + " or a model of type " + CBIMultiCurveModelForNegativeRates.class + ".");
		}
	}
	
	@Override
	public String getUnderlying() {
		return (this.tenor).getTenorName();
	}
	
	public String getTenorName() {
		return (this.tenor).getTenorName();
	}
	
	public double getTenorLength() {
    	return (this.tenor).getTenorLength();
    }
	
	public AnalyticModel getAnalyticModel() {
		return this.curves;
	}
	
	public Map<String,Curve> getInitialCurves() {
		return (this.curves).getCurves();
	}
	
	public DiscountCurve getDiscountCurve() {
		return this.initialDC;
	}
	
	public ForwardCurve getForwardCurve() {
		return this.initialFC;
	}
	
	@Override
	public double getMaturity() {
		return this.maturity;
	}
	
	@Override
	public int getLevel() {
		return this.level;
	}
	
	@Override
	public MultivariateCalibrableProcessInterface getUnderlyingModel() {
		return this.model;
	}

	@Override
	public double[] getQuantizationGrid() {
		return this.quantizer;
	}

	@Override
	public RandomVariable getVoronoiQuantization() {
		return new RandomVariableFromDoubleArray(this.maturity, this.getQuantizationGrid());
	}

	@Override
	public RandomVariable getCompanionWeights() {
		return new RandomVariableFromDoubleArray(this.maturity, this.companionWeights);
	}

	@Override
	public QuantizableMultiDimProcessInterface getCloneForModifiedParameters(double[] parameters) {
		return new QuantizableCBIDrivenMultiCurveModel(this.maturity, this.level, (this.model).getCloneForModifiedParameters(parameters), this.tenor, this.quantizer);
	}
	
	private void generateQuantizer() {
		
		TrapezoidalRealIntegrator mc = new TrapezoidalRealIntegrator(0.01, 100, 50);
		
		/* Initial guess for the Newton-Raphson algorithm. */
		double[] x0 = new double[this.level];
		if(this.quantizer == null) {
			
			double e = Math.abs(this.model.apply(this.maturity, this.getTenorName()).apply(Complex.I.negate()).getReal())*(1.0 / (this.getDiscountCurve()).getDiscountFactor(this.maturity+this.getTenorLength())); // Expectation of the state variable.
		
			if(Double.isNaN(e)) {
				e = 1 + this.getTenorLength()*this.getForwardCurve().getForward(this.getAnalyticModel(), this.getMaturity());
			}
			
	        double dev = Math.sqrt(Math.abs(Math.abs(this.model.apply(this.maturity, this.getTenorName()).apply(Complex.I.negate().multiply(2)).getReal())*(1.0 / (this.getDiscountCurve()).getDiscountFactor(this.maturity+this.getTenorLength())) - e*e)); // Standard deviation.
	        
	        if((dev >= 0.01)&&(dev < e)) {
	        	double h = 2*dev / this.level;
	            x0[0] = e - 1*dev;
	            this.lowerThreshold = x0[0];
	            x0[this.level-1] = e + 1*dev;
	            this.upperThreshold = x0[this.level-1];
	            for(int j = 1; j < this.level-1; j++) {
	            	x0[j] = x0[0] + j*h;
	            }
	        } else if((dev >= e)||(dev < 0.01)||(Double.isNaN(dev))) {
	        	dev = 0.1*e;
	        	double h = 2*dev / this.level;
	            x0[0] = e - 1*dev;
	            this.lowerThreshold = x0[0];
	            x0[this.level-1] = e + 1*dev;
	            this.upperThreshold = x0[this.level-1];
	            for(int j = 1; j < this.level-1; j++) {
	            	x0[j] = x0[0] + j*h;
	            }
	        }
		} else {
			x0 = this.quantizer;
			this.lowerThreshold = x0[0];
			this.upperThreshold = x0[this.level-1];
		}
	     
		double[] v = x0;
		
		/* Start of the Newton-Raphson algorithm. */
		int maxIterations = 10;
		double xTolerance = 0.01;
		for(int l = 0; l < maxIterations; l++) {
					
			double[] w = v;
			
			double[][] d = new double[this.level][this.level];
			double[] g = new double[this.level];
		
			g[0] = (2.0 / Math.PI)*v[0]*mc.integrate(u -> ((((((Complex.I).negate()).multiply(u*Math.log(w[0]))).exp()).multiply( ((model.apply(maturity, tenor.getTenorName())).apply(new Complex(u,0.001)))  )).multiply( (ComplexSpecialFunctions.beta(new Complex(0.001,-u), new Complex(2,0.001)).multiply(new Complex(1,0).subtract(new Complex(lowerThreshold / w[0],0).pow(new Complex(0,-u))).subtract(new Complex(lowerThreshold / w[0],0).pow(new Complex(0,-u)).multiply(new Complex(1-(lowerThreshold / w[0]),0)).divide(ComplexSpecialFunctions.beta(new Complex(0.001,-u), new Complex(1,0.001))))    )  ).subtract(ComplexSpecialFunctions.beta(new Complex(-1,u), new Complex(2,0.001)).multiply(new Complex(1,0).subtract(new Complex(2*w[0] / (w[0]+w[1]),0).pow(new Complex(-1,u))).subtract(new Complex(2*w[0] / (w[0]+w[1]),0).pow(new Complex(-1,u)).multiply(new Complex(1-(2*w[0] / (w[0]+w[1])),0)).divide(ComplexSpecialFunctions.beta(new Complex(-1,u), new Complex(1,0.001))))    )          )) ).getReal()   );
			g[this.level-1] = (2.0 / Math.PI)*v[this.level-1]*mc.integrate(u -> ((((((Complex.I).negate()).multiply(u*Math.log(w[level-1]))).exp()).multiply( ((model.apply(maturity, tenor.getTenorName())).apply(new Complex(u,0.001)))  )).multiply( (ComplexSpecialFunctions.beta(new Complex(0.001,-u), new Complex(2,0.001)).multiply(new Complex(1,0).subtract(new Complex((w[level-2]+w[level-1]) / 2*w[level-1],0).pow(new Complex(0,-u))).subtract(new Complex((w[level-2]+w[level-1]) / 2*w[level-1],0).pow(new Complex(0,-u)).multiply(new Complex(1-((w[level-2]+w[level-1]) / 2*w[level-1]),0)).divide(ComplexSpecialFunctions.beta(new Complex(0.001,-u), new Complex(1,0.001))))    )  ).subtract(ComplexSpecialFunctions.beta(new Complex(-1,u), new Complex(2,0.001)).multiply(new Complex(1,0).subtract(new Complex(w[level-1] / upperThreshold,0).pow(new Complex(-1,u))).subtract(new Complex(w[level-1] / upperThreshold,0).pow(new Complex(-1,u)).multiply(new Complex(1-(w[level-1] / upperThreshold),0)).divide(ComplexSpecialFunctions.beta(new Complex(-1,u), new Complex(2,0.001))))    ) )) ).getReal()   );
	
			d[0][0] = (2.0 / Math.PI)*mc.integrate(u -> ((((((Complex.I).negate()).multiply(u*Math.log(w[0]))).exp()).multiply( ((model.apply(maturity, tenor.getTenorName())).apply(new Complex(u,0.001)))  )).multiply( (ComplexSpecialFunctions.beta(new Complex(0.001,-u), new Complex(1,0.001)).multiply(new Complex(1,0).subtract(new Complex(lowerThreshold / w[0],0).pow(new Complex(0,-u)))   ) ).add( ComplexSpecialFunctions.beta(new Complex(0.001,u), new Complex(1,0.001)).multiply(new Complex(1,0).subtract(new Complex(2*w[0] / (w[1]+w[0]),0).pow(new Complex(0,u)))   )      )) ).getReal()   ) + (-2.0 / ((v[0]+v[1])*Math.PI))*(v[1]-v[0])*0.5*mc.integrate(u -> ((((Complex.I).negate()).multiply(u*Math.log((w[0]+w[1])*0.5)).exp()).multiply( ((model.apply(maturity, tenor.getTenorName())).apply(new Complex(u,0.001)))  )).getReal()   );
			d[this.level-1][this.level-1] = (2.0 / Math.PI)*mc.integrate(u -> ((((((Complex.I).negate()).multiply(u*Math.log(w[level-1]))).exp()).multiply( ((model.apply(maturity, tenor.getTenorName())).apply(new Complex(u,0.001)))  )).multiply( (ComplexSpecialFunctions.beta(new Complex(0.001,-u), new Complex(1,0.001)).multiply(new Complex(1,0).subtract(new Complex((w[level-2]+w[level-1]) / 2*w[level-1],0).pow(new Complex(0,-u)))   ) ).add(ComplexSpecialFunctions.beta(new Complex(0.001,u), new Complex(1,0.001)).multiply(new Complex(1,0).subtract(new Complex(w[level-1] / upperThreshold,0).pow(new Complex(0,u)))   )      )) ).getReal()   ) + (-2.0 / ((v[level-2]+v[level-1])*Math.PI))*(v[level-1]-v[level-2])*0.5*mc.integrate(u -> ((((Complex.I).negate()).multiply(u*Math.log((w[level-2]+w[level-1])*0.5)).exp()).multiply( ((model.apply(maturity, tenor.getTenorName())).apply(new Complex(u,0.001))) )).getReal()   );
			d[0][1] = (-2.0 / ((v[0]+v[1])*Math.PI))*(v[1]-v[0])*0.5*mc.integrate(u -> ((((Complex.I).negate()).multiply(u*Math.log((w[0]+w[1])*0.5)).exp()).multiply( ((model.apply(maturity, tenor.getTenorName())).apply(new Complex(u,0.001)))  )).getReal()   );
			d[this.level-1][this.level-2] = (-2.0 / ((v[this.level-2]+v[this.level-1])*Math.PI))*(v[this.level-1]-v[this.level-2])*0.5*mc.integrate(u -> ((((Complex.I).negate()).multiply(u*Math.log((w[level-2]+w[level-1])*0.5)).exp()).multiply( ((model.apply(maturity, tenor.getTenorName())).apply(new Complex(u,0.001)))  )).getReal()   );
			
			for(int i = 1; i < this.level-1; i++) {
				
				int j = i;
				g[j] = (2.0 / Math.PI)*v[j]*mc.integrate(u -> ((((((Complex.I).negate()).multiply(u*Math.log(w[j]))).exp()).multiply(((model.apply(maturity, tenor.getTenorName())).apply(new Complex(u,0.001)))  )).multiply( (ComplexSpecialFunctions.beta(new Complex(0.001,-u), new Complex(2,0.001)).multiply(new Complex(1,0).subtract(new Complex((w[j]+w[j-1])/2*w[j],0).pow(new Complex(0,-u))).subtract(new Complex((w[j]+w[j-1])/2*w[j],0).pow(new Complex(0,-u)).multiply(new Complex(1-((w[j]+w[j-1])/2*w[j]),0)).divide(ComplexSpecialFunctions.beta(new Complex(0.001,-u), new Complex(1,0.001))))    )  ).subtract(ComplexSpecialFunctions.beta(new Complex(-1,u), new Complex(2,0.001)).multiply(new Complex(1,0).subtract(new Complex(2*w[j] / (w[j]+w[j+1]),0).pow(new Complex(-1,u))).subtract(new Complex(2*w[j] / (w[j]+w[j+1]),0).pow(new Complex(-1,u)).multiply(new Complex(1-(2*w[j] / (w[j]+w[j+1])),0)).divide(ComplexSpecialFunctions.beta(new Complex(-1,u), new Complex(2,0.001))))    ) ) )).getReal()   );
				d[j][j-1] = (-2.0 / ((v[j]+v[j-1])*Math.PI))*(v[j]-v[j-1])*0.5*mc.integrate(u -> ((((Complex.I).negate()).multiply(u*Math.log((w[j-1]+w[j])*0.5)).exp()).multiply(((model.apply(maturity, tenor.getTenorName())).apply(new Complex(u,0.001)))  )).getReal()   );
				d[j][j+1] = (-2.0 / ((v[j+1]+v[j])*Math.PI))*(v[j+1]-v[j])*0.5*mc.integrate(u -> ((((Complex.I).negate()).multiply(u*Math.log((w[j+1]+w[j])*0.5)).exp()).multiply(((model.apply(maturity, tenor.getTenorName())).apply(new Complex(u,0.001)))  )).getReal()   );
				d[j][j] = (2.0 / Math.PI)*mc.integrate(u -> ((((((Complex.I).negate()).multiply(u*Math.log(w[j]))).exp()).multiply(((model.apply(maturity, tenor.getTenorName())).apply(new Complex(u,0.001))) )).multiply( (ComplexSpecialFunctions.beta(new Complex(0.001,-u), new Complex(1,0.001)).multiply(new Complex(1,0).subtract(new Complex((w[j-1]+w[j]) / 2*w[j],0).pow(new Complex(0,-u)))   ) ).add(ComplexSpecialFunctions.beta(new Complex(0.001,u), new Complex(1,0.001)).multiply(new Complex(1,0).subtract(new Complex(2*w[j] / (w[j+1]+w[j]),0).pow(new Complex(0,u)))   )      )) ).getReal()  ) + d[j][j+1] + d[j][j-1];
				
			}
			
			/* Computation of the determinant of the Jacobi tridiagonal matrix d : */
			double[] theta = new double[this.level];
			theta[0] = d[0][0];
			theta[1] = d[0][0]*d[1][1] - d[0][1]*d[1][0];
			for(int n = 2; n < this.level; n++) {
				theta[n] = theta[n-1]*d[n][n] - theta[n-2]*d[n][n-1]*d[n-1][n];
			}
		
			double[] phi = new double[this.level+2];
			phi[this.level+1] = 0;
			phi[this.level] = 1;
			phi[this.level-1] = d[this.level-1][this.level-1];
			for(int n = this.level-2; n >= 0; n--) {
				phi[n] = phi[n+1]*d[n][n] - phi[n+2]*d[n][n+1]*d[n+1][n];
			}
	
			/* Computation of the inverse matrix of the matrix d : */
			double[][] m = new double[this.level][this.level];
			for(int i = 0; i < this.level; i++) {
				if(i == 0) {
					m[0][0] = phi[1] / theta[this.level-1];
				} else {
					for(int j = 0; j <= i; j++) {
						if(j == i) {
							m[i][j] = theta[i-1]*phi[i+1] / theta[this.level-1];
						} else if(j == 0) {
							double p = 1;
							for(int k = j+1; k < i+1; k++) {
								p = p*d[k][k-1];
							}
							m[i][j] = Math.pow(-1, i+j)*p*(phi[i+1] / theta[this.level-1]);
						} else {
							double x = 1;
							for(int k = j+1; k < i+1; k++) {
								x = x*d[k][k-1];
							}
							m[i][j] = Math.pow(-1, i+j)*x*(theta[j-1]*phi[i+1] / theta[this.level-1]);
						}
					}	 
				}
			}
			for(int j = 0; j < this.level ; j++) {
				for(int i = 0; i < j; i++) {
					m[i][j] = m[j][i];
				}
			}
			
			/* Newton-Raphson iteration. */
			double[] r = new double[this.level];
			
			for(int i = 0; i < this.level; i++) {
				
				double sum = 0;
				for(int j = 0; j < this.level; j++) {
					sum = sum + m[i][j]*g[j];
				}
		
				r[i] = Math.abs(v[i]-sum); // We make sure that we get a quantization grid with only positive components after iteration.
				
				/* We substitute the abnormal components (too high or too low) of the grid for better ones (for instance the ones of the initialization point of the algorithm). */
				if(r[i] >= this.upperThreshold || r[i] <= this.lowerThreshold) {
					r[i] = x0[i];
				}
				
			}
			
			Arrays.sort(r); // Every grid has to be sorted and contains only positive numbers before using the algorithm.
			
			RealVector previous = new ArrayRealVector(v);
			RealVector current = new ArrayRealVector(r);
			double distance = previous.subtract(current).getNorm();

			v = r;

			if(distance < xTolerance) {
			break;
			}
				
		} 
				
		this.quantizer = v;
		
	}
	
	private void generateCompanionWeights() {
		
		TrapezoidalRealIntegrator mc = new TrapezoidalRealIntegrator(0.01, 100, 50);
		double sum = 0;
		
		this.companionWeights[0] = (1.0 / (this.getDiscountCurve()).getDiscountFactor(this.maturity+this.getTenorLength()))*Math.abs((this.getDiscountCurve()).getDiscountFactor(this.maturity+this.getTenorLength())*0.5-((1.0 / Math.PI)*mc.integrate(u -> (((((Complex.I).negate()).multiply(u*Math.log((quantizer[0]+quantizer[1])*0.5)).exp()).multiply( ((model.apply(maturity, tenor.getTenorName())).apply(new Complex(u,0.001)))  )).divide(new Complex(0,u))).getReal()   )));
		this.companionWeights[this.level-1] = (1.0 / (this.getDiscountCurve()).getDiscountFactor(this.maturity+this.getTenorLength()))*Math.abs((this.getDiscountCurve()).getDiscountFactor(this.maturity+this.getTenorLength())*0.5+((1.0 / Math.PI)*mc.integrate(u -> (((((Complex.I).negate()).multiply(u*Math.log((quantizer[level-1]+quantizer[level-2])*0.5)).exp()).multiply( ((model.apply(maturity, tenor.getTenorName())).apply(new Complex(u,0.001))) )).divide(new Complex(0,u))).getReal()  )));
		sum = this.companionWeights[0];
		
		/* Use of the Cumulative distribution function of the variable to quantize in order to compute the companion weights. */
		for(int j = 1; j < this.level-1; j++) {
			
			int i = j;
			this.companionWeights[i] = Math.abs(mc.integrate(u -> (((((Complex.I).negate()).multiply(u*Math.log((quantizer[i]+quantizer[i-1])*0.5)).exp()).multiply( ((model.apply(maturity, tenor.getTenorName())).apply(new Complex(u,0.001))) )).divide(new Complex(0,u))).getReal() ) - mc.integrate(u -> (((((Complex.I).negate()).multiply(u*Math.log((this.quantizer[i]+this.quantizer[i+1])*0.5)).exp()).multiply( ((model.apply(maturity, tenor.getTenorName())).apply(new Complex(u,0.001)))   )).divide(new Complex(0,u))).getReal()   ))*(1.0 / Math.PI )*(1.0 / (this.getDiscountCurve()).getDiscountFactor(this.maturity+this.getTenorLength()));
			sum = sum + this.companionWeights[i];
			
		}
		
		sum = sum + this.companionWeights[this.level-1] ;

		for(int i = 0; i < level; i++) {
			this.companionWeights[i] = this.companionWeights[i] / sum;
		}

	}
	
	private double underlyingToTenor(String underlyingName) {
		if(underlyingName.contains("3M")) {
			return 0.25;
		}else if(underlyingName.contains("6M")) {
			return 0.5;
		}else if(underlyingName.contains("12M")) {
			return 1.0;
			
		}else if(underlyingName.contains("1Y")) {
			return 1.0;
			
		}else {
			throw new IllegalArgumentException("Tenor not recognized.");
		}
	}

	@Override
	public double[] getParameterUpperBounds() {
		return (this.model).getParameterUpperBounds();
	}

	@Override
	public double[] getParameterLowerBounds() {
		return (this.model).getParameterLowerBounds();
	}

	@Override
	public double[] getParameters() {
		return (this.model).getParameters();
	}
	
}
