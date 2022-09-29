package org.calibrationframework.fouriermethod.calibration;

import java.util.*;
import java.util.Map.Entry;

import org.apache.commons.lang3.ArrayUtils;

import org.calibrationframework.exception.CalculationException;
import org.calibrationframework.fouriermethod.calibration.models.MultivariateCalibrableProcessInterface;
import org.calibrationframework.fouriermethod.models.MultivariateProcessCharacteristicFunctionInterface;
import org.calibrationframework.fouriermethod.products.EuropeanOptionSmileMultiAsset;
import org.calibrationframework.marketdata.model.volatilities.*;
import org.calibrationframework.marketdata.model.volatilities.VolatilitySurfaceInterface.QuotingConvention;
import org.calibrationframework.optimizer.*;
import org.nd4j.common.primitives.Pair;

public class CurrencyOptionCalibrationProblem {

	private final LinkedHashMap<String, OptionSurfaceData> surfaces;
	private final MultivariateCalibrableProcessInterface model; //Pricing model
	private final OptimizerFactoryInterface optimizerFactory; //construct the instance of the optimization algorithm inside the class.
	private final EuropeanOptionSmileMultiAsset pricer; //How do we compute prices: Carr Madan, Cos, Conv, Lewis...
	
	private final Pair< double[], MultivariateCalibrableProcessInterface > samplePair;

	//Optimizer parameters
	private final double[] initialParameters;
	private final double[] lowerBound;
	private final double[] upperBound;
	private final double[] parameterStep;
	
	public CurrencyOptionCalibrationProblem(LinkedHashMap<String, OptionSurfaceData> surfaces, MultivariateCalibrableProcessInterface model, 
			OptimizerFactoryInterface optimizerFactory, EuropeanOptionSmileMultiAsset pricer) throws IllegalArgumentException {
		
		this.surfaces = surfaces;
		this.model = model;
		this.optimizerFactory = optimizerFactory;
		this.pricer = pricer;
		
		this.samplePair = generateSamplePair();
		
		this.initialParameters = this.samplePair.getFirst();
			
		this.lowerBound = this.samplePair.getSecond().getParameterLowerBounds();
		this.upperBound = this.samplePair.getSecond().getParameterUpperBounds();
			
		this.parameterStep = new double[this.initialParameters.length];
		for(int i = 0; i < this.initialParameters.length; i++) 
			this.parameterStep[i] = 0.01;
	
	}
	
	public CurrencyOptionCalibrationProblem(LinkedHashMap<String, OptionSurfaceData> surfaces, MultivariateCalibrableProcessInterface model, 
			OptimizerFactoryInterface optimizerFactory, EuropeanOptionSmileMultiAsset pricer,
			double[] initialParameters, double[] parameterStep) throws IllegalArgumentException {
		
		this.surfaces = surfaces;
		this.model = model;
		this.optimizerFactory = optimizerFactory;
		this.pricer = pricer;
		
		this.samplePair = generateSamplePair();
		
		this.initialParameters = initialParameters;
					
		this.lowerBound = this.samplePair.getSecond().getParameterLowerBounds();
		this.upperBound = this.samplePair.getSecond().getParameterUpperBounds();
					
		this.parameterStep = parameterStep;
	
	}

	public OptimizationResult runCalibration() throws SolverException {
		OptimizerInterface.ObjectiveFunction objectiveFunction = new OptimizerInterface.ObjectiveFunction() {

			@Override
			public void setValues(double[] parameters, double[] values) throws SolverException {
				
				long tstart = System.currentTimeMillis();
				//We change the parameters of the model
				MultivariateCalibrableProcessInterface newModel = model.getCloneForModifiedParameters(parameters);
				MultivariateProcessCharacteristicFunctionInterface newModelFourier = newModel.getCharacteristiFunction();

				Iterator<Entry<String, OptionSurfaceData>> iterator = surfaces.entrySet().iterator();
				final ArrayList<Double> vals = new ArrayList<>();
				
				//loop over different surfaces
				while (iterator.hasNext()) {
					Map.Entry<String, OptionSurfaceData> pair = (Entry<String, OptionSurfaceData>) iterator.next();
					String underlying = pair.getKey();
					OptionSurfaceData surface = pair.getValue();
					
					int numberOfMaturities = surface.getMaturities().length;
					double mats[] = surface.getMaturities();
					
					QuotingConvention targetConvention = surface.getQuotingConvention();
					
					for(int t = 0; t<numberOfMaturities; t++) {
						
						double[] currentStrikes = surface.getSmile(mats[t]).getStrikes();
						
						EuropeanOptionSmileMultiAsset newPricer = pricer.getCloneWithModifiedParameters(underlying, mats[t],currentStrikes);
						
						try {
							
							Map<Double, Double> currentModelPrices = newPricer.getValue(newModelFourier);
							
							for(int i = 0; i<currentStrikes.length;i++) {
								
								if(targetConvention.equals(QuotingConvention.VOLATILITYLOGNORMAL)) {
									
									//we convert prices into lognormal volatilities
									double forward = surface.getEquityForwardCurve().getDiscountFactor(mats[t]);
									double optionMaturity =mats[t];
									double optionStrike = currentStrikes[i];
									double payoffUnit = surface.getDiscountCurve().getDiscountFactor(mats[t]);
									double optionValue = currentModelPrices.get(currentStrikes[i]);
									double implVol = Math.abs(net.finmath.functions.AnalyticFormulas.blackScholesOptionImpliedVolatility(forward, optionMaturity, optionStrike, payoffUnit, optionValue));
									if(implVol < 2.0) {
										vals.add(implVol);
									} else {
										vals.add( surface.getSurface().get(optionMaturity).getSmile().get(optionStrike).getValue() - 1E-4 );
									}
									
								
								} else if(targetConvention.equals(QuotingConvention.VOLATILITYNORMAL)) {
									
									//we convert prices into normal volatilities
									double forward = surface.getEquityForwardCurve().getDiscountFactor(mats[t]);
									double optionMaturity =mats[t];
									double optionStrike = currentStrikes[i];
									double payoffUnit = surface.getDiscountCurve().getDiscountFactor(mats[t]);
									double optionValue = currentModelPrices.get(currentStrikes[i]);
									double implVol = Math.abs(net.finmath.functions.AnalyticFormulas.bachelierOptionImpliedVolatility(forward, optionMaturity, optionStrike, payoffUnit, optionValue));
									if(implVol < 2.0) {
										vals.add(implVol);
									} else {
										vals.add( surface.getSurface().get(optionMaturity).getSmile().get(optionStrike).getValue() - 1E-4 );
									}
									
								} else {
									
									//just output the prices
									vals.add(currentModelPrices.get(currentStrikes[i]));
									
								}						
								
							}		
							
						} catch (CalculationException e) {
							e.printStackTrace();
						}	
					}
					
				}
				
				for(int i = 0; i<values.length; i++) {values[i] = vals.get(i);}
				
				long tend = System.currentTimeMillis();
				System.out.println("Iteration required: " + (tend-tstart)/1000.0 + "seconds");
				
			}
			
		};
		
		OptimizerInterface optimizer = optimizerFactory.getOptimizer(
				objectiveFunction,
				initialParameters,
				lowerBound,
				upperBound,
				parameterStep,
				formatTargetValuesForOptimizer() /* targetValues */);
		
		
		optimizer.run();
		
		ArrayList<String> calibrationOutput = outputCalibrationResult(optimizer.getBestFitParameters()); 
		
		MultivariateCalibrableProcessInterface calibratedModel = samplePair.getSecond().getCloneForModifiedParameters(optimizer.getBestFitParameters());
		
		return new OptimizationResult(calibratedModel, optimizer.getBestFitParameters(),
				optimizer.getIterations(), optimizer.getRootMeanSquaredError(), calibrationOutput);
		
	}


	/**
	 * This is a service method that takes care of putting all the target values in a single array.
	 * @return
	 */
	private double[] formatTargetValuesForOptimizer() {
	
		final ArrayList<Double> vals = new ArrayList<>();
		
		for(String key : surfaces.keySet()) {
			OptionSurfaceData surface = surfaces.get(key);
						
			final int numberOfMaturities = surface.getMaturities().length;
			final double[] mats = surface.getMaturities();
			
			for(int t = 0; t<numberOfMaturities; t++) {
				final double mat = mats[t];
				final double[] myStrikes = surface.getSurface().get(mat).getStrikes();

				final OptionSmileData smileOfInterest = surface.getSurface().get(mat);

				for(int k = 0; k < myStrikes.length; k++) {
					vals.add(smileOfInterest.getSmile().get(myStrikes[k]).getValue());
				}
			}
					
		}
		final Double[] targetVals = new Double[vals.size()];
		return ArrayUtils.toPrimitive(vals.toArray(targetVals));		
	}

	/**
	 * When the calibration is over this method is called to produce a table
	 * @param parameters
	 */
	private ArrayList<String> outputCalibrationResult(double[] parameters){
		
		ArrayList<String> calibrationOutput = new ArrayList<String>();
		
		//We change the parameters of the model
		final MultivariateCalibrableProcessInterface newModel = model.getCloneForModifiedParameters(parameters);
		final MultivariateProcessCharacteristicFunctionInterface newModelFourier = newModel.getCharacteristiFunction();

		Iterator<Entry<String, OptionSurfaceData>> iterator = surfaces.entrySet().iterator();
		
		//loop over different surfaces
		while (iterator.hasNext()) {
			final Map.Entry<String, OptionSurfaceData> pair = (Entry<String, OptionSurfaceData>) iterator.next();
			final String underlying = pair.getKey();
			final OptionSurfaceData surface = pair.getValue();
			
			final int numberOfMaturities = surface.getMaturities().length;
			final double mats[] = surface.getMaturities();
			
			final QuotingConvention targetConvention = surface.getQuotingConvention();
			
			
			calibrationOutput.add("Calibration results for " + underlying);
			calibrationOutput.add("Strike"+ "\t" + "Maturity"+ "\t" + "Market Value" + "\t" + "Model Value" + "\t" + "Squared Error");
			
			for(int t = 0; t<numberOfMaturities; t++) {
				final double T = mats[t];
				final OptionSmileData currentSmile = surface.getSmile(mats[t]);
				final double[] currentStrikes = currentSmile.getStrikes();
				
				EuropeanOptionSmileMultiAsset newPricer = pricer.getCloneWithModifiedParameters(underlying, mats[t],currentStrikes);
				
				try {
					
					Map<Double, Double> currentModelPrices = newPricer.getValue(newModelFourier);
					
					for(int i = 0; i<currentStrikes.length;i++) {
						
						final double K = currentStrikes[i];
						final double targetValue = currentSmile.getOption(currentStrikes[i]).getValue();
						final double value;
						
						if(targetConvention.equals(QuotingConvention.VOLATILITYLOGNORMAL)) {
							
							//we convert prices into lognormal volatilities
							final double forward = surface.getEquityForwardCurve().getDiscountFactor(mats[t]);
							final double optionMaturity =mats[t];
							final double optionStrike = currentStrikes[i];
							final double payoffUnit = surface.getDiscountCurve().getDiscountFactor(mats[t]);
							final double optionValue = currentModelPrices.get(currentStrikes[i]);
							double implVol = Math.abs(net.finmath.functions.AnalyticFormulas.blackScholesOptionImpliedVolatility(forward, optionMaturity, optionStrike, payoffUnit, optionValue));
							if(implVol < 2.0) {
								value = implVol;
							} else {
								value = surface.getSurface().get(optionMaturity).getSmile().get(optionStrike).getValue() - 1E-4;
							}
							
						} else if(targetConvention.equals(QuotingConvention.VOLATILITYNORMAL)) {
							
							//we convert prices into normal volatilities
							final double forward = surface.getEquityForwardCurve().getDiscountFactor(mats[t]);
							final double optionMaturity =mats[t];
							final double optionStrike = currentStrikes[i];
							final double payoffUnit = surface.getDiscountCurve().getDiscountFactor(mats[t]);
							final double optionValue = currentModelPrices.get(currentStrikes[i]);
							double implVol = Math.abs(net.finmath.functions.AnalyticFormulas.bachelierOptionImpliedVolatility(forward, optionMaturity, optionStrike, payoffUnit, optionValue));
							if(implVol < 2.0) {
								value = implVol;
							} else {
								value =  surface.getSurface().get(optionMaturity).getSmile().get(optionStrike).getValue() - 1E-4;
							}
							
						} else {
							
							//just output the prices
							value = currentModelPrices.get(currentStrikes[i]);
							
						}	
						
						calibrationOutput.add(K+ "\t" + T + "\t" + targetValue + "\t" + value+ "\t" + Math.pow(targetValue-value,2));
						
					}		
					
				} catch (CalculationException e) {
					e.printStackTrace();
				}	
				
			}
			
		}
		
		return calibrationOutput;
	}

	
	/**
	 * Helper class for calibration results.
	 * @author Alessandro Gnoatto
	 *
	 */
	public class OptimizationResult{
		private final MultivariateCalibrableProcessInterface model; //the calibrated model
		private final double[] bestFitParameters;
		private final int iterations;
		private final double rootMeanSquaredError;
		private final ArrayList<String> calibrationOutput;
		
		public OptimizationResult(MultivariateCalibrableProcessInterface model, double[] bestFitParameters,
				int iterations, double rootMeanSquaredError, ArrayList<String> calibrationOutput) {
			super();
			this.model = model;
			this.bestFitParameters = bestFitParameters;
			this.iterations = iterations;
			this.rootMeanSquaredError = rootMeanSquaredError;
			this.calibrationOutput = calibrationOutput;
		}

		public MultivariateCalibrableProcessInterface getModel() {
			return model;
		}
		
		public double[] getBestFitParameters() {
			return bestFitParameters;
		}

		public int getIterations() {
			return iterations;
		}

		public double getRootMeanSquaredError() {
			return rootMeanSquaredError;
		}
		
		public ArrayList<String> getCalibrationOutput(){
			return this.calibrationOutput;
		}
		
	}
	
	private Pair< double[], MultivariateCalibrableProcessInterface > generateSamplePair() throws IllegalArgumentException {
		
		Pair< List<Double>, MultivariateCalibrableProcessInterface > pair = this.model.generateSamplePair();
			
		List<Double> newParams = pair.getFirst();
			
		double[] params = ArrayUtils.toPrimitive( newParams.toArray(new Double[newParams.size()]) );
					
		MultivariateCalibrableProcessInterface newModel = pair.getSecond();
			
		return Pair.< double[], MultivariateCalibrableProcessInterface >create(params, newModel);
		
	}

}
