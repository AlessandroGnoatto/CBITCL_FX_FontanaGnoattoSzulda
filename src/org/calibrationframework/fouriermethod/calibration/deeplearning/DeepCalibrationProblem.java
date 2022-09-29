package org.calibrationframework.fouriermethod.calibration.deeplearning;

import java.io.*;
import java.util.*;

import org.apache.commons.lang3.ArrayUtils;
import org.calibrationframework.fouriermethod.calibration.models.MultivariateCalibrableProcessInterface;
import org.calibrationframework.marketdata.model.volatilities.*;
import org.calibrationframework.optimizer.*;

import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.deeplearning4j.util.ModelSerializer;
import org.nd4j.common.primitives.Pair;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.api.preprocessor.DataNormalization;
import org.nd4j.linalg.factory.Nd4j;

/**
 * This class enables to calibrate the underlying pricing model to the market surfaces via optimization.
 * 
 * It differs from standard calibration procedures in the sense that the underlying pricing function(s) is/are replaced by the previously trained and tested neural-network approximation(s),
 * in the calibration criterion and then within the iterative optimization algorithm, thus making the procedure of deterministic nature and then faster than others.
 * 
 * To proceed, its main method runCalibration() has to be launched after runDeepLearningApproximation() of DeepApproximation, 
 * thus allowing for the restoration of the trained neural network(s) and the associated data normalizer(s) that lie inside the corresponding .zip file(s).
 * 
 * At the very end of the optimization, the outcome will be displayed similarly to the other calibration classes, comparing the markets prices/impl vols to the neural network ones for the best fit parameter set.
 * 
 * @author Szulda Guillaume
 *
 */
public class DeepCalibrationProblem {
	
	private final MultivariateCalibrableProcessInterface model;
	private final LinkedHashMap<String, WeightedOptionSurfaceData> surfaces;
	private final OptimizerFactoryInterface optimizerFactory; 
	
	private final Pair< double[], MultivariateCalibrableProcessInterface > samplePair;
	
	private final DeepApproximation.Data dataType;
	
	private MultiLayerNetwork neuralNetwork;
	private LinkedHashMap<String, MultiLayerNetwork> neuralNetworks;
	
	private DataNormalization dataNormalizer;
	private LinkedHashMap<String, DataNormalization> dataNormalizers;

	//Optimizer parameters
	private final double[] initialParameters;
	private final double[] lowerBound;
	private final double[] upperBound;
	private final double[] parameterStep;
	
	/**
	 * Constructor for the deep calibration problem.
	 * 
	 * The input parameter dataType gives information about the data that has been previously generated and the associated trained and tested neural network(s).
	 * 
	 * Here, the initial parameters for the iterative optimization algorithm are drawn randomly and the associated steps are all set to 0.01.
	 * 
	 * @param model
	 * @param surfaces
	 * @param optimizerFactory
	 * @param dataType
	 * @throws IllegalArgumentException
	 * @throws IOException
	 */
	public DeepCalibrationProblem(MultivariateCalibrableProcessInterface model, 
			LinkedHashMap<String, WeightedOptionSurfaceData> surfaces, OptimizerFactoryInterface optimizerFactory,
			DeepApproximation.Data dataType) throws IllegalArgumentException, IOException {
		
		this.model = model;
		this.surfaces = surfaces;
		this.optimizerFactory = optimizerFactory;
		
		this.samplePair = generateSamplePair();
		
		this.initialParameters = this.samplePair.getFirst();
			
		this.lowerBound = this.samplePair.getSecond().getParameterLowerBounds();
		this.upperBound = this.samplePair.getSecond().getParameterUpperBounds();
			
		this.parameterStep = new double[this.initialParameters.length];
		for(int i = 0; i < this.initialParameters.length; i++) 
			this.parameterStep[i] = 0.01;
		
		this.dataType = dataType;
		
		this.neuralNetwork = null;
		this.neuralNetworks = null;
		
		this.dataNormalizer = null;
		this.dataNormalizers = null;
		
		if(this.dataType==DeepApproximation.Data.ALLSURFACES) {
			
			//Restoration of the neural network model plus its data normalizer
			File modelFile = new File("modelAllSurfaces.zip");
			
			if(modelFile.exists()) {
				this.neuralNetwork = ModelSerializer.restoreMultiLayerNetwork(modelFile);
				this.dataNormalizer = ModelSerializer.restoreNormalizerFromFile(modelFile);
			} else {
				throw new IllegalArgumentException("runDeepLearningApproximation() of DeepApproximation must be carried out before calibration");
			}
			
		} else if(this.dataType==DeepApproximation.Data.EACHSURFACE) {
			
			this.neuralNetworks = new LinkedHashMap<String, MultiLayerNetwork>(this.surfaces.size());
			this.dataNormalizers = new LinkedHashMap<String, DataNormalization>(this.surfaces.size());
			
			Iterator<Map.Entry<String, WeightedOptionSurfaceData>> surfaceIterator = this.surfaces.entrySet().iterator();
			
			//Loop over the surfaces
			while (surfaceIterator.hasNext()) {
				
				Map.Entry<String, WeightedOptionSurfaceData> pair = 
						(Map.Entry<String, WeightedOptionSurfaceData>) surfaceIterator.next();
				
				String surfaceName = pair.getKey();
				
				//Restoration of the neural network model plus its data normalizer
				File modelFile = new File("modelSurface" + surfaceName + ".zip");
				
				if(modelFile.exists()) {
					this.neuralNetworks.put(surfaceName, ModelSerializer.restoreMultiLayerNetwork(modelFile));
					this.dataNormalizers.put(surfaceName, ModelSerializer.restoreNormalizerFromFile(modelFile));
				} else {
					throw new IllegalArgumentException("runDeepLearningApproximation() of DeepApproximation must be carried out before calibration");
				}
				
			}
		
		} else if(this.dataType==DeepApproximation.Data.PRICES) {
			
			this.neuralNetworks = new LinkedHashMap<String, MultiLayerNetwork>(this.surfaces.size());
			this.dataNormalizers = new LinkedHashMap<String, DataNormalization>(this.surfaces.size());
			
			Iterator<Map.Entry<String, WeightedOptionSurfaceData>> surfaceIterator = this.surfaces.entrySet().iterator();
			
			//Loop over the surfaces
			while (surfaceIterator.hasNext()) {
				
				Map.Entry<String, WeightedOptionSurfaceData> pair = 
						(Map.Entry<String, WeightedOptionSurfaceData>) surfaceIterator.next();
				
				String surfaceName = pair.getKey();
				
				//Restoration of the neural network model plus its data normalizer
				File modelFile = new File("modelPrices" + surfaceName + ".zip");
				
				if(modelFile.exists()) {
					this.neuralNetworks.put(surfaceName, ModelSerializer.restoreMultiLayerNetwork(modelFile));
					this.dataNormalizers.put(surfaceName, ModelSerializer.restoreNormalizerFromFile(modelFile));
				} else {
					throw new IllegalArgumentException("runDeepLearningApproximation() of DeepApproximation must be carried out before calibration");
				}
				
			}
			
			
		} else if(this.dataType==DeepApproximation.Data.SMILES) {
			
			this.neuralNetworks = new LinkedHashMap<String, MultiLayerNetwork>(this.surfaces.size());
			this.dataNormalizers = new LinkedHashMap<String, DataNormalization>(this.surfaces.size());
			
			Iterator<Map.Entry<String, WeightedOptionSurfaceData>> surfaceIterator = this.surfaces.entrySet().iterator();
			
			//Loop over the surfaces
			while (surfaceIterator.hasNext()) {
				
				Map.Entry<String, WeightedOptionSurfaceData> pair = 
						(Map.Entry<String, WeightedOptionSurfaceData>) surfaceIterator.next();
				
				String surfaceName = pair.getKey();
				
				//Restoration of the neural network model plus its data normalizer
				File modelFile = new File("modelSmiles" + surfaceName + ".zip");
				
				if(modelFile.exists()) {
					this.neuralNetworks.put(surfaceName, ModelSerializer.restoreMultiLayerNetwork(modelFile));
					this.dataNormalizers.put(surfaceName, ModelSerializer.restoreNormalizerFromFile(modelFile));
				} else {
					throw new IllegalArgumentException("runDeepLearningApproximation() of DeepApproximation must be carried out before calibration");
				}
				
			}
			
		}
		
	}
	
	/**
	 * Constructor for the deep calibration problem.
	 * 
	 * The input parameter dataType gives information about the data that has been previously generated along with the associated trained and tested neural network(s).
	 * 
	 * Here, the initial parameters for the iterative optimization algorithm and the associated steps are chosen by the user.
	 * 
	 * @param model
	 * @param surfaces
	 * @param optimizerFactory
	 * @param dataType
	 * @param initialParameters
	 * @param parameterStep
	 * @throws IllegalArgumentException
	 * @throws IOException
	 */
	public DeepCalibrationProblem(MultivariateCalibrableProcessInterface model, 
			LinkedHashMap<String, WeightedOptionSurfaceData> surfaces, OptimizerFactoryInterface optimizerFactory, 
			DeepApproximation.Data dataType, double[] initialParameters, double[] parameterStep) 
					throws IllegalArgumentException, IOException {
		
		this.model = model;
		this.surfaces = surfaces;
		this.optimizerFactory = optimizerFactory;
					
		this.samplePair = generateSamplePair();
				
		this.initialParameters = initialParameters;
					
		this.lowerBound = this.samplePair.getSecond().getParameterLowerBounds();
		this.upperBound = this.samplePair.getSecond().getParameterUpperBounds();
					
		this.parameterStep = parameterStep;
		
		this.dataType = dataType;
		
		this.neuralNetwork = null;
		this.neuralNetworks = null;
		
		this.dataNormalizer = null;
		this.dataNormalizers = null;
		
		if(this.dataType==DeepApproximation.Data.ALLSURFACES) {
			
			//Restoration of the neural network model plus its data normalizer
			File modelFile = new File("modelAllSurfaces.zip");
			
			if(modelFile.exists()) {
				this.neuralNetwork = ModelSerializer.restoreMultiLayerNetwork(modelFile);
				this.dataNormalizer = ModelSerializer.restoreNormalizerFromFile(modelFile);
			} else {
				throw new IllegalArgumentException("runDeepLearningApproximation() of DeepApproximation must be carried out before calibration");
			}
			
		} else if(this.dataType==DeepApproximation.Data.EACHSURFACE) {
			
			this.neuralNetworks = new LinkedHashMap<String, MultiLayerNetwork>(this.surfaces.size());
			this.dataNormalizers = new LinkedHashMap<String, DataNormalization>(this.surfaces.size());
			
			Iterator<Map.Entry<String, WeightedOptionSurfaceData>> surfaceIterator = this.surfaces.entrySet().iterator();
			
			//Loop over the surfaces
			while (surfaceIterator.hasNext()) {
				
				Map.Entry<String, WeightedOptionSurfaceData> pair = 
						(Map.Entry<String, WeightedOptionSurfaceData>) surfaceIterator.next();
				
				String surfaceName = pair.getKey();
				
				//Restoration of the neural network model plus its data normalizer
				File modelFile = new File("modelSurface" + surfaceName + ".zip");
				
				if(modelFile.exists()) {
					this.neuralNetworks.put(surfaceName, ModelSerializer.restoreMultiLayerNetwork(modelFile));
					this.dataNormalizers.put(surfaceName, ModelSerializer.restoreNormalizerFromFile(modelFile));
				} else {
					throw new IllegalArgumentException("runDeepLearningApproximation() of DeepApproximation must be carried out before calibration");
				}
				
			}
		
		} else if(this.dataType==DeepApproximation.Data.PRICES) {
			
			this.neuralNetworks = new LinkedHashMap<String, MultiLayerNetwork>(this.surfaces.size());
			this.dataNormalizers = new LinkedHashMap<String, DataNormalization>(this.surfaces.size());
			
			Iterator<Map.Entry<String, WeightedOptionSurfaceData>> surfaceIterator = this.surfaces.entrySet().iterator();
			
			//Loop over the surfaces
			while (surfaceIterator.hasNext()) {
				
				Map.Entry<String, WeightedOptionSurfaceData> pair = 
						(Map.Entry<String, WeightedOptionSurfaceData>) surfaceIterator.next();
				
				String surfaceName = pair.getKey();
				
				//Restoration of the neural network model plus its data normalizer
				File modelFile = new File("modelPrices" + surfaceName + ".zip");
				
				if(modelFile.exists()) {
					this.neuralNetworks.put(surfaceName, ModelSerializer.restoreMultiLayerNetwork(modelFile));
					this.dataNormalizers.put(surfaceName, ModelSerializer.restoreNormalizerFromFile(modelFile));
				} else {
					throw new IllegalArgumentException("runDeepLearningApproximation() of DeepApproximation must be carried out before calibration");
				}
				
			}
			
		} else if(this.dataType==DeepApproximation.Data.SMILES) {
			
			this.neuralNetworks = new LinkedHashMap<String, MultiLayerNetwork>(this.surfaces.size());
			this.dataNormalizers = new LinkedHashMap<String, DataNormalization>(this.surfaces.size());
			
			Iterator<Map.Entry<String, WeightedOptionSurfaceData>> surfaceIterator = this.surfaces.entrySet().iterator();
			
			//Loop over the surfaces
			while (surfaceIterator.hasNext()) {
				
				Map.Entry<String, WeightedOptionSurfaceData> pair = 
						(Map.Entry<String, WeightedOptionSurfaceData>) surfaceIterator.next();
				
				String surfaceName = pair.getKey();
				
				//Restoration of the neural network model plus its data normalizer
				File modelFile = new File("modelSmiles" + surfaceName + ".zip");
				
				if(modelFile.exists()) {
					this.neuralNetworks.put(surfaceName, ModelSerializer.restoreMultiLayerNetwork(modelFile));
					this.dataNormalizers.put(surfaceName, ModelSerializer.restoreNormalizerFromFile(modelFile));
				} else {
					throw new IllegalArgumentException("runDeepLearningApproximation() of DeepApproximation must be carried out before calibration");
				}
				
			}
			
		}
		
	}
	
	public OptimizationResult runCalibration() throws SolverException {
		
		OptimizerInterface.ObjectiveFunction objectiveFunction = null;
		
		if(this.dataType==DeepApproximation.Data.ALLSURFACES) {
			
			objectiveFunction = new OptimizerInterface.ObjectiveFunction() {

				@Override
				public void setValues(double[] parameters, double[] values) throws SolverException {
					
					//Creation of the input for the neural network
					INDArray input = Nd4j.create(parameters, new int[] {1, parameters.length});
					
					//Normalization of the input by means of the normalizer of the training/testing data set
					dataNormalizer.transform(input);
					
					//Determination of the corresponding output through the neural network 
					INDArray output = neuralNetwork.output(input);
					
					//Normalization of the output by the normalizer
					//dataNormalizer.transformLabel(output);
					//dataNormalizer.revertLabels(output);
					
					double[] myValues = output.toDoubleVector();
					
					for(int i = 0; i < values.length; i++) {
						values[i] = myValues[i];
					}
					
				}
				
			};
			
		} else if(this.dataType==DeepApproximation.Data.EACHSURFACE) {
			
			objectiveFunction = new OptimizerInterface.ObjectiveFunction() {

				@Override
				public void setValues(double[] parameters, double[] values) throws SolverException {
					
					List<Double> vals = new ArrayList<Double>();
					
					//Iterator pointing at the surfaces
					Iterator<Map.Entry<String, WeightedOptionSurfaceData>> surfaceIterator = surfaces.entrySet().iterator();
					
					//Loop over the surfaces
					while (surfaceIterator.hasNext()) {
						
						Map.Entry<String, WeightedOptionSurfaceData> pair = 
								(Map.Entry<String, WeightedOptionSurfaceData>) surfaceIterator.next();
						
						String surfaceName = pair.getKey();
						
						//Creation of the input for this surface
						INDArray currentInput = Nd4j.create(parameters, new int[] {1, parameters.length});
						
						//Normalization of the input by means of the normalizer of the training/testing data set for the current surface
						DataNormalization currentDataNormalizer = dataNormalizers.get(surfaceName);
						currentDataNormalizer.transform(currentInput);
						
						//Determination of the corresponding output 
						MultiLayerNetwork currentNeuralNetwork = neuralNetworks.get(surfaceName);
						INDArray currentOutput = currentNeuralNetwork.output(currentInput);
						double[] currentOutputArray = currentOutput.toDoubleVector();
						for(double x : currentOutputArray)
							vals.add(x);
						
					}
					
					double[] myValues = ArrayUtils.toPrimitive( vals.toArray(new Double[vals.size()]) );
					
					for(int i = 0; i < values.length; i++) {
						values[i] = myValues[i];
					}
					
				}
				
			};
			
		
		} else if(this.dataType==DeepApproximation.Data.PRICES) {
		
			objectiveFunction = new OptimizerInterface.ObjectiveFunction() {

				@Override
				public void setValues(double[] parameters, double[] values) throws SolverException {
					
					List<Double> vals = new ArrayList<Double>();
					
					//Iterator pointing at the surfaces
					Iterator<Map.Entry<String, WeightedOptionSurfaceData>> surfaceIterator = surfaces.entrySet().iterator();
					
					//Loop over the surfaces
					while (surfaceIterator.hasNext()) {
						
						Map.Entry<String, WeightedOptionSurfaceData> pair = 
								(Map.Entry<String, WeightedOptionSurfaceData>) surfaceIterator.next();
						
						String surfaceName = pair.getKey();
						
						DataNormalization currentDataNormalizer = dataNormalizers.get(surfaceName);
						MultiLayerNetwork currentNeuralNetwork = neuralNetworks.get(surfaceName);
						
						OptionSurfaceData currentSurface = pair.getValue();
						
						double mats[] = currentSurface.getMaturities();
						
						for(int t = 0; t< mats.length; t++) {
							
							double maturity = mats[t];
									
							double[] strikes = currentSurface.getSmile(maturity).getStrikes();
									
							for(int i = 0; i < strikes.length; i++) {
								
								double strike = strikes[i];
								
								List<Double> newParams = new ArrayList<Double>(parameters.length+2);
								for(double x : parameters)
									newParams.add(x);
								newParams.add(maturity);
								newParams.add(strike);
								double[] newParamsArray = ArrayUtils.toPrimitive( newParams.toArray(new Double[newParams.size()]) );
								
								//Creation of the current input
								INDArray currentInput = Nd4j.create(newParamsArray, new int[] {1, newParamsArray.length});
								
								//Normalization of the input by means of the normalizer of the training/testing data set for the current surface
								currentDataNormalizer.transform(currentInput);
								
								//Determination of the corresponding output 
								INDArray currentOutput = currentNeuralNetwork.output(currentInput);
								double[] currentOutputArray = currentOutput.toDoubleVector();
								for(double x : currentOutputArray)
									vals.add(x);
								
							}
							
						}
					
					}
					
					double[] myValues = ArrayUtils.toPrimitive( vals.toArray(new Double[vals.size()]) );
					
					for(int i = 0; i < values.length; i++) {
						values[i] = myValues[i];
					}
					
				}
				
			};
			
		} else if(this.dataType==DeepApproximation.Data.SMILES) {
			
			objectiveFunction = new OptimizerInterface.ObjectiveFunction() {

				@Override
				public void setValues(double[] parameters, double[] values) throws SolverException {
					
					List<Double> vals = new ArrayList<Double>();
					
					//Iterator pointing at the surfaces
					Iterator<Map.Entry<String, WeightedOptionSurfaceData>> surfaceIterator = surfaces.entrySet().iterator();
					
					//Loop over the surfaces
					while (surfaceIterator.hasNext()) {
						
						Map.Entry<String, WeightedOptionSurfaceData> pair = 
								(Map.Entry<String, WeightedOptionSurfaceData>) surfaceIterator.next();
						
						String surfaceName = pair.getKey();
						
						DataNormalization currentDataNormalizer = dataNormalizers.get(surfaceName);
						MultiLayerNetwork currentNeuralNetwork = neuralNetworks.get(surfaceName);
						
						OptionSurfaceData currentSurface = pair.getValue();
						
						double mats[] = currentSurface.getMaturities();
						
						for(int t = 0; t< mats.length; t++) {
							
							double maturity = mats[t];
									
							double[] strikes = currentSurface.getSmile(maturity).getStrikes();
							
							List<Double> newParams = new ArrayList<Double>(parameters.length+1+strikes.length);
							for(double x : parameters)
								newParams.add(x);
							newParams.add(maturity);
							for(double y : strikes)
								newParams.add(y);
							double[] newParamsArray = ArrayUtils.toPrimitive( newParams.toArray(new Double[newParams.size()]) );
								
							//Creation of the current input
							INDArray currentInput = Nd4j.create(newParamsArray, new int[] {1, newParamsArray.length});
								
							//Normalization of the input by means of the normalizer of the training/testing data set for the current surface
							currentDataNormalizer.transform(currentInput);
								
							//Determination of the corresponding output 
							INDArray currentOutput = currentNeuralNetwork.output(currentInput);
							double[] currentOutputArray = currentOutput.toDoubleVector();
							for(double x : currentOutputArray)
								vals.add(x);
								
						}
							
					}
					
					double[] myValues = ArrayUtils.toPrimitive( vals.toArray(new Double[vals.size()]) );
					
					for(int i = 0; i < values.length; i++) 
						values[i] = myValues[i];
					
				}
				
			};
			
		}
		
		OptimizerInterface optimizer = optimizerFactory.getOptimizer(
				objectiveFunction,
				initialParameters,
				lowerBound,
				upperBound,
				parameterStep,
				formatTargetValuesForOptimizer() /* targetValues */);
		
		optimizer.run();
		
		ArrayList<String> comparisonMarketNeuralNetwork = comparisonMarketNeuralNetwork(optimizer.getBestFitParameters()); 
		
		MultivariateCalibrableProcessInterface calibratedModel = samplePair.getSecond().getCloneForModifiedParameters(optimizer.getBestFitParameters());
		
		return new OptimizationResult(calibratedModel, optimizer.getBestFitParameters(), optimizer.getIterations(), 
				optimizer.getRootMeanSquaredError(), comparisonMarketNeuralNetwork);
		
	}

	/**
	 * This is a service method that takes care of putting all the target values in a single array.
	 * @return
	 */
	private double[] formatTargetValuesForOptimizer() {
	
		final ArrayList<Double> vals = new ArrayList<>();
		
		for(String key : this.surfaces.keySet()) {
			OptionSurfaceData surface = this.surfaces.get(key);
						
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
	 * When the calibration is over this method is called to produce a comparison table,
	 * between market and neural network prices/impl vols.
	 * 
	 * @param parameters
	 */
	private ArrayList<String> comparisonMarketNeuralNetwork(double[] parameters) {
		
		ArrayList<String> comparisonMarketNeuralNetwork = new ArrayList<String>();
		
		comparisonMarketNeuralNetwork.add("Calibration results: comparison table between market and neural network prices/impl vols");
		
		if(this.dataType==DeepApproximation.Data.ALLSURFACES) {
			
			//Creation of the input for the neural network
			INDArray input = Nd4j.create(parameters, new int[] {1, parameters.length});
			
			//Normalization of the input by means of the normalizer of the training/testing data set
			this.dataNormalizer.transform(input);
			
			//Determination of the corresponding output through the neural network 
			INDArray output = this.neuralNetwork.output(input);
			
			//Normalization of the output by the normalizer
			//this.dataNormalizer.transformLabel(output);
			//this.dataNormalizer.revertLabels(output);
			
			//Transformation of the output format, from INDArray to standard Array
			double[] outputArray = output.toDoubleVector();
			
			int surfaceCapacity = 0;

			//Iterator pointing at the surfaces
			Iterator<Map.Entry<String, WeightedOptionSurfaceData>> iterator = this.surfaces.entrySet().iterator();
			
			//loop over different surfaces
			while (iterator.hasNext()) {
				
				Map.Entry<String, WeightedOptionSurfaceData> pair = 
						(Map.Entry<String, WeightedOptionSurfaceData>) iterator.next();
				
				String underlying = pair.getKey();
				OptionSurfaceData surface = pair.getValue();
				
				int numberOfMaturities = surface.getMaturities().length;
				double mats[] = surface.getMaturities();
				
				comparisonMarketNeuralNetwork.add("Calibration results for " + underlying);
				comparisonMarketNeuralNetwork.add("Strike"+ "\t" + "Maturity"+ "\t" + "Market Value" + "\t" + "Neural Network Value" + "\t" + "Squared Error");
				
				for(int t = 0; t<numberOfMaturities; t++) {
					
					double T = mats[t];
					
					OptionSmileData currentSmile = surface.getSmile(T);
					
					double[] currentStrikes = currentSmile.getStrikes();
						
					for(int i = 0; i<currentStrikes.length;i++) {
						
						double K = currentStrikes[i];
						
						double targetValue = currentSmile.getOption(K).getValue();
						
						double value = Math.abs(outputArray[surfaceCapacity + i]);
						
						comparisonMarketNeuralNetwork.add(K+ "\t" + T + "\t" + targetValue + "\t" + value+ "\t" + Math.pow(targetValue-value,2));
				
					}
					
					surfaceCapacity = surfaceCapacity + currentStrikes.length;
						
				}
				
			}
			
		} else if(this.dataType==DeepApproximation.Data.EACHSURFACE) {

			//Iterator pointing at the surfaces
			Iterator<Map.Entry<String, WeightedOptionSurfaceData>> iterator = this.surfaces.entrySet().iterator();
			
			//loop over different surfaces
			while (iterator.hasNext()) {
				
				//Creation of the current input
				INDArray currentInput = Nd4j.create(parameters, new int[] {1, parameters.length});
				
				Map.Entry<String, WeightedOptionSurfaceData> pair = 
						(Map.Entry<String, WeightedOptionSurfaceData>) iterator.next();
				
				String underlying = pair.getKey();
				
				//Normalization of the input by means of the normalizer of the training/testing data set
				DataNormalization currentDataNormalizer = this.dataNormalizers.get(underlying);
				currentDataNormalizer.transform(currentInput);
				
				//Determination of the corresponding output 
				MultiLayerNetwork currentNeuralNetwork = this.neuralNetworks.get(underlying);
				INDArray currentOutput = currentNeuralNetwork.output(currentInput);
				double[] currentOutputArray = currentOutput.toDoubleVector();
				
				OptionSurfaceData surface = pair.getValue();
				
				int numberOfMaturities = surface.getMaturities().length;
				double mats[] = surface.getMaturities();
				
				comparisonMarketNeuralNetwork.add("Calibration results for " + underlying);
				comparisonMarketNeuralNetwork.add("Strike"+ "\t" + "Maturity"+ "\t" + "Market Value" + "\t" + "Neural Network Value" + "\t" + "Squared Error");
				
				int surfaceCapacity = 0;
				
				for(int t = 0; t<numberOfMaturities; t++) {
					
					double T = mats[t];
					
					OptionSmileData currentSmile = surface.getSmile(T);
					
					double[] currentStrikes = currentSmile.getStrikes();
						
					for(int i = 0; i<currentStrikes.length;i++) {
						
						double K = currentStrikes[i];
						
						double targetValue = currentSmile.getOption(K).getValue();
						
						double value = Math.abs(currentOutputArray[surfaceCapacity + i]);
						
						comparisonMarketNeuralNetwork.add(K+ "\t" + T + "\t" + targetValue + "\t" + value+ "\t" + Math.pow(targetValue-value,2));
				
					}
					
					surfaceCapacity = surfaceCapacity + currentStrikes.length;
						
				}
				
			}
			
			
		} else if(this.dataType==DeepApproximation.Data.PRICES) {
			
			//Iterator pointing at the surfaces
			Iterator<Map.Entry<String, WeightedOptionSurfaceData>> iterator = this.surfaces.entrySet().iterator();
			
			//loop over different surfaces
			while (iterator.hasNext()) {
				
				Map.Entry<String, WeightedOptionSurfaceData> pair = 
						(Map.Entry<String, WeightedOptionSurfaceData>) iterator.next();
	
				String underlying = pair.getKey();
				
				DataNormalization currentDataNormalizer = this.dataNormalizers.get(underlying);
				MultiLayerNetwork currentNeuralNetwork = this.neuralNetworks.get(underlying);
				
				OptionSurfaceData surface = pair.getValue();
				
				int numberOfMaturities = surface.getMaturities().length;
				double mats[] = surface.getMaturities();
				
				comparisonMarketNeuralNetwork.add("Calibration results for " + underlying);
				comparisonMarketNeuralNetwork.add("Strike"+ "\t" + "Maturity"+ "\t" + "Market Value" + "\t" + "Neural Network Value" + "\t" + "Squared Error");
				
				for(int t = 0; t<numberOfMaturities; t++) {
					
					double T = mats[t];
					
					OptionSmileData currentSmile = surface.getSmile(T);
					
					double[] currentStrikes = currentSmile.getStrikes();
						
					for(int i = 0; i<currentStrikes.length;i++) {
						
						double K = currentStrikes[i];
						
						double targetValue = currentSmile.getOption(K).getValue();
						
						//Creation of the current input
						List<Double> newParams = new ArrayList<Double>(parameters.length+2);
						for(double x : parameters)
							newParams.add(x);
						newParams.add(T);
						newParams.add(K);
						double[] newParamsArray = ArrayUtils.toPrimitive( newParams.toArray(new Double[newParams.size()]) );
						INDArray currentInput = Nd4j.create(newParamsArray, new int[] {1, newParamsArray.length});
					
						//Normalization of the input by means of the normalizer of the training/testing data set
						currentDataNormalizer.transform(currentInput);
						
						//Determination of the corresponding output
						INDArray currentOutput = currentNeuralNetwork.output(currentInput);
						double[] currentOutputArray = currentOutput.toDoubleVector();
						double value = Math.abs(currentOutputArray[0]);
						
						comparisonMarketNeuralNetwork.add(K+ "\t" + T + "\t" + targetValue + "\t" + value+ "\t" + Math.pow(targetValue-value,2));
				
					}
						
				}
				
			}
			
		} else if(this.dataType==DeepApproximation.Data.SMILES) {
			
			//Iterator pointing at the surfaces
			Iterator<Map.Entry<String, WeightedOptionSurfaceData>> iterator = this.surfaces.entrySet().iterator();
			
			//loop over different surfaces
			while (iterator.hasNext()) {
				
				Map.Entry<String, WeightedOptionSurfaceData> pair = 
						(Map.Entry<String, WeightedOptionSurfaceData>) iterator.next();
	
				String underlying = pair.getKey();
				
				DataNormalization currentDataNormalizer = this.dataNormalizers.get(underlying);
				MultiLayerNetwork currentNeuralNetwork = this.neuralNetworks.get(underlying);
				
				OptionSurfaceData surface = pair.getValue();
				
				int numberOfMaturities = surface.getMaturities().length;
				double mats[] = surface.getMaturities();
				
				comparisonMarketNeuralNetwork.add("Calibration results for " + underlying);
				comparisonMarketNeuralNetwork.add("Strike"+ "\t" + "Maturity"+ "\t" + "Market Value" + "\t" + "Neural Network Value" + "\t" + "Squared Error");
				
				for(int t = 0; t<numberOfMaturities; t++) {
					
					double T = mats[t];
					
					OptionSmileData currentSmile = surface.getSmile(T);
					
					double[] strikes = currentSmile.getStrikes();
					
					//Creation of the current input
					List<Double> newParams = new ArrayList<Double>(parameters.length+1+strikes.length);
					for(double x : parameters)
						newParams.add(x);
					newParams.add(T);
					for(double y : strikes)
						newParams.add(y);
					double[] newParamsArray = ArrayUtils.toPrimitive( newParams.toArray(new Double[newParams.size()]) );
					INDArray currentInput = Nd4j.create(newParamsArray, new int[] {1, newParamsArray.length});
						
					//Normalization of the input by means of the normalizer of the training/testing data set for the current surface
					currentDataNormalizer.transform(currentInput);
						
					//Determination of the corresponding output 
					INDArray currentOutput = currentNeuralNetwork.output(currentInput);
					double[] currentOutputArray = currentOutput.toDoubleVector();
						
					for(int i = 0; i < strikes.length; i++) {
						
						double K = strikes[i];
						
						double targetValue = currentSmile.getOption(K).getValue();
						
						double value = Math.abs(currentOutputArray[i]);
						
						comparisonMarketNeuralNetwork.add(K+ "\t" + T + "\t" + targetValue + "\t" + value+ "\t" + Math.pow(targetValue-value,2));
				
					}
						
				}
				
			}
			
		}
		
		return comparisonMarketNeuralNetwork;
		
	}
	
	/**
	 * Helper class for calibration results.
	 * @author Alessandro Gnoatto
	 */
	public class OptimizationResult {
		
		private final MultivariateCalibrableProcessInterface model; //The calibrated model
		private final double[] bestFitParameters; //The best-fit parameters returned by the optimizer at the end of the calibration routine
		private final int iterations;
		private final double rootMeanSquaredError;
		private final ArrayList<String> comparisonMarketNeuralNetwork;
		
		public OptimizationResult(MultivariateCalibrableProcessInterface model, double[] bestFitParameters, int iterations, 
				double rootMeanSquaredError, ArrayList<String> comparisonMarketNeuralNetwork) {
			this.model = model;
			this.bestFitParameters = bestFitParameters;
			this.iterations = iterations;
			this.rootMeanSquaredError = rootMeanSquaredError;
			this.comparisonMarketNeuralNetwork = comparisonMarketNeuralNetwork;
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
		
		public ArrayList<String> getComparisonMarketNeuralNetwork() {
			return comparisonMarketNeuralNetwork;
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
