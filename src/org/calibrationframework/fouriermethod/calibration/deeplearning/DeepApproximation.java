package org.calibrationframework.fouriermethod.calibration.deeplearning;

import java.io.*;
import java.util.*;

import org.calibrationframework.fouriermethod.calibration.models.*;
import org.calibrationframework.marketdata.model.volatilities.*;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.deeplearning4j.nn.weights.WeightInit;
import org.deeplearning4j.nn.api.OptimizationAlgorithm;
import org.deeplearning4j.nn.conf.*;
import org.deeplearning4j.nn.conf.layers.*;
import org.deeplearning4j.optimize.listeners.*;
import org.deeplearning4j.util.ModelSerializer;

import org.nd4j.evaluation.regression.RegressionEvaluation;
import org.nd4j.linalg.dataset.*;
import org.nd4j.linalg.dataset.api.preprocessor.*;
import org.nd4j.linalg.activations.Activation;
import org.nd4j.linalg.dataset.api.iterator.*;
import org.nd4j.linalg.learning.config.Adam;
import org.nd4j.linalg.lossfunctions.LossFunctions.LossFunction;

/**
 * This class deals with the training and testing of the neural network(s) approximating pricing function(s).
 * 
 * To proceed:
 * 
 * 			- One of the main methods of the class DataGenerator has to be launched first, which will generate one or several data set(s);
 * 
 * 			- Then, the main method of the present class runDeepLearningApproximation() will train and test the neural network(s).
 * 
 * At the end of this procedure, for each neural network, the resulting metrics of the evaluation/testing phase will be displayed.
 * 
 * Finally, the progress and status of each training phase can be monitored by means of a UI server,
 * 
 * To access it, open your browser and go to http://localhost:9000/train/overview
 * 
 * @author Szulda Guillaume
 * 
 */
public class DeepApproximation {
	
	private final MultivariateCalibrableProcessInterface model;
	private final LinkedHashMap<String, WeightedOptionSurfaceData> surfaces;
	
	private final DeepApproximation.Data dataType;
	
	private final DataNormalization dataNormalizer;
	private LinkedHashMap<String, DataNormalization> dataNormalizers;
	
	private int[] nIns; //Input size(s)
	private int[] nOuts; //Output size(s)
	
	/**
	 * Constructor of the class dealing with the training and testing of the neural network(s):
	 * 
	 * 		- The input parameter dataType provides information about the data set(s) previously generated via the class DataGenerator:
	 * 				
	 * 				* ALLSURFACES for generateDataAllSurfaces(), then one data set;
	 * 
	 * 				* EACHSURFACE for generateDataEachSurface(), then as many data sets as there are surfaces;
	 * 
	 * 				* PRICES for generateDataPrices(), then as many data sets as there are surfaces;
	 * 
	 * 				* SMILES for generateDataSmiles(), then as many data sets as there are surfaces;
	 * 
	 * 		- The input model allows to retrieve the number of parameters to calibrate,
	 * 			which represent a part of the input(s) of the neural network(s);
	 * 
	 * 		- The input surfaces allow to retrieve the output size(s) of the neural network(s);
	 * 
	 * 		- The input data normalizer takes care of the normalization of the features (and labels) of the data set(s) before both training and testing.
	 * 
	 * @param model
	 * @param surfaces
	 * @param dataNormalizer
	 * @param dataType
	 * @throws IllegalArgumentException
	 */
	public DeepApproximation(MultivariateCalibrableProcessInterface model, 
			LinkedHashMap<String, WeightedOptionSurfaceData> surfaces, DataNormalization dataNormalizer,
			DeepApproximation.Data dataType) throws IllegalArgumentException {
		
		this.dataType = dataType;
		
		this.model = model;
		this.surfaces = surfaces;
		
		this.dataNormalizer = dataNormalizer;
		
		this.dataNormalizers = null;
		this.nIns = null;
		this.nOuts = null;
		
		if(this.dataType==DeepApproximation.Data.ALLSURFACES) {
			
			//Determination of the input size (just one here)
			this.nIns = new int[] {this.model.getParameters().length};
			
			//Determination of the output sizes
			int nOut = 0;
			
			Iterator<Map.Entry<String, WeightedOptionSurfaceData>> surfaceIterator = this.surfaces.entrySet().iterator();
			
			//Loop over the surfaces
			while (surfaceIterator.hasNext()) {
				
				Map.Entry<String, WeightedOptionSurfaceData> pair = 
						(Map.Entry<String, WeightedOptionSurfaceData>) surfaceIterator.next();
				
				OptionSurfaceData surface = pair.getValue();
				
				nOut = nOut + surface.getSurfaceSize();
			
			}
			
			//Only one output (all the prices/impl vols of all the surfaces stored in a 1D array)
			this.nOuts = new int[] {nOut};
			
		} else if(this.dataType==DeepApproximation.Data.EACHSURFACE) {
			
			//Determination of the input sizes
			this.nIns = new int[this.surfaces.size()];
			for(int i = 0; i < this.nIns.length; i++) 
				this.nIns[i] = this.model.getParameters().length;
			
			//Determination of the output sizes 
			this.nOuts = new int[this.surfaces.size()];
			int i = 0;
			
			this.dataNormalizers = new LinkedHashMap<String, DataNormalization>(this.surfaces.size());
			
			Iterator<Map.Entry<String, WeightedOptionSurfaceData>> surfaceIterator = this.surfaces.entrySet().iterator();
			
			//Loop over the surfaces
			while (surfaceIterator.hasNext()) {
				
				Map.Entry<String, WeightedOptionSurfaceData> pair = 
						(Map.Entry<String, WeightedOptionSurfaceData>) surfaceIterator.next();
				
				OptionSurfaceData surface = pair.getValue();
				
				this.nOuts[i] = surface.getSurfaceSize();
				
				String surfaceName = pair.getKey();
				
				this.dataNormalizers.put(surfaceName, this.dataNormalizer);
				
				i++;
			
			}	
		
		} else if(this.dataType==DeepApproximation.Data.PRICES) {
			
			//Determination of the input sizes taking into account (T,K) 
			this.nIns = new int[this.surfaces.size()];
			for(int i = 0; i < this.nIns.length; i++)
				this.nIns[i] = this.model.getParameters().length + 2;
			
			//Determination of the output sizes (nOut = 1 for each surface, just one price/impl vol)
			this.nOuts = new int[this.surfaces.size()];
			for(int i = 0; i < this.nOuts.length; i++)
				this.nOuts[i] = 1;
			
			this.dataNormalizers = new LinkedHashMap<String, DataNormalization>(this.surfaces.size());
			
			Iterator<Map.Entry<String, WeightedOptionSurfaceData>> surfaceIterator = this.surfaces.entrySet().iterator();
			
			//Loop over the surfaces
			while (surfaceIterator.hasNext()) {
				
				Map.Entry<String, WeightedOptionSurfaceData> pair = 
						(Map.Entry<String, WeightedOptionSurfaceData>) surfaceIterator.next();
				
				String surfaceName = pair.getKey();
				
				this.dataNormalizers.put(surfaceName, this.dataNormalizer);
				
			}
			
		} else if(this.dataType==DeepApproximation.Data.SMILES) {
			
			//Determination of the input sizes
			this.nIns = new int[this.surfaces.size()];
			//Determination of the output sizes 
			this.nOuts = new int[this.surfaces.size()];
			int i = 0;
			
			this.dataNormalizers = new LinkedHashMap<String, DataNormalization>(this.surfaces.size());
			
			Iterator<Map.Entry<String, WeightedOptionSurfaceData>> surfaceIterator = this.surfaces.entrySet().iterator();
			
			//Loop over the surfaces
			while (surfaceIterator.hasNext()) {
				
				Map.Entry<String, WeightedOptionSurfaceData> pair = 
						(Map.Entry<String, WeightedOptionSurfaceData>) surfaceIterator.next();
				
				String surfaceName = pair.getKey();
				
				this.dataNormalizers.put(surfaceName, this.dataNormalizer);
				
				OptionSurfaceData surface = pair.getValue();
				
				double[] maturities = surface.getMaturities();
				
				int smileSize = surface.getSmileSize(maturities[0]);
				
				for(int j = 1; j < maturities.length; j++) {
					if(smileSize!=surface.getSmileSize(maturities[j])) 
						throw new IllegalArgumentException("all smiles of the surface must share the same number of option quotes (same size)");
				}
				
				this.nIns[i] = this.model.getParameters().length + smileSize + 1;
				this.nOuts[i] = smileSize;
				i++;
				
			}
			
		}
			
	}
	
	public MultivariateCalibrableProcessInterface getModel() {
		return this.model;
	}
	
	public LinkedHashMap<String, WeightedOptionSurfaceData> getSurfaces() {
		return this.surfaces;
	}
	
	public DeepApproximation.Data getDataType() {
		return this.dataType;
	}
	
	public DataNormalization getDataNormalizer() {
		return this.dataNormalizer;
	}
	
	public LinkedHashMap<String, DataNormalization> getDataNormalizers() {
		return this.dataNormalizers;
	}
	
	public int[] getInputSizes() {
		return this.nIns;
	}
	
	public int[] getOutputSizes() {
		return this.nOuts;
	}
	
	/**
	 * runDeepLearningApproximation() is the main method of the class DeepApproximation.
	 * 
	 * It takes care of the training and testing of the neural network(s) with respect to the previously generated data set(s).
	 * 
	 * At the very end, the resulting neural network(s) will be stored according to the following:
	 * 
	 * 		- generateDataAllSurfaces(): one data set, then one neural network stored in modelAllSurfaces.zip;
	 * 
	 * 		- generateDataEachSurface(): several data sets, one per surface, then each neural network stored in modelSurface"name of the surface".zip;
	 * 
	 * 		- generateDataPrices(): several data sets, one per surface, then each neural network stored in modelPrices"name of the surface".zip;
	 * 
	 * 		- generateDataSmiles(): several data sets, one per surface, then each neural network stored in modelSmiles"name of the surface".zip.
	 * 
	 * The input parameters stand for the hyper-parameters controlling the training/testing architecture. 
	 * 
	 * Note that if there are several neural networks to train and test, 
	 * the corresponding architectures will share the same hyper-parameters:
	 * 
	 * @param miniBatchSizeTrainingSet //The size of a mini-batch for the training set;
	 * @param miniBatchSizeTestingSet //The size of a mini-batch for the testing set;
	 * @param ratioTrainTest //Ratio for the splitting of each data set into training and testing sets;
	 * @param numEpochs //The number of epochs;
	 * @param numHiddenLayers //The number of hidden layers;
	 * @param labelNormalization //Whether the labels has to be normalized or not.
	 * @throws IllegalArgumentException
	 * @throws IOException
	 */
	public void runDeepLearningApproximation(int miniBatchSizeTrainingSet, int miniBatchSizeTestingSet,
			double ratioTrainTest, int numEpochs, int numHiddenLayers, boolean labelNormalization) throws IllegalArgumentException, IOException {
		
		if(this.dataType==DeepApproximation.Data.ALLSURFACES) {
			
			//Retrieving of the previously generated data
			File dataFile = new File("surfaces.txt");
			DataSet dataSet = new DataSet();
			
			if(dataFile.exists()) {
				dataSet.load(dataFile);
			} else {
				throw new IllegalArgumentException("generateDataAllSurfaces() of DataGenerator must be performed before training");
			}
			
			//Splitting the data set into training and testing parts
			SplitTestAndTrain splitter = dataSet.splitTestAndTrain(ratioTrainTest);
			DataSet trainingSet = splitter.getTrain();
			DataSet testingSet = splitter.getTest();
			
			//Creation of the data set iterators for training and testing
			DataSetIterator trainIterator = new org.nd4j.linalg.dataset.api.iterator.SamplingDataSetIterator(trainingSet,
					miniBatchSizeTrainingSet, trainingSet.numExamples());
			DataSetIterator testIterator = new org.nd4j.linalg.dataset.api.iterator.SamplingDataSetIterator(testingSet, 
					miniBatchSizeTestingSet, testingSet.numExamples());
			
			//Normalization of the training data
			this.dataNormalizer.fitLabel(labelNormalization);
			this.dataNormalizer.fit(trainIterator);
			trainIterator.setPreProcessor(this.dataNormalizer);
				
			//Setting up the neural network configuration
			NeuralNetConfiguration.Builder neuralNetConf = new NeuralNetConfiguration.Builder();
				
			//Setting the seed of the stochastic optimization algorithm used for training the neural network
			neuralNetConf = neuralNetConf.seed(System.currentTimeMillis());
				
			//Initialization of the weights of the neural network (same scheme used for all layers)
			neuralNetConf = neuralNetConf.weightInit(WeightInit.XAVIER);
				
			//Choice of the stochastic optimization algorithm that will be used for training
			neuralNetConf = neuralNetConf.optimizationAlgo(OptimizationAlgorithm.STOCHASTIC_GRADIENT_DESCENT);
				
			//Specification of the stochastic gradient descent in use, here we choose the Adam scheme introduced by Kingma and Ba (2017):
			//Adam: A Method for Stochastic Optimization, https://arxiv.org/abs/1412.6980
			neuralNetConf = neuralNetConf.updater(new Adam());
				
			//L1 regularization and plus one of L2 kind (weight decay)
			neuralNetConf = neuralNetConf.l1(0.01).l2(0.01);
				
			//Adding layers
			NeuralNetConfiguration.ListBuilder neuralNetLayerList = neuralNetConf.list();
				
			//Input layer
			neuralNetLayerList = neuralNetLayerList.layer(0, new DenseLayer.Builder()
							.nIn(this.nIns[0])
							.nOut(30)
							.activation(Activation.ELU)
							//.dropOut()
							.build());
						
			//Hidden layers
			for(int i = 1; i < numHiddenLayers + 1; i++) {

				neuralNetLayerList = neuralNetLayerList.layer(i, new DenseLayer.Builder()
									.nIn(30)
									.nOut(30)
									.activation(Activation.ELU)
									//.dropOut()
									.build());

			}
						
			//Output layer
			neuralNetLayerList = neuralNetLayerList.layer(numHiddenLayers + 1, new OutputLayer.Builder(LossFunction.L2)
							.nIn(30)
							.nOut(this.nOuts[0])
							.activation(Activation.SIGMOID)
							//.dropOut()
							.build());
				
			//Building of the neural configuration
			MultiLayerConfiguration conf = neuralNetLayerList.build();
				
			//Creation of the neural network for the freshly generated configuration
			MultiLayerNetwork neuralNetwork = new MultiLayerNetwork(conf);
				
			//Initialization of the neural network
			neuralNetwork.init();
				
			//Setting up the user interface
			org.deeplearning4j.ui.api.UIServer uiServer = org.deeplearning4j.ui.api.UIServer.getInstance();
			org.deeplearning4j.core.storage.StatsStorage statsStorage = new org.deeplearning4j.ui.model.storage.InMemoryStatsStorage(); 
			uiServer.attach(statsStorage);
				
			//Printing the score of the loss function during training every 10 iterations
			//Storing the score of the loss function during training into a file at every iteration
			neuralNetwork.setListeners(new org.deeplearning4j.ui.model.stats.StatsListener(statsStorage),
					new ScoreIterationListener(), new CollectScoresIterationListener());
				
			//Training of the neural network
			neuralNetwork.fit(trainIterator, numEpochs);
			
			//Normalization of the testing data before evaluation
			this.dataNormalizer.fitLabel(labelNormalization);
			this.dataNormalizer.fit(testIterator);
			testIterator.setPreProcessor(this.dataNormalizer);
			
			//Evaluation of the model and then printing of the metrics
			RegressionEvaluation eval = neuralNetwork.evaluateRegression(testIterator);
			System.out.println(eval.stats());
			
			//Saving of the neural network model along with its data normalizer
			File modelFile = new File("modelAllSurfaces.zip");
					
			if(modelFile.createNewFile()) {
				ModelSerializer.writeModel(neuralNetwork, modelFile, true);
				ModelSerializer.addNormalizerToModel(modelFile, this.dataNormalizer);
			} else {
				modelFile.delete();
				modelFile.createNewFile();
				ModelSerializer.writeModel(neuralNetwork, modelFile, true);
				ModelSerializer.addNormalizerToModel(modelFile, this.dataNormalizer);
			}
			
		} else if(this.dataType==DeepApproximation.Data.EACHSURFACE) {
			
			//Iterator pointing at all the surfaces
			Iterator<Map.Entry<String, WeightedOptionSurfaceData>> iterator = this.surfaces.entrySet().iterator();
			int l = 0;
					
			//Loop over all the surfaces
			while (iterator.hasNext()) {
						
				Map.Entry<String, WeightedOptionSurfaceData> pair = 
						(Map.Entry<String, WeightedOptionSurfaceData>) iterator.next();
								
				String underlying = pair.getKey();
			
				//Retrieving of the previously generated data
				File dataFile = new File("surface" + underlying + ".txt");
				DataSet dataSet = new DataSet();
			
				if(dataFile.exists()) {
					dataSet.load(dataFile);
				} else {
					throw new IllegalArgumentException("generateDataEachSurface() of DataGenerator must be performed before training");
				}
			
				//Splitting the data set into training and testing parts
				SplitTestAndTrain splitter = dataSet.splitTestAndTrain(ratioTrainTest);
				DataSet trainingSet = splitter.getTrain();
				DataSet testingSet = splitter.getTest();
			
				//Creation of the data set iterators for training and testing
				DataSetIterator trainIterator = new org.nd4j.linalg.dataset.api.iterator.SamplingDataSetIterator(trainingSet,
						miniBatchSizeTrainingSet, trainingSet.numExamples());
				DataSetIterator testIterator = new org.nd4j.linalg.dataset.api.iterator.SamplingDataSetIterator(testingSet, 
						miniBatchSizeTestingSet, testingSet.numExamples());
				
				//Normalization of the training data
				DataNormalization dataNormalizer = this.dataNormalizers.get(underlying);
				dataNormalizer.fitLabel(labelNormalization);
				dataNormalizer.fit(trainIterator);
				trainIterator.setPreProcessor(dataNormalizer);
				
				//Setting up the neural network configuration
				NeuralNetConfiguration.Builder neuralNetConf = new NeuralNetConfiguration.Builder();
				
				//Setting the seed of the stochastic optimization algorithm used for training the neural network
				neuralNetConf = neuralNetConf.seed(System.currentTimeMillis());
				
				//Initialization of the weights of the neural network (same scheme used for all layers)
				neuralNetConf = neuralNetConf.weightInit(WeightInit.XAVIER);
				
				//Choice of the stochastic optimization algorithm that will be used for training
				neuralNetConf = neuralNetConf.optimizationAlgo(OptimizationAlgorithm.STOCHASTIC_GRADIENT_DESCENT);
				
				//Specification of the stochastic gradient descent in use, here we choose the Adam scheme introduced by Kingma and Ba (2017):
				//Adam: A Method for Stochastic Optimization, https://arxiv.org/abs/1412.6980
				neuralNetConf = neuralNetConf.updater(new Adam());
				
				//L1 regularization and plus one of L2 kind (weight decay)
				neuralNetConf = neuralNetConf.l1(0.01).l2(0.01);
				
				//Adding layers
				NeuralNetConfiguration.ListBuilder neuralNetLayerList = neuralNetConf.list();
				
				//Input layer
				neuralNetLayerList = neuralNetLayerList.layer(0, new DenseLayer.Builder()
								.nIn(this.nIns[l])
								.nOut(30)
								.activation(Activation.ELU)
								//.dropOut()
								.build());
						
				//Hidden layers
				for(int i = 1; i < numHiddenLayers + 1; i++) {

					neuralNetLayerList = neuralNetLayerList.layer(i, new DenseLayer.Builder()
										.nIn(30)
										.nOut(30)
										.activation(Activation.ELU)
										//.dropOut()
										.build());

				}
						
				//Output layer
				neuralNetLayerList = neuralNetLayerList.layer(numHiddenLayers + 1, new OutputLayer.Builder(LossFunction.L2)
								.nIn(30)
								.nOut(this.nOuts[l])
								.activation(Activation.SIGMOID)
								//.dropOut()
								.build());
				
				//Building of the neural configuration
				MultiLayerConfiguration conf = neuralNetLayerList.build();
				
				//Creation of the neural network for the freshly generated configuration
				MultiLayerNetwork neuralNetwork = new MultiLayerNetwork(conf);
				
				//Initialization of the neural network
				neuralNetwork.init();
				
				//Setting up the user interface
				org.deeplearning4j.ui.api.UIServer uiServer = org.deeplearning4j.ui.api.UIServer.getInstance();
				org.deeplearning4j.core.storage.StatsStorage statsStorage = new org.deeplearning4j.ui.model.storage.InMemoryStatsStorage(); 
				uiServer.attach(statsStorage);
				
				//Printing the score of the loss function during training every 10 iterations
				//Storing the score of the loss function during training into a file at every iteration
				neuralNetwork.setListeners(new org.deeplearning4j.ui.model.stats.StatsListener(statsStorage),
						new ScoreIterationListener(), new CollectScoresIterationListener());
				
				//Training of the neural network
				neuralNetwork.fit(trainIterator, numEpochs);
				
				//Normalization of the testing data before evaluation
				dataNormalizer.fitLabel(labelNormalization);
				dataNormalizer.fit(testIterator);
				testIterator.setPreProcessor(dataNormalizer);
			
				//Evaluation of the model and then printing of the metrics
				RegressionEvaluation eval = neuralNetwork.evaluateRegression(testIterator);
				System.out.println(eval.stats());
			
				//Saving of the neural network model along with its data normalizer
				File modelFile = new File("modelSurface" + underlying + ".zip");
					
				if(modelFile.createNewFile()) {
					ModelSerializer.writeModel(neuralNetwork, modelFile, true);
					ModelSerializer.addNormalizerToModel(modelFile, dataNormalizer);
				} else {
					modelFile.delete();
					modelFile.createNewFile();
					ModelSerializer.writeModel(neuralNetwork, modelFile, true);
					ModelSerializer.addNormalizerToModel(modelFile, dataNormalizer);
				}
				
				l++;

			}
		
		} else if(this.dataType==DeepApproximation.Data.PRICES) {
			
			//Iterator pointing at all the surfaces
			Iterator<Map.Entry<String, WeightedOptionSurfaceData>> iterator = this.surfaces.entrySet().iterator();
			int l = 0;
					
			//Loop over all surfaces and retrieving of all above data needed for random generation
			while (iterator.hasNext()) {
						
				Map.Entry<String, WeightedOptionSurfaceData> pair = 
						(Map.Entry<String, WeightedOptionSurfaceData>) iterator.next();
								
				String underlying = pair.getKey();
			
				//Retrieving of the previously generated data
				File dataFile = new File("prices" + underlying + ".txt");
				DataSet dataSet = new DataSet();
			
				if(dataFile.exists()) {
					dataSet.load(dataFile);
				} else {
					throw new IllegalArgumentException("generateDataPrices() of DataGenerator must be performed before training");
				}
			
				//Splitting the data set into training and testing parts
				SplitTestAndTrain splitter = dataSet.splitTestAndTrain(ratioTrainTest);
				DataSet trainingSet = splitter.getTrain();
				DataSet testingSet = splitter.getTest();
			
				//Creation of the data set iterators for training and testing
				DataSetIterator trainIterator = new org.nd4j.linalg.dataset.api.iterator.SamplingDataSetIterator(trainingSet,
						miniBatchSizeTrainingSet, trainingSet.numExamples());
				DataSetIterator testIterator = new org.nd4j.linalg.dataset.api.iterator.SamplingDataSetIterator(testingSet, 
						miniBatchSizeTestingSet, testingSet.numExamples());
				
				//Normalization of the training data
				DataNormalization dataNormalizer = this.dataNormalizers.get(underlying);
				dataNormalizer.fitLabel(labelNormalization);
				dataNormalizer.fit(trainIterator);
				trainIterator.setPreProcessor(dataNormalizer);
				
				//Setting up the neural network configuration
				NeuralNetConfiguration.Builder neuralNetConf = new NeuralNetConfiguration.Builder();
				
				//Setting the seed of the stochastic optimization algorithm used for training the neural network
				neuralNetConf = neuralNetConf.seed(System.currentTimeMillis());
				
				//Initialization of the weights of the neural network (same scheme used for all layers)
				neuralNetConf = neuralNetConf.weightInit(WeightInit.XAVIER);
				
				//Choice of the stochastic optimization algorithm that will be used for training
				neuralNetConf = neuralNetConf.optimizationAlgo(OptimizationAlgorithm.STOCHASTIC_GRADIENT_DESCENT);
				
				//Specification of the stochastic gradient descent in use, here we choose the Adam scheme introduced by Kingma and Ba (2017):
				//Adam: A Method for Stochastic Optimization, https://arxiv.org/abs/1412.6980
				neuralNetConf = neuralNetConf.updater(new Adam());
				
				//L1 regularization and plus one of L2 kind (weight decay)
				neuralNetConf = neuralNetConf.l1(0.01).l2(0.01);
				
				//Adding layers
				NeuralNetConfiguration.ListBuilder neuralNetLayerList = neuralNetConf.list();
				
				//Input layer
				neuralNetLayerList = neuralNetLayerList.layer(0, new DenseLayer.Builder()
								.nIn(this.nIns[l])
								.nOut(50)
								.activation(Activation.ELU)
								//.dropOut()
								.build());
						
				//Hidden layers
				for(int i = 1; i < numHiddenLayers + 1; i++) {

					neuralNetLayerList = neuralNetLayerList.layer(i, new DenseLayer.Builder()
										.nIn(50)
										.nOut(50)
										.activation(Activation.ELU)
										//.dropOut()
										.build());

				}
						
				//Output layer
				neuralNetLayerList = neuralNetLayerList.layer(numHiddenLayers + 1, new OutputLayer.Builder(LossFunction.L2)
								.nIn(50)
								.nOut(this.nOuts[l])
								.activation(Activation.SIGMOID)
								//.dropOut()
								.build());
				
				//Building of the neural configuration
				MultiLayerConfiguration conf = neuralNetLayerList.build();
				
				//Creation of the neural network for the freshly generated configuration
				MultiLayerNetwork neuralNetwork = new MultiLayerNetwork(conf);
				
				//Initialization of the neural network
				neuralNetwork.init();
				
				//Setting up the user interface
				org.deeplearning4j.ui.api.UIServer uiServer = org.deeplearning4j.ui.api.UIServer.getInstance();
				org.deeplearning4j.core.storage.StatsStorage statsStorage = new org.deeplearning4j.ui.model.storage.InMemoryStatsStorage(); 
				uiServer.attach(statsStorage);
				
				//Printing the score of the loss function during training every 10 iterations
				//Storing the score of the loss function during training into a file at every iteration
				neuralNetwork.setListeners(new org.deeplearning4j.ui.model.stats.StatsListener(statsStorage),
						new ScoreIterationListener(), new CollectScoresIterationListener());
				
				//Training of the neural network
				neuralNetwork.fit(trainIterator, numEpochs);
				
				//Normalization of the testing data before evaluation
				dataNormalizer.fitLabel(labelNormalization);
				dataNormalizer.fit(testIterator);
				testIterator.setPreProcessor(dataNormalizer);
			
				//Evaluation of the model and then printing of the metrics
				RegressionEvaluation eval = neuralNetwork.evaluateRegression(testIterator);
				System.out.println(eval.stats());
			
				//Saving of the neural network model along with its data normalizer
				File modelFile = new File("modelPrices" + underlying + ".zip");
					
				if(modelFile.createNewFile()) {
					ModelSerializer.writeModel(neuralNetwork, modelFile, true);
					ModelSerializer.addNormalizerToModel(modelFile, dataNormalizer);
				} else {
					modelFile.delete();
					modelFile.createNewFile();
					ModelSerializer.writeModel(neuralNetwork, modelFile, true);
					ModelSerializer.addNormalizerToModel(modelFile, dataNormalizer);
				}
				
				l++;

			}
			
		} else if(this.dataType==DeepApproximation.Data.SMILES) {
			
			//Iterator pointing at all the surfaces
			Iterator<Map.Entry<String, WeightedOptionSurfaceData>> iterator = this.surfaces.entrySet().iterator();
			int l = 0;
					
			//Loop over all surfaces and retrieving of all above data needed for random generation
			while (iterator.hasNext()) {
						
				Map.Entry<String, WeightedOptionSurfaceData> pair = 
						(Map.Entry<String, WeightedOptionSurfaceData>) iterator.next();
								
				String underlying = pair.getKey();
			
				//Retrieving of the previously generated data
				File dataFile = new File("smiles" + underlying + ".txt");
				DataSet dataSet = new DataSet();
			
				if(dataFile.exists()) {
					dataSet.load(dataFile);
				} else {
					throw new IllegalArgumentException("generateDataSmiles() of DataGenerator must be performed before training");
				}
			
				//Splitting the data set into training and testing parts
				SplitTestAndTrain splitter = dataSet.splitTestAndTrain(ratioTrainTest);
				DataSet trainingSet = splitter.getTrain();
				DataSet testingSet = splitter.getTest();
			
				//Creation of the data set iterators for training and testing
				DataSetIterator trainIterator = new org.nd4j.linalg.dataset.api.iterator.SamplingDataSetIterator(trainingSet,
						miniBatchSizeTrainingSet, trainingSet.numExamples());
				DataSetIterator testIterator = new org.nd4j.linalg.dataset.api.iterator.SamplingDataSetIterator(testingSet, 
						miniBatchSizeTestingSet, testingSet.numExamples());
				
				//Normalization of the training data
				DataNormalization dataNormalizer = this.dataNormalizers.get(underlying);
				dataNormalizer.fitLabel(labelNormalization);
				dataNormalizer.fit(trainIterator);
				trainIterator.setPreProcessor(dataNormalizer);
				
				//Setting up the neural network configuration
				NeuralNetConfiguration.Builder neuralNetConf = new NeuralNetConfiguration.Builder();
				
				//Setting the seed of the stochastic optimization algorithm used for training the neural network
				neuralNetConf = neuralNetConf.seed(System.currentTimeMillis());
				
				//Initialization of the weights of the neural network (same scheme used for all layers)
				neuralNetConf = neuralNetConf.weightInit(WeightInit.XAVIER);
				
				//Choice of the stochastic optimization algorithm that will be used for training
				neuralNetConf = neuralNetConf.optimizationAlgo(OptimizationAlgorithm.STOCHASTIC_GRADIENT_DESCENT);
				
				//Specification of the stochastic gradient descent in use, here we choose the Adam scheme introduced by Kingma and Ba (2017):
				//Adam: A Method for Stochastic Optimization, https://arxiv.org/abs/1412.6980
				neuralNetConf = neuralNetConf.updater(new Adam());
				
				//L1 regularization and plus one of L2 kind (weight decay)
				neuralNetConf = neuralNetConf.l1(0.01).l2(0.01);
				
				//Adding layers
				NeuralNetConfiguration.ListBuilder neuralNetLayerList = neuralNetConf.list();
				
				//Input layer
				neuralNetLayerList = neuralNetLayerList.layer(0, new DenseLayer.Builder()
								.nIn(this.nIns[l])
								.nOut(60)
								.activation(Activation.ELU)
								//.dropOut()
								.build());
						
				//Hidden layers
				for(int i = 1; i < numHiddenLayers + 1; i++) {

					neuralNetLayerList = neuralNetLayerList.layer(i, new DenseLayer.Builder()
										.nIn(60)
										.nOut(60)
										.activation(Activation.ELU)
										//.dropOut()
										.build());

				}
						
				//Output layer
				neuralNetLayerList = neuralNetLayerList.layer(numHiddenLayers + 1, new OutputLayer.Builder(LossFunction.L2)
								.nIn(60)
								.nOut(this.nOuts[l])
								.activation(Activation.SIGMOID)
								//.dropOut()
								.build());
				
				//Building of the neural configuration
				MultiLayerConfiguration conf = neuralNetLayerList.build();
				
				//Creation of the neural network for the freshly generated configuration
				MultiLayerNetwork neuralNetwork = new MultiLayerNetwork(conf);
				
				//Initialization of the neural network
				neuralNetwork.init();
				
				//Setting up the user interface
				org.deeplearning4j.ui.api.UIServer uiServer = org.deeplearning4j.ui.api.UIServer.getInstance();
				org.deeplearning4j.core.storage.StatsStorage statsStorage = new org.deeplearning4j.ui.model.storage.InMemoryStatsStorage(); 
				uiServer.attach(statsStorage);
				
				//Printing the score of the loss function during training every 10 iterations
				//Storing the score of the loss function during training into a file at every iteration
				neuralNetwork.setListeners(new org.deeplearning4j.ui.model.stats.StatsListener(statsStorage),
						new ScoreIterationListener(), new CollectScoresIterationListener());
				
				//Training of the neural network
				neuralNetwork.fit(trainIterator, numEpochs);
				
				//Normalization of the testing data before evaluation
				dataNormalizer.fitLabel(labelNormalization);
				dataNormalizer.fit(testIterator);
				testIterator.setPreProcessor(dataNormalizer);
			
				//Evaluation of the model and then printing of the metrics
				RegressionEvaluation eval = neuralNetwork.evaluateRegression(testIterator);
				System.out.println(eval.stats());
			
				//Saving of the neural network model along with its data normalizer
				File modelFile = new File("modelSmiles" + underlying + ".zip");
					
				if(modelFile.createNewFile()) {
					ModelSerializer.writeModel(neuralNetwork, modelFile, true);
					ModelSerializer.addNormalizerToModel(modelFile, dataNormalizer);
				} else {
					modelFile.delete();
					modelFile.createNewFile();
					ModelSerializer.writeModel(neuralNetwork, modelFile, true);
					ModelSerializer.addNormalizerToModel(modelFile, dataNormalizer);
				}
				
				l++;

			}
			
		}
			
	}
	
	public enum Data {
		ALLSURFACES,
		EACHSURFACE,
		PRICES,
		SMILES
	}

}
