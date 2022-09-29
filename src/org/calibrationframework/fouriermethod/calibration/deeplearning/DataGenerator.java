package org.calibrationframework.fouriermethod.calibration.deeplearning;

import java.io.*;
import java.util.*;

import org.apache.commons.lang3.ArrayUtils;

import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.distribution.*;

import org.calibrationframework.exception.*;
import org.calibrationframework.fouriermethod.calibration.models.*;
import org.calibrationframework.fouriermethod.products.*;
import org.calibrationframework.marketdata.model.volatilities.*;
import org.calibrationframework.marketdata.model.volatilities.VolatilitySurfaceInterface.QuotingConvention;

import org.nd4j.linalg.dataset.DataSet;
import org.deeplearning4j.datasets.iterator.DoublesDataSetIterator;
import org.nd4j.common.primitives.Pair;
import org.nd4j.linalg.dataset.api.iterator.DataSetIterator;

/**
 * This class consists in the random generation of training set(s) for the deep-learning approximation of pricing function(s) (via neural networks).
 * 
 * There exist in the literature different approaches to the deep calibration problem. Most of them have been implemented in the present class.
 * 
 * @author Szulda Guillaume
 *
 */
public class DataGenerator {
	
	private final MultivariateCalibrableProcessInterface model;
	private final EuropeanOptionSmileMultiAsset pricer;
	private final LinkedHashMap<String, WeightedOptionSurfaceData> surfaces;
	
	/**
	 * Constructor for the generation of training set(s) for the deep-learning approximation of pricing function(s).
	 * 
	 * The surfaces are needed for:
	 * 		- the retrieving of the fixed grids of maturities and strikes regarding generateDataAllSurfaces() and generateDataEachSurface();
	 * 		- the estimation of the density of (T,K) that has to be consistent with the market data for generateDataPrices();
	 * 		- the estimation of the density of (T, K1, ..., Kn) that has to be consistent with the market data for generateDataSmiles().
	 * 
	 * @param model
	 * @param pricer
	 * @param surfaces
	 */
	public DataGenerator(MultivariateCalibrableProcessInterface model, EuropeanOptionSmileMultiAsset pricer, 
			LinkedHashMap<String, WeightedOptionSurfaceData> surfaces) {
		
		this.model = model;
		this.pricer = pricer;
		this.surfaces = surfaces;
		
	}
	
	/**
	 * generateDataAllSurfaces() follows Horvath et al. "Deep Learning Volatility" (2021).
	 * 
	 * It allows to randomly generate a data set of chosen size, 
	 * consisting of couples of the form: array of parameter - model prices/impl vols of all considered surfaces, all of which are stored in a 1D array.
	 * 
	 * For each couple, the parameter array, generated randomly, stands for the parameters of the underlying model.
	 * It is coupled with a 1D array, representing all the prices/impl vols of all considered surfaces, computed through the pricing model for this set of parameters.
	 * 
	 * This approach then turns out to consider fixed grids of maturities and strikes throughout the generation.
	 * These fixed grids have to be the same as the market ones over which the model is wished to be calibrated.
	 * 
	 * The resulting data set, which will be saved inside a file named surfaces.txt, will be used for determining the neural network approximation of the underlying pricing function,
	 * the one that takes the parameter set as input and returns, as a 1D output array, all the surfaces of model prices/impl vols.
  	 *
	 * @param sizeDataSet
	 * @throws IOException
	 */
	public void generateDataAllSurfaces(int sizeDataSet) throws IOException {
		
		//Creation of the new data list that will be used for training the neural network
		List< Pair< double[], double[] > > dataList = new ArrayList< Pair< double[], double[] > >(sizeDataSet);
				
		//Generation of all data sample couples: parameter set - model prices/impl vols of all surfaces in the form of 1D array
		for(int k = 0; k < sizeDataSet; k++) {
			
			long startMillis = System.currentTimeMillis();
					
			//Retrieving of the k-th couple: parameter set - corresponding MultivariateCalibrableProcessInterface model
			Pair< List<Double>, MultivariateCalibrableProcessInterface > currentPair = this.model.generateSamplePair();
					
			List<Double> newParams = currentPair.getFirst();
					
			double[] params = ArrayUtils.toPrimitive( newParams.toArray(new Double[newParams.size()]) );

			MultivariateCalibrableProcessInterface newModel = currentPair.getSecond();
					
			//Starting point of the computation of the k-th output of the pricing function
			//Iterator pointing at the target surfaces
			Iterator<Map.Entry<String, WeightedOptionSurfaceData>> iterator = this.surfaces.entrySet().iterator();
					
			//The k-th output of the pricing function
			List<Double> vals = new ArrayList<Double>();
					
			//Loop over the surfaces
			while (iterator.hasNext()) {
				
				Map.Entry<String, WeightedOptionSurfaceData> pair = 
						(Map.Entry<String, WeightedOptionSurfaceData>) iterator.next();
						
				String underlying = pair.getKey();
						
				OptionSurfaceData surface = pair.getValue();
						
				int numberOfMaturities = surface.getMaturities().length;
						
				double mats[] = surface.getMaturities();
						
				QuotingConvention targetConvention = surface.getQuotingConvention();
						
				for(int t = 0; t<numberOfMaturities; t++) {
							
					double[] currentStrikes = surface.getSmile(mats[t]).getStrikes();
							
					EuropeanOptionSmileMultiAsset newPricer = this.pricer.getCloneWithModifiedParameters(underlying, mats[t],currentStrikes);
							
					try {
								
						Map<Double, Double> currentModelPrices = newPricer.getValue(newModel);
								
						for(int i = 0; i < currentStrikes.length; i++) {
						
							if(targetConvention.equals(QuotingConvention.VOLATILITYLOGNORMAL)) {
										
								//we convert prices into lognormal volatilities
								double forward = surface.getEquityForwardCurve().getDiscountFactor(mats[t]);
								double optionMaturity = mats[t];
								double optionStrike = currentStrikes[i];
								double payoffUnit = surface.getDiscountCurve().getDiscountFactor(mats[t]);
								double optionValue = currentModelPrices.get(currentStrikes[i]);
								double implVol = Math.abs(net.finmath.functions.AnalyticFormulas.blackScholesOptionImpliedVolatility(forward, optionMaturity, optionStrike, payoffUnit, optionValue));
								if(implVol < 1.0) {
									vals.add(implVol);
								} else {
									vals.add( surface.getSurface().get(optionMaturity).getSmile().get(optionStrike).getValue() );
								}
								
							} else if(targetConvention.equals(QuotingConvention.VOLATILITYNORMAL)) {
										
								//we convert prices into normal volatilities
								double forward = surface.getEquityForwardCurve().getDiscountFactor(mats[t]);
								double optionMaturity =mats[t];
								double optionStrike = currentStrikes[i];
								double payoffUnit = surface.getDiscountCurve().getDiscountFactor(mats[t]);
								double optionValue = currentModelPrices.get(currentStrikes[i]);
								double implVol = Math.abs(net.finmath.functions.AnalyticFormulas.bachelierOptionImpliedVolatility(forward, optionMaturity, optionStrike, payoffUnit, optionValue));
								if(implVol < 1.0) {
									vals.add(implVol);
								} else {
									vals.add( surface.getSurface().get(optionMaturity).getSmile().get(optionStrike).getValue() );
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
			System.out.println("Paramters:");
			System.out.println(newParams);
			System.out.println("Implied volatilites:");
			System.out.println(vals);
			//k-th output of the pricing function in the form of an array
			double[] values = ArrayUtils.toPrimitive( vals.toArray(new Double[vals.size()]) );
					
			dataList.add( Pair.< double[], double[] >create(params, values) );
			
			long endMillis = System.currentTimeMillis();
			
			double calculationTime = ((endMillis-startMillis)/1000.0);
			
			System.out.println(k + "-th data sample couple retrieved in: " + calculationTime + " seconds");
					
		}
				
		//Creation of the data set 
		DataSetIterator iterator = new DoublesDataSetIterator(dataList, sizeDataSet);
		DataSet dataSet = iterator.next();
		
		//Creation of the file in which the data set will be saved
		File dataFile = new File("surfaces.txt");
		
		//Check whether data.txt already exists or not
		if(dataFile.createNewFile()) {
			dataSet.save(dataFile);
		} else {
			dataFile.delete();
			dataFile.createNewFile();
			dataSet.save(dataFile);
		}
		
	}
	
	/**
	 * generateDataEachSurface() follows Horvath et al. "Deep Learning Volatility" (2021).
	 * 
	 * It randomly generates as many data sets as there are surfaces over which the model has to be calibrated, all data sets of the same chosen size.
	 * Each of them consists of couples of the form: array of parameter - model prices/impl vols of one surface, all of which are stored in a 1D array.
	 * 
	 * In the same way as above, the considered grids of maturities and strikes will be fixed throughout the generation.
	 * These fixed grids have to be the same as the market ones over which the model is wished to be calibrated.
	 * 
	 * Each resulting data set, one per surface, will be saved inside a file named surface"name of the surface".txt, 
	 * and restored later for the training of the neural network approximating the pricing function of the surface.	
  	 *
	 * @param sizeDataSet
	 * @throws IOException
	 */
	public void generateDataEachSurface(int sizeDataSet) throws IOException {
		
		//Storing of the names of the surfaces in a list
		//Ordered like LinkedHashMap<String, WeightedOptionSurfaceData> surfaces
		List< String > names = new ArrayList< String >(this.surfaces.size());
				
		//Creation of the data lists that will be used for training the neural networks
		//As many data lists/neural networks as there are surfaces
		//In the same order as LinkedHashMap<String, WeightedOptionSurfaceData> surfaces
		List< List< Pair< double[], double[] > > > dataLists = 
				new ArrayList< List< Pair< double[], double[] > > >(this.surfaces.size());
										
		//Iterator pointing at all the surfaces
		Iterator<Map.Entry<String, WeightedOptionSurfaceData>> surfacesIterator = this.surfaces.entrySet().iterator();
				
		//Loop over all surfaces and retrieving of all above data needed for generation
		while (surfacesIterator.hasNext()) {
					
			Map.Entry<String, WeightedOptionSurfaceData> pair = 
					(Map.Entry<String, WeightedOptionSurfaceData>) surfacesIterator.next();
							
			String underlying = pair.getKey();
					
			//Storing of the name of the surface
			names.add(underlying);
					
			List< Pair< double[], double[] > > dataList = new ArrayList< Pair< double[], double[] > >(sizeDataSet);
					
			//Declaration of the data list for the current surface
			dataLists.add(dataList);
					
		}
				
		//For each surface, generation of all data sample couples: parameter set - surface of model prices/impl vols in the form of 1D array
		for(int k = 0; k < sizeDataSet; k++) {
			
			long startMillis = System.currentTimeMillis();
					
			//Retrieving of the k-th couple: parameter set - corresponding MultivariateCalibrableProcessInterface model
			Pair< List<Double>, MultivariateCalibrableProcessInterface > currentPair = this.model.generateSamplePair();
					
			List<Double> newParams = currentPair.getFirst();
					
			double[] params = ArrayUtils.toPrimitive( newParams.toArray(new Double[newParams.size()]) );

			MultivariateCalibrableProcessInterface newModel = currentPair.getSecond();
					
			//Loop going through all surfaces
			for(int l = 0; l < this.surfaces.size(); l++) {
				
				//The k-th output of the pricing function of the l-th surface
				List<Double> vals = new ArrayList<Double>();
				
				String surfaceName = names.get(l);
				
				OptionSurfaceData surface = this.surfaces.get(surfaceName);
				
				int numberOfMaturities = surface.getMaturities().length;
				
				double mats[] = surface.getMaturities();
						
				QuotingConvention targetConvention = surface.getQuotingConvention();
						
				for(int t = 0; t<numberOfMaturities; t++) {
							
					double[] currentStrikes = surface.getSmile(mats[t]).getStrikes();
							
					EuropeanOptionSmileMultiAsset newPricer = this.pricer.getCloneWithModifiedParameters(surfaceName,
							mats[t],currentStrikes);
							
					try {
								
						Map<Double, Double> currentModelPrices = newPricer.getValue(newModel);
								
						for(int i = 0; i < currentStrikes.length; i++) {
						
							if(targetConvention.equals(QuotingConvention.VOLATILITYLOGNORMAL)) {
										
								//we convert prices into lognormal volatilities
								double forward = surface.getEquityForwardCurve().getDiscountFactor(mats[t]);
								double optionMaturity = mats[t];
								double optionStrike = currentStrikes[i];
								double payoffUnit = surface.getDiscountCurve().getDiscountFactor(mats[t]);
								double optionValue = currentModelPrices.get(currentStrikes[i]);
								double implVol = Math.abs(net.finmath.functions.AnalyticFormulas.blackScholesOptionImpliedVolatility(forward, optionMaturity, optionStrike, payoffUnit, optionValue));
								if(implVol < 10.0) {
									vals.add(implVol);
								} else {
									vals.add( surface.getSurface().get(optionMaturity).getSmile().get(optionStrike).getValue() );
								}
								
							} else if(targetConvention.equals(QuotingConvention.VOLATILITYNORMAL)) {
										
								//we convert prices into normal volatilities
								double forward = surface.getEquityForwardCurve().getDiscountFactor(mats[t]);
								double optionMaturity =mats[t];
								double optionStrike = currentStrikes[i];
								double payoffUnit = surface.getDiscountCurve().getDiscountFactor(mats[t]);
								double optionValue = currentModelPrices.get(currentStrikes[i]);
								double implVol = Math.abs(net.finmath.functions.AnalyticFormulas.bachelierOptionImpliedVolatility(forward, optionMaturity, optionStrike, payoffUnit, optionValue));
								if(implVol < 10.0) {
									vals.add(implVol);
								} else {
									vals.add( surface.getSurface().get(optionMaturity).getSmile().get(optionStrike).getValue() );
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
				
				//k-th output of the l-th pricing function in the form of an array (l-th surface)
				double[] values = ArrayUtils.toPrimitive( vals.toArray(new Double[vals.size()]) );
				
				dataLists.get(l).add( Pair.< double[], double[] >create(params, values) );
						
			}
			
			long endMillis = System.currentTimeMillis();
			
			double calculationTime = ((endMillis-startMillis)/1000.0);
			
			System.out.println(k + "-th data sample couples for all surfaces retrieved in: " + calculationTime + " seconds");
					
		}
		
		for(int l = 0; l < this.surfaces.size(); l++) {
			
			//Creation of the data set 
			DataSetIterator iterator = new DoublesDataSetIterator(dataLists.get(l), sizeDataSet);
			DataSet dataSet = iterator.next();
			
			//Creation of the file in which the data set will be saved
			File dataFile = new File("surface" + names.get(l) + ".txt");
			
			//Check whether data.txt already exists or not
			if(dataFile.createNewFile()) {
				dataSet.save(dataFile);
			} else {
				dataFile.delete();
				dataFile.createNewFile();
				dataSet.save(dataFile);
			}
			
		}
		
	}
	
    /**
     * generateDataPrices() follows Bayer et Stemper "Deep calibration of rough stochastic volatility models" (2018).
     * 
     * Similarly to generateDataEachSurface(), it randomly generates as many data sets as there are price/impl vol surfaces. One data set for each surface.
     * However, here for each surface, the one-output pricing function is approximated instead of the whole pricing function returning the complete surface of model prices/impl vols.
     * This means that for each generated data set, couples are of the form: array of parameter plus maturity and strike (input size = numOfParams + 2) - corresponding price/impl vol (output size = 1).
     * 
     * For each data set/surface, maturities and strikes are generated randomly from a joint distribution, the latter of which is estimated via weighted (bivariate Gaussian) kernel density estimation (wKDE).
     * Correlation is assumed between maturities and strikes, meaning that the bandwidth matrix reveals to be full.
     * The latter will be computed by means of Scott and Silverman's rule of thumb, where the co-variance matrix,
     * will be determined empirically through the market data of the surface in question, while taking into account the liquidity weight of every single option quote.
     *          
     * Finally, every single data set will be saved inside a file named prices"name of the surface".txt, 
     * which will be used for the determination of the neural network approximating the corresponding one-output pricing function (one for each surface).
     * We will then end up with as many neural network as there are surfaces of prices/impl vols.
     *
     * @param sizeDataSet
     * @throws IOException 
     */
	public void generateDataPrices(int sizeDataSet) throws IOException, IllegalArgumentException {
		
		//Storing of the names of the surfaces in a list
		//Ordered like LinkedHashMap<String, WeightedOptionSurfaceData> surfaces
		List< String > names = new ArrayList< String >(this.surfaces.size());
		
		//Creation of the data lists that will be used for training the neural networks
		//As many data lists/neural networks as there are surfaces
		//In the same order as LinkedHashMap<String, WeightedOptionSurfaceData> surfaces
		List< List< Pair< double[], double[] > > > dataLists = 
				new ArrayList< List< Pair< double[], double[] > > >(this.surfaces.size());
		
		//Determination of the density of (T,K) for every single surface of LinkedHashMap<String, WeightedOptionSurfaceData> surfaces
		//Performed through weighted kernel density estimation (wKDE) with Gaussian kernel and Scott and Silverman's rule of thumb for the selection of the full bandwidth matrix.
		//Again in the same order as LinkedHashMap<String, WeightedOptionSurfaceData> surfaces
		List< MixtureMultivariateNormalDistribution > densities = new ArrayList< MixtureMultivariateNormalDistribution >(this.surfaces.size());
								
		//Iterator pointing at all the surfaces
		Iterator<Map.Entry<String, WeightedOptionSurfaceData>> surfacesIterator = this.surfaces.entrySet().iterator();
		
		//Loop over all surfaces and retrieving of all above data needed for random generation
		while (surfacesIterator.hasNext()) {
			
			Map.Entry<String, WeightedOptionSurfaceData> pair = 
					(Map.Entry<String, WeightedOptionSurfaceData>) surfacesIterator.next();
					
			String underlying = pair.getKey();
			
			//Storing of the name of the surface
			names.add(underlying);
			
			List< Pair< double[], double[] > > dataList = new ArrayList< Pair< double[], double[] > >(sizeDataSet);
			
			//Declaration of the data list for the current surface
			dataLists.add(dataList);
			
			WeightedOptionSurfaceData weightedSurface = pair.getValue();
			
			//Computation of the wKDE density of (T,K) for this surface
			densities.add( bivariatewKDE(weightedSurface) );
			
		}
		
		//Generation of all data sample couples for all surfaces
		//Every single couple of the form: parameter set + maturity and strike - price/impl vol 
		for(int k = 0; k < sizeDataSet; k++) {
					
			long startMillis = System.currentTimeMillis();
							
			//Retrieving of the k-th couple: parameter set - corresponding MultivariateCalibrableProcessInterface model
			Pair< List<Double>, MultivariateCalibrableProcessInterface > currentPair = this.model.generateSamplePair();
			
			//Retrieving of the model
			MultivariateCalibrableProcessInterface newModel = currentPair.getSecond();
			
			//Loop over the different surfaces
			for(int l = 0; l < this.surfaces.size(); l++) {
				
				//Retrieving of the name
				String surfaceName = names.get(l);
				
				//Retrieving of the surface in question to get:
				//- Quoting convention;
				//- Discount curve;
				//- Forward curve.
				WeightedOptionSurfaceData weightedSurface = this.surfaces.get(surfaceName);
				
				//Retrieving of the sample parameter set of the model
				@SuppressWarnings("unchecked")
				List<Double> newParams = (List<Double>)((ArrayList<Double>)currentPair.getFirst()).clone();
				
				//Generation of a sample couple of the wKDE distribution of (T,K) for this surface
				double[] coupleTK = densities.get(l).sample();
				double maturity = Math.abs(coupleTK[0]);
				double strike = Math.abs(coupleTK[1]);
				while(maturity > newModel.getTimeHorizon()) {
					coupleTK = densities.get(l).sample();
					maturity = Math.abs(coupleTK[0]);
					strike = Math.abs(coupleTK[1]);
				}
				
				//Addition of the (T,K) sample couple to the sample parameter set of the model
				newParams.add(maturity);
				newParams.add(strike);
				
				//Declaration of the pricer
				EuropeanOptionSmileMultiAsset newPricer = this.pricer.getCloneWithModifiedParameters(
						surfaceName, maturity, new double[] {strike} );
				
				//This value will contain the output price/impl vol returned by the pricing function
				double value = 0;
				
				try {
					
					Map<Double, Double> modelPricesImplVols = newPricer.getValue(newModel);
					
					if(weightedSurface.getQuotingConvention().equals(QuotingConvention.VOLATILITYLOGNORMAL)) {
									
						double forward = weightedSurface.getEquityForwardCurve().getDiscountFactor(maturity);
						double optionMaturity = maturity;
						double optionStrike = strike;
						double payoffUnit = weightedSurface.getDiscountCurve().getDiscountFactor(maturity);
						double optionValue = modelPricesImplVols.get(strike);
						double implVol = Math.abs(net.finmath.functions.AnalyticFormulas.blackScholesOptionImpliedVolatility(forward, optionMaturity, optionStrike, payoffUnit, optionValue));
						if(implVol < 10.0) {
							value = implVol;
						} else {
							double distance1 = Math.abs(weightedSurface.getMaturities()[0] - optionMaturity);
							int idx1 = 0;
							for(int c = 1; c < weightedSurface.getMaturities().length; c++) {
								double cdistance = Math.abs(weightedSurface.getMaturities()[c] - optionMaturity);
								if(cdistance < distance1){
									idx1 = c;
									distance1 = cdistance;
								}
							}
							double[] strikes = weightedSurface.getSmile(weightedSurface.getMaturities()[idx1]).getStrikes();
							double distance2 = Math.abs(strikes[0] - optionStrike);
							int idx2 = 0;
							for(int c = 1; c < strikes.length; c++) {
								double cdistance = Math.abs(strikes[c] - optionStrike);
								if(cdistance < distance2){
									idx2 = c;
									distance2 = cdistance;
								}
							}
							value = weightedSurface.getSmile(weightedSurface.getMaturities()[idx1]).getOption(strikes[idx2]).getValue();	
						}
							
					} else if(weightedSurface.getQuotingConvention().equals(QuotingConvention.VOLATILITYNORMAL)) {
									
						double forward = weightedSurface.getEquityForwardCurve().getDiscountFactor(maturity);
						double optionMaturity = maturity;
						double optionStrike = strike;
						double payoffUnit = weightedSurface.getDiscountCurve().getDiscountFactor(maturity);
						double optionValue = modelPricesImplVols.get(strike);
						double implVol = Math.abs(net.finmath.functions.AnalyticFormulas.bachelierOptionImpliedVolatility(forward, optionMaturity, optionStrike, payoffUnit, optionValue));
						if(implVol < 10.0) {
							value = implVol;
						} else {
							double distance1 = Math.abs(weightedSurface.getMaturities()[0] - optionMaturity);
							int idx1 = 0;
							for(int c = 1; c < weightedSurface.getMaturities().length; c++) {
								double cdistance = Math.abs(weightedSurface.getMaturities()[c] - optionMaturity);
								if(cdistance < distance1){
									idx1 = c;
									distance1 = cdistance;
								}
							}
							double[] strikes = weightedSurface.getSmile(weightedSurface.getMaturities()[idx1]).getStrikes();
							double distance2 = Math.abs(strikes[0] - optionStrike);
							int idx2 = 0;
							for(int c = 1; c < strikes.length; c++) {
								double cdistance = Math.abs(strikes[c] - optionStrike);
								if(cdistance < distance2){
									idx2 = c;
									distance2 = cdistance;
								}
							}
							value = weightedSurface.getSmile(weightedSurface.getMaturities()[idx1]).getOption(strikes[idx2]).getValue();	
						}
									
					} else {
							
						value = modelPricesImplVols.get(strike);
									
					}						
							
				} catch (CalculationException e) {
					e.printStackTrace();
				}
				
				double[] params = ArrayUtils.toPrimitive( newParams.toArray(new Double[newParams.size()]) );
				double[] values = {value};
				
				dataLists.get(l).add( Pair.< double[], double[] >create(params, values) );
				
			}
							
			long endMillis = System.currentTimeMillis();
			
			double calculationTime = ((endMillis-startMillis)/1000.0);
			
			System.out.println(k + "-th data sample couples for all surfaces retrieved in: " + calculationTime + " seconds");
			
		}
		
		//Loop over the different surfaces
		for(int p = 0; p < this.surfaces.size(); p++) {
			
			//Creation of the data set 
			DataSetIterator iterator = new DoublesDataSetIterator(dataLists.get(p), sizeDataSet);
			DataSet dataSet = iterator.next();
			
			//Creation of the file in which the data set will be saved
			File dataFile = new File("prices" + names.get(p) + ".txt");
			
			//Check whether data.txt already exists or not
			if(dataFile.createNewFile()) {
				dataSet.save(dataFile);
			} else {
				dataFile.delete();
				dataFile.createNewFile();
				dataSet.save(dataFile);
			}
			
		}
		
	}
	
	/**
	 * This private method returns the bivariate weighted Kernel Density Estimation (wKDE) of the density of (T,K) (couple maturity - strike) for the input weighted surface.
	 * The adjective weighted here means that the standard KDE formula is computed by taking into account the weight of every single (T,K),
	 * which measures the liquidity on the market of the corresponding option.
	 * The latter is observable in the formula where it replaces (1/n) for every single summand.
	 * Finally, Gaussian Kernel will be used for computations along with Scott and Silverman's rule of thumb for bandwidth selection.
	 * 
	 * @param weightedSurface
	 * @return
	 */
	private MixtureMultivariateNormalDistribution bivariatewKDE(WeightedOptionSurfaceData weightedSurface) {
		
		//Preparation of the inputs required for the creation of the wKDE density of (T,K) for this surface
		List< org.apache.commons.math3.util.Pair< Double, MultivariateNormalDistribution > > components = 
				new ArrayList< org.apache.commons.math3.util.Pair< Double, MultivariateNormalDistribution > >(weightedSurface.getSurfaceSize());
		
		//Computation of the co-variance matrix of (T,K) for the surface
		double[][] coVarianceMatrix = twoDimensionalCoVarianceMatrix(weightedSurface);
		
		//Bandwidth selection by means of Scott and Silverman's rule of thumb
		RealMatrix bandwidthMatrix = new Array2DRowRealMatrix(coVarianceMatrix);
		bandwidthMatrix = bandwidthMatrix.scalarMultiply(
				Math.pow(weightedSurface.getSurfaceSize(), -1.0/3.0) );
		double[][] bandwidth2DArray = bandwidthMatrix.getData();
		
		double mats[] = weightedSurface.getMaturities();
				
		for(int t = 0; t< mats.length; t++) {
			
			double maturity = mats[t];
			
			double[] strikes = weightedSurface.getSmile(maturity).getStrikes();
					
			for(int i = 0; i < strikes.length; i++) {
				
				double strike = strikes[i];
				
				double weight = weightedSurface.getWeight(maturity, strike);
				
				double[] dataCouple = {maturity, strike};
				
				MultivariateNormalDistribution gaussianKernel = new MultivariateNormalDistribution(dataCouple, bandwidth2DArray);
				
				components.add( new org.apache.commons.math3.util.Pair< Double, MultivariateNormalDistribution >(
								weight, gaussianKernel) );
		
			}
			
		}
		
		return new MixtureMultivariateNormalDistribution(components);
		
	}
	
	/**
	 * This private method computes the empirically weighted two-dimensional co-variance matrix of (T,K) (couple maturity - strike),
	 * through the data of the input weighted option surface. The adjective weighted here means that,
	 * empirical means, variances and co-variances are computed through the standard formulas,
	 * while taking into account the weight of every single (T,K) that measures the liquidity on the market of the corresponding option.
	 * 
	 * @param weightedSurface
	 * @return
	 */
	private double[][] twoDimensionalCoVarianceMatrix(WeightedOptionSurfaceData weightedSurface) {
		
		double maturityMean = 0;
		double strikeMean = 0;
		
		double maturityVariance = 0;
		double strikeVariance = 0;
		
		double maturityStrikeCovariance = 0;
		
		double mats[] = weightedSurface.getMaturities();
		
		for(int t = 0; t< mats.length; t++) {
			
			double maturity = mats[t];
					
			double[] strikes = weightedSurface.getSmile(maturity).getStrikes();
					
			for(int i = 0; i < strikes.length; i++) {
				
				double strike = strikes[i];
				
				double weight = weightedSurface.getWeight(maturity, strike);
				
				maturityMean = maturityMean + weight*maturity;
				strikeMean = strikeMean + weight*strike;
				
				maturityVariance = maturityVariance + weight*maturity*maturity;
				strikeVariance = strikeVariance + weight*strike*strike;
				
				maturityStrikeCovariance = maturityStrikeCovariance + weight*maturity*strike;
		
			}
			
		}
		
		maturityVariance = maturityVariance - maturityMean*maturityMean;
		strikeVariance = strikeVariance - strikeMean*strikeMean;
		
		maturityStrikeCovariance = maturityStrikeCovariance - maturityMean*strikeMean; 
		
		return new double[][] {{maturityVariance, maturityStrikeCovariance}, {maturityStrikeCovariance, strikeVariance}};
		
	}
	
	/**
	 * generateDataSmiles() follows McGhee "An Artificial Neural Network Representation of the SABR Stochastic Volatility Model" (2018).
	 *
	 * Similarly to generateDataEachSurface() and generateDataPrices(), it randomly generates as many data sets as there are surfaces over which the model has to be calibrated, 
	 * except that in each case, the smile pricing function is approximated instead of the one-output pricing function or the surface pricing function.
	 * 
	 * This means that for each generated data set, couples are of the form: array of parameter plus T, K1, ..., Kn (input size = numOfParams + n + 1) - array of prices/impl vols (output size = n),
	 * where "n" stands for the size of the returned smile (number of strikes/option quotes in it), which depends on the surface in question.
	 * 
	 * Then, using this approach imposes the restriction that for each surface, every smile/maturity has the same number of option quotes.
	 * If the user lies in presence of surfaces with smiles that do not share the same number of option quotes, it is preferable to use directly generateDataEachSurface() or generateDataPrices() instead.
	 *  
	 * For each data set, the density of (T, K1, ..., Kn) is estimated through weighted (multivariate Gaussian) kernel density estimation (wKDE),
	 * where correlation is assumed between maturities and strikes, meaning that the bandwidth (n+1)x(n+1) matrix reveals to be full.
	 * The latter will be computed by means of either Scott's or Silverman's rule of thumb (they are equivalent only for the bivariate case), where the co-variance matrix,
	 * will be determined empirically through the market data of the surface in question, while taking into account the liquidity weights.
	 * Note that the total liquidity weight of a smile equals the sum of all the weights of its options quotes.
	 *          
	 * Finally, every single data set will be saved inside a file named smiles"name of the surface".txt, 
	 * which will be used for the determination of the neural network approximating the corresponding smile pricing function (one for each surface).
	 * We will then end up with as many neural network as there are surfaces of prices/impl vols.
	 *
     * @param sizeDataSet
     * @throws IOException, IllegalArgumentException
     */
	public void generateDataSmiles(int sizeDataSet) throws IOException, IllegalArgumentException {
		
		//Storing of the names of the surfaces in a list
		//Ordered like LinkedHashMap<String, WeightedOptionSurfaceData> surfaces
		List< String > names = new ArrayList< String >(this.surfaces.size());
		
		//Creation of the data lists that will be used for training the neural networks
		//As many data lists/neural networks as there are surfaces
		//In the same order as LinkedHashMap<String, WeightedOptionSurfaceData> surfaces
		List< List< Pair< double[], double[] > > > dataLists = 
				new ArrayList< List< Pair< double[], double[] > > >(this.surfaces.size());
		
		//Determination of the density of the smile (T, K1, ..., Kn) for every single surface of LinkedHashMap<String, WeightedOptionSurfaceData> surfaces
		//Performed through multivariate weighted kernel density estimation (wKDE) with Gaussian kernel and either Scott's or Silverman's rule of thumb for the selection of the full bandwidth matrix.
		//Again in the same order as LinkedHashMap<String, WeightedOptionSurfaceData> surfaces
		List< MixtureMultivariateNormalDistribution > densities = new ArrayList< MixtureMultivariateNormalDistribution >(this.surfaces.size());
								
		//Iterator pointing at all the surfaces
		Iterator<Map.Entry<String, WeightedOptionSurfaceData>> surfacesIterator = this.surfaces.entrySet().iterator();
		
		//Loop over all surfaces and retrieving of all above data needed for random generation
		while(surfacesIterator.hasNext()) {
			
			Map.Entry<String, WeightedOptionSurfaceData> pair = 
					(Map.Entry<String, WeightedOptionSurfaceData>) surfacesIterator.next();
					
			String underlying = pair.getKey();
			
			//Storing of the name of the surface
			names.add(underlying);
			
			List< Pair< double[], double[] > > dataList = new ArrayList< Pair< double[], double[] > >(sizeDataSet);
			
			//Declaration of the data list for the current surface
			dataLists.add(dataList);
			
			WeightedOptionSurfaceData weightedSurface = pair.getValue();
			
			//Computation of the wKDE density of the smile (T, K1, ..., Kn) for this surface
			/* This also returns an exception in the case that for the input surface,
			 * one of its smiles has a different number of option quotes (smile size) from the other smiles of the surface
			 */
			densities.add( multivariatewKDE(weightedSurface) );
			
		}
		
		//Generation of all data sample couples for all surfaces
		//Every single couple of the form: parameter set + maturity + strike1 + ... + strikeN - price1/impl vol1, ..., priceN/impl volN 
		for(int k = 0; k < sizeDataSet; k++) {
					
			long startMillis = System.currentTimeMillis();
							
			//Retrieving of the k-th couple: parameter set - corresponding MultivariateCalibrableProcessInterface model
			Pair< List<Double>, MultivariateCalibrableProcessInterface > currentPair = this.model.generateSamplePair();
			
			//Retrieving of the model
			MultivariateCalibrableProcessInterface newModel = currentPair.getSecond();
			
			//Loop over the different surfaces
			for(int l = 0; l < this.surfaces.size(); l++) {
				
				//Retrieving of the name
				String surfaceName = names.get(l);
				
				//Retrieving of the surface in question to get:
				//- Quoting convention;
				//- Discount curve;
				//- Forward curve.
				WeightedOptionSurfaceData weightedSurface = this.surfaces.get(surfaceName);
				
				//Retrieving of the sample parameter set of the model
				@SuppressWarnings("unchecked")
				List<Double> newParams = (List<Double>)((ArrayList<Double>)currentPair.getFirst()).clone();
				
				//Generation of a sample smile (T, K1, ..., Kn) for this surface
				double[] sampleSmile = densities.get(l).sample();
				double maturity = Math.abs(sampleSmile[0]);
				double[] strikes = new double[sampleSmile.length-1];
				for(int i = 0; i < strikes.length; i++) 
					strikes[i] = sampleSmile[i+1];
				while(maturity > newModel.getTimeHorizon()) {
					sampleSmile = densities.get(l).sample();
					maturity = Math.abs(sampleSmile[0]);
					for(int i = 0; i < strikes.length; i++) 
						strikes[i] = sampleSmile[i+1];
				}
				
				//Addition of the sample smile (T, K1, ..., Kn) to the sample parameter set of the model
				newParams.add(maturity);
				for(double x : strikes) 
					newParams.add(x);
				
				//Declaration of the pricer 
				EuropeanOptionSmileMultiAsset newPricer = this.pricer.getCloneWithModifiedParameters(
						surfaceName, maturity, strikes);
				
				//The k-th output of the pricing function of the l-th surface
				List<Double> vals = new ArrayList<Double>();
				
				try {
					
					Map<Double, Double> modelPricesImplVols = newPricer.getValue(newModel);
					
					for(int i = 0; i < strikes.length; i++) {
						
						if(weightedSurface.getQuotingConvention().equals(QuotingConvention.VOLATILITYLOGNORMAL)) {
							
							double forward = weightedSurface.getEquityForwardCurve().getDiscountFactor(maturity);
							double optionMaturity = maturity;
							double optionStrike = strikes[i];
							double payoffUnit = weightedSurface.getDiscountCurve().getDiscountFactor(maturity);
							double optionValue = modelPricesImplVols.get(optionStrike);
							double implVol = Math.abs(net.finmath.functions.AnalyticFormulas.blackScholesOptionImpliedVolatility(forward, optionMaturity, optionStrike, payoffUnit, optionValue));
							if(implVol < 10.0) {
								vals.add(implVol);
							} else {
								double distance1 = Math.abs(weightedSurface.getMaturities()[0] - optionMaturity);
								int idx1 = 0;
								for(int c = 1; c < weightedSurface.getMaturities().length; c++) {
									double cdistance = Math.abs(weightedSurface.getMaturities()[c] - optionMaturity);
									if(cdistance < distance1){
										idx1 = c;
										distance1 = cdistance;
									}
								}
								double[] optStrikes = weightedSurface.getSmile(weightedSurface.getMaturities()[idx1]).getStrikes();
								double distance2 = Math.abs(optStrikes[0] - optionStrike);
								int idx2 = 0;
								for(int c = 1; c < optStrikes.length; c++) {
									double cdistance = Math.abs(optStrikes[c] - optionStrike);
									if(cdistance < distance2){
										idx2 = c;
										distance2 = cdistance;
									}
								}
								vals.add(weightedSurface.getSmile(weightedSurface.getMaturities()[idx1]).getOption(optStrikes[idx2]).getValue());	
							}
								
						} else if(weightedSurface.getQuotingConvention().equals(QuotingConvention.VOLATILITYNORMAL)) {
										
							double forward = weightedSurface.getEquityForwardCurve().getDiscountFactor(maturity);
							double optionMaturity = maturity;
							double optionStrike = strikes[i];
							double payoffUnit = weightedSurface.getDiscountCurve().getDiscountFactor(maturity);
							double optionValue = modelPricesImplVols.get(optionStrike);
							double implVol = Math.abs(net.finmath.functions.AnalyticFormulas.bachelierOptionImpliedVolatility(forward, optionMaturity, optionStrike, payoffUnit, optionValue));
							if(implVol < 10.0) {
								vals.add(implVol);
							} else {
								double distance1 = Math.abs(weightedSurface.getMaturities()[0] - optionMaturity);
								int idx1 = 0;
								for(int c = 1; c < weightedSurface.getMaturities().length; c++) {
									double cdistance = Math.abs(weightedSurface.getMaturities()[c] - optionMaturity);
									if(cdistance < distance1){
										idx1 = c;
										distance1 = cdistance;
									}
								}
								double[] optStrikes = weightedSurface.getSmile(weightedSurface.getMaturities()[idx1]).getStrikes();
								double distance2 = Math.abs(optStrikes[0] - optionStrike);
								int idx2 = 0;
								for(int c = 1; c < optStrikes.length; c++) {
									double cdistance = Math.abs(optStrikes[c] - optionStrike);
									if(cdistance < distance2){
										idx2 = c;
										distance2 = cdistance;
									}
								}
								vals.add(weightedSurface.getSmile(weightedSurface.getMaturities()[idx1]).getOption(optStrikes[idx2]).getValue());	
							}
										
						} else {
								
							vals.add(modelPricesImplVols.get(strikes[i]));
										
						}	
						
					}					
							
				} catch (CalculationException e) {
					e.printStackTrace();
				}
				
				double[] params = ArrayUtils.toPrimitive( newParams.toArray(new Double[newParams.size()]) );
				double[] values = ArrayUtils.toPrimitive( vals.toArray(new Double[vals.size()]) );;
				
				dataLists.get(l).add( Pair.< double[], double[] >create(params, values) );
				
			}
							
			long endMillis = System.currentTimeMillis();
			
			double calculationTime = ((endMillis-startMillis)/1000.0);
			
			System.out.println(k + "-th data sample couples for all surfaces retrieved in: " + calculationTime + " seconds");
			
		}
		
		//Loop over the different surfaces
		for(int p = 0; p < this.surfaces.size(); p++) {
			
			//Creation of the data set 
			DataSetIterator iterator = new DoublesDataSetIterator(dataLists.get(p), sizeDataSet);
			DataSet dataSet = iterator.next();
			
			//Creation of the file in which the data set will be saved
			File dataFile = new File("smiles" + names.get(p) + ".txt");
			
			//Check whether data.txt already exists or not
			if(dataFile.createNewFile()) {
				dataSet.save(dataFile);
			} else {
				dataFile.delete();
				dataFile.createNewFile();
				dataSet.save(dataFile);
			}
			
		}
		
	}
	
	/**
	 * This private method returns the multivariate weighted Kernel Density Estimation (wKDE) of the density of the smile (T, K1, ..., Kn) for the input weighted surface.
	 * The adjective weighted here means that the standard KDE formula is computed by taking into account the weight of every smile (T, K1, ..., Kn),
	 * which measures the liquidity on the market of the corresponding smile (the smile weight equals the sum of the weights of all its option quotes).
	 * The latter is observable in the formula where it replaces (1/n) for every single summand.
	 * Finally, Gaussian Kernel will be used for computations along with either Scott's or Silverman's rule of thumb for bandwidth selection (they are equal only for bivariate kernel density estimation).
	 * 
	 * Note also that an exception will be thrown if in the input surface, all the smiles do not share the same number of options quotes (smile size).
	 * 
	 * @param weightedSurface
	 * @throws IllegalArgumentException
	 * @return
	 */
	private MixtureMultivariateNormalDistribution multivariatewKDE(WeightedOptionSurfaceData weightedSurface) throws IllegalArgumentException {
		
		//Preparation of the inputs required for the creation of the wKDE density of the smile (T, K1, ..., Kn) for this surface
		List< org.apache.commons.math3.util.Pair< Double, MultivariateNormalDistribution > > components = 
				new ArrayList< org.apache.commons.math3.util.Pair< Double, MultivariateNormalDistribution > >(weightedSurface.getNumberOfSmiles());
		
		//Computation of the multi-dimnesional co-variance matrix of the smile (T, K1, ..., Kn) for the surface (size (n+1)x(n+1)):
		/* If one of the smiles of the surface has a different size (number of option quotes) from the other smiles,
		 * an exception will be thrown at this point.
		 */
		RealMatrix coVarianceMatrix = multiDimensionalCoVarianceMatrix(weightedSurface);
		int totalSmileSize = coVarianceMatrix.getRowDimension(); //Size of (T, K1, ..., Kn), which is n + 1
		
		//Bandwidth selection by means of Scott's rule of thumb (more stable than Silverman's as totalSmileSize rises)
		RealMatrix bandwidthMatrix = coVarianceMatrix.scalarMultiply(
				Math.pow( weightedSurface.getNumberOfSmiles(), -2.0/(totalSmileSize+4.0) ) );
		double[][] bandwidth2DArray = bandwidthMatrix.getData();
		
		double mats[] = weightedSurface.getMaturities();
				
		for(int t = 0; t< mats.length; t++) {
			
			double[] maturity = {mats[t]};
			
			double[] strikes = weightedSurface.getSmile(maturity[0]).getStrikes();
			
			double[] smile = ArrayUtils.addAll(maturity, strikes);
					
			MultivariateNormalDistribution gaussianKernel = new MultivariateNormalDistribution(smile, bandwidth2DArray);
				
			components.add( new org.apache.commons.math3.util.Pair< Double, MultivariateNormalDistribution >(
					weightedSurface.getWeightedSmile(maturity[0]).getSmileWeight(), gaussianKernel) );
		
		}
		
		return new MixtureMultivariateNormalDistribution(components);
		
	}
	
	/**
	 * This private method computes the empirically weighted multi-dimensional co-variance matrix of the smile (T, K1, ..., Kn) (matrix size (n+1)x(n+1)),
	 * through the data of the input weighted option surface. The adjective weighted here means that,
	 * empirical means, variances and co-variances are computed through the standard formulas,
	 * while taking into account the weight of every single smile (T, K1, ..., Kn) that measures its market liquidity.
	 * Note that for each smile, the weight equals the sum of all the weights of the option quotes that it contains.
	 * Finally, an exception will be thrown by this method if one of the smiles of the input surface has a size (number of quotes),
	 * which is different from the other smiles.
	 * 
	 * @param weightedSurface
	 * @throws IllegalArgumentException
	 * @return
	 */
	private RealMatrix multiDimensionalCoVarianceMatrix(WeightedOptionSurfaceData weightedSurface) throws IllegalArgumentException {
		
		double[] maturities = weightedSurface.getMaturities();
		int smileSize = weightedSurface.getSmileSize(maturities[0]);
		for(int i = 1; i < maturities.length; i++) {
			if(smileSize!=weightedSurface.getSmileSize(maturities[i])) 
				throw new IllegalArgumentException("all smiles of the surface must share the same number of option quotes (same size)");
		}
		
		double maturityMean = 0;
		double maturityVariance = 0;
		double[] strikesMean = new double[smileSize];
		double[] maturityStrikesCovariance = new double[smileSize];
		double[][] strikesCovariance = new double[smileSize][smileSize];
		
		for(int i = 0; i < smileSize; i++) {
			strikesMean[i] = 0;
			maturityStrikesCovariance[i] = 0;
			for(int j = 0; j < smileSize; j++)
				strikesCovariance[i][j] = 0;
		}
		
		double mats[] = weightedSurface.getMaturities();
		
		for(int t = 0; t < mats.length; t++) {
			
			double maturity = mats[t];
			
			double weight = weightedSurface.getWeightedSmile(maturity).getSmileWeight();
			
			maturityMean = maturityMean + weight*maturity;
			
			maturityVariance = maturityVariance + weight*maturity*maturity;
					
			double[] strikes = weightedSurface.getSmile(maturity).getStrikes();
					
			for(int i = 0; i < strikes.length; i++) {
				
				strikesMean[i] = strikesMean[i] + weight*strikes[i];
				
				maturityStrikesCovariance[i] = maturityStrikesCovariance[i] + weight*maturity*strikes[i];
				
				for(int j = 0; j < strikes.length; j++) 
					strikesCovariance[i][j] = strikesCovariance[i][j] + weight*strikes[i]*strikes[j];
		
			}
			
		}
		
		maturityVariance = maturityVariance - maturityMean*maturityMean;
		for(int i = 0; i < smileSize; i++) {
			maturityStrikesCovariance[i] = maturityStrikesCovariance[i] - strikesMean[i]*maturityMean;
			for(int j = 0; j < smileSize; j++)
				strikesCovariance[i][j] = strikesCovariance[i][j] - strikesMean[i]*strikesMean[j];
		}
		
		double[][] coVarianceMatrix = new double[smileSize+1][smileSize+1];
		coVarianceMatrix[0][0] = maturityVariance;
		for(int i = 1; i < smileSize+1; i++) {
			coVarianceMatrix[0][i] = maturityStrikesCovariance[i-1];
			coVarianceMatrix[i][0] = maturityStrikesCovariance[i-1];
			for(int j = 1; j < smileSize+1; j++)
				coVarianceMatrix[i][j] = strikesCovariance[i-1][j-1];
		}
		
		return new Array2DRowRealMatrix(coVarianceMatrix);
		
	}

}
