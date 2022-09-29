package org.calibrationframework.fouriermethod.calibration.models;

import org.calibrationframework.stochastic.*;
import org.calibrationframework.fouriermethod.CharacteristicFunctionInterface;
import org.calibrationframework.fouriermethod.calibration.constraints.*;
import org.calibrationframework.fouriermethod.models.MultivariateProcessCharacteristicFunctionInterface;

import java.util.*;
import java.util.function.UnaryOperator;

import org.apache.commons.lang3.ArrayUtils;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.random.*;

import org.nd4j.common.primitives.Pair;

/**
 * This class represents the multi-currency CBITCL modeling framework where each CBICTL process is taken to be tempered stable.
 * 
 * The setting is directly reduced to three currencies (JPY-1, USD-2, EUR-3) in view of a calibration of the model to such a currency triangle.
 * 
 * The main feature of the class is given by both methods apply and applyForTwoStrings. For a couple (i,j) denoting the currency pair, 
 * they provide the characteristic function of $\log S^{i,j}$ under $\QQ_i$, thus paving the way to pricing and calibration.
 * 
 * @author Szulda Guillaume 
 */
public class TemperedStableCBITCLMultiCurrencyModel implements MultivariateCalibrableProcessInterface {
	
	private final double minP; //Defines the extended domain of the characteristic function of the model.
	
	private final double timeHorizon;
	private final int numberOfTimeSteps;
	
	private final int dimension; //Size of temperedStableCBITCLs.
	
	private TemperedStableCBITCLProcess[] temperedStableCBITCLs; //Tempered stable CBITCL processes under $\QQ_0$.
	
	private TemperedStableCBITCLProcess[] temperedStableCBITCLsJPY; //Tempered stable CBITCL processes under $\QQ_{JPY}$.
	private TemperedStableCBITCLProcess[] temperedStableCBITCLsUSD; //Tempered stable CBITCL processes under $\QQ_{USD}$.
	private TemperedStableCBITCLProcess[] temperedStableCBITCLsEUR; //Tempered stable CBITCL processes under $\QQ_{EUR}$.
	
	private double interestRateJPY;
	private double interestRateUSD;
	private double interestRateEUR;
	
	private double initialValueJPYtoUSD;
	private double initialValueJPYtoEUR;
	private double initialValueUSDtoEUR;
	
	private double[] zetasJPY;
	private double[] zetasUSD;
	private double[] zetasEUR;
	
	private double[] lambdasJPY;
	private double[] lambdasUSD;
	private double[] lambdasEUR;
	
	private ScalarParameterInformationInterface[] zetasJPYInfo;
	private ScalarParameterInformationInterface[] zetasUSDInfo;
	private ScalarParameterInformationInterface[] zetasEURInfo;
	
	private ScalarParameterInformationInterface[] lambdasJPYInfo;
	private ScalarParameterInformationInterface[] lambdasUSDInfo;
	private ScalarParameterInformationInterface[] lambdasEURInfo;
	
	/*
	 * Upper and lower bounds are collected here for convenience:
	 * such vectors are then passed to the factory of the optimization algorithm.
	 * In this way we guarantee consistency between the constraints in the model
	 * and the constraints in the optimizer factory.
	 */
	private double[] parameterLowerBounds;
	private double[] parameterUpperBounds;
	
	/**
	 * First constructor creates a textbook specification of the multi-currency tempered stable CBITCL modeling framework.
	 * 
	 * @param temperedStableCBITCLs
	 * @param interestRateJPY
	 * @param interestRateUSD
	 * @param interestRateEUR
	 * @param initialValueJPYtoUSD
	 * @param initialValueJPYtoEUR
	 * @param initialValueUSDtoEUR
	 * @param zetasJPY
	 * @param zetasUSD
	 * @param zetasEUR
	 * @param lambdasJPY
	 * @param lambdasUSD
	 * @param lambdasEUR
	 * @throws IllegalArgumentException
	 */
	public TemperedStableCBITCLMultiCurrencyModel(TemperedStableCBITCLProcess[] temperedStableCBITCLs, 
			double interestRateJPY, double interestRateUSD, double interestRateEUR, 
			double initialValueJPYtoUSD, double initialValueJPYtoEUR, double initialValueUSDtoEUR,
			double[] zetasJPY, double[] zetasUSD, double[] zetasEUR, 
			double[] lambdasJPY, double[] lambdasUSD, double[] lambdasEUR ) throws IllegalArgumentException {
		
		if(temperedStableCBITCLs.length != zetasJPY.length || temperedStableCBITCLs.length != zetasUSD.length || temperedStableCBITCLs.length != zetasEUR.length ||
				temperedStableCBITCLs.length != lambdasJPY.length || temperedStableCBITCLs.length != lambdasUSD.length || temperedStableCBITCLs.length != lambdasEUR.length) {
			
			throw new IllegalArgumentException("temperedStableCBITCLs, zetasJPY/USD/EUR and lambdasJPY/USD/EUR must all have the same size");
		
		} else {
			
			this.dimension = temperedStableCBITCLs.length;
			
			for(int i = 0; i < this.dimension; i++) {
				
				if(lambdasJPY[i] <= 0.0 || lambdasUSD[i] <= 0.0 || lambdasEUR[i] <= 0.0 || 
						lambdasJPY[i] >= 1.0 || lambdasUSD[i] >= 1.0 || lambdasEUR[i] >= 1.0 ||	
						zetasJPY[i] <= 0.0 || zetasUSD[i] <= 0.0 || zetasEUR[i] <= 0.0 || 
						zetasJPY[i] >= 1.0 || zetasUSD[i] >= 1.0 || zetasEUR[i] >= 1.0) {
					
					throw new IllegalArgumentException("All lambdasJPY/USD/EUR and zetasJPY/USD/EUR must be included in (0, 1)");
					
				}
				
			}
			
			double timeHorizon = temperedStableCBITCLs[0].getTimeHorizon();
			int numberOfTimeSteps = temperedStableCBITCLs[0].getNumberOfTimeSteps();
			
			for(int i = 1; i < this.dimension; i++) {
				
				if(timeHorizon != temperedStableCBITCLs[i].getTimeHorizon() || numberOfTimeSteps != temperedStableCBITCLs[i].getNumberOfTimeSteps()) {
					
					throw new IllegalArgumentException("All tempered stable CBITCL processes must have the same time horizon and "
							+ "the same number of time steps");
				
				} 
				
			}
				
			this.timeHorizon = timeHorizon;
			this.numberOfTimeSteps = numberOfTimeSteps;
			
			this.temperedStableCBITCLs = temperedStableCBITCLs;
			
			this.interestRateJPY = interestRateJPY;
			this.interestRateUSD = interestRateUSD;
			this.interestRateEUR = interestRateEUR;
				
			this.initialValueJPYtoUSD = initialValueJPYtoUSD;
			this.initialValueJPYtoEUR = initialValueJPYtoEUR;
			this.initialValueUSDtoEUR = initialValueUSDtoEUR;
				
			this.zetasJPY = zetasJPY;
			this.zetasUSD = zetasUSD;
			this.zetasEUR = zetasEUR;
				
			this.lambdasJPY = lambdasJPY;
			this.lambdasUSD = lambdasUSD;
			this.lambdasEUR = lambdasEUR;
				
			this.zetasJPYInfo = new ScalarParameterInformationInterface[this.dimension]; 
			this.zetasUSDInfo = new ScalarParameterInformationInterface[this.dimension]; 
			this.zetasEURInfo = new ScalarParameterInformationInterface[this.dimension]; 
			
			this.lambdasJPYInfo = new ScalarParameterInformationInterface[this.dimension]; 
			this.lambdasUSDInfo = new ScalarParameterInformationInterface[this.dimension]; 
			this.lambdasEURInfo = new ScalarParameterInformationInterface[this.dimension]; 
				
			/* The next one will determine the information on the parameters that was predetermined for this specification of the model that will be used for calibration,
			 * which consists of the constraints that make the model admissible and whether one parameter has to be calibrated or not. 
			 */
			generateInformationForParameters();
			this.parameterLowerBounds = extractLowerBounds();
			this.parameterUpperBounds = extractUpperBounds();
			
			/* The following loop forces the martingale measure assumption of the multi-currency tempered stable CBITCL modeling framework to be satisfied.
			 * It also includes the model specification using tempered stable CBITCL processes of CGMY type.
			 */
			for(int k = 0; k < this.dimension; k++) 
				martingaleMeasureAssumption(k);
			
			this.temperedStableCBITCLsJPY = new TemperedStableCBITCLProcess[this.dimension];
			this.temperedStableCBITCLsUSD = new TemperedStableCBITCLProcess[this.dimension];
			this.temperedStableCBITCLsEUR = new TemperedStableCBITCLProcess[this.dimension];
			
			/* This generates the $\QQ_{JPY}$,$\QQ_{USD}$ and $\QQ_{EUR}$-versions of the tempered stable CBITCL processes,
			 * also including tempered stable CBITCL processes of CGMY type. This is permitted thanks to 
			 * the validity of the martingale measure assumption of the multi-currency tempered stable CBITCL modeling framework.
			 */
			generateRiskNeutralTemperedStableCBITCLs();
			
			//Determination of the extended domain of the characteristic function of the model.
			this.minP = computationOfMinP();
			
		}
		
	}
	
	/**
	 * Second constructor provides a more general specification of the multi-currency tempered stable CBITCL modeling framework.
	 * In particular, we are allowed to have more general information on parameters.
	 * 
	 * @param temperedStableCBITCLs
	 * @param interestRateJPY
	 * @param interestRateUSD
	 * @param interestRateEUR
	 * @param initialValueJPYtoUSD
	 * @param initialValueJPYtoEUR
	 * @param initialValueUSDtoEUR
	 * @param zetasJPY
	 * @param zetasJPYInfo
	 * @param zetasUSD
	 * @param zetasUSDInfo
	 * @param zetasEUR
	 * @param zetasEURInfo
	 * @param lambdasJPY
	 * @param lambdasJPYInfo
	 * @param lambdasUSD
	 * @param lambdasUSDInfo
	 * @param lambdasEUR
	 * @param lambdasEURInfo
	 * @throws IllegalArgumentException
	 */
	public TemperedStableCBITCLMultiCurrencyModel(TemperedStableCBITCLProcess[] temperedStableCBITCLs, double interestRateJPY, double interestRateUSD, 
			double interestRateEUR, double initialValueJPYtoUSD, double initialValueJPYtoEUR, double initialValueUSDtoEUR,
			double[] zetasJPY, ScalarParameterInformationInterface[] zetasJPYInfo,  double[] zetasUSD, ScalarParameterInformationInterface[] zetasUSDInfo, 
			double[] zetasEUR, ScalarParameterInformationInterface[] zetasEURInfo, double[] lambdasJPY, ScalarParameterInformationInterface[] lambdasJPYInfo, 
			double[] lambdasUSD, ScalarParameterInformationInterface[] lambdasUSDInfo,  
			double[] lambdasEUR, ScalarParameterInformationInterface[] lambdasEURInfo) throws IllegalArgumentException {
		
		if(temperedStableCBITCLs.length != zetasJPY.length || temperedStableCBITCLs.length != zetasUSD.length || temperedStableCBITCLs.length != zetasEUR.length ||
				temperedStableCBITCLs.length != lambdasJPY.length || temperedStableCBITCLs.length != zetasUSD.length || temperedStableCBITCLs.length != lambdasEUR.length ||
				temperedStableCBITCLs.length != zetasJPYInfo.length || temperedStableCBITCLs.length != zetasUSDInfo.length || temperedStableCBITCLs.length != zetasEURInfo.length ||
				temperedStableCBITCLs.length != lambdasJPYInfo.length || temperedStableCBITCLs.length != lambdasUSDInfo.length || temperedStableCBITCLs.length != lambdasEURInfo.length) {
			
			throw new IllegalArgumentException("temperedStableCBITCLs, zetasJPY/USD/EUR, zetasInfoJPY/USD/EUR, lambdasJPY/USD/EUR and lambdasInfoJPY/USD/EUR must all have the same size");
		
		} else {
			
			this.dimension = temperedStableCBITCLs.length;
			
			for(int i = 0; i < this.dimension; i++) {
				
				if(lambdasJPY[i] <= 0.0 || lambdasUSD[i] <= 0.0 || lambdasEUR[i] <= 0.0 || 
						lambdasJPY[i] >= 1.0 || lambdasUSD[i] >= 1.0 || lambdasEUR[i] >= 1.0 ||	
						zetasJPY[i] <= 0.0 || zetasUSD[i] <= 0.0 || zetasEUR[i] <= 0.0 || 
						zetasJPY[i] >= 1.0 || zetasUSD[i] >= 1.0 || zetasEUR[i] >= 1.0) {
					
					throw new IllegalArgumentException("All lambdasJPY/USD/EUR and zetasJPY/USD/EUR must be included in (0, 1)");
					
				}
				
			}
			
			double timeHorizon = temperedStableCBITCLs[0].getTimeHorizon();
			int numberOfTimeSteps = temperedStableCBITCLs[0].getNumberOfTimeSteps();
			
			for(int i = 1; i < this.dimension; i++) {
				
				if(timeHorizon != temperedStableCBITCLs[i].getTimeHorizon() || numberOfTimeSteps != temperedStableCBITCLs[i].getNumberOfTimeSteps()) {
					
					throw new IllegalArgumentException("All tempered stable CBITCL processes must have the same time horizon "
							+ "and the same number of time steps");
				
				} 
				
			}
				
			this.timeHorizon = timeHorizon;
			this.numberOfTimeSteps = numberOfTimeSteps;
			
			this.temperedStableCBITCLs = temperedStableCBITCLs;
			
			this.interestRateJPY = interestRateJPY;
			this.interestRateUSD = interestRateUSD;
			this.interestRateEUR = interestRateEUR;
				
			this.initialValueJPYtoUSD = initialValueJPYtoUSD;
			this.initialValueJPYtoEUR = initialValueJPYtoEUR;
			this.initialValueUSDtoEUR = initialValueUSDtoEUR;
				
			this.zetasJPY = zetasJPY;
			this.zetasUSD = zetasUSD;
			this.zetasEUR = zetasEUR;
				
			this.lambdasJPY = lambdasJPY;
			this.lambdasUSD = lambdasUSD;
			this.lambdasEUR = lambdasEUR;
				
			this.zetasJPYInfo = zetasJPYInfo;
			this.zetasUSDInfo = zetasUSDInfo;
			this.zetasEURInfo = zetasEURInfo;
			
			this.lambdasJPYInfo = lambdasJPYInfo;
			this.lambdasUSDInfo = lambdasUSDInfo;
			this.lambdasEURInfo = lambdasEURInfo;
				
			this.parameterLowerBounds = extractLowerBounds();
			this.parameterUpperBounds = extractUpperBounds();
			
			/* The following loop forces the martingale measure assumption of the multi-currency tempered stable CBITCL modeling framework to be satisfied.
			 * It also includes the model specification using tempered stable CBITCL processes of CGMY type.
			 */
			for(int k = 0; k < this.dimension; k++) 
				martingaleMeasureAssumption(k);
			
			this.temperedStableCBITCLsJPY = new TemperedStableCBITCLProcess[this.dimension];
			this.temperedStableCBITCLsUSD = new TemperedStableCBITCLProcess[this.dimension];
			this.temperedStableCBITCLsEUR = new TemperedStableCBITCLProcess[this.dimension];
			
			/* This generates the $\QQ_{JPY}$,$\QQ_{USD}$ and $\QQ_{EUR}$-versions of the tempered stable CBITCL processes,
			 * also including tempered stable CBITCL processes of CGMY type. This is permitted thanks to 
			 * the validity of the martingale measure assumption of the multi-currency tempered stable CBITCL modeling framework.
			 */
			generateRiskNeutralTemperedStableCBITCLs();
			
			//Determination of the extended domain of the characteristic function of the model.
			this.minP = computationOfMinP();
			
		}
		
	}
	
	public double getMinP() {
		return this.minP;
	}
	
	@Override
	public double getTimeHorizon() {
		return this.timeHorizon;
	}
	
	@Override
	public int getNumberOfTimeSteps() {
		return this.numberOfTimeSteps;
	}
	
	public int getNumberOfCurrencies() {
		return 3;
	}
	
	public int getDimension() {
		return this.dimension;
	}
	
	public int getCurrencyToIndex(String underlyingName) throws IllegalArgumentException {
		if(underlyingName.contains("JPY")) {
			return 1;
		} else if(underlyingName.contains("USD")) {
			return 2;
		} else if(underlyingName.contains("EUR")) {
			return 3;
		} else {
			throw new IllegalArgumentException("Currency not recognized");
		}
	}
	
	public double getInitialValue(String underlyingName) {
		String underlyingNameForeign = underlyingName.subSequence(0,3).toString();
		String underlyingNameDomestic = underlyingName.subSequence(3,6).toString();
		return getInitialValue(underlyingNameDomestic, underlyingNameForeign);
	}
	
	public double getInitialValue(String underlyingNameDomestic, String underlyingNameForeign) throws IllegalArgumentException {
		return this.getInitialValue(this.getCurrencyToIndex(underlyingNameDomestic), this.getCurrencyToIndex(underlyingNameForeign));
	}
	
	public double getInitialValue(int i, int j) throws IllegalArgumentException  {
		if(i == 1) {
			if(j == 2) {
				return this.initialValueJPYtoUSD;
			} else if(j == 3) {
				return this.initialValueJPYtoEUR;
			} else {
				throw new IllegalArgumentException("Index non valid");
			}
		} else if(i == 2) {
			if(j == 3) {
				return this.initialValueUSDtoEUR;
			} else {
				throw new IllegalArgumentException("Index non valid");
			}
		} else {
			throw new IllegalArgumentException("Index non valid");
		}
	}
	
	public double getZeta(String underlyingName, int k) throws IllegalArgumentException {
		return this.getZeta(this.getCurrencyToIndex(underlyingName), k);
	}
	
	public double getZeta(int i, int k) throws IllegalArgumentException {
		if(i == 1) {
			return zetasJPY[k];
		} else if(i == 2) {
			return zetasUSD[k];
		} else if(i == 3) {
			return zetasEUR[k];
		} else {
			throw new IllegalArgumentException("Index non valid");
		}
	}
	
	public double getLambda(String underlyingName, int k) throws IllegalArgumentException {
		return this.getLambda(this.getCurrencyToIndex(underlyingName), k);
	}
	
	public double getLambda(int i, int k) throws IllegalArgumentException {
		if(i == 1) {
			return lambdasJPY[k];
		} else if(i == 2) {
			return lambdasUSD[k];
		} else if(i == 3) {
			return lambdasEUR[k];
		} else {
			throw new IllegalArgumentException("Index non valid");
		}
	}
	
	public double getInterestRate(String underlyingName) throws IllegalArgumentException {
		return this.getInterestRate(this.getCurrencyToIndex(underlyingName));
	}
	
	public double getInterestRate(int i) throws IllegalArgumentException {
		if(i == 1) {
			return interestRateJPY;
		} else if(i == 2) {
			return interestRateUSD;
		} else if(i == 3) {
			return interestRateEUR;
		} else {
			throw new IllegalArgumentException("Index non valid");
		}
	}
	
	@Override
	public MultivariateProcessCharacteristicFunctionInterface getCharacteristiFunction() {
		return this;
	}
	
	@Override
	public CharacteristicFunctionInterface apply(double time, String underlyingName) throws IllegalArgumentException  {
		String underlyingNameForeign = underlyingName.subSequence(0,3).toString();
		String underlyingNameDomestic = underlyingName.subSequence(3,6).toString();
		return applyForTwoStrings(time, underlyingNameDomestic, underlyingNameForeign);
	}
	
	/**
	 * This method represents the main feature of the present class.
	 * 
	 * For a couple (i,j) denoting the currency pair,
	 * where int i = this.getCurrencyToIndex(underlyingNameDomestic),
	 * and int j = this.getCurrencyToIndex(underlyingNameForeign),
	 * it provides the characteristic function of $\log S^{i,j}$ under $\QQ_i$.
	 * 
	 * The latter can be computed for all $w \in \C_+$ such that $0 \leq \Re(w) < \min P$, 
	 * where $\min P$ stands for the minimal value of all the P's of the tempered stable CBITCL processes (of CGMY type).
	 * 
	 * @param time
	 * @param underlyingNameDomestic
	 * @param underlyingNameForeign
	 * @throws IllegalArgumentException
	 */
	public CharacteristicFunctionInterface applyForTwoStrings(double time, String underlyingNameDomestic, String underlyingNameForeign) 
			throws IllegalArgumentException  {
		
		if(time <= this.getTimeHorizon()) {
			
			int i = this.getCurrencyToIndex(underlyingNameDomestic);
			int j = this.getCurrencyToIndex(underlyingNameForeign);
			
			if(i == 1) {
				if(j == 2) {
					
					return new CharacteristicFunctionInterface() {

						@Override
						public Complex apply(Complex w) {
							
							if(w.getReal() >= 0.0 && w.getReal() < minP) {
							
							int d = getDimension();
							
							Complex m = new Complex(1.0);
							Complex n = new Complex(1.0);
							for(int k = 0; k < d; k++) {
								
								UnaryOperator<Complex> branchingMechanism = temperedStableCBITCLs[k].getBranchingMechanism();
								UnaryOperator<Complex> immigrationRate = temperedStableCBITCLs[k].getImmigrationRate();
								UnaryOperator<Complex> levyExponent = temperedStableCBITCLs[k].getLevyExponent();
								
								Complex mk = ( w.multiply(time).multiply( 
										immigrationRate.apply(new Complex(zetasJPY[k])).subtract( immigrationRate.apply(new Complex(zetasUSD[k])) ) ) ).exp();
								
								m = m.multiply(mk);
								
								Complex u1k = w.multiply( zetasUSD[k] - zetasJPY[k] );
								Complex u2k = w.multiply( branchingMechanism.apply(new Complex(zetasJPY[k])).add(levyExponent.apply(new Complex(lambdasJPY[k])))
										.subtract(branchingMechanism.apply(new Complex(zetasUSD[k]))).subtract(levyExponent.apply(new Complex(lambdasUSD[k]))) );
								Complex u3k = w.multiply( lambdasUSD[k] - lambdasJPY[k] );
					
								Complex nk = temperedStableCBITCLsJPY[k].getLaplaceFourierTransform(time, u1k, u2k, u3k);
								
								n = n.multiply(nk);
								
							}
							
							return m.multiply(n).multiply( ( w.multiply( Math.log(initialValueJPYtoUSD) + time*(interestRateJPY - interestRateUSD) ) ).exp() );
							
							} else {
							
							throw new IllegalArgumentException("The real part of the argument must be included in [0, minP)");
							
							}
							
						}
						
					};
					
				} else if(j == 3) {
					
					return new CharacteristicFunctionInterface() {

						@Override
						public Complex apply(Complex w) {
							
							if(w.getReal() >= 0.0 && w.getReal() < minP) {
							
							int d = getDimension();
							
							Complex m = new Complex(1.0);
							Complex n = new Complex(1.0);
							for(int k = 0; k < d; k++) {
								
								UnaryOperator<Complex> branchingMechanism = temperedStableCBITCLs[k].getBranchingMechanism();
								UnaryOperator<Complex> immigrationRate = temperedStableCBITCLs[k].getImmigrationRate();
								UnaryOperator<Complex> levyExponent = temperedStableCBITCLs[k].getLevyExponent();
								
								Complex mk = ( w.multiply(time).multiply( 
										immigrationRate.apply(new Complex(zetasJPY[k])).subtract( immigrationRate.apply(new Complex(zetasEUR[k])) ) ) ).exp();
								
								m = m.multiply(mk);
								
								Complex u1k = w.multiply( zetasEUR[k] - zetasJPY[k] );
								Complex u2k = w.multiply( branchingMechanism.apply(new Complex(zetasJPY[k])).add(levyExponent.apply(new Complex(lambdasJPY[k])))
										.subtract(branchingMechanism.apply(new Complex(zetasEUR[k]))).subtract(levyExponent.apply(new Complex(lambdasEUR[k]))) );
								Complex u3k = w.multiply( lambdasEUR[k] - lambdasJPY[k] );
								
								Complex nk = temperedStableCBITCLsJPY[k].getLaplaceFourierTransform(time, u1k, u2k, u3k);
								
								n = n.multiply(nk);
								
							}
							
							return m.multiply(n).multiply( ( w.multiply( Math.log(initialValueJPYtoEUR) + time*(interestRateJPY - interestRateEUR) ) ).exp() );
							
							} else {
								
								throw new IllegalArgumentException("The real part of the argument must be included in [0, minP)");
								
							}
							
						}
						
					};
					
				} else {
					throw new IllegalArgumentException("Index non valid");
				}
			} else if(i == 2) {
				if(j == 3) {
					
					return new CharacteristicFunctionInterface() {

						@Override
						public Complex apply(Complex w) {
							
							if(w.getReal() >= 0.0 && w.getReal() < minP) {
							
							int d = getDimension();
							
							Complex m = new Complex(1.0);
							Complex n = new Complex(1.0);
							for(int k = 0; k < d; k++) {
								
								UnaryOperator<Complex> branchingMechanism = temperedStableCBITCLs[k].getBranchingMechanism();
								UnaryOperator<Complex> immigrationRate = temperedStableCBITCLs[k].getImmigrationRate();
								UnaryOperator<Complex> levyExponent = temperedStableCBITCLs[k].getLevyExponent();
								
								Complex mk = ( w.multiply(time).multiply( 
										immigrationRate.apply(new Complex(zetasUSD[k])).subtract( immigrationRate.apply(new Complex(zetasEUR[k])) ) ) ).exp();
								
								m = m.multiply(mk);
								
								Complex u1k = w.multiply( zetasEUR[k] - zetasUSD[k] );
							
								Complex u2k = w.multiply( branchingMechanism.apply(new Complex(zetasUSD[k])).add(levyExponent.apply(new Complex(lambdasUSD[k])))
										.subtract(branchingMechanism.apply(new Complex(zetasEUR[k]))).subtract(levyExponent.apply(new Complex(lambdasEUR[k]))) );
								
								Complex u3k = w.multiply( lambdasEUR[k] - lambdasUSD[k] );
						
								Complex nk = temperedStableCBITCLsUSD[k].getLaplaceFourierTransform(time, u1k, u2k, u3k);
								
								n = n.multiply(nk);
								
							}
							
							return m.multiply(n).multiply( ( w.multiply( Math.log(initialValueUSDtoEUR) + time*(interestRateUSD - interestRateEUR) ) ).exp() );
							
							} else {
								
								throw new IllegalArgumentException("The real part of the argument must be included in [0, minP)");
								
							}
							
						}
						
					};
					
				} else {
					throw new IllegalArgumentException("Index non valid");
				}
			} else {
				throw new IllegalArgumentException("Index non valid");
			}
			
		} else {
			
			throw new IllegalArgumentException("The time at which the characteristic function is considered must be lower than the time horizon");
			
		}
		
	}

	@Override
	public TemperedStableCBITCLMultiCurrencyModel getCloneForModifiedParameters(double[] parameters) {
		
		/* For each parameter, we check whether it has to be calibrated or not. If so, we replace it for the new one to which we apply the corresponding constraint.
		 * If not, this parameter is not modified.
		 */
		int d = this.getDimension();
		TemperedStableCBITCLProcess[] newTemperedStableCBITCLs = new TemperedStableCBITCLProcess[d];
		int p = 0;
		
		for(int k = 0; k < d; k++) {
			
			if(this.temperedStableCBITCLs[k] instanceof TemperedStableCBITCBrownian) {
				
				p = this.temperedStableCBITCLs[k].getNumberOfParameters();
				double[] params = new double[p];
				for(int l = 0; l < p; l++) {
					params[l] = parameters[k*p + l];
				}
				newTemperedStableCBITCLs[k] = ((TemperedStableCBITCBrownian)this.temperedStableCBITCLs[k]).getCloneForModifiedParameters(params);
			
			} else if(this.temperedStableCBITCLs[k] instanceof TemperedStableCBITCLofCGMYtype) {
				
				p = ((TemperedStableCBITCLofCGMYtype)this.temperedStableCBITCLs[k]).getNumberOfParameters();
				double[] params = new double[p];
				for(int l = 0; l < p; l++) {
					params[l] = parameters[k*p + l];
				}
				newTemperedStableCBITCLs[k] = ((TemperedStableCBITCLofCGMYtype)this.temperedStableCBITCLs[k]).getCloneForModifiedParameters(params);
				
			}
			
		}
		
		double[] newZetasJPY = new double[d];
		double[] newZetasUSD = new double[d];
		double[] newZetasEUR = new double[d];
		
		for(int k = 0; k < d; k++) {
			
			newZetasJPY[k] = this.zetasJPYInfo[k].getIsParameterToCalibrate() == true ? this.zetasJPYInfo[k].getConstraint().applyConstraint(parameters[ 3*k + p*d ]) : this.zetasJPY[k];
			newZetasUSD[k] = this.zetasUSDInfo[k].getIsParameterToCalibrate() == true ? this.zetasUSDInfo[k].getConstraint().applyConstraint(parameters[ 3*k + p*d + 1 ]) : this.zetasUSD[k];
			newZetasEUR[k] = this.zetasEURInfo[k].getIsParameterToCalibrate() == true ? this.zetasEURInfo[k].getConstraint().applyConstraint(parameters[ 3*k + p*d + 2 ]) : this.zetasEUR[k];
			
		}
		
		double[] newLambdasJPY = new double[d];
		double[] newLambdasUSD = new double[d];
		double[] newLambdasEUR = new double[d];
		
		for(int k = 0; k < d; k++) {
			
			newLambdasJPY[k] = this.lambdasJPYInfo[k].getIsParameterToCalibrate() == true ? this.lambdasJPYInfo[k].getConstraint().applyConstraint(parameters[ 3*k + d*(p+3) ]) : this.lambdasJPY[k];
			newLambdasUSD[k] = this.lambdasUSDInfo[k].getIsParameterToCalibrate() == true ? this.lambdasUSDInfo[k].getConstraint().applyConstraint(parameters[ 3*k + d*(p+3) + 1 ]) : this.lambdasUSD[k];
			newLambdasEUR[k] = this.lambdasEURInfo[k].getIsParameterToCalibrate() == true ? this.lambdasEURInfo[k].getConstraint().applyConstraint(parameters[ 3*k + d*(p+3) + 2 ]) : this.lambdasEUR[k];
			
		}
		
		return new TemperedStableCBITCLMultiCurrencyModel(newTemperedStableCBITCLs, this.interestRateJPY, this.interestRateUSD, this.interestRateEUR, 
				this.initialValueJPYtoUSD, this.initialValueJPYtoEUR, this.initialValueUSDtoEUR,
				newZetasJPY, this.zetasJPYInfo,  newZetasUSD, this.zetasUSDInfo, newZetasEUR, this.zetasEURInfo, 
				newLambdasJPY, this.lambdasJPYInfo, newLambdasUSD, this.lambdasUSDInfo, newLambdasEUR, this.lambdasEURInfo);
	
	}

	@Override
	public double[] getParameterLowerBounds() {
		return parameterLowerBounds;
	}
	
	@Override
	public double[] getParameterUpperBounds() {
		return parameterUpperBounds;
	}
	
	/**
	 * This method allows to impose the usual constraints on the parameters of the model.
	 * It also informs whether one parameter or another has to be calibrated or not.
	 **/
	private void generateInformationForParameters() {
		
		for(int k = 0; k < this.getDimension(); k++) {
			
			this.zetasJPYInfo[k] = new ScalarParameterInformation(true, new BoundConstraint(0, 1));
			this.zetasUSDInfo[k] = new ScalarParameterInformation(true, new BoundConstraint(0, 1));
			this.zetasEURInfo[k] = new ScalarParameterInformation(true, new BoundConstraint(0, 1));
	
			this.lambdasJPYInfo[k] = new ScalarParameterInformation(true, new BoundConstraint(0, 1));
			this.lambdasUSDInfo[k] = new ScalarParameterInformation(true, new BoundConstraint(0, 1));
			this.lambdasEURInfo[k] = new ScalarParameterInformation(true, new BoundConstraint(0, 1));
			
		}
	
	}
	
	private double[] extractUpperBounds() {
		
		List<Double> upperBounds = new ArrayList<Double>();
		  
		for(int k = 0; k < this.getDimension(); k++) {
			
			if(this.temperedStableCBITCLs[k] instanceof TemperedStableCBITCBrownian) {
				
				int numParam = this.temperedStableCBITCLs[k].getNumberOfParameters();
				double[] paramUpperBounds = this.temperedStableCBITCLs[k].getParameterUpperBounds();
				for(int l = 0; l < numParam; l++) 
					upperBounds.add( paramUpperBounds[l] );
			
			} else if(this.temperedStableCBITCLs[k] instanceof TemperedStableCBITCLofCGMYtype) {
				
				int numParam = ((TemperedStableCBITCLofCGMYtype)this.temperedStableCBITCLs[k]).getNumberOfParameters();
				double[] paramUpperBounds = ((TemperedStableCBITCLofCGMYtype)this.temperedStableCBITCLs[k]).getParameterUpperBounds();
				for(int l = 0; l < numParam; l++) 
					upperBounds.add( paramUpperBounds[l] );
				
			}
			
		}
				
		double threshold = 1E6;
		
		for(int k = 0; k < this.getDimension(); k++) {
			upperBounds.add( this.zetasJPYInfo[k].getConstraint().getUpperBound() > threshold ? threshold : this.zetasJPYInfo[k].getConstraint().getUpperBound() );
			upperBounds.add( this.zetasUSDInfo[k].getConstraint().getUpperBound() > threshold ? threshold : this.zetasUSDInfo[k].getConstraint().getUpperBound() );
			upperBounds.add( this.zetasEURInfo[k].getConstraint().getUpperBound() > threshold ? threshold : this.zetasEURInfo[k].getConstraint().getUpperBound() );
		}
		
		for(int k = 0; k < this.getDimension(); k ++) {
			upperBounds.add(this.lambdasJPYInfo[k].getConstraint().getUpperBound() > threshold ? threshold : this.lambdasJPYInfo[k].getConstraint().getUpperBound());
			upperBounds.add(this.lambdasUSDInfo[k].getConstraint().getUpperBound() > threshold ? threshold : this.lambdasUSDInfo[k].getConstraint().getUpperBound());
			upperBounds.add(this.lambdasEURInfo[k].getConstraint().getUpperBound() > threshold ? threshold : this.lambdasEURInfo[k].getConstraint().getUpperBound());
		}

		Double[] upperBoundsArray = upperBounds.toArray(new Double[upperBounds.size()]);
		  
		return ArrayUtils.toPrimitive(upperBoundsArray);
		
	}

	private double[] extractLowerBounds() {
		
		List<Double> lowerBounds = new ArrayList<Double>();
		 
		for(int k = 0; k < this.getDimension(); k++) {
			
			if(this.temperedStableCBITCLs[k] instanceof TemperedStableCBITCBrownian) {
				
				int numParam = this.temperedStableCBITCLs[k].getNumberOfParameters();
				double[] paramLowerBounds = this.temperedStableCBITCLs[k].getParameterLowerBounds();
				for(int l = 0; l < numParam; l++) 
					lowerBounds.add( paramLowerBounds[l] );
			
			} else if(this.temperedStableCBITCLs[k] instanceof TemperedStableCBITCLofCGMYtype) {
				
				int numParam = ((TemperedStableCBITCLofCGMYtype)this.temperedStableCBITCLs[k]).getNumberOfParameters();
				double[] paramLowerBounds = ((TemperedStableCBITCLofCGMYtype)this.temperedStableCBITCLs[k]).getParameterLowerBounds();
				for(int l = 0; l < numParam; l++) 
					lowerBounds.add( paramLowerBounds[l] );
				
			}
			
		}

		double threshold = -1E6;
		
		for(int k = 0; k < this.getDimension(); k ++) {
			
			lowerBounds.add(this.zetasJPYInfo[k].getConstraint().getLowerBound() < threshold ? threshold : this.zetasJPYInfo[k].getConstraint().getLowerBound());
			lowerBounds.add(this.zetasUSDInfo[k].getConstraint().getLowerBound() < threshold ? threshold : this.zetasUSDInfo[k].getConstraint().getLowerBound());
			lowerBounds.add(this.zetasEURInfo[k].getConstraint().getLowerBound() < threshold ? threshold : this.zetasEURInfo[k].getConstraint().getLowerBound());
			
		}
		
		for(int k = 0; k < this.getDimension(); k ++) {
			
			lowerBounds.add(this.lambdasJPYInfo[k].getConstraint().getLowerBound() < threshold ? threshold : this.lambdasJPYInfo[k].getConstraint().getLowerBound());
			lowerBounds.add(this.lambdasUSDInfo[k].getConstraint().getLowerBound() < threshold ? threshold : this.lambdasUSDInfo[k].getConstraint().getLowerBound());
			lowerBounds.add(this.lambdasEURInfo[k].getConstraint().getLowerBound() < threshold ? threshold : this.lambdasEURInfo[k].getConstraint().getLowerBound());
			
		}
		
		Double[] lowerBoundsArray = lowerBounds.toArray(new Double[lowerBounds.size()]);
		  
		return ArrayUtils.toPrimitive(lowerBoundsArray);
		
	}

	@Override
	public double[] getParameters() {
		  
		  List<Double> params = new ArrayList<Double>();
		  
		  for(int k = 0; k < this.getDimension(); k++) {
			  
			  	if(this.temperedStableCBITCLs[k] instanceof TemperedStableCBITCBrownian) {
			  		
			  		int numParam = this.temperedStableCBITCLs[k].getNumberOfParameters();
			  		double[] parameters = this.temperedStableCBITCLs[k].getParameters();
					for(int l = 0; l < numParam; l++) 
						params.add( parameters[l] );
				
		  		} else if(this.temperedStableCBITCLs[k] instanceof TemperedStableCBITCLofCGMYtype) {
		  			
		  			int numParam = ((TemperedStableCBITCLofCGMYtype)this.temperedStableCBITCLs[k]).getNumberOfParameters();
			  		double[] parameters = ((TemperedStableCBITCLofCGMYtype)this.temperedStableCBITCLs[k]).getParameters();
					for(int l = 0; l < numParam; l++) 
						params.add( parameters[l] );
					
		  		}
				
		  }
		  
		  for(int k = 0; k < this.getDimension(); k++) {
				  params.add(this.zetasJPY[k]);
				  params.add(this.zetasUSD[k]);
				  params.add(this.zetasEUR[k]);
			  }
		  
		  for(int k = 0; k < this.getDimension(); k++) {
				  params.add(this.lambdasJPY[k]);
				  params.add(this.lambdasUSD[k]);
				  params.add(this.lambdasEUR[k]);
			  }

		  
		  Double[] array = params.toArray(new Double[params.size()]);
		  
		  return ArrayUtils.toPrimitive(array);
		  
	}
	
	@Override
	public Pair< List<Double>, MultivariateCalibrableProcessInterface > generateSamplePair() {
		
		//Creation of the Mersenne Twister pseudo-random number generator
		RandomGenerator rng = new MersenneTwister();
		RandomDataGenerator rndNumGen = new RandomDataGenerator(rng);
		
		//Retrieving of data for each tempered stable CBITCL process
		int d = this.getDimension();
		
		@SuppressWarnings("unchecked")
		Pair< List<Double>, TemperedStableCBITCLProcess >[] temperedStableCBITCLsData = ( Pair< List<Double>, TemperedStableCBITCLProcess >[] )new Pair[d];
		
		for(int k = 0; k < d; k++)
			temperedStableCBITCLsData[k] = this.temperedStableCBITCLs[k].generateSamplePair();
			
		//Creation of the parameter list
		List<Double> params = new ArrayList<Double>(this.getParameters().length);
		TemperedStableCBITCLProcess[] newTemperedStableCBITCLs = new TemperedStableCBITCLProcess[d];
			
		double[] newZetasJPY = new double[d];
		double[] newZetasUSD = new double[d];
		double[] newZetasEUR = new double[d];
			
		double[] newLambdasJPY = new double[d];
		double[] newLambdasUSD = new double[d];
		double[] newLambdasEUR = new double[d];
			
		for(int k = 0; k < d; k++) {
				
			//Addition of the randomly-generated parameter set of temperedStableCBITCLs[k]
			params.addAll(  temperedStableCBITCLsData[k].getFirst() );
				
			//Retrieving of the sample specification of temperedStableCBITCLs[k]
			newTemperedStableCBITCLs[k] = temperedStableCBITCLsData[k].getSecond();
				
			//Random generation of the zetas
			newZetasJPY[k] = rndNumGen.nextUniform(0.0, 1.0);
			newZetasUSD[k] = rndNumGen.nextUniform(0.0, 1.0);
			newZetasEUR[k] = rndNumGen.nextUniform(0.0, 1.0);
			
			//Random generation of the lambdas
			newLambdasJPY[k] = rndNumGen.nextUniform(0.0, 1.0);
			newLambdasUSD[k] = rndNumGen.nextUniform(0.0, 1.0);
			newLambdasEUR[k] = rndNumGen.nextUniform(0.0, 1.0);
				
		}
			
		//Creation of the tempered stable CBITCL multi-currency model
		TemperedStableCBITCLMultiCurrencyModel model = new  TemperedStableCBITCLMultiCurrencyModel(newTemperedStableCBITCLs, 
				this.interestRateJPY, this.interestRateUSD, this.interestRateEUR, this.initialValueJPYtoUSD, this.initialValueJPYtoEUR, this.initialValueUSDtoEUR,
				newZetasJPY, newZetasUSD, newZetasEUR, newLambdasJPY, newLambdasUSD, newLambdasEUR); 
			
		//Addition of the zetas to the parameter set
		for(int k = 0; k < d; k++) {
			params.add(newZetasJPY[k]);
			params.add(newZetasUSD[k]);
			params.add(newZetasEUR[k]);
		}
		  
		//Addition of the lambdas to the parameter set
		for(int k = 0; k < d; k++) {
			params.add(newLambdasJPY[k]);
			params.add(newLambdasUSD[k]);
			params.add(newLambdasEUR[k]);
		}
			
		return Pair.< List<Double>, MultivariateCalibrableProcessInterface >create(params, model);
		
	}
	
	/* First part: $0 < \zeta_3 \leq \zeta_2 \leq \zeta_1 < \theta$, 
	 * Second part (provided that the tempered stable CBITCL processes are of CGMY type): $-G_L < 0 < \lambda_3 \leq \lambda_2 \leq \lambda_1 < M_L$,
	 * for 1 \leq k \leq d$ where $d$ is the size of temperedStableCBITCLs and where we recall that JPY-1, USD-2, EUR-3.
	 */
	private void martingaleMeasureAssumption(int k) {
		
		//First part:
		if(this.zetasJPY[k] >= this.temperedStableCBITCLs[k].getTheta()) 
			this.zetasJPY[k] = this.temperedStableCBITCLs[k].getTheta()*0.5;	
		
		if(this.zetasUSD[k] > this.zetasJPY[k])
			this.zetasUSD[k] = this.zetasJPY[k]*0.5;
		
		if(this.zetasEUR[k] > this.zetasUSD[k]) 
			this.zetasEUR[k] = this.zetasUSD[k]*0.5; 
		
		//Second part:
		if(this.temperedStableCBITCLs[k] instanceof TemperedStableCBITCLofCGMYtype) {
			
			if(this.lambdasJPY[k] >= ((TemperedStableCBITCLofCGMYtype)this.temperedStableCBITCLs[k]).getM() ) 
				this.lambdasJPY[k] = ((TemperedStableCBITCLofCGMYtype)this.temperedStableCBITCLs[k]).getM()*0.5;
			
			if(this.lambdasUSD[k] > this.lambdasJPY[k]) 
				this.lambdasUSD[k] = this.lambdasJPY[k]*0.5;
			
			if(this.lambdasEUR[k] > this.lambdasUSD[k]) 
				this.lambdasEUR[k] = this.lambdasUSD[k]*0.5;
		
		}
		
	}
	
	/* This generates the $\QQ_{JPY}$,$\QQ_{USD}$ and $\QQ_{EUR}$-versions of the tempered stable CBITCL processes, 
	 * also including tempered stable CBITCL processes of CGMY type. This is permitted thanks to 
	 * the validity of the martingale measure assumption of the multi-currency tempered stable CBITCL modeling framework.
	 */
	private void generateRiskNeutralTemperedStableCBITCLs() {
		
		for(int k = 0; k < this.getDimension(); k++) {
			
			TemperedStableCBITCLProcess temperedStableCBITCLproc = this.temperedStableCBITCLs[k];
			
			double v0JPY = temperedStableCBITCLproc.getInitialValueOfCBI(); 
			double v0USD = v0JPY;
			double v0EUR = v0JPY;
			
			double betaJPY = temperedStableCBITCLproc.getBeta();
			double betaUSD = betaJPY;
			double betaEUR = betaJPY;
			
			double sigmaJPY = temperedStableCBITCLproc.getSigma(); 
			double sigmaUSD = sigmaJPY;
			double sigmaEUR = sigmaJPY;
			
			double etaJPY = temperedStableCBITCLproc.getEta();
			double etaUSD = etaJPY;
			double etaEUR = etaJPY;
			
			double thetaJPY = temperedStableCBITCLproc.getTheta() - this.zetasJPY[k];
			double thetaUSD = temperedStableCBITCLproc.getTheta() - this.zetasUSD[k];
			double thetaEUR = temperedStableCBITCLproc.getTheta() - this.zetasEUR[k];
			
			double alphaJPY = temperedStableCBITCLproc.getAlpha();
			double alphaUSD = alphaJPY;
			double alphaEUR = alphaJPY;
		
			double bJPY = temperedStableCBITCLproc.getb() - 
					sigmaJPY*sigmaJPY*this.zetasJPY[k] - alphaJPY*Math.pow(etaJPY, alphaJPY)*(
							Math.pow(temperedStableCBITCLproc.getTheta(), alphaJPY - 1) - Math.pow(thetaJPY, alphaJPY - 1) );
			
			double bUSD = temperedStableCBITCLproc.getb() - 
					sigmaUSD*sigmaUSD*this.zetasUSD[k] - alphaUSD*Math.pow(etaUSD, alphaUSD)*(
							Math.pow(temperedStableCBITCLproc.getTheta(), alphaUSD - 1) - Math.pow(thetaUSD, alphaUSD - 1) );
			
			double bEUR = temperedStableCBITCLproc.getb() - 
					sigmaEUR*sigmaEUR*this.zetasEUR[k] - alphaEUR*Math.pow(etaEUR, alphaEUR)*(
							Math.pow(temperedStableCBITCLproc.getTheta(), alphaEUR - 1) - Math.pow(thetaEUR, alphaEUR - 1) );
			
			if(temperedStableCBITCLproc instanceof TemperedStableCBITCBrownian) {
			
				this.temperedStableCBITCLsJPY[k] = new TemperedStableCBITCBrownian(this.timeHorizon, this.numberOfTimeSteps, 
						v0JPY, betaJPY, bJPY, sigmaJPY, etaJPY, thetaJPY, alphaJPY);
				
				this.temperedStableCBITCLsUSD[k] = new TemperedStableCBITCBrownian(this.timeHorizon, this.numberOfTimeSteps, 
						v0USD, betaUSD, bUSD, sigmaUSD, etaUSD, thetaUSD, alphaUSD);
						
				this.temperedStableCBITCLsEUR[k] = new TemperedStableCBITCBrownian(this.timeHorizon, this.numberOfTimeSteps, 
						v0EUR, betaEUR, bEUR, sigmaEUR, etaEUR, thetaEUR, alphaEUR);
						
			} else if(temperedStableCBITCLproc instanceof TemperedStableCBITCLofCGMYtype) {
				
				double P = ((TemperedStableCBITCLofCGMYtype)temperedStableCBITCLproc).getP();
				
				double gJPY = ((TemperedStableCBITCLofCGMYtype)temperedStableCBITCLproc).getG() + this.lambdasJPY[k];
				double gUSD = ((TemperedStableCBITCLofCGMYtype)temperedStableCBITCLproc).getG() + this.lambdasUSD[k];
				double gEUR = ((TemperedStableCBITCLofCGMYtype)temperedStableCBITCLproc).getG() + this.lambdasEUR[k];
				
				double mJPY = ((TemperedStableCBITCLofCGMYtype)temperedStableCBITCLproc).getM() - this.lambdasJPY[k];
				double mUSD = ((TemperedStableCBITCLofCGMYtype)temperedStableCBITCLproc).getM() - this.lambdasUSD[k];
				double mEUR = ((TemperedStableCBITCLofCGMYtype)temperedStableCBITCLproc).getM() - this.lambdasEUR[k];
				
				double yJPY = ((TemperedStableCBITCLofCGMYtype)temperedStableCBITCLproc).getY();
				double yUSD = yJPY;
				double yEUR = yJPY;
				
				double betaLJPY = ((TemperedStableCBITCLofCGMYtype)temperedStableCBITCLproc).getBetaL() + yJPY*(
						Math.pow(((TemperedStableCBITCLofCGMYtype)temperedStableCBITCLproc).getM(), yJPY - 1) - Math.pow(mJPY, yJPY - 1) 
						+ Math.pow(((TemperedStableCBITCLofCGMYtype)temperedStableCBITCLproc).getG(), yJPY - 1) - Math.pow(gJPY, yJPY - 1) );
				
				double betaLUSD = ((TemperedStableCBITCLofCGMYtype)temperedStableCBITCLproc).getBetaL() + yUSD*(
						Math.pow(((TemperedStableCBITCLofCGMYtype)temperedStableCBITCLproc).getM(), yUSD - 1) - Math.pow(mUSD, yUSD - 1) 
						+ Math.pow(((TemperedStableCBITCLofCGMYtype)temperedStableCBITCLproc).getG(), yUSD - 1) - Math.pow(gUSD, yUSD - 1) );
				
				double betaLEUR = ((TemperedStableCBITCLofCGMYtype)temperedStableCBITCLproc).getBetaL() + yEUR*(
						Math.pow(((TemperedStableCBITCLofCGMYtype)temperedStableCBITCLproc).getM(), yEUR - 1) - Math.pow(mEUR, yEUR - 1) 
						+ Math.pow(((TemperedStableCBITCLofCGMYtype)temperedStableCBITCLproc).getG(), yEUR - 1) - Math.pow(gEUR, yEUR - 1) );
				
				this.temperedStableCBITCLsJPY[k] = new TemperedStableCBITCLofCGMYtype(P, this.timeHorizon, this.numberOfTimeSteps, 
						v0JPY, betaJPY, bJPY, sigmaJPY, etaJPY, thetaJPY, alphaJPY, betaLJPY, gJPY, mJPY, yJPY);
				
				this.temperedStableCBITCLsUSD[k] = new TemperedStableCBITCLofCGMYtype(P, this.timeHorizon, this.numberOfTimeSteps, 
						v0USD, betaUSD, bUSD, sigmaUSD, etaUSD, thetaUSD, alphaUSD, betaLUSD, gUSD, mUSD, yUSD);
						
				this.temperedStableCBITCLsEUR[k] = new TemperedStableCBITCLofCGMYtype(P, this.timeHorizon, this.numberOfTimeSteps, 
						v0EUR, betaEUR, bEUR, sigmaEUR, etaEUR, thetaEUR, alphaEUR, betaLEUR, gEUR, mEUR, yEUR);
					
			}
			
		}
		
	}
	
	/* The following returns the constant minP needed for the determination of the extended domain of the characteristic function.
	 * The latter corresponds to the minimum value of all the P's of the tempered stable CBITCL processes (of CGMY type).
	 */
	private double computationOfMinP() {
		
		double minP = 1.5;//> 1, value if for every $k$, temperedStableCBITCLs[k] only defines a tempered stable CBI-time-changed Brownian motion.
		
		if(this.temperedStableCBITCLs[0] instanceof TemperedStableCBITCLofCGMYtype) {
			
			minP = ((TemperedStableCBITCLofCGMYtype)this.temperedStableCBITCLs[0]).getP();
			
			for(int k = 1; k < this.dimension; k++) {
				
				if( ((TemperedStableCBITCLofCGMYtype)this.temperedStableCBITCLs[k]).getP() <= minP ) 
					minP = ((TemperedStableCBITCLofCGMYtype)this.temperedStableCBITCLs[k]).getP();
					
			}
			
			return minP;
			
		} else {
			
			return minP;
			
		}
		
	}

}
