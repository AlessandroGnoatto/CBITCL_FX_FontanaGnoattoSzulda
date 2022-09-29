package org.calibrationframework.currencyoptioncalibration;

import org.calibrationframework.exception.CalculationException;
import org.calibrationframework.fouriermethod.calibration.models.*;
import org.calibrationframework.fouriermethod.products.*;
import org.calibrationframework.marketdata.model.volatilities.OptionSmileData;
import org.calibrationframework.marketdata.model.volatilities.OptionSurfaceData;
import org.calibrationframework.marketdata.model.volatilities.WeightedOptionSmileData;
import org.calibrationframework.marketdata.model.volatilities.WeightedOptionSurfaceData;
import org.calibrationframework.marketdata.model.volatilities.VolatilitySurfaceInterface.QuotingConvention;

import net.finmath.marketdata.model.curves.*;
import net.finmath.marketdata.model.curves.CurveInterpolation.ExtrapolationMethod;
import net.finmath.marketdata.model.curves.CurveInterpolation.InterpolationEntity;
import net.finmath.marketdata.model.curves.CurveInterpolation.InterpolationMethod;
import org.calibrationframework.stochastic.*;
import org.nd4j.common.primitives.Pair;

import java.time.LocalDate;
import java.util.*;

public class TestTemperedStableCBITCLMultiCurrencyModel {
	
	/**
	 * Get the discount curve using the riskFreeRate.
	 * 
	 * @param name Name of the curve
	 * @param referenceDate Date corresponding to t=0.
	 * @param riskFreeRate Constant continuously compounded rate
	 * 
	 * @return the discount curve using the riskFreeRate.
	 */
	private static DiscountCurve getDiscountCurve(String name, LocalDate referenceDate, double riskFreeRate) {
		double[] times = new double[] { 1.0 };
		double[] givenAnnualizedZeroRates = new double[] { riskFreeRate };
		InterpolationMethod interpolationMethod = InterpolationMethod.LINEAR;
		InterpolationEntity interpolationEntity = InterpolationEntity.LOG_OF_VALUE_PER_TIME;
		ExtrapolationMethod extrapolationMethod = ExtrapolationMethod.CONSTANT;
		DiscountCurve discountCurve = DiscountCurveInterpolation.createDiscountCurveFromAnnualizedZeroRates(name, referenceDate, times, givenAnnualizedZeroRates, interpolationMethod, extrapolationMethod, interpolationEntity);
		return discountCurve;
	}


	public static void main(String[] args) throws CalculationException {
		
		LocalDate referenceDate = LocalDate.of(2020, 4, 15); //t_0
		
		QuotingConvention convention = QuotingConvention.VOLATILITYLOGNORMAL;
		
		double[] maturities = {/*1.0/365,*/ 7.0/365, 14.0/365, 1.0/12, 3.0/12, 6.0/12, 1.0/*,2.0*/};
		
		/*
		 * EURUSD
		 */
		/*double[] volEURUSD0 = {0.119675, 0.11001250000000001, 0.10587500000000001, 0.1077875, 0.116675};*/
		double[] volEURUSD1 = {0.09563750000000001, 0.08641249999999999, 0.08145, 0.0830875/*, 0.08946249999999999*/};
		double[] volEURUSD2 = {0.0957875, 0.08443749999999998, 0.078175, 0.08011249999999999/*, 0.0877625*/};
		double[] volEURUSD3 = {0.0984, 0.085725, 0.07885, 0.08057500000000001/*, 0.0891*/};
		double[] volEURUSD4 = {0.09882500000000001, 0.08462500000000002, 0.07690000000000001, 0.07867500000000001/*, 0.088075*/};
		double[] volEURUSD5 = {0.10085000000000001, 0.08511250000000001, 0.07655, 0.07793750000000002/*, 0.08755*/};
		double[] volEURUSD6 = {0.10387500000000001, 0.0870625, 0.078, 0.0792875/*, 0.089325*/};
		/*double[] volEURUSD7 = {0.10718749999999999, 0.090375, 0.08145, 0.083025, 0.09311249999999999};*/

		/*double[] strikesEURUSD0 = {1.0830916533264612, 1.0875857847425179, 1.091816765264554, 1.0959800649113822, 1.1003990162806785};*/
		double[] strikesEURUSD1 = {1.0736665485947974, 1.0832493552102738, 1.092019466260985, 1.1005303661154329/*, 1.1095108545186112*/};
		double[] strikesEURUSD2 = {1.0901695159533655, 1.1042668007956644, 1.1166308654893506, 1.128517131608227/*, 1.1415350335611791*/};
		double[] strikesEURUSD3 = {1.0924070710388567, 1.114045792484283, 1.1327434049550982, 1.1506679179112815/*, 1.170787899495015*/};
		double[] strikesEURUSD4 = {1.0276082002234361, 1.0636364621501644, 1.0942585784210839, 1.1237196462929744/*, 1.1580572459547283*/};
		double[] strikesEURUSD5 = {1.0013144785277734, 1.052721141626843, 1.0959543709704656, 1.1375198977428984/*, 1.1869821738969657*/};
		double[] strikesEURUSD6 = {0.9653754927945374, 1.038269302786766, 1.100241850177062, 1.1608005301886808/*, 1.2348531928649114*/};
		/*double[] strikesEURUSD7 = {0.9183740977750597, 1.0197801730846652, 1.1098887457916413, 1.2016722135780769, 1.3165985533837115};*/
		
		double[] weightEURUSD1 = {0.0129453183, 0.0147946494, 0.0184933118, 0.0147946494/*, 0.0129453183*/};
		double[] weightEURUSD2 = {0.0159326994, 0.0182087993, 0.0227609991, 0.0182087993/*, 0.0159326994*/};
		double[] weightEURUSD3 = {0.0207125092, 0.0236714391, 0.0295892989, 0.0236714391/*, 0.0207125092*/};
		double[] weightEURUSD4 = {0.0258906365, 0.0295892989, 0.0369866236, 0.0295892989/*, 0.0258906365*/};
		double[] weightEURUSD5 = {0.0295892989, 0.0338163416, 0.042270427, 0.0338163416/*, 0.0295892989*/};
		double[] weightEURUSD6 = {0.0318653988, 0.0364175986, 0.0455219983, 0.0364175986/*, 0.0318653988*/};
		/*double[] weightEURUSD7 = {0.0295892989, 0.0338163416, 0.042270427, 0.0338163416, 0.0295892989};*/
		
		double[][] volEURUSD = {/*volEURUSD0,*/ volEURUSD1, volEURUSD2, volEURUSD3, volEURUSD4,volEURUSD5,volEURUSD6/*,volEURUSD7*/};		
		double[][] strikesEURUSD = {/*strikesEURUSD0,*/strikesEURUSD1,strikesEURUSD2,strikesEURUSD3,strikesEURUSD4,strikesEURUSD5,strikesEURUSD6/*,strikesEURUSD7*/};
		double[][] weightEURUSD = {weightEURUSD1,weightEURUSD2,weightEURUSD3,weightEURUSD4,weightEURUSD5,weightEURUSD6/*,weightEURUSD7*/};
		
		double[] forwardsEURUSD = {/*1.0918,*/ 1.09195, 1.1165, 1.13245,1.09345, 1.09435, 1.0969/*, 1.10255*/};
				
		DiscountCurve discountCurveUSD = getDiscountCurve("discountCurve", referenceDate, 0.0);
		
		DiscountCurve forwardCurveEURUSD = DiscountCurveInterpolation.createDiscountCurveFromDiscountFactors(
				"EURUSDforwardCurve"								/* name */,
				maturities,
				forwardsEURUSD);
		
		/*
		 *EURJPY 
		 */
		/*double[] volEURJPY0 = {0.171075, 0.132825, 0.114375, 0.108325, 0.12357500000000002};*/
		double[] volEURJPY1 = {0.150925, 0.11357499999999998, 0.09392499999999998, 0.08607499999999998/*, 0.099525*/};
		double[] volEURJPY2 = {0.1442125, 0.1153875, 0.094375, 0.08476249999999999/*, 0.0898375*/};
		double[] volEURJPY3 = {0.145675, 0.1183125, 0.0951, 0.0838375/*, 0.08252500000000002*/};
		double[] volEURJPY4 = {0.154675, 0.12383749999999999, 0.09705, 0.0839125/*, 0.081075*/};
		double[] volEURJPY5 = {0.1595625, 0.1267625, 0.09795, 0.08378749999999999/*, 0.0806375*/};
		double[] volEURJPY6 = {0.1649, 0.12983750000000002, 0.09897500000000001, 0.08416250000000002/*, 0.08055*/};
		/*double[] volEURJPY7 = {0.15788775, 0.12465000000000001, 0.096225, 0.08285000000000001, 0.08061225000000001};*/

		/*double[] strikesEURJPY0 = {115.97119608544445, 116.75903024478252, 117.3071021309632, 117.75636713701685, 118.28389744718632};*/
		double[] strikesEURJPY1 = {114.22963028417423, 116.08153529139031, 117.31497367803952, 118.26038325033285/*, 119.40678776105385*/};
		double[] strikesEURJPY2 = {113.18017202328735, 115.56010762444403, 117.32508885513762, 118.64222676330878/*, 119.99868210078732*/};
		double[] strikesEURJPY3 = {111.24886090978765, 114.70057494889274, 117.34931285065109, 119.27061222617105/*, 120.97601384481221*/};
		double[] strikesEURJPY4 = {106.5544296105295, 112.72269589869013, 117.44318894366741, 120.77829648932939/*, 123.66173835127871*/};
		double[] strikesEURJPY5 = {102.16056564352186, 110.86738944725796, 117.58644908112925, 122.30178450902785/*, 126.4028531242797*/};
		double[] strikesEURJPY6 = {96.25913015092581, 108.37858355444538, 117.88077128916112, 124.59694434569336/*, 130.48376402131106*/};
		/*double[] strikesEURJPY7 = {90.33725909078565, 105.78541528153434, 118.39594812217175, 127.82571221721972, 136.64334728497815};*/

		double[] weightEURJPY1 = {0.0088356809, 0.010097921, 0.0126224013, 0.010097921/*, 0.0088356809*/};
		double[] weightEURJPY2 = {0.0138057514, 0.0157780016, 0.019722502, 0.0157780016/*, 0.0138057514*/};
		double[] weightEURJPY3 = {0.0232517918, 0.0265734764, 0.0332168454, 0.0265734764/*, 0.0232517918*/};
		double[] weightEURJPY4 = {0.0266134967, 0.0304154247, 0.0380192809, 0.0304154247/*, 0.0266134967*/};
		double[] weightEURJPY5 = {0.0334684882, 0.0382497008, 0.047812126, 0.0382497008/*, 0.0334684882*/};
		double[] weightEURJPY6 = {0.0356277455, 0.0407174234, 0.0508967793, 0.0407174234/*, 0.0356277455*/};
		/*double[] weightEURJPY7 = {0.0276115028, 0.0315560032, 0.039445004, 0.0315560032, 0.0276115028};*/
		
		double[][] volEURJPY = {/*volEURJPY0,*/ volEURJPY1, volEURJPY2, volEURJPY3, volEURJPY4, volEURJPY5, volEURJPY6/*, volEURJPY7*/};
		double[][] strikesEURJPY = {/*strikesEURJPY0,*/strikesEURJPY1,strikesEURJPY2,strikesEURJPY3,strikesEURJPY4,strikesEURJPY5,strikesEURJPY6/*,strikesEURJPY7*/};
		double[][] weightEURJPY = {weightEURJPY1,weightEURJPY2,weightEURJPY3,weightEURJPY4,weightEURJPY5,weightEURJPY6/*,weightEURJPY7*/};
		
		double[] forwardsEURJPY = {/*117.305,*/ 117.30505, 117.30505, 117.3051, 117.305, 117.30475, 117.3048/*, 117.30475*/};
		
		DiscountCurve discountCurveJPY = getDiscountCurve("discountCurve", referenceDate, 0.0);
		
		DiscountCurve forwardCurveEURJPY = DiscountCurveInterpolation.createDiscountCurveFromDiscountFactors(
				"EURUSDforwardCurve"								/* name */,
				maturities,
				forwardsEURJPY);
		
		/*
		 * USDJPY
		 */
		/*double[] volUSDJPY0 = {0.11852500000000002, 0.10423750000000001, 0.09265000000000001, 0.08656250000000001, 0.08637500000000001};*/
		double[] volUSDJPY1 = {0.11247499999999999, 0.0960375, 0.08199999999999999, 0.0752125/*, 0.074625*/};
		double[] volUSDJPY2 = {0.116425, 0.09887499999999999, 0.08355, 0.07502500000000001/*, 0.072875*/};
		double[] volUSDJPY3 = {0.125275, 0.10513750000000001, 0.08715, 0.07691250000000001/*, 0.073825*/};
		double[] volUSDJPY4 = {0.13245, 0.10792500000000001, 0.0868, 0.074975/*, 0.07145*/};
		double[] volUSDJPY5 = {0.137275, 0.10953750000000001, 0.08615, 0.0731125/*, 0.06992499999999999*/};
		double[] volUSDJPY6 = {0.14104999999999998, 0.11109999999999999, 0.08654999999999999, 0.07309999999999998/*, 0.07094999999999999*/};
		/*double[] volUSDJPY7 = {0.13604999999999998, 0.1061375, 0.08324999999999999, 0.0711125, 0.07135};*/

		/*double[] strikesUSDJPY0 = {106.56642359493159, 107.02202754497327, 107.41626309342197, 107.7448723378179, 108.04027203309737};*/
		double[] strikesUSDJPY1 = {105.30476233695445, 106.46505440574764, 107.42182598288703, 108.17805214717387/*, 108.85279052133608*/};
		double[] strikesUSDJPY2 = {104.348487456532, 106.04079825876437, 107.42913106591406, 108.49629613950225/*, 109.40867035392269*/};
		double[] strikesUSDJPY3 = {102.616733972536, 105.28619167006954, 107.44824806207954, 109.06181664671308/*, 110.41341677644552*/};
		double[] strikesUSDJPY4 = {98.88902086988593, 103.7241446335822, 107.51365654683141, 110.24042781110688/*, 112.51625101015048*/};
		double[] strikesUSDJPY5 = {95.29380256622512, 102.2487565742581, 107.60927912269112, 111.36997740160989/*, 114.57611977337466*/};
		double[] strikesUSDJPY6 = {90.539424617409, 100.26698495302524, 107.80738205426069, 113.13453545550729/*, 117.92439736975528*/};
		/*double[] strikesUSDJPY7 = {85.49382664367916, 98.15358270072305, 108.14149018241402, 115.51485668769966, 122.84391669499219};*/

		double[] weightUSDJPY1 = {0.0088356809, 0.010097921, 0.0126224013, 0.010097921/*, 0.0088356809*/};
		double[] weightUSDJPY2 = {0.0138057514, 0.0157780016, 0.019722502, 0.0157780016/*, 0.0138057514*/};
		double[] weightUSDJPY3 = {0.0232517918, 0.0265734764, 0.0332168454, 0.0265734764/*, 0.0232517918*/};
		double[] weightUSDJPY4 = {0.0266134967, 0.0304154247, 0.0380192809, 0.0304154247/*, 0.0266134967*/};
		double[] weightUSDJPY5 = {0.0334684882, 0.0382497008, 0.047812126, 0.0382497008/*, 0.0334684882*/};
		double[] weightUSDJPY6 = {0.0356277455, 0.0407174234, 0.0508967793, 0.0407174234/*, 0.0356277455*/};
		/*double[] weightUSDJPY7 = {0.0276115028, 0.0315560032, 0.039445004, 0.0315560032, 0.0276115028};*/
		
		double[][] volUSDJPY = {/*volUSDJPY0,*/ volUSDJPY1, volUSDJPY2, volUSDJPY3, volUSDJPY4, volUSDJPY5, volUSDJPY6/*, volUSDJPY7*/};
		double[][] strikesUSDJPY = {/*strikesUSDJPY0,*/strikesUSDJPY1,strikesUSDJPY2, strikesUSDJPY3,strikesUSDJPY4,strikesUSDJPY5,strikesUSDJPY6/*,strikesUSDJPY7*/};
		double[][] weightUSDJPY = {weightUSDJPY1,weightUSDJPY2,weightUSDJPY3,weightUSDJPY4,weightUSDJPY5,weightUSDJPY6/*,weightUSDJPY7*/};
		
		double[] forwardsUSDJPY = {/*107.415,*/ 107.4149, 107.41475, 107.41425,107.41245, 107.4098, 107.40435/*, 107.3946*/};
		
		DiscountCurve forwardCurveUSDJPY = DiscountCurveInterpolation.createDiscountCurveFromDiscountFactors(
				"USDJPYforwardCurve"								/* name */,
				maturities,
				forwardsUSDJPY);				
		
		OptionSmileData[] smilesEURUSD = new OptionSmileData[maturities.length];
		OptionSmileData[] smilesEURJPY = new OptionSmileData[maturities.length];
		OptionSmileData[] smilesUSDJPY = new OptionSmileData[maturities.length];
				
		double maturity;
		
		for(int i = 0; i< maturities.length; i++) {
			
			maturity = maturities[i];
			
			smilesEURUSD[i] = new OptionSmileData("EURUSD", referenceDate, strikesEURUSD[i], maturity, volEURUSD[i], convention);
			smilesEURJPY[i] = new OptionSmileData("EURJPY", referenceDate, strikesEURJPY[i], maturity, volEURJPY[i], convention);
			smilesUSDJPY[i] = new OptionSmileData("USDJPY", referenceDate, strikesUSDJPY[i], maturity, volUSDJPY[i], convention);
			
		}
		
		OptionSurfaceData surfaceEURUSD = new OptionSurfaceData(smilesEURUSD, discountCurveUSD, forwardCurveEURUSD);
		OptionSurfaceData surfaceEURJPY = new OptionSurfaceData(smilesEURJPY, discountCurveJPY, forwardCurveEURJPY);
		OptionSurfaceData surfaceUSDJPY = new OptionSurfaceData(smilesUSDJPY, discountCurveJPY, forwardCurveUSDJPY);
		
		LinkedHashMap<String, OptionSurfaceData> surfaces = new LinkedHashMap<String, OptionSurfaceData>();
		surfaces.put("EURUSD", surfaceEURUSD);
		surfaces.put("EURJPY", surfaceEURJPY);
		surfaces.put("USDJPY", surfaceUSDJPY);
		
		/*
		 * Initial Model Set-up
		 */
		double timeHorizon = 1.2;
		int numberOfTimeSteps = 200;
		
		double P = 1.5;
		
		int d = 2; //Dimension of the underlying process.
		
		double interestRateJPY = 0.00; 
		double interestRateUSD = 0.00;
		double interestRateEUR = 0.00;
		
		double initialValueJPYtoUSD = 107;
		double initialValueJPYtoEUR = 117;
		double initialValueUSDtoEUR = 1.09; 
		
		double[] zetasJPY = {0.6, 0.7};
		double[] zetasUSD = {0.9, 0.3};
		double[] zetasEUR = {0.8, 0.1};
		
		double[] lambdasJPY = {0.6, 0.7};
		double[] lambdasUSD = {0.9, 0.6};
		double[] lambdasEUR = {0.8, 0.9};
		
		double[] v0 = {2.1, 3.3, 1.2};
		
		double[] beta = {2.1, 3.7, 1.3};
		
		double[] b = {-0.3, 0.5, 1.1};
		double[] sigma = {3.1, 3.3, 0.3};
		double[] eta = {0.1, 0.3, 1.4};
		double[] theta = {3.3, 3.6, 2.3};
		double[] alpha = {1.5, 1.3, 1.1};
		
		double[] betaL = {-0.1, 1.2};
		double[] G = {1.2, 1.3};
		double[] M = {1.1, 1.9};
		double[] Y = {1.3, 1.8};
			
		TemperedStableCBITCLofCGMYtype[] temperedStableCBITCLs = new TemperedStableCBITCLofCGMYtype[d];
		
		for(int k = 0; k < d; k++) {
			temperedStableCBITCLs[k] =  new TemperedStableCBITCLofCGMYtype(P, timeHorizon, numberOfTimeSteps, 
					v0[k], beta[k], b[k], sigma[k], eta[k], theta[k], alpha[k],
					betaL[k], G[k], M[k], Y[k]);
		}
		
		TemperedStableCBITCLMultiCurrencyModel model = new TemperedStableCBITCLMultiCurrencyModel(temperedStableCBITCLs, interestRateJPY, interestRateUSD, interestRateEUR, 
				initialValueJPYtoUSD, initialValueJPYtoEUR, initialValueUSDtoEUR,
				zetasJPY, zetasUSD, zetasEUR, lambdasJPY, lambdasUSD, lambdasEUR );
		
		//Best-fit parameter set:
		double[] bestFitParams =  new double[] {
				1.2561819794656286, 
				-0.2194282411176097, 
				1.8386482622023261,
				0.6250931769322416,
				2.7658286437245834,
				4.138611331841888,
				1.113904250715546,
				0.4845411702381398,
				1.0442417787840679,
				2.8998252841886254,
				4.134687760246713,
				-0.7437580013804046,
				-0.9637809864784093,
				-0.8120698400968506,
				-0.2607496775743369,
				-1.4251887582654499,
				1.181714659783091,
				2.063967576240257,
				0.7046734515636803,
				1.495726665130293,
				-0.47803882217384946,
				1.165432621761013,
				0.7855777400208042,
				0.1668331754294542,
				1.6802597420093992,
				1.7915051355902205,
				0.1425805625718937,
				1.8708607080100648,
				2.398632857679031,
				0.27339988214145244,
				-1.5030473637101807,
				-0.9052326371672694,
				2.464823373221091,
				-1.8114347405955225
		};
		
		TemperedStableCBITCLMultiCurrencyModel calibratedModel = model.getCloneForModifiedParameters(bestFitParams);
		
		double[] calibratedParams = calibratedModel.getParameters();
		
		for(int i = 0; i < calibratedParams.length; i++)
			System.out.println(calibratedParams[i]);

		/*
		//Generation of the model to be tested
		Pair< List<Double>, MultivariateCalibrableProcessInterface > modelPair = model.generateSamplePair();
		
		//Retrieving of the parameters
		List<Double> params = modelPair.getLeft();
		
		//Retrieving of the model
		MultivariateCalibrableProcessInterface newModel = modelPair.getRight();
		
		//Creation of the pricer
		double myTime = maturities[0];
		
		double[] myStrikes = strikesEURUSD1;
		double myStrike = myStrikes[0];
		
		double lineOfIntegration = 0.5*(1 + model.getMinP());
		
		EuropeanOptionSmileMultiAsset pricer = /*new CurrencyOptionByCOSmethod("EURUSD", myTime, new double[] {myStrike},
				discountCurveUSD, forwardCurveEURUSD, discountCurveJPY, forwardCurveEURJPY, forwardCurveUSDJPY, maturities);*/
				/*new CurrencyOptionByCarrMadan("EURUSD", myTime, new double[] {myStrike}, lineOfIntegration);

		Iterator<Map.Entry<String, OptionSurfaceData>> iterator = surfaces.entrySet().iterator();
				
		//Loop over the surfaces
		while (iterator.hasNext()) {
			
			Map.Entry<String, OptionSurfaceData> pair = (Map.Entry<String, OptionSurfaceData>) iterator.next();
					
			String underlying = pair.getKey();
			
			System.out.println(underlying);
					
			OptionSurfaceData surface = pair.getValue(); 
					
			int numberOfMaturities = surface.getMaturities().length;
					
			double mats[] = surface.getMaturities();
					
			QuotingConvention targetConvention = surface.getQuotingConvention();
					
			for(int t = 0; t<numberOfMaturities; t++) {
				
				System.out.println("T = " + mats[t]);
						
				double[] currentStrikes = surface.getSmile(mats[t]).getStrikes();
						
				EuropeanOptionSmileMultiAsset newPricer = pricer.getCloneWithModifiedParameters(underlying, mats[t],currentStrikes);
						
				try {
					
					Map<Double, Double> currentModelPrices = newPricer.getValue(newModel);
							
					for(int i = 0; i < currentStrikes.length; i++) {
						
						System.out.println("K = " + currentStrikes[i]);
					
						if(targetConvention.equals(QuotingConvention.VOLATILITYLOGNORMAL)) {
									
							//we convert prices into lognormal volatilities
							double forward = surface.getEquityForwardCurve().getDiscountFactor(mats[t]);
							double optionMaturity = mats[t];
							double optionStrike = currentStrikes[i];
							double payoffUnit = surface.getDiscountCurve().getDiscountFactor(mats[t]);
							double optionValue = currentModelPrices.get(currentStrikes[i]);
							System.out.println(optionValue);
							System.out.println( net.finmath.functions.AnalyticFormulas.blackScholesOptionImpliedVolatility(forward, optionMaturity, optionStrike, payoffUnit, optionValue) );
								
						} else if(targetConvention.equals(QuotingConvention.VOLATILITYNORMAL)) {
									
							//we convert prices into normal volatilities
							double forward = surface.getEquityForwardCurve().getDiscountFactor(mats[t]);
							double optionMaturity =mats[t];
							double optionStrike = currentStrikes[i];
							double payoffUnit = surface.getDiscountCurve().getDiscountFactor(mats[t]);
							double optionValue = currentModelPrices.get(currentStrikes[i]);
							System.out.println(optionValue);
							System.out.println( net.finmath.functions.AnalyticFormulas.bachelierOptionImpliedVolatility(forward, optionMaturity, optionStrike, payoffUnit, optionValue) );
									
						} else {
							
							//just output the prices
							System.out.println( currentModelPrices.get(currentStrikes[i]) );
									
						}						
								
					}	
							
				} catch (CalculationException e) {
					e.printStackTrace();
				}						
							
			}
					
		}*/
		
	}
	
}
