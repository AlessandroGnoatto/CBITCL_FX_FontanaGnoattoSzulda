package org.calibrationframework.capletcalibration;

import java.time.LocalDate;
import java.util.*;
import java.util.function.*;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.*;

import org.calibrationframework.exception.CalculationException;
import org.calibrationframework.fouriermethod.calibration.models.CBIDrivenMultiCurveModel;
import org.calibrationframework.fouriermethod.products.*;
import org.calibrationframework.fouriermethod.quantization.QuantizableCBIDrivenMultiCurveModel;
import org.calibrationframework.fouriermethod.quantization.QuantizationMultiCurveCapletPricer;
import org.calibrationframework.quantization.*;

import net.finmath.marketdata.calibration.CalibratedCurves;
import net.finmath.marketdata.calibration.CalibratedCurves.CalibrationSpec;
import net.finmath.marketdata.model.*;
import net.finmath.marketdata.model.curves.*;
import net.finmath.optimizer.SolverException;

import org.calibrationframework.stochastic.*;

import net.finmath.marketdata.model.curves.CurveInterpolation.*;
import net.finmath.time.*;
import net.finmath.time.businessdaycalendar.*;


public class TestCBIProcess {

	public static void main(String[] args) throws SolverException, CloneNotSupportedException, CalculationException {	
		
		/*
		 * Calibration of a single curve - OIS curve - self disocunted curve, from a set of calibration products.
		 */
		LocalDate referenceDate = LocalDate.of(2018,9,24);

		/*
		 * Define the calibration spec generators for our calibration products
		 */
		Function<String,String> frequencyForTenor = (tenor) -> {
			switch(tenor) {
			case "3M":
				return "quarterly";
			case "6M":
				return "semiannual";
			}
			throw new IllegalArgumentException("Unkown tenor " + tenor);
		};

		BiFunction<String, Double, CalibrationSpec> deposit = (maturity, rate) -> {
			Schedule scheduleInterfaceRec = ScheduleGenerator.createScheduleFromConventions(referenceDate, 2, "0D", maturity, "tenor", "act/360", "first", "following", new BusinessdayCalendarExcludingTARGETHolidays(), 0, 0);
			Schedule scheduleInterfacePay = null;
			double calibrationTime = scheduleInterfaceRec.getPayment(scheduleInterfaceRec.getNumberOfPeriods()-1);
			CalibrationSpec calibrationSpec = new CalibratedCurves.CalibrationSpec("EUR-OIS-" + maturity, "Deposit", scheduleInterfaceRec, "", rate, "discount-EUR-OIS", scheduleInterfacePay, null, 0.0, null, "discount-EUR-OIS", calibrationTime);
			return calibrationSpec;
		};

		BiFunction<String, Double, CalibrationSpec> swapSingleCurve = (maturity, rate) -> {
			Schedule scheduleInterfaceRec = ScheduleGenerator.createScheduleFromConventions(referenceDate, 2, "0D", maturity, "annual", "act/360", "first", "modified_following", new BusinessdayCalendarExcludingTARGETHolidays(), 0, 1);
			Schedule scheduleInterfacePay = ScheduleGenerator.createScheduleFromConventions(referenceDate, 2, "0D", maturity, "annual", "act/360", "first", "modified_following", new BusinessdayCalendarExcludingTARGETHolidays(), 0, 1);
			double calibrationTime = scheduleInterfaceRec.getPayment(scheduleInterfaceRec.getNumberOfPeriods() - 1);
			CalibrationSpec calibrationSpec = new CalibratedCurves.CalibrationSpec("EUR-OIS-" + maturity, "Swap", scheduleInterfaceRec, "forward-EUR-OIS", 0.0, "discount-EUR-OIS", scheduleInterfacePay, "", rate, "discount-EUR-OIS", "discount-EUR-OIS", calibrationTime);
			return calibrationSpec;
		};

		Function<String,BiFunction<String, Double, CalibrationSpec>> fra = (tenor) -> {
			return (fixing, rate) -> {
				Schedule scheduleInterfaceRec = ScheduleGenerator.createScheduleFromConventions(referenceDate, 2, fixing, tenor, "tenor", "act/360", "first", "modified_following", new BusinessdayCalendarExcludingTARGETHolidays(), 0, 0);
				double calibrationTime = scheduleInterfaceRec.getFixing(scheduleInterfaceRec.getNumberOfPeriods() - 1);
				String curveName = "forward-EUR-" + tenor;
				CalibrationSpec calibrationSpec = new CalibratedCurves.CalibrationSpec("EUR-" + tenor + "-" + fixing, "FRA", scheduleInterfaceRec, curveName, rate, "discount-EUR-OIS", null, null, 0.0, null, curveName, calibrationTime);
				return calibrationSpec;
			};
		};

		Function<String,BiFunction<String, Double, CalibrationSpec>> swap = (tenor) -> {
			return (maturity, rate) -> {
				String frequencyRec = frequencyForTenor.apply(tenor);

				Schedule scheduleInterfaceRec = ScheduleGenerator.createScheduleFromConventions(referenceDate, 2, "0D", maturity, frequencyRec, "act/360", "first", "following", new BusinessdayCalendarExcludingTARGETHolidays(), 0, 0);
				Schedule scheduleInterfacePay = ScheduleGenerator.createScheduleFromConventions(referenceDate, 2, "0D", maturity, "annual", "E30/360", "first", "following", new BusinessdayCalendarExcludingTARGETHolidays(), 0, 0);
				double calibrationTime = scheduleInterfaceRec.getFixing(scheduleInterfaceRec.getNumberOfPeriods() - 1);
				String curveName = "forward-EUR-" + tenor;
				CalibrationSpec calibrationSpec = new CalibratedCurves.CalibrationSpec("EUR-" + tenor + maturity, "Swap", scheduleInterfaceRec, curveName, 0.0, "discount-EUR-OIS", scheduleInterfacePay, "", rate, "discount-EUR-OIS", curveName, calibrationTime);
				return calibrationSpec;
			};
		};

		BiFunction<String,String,BiFunction<String, Double, CalibrationSpec>> swapBasis = (tenorRec,tenorPay) -> {
			return (maturity, rate) -> {
				String curveNameRec = "forward-EUR-" + tenorRec;
				String curveNamePay = "forward-EUR-" + tenorPay;

				String frequencyRec = frequencyForTenor.apply(tenorRec);
				String frequencyPay = frequencyForTenor.apply(tenorPay);

				Schedule scheduleInterfaceRec = ScheduleGenerator.createScheduleFromConventions(referenceDate, 2, "0D", maturity, frequencyRec, "act/360", "first", "following", new BusinessdayCalendarExcludingTARGETHolidays(), 0, 0);
				Schedule scheduleInterfacePay = ScheduleGenerator.createScheduleFromConventions(referenceDate, 2, "0D", maturity, frequencyPay, "act/360", "first", "following", new BusinessdayCalendarExcludingTARGETHolidays(), 0, 0);
				double calibrationTime = scheduleInterfaceRec.getFixing(scheduleInterfaceRec.getNumberOfPeriods() - 1);

				CalibrationSpec calibrationSpec = new CalibratedCurves.CalibrationSpec("EUR-" + tenorRec + "-" + tenorPay + maturity, "Swap", scheduleInterfaceRec, curveNameRec, 0.0, "discount-EUR-OIS", scheduleInterfacePay, curveNamePay, rate, "discount-EUR-OIS", curveNameRec, calibrationTime);
				return calibrationSpec;
			};
		};

		/*
		 * Generate empty curve template (for cloning during calibration)
		 */
		double[] times = { 0.0 };
		double[] discountFactors = { 1.0 };
		boolean[] isParameter = { false };

		DiscountCurve discountCurveOIS = DiscountCurveInterpolation.createDiscountCurveFromDiscountFactors("discount-EUR-OIS", referenceDate, times, discountFactors, isParameter, InterpolationMethod.LINEAR, ExtrapolationMethod.CONSTANT, InterpolationEntity.LOG_OF_VALUE);
		ForwardCurve forwardCurveOIS = new ForwardCurveFromDiscountCurve("forward-EUR-OIS", "discount-EUR-OIS", referenceDate, "3M");
		ForwardCurve forwardCurve3M = new ForwardCurveInterpolation("forward-EUR-3M", referenceDate, "3M", new BusinessdayCalendarExcludingTARGETHolidays(), BusinessdayCalendar.DateRollConvention.FOLLOWING, CurveInterpolation.InterpolationMethod.LINEAR, CurveInterpolation.ExtrapolationMethod.CONSTANT, CurveInterpolation.InterpolationEntity.VALUE,ForwardCurveInterpolation.InterpolationEntityForward.FORWARD, "discount-EUR-OIS");
		ForwardCurve forwardCurve6M = new ForwardCurveInterpolation("forward-EUR-6M", referenceDate, "6M", new BusinessdayCalendarExcludingTARGETHolidays(), BusinessdayCalendar.DateRollConvention.FOLLOWING, CurveInterpolation.InterpolationMethod.LINEAR, CurveInterpolation.ExtrapolationMethod.CONSTANT, CurveInterpolation.InterpolationEntity.VALUE,ForwardCurveInterpolation.InterpolationEntityForward.FORWARD, "discount-EUR-OIS");

		AnalyticModelFromCurvesAndVols forwardCurveModel = new AnalyticModelFromCurvesAndVols(new Curve[] { discountCurveOIS, forwardCurveOIS, forwardCurve3M, forwardCurve6M });

		List<CalibrationSpec> calibrationSpecs = new LinkedList<>();

		/*
		 * Calibration products for OIS curve: Deposits
		 */
		calibrationSpecs.add(deposit.apply("1D", 0.202 / 100.0));
		calibrationSpecs.add(deposit.apply("1W", 0.195 / 100.0));
		calibrationSpecs.add(deposit.apply("2W", 0.193 / 100.0));
		calibrationSpecs.add(deposit.apply("3W", 0.193 / 100.0));
		calibrationSpecs.add(deposit.apply("1M", 0.191 / 100.0));
		calibrationSpecs.add(deposit.apply("2M", 0.185 / 100.0));
		calibrationSpecs.add(deposit.apply("3M", 0.180 / 100.0));
		calibrationSpecs.add(deposit.apply("4M", 0.170 / 100.0));
		calibrationSpecs.add(deposit.apply("5M", 0.162 / 100.0));
		calibrationSpecs.add(deposit.apply("6M", 0.156 / 100.0));
		calibrationSpecs.add(deposit.apply("7M", 0.150 / 100.0));
		calibrationSpecs.add(deposit.apply("8M", 0.145 / 100.0));
		calibrationSpecs.add(deposit.apply("9M", 0.141 / 100.0));
		calibrationSpecs.add(deposit.apply("10M", 0.136 / 100.0));
		calibrationSpecs.add(deposit.apply("11M", 0.133 / 100.0));
		calibrationSpecs.add(deposit.apply("12M", 0.129 / 100.0));

		/*
		 * Calibration products for OIS curve: Swaps
		 */
		calibrationSpecs.add(swapSingleCurve.apply("15M", 0.118 / 100.0));
		calibrationSpecs.add(swapSingleCurve.apply("18M", 0.108 / 100.0));
		calibrationSpecs.add(swapSingleCurve.apply("21M", 0.101 / 100.0));
		calibrationSpecs.add(swapSingleCurve.apply("2Y", 0.101 / 100.0));
		calibrationSpecs.add(swapSingleCurve.apply("3Y", 0.194 / 100.0));
		calibrationSpecs.add(swapSingleCurve.apply("4Y", 0.346 / 100.0));
		calibrationSpecs.add(swapSingleCurve.apply("5Y", 0.534 / 100.0));
		calibrationSpecs.add(swapSingleCurve.apply("6Y", 0.723 / 100.0));
		calibrationSpecs.add(swapSingleCurve.apply("7Y", 0.895 / 100.0));
		calibrationSpecs.add(swapSingleCurve.apply("8Y", 1.054 / 100.0));
		calibrationSpecs.add(swapSingleCurve.apply("9Y", 1.189 / 100.0));
		calibrationSpecs.add(swapSingleCurve.apply("10Y", 1.310 / 100.0));
		calibrationSpecs.add(swapSingleCurve.apply("11Y", 1.423 / 100.0));
		calibrationSpecs.add(swapSingleCurve.apply("12Y", 1.520 / 100.0));
		calibrationSpecs.add(swapSingleCurve.apply("15Y", 1.723 / 100.0));
		calibrationSpecs.add(swapSingleCurve.apply("20Y", 1.826 / 100.0));
		calibrationSpecs.add(swapSingleCurve.apply("25Y", 1.877 / 100.0));
		calibrationSpecs.add(swapSingleCurve.apply("30Y", 1.910 / 100.0));
		calibrationSpecs.add(swapSingleCurve.apply("40Y", 2.025 / 100.0));
		calibrationSpecs.add(swapSingleCurve.apply("50Y", 2.101 / 100.0));

		/*
		 * Calibration products for 3M curve: FRAs
		 */
		calibrationSpecs.add(fra.apply("3M").apply("0D", 0.322 / 100.0));
		calibrationSpecs.add(fra.apply("3M").apply("1M", 0.329 / 100.0));
		calibrationSpecs.add(fra.apply("3M").apply("2M", 0.328 / 100.0));
		calibrationSpecs.add(fra.apply("3M").apply("3M", 0.326 / 100.0));
		calibrationSpecs.add(fra.apply("3M").apply("6M", 0.323 / 100.0));
		calibrationSpecs.add(fra.apply("3M").apply("9M", 0.316 / 100.0));
		calibrationSpecs.add(fra.apply("3M").apply("12M", 0.360 / 100.0));
		calibrationSpecs.add(fra.apply("3M").apply("15M", 0.390 / 100.0));

		/*
		 * Calibration products for 3M curve: swaps
		 */
		calibrationSpecs.add(swap.apply("3M").apply("2Y", 0.380 / 100.0));
		calibrationSpecs.add(swap.apply("3M").apply("3Y", 0.485 / 100.0));
		calibrationSpecs.add(swap.apply("3M").apply("4Y", 0.628 / 100.0));
		calibrationSpecs.add(swap.apply("3M").apply("5Y", 0.812 / 100.0));
		calibrationSpecs.add(swap.apply("3M").apply("6Y", 0.998 / 100.0));
		calibrationSpecs.add(swap.apply("3M").apply("7Y", 1.168 / 100.0));
		calibrationSpecs.add(swap.apply("3M").apply("8Y", 1.316 / 100.0));
		calibrationSpecs.add(swap.apply("3M").apply("9Y", 1.442 / 100.0));
		calibrationSpecs.add(swap.apply("3M").apply("10Y", 1.557 / 100.0));
		calibrationSpecs.add(swap.apply("3M").apply("12Y", 1.752 / 100.0));
		calibrationSpecs.add(swap.apply("3M").apply("15Y", 1.942 / 100.0));
		calibrationSpecs.add(swap.apply("3M").apply("20Y", 2.029 / 100.0));
		calibrationSpecs.add(swap.apply("3M").apply("25Y", 2.045 / 100.0));
		calibrationSpecs.add(swap.apply("3M").apply("30Y", 2.097 / 100.0));
		calibrationSpecs.add(swap.apply("3M").apply("40Y", 2.208 / 100.0));
		calibrationSpecs.add(swap.apply("3M").apply("50Y", 2.286 / 100.0));

		/*
		 * Calibration products for 6M curve: FRAs
		 */

		calibrationSpecs.add(fra.apply("6M").apply("0D", 0.590 / 100.0));
		calibrationSpecs.add(fra.apply("6M").apply("1M", 0.597 / 100.0));
		calibrationSpecs.add(fra.apply("6M").apply("2M", 0.596 / 100.0));
		calibrationSpecs.add(fra.apply("6M").apply("3M", 0.594 / 100.0));
		calibrationSpecs.add(fra.apply("6M").apply("6M", 0.591 / 100.0));
		calibrationSpecs.add(fra.apply("6M").apply("9M", 0.584 / 100.0));
		calibrationSpecs.add(fra.apply("6M").apply("12M", 0.584 / 100.0));

		/*
		 * Calibration products for 6M curve: tenor basis swaps
		 * Note: the fixed bases is added to the second argument tenor (here 3M).
		 */
		calibrationSpecs.add(swapBasis.apply("6M","3M").apply("2Y", 0.255 / 100.0));
		calibrationSpecs.add(swapBasis.apply("6M","3M").apply("3Y", 0.245 / 100.0));
		calibrationSpecs.add(swapBasis.apply("6M","3M").apply("4Y", 0.227 / 100.0));
		calibrationSpecs.add(swapBasis.apply("6M","3M").apply("5Y", 0.210 / 100.0));
		calibrationSpecs.add(swapBasis.apply("6M","3M").apply("6Y", 0.199 / 100.0));
		calibrationSpecs.add(swapBasis.apply("6M","3M").apply("7Y", 0.189 / 100.0));
		calibrationSpecs.add(swapBasis.apply("6M","3M").apply("8Y", 0.177 / 100.0));
		calibrationSpecs.add(swapBasis.apply("6M","3M").apply("9Y", 0.170 / 100.0));
		calibrationSpecs.add(swapBasis.apply("6M","3M").apply("10Y", 0.164 / 100.0));
		calibrationSpecs.add(swapBasis.apply("6M","3M").apply("12Y", 0.156 / 100.0));
		calibrationSpecs.add(swapBasis.apply("6M","3M").apply("15Y", 0.135 / 100.0));
		calibrationSpecs.add(swapBasis.apply("6M","3M").apply("20Y", 0.125 / 100.0));
		calibrationSpecs.add(swapBasis.apply("6M","3M").apply("25Y", 0.117 / 100.0));
		calibrationSpecs.add(swapBasis.apply("6M","3M").apply("30Y", 0.107 / 100.0));
		calibrationSpecs.add(swapBasis.apply("6M","3M").apply("40Y", 0.095 / 100.0));
		calibrationSpecs.add(swapBasis.apply("6M","3M").apply("50Y", 0.088 / 100.0));

		/*
		 * Calibrate
		 */
		CalibratedCurves calibratedCurves = new CalibratedCurves(calibrationSpecs.toArray(new CalibrationSpec[calibrationSpecs.size()]), forwardCurveModel, 1E-15);

		/*
		 * Get the calibrated model
		 */
		AnalyticModel calibratedModel = calibratedCurves.getModel();
		
		//Params of the FFT pricer:
		int numberOfPoints = 1024;
		double gridSpacing = 0.05;
		net.finmath.interpolation.RationalFunctionInterpolation.InterpolationMethod intMethod = net.finmath.interpolation.RationalFunctionInterpolation.InterpolationMethod.LINEAR;
		net.finmath.interpolation.RationalFunctionInterpolation.ExtrapolationMethod extMethod = net.finmath.interpolation.RationalFunctionInterpolation.ExtrapolationMethod.CONSTANT;
		
		/*
		//Flow of CBI params:
		double b = 0.05769718466962348;
		double sigma = 0.015789007994048317;
		double eta = 0.05918343595932772;
		double zeta = 0.06918343595932772;
		double alpha = 1.6979723707482892;
		double[] initialValues = {0.0012, 0.0013};
		double[] immigrationRates = {0.001, 0.005};
		double[] lambda = {1.3937638262955285, 1.0};
		*/
		
		//Flow of CBI params:
		double b = 0.22926988700762607;
		double sigma = 0.03351058076450054;
		double eta =  0.09734037637289966;
		double zeta = 0.507210173470597;
		double alpha = 1.3220545188973394;
		double[] initialValues = {0.0012256627870067618, 1.0E-4};
		double[] immigrationRates = {0.9916429538266227, 0.1161235882057263};
		double[] lambda = {1.3937638262955285, 1.0};
		
		/*
		double b = 0.05353549346164644;
		double sigma = 0.005827989181896429;
		double eta =  0.04070169217539017;
		double zeta = 0.050701692175390174;
		double alpha = 1.3175352727830814;
		double[] initialValues = {0.004953850642168643, 0.005076590407389615};
		double[] immigrationRates = {0.0001, 0.0034047019048037384};
		double[] lambda = {1.4999999999998428, 1.0000000348864304};
		*/
		
		MultiCurveTenor threeMonth = new MultiCurveTenor(0.25, "3M");
		MultiCurveTenor sixMonth = new MultiCurveTenor(0.5, "6M");
		MultiCurveTenor[] tenors = {threeMonth,sixMonth};
		
		double timeHorizon = 10.0;
		int numberOfTimeSteps = 500;
		
		CBIDrivenMultiCurveModel model = new CBIDrivenMultiCurveModel(timeHorizon, numberOfTimeSteps, calibratedModel, tenors, initialValues, immigrationRates, b, sigma, eta, zeta, alpha, lambda);
		
		int N = 10;
		
		double[] strikes = {-0.05};
		
		MultiCurveTenor tenor = new MultiCurveTenor(0.5, "forward-EUR-6M");
	
		System.out.println("Strike: " + strikes[0] + ", Tenor: " + tenor.getTenorName() + ", Level of Quantization: " + N);
		
		System.out.println("FFT price             " + "Quantization price    " + "Maturity   ");
		
		double[] initialGuess = new double[N];
		
		double[] previous = new double[N];
		
		double[] next = new double[N];
	
		initialGuess = null;
	
		for(int i = 20; i < 40; i++) {
			
			double maturity = 0.1 + i*0.1;
			
			previous = initialGuess;
			
			CapletByCarrMadan smile = new CapletByCarrMadan(tenor.getTenorName(), maturity, strikes, numberOfPoints, gridSpacing,  intMethod,  extMethod);
			
			Map<Double, Double> prices = smile.getValue(model);
			
			QuantizableCBIDrivenMultiCurveModel quantizedModel = new QuantizableCBIDrivenMultiCurveModel(maturity, N, model, tenor, initialGuess);
			
			QuantizationMultiCurveCapletPricer myCapletQunt = new QuantizationMultiCurveCapletPricer(strikes[0], maturity, tenor);
			
			next = quantizedModel.getQuantizationGrid();
			
			if(previous != null) {
				RealVector previousVect = new ArrayRealVector(previous);
				RealVector nextVect = new ArrayRealVector(next);
				double distance = previousVect.subtract(nextVect).getNorm();
				if(distance >= 0.1) {
					initialGuess = null;
					quantizedModel = new QuantizableCBIDrivenMultiCurveModel(maturity, N, model, tenor, initialGuess);
					initialGuess = quantizedModel.getQuantizationGrid();
				} else {
					initialGuess = quantizedModel.getQuantizationGrid();
				}
			} else {
				initialGuess = quantizedModel.getQuantizationGrid();
			}
					
			System.out.println(prices.get(strikes[0]) + "  " + myCapletQunt.getValue(quantizedModel) + "  " + maturity);
			
		}	

	}

}