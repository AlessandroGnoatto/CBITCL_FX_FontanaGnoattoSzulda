package org.calibrationframework.capletcalibration;

import java.time.LocalDate;
import java.util.*;
import java.util.function.*;

import org.calibrationframework.fouriermethod.calibration.*;
import org.calibrationframework.fouriermethod.calibration.constraints.BoundConstraint;
import org.calibrationframework.fouriermethod.calibration.constraints.PositivityConstraint;
import org.calibrationframework.fouriermethod.calibration.constraints.ScalarParameterInformation;
import org.calibrationframework.fouriermethod.calibration.constraints.ScalarParameterInformationInterface;
import org.calibrationframework.fouriermethod.calibration.models.*;
import org.calibrationframework.fouriermethod.quantization.QuantizationMultiCurveSmileOfCaplets;

import net.finmath.marketdata.calibration.CalibratedCurves;
import net.finmath.marketdata.calibration.CalibratedCurves.CalibrationSpec;
import net.finmath.marketdata.model.*;
import net.finmath.marketdata.model.curves.*;
import net.finmath.marketdata.model.curves.CurveInterpolation.*;
import org.calibrationframework.marketdata.model.volatilities.*;
import org.calibrationframework.marketdata.model.volatilities.VolatilitySurfaceInterface.QuotingConvention;
import org.calibrationframework.optimizer.*;
import org.calibrationframework.quantization.*;
import org.calibrationframework.stochastic.*;
import net.finmath.time.*;
import net.finmath.time.businessdaycalendar.*;
import net.finmath.optimizer.SolverException;

public class CapletCalibrationQuantizationTest {
	
	public static void main(String[] args) throws SolverException, CloneNotSupportedException, org.calibrationframework.optimizer.SolverException, SolverException {
		
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
		//calibrationSpecs.add(deposit.apply("1D", 0.202 / 100.0));
		calibrationSpecs.add(deposit.apply("1W", -0.353 / 100.0));
		calibrationSpecs.add(deposit.apply("2W", -0.357 / 100.0));
		calibrationSpecs.add(deposit.apply("3W", -0.358 / 100.0));
		calibrationSpecs.add(deposit.apply("1M", -0.359 / 100.0));
		calibrationSpecs.add(deposit.apply("2M", -0.359/ 100.0));
		calibrationSpecs.add(deposit.apply("3M", -0.359/ 100.0));
		calibrationSpecs.add(deposit.apply("4M", -0.359/ 100.0));
		calibrationSpecs.add(deposit.apply("5M", -0.359/ 100.0));
		calibrationSpecs.add(deposit.apply("6M", -0.358/ 100.0));
		calibrationSpecs.add(deposit.apply("7M", -0.358/ 100.0));
		calibrationSpecs.add(deposit.apply("8M", -0.357/ 100.0));
		calibrationSpecs.add(deposit.apply("9M", -0.357/ 100.0));
		calibrationSpecs.add(deposit.apply("10M",-0.356 / 100.0));
		calibrationSpecs.add(deposit.apply("11M",-0.356 / 100.0));
		calibrationSpecs.add(deposit.apply("12M",-0.355 / 100.0));

		/*
		 * Calibration products for OIS curve: Swaps
		 */
		calibrationSpecs.add(swapSingleCurve.apply("15M", -0.351 / 100.0));
		calibrationSpecs.add(swapSingleCurve.apply("18M", -0.340 / 100.0));
		//calibrationSpecs.add(swapSingleCurve.apply("21M", 0.101 / 100.0));
		calibrationSpecs.add(swapSingleCurve.apply("2Y", -0.323 / 100.0));
		calibrationSpecs.add(swapSingleCurve.apply("3Y", -0.185/ 100.0));
		calibrationSpecs.add(swapSingleCurve.apply("4Y", -0.049/ 100.0));
		calibrationSpecs.add(swapSingleCurve.apply("5Y", 0.088/ 100.0));
		calibrationSpecs.add(swapSingleCurve.apply("6Y", 0.223 / 100.0));
		calibrationSpecs.add(swapSingleCurve.apply("7Y", 0.352 / 100.0));
		calibrationSpecs.add(swapSingleCurve.apply("8Y", 0.475 / 100.0));
		calibrationSpecs.add(swapSingleCurve.apply("9Y", 0.590 / 100.0));
		calibrationSpecs.add(swapSingleCurve.apply("10Y", 0.694 / 100.0));
		calibrationSpecs.add(swapSingleCurve.apply("11Y", 0.788 / 100.0));
		calibrationSpecs.add(swapSingleCurve.apply("12Y", 0.873 / 100.0));
		calibrationSpecs.add(swapSingleCurve.apply("15Y", 1.065 / 100.0));
		calibrationSpecs.add(swapSingleCurve.apply("20Y", 1.250 / 100.0));
		calibrationSpecs.add(swapSingleCurve.apply("25Y", 1.320 / 100.0));
		calibrationSpecs.add(swapSingleCurve.apply("30Y", 1.346 / 100.0));
		calibrationSpecs.add(swapSingleCurve.apply("40Y", 1.364 / 100.0));
		calibrationSpecs.add(swapSingleCurve.apply("50Y", 1.353 / 100.0));

		/*
		 * Calibration products for 3M curve: FRAs
		 */
		calibrationSpecs.add(fra.apply("3M").apply("0D", -0.324 / 100.0));
		calibrationSpecs.add(fra.apply("3M").apply("1M", -0.313 / 100.0));
		calibrationSpecs.add(fra.apply("3M").apply("2M", -0.305 / 100.0));
		calibrationSpecs.add(fra.apply("3M").apply("3M", -0.297 / 100.0));
		//calibrationSpecs.add(fra.apply("3M").apply("6M", 0.323 / 100.0));
		//calibrationSpecs.add(fra.apply("3M").apply("9M", 0.316 / 100.0));
		//calibrationSpecs.add(fra.apply("3M").apply("12M", 0.360 / 100.0));
		//calibrationSpecs.add(fra.apply("3M").apply("15M", 0.390 / 100.0));

		/*
		 * Calibration products for 3M curve: swaps
		 */
		calibrationSpecs.add(swap.apply("3M").apply("2Y", -0.231 / 100.0));
		calibrationSpecs.add(swap.apply("3M").apply("3Y", -0.102 / 100.0));
		calibrationSpecs.add(swap.apply("3M").apply("4Y", 0.046 / 100.0));
		calibrationSpecs.add(swap.apply("3M").apply("5Y", 0.194 / 100.0));
		calibrationSpecs.add(swap.apply("3M").apply("6Y", 0.336 / 100.0));
		calibrationSpecs.add(swap.apply("3M").apply("7Y", 0.471 / 100.0));
		calibrationSpecs.add(swap.apply("3M").apply("8Y", 0.598 / 100.0));
		calibrationSpecs.add(swap.apply("3M").apply("9Y", 0.715 / 100.0));
		calibrationSpecs.add(swap.apply("3M").apply("10Y", 0.822 / 100.0));
		calibrationSpecs.add(swap.apply("3M").apply("12Y", 1.003 / 100.0));
		calibrationSpecs.add(swap.apply("3M").apply("15Y", 1.198 / 100.0));
		calibrationSpecs.add(swap.apply("3M").apply("20Y", 1.380 / 100.0));
		calibrationSpecs.add(swap.apply("3M").apply("25Y", 1.447 / 100.0));
		calibrationSpecs.add(swap.apply("3M").apply("30Y", 1.471 / 100.0));
		calibrationSpecs.add(swap.apply("3M").apply("40Y", 1.481 / 100.0));
		calibrationSpecs.add(swap.apply("3M").apply("50Y", 1.463 / 100.0));

		/*
		 * Calibration products for 6M curve: FRAs
		 */

		calibrationSpecs.add(fra.apply("6M").apply("0D", -0.260 / 100.0));
		calibrationSpecs.add(fra.apply("6M").apply("1M", -0.260 / 100.0));
		calibrationSpecs.add(fra.apply("6M").apply("2M", -0.252 / 100.0));
		calibrationSpecs.add(fra.apply("6M").apply("3M", -0.243 / 100.0));
		calibrationSpecs.add(fra.apply("6M").apply("6M", -0.223 / 100.0));
		calibrationSpecs.add(fra.apply("6M").apply("9M", -0.208 / 100.0));
		calibrationSpecs.add(fra.apply("6M").apply("12M",-0.166 / 100.0));

		/*
		 * Calibration products for 6M curve:  swaps
		 */
		calibrationSpecs.add(swap.apply("6M").apply("2Y", -0.173 / 100.0));
		calibrationSpecs.add(swap.apply("6M").apply("3Y", -0.041 / 100.0));
		calibrationSpecs.add(swap.apply("6M").apply("4Y", 0.110 / 100.0));
		calibrationSpecs.add(swap.apply("6M").apply("5Y", 0.259 / 100.0));
		calibrationSpecs.add(swap.apply("6M").apply("6Y", 0.404 / 100.0));
		calibrationSpecs.add(swap.apply("6M").apply("7Y", 0.539 / 100.0));
		calibrationSpecs.add(swap.apply("6M").apply("8Y", 0.667 / 100.0));
		calibrationSpecs.add(swap.apply("6M").apply("9Y", 0.784 / 100.0));
		calibrationSpecs.add(swap.apply("6M").apply("10Y", 0.890 / 100.0));
		calibrationSpecs.add(swap.apply("6M").apply("12Y", 1.072 / 100.0));
		calibrationSpecs.add(swap.apply("6M").apply("15Y", 1.259 / 100.0));
		calibrationSpecs.add(swap.apply("6M").apply("20Y", 1.431 / 100.0));
		calibrationSpecs.add(swap.apply("6M").apply("25Y", 1.492 / 100.0));
		calibrationSpecs.add(swap.apply("6M").apply("30Y", 1.511 / 100.0));
		calibrationSpecs.add(swap.apply("6M").apply("40Y", 1.516 / 100.0));
		calibrationSpecs.add(swap.apply("6M").apply("50Y", 1.511 / 100.0));

		/*
		 * Calibrate
		 */
		CalibratedCurves calibratedCurves = new CalibratedCurves(calibrationSpecs.toArray(new CalibrationSpec[calibrationSpecs.size()]), forwardCurveModel, 1E-15);
		
		
		/*System.out.println(calibratedCurves.getCurve("discount-EUR-OIS").toString());
		System.out.println(calibratedCurves.getCurve("forward-EUR-3M").toString());
		System.out.println(calibratedCurves.getCurve("forward-EUR-6M").toString());*/
		int threads = Runtime.getRuntime().availableProcessors();
		System.out.println("Available processors to the Java Virtual Machine " +threads);
		/*
		 * Get the calibrated model
		 */
		AnalyticModel calibratedModel = calibratedCurves.getModel();
		
		double[] strikes = {-0.005,
				-0.0025,
				-0.0013,
				0.0,
				0.0025,
				0.005,
				0.01,
				0.015,
				0.02,
				0.03,
				0.05};
		
	
		double[][] volatilities = 
			{{0.0025,0.0024,0.0027,0.0031,0.0037,0.0043,0.0052,0.0061,0.0069,0.0085,0.0114},
					{0.0025,0.0024,0.0027,0.0031,0.0037,0.0043,0.0052,0.0061,0.0069,0.0085,0.0114},
					{0.0025,0.0024,0.0028,0.0031,0.0038,0.0043,0.0052,0.0061,0.0069,0.0085,0.0114},
					{0.0025,0.0024,0.0027,0.0031,0.0037,0.0043,0.0052,0.0061,0.0069,0.0085,0.0114},
					{0.0025,0.0024,0.0027,0.0031,0.0038,0.0043,0.0052,0.0061,0.0069,0.0085,0.0114},
					{0.0040,0.0043,0.0045,0.0047,0.0049,0.0052,0.0058,0.0065,0.0072,0.0087,0.0114},
					{0.0047,0.0049,0.0051,0.0053,0.0055,0.0057,0.0062,0.0068,0.0075,0.0089,0.0115},
					{0.0050,0.0052,0.0053,0.0054,0.0056,0.0058,0.0063,0.0069,0.0075,0.0088,0.0113},
					{0.0055,0.0056,0.0057,0.0058,0.0060,0.0061,0.0065,0.0070,0.0076,0.0088,0.0113},
					{0.0055,0.0056,0.0057,0.0058,0.0060,0.0061,0.0064,0.0069,0.0074,0.0084,0.0107},
					{0.0058,0.0059,0.0061,0.0061,0.0062,0.0064,0.0066,0.0070,0.0074,0.0083,0.0103},
					{0.0058,0.0059,0.0059,0.0060,0.0062,0.0063,0.0066,0.0069,0.0073,0.0081,0.0099},
					{0.0060,0.0061,0.0062,0.0063,0.0064,0.0065,0.0067,0.0069,0.0072,0.0079,0.0095},
					{0.0059,0.0060,0.0061,0.0061,0.0063,0.0063,0.0066,0.0068,0.0071,0.0079,0.0095},
					{0.0061,0.0062,0.0062,0.0062,0.0064,0.0065,0.0066,0.0068,0.0071,0.0077,0.0092},
					{0.0059,0.0060,0.0061,0.0061,0.0062,0.0063,0.0065,0.0067,0.0070,0.0077,0.0091},
					{0.0060,0.0061,0.0062,0.0062,0.0063,0.0063,0.0066,0.0067,0.0070,0.0075,0.0089},
					{0.0060,0.0060,0.0060,0.0061,0.0062,0.0063,0.0064,0.0067,0.0069,0.0075,0.0088},
					{0.0060,0.0061,0.0061,0.0062,0.0062,0.0063,0.0064,0.0067,0.0068,0.0073,0.0086},
					{0.0058,0.0058,0.0059,0.0059,0.0060,0.0061,0.0063,0.0065,0.0067,0.0073,0.0085},
					{0.0058,0.0059,0.0060,0.0060,0.0060,0.0061,0.0063,0.0065,0.0067,0.0071,0.0083},
					{0.0059,0.0059,0.0060,0.0060,0.0061,0.0061,0.0063,0.0064,0.0066,0.0070,0.0081},
					{0.0059,0.0060,0.0061,0.0061,0.0061,0.0062,0.0063,0.0064,0.0065,0.0069,0.0079},
					{0.0056,0.0057,0.0057,0.0057,0.0058,0.0059,0.0060,0.0062,0.0065,0.0069,0.0081},
					{0.0056,0.0057,0.0057,0.0057,0.0058,0.0059,0.0060,0.0062,0.0064,0.0068,0.0079},
					{0.0057,0.0057,0.0057,0.0058,0.0058,0.0059,0.0060,0.0061,0.0063,0.0067,0.0077},
					{0.0057,0.0057,0.0057,0.0058,0.0058,0.0059,0.0060,0.0061,0.0063,0.0066,0.0076},
					{0.0057,0.0057,0.0058,0.0058,0.0058,0.0059,0.0059,0.0060,0.0062,0.0065,0.0074},
					{0.0057,0.0058,0.0058,0.0058,0.0059,0.0059,0.0059,0.0060,0.0061,0.0064,0.0073},
					{0.0054,0.0054,0.0054,0.0055,0.0055,0.0056,0.0058,0.0059,0.0061,0.0065,0.0076},
					{0.0054,0.0054,0.0054,0.0054,0.0055,0.0056,0.0057,0.0059,0.0060,0.0065,0.0075},
					{0.0054,0.0053,0.0054,0.0054,0.0055,0.0055,0.0057,0.0058,0.0060,0.0064,0.0074},
					{0.0053,0.0054,0.0054,0.0054,0.0055,0.0055,0.0056,0.0058,0.0059,0.0063,0.0073},
					{0.0053,0.0054,0.0054,0.0054,0.0054,0.0055,0.0056,0.0057,0.0059,0.0062,0.0071},
					{0.0053,0.0054,0.0054,0.0054,0.0054,0.0055,0.0056,0.0057,0.0058,0.0061,0.0070},
					{0.0053,0.0053,0.0054,0.0054,0.0054,0.0054,0.0055,0.0056,0.0057,0.0061,0.0069},
					{0.0053,0.0054,0.0054,0.0054,0.0054,0.0054,0.0055,0.0056,0.0057,0.0060,0.0068},
					{0.0053,0.0054,0.0054,0.0054,0.0054,0.0054,0.0054,0.0055,0.0056,0.0059,0.0067},
					{0.0053,0.0053,0.0054,0.0054,0.0054,0.0054,0.0054,0.0055,0.0055,0.0058,0.0066},
					{0.0051,0.0051,0.0051,0.0051,0.0052,0.0053,0.0054,0.0055,0.0057,0.0060,0.0070},
					{0.0051,0.0051,0.0051,0.0051,0.0052,0.0052,0.0053,0.0055,0.0056,0.0060,0.0069},
					{0.0051,0.0051,0.0051,0.0051,0.0052,0.0052,0.0053,0.0054,0.0056,0.0059,0.0068},
					{0.0051,0.0051,0.0051,0.0051,0.0051,0.0052,0.0052,0.0054,0.0055,0.0058,0.0067},
					{0.0051,0.0051,0.0051,0.0051,0.0051,0.0051,0.0052,0.0053,0.0055,0.0058,0.0066},
					{0.0050,0.0051,0.0051,0.0050,0.0051,0.0051,0.0052,0.0053,0.0054,0.0057,0.0066},
					{0.0050,0.0050,0.0051,0.0050,0.0051,0.0051,0.0051,0.0052,0.0054,0.0056,0.0065},
					{0.0050,0.0050,0.0050,0.0050,0.0051,0.0051,0.0051,0.0052,0.0053,0.0056,0.0064},
					{0.0050,0.0050,0.0050,0.0050,0.0050,0.0050,0.0050,0.0051,0.0053,0.0055,0.0063},
					{0.0050,0.0050,0.0050,0.0050,0.0050,0.0050,0.0050,0.0051,0.0052,0.0054,0.0063},
					{0.0050,0.0049,0.0049,0.0050,0.0050,0.0051,0.0051,0.0052,0.0053,0.0057,0.0065},
					{0.0049,0.0049,0.0049,0.0050,0.0050,0.0050,0.0051,0.0052,0.0053,0.0057,0.0065},
					{0.0049,0.0049,0.0049,0.0049,0.0050,0.0050,0.0051,0.0052,0.0053,0.0056,0.0064},
					{0.0049,0.0049,0.0049,0.0049,0.0049,0.0050,0.0050,0.0051,0.0052,0.0056,0.0064},
					{0.0049,0.0048,0.0049,0.0049,0.0049,0.0050,0.0050,0.0051,0.0052,0.0056,0.0063},
					{0.0048,0.0048,0.0048,0.0049,0.0049,0.0049,0.0050,0.0050,0.0051,0.0055,0.0063},
					{0.0048,0.0048,0.0048,0.0048,0.0049,0.0049,0.0049,0.0050,0.0051,0.0055,0.0062},
					{0.0048,0.0048,0.0048,0.0048,0.0048,0.0049,0.0049,0.0050,0.0051,0.0054,0.0062},
					{0.0048,0.0047,0.0048,0.0048,0.0048,0.0048,0.0048,0.0049,0.0050,0.0054,0.0061},
					{0.0047,0.0047,0.0047,0.0048,0.0048,0.0048,0.0048,0.0049,0.0050,0.0053,0.0061}};
		
		double[] maturities = {0.5000,1.0000,1.5000,2.0000,2.5000,3.0000,3.5000,4.0000,4.5000,5.0000,
				5.5000,6.0000,6.5000,7.0000,7.5000,8.0000,8.5000,9.0000,9.5000,10.0000,
				10.5000,11.0000,11.5000,12.0000,12.5000,13.0000,13.5000,14.0000,14.5000,
				15.0000,15.5000,16.0000,16.5000,17.0000,17.5000,18.0000,18.5000,19.0000,
				19.5000,20.0000,20.5000,21.0000,21.5000,22.0000,22.5000,23.0000,23.5000,
				24.0000,24.5000,25.0000,25.5000,26.0000,26.5000,27.0000,27.5000,28.0000,
				28.5000,29.0000,29.5000};
		
		int[] indeces = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}; //chose the index of the maturities we want to calibrate
		int[] strikeIndeces = {1,2,3,4,5,6,7,8}; //choose the index of the strikes we want to calibrate
		
		double[] calibratedStrikes = new double[strikeIndeces.length];
		
		for(int j = 0; j <strikeIndeces.length; j++) {
			calibratedStrikes[j] = strikes[strikeIndeces[j]];
		}
		
		
		//Target caplet data
		QuotingConvention convention = QuotingConvention.VOLATILITYNORMAL;
		
		ArrayList<CapletSmileData> smiles = new ArrayList<CapletSmileData>();
		
		for(int i = 0; i<indeces.length; i++) {
			int index = indeces[i];
			double[] values = new double[strikeIndeces.length];
			
			for(int j = 0; j <strikeIndeces.length; j++) {
				values[j] = volatilities[index][strikeIndeces[j]];
			}
			
			double time = maturities[index];
			
			String forwardCurve;
			
			if(time < 2.0) {
				forwardCurve = "forward-EUR-3M";
			} else {
				forwardCurve = "forward-EUR-6M";
			}
			
			CapletSmileData ithSmile = new CapletSmileData(forwardCurve,"discount-EUR-OIS",referenceDate,calibratedStrikes,time, values, convention);
			smiles.add(ithSmile);
		}
		
		CapletSmileData[] smileArray = new CapletSmileData[smiles.size()];
		smileArray = smiles.toArray(new CapletSmileData[smiles.size()]);
		
		CapletSurfaceData surface =  new CapletSurfaceData(smileArray, calibratedModel);
		
		OptimizerFactoryInterface optimizerFactory = new OptimizerFactoryLevenbergMarquardt(400 /* maxIterations */, 8 /* maxThreads */);

		
		double[] initialParameters = new double[] {0.1, 0.2,0.1, 0.3,0.1,0.7,0.4,0.6,0.1,0.2,0.1} /* initialParameters */;
		double[] parameterStep = new double[] {0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01, 0.01} /* parameterStep */;
		
		QuantizationMultiCurveSmileOfCaplets pricer = new QuantizationMultiCurveSmileOfCaplets("forward-EUR-6M", 5.0, strikes);
		
		/*
		 * Parameters here are in the same order
		 * as returned by getParameters in FlowOfTemperedStableCBIProcess
		 */
		double b = 0.3;
		double sigma = 0.5;
		double eta =  0.5;
		double zeta = 0.2;
		double alpha = 1.8;
		double[] initialValues = {0.0012256627870067618, 0.00148606411515568};
		double[] immigrationRates = {0.02, 0.04};
		double[] lambda = {1.0, 1.0};
		
		double[] tenorLengths = {0.25,0.5};
		String[] tenorNames = {"3M", "6M"};
		
		double timeHorizon = 10.0;		
		int numberOfTimeSteps = 500;
		
		ScalarParameterInformationInterface[] lambdaInfo = 
			{new ScalarParameterInformation(false, new PositivityConstraint()),
					new ScalarParameterInformation(false, new PositivityConstraint())};
		ScalarParameterInformationInterface[] immigrationRatesInfo =
			{new ScalarParameterInformation(true, new BoundConstraint(1E-5, 1E-3)),
					new ScalarParameterInformation(true, new BoundConstraint(1E-4, 5*1E-3))};
		
		ScalarParameterInformationInterface bInfo = new ScalarParameterInformation(true, new PositivityConstraint());
		ScalarParameterInformationInterface sigmaInfo = new ScalarParameterInformation(true, new PositivityConstraint());
		ScalarParameterInformationInterface etaInfo = new ScalarParameterInformation(true, new PositivityConstraint());
		ScalarParameterInformationInterface zetaInfo = new ScalarParameterInformation(true, new PositivityConstraint());
		ScalarParameterInformationInterface alphaInfo = new ScalarParameterInformation(true, new BoundConstraint(1.0, 2.0));
		
		ScalarParameterInformationInterface[] initialValuesInfo = 
			{new ScalarParameterInformation(true, new BoundConstraint(1E-5, 1E-3)),
					new ScalarParameterInformation(true, new BoundConstraint(1E-4, 3* 1E-3))};
		
		boolean functionVConstraint = false;
		boolean expMomentConstraint = true;
		
		CBIProcessInterface cbiProcess = new FlowOfTemperedAlphaStableCBIprocess(timeHorizon, numberOfTimeSteps, initialValues, immigrationRates,  b, sigma, 
				 eta,  zeta,  alpha,lambda, 
				lambdaInfo, immigrationRatesInfo, 
				bInfo, sigmaInfo, etaInfo, 
				zetaInfo, alphaInfo, initialValuesInfo,
				functionVConstraint, expMomentConstraint);
		
		CBIDrivenMultiCurveModel model = new CBIDrivenMultiCurveModel(calibratedModel, cbiProcess, tenorLengths, tenorNames);
		
		int levelOfQuantization = 10;
		
		CapletCalibrationProblemForQuantization problem = new CapletCalibrationProblemForQuantization(surface, model, optimizerFactory, pricer, levelOfQuantization, initialParameters, parameterStep,CapletCalibrationProblemForQuantization.Error.RMSE);
		
		System.out.println("Calibration started");
		
		long startMillis	= System.currentTimeMillis();
		CapletCalibrationProblemForQuantization.OptimizationResult result = problem.runCalibration();
		long endMillis		= System.currentTimeMillis();
		
		double calculationTime = ((endMillis-startMillis)/1000.0);
		
		System.out.println("Calibration completed in: " +calculationTime + " seconds");
		
		System.out.println("The solver required " + result.getIterations() + " iterations.");
		System.out.println("RMSQE " +result.getRootMeanSquaredError());
		
		double[] parameters = result.getModel().getParameters();
		for(int i =0; i<parameters.length; i++) {
			System.out.println(parameters[i]);
		}
		
		ArrayList<String> errorsOverview = result.getCalibrationOutput();
		
		for(String myString : errorsOverview)
			System.out.println(myString);
		
		
		System.out.println("Finished.");
		
	}

}
