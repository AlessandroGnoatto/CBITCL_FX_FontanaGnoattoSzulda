package org.calibrationframework.timeseries;

import java.util.Iterator;
import java.util.function.UnaryOperator;

import org.apache.commons.math3.complex.Complex;

import net.finmath.time.*;

public class FunctionV implements TimeSeriesInterface {
	
	private final Complex[] values;
	private final TimeDiscretization timeGrid;
	private final UnaryOperator<Complex> branchingMechanism;
	private final UnaryOperator<Complex> immigrationRate;
	private final UnaryOperator<Complex> levyExponent;
	private final double upperBoundD1;
	
	public FunctionV(double timeHorizon, int numberOfTimeSteps, UnaryOperator<Complex> branchingMechanism, UnaryOperator<Complex> immigrationRate, UnaryOperator<Complex> levyExponent, Complex u1, Complex u2, Complex u3, double upperBoundD1) {
		this.branchingMechanism = branchingMechanism;
		this.immigrationRate = immigrationRate;
		this.levyExponent = levyExponent;
		this.upperBoundD1 = upperBoundD1;
		double deltaT = timeHorizon / (double)(numberOfTimeSteps);
		this.timeGrid = new TimeDiscretizationFromArray(0.0, numberOfTimeSteps, deltaT);
		this.values = new Complex[numberOfTimeSteps+1];
		values[0] = u1;
		int j = 0;
		for(int i = 1; i < numberOfTimeSteps+1; i++) {
			Complex k1 = (levyExponent.apply(u3).add(u2).add(branchingMechanism.apply(values[i-1]))).multiply(deltaT);
			Complex k2 = (levyExponent.apply(u3).add(u2).add(branchingMechanism.apply(values[i-1].add(k1.multiply(deltaT*0.5))))).multiply(deltaT);
			Complex k3 = (levyExponent.apply(u3).add(u2).add(branchingMechanism.apply(values[i-1].add(k2.multiply(deltaT*0.5))))).multiply(deltaT);
			Complex k4 = (levyExponent.apply(u3).add(u2).add(branchingMechanism.apply(values[i-1].add(k3.multiply(deltaT))))).multiply(deltaT);
			values[i] = ((k1.add(k2.multiply(2)).add(k3.multiply(2)).add(k4)).divide(6.0)).add(values[i-1]);
			if(values[i].getReal() >= this.upperBoundD1) {
				values[i] = new Complex(0.0); //new Complex(this.upperBoundD1);
				j = i;
				break;
			}
		}
		for(int i = j+1; i < numberOfTimeSteps+1; i++) {
			values[i] = values[j];
		}
	}
	
	@Override
	public double getTime(int index) {
		return this.timeGrid.getTime(index);
	}

	@Override
	public double getValue(int index) {
		return this.values[index].getReal();
	}

	@Override
	public int getNumberOfTimePoints() {
		return this.values.length;
	}
	
	@Override
	public Iterable<Double> getValues() {
		return new Iterable<Double>() {
			private int index = 0;

			@Override
			public Iterator<Double> iterator() {
				return new Iterator<Double>() {
					@Override
					public boolean hasNext() {
						return index < FunctionV.this.getNumberOfTimePoints();
					}

					@Override
					public Double next() {
						return FunctionV.this.getValue(index++);
					}
				};
			}

		};
	}
	
	public Complex getValue(double time) {
		if(this.timeGrid.getTimeIndex(time) < 0) {
			return this.values[this.timeGrid.getTimeIndexNearestLessOrEqual(time)];
		} else {
			return this.values[this.timeGrid.getTimeIndex(time)];
		}
	}
	
	/**
	 * Returns function U
	 * @param firstTime
	 * @param lastTime
	 */
	public Complex getIntegralOfImmigrationRate(double firstTime, double lastTime) {
		int firstIndex = 0;
		int lastIndex = 0;
		if(this.timeGrid.getTimeIndex(firstTime) < 0) {
			firstIndex = this.timeGrid.getTimeIndexNearestLessOrEqual(firstTime);
		} else {
			firstIndex = this.timeGrid.getTimeIndex(firstTime);
		}
		if(this.timeGrid.getTimeIndex(lastTime) < 0) {
			lastIndex = this.timeGrid.getTimeIndexNearestLessOrEqual(lastTime);
		} else {
			lastIndex = this.timeGrid.getTimeIndex(lastTime);
		}
		Complex sum = new Complex(0,0);
		for(int i = firstIndex; i < lastIndex + 1; i++) {
			sum = sum.add( this.getImmigrationRate().apply( this.values[i] ).multiply( this.timeGrid.getTimeStep(i) ) );
		}
		
		return sum;
	}
	
	public UnaryOperator<Complex> getBranchingMechanism() {
		return this.branchingMechanism;
	}
	
	public UnaryOperator<Complex> getImmigrationRate() {
		return this.immigrationRate;
	}
	
	public UnaryOperator<Complex> getLevyExponent() {
		return this.levyExponent;
	}
	
}
