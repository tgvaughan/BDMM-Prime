package math;

import beast.core.Description;

@Description("Set of initial conditions for a system of ODEs, all increased by the same factor to prevent underflowing (which appears with numbers too small to be distinguished from zero, when using the type 'double')")
public class ScaledNumbers {

	// scale factor(s)
	private int[] factor;

	// set(s) of initial conditions on which the scale factor(s) was/were applied
	private double[] equation;

	public ScaledNumbers(int[] f, double[] e){
		this.factor = f;
		this.equation = e;
	}

	public ScaledNumbers(){
		this.factor = new int[0];
		this.equation = new double[0];
	}

	public double[] getEquation(){
		return this.equation;
	}

	public int[] getScalingFactor(){
		return this.factor;
	}
}
