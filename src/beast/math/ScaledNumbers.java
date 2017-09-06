package beast.math;

import beast.core.Description;

/**
 * Created by Jeremie Scire (jscire) on 24.06.16.
 */

@Description("Set of initial conditions for a system of ODEs, all increased by the same factor to prevent underflowing (which appears with numbers too small to be distinguished from zero, when using the type 'double')")
public class ScaledNumbers {

	// scale factor(s)
	private int factor;

	// set(s) of initial conditions on which the scale factor(s) was/were applied
	private double[] equation;

	public ScaledNumbers(int f, double[] e){
		this.factor = f;
		this.equation = e;
	}

	public ScaledNumbers(){
		this.factor = 0;
		this.equation = new double[0];
	}

	public double[] getEquation(){
		return this.equation;
	}
	
	public void setEquation(double[] newEquation){
		this.equation = newEquation;
	}

	public int getScalingFactor(){
		return this.factor;
	}
	
	public void setScalingFactor(int newFactor){
		this.factor = newFactor; 
	}
	
	public void augmentFactor(int increase) {
		factor += increase;
	}
}
