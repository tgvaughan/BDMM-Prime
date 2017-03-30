package beast.math;

import beast.core.Description;

/**
 * Created by Jeremie Scire (jscire)
 */

@Description("This type contains both sets of initial conditions for both equation types: p and ge")
public class p0ge_InitialConditions {
	
	int dimension;
	public SmallNumber[] conditionsOnG;
	public double[] conditionsOnP;
	
	public p0ge_InitialConditions(double[] pcond, SmallNumber[] gcond) {
		if(pcond.length != gcond.length) {
			throw new RuntimeException("Incorrect initialization: difference of size between conditionsOnG and conditionsOnP");
		}
		dimension = pcond.length;
		conditionsOnP = pcond;
		conditionsOnG = gcond;
	}
	
	public p0ge_InitialConditions() {		
		dimension = 1;
		conditionsOnP = new double[] {0};
		conditionsOnG = new SmallNumber[] {new SmallNumber()};
	}
	
	public double[] getConditionsOnP(){
		return this.conditionsOnP;
	}
	
	public SmallNumber[] getConditionsOnG(){
		return this.conditionsOnG;
	}

}
