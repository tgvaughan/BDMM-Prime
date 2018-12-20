package bdmm.distributions;

import beast.core.Description;

/**
 * Created by Jeremie Scire (jscire)
 */

@Description("This type contains both sets of initial conditions for both equation types: p and ge")
public class P0GeState {
	
	int dimension;
	public SmallNumber[] conditionsOnG;
	public double[] conditionsOnP;
	
	public P0GeState(double[] pcond, SmallNumber[] gcond) {
		if(pcond.length != gcond.length) {
			throw new RuntimeException("Incorrect initialization: difference of size between conditionsOnG and conditionsOnP");
		}
		dimension = pcond.length;
		conditionsOnP = pcond;
		conditionsOnG = gcond;
	}
	
	public P0GeState() {
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
