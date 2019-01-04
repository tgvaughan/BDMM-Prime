package bdmm.distributions;

/**
 * Created by Jeremie Scire (jscire)
 */

/**
 * Class containing the values of P0 and Ge.
 */
public class P0GeState extends P0State {
	
	int dimension;
	public SmallNumber[] ge;
	public double[] p0;

	public P0GeState(int nTypes) {
	    super(nTypes);
		ge = new SmallNumber[nTypes];
		for (int i = 0; i<nTypes; i++)
		    ge[i] = new SmallNumber();

    }
	
	public P0GeState(double[] p0, SmallNumber[] ge) {
        super(p0);
		if(p0.length != ge.length) {
			throw new RuntimeException("Incorrect initialization: difference of size between ge and p0");
		}
		this.ge = ge;
	}
	
	public P0GeState() {
	    super();
		ge = new SmallNumber[] {new SmallNumber()};
	}

    @Override
    public String toString() {
	    StringBuilder sb = new StringBuilder();

        for (int type=0; type<dimension; type++) {
            if (type>0)
                sb.append(" ");

            sb.append("p0[").append(type).append("]=").append(p0[type]);
            sb.append(" ge[").append(type).append("]=").append(ge[type]);
        }

        return sb.toString();
    }

    /**
	 * TODO: change the comments here
	 * The maximal value a double can take is 1.80E308, the minimal value is 4.9E-324.
	 * With a safety margin, two doubles can be dealt with by Java when there is 630 orders of magnitude (2040 in base 2) between them.
	 */
	final static int exponentMaxValueDouble = 1023;
	final static int exponentMinValueDouble = -1022;
	final static int safeGapMinMaxDouble = 2040;

	/**
	 * Determine the appropriate scale factor(sampling) and perform the multiplication by said scale factor.
	 * This "factor" is defined as the magnitude in base 2 by which the original values are increased.
	 * They are increased so as to fit into the window of values accepted as doubles by java.
	 * A parameter 'eqtype' is introduced, so that different rules for choosing a scale factor can be introduced for each type of equations dealt with.
	 * For now, whether we are dealing with only p equations or g and p equations, the rules are basically the same for each ODE system.
	 *
	 * WARNING: If the range of values in input is bigger than the range between the maximal and minimal values (in absolute value) represented in the java 'double' type, an approximation is needed.
	 * 			The smallest values are set to zero until all values left can fit in the window allowed by the type 'double'.
	 *
	 * @return
	 */
	public ScaledNumbers getScaledState() {
		int n = p0.length;
		SmallNumber[] geConditions = ge;

		// scalingFactors will store the scaling factor chosen for the ge equations in the array 'geConditions'.
		int scalingFactor = 0;
		double[] scaledEquation = new double[2*n];

		// the first half of 'scaledEquations' contains the initial conditions for p equations, no scaling process is needed there
		for (int i=0; i<n; i++){
			scaledEquation[i] = p0[i];
		}

		if (geConditions.length > 0) {

			// initialization of minExponent and maxExponent with geConditions[idx] =0 would cause issues with further determination of the scale factor. So, we go look at the first value that is not zero, if it exists.
			int idx = 0;
			while(idx < (n-1)){
				if (geConditions[idx].getMantissa()==0) idx++;
				else break;
			}
			int maxExponent = geConditions[idx].getExponent();
			int minExponent = geConditions[idx].getExponent();

			// look for the highest and lowest orders of magnitude for values in 'equation'
			if (n > idx){
				for (int i=idx; i< n; i++) {

					// only non-zero numbers are taken into account
					if (geConditions[i].getMantissa()!=0) {
						if (geConditions[i].getExponent() > maxExponent) {
							maxExponent = geConditions[i].getExponent();
						} else if (geConditions[i].getExponent() < minExponent) {
							minExponent = geConditions[i].getExponent();
						}
					}
				}
			}

			// if the range of values in the initial conditions - 'geConditions' - does not exceed the size of the window of values authorized by 'double' type,
			// a scale factor is chosen and all input values will see their order of magnitude increased by this factor.
			if ((maxExponent - minExponent)< safeGapMinMaxDouble) {

				// if possible (if the highest value is not too close to the highest value allowed in a 'double'), for simplicity, the scale factor chosen is the opposite of the exponent of the lowest value of the array.
				// in most cases, this is enough.
				if ((maxExponent - minExponent) < exponentMaxValueDouble) scalingFactor = - minExponent;

				// else, the scale factor chosen will center the range of values of the array 'equation' in the window of values authorized by java for numbers of type 'double'.
				else scalingFactor = exponentMinValueDouble + (safeGapMinMaxDouble - (maxExponent - minExponent))/2 -minExponent;

				// finally, store in scaledEquation the corresponding numbers, increased by scalingFactor orders of magnitude.
				// scaledEquation[] is of type double[]
				for (int i=0; i<n;i++)
					scaledEquation[i+n] = multiplyByPowerOfTwo(geConditions[i].getMantissa(), geConditions[i].getExponent() + scalingFactor);



				// else, then it is impossible to fit all values from the array in the window allowed for 'double' type.
				// As an approximation, smallest values are then set to zero.
				// NOTE: The implementation below is computationally expensive and naive. However, this part deals with an extreme case, that would hardly ever occur.
				// Keeping it naive and expensive allows for simplification of the much more frequent case (above).
			} else {
				// work on a copy of geConditions[]
				SmallNumber[] eqcopy = new SmallNumber[n];

				for (int i = 0; i< n; i++ ){
					eqcopy[i] = geConditions[i];
				}
				while ((maxExponent - minExponent) >= safeGapMinMaxDouble) {

					// set smallest values to zero
					for (int i=0; i< n; i++) {
						if (eqcopy[i].getExponent() == minExponent) {
							eqcopy[i] = new SmallNumber(0);
						}
					}

					// re-initialize the value of minExponent
					idx = 0;
					while(eqcopy[idx].getMantissa()==0 && (idx < n)){
						idx++;
					}
					minExponent = eqcopy[idx].getExponent();

					for (int i=idx; i< n; i++) {

						// only non-zero numbers are taken into account
						if (eqcopy[i].getMantissa()!=0) {
							if (eqcopy[i].getExponent() < minExponent) {
								minExponent = eqcopy[i].getExponent();
							}
						}
					}
				}

				// then, choose the proper scale factor
				if ((maxExponent - minExponent) < exponentMaxValueDouble) scalingFactor = - minExponent;
				else scalingFactor = exponentMinValueDouble + (safeGapMinMaxDouble - (maxExponent - minExponent))/2 -minExponent;

				// finally, store in scaledEquation the corresponding numbers, increased by scalingFactor orders of magnitude.
				// scaledEquation[] is of type double[]
				for (int i=0; i<n;i++)
					scaledEquation[i+n] = multiplyByPowerOfTwo(eqcopy[i].getMantissa(), eqcopy[i].getExponent() + scalingFactor);

			}
		}

		return new ScaledNumbers(scalingFactor, scaledEquation);
	}

	/**
	 * Retrieve values of accurate magnitude from the 'scaled' ones.
	 * @param numbers
	 * @param factor
	 */
	public void setFromScaledState(double[] numbers, int factor){

		if (numbers.length % 2 == 1 )
			throw new RuntimeException("input should be of even length");


		System.arraycopy(numbers, 0, p0, 0, p0.length);

		for (int i = 0; i < p0.length; i++){
			ge[i] = new SmallNumber(numbers[i+p0.length]);
			ge[i].incrementExponent(-factor);
		}
	}

	/**
	 * Helper method that multiplies a double x by 2^n
	 * When n is 'small', especially when less than 30, this method is significantly faster than x*Math.pow(2, n)
	 * @param x
	 * @param n
	 * @return x * 2^n
	 */
	private double multiplyByPowerOfTwo(double x, int n){
		if(x !=0) {
			if(n>180 || n<0) // the threshold of 180 was chosen empirically to maximize the method's sampling speed
				return x*Math.pow(2, n);
			else {
				while(n>30) { //the number 30 comes from the max number of bits that can be left-shifted on an int and still be equivalent to a 2 to the power of n operation
					x = x*(1<<30);
					n-=30;
				}
				return x*(1<<n); // the remaining part of the original power of 2 with which the mantissa must be multiplied
			}
		}
		return 0;
	}
}
