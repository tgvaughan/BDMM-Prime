package beast.math;

import beast.core.Description;

@Description("SmallNumberScaler contains methods to perform scaling/unscaling on a set of Small Numbers.")
public class SmallNumberScaler {

	/**
	 * TODO: change the comments here
	 * The maximal value a double can take is 1.80E308, the minimal value is 4.9E-324.
	 * With a safety margin, two doubles can be dealt with by Java when there is 630 orders of magnitude between them.
	 */
	final static int exponentMaxValueDouble = 1023;
	final static int exponentMinValueDouble = -1022;
	final static int safeGapMinMaxDouble = 2040;

	/**
	 * Determine the appropriate scale factor(s) and perform the multiplication by said scale factor.
	 * This "factor" is defined as the magnitude in base 10 by which the original values are increased.
	 * They are increased so as to fit into the window of values accepted as doubles by java.
	 * A parameter 'eqtype' is introduced, so that different rules for choosing a scale factor can be introduced for each type of equations dealt with.
	 * For now, whether we are dealing with only p equations or g and p equations, the rules are basically the same for each ODE system.
	 * 
	 * WARNING: If the range of values in input is bigger than the range between the maximal and minimal values (in absolute value) represented in the java 'double' type, an approximation is needed.
	 * 			The smallest values are set to zero until all values left can fit in the window allowed by the type 'double'.
	 * 
	 * @param conditions
	 * @return
	 */
	public static ScaledNumbers scale(p0ge_InitialConditions conditions) {
		if (conditions == null)
			throw new RuntimeException("Incorrect input (null) in method scale");

		int n = conditions.getConditionsOnP().length;
		SmallNumber[] geConditions = conditions.getConditionsOnG();

		// scalingFactors will store the scaling factor chosen for the ge equations in the array 'geConditions'.
		int scalingFactor = 0;
		double[] scaledEquation = new double[2*n];

		// the first half of 'scaledEquations' contains the initial conditions for p equations, no scaling process is needed there
		for (int i=0; i<n; i++){
			scaledEquation[i] = conditions.getConditionsOnP()[i];
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
				for (int i=0; i<n;i++){
					double a = geConditions[i].getMantissa();
					int b = geConditions[i].getExponent() + scalingFactor;
					if(a!=0) {
						if(b>180 || b<0)
							scaledEquation[i+n] = a*Math.pow(2, b);
						else {
							while(b>30) { //the number 30 comes from the max number of bits that can be left-shifted on an int and still be equivalent to a 2 to the power of n operation
								a = a*(1<<30);
								b-=30;
							}
							a = a*(1<<b); // the remaining part of the original power of 2 with which the mantissa must be multiplied
							scaledEquation[i+n] = a;
						}
					} else {
						scaledEquation[i+n]=0;
					}
				}


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
				for (int i=0; i<n;i++){
						double a = eqcopy[i].getMantissa();
						int b = eqcopy[i].getExponent() + scalingFactor;
						if(a!=0) {
							//scaledEquation[i+n] = eqcopy[i].getMantissa()*Math.pow(2, (eqcopy[i].getExponent() + scalingFactor)); TO DO REMOVE LINE
							if(b>180 || b<0)
								scaledEquation[i+n] = a*Math.pow(2,b);
							else {
								while(b>30) { //the number 30 comes from the max number of bits that can be left-shifted on an int and still be equivalent to a 2 to the power of n operation
									a = a*(1<<30);
									b-=30;
								}
								a = a*(1<<b); // the remaining part of the original power of 2 with which the mantissa must be multiplied
								scaledEquation[i+n] = a;
							}
						} else {
							scaledEquation[i+n]=0;
						}
					}
				}
			}

			return new ScaledNumbers(scalingFactor, scaledEquation);
		}

		/**
		 * Retrieve values of accurate magnitude from the 'scaled' ones.
		 * @param numbers
		 * @param factor
		 * @return an instance of p0ge_InitialConditions containing an array of SmallNumbers (ge equations) and an array of doubles (p equations)
		 */
		public static p0ge_InitialConditions unscale(double[] numbers, int factor){

			if (numbers.length % 2 == 1 )
				throw new RuntimeException("input should be of even length");

			int dim = numbers.length/2;

			double[] pConditions = new double[dim];

			System.arraycopy(numbers, 0, pConditions, 0, dim);

			SmallNumber[] unscaledGe = new SmallNumber[dim];

			for (int i = 0; i < dim; i++){
				unscaledGe[i] = new SmallNumber(numbers[i+dim]);
				unscaledGe[i].addExponent(-factor);
			}

			return new p0ge_InitialConditions(pConditions, unscaledGe);
		}
	}

