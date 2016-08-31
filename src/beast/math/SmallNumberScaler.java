package beast.math;

import org.apache.commons.lang3.*;

import beast.core.Description;
import beast.evolution.speciation.PiecewiseBirthDeathMigrationDistribution;

@Description("SmallNumberScaler contains methods to perform scaling/unscaling on a set of Small Numbers.")
public class SmallNumberScaler {

	/**
	 * The maximal value a double can take is 1.80E308, the minimal value is 4.9E-324.
	 * With a safety margin, two doubles can be dealt with by Java when there is 630 orders of magnitude between them.
	 */
	final static int exponentMaxValueDouble = 308;
	final static int exponentMinValueDouble = -324;
	final static int safeGapMinMaxDouble = 630;

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
	 * @param equation
	 * @param eqtype
	 * @return
	 */
	public static ScaledNumbers scale(SmallNumber[] equation, EquationType eqtype) {
		if (equation == null) 
			throw new RuntimeException("Incorrect input (null) in method scale");

		int n = equation.length;

		// scalingFactors will store the scaling factor chosen for the array 'equation'. It is an array to allow for choosing multiple scale factors, for instance when 'equation' represents initial conditions for more than one system of ODEs.
		int[] scalingFactors = null;
		double[] scaledEquation = new double[n];

		// Depending on the type of equation, different rules of scaling can be applied
		switch (eqtype)
		{
		case EquationOnP:


			// initialization of minExponent and maxExponent with equation[idx] =0 would cause issues with further determination of the scale factor. So, we go look at the first value that is not zero, if it exists.
			int idx = 0;
			while(idx < (n-1)){
				if (equation[idx].getRoot()==0) idx++;
				else break;
			}
			int maxExponent = equation[idx].getExponent();
			int minExponent = equation[idx].getExponent();

			// look for the highest and lowest orders of magnitude for values in 'equation'
			if (n > idx){
				for (int i=idx; i< n; i++) {

					// only non-zero numbers are taken into account
					if (equation[i].getRoot()!=0) {
						if (equation[i].getExponent() > maxExponent) {
							maxExponent = equation[i].getExponent();
						} else if (equation[i].getExponent() < minExponent) {
							minExponent = equation[i].getExponent();
						}
					}
				}
			}

			// if the range of values in the initial conditions - 'equation' - does not exceed the size of the window of values authorized by 'double' type,
			// a scale factor is chosen and all input values will see their order of magnitude increased by this factor.
			if ((maxExponent - minExponent)< safeGapMinMaxDouble) {

				// if possible (if the highest value is not too close from the highest value allowed in a 'double'), for simplicity, the scale factor chosen is the opposite of the exponent of the lowest value of the array.
				// in most cases, this is enough.
				if ((maxExponent - minExponent) < exponentMaxValueDouble) scalingFactors = new int[]{ - minExponent};

				// else, the scale factor chosen will center the range of values of the array 'equation' in the window of values authorized by java for numbers of type 'double'.
				else scalingFactors = new int[] { exponentMinValueDouble + (safeGapMinMaxDouble - (maxExponent - minExponent))/2 -minExponent};

				// finally, store in scaledEquation the corresponding numbers, increased by scalingFactor orders of magnitude.
				// scaledEquation[] is of type double[]
				for (int i=0; i<n;i++){
					if(equation[i].getRoot()!=0) {
						scaledEquation[i] = equation[i].getRoot()*Math.exp(Math.log(10)*(equation[i].getExponent() + scalingFactors[0]));
					} else {
						scaledEquation[i]=0;
					}
				}


				// else, then it is impossible to fit all values from the array in the window allowed for 'double' type.
				// As an approximation, smallest values are then set to zero.
				// NOTE: The implementation below is computationally expensive and naive. However, this part deals with an extreme case, that would hardly ever occur.
				// Keeping it naive and expensive allows for simplification of the much more frequent case (above).
			} else {
				// work on a copy of equation[]
				SmallNumber[] eqcopy = new SmallNumber[n];

				for (int i = 0; i< n; i++ ){
					eqcopy[i] = equation[i];
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
					while(eqcopy[idx].getRoot()==0 && (idx < n)){
						idx++;
					}
					minExponent = eqcopy[idx].getExponent();

					for (int i=idx; i< n; i++) {

						// only non-zero numbers are taken into account
						if (eqcopy[i].getRoot()!=0) {
							if (eqcopy[i].getExponent() < minExponent) {
								minExponent = eqcopy[i].getExponent();
							}
						}
					}
				}

				// then, choose the proper scale factor
				if ((maxExponent - minExponent) < exponentMaxValueDouble) scalingFactors = new int[]{ - minExponent};
				else scalingFactors = new int[] { exponentMinValueDouble + (safeGapMinMaxDouble - (maxExponent - minExponent))/2 -minExponent};

				// finally, store in scaledEquation the corresponding numbers, increased by scalingFactor orders of magnitude.
				// scaledEquation[] is of type double[]
				for (int i=0; i<n;i++){
					if(eqcopy[i].getRoot()!=0) {
						scaledEquation[i] = eqcopy[i].getRoot()*Math.exp(Math.log(10)*(eqcopy[i].getExponent() + scalingFactors[0]));
					} else {
						scaledEquation[i]=0;
					}
				}


			}



			break;


		case EquationOnGe:

			if (n % 2 == 1 )
				throw new RuntimeException("equation[] should be of even length if EquationType is EquationOnGe");
			int dim = n/2;

			// Call scale for  p equations: (0 .. dim-1)
			SmallNumber[] pequation = new SmallNumber[dim];
			for (int i = 0; i<dim; i++){
				pequation[i]=equation[i];
			}

			ScaledNumbers pscaled = SmallNumberScaler.scale(pequation, EquationType.EquationOnP);

			// Call scale with type EquationOnP also for ge equations: (dim .. 2*dim-1)
			SmallNumber[] geequation = new SmallNumber[dim];
			for (int i = 0; i<dim; i++){
				geequation[i]=equation[i+dim];
			}

			ScaledNumbers gescaled = SmallNumberScaler.scale(geequation, EquationType.EquationOnP);

			scalingFactors = new int[] {pscaled.getScalingFactor()[0], gescaled.getScalingFactor()[0]};

			scaledEquation = ArrayUtils.addAll(pscaled.getEquation(), gescaled.getEquation());
			break;
		}

		return new ScaledNumbers(scalingFactors, scaledEquation);
	}

	/**
	 * Retrieve values of accurate magnitude from the 'scaled' ones.
	 * @param numbers
	 * @param factors
	 * @param eqType
	 * @return an array of SmallNumbers
	 */
	public static SmallNumber[] unscale(double[] numbers, int[] factors, EquationType eqType){

		SmallNumber[] unscaledNumbers = new SmallNumber[numbers.length];
		

		switch(eqType){

		// Simplest case: all values were increased by the same order of magnitude.
		case EquationOnP:
			for (int i =0; i<numbers.length; i++){
				unscaledNumbers[i] = new SmallNumber(numbers[i]);
				unscaledNumbers[i].addExponent(-factors[0]);
			}
			break;

			// In this case, two different scale factors were used to deal with equations p and equations g.
		case EquationOnGe:
			if (numbers.length % 2 == 1 )
				throw new RuntimeException("numbers should be of even length if EquationType is EquationOnGe");

			if (factors.length != 2)
				throw new RuntimeException("Array factors[] is of incorrect size. Check eqtype maybe");

			int dim = numbers.length/2;
			for (int i = 0; i < dim; i++){
				unscaledNumbers[i] = new SmallNumber(numbers[i]);
				unscaledNumbers[i].addExponent(-factors[0]);
				
				unscaledNumbers[i+dim] = new SmallNumber(numbers[i+dim]);
				unscaledNumbers[i+dim].addExponent(-factors[1]);
			}
				
			break;
		}

		return unscaledNumbers;
	}
}

