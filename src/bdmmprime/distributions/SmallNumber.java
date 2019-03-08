package bdmmprime.distributions;

import java.text.DecimalFormat;

import beast.core.Description;

/**
 * Created by Jeremie Scire (jscire) on 24.06.16.
 */

@Description("This class contains the tools needed to represent and do basic " +
        "calculations with numbers in extended floating point representation (EPFP).")
public class SmallNumber {

	private double mantissa;
	private int exponent;
	public static final int numbersPrinted = 12;
	public static final int threshold = -4;

	// number of orders of magnitude between two Small Numbers needed to consider that the lower one is negligible compared to the higher one
	public static final int approximationThreshold = 53;

	final static int exponentMaxValueDouble = 1023;
	final static int exponentMinValueDouble = -1022;

	public SmallNumber(){
		this.mantissa = 0;
		this.exponent = 0;
	}

	public SmallNumber(double r, int exp){
		if (Double.isInfinite(r))
			throw new RuntimeException("Unauthorized number (Infinity) used for conversion into SmallNumber");
		if (r == 0) {
			this.mantissa = 0;
			this.exponent = 0;
		} else {
			this.mantissa = r;
			this.exponent = exp;
		}
		this.update();
	}

	/**
	 * creates a SmallNumber from a double
	 * @param num
	 */
	public SmallNumber(double num){
		if (Double.isInfinite(num))
			throw new RuntimeException("Unauthorized number (Infinity) used for conversion into SmallNumber");
		if (num == 0){
			this.mantissa = 0;
			this.exponent = 0;
		} else {
			int numExponent = Math.getExponent(num);


			if(numExponent > 0) {
				// instead of doing nothing in case of positive exponent, we decide to scale it as well
				this.mantissa = num *Math.pow(2, -numExponent);
				this.exponent = numExponent;
			} else {
				
				// similar approach as in method SmallNumberScaler.multiplyByPowerOfTwo()
				if(-numExponent>180)
					this.mantissa = num *Math.pow(2, -numExponent);
				else {
					int temporaryExponent = numExponent;
					while(-temporaryExponent>30) {
						num = num*(1<<30);
						temporaryExponent+=30;
					}
					num = num*(1<<(-temporaryExponent));
					this.mantissa = num;
				}
				
				this.exponent = numExponent;

			}
		}		
	}

	/**
	 * Make sure a number is in proper SmallNumber representation (mantissa not too large, not too small), if not make it so.
	 */
	public void update(){
		if (Double.isInfinite(this.mantissa))
			throw new RuntimeException("Unauthorized number (Infinity) used for conversion in SmallNumber");
		if (this.mantissa == 0) {
			this.exponent = 0;

			// if the exponent is too low now, apply 
		} else {
			int tempExp = Math.getExponent(this.mantissa);

			if(Math.abs(tempExp) > 200){ //arbitrary threshold set at 200 to avoid underflow and overflow
				this.mantissa  *= Math.pow(2, -tempExp);
				this.exponent += tempExp;
			}
		}
	}

	public int getExponent(){
		return this.exponent;
	}

	public double getMantissa(){
		return this.mantissa;
	}

    /**
     * Return a new SmallNumber which is the result of multiplying this number
     * by multiplicand.
     *
     * @param other number to multiply this number by.
     * @return result of multiplication.
     */
	public SmallNumber multiplyBy(SmallNumber other){
		if (mantissa == 0 || other.mantissa == 0){
			return new SmallNumber(0);
		}
		SmallNumber result = new SmallNumber(mantissa * other.mantissa, exponent + other.exponent);

		// use update to make sure the result is in the correct representation
		result.update();
		return result;
	}
	
	/**
	 * Multiply a SmallNumber with a double
	 * @param lambda
	 * @return a new SmallNumber
	 */
	public SmallNumber scalarMultiplyBy(double lambda){
		if (Double.isInfinite(lambda))
			throw new RuntimeException("Unauthorized number (Infinity) used for multiplication with a SmallNumber");
		SmallNumber res = new SmallNumber(mantissa * lambda, exponent);
		res.update();
		return res;
	}

	/**
	 * Increase the value of a SmallNumber by 'exp' orders of magnitude (in base 2)
	 * @param exp
	 */
	public void incrementExponent(int exp){
		this.exponent += exp;
		this.update();
	}

	/**
	 * Return the negative of a SmallNumber
	 * @return a new SmallNumber
	 */
	public SmallNumber negate(){
		return new SmallNumber(-this.mantissa, this.exponent);
	}

    /**
     * Return a new SmallNumber which is the result of adding this
     * small number to otherNumber.
     *
     * @param other number to add to this one.
     * @return result of addition.
     */
	public SmallNumber addTo(SmallNumber other){

		if (mantissa == 0 || ((other.exponent - exponent)> approximationThreshold && other.mantissa !=0)) {
			return other;
		} else if (other.mantissa == 0 || (exponent - other.exponent) > approximationThreshold) {
			return this;
		} else {
			if (exponent > other.exponent) {
				return new SmallNumber(mantissa + other.mantissa * Math.pow(2, other.exponent-exponent), exponent);
			} else {
				return new SmallNumber(other.mantissa + mantissa * Math.pow(2, exponent-other.exponent), other.exponent);
			}	
		}
	}

	/**
	 * If possible, return the number of interest in a 'double'
	 * WARNING: if the number stored in this SmallNumber is too low to be accurately represented in a double, revert will return 0.
	 * @return a double
	 */
	public double revert(){
		if(this.exponent < exponentMinValueDouble)
			return 0;

		return this.mantissa*Math.pow(2, this.exponent);
	}

	public SmallNumber getScientificRepresentation() {

		if (mantissa == 0) {
			return this;
		}
		SmallNumber a10 = new SmallNumber(0);
		// Transformation of a * 2^(birth) in alpha * c * 10^(beta+z) where a = c * 10^z and 2^birth = alpha * 10^beta
		double exponentBase10  = exponent * Math.log(2)/Math.log(10);

		int beta = (int) Math.floor(exponentBase10);
		double alpha = Math.pow(10, exponentBase10 - beta);

		double sign = Math.signum(mantissa);
		double log = Math.log10(Math.abs(mantissa));
		int z = (int)Math.floor(log);
		double c = sign * Math.pow(10, log - z);

		a10.mantissa = alpha * c;
		a10.exponent = beta+z;
		return a10;
	}


	/**
	 * Convert an array of doubles into an array of SmallNumbers 
	 * @param numbers
	 * @return
	 */
	public static SmallNumber[] getSmallNumbers(double[] numbers){
		SmallNumber[] smallNums = new SmallNumber[numbers.length];
		for (int i=0; i<numbers.length; i++){
			smallNums[i] = new SmallNumber(numbers[i]);
		}
		return smallNums;
	}

	public String toString() {

		SmallNumber num10 = getScientificRepresentation();
		String pattern = "0.";
		for (int i=1; i<numbersPrinted; i++) {
			pattern += "#";
		}
		DecimalFormat dF  = new DecimalFormat(pattern);
		return "" + dF.format(num10.mantissa) + "E" + num10.exponent;
	}

	/**
	 * toString method for an array of SmallNumbers
	 * @param nums
	 * @return
	 */
	public static String toString(SmallNumber[] nums) {

		String result = "";
		for (SmallNumber sn: nums){
			result = result + sn.toString() + "\t";
		}
		return result;
	}

	/**
	 * Convert an array of SmallNumbers into an array of doubles.
	 * To be used when the risk of underflowing is negligible.
	 * WARNING: if the Small Number number stored in one of the boxes is too low to be accurately represented in a double, it will return as 0.
	 * @param nums
	 * @return
	 */
	public static double[] getDoubles(SmallNumber[] nums){
		double[] res = new double[nums.length];
		for (int i =0; i<nums.length; i++){
			res[i] = nums[i].revert();
		}
		return res;
	}

	/**
	 * @return the log of a Small Number
	 */
	public double log(){
		if (this.mantissa <= 0) 
			return Double.NEGATIVE_INFINITY;
		return Math.log(this.mantissa) + this.exponent*Math.log(2);
	}

	public static boolean isInfinite(SmallNumber a){
		return (Double.isInfinite(a.mantissa));
	}

	/**
	 * Average exponent of an array of Small Numbers. Can be useful for testing properties of the likelihood input and results.
	 * @return
	 */
	public static double averageExponent(SmallNumber[] nums){
		double res = 0;
		int n = nums.length;

		for (int i=0; i<n; i++) {
			res = res + (nums[i].exponent*1.0)/n;
		}
		return res;
	}

	/**
	 * Compare the minimal exponent of an array of Small Numbers with an arbitrary threshold. 
	 * @param nums
	 * @return 1 if the minimal exponent went over the threshold (in absolute value), 0 if it did not
	 */
	public static int compareExponent(SmallNumber[] nums){
		int min = nums[0].exponent;
		for (int i=1; i<nums.length; i++) {
			min = min>nums[i].exponent? nums[i].exponent: min;
		}
		int res = min < SmallNumber.threshold ? 1: 0;
		return res;
	}

	public static void main(String args[]) {

		//Simple test 
		
		double testedA =  2.2914985084252684e90;
		SmallNumber snA = new SmallNumber(testedA);
		System.out.println(snA);

		double aOld = 1.2364544564e-75;
		SmallNumber a = new SmallNumber(aOld);
		System.out.println("The value of a is " + a.toString());
		SmallNumber f = new SmallNumber(Math.pow(Math.E, -50));
		System.out.println("The value of log f is " + f.log());
		SmallNumber c = new SmallNumber(0);
		double abis = a.getMantissa()*Math.pow(2, a.getExponent());
		System.out.println("The value of a is " + a.toString() + "\t vs 1,2364544564E-75");
		System.out.println("The value of c is " + c.toString());
		SmallNumber b = a
                .multiplyBy(a)
                .multiplyBy(a)
                .multiplyBy(a)
                .multiplyBy(a)
                .multiplyBy(a)
                .multiplyBy(a)
                .multiplyBy(a);
		double bOld = aOld*aOld*aOld*aOld*aOld*aOld*aOld*aOld;
		System.out.println("The value of birth is " + b.toString());
		System.out.println("The value of birth using double would have been " + bOld);


		//		SmallNumber sampling = new SmallNumber();
		//		System.out.println(sampling.toString());
		//		SmallNumber[] emptyTable = new SmallNumber[4];
		//		double emptyDouble[] = new double[4];
		//		for (int i=0; i<4; i++) emptyTable[i] = new SmallNumber();
		//		System.out.println(SmallNumber.toString(emptyTable));
		//		System.out.println(emptyDouble[2]);


		// Tests on basic operations
//				double aOld = 0.45643453;
//				double bOld = 8900.;
//				//				double bOld = 0;
//				SmallNumber a = new SmallNumber(aOld);
//				SmallNumber birth = new SmallNumber(bOld);
//				SmallNumber c = SmallNumber.multiply(a, birth);
//				System.out.println("The value of c is " + c.toString()+" vs " + aOld*bOld);
//				double lambda = 0.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000005;
//				c = c.scalarMultiply(lambda);
//				System.out.println("The value of c is " + c.toString()+" vs " + aOld*bOld*lambda);
//				SmallNumber death = SmallNumber.add(a, birth);
//				System.out.println("With SmallNumber implementation: " + death.toString() + " vs "+ (aOld+bOld));
		
		//		double aOld = 4564.3453;
		//		double bOld = 89;
		//		double cOld = aOld*bOld*lambda;
		//		double dOld = aOld+bOld;
		//		System.out.println("With classic double implementation: " + dOld);

		// Test on scaledNumbers
		double[] eqp = {0, 1, 0.5, 0.8, 0.9, 1.0, 0.6};
		SmallNumber[] eq = {new SmallNumber(0), new SmallNumber(0), new SmallNumber(1.5), new SmallNumber(0), new SmallNumber(1., 400), new SmallNumber(1., -200), new SmallNumber(1., -1000)};
		double m = SmallNumber.averageExponent(eq);

//		ScaledNumbers scaeq = SmallNumberScaler.scale(new P0GeState(eqp, eq));
//		System.out.println(SmallNumber.toString(eq) +  "with an average exponent of: " + m + "\t and minimal exponent compared to the set threshold of: " + SmallNumber.compareExponent(eq));
//		System.out.println(scaeq.getScalingFactor());
//		System.out.println("\n" + scaeq.getEquation()[0] + " " + scaeq.getEquation()[1] + " " + scaeq.getEquation()[2] + " " + scaeq.getEquation()[3] + " " + scaeq.getEquation()[4] + " " + scaeq.getEquation()[5] + " " + scaeq.getEquation()[6]);


		//		double res = 0*Math.exp(Math.log(10)*(389));
		//		System.out.println(Math.exp(2*389));
		//
		//		SmallNumber ka = new SmallNumber(Double.POSITIVE_INFINITY);
		//		System.out.println(SmallNumber.isInfinite(ka));
		//		System.out.println(ka.toString());
		//		SmallNumber kb = new SmallNumber(0);
		//		System.out.println(SmallNumber.isInfinite(kb));

	}

}
