package beast.math;

import java.text.DecimalFormat;

import beast.core.Description;

@Description("This class contains the tools needed to represent and do basic calculations with numbers in scientific representation."
		+ " For instance, 0.0523 would be written as 5.23E-2 in scientific representation. In this implementation, the attribute 'root' would be 5.23 and 'exponent' would be -2.")
public class SmallNumber {

	private double root;
	private int exponent;
	public static final int numbersPrinted = 12;
	public static final int threshold = -4;

	// number of orders of magnitude between two Small Numbers needed to consider that the lower one is negligible compared to the higher one
	public static final int approximationThreshold = 18;

	final static int exponentMaxValueDouble = 308;
	final static int exponentMinValueDouble = -324;

	public SmallNumber(){
		this.root = 0;
		this.exponent = 0;
	}

	public SmallNumber(double r, int exp){
		if (Double.isInfinite(r))
			throw new RuntimeException("Unauthorized number (Infinity) used for conversion into SmallNumber");
		if (r == 0) {
			this.root = 0;
			this.exponent = 0;
		} else {
			this.root = r;
			this.exponent = exp;
		}
		this.update();
	}

	/**
	 * creates a SmallNumber (number in scientific representation) from a double
	 * @param num
	 */
	public SmallNumber(double num){
		if (Double.isInfinite(num))
			throw new RuntimeException("Unauthorized number (Infinity) used for conversion into SmallNumber");
		if (num == 0){
			this.root = 0;
			this.exponent = 0;
		} else {
			double sign = Math.signum(num);
			double log = Math.log10(Math.abs(num));
			int exponent = (int)Math.floor(log);
			this.root = sign * Math.pow(10, log - exponent);
			this.exponent = exponent;
		}
	}

	/**
	 * Make sure a number is in scientific representation, if not make it so.
	 */
	public void update(){
		if (Double.isInfinite(this.root))
			throw new RuntimeException("Unauthorized number (Infinity) used for conversion in SmallNumber");
		if (this.root == 0) {
			this.exponent = 0;

			// if root is lower than 1 or higher than 10 in absolute value, the exponent and root are changed to fit the scientific representation format
		} else if (Math.abs(this.root)<1 || Math.abs(this.root) >= 10){
			double sign = Math.signum(this.root);
			double log = Math.log10(Math.abs(this.root));
			int exp = (int)Math.floor(log);
			this.root = sign * Math.pow(10, log - exp);
			this.exponent += exp;
		}
	}

	public int getExponent(){
		return this.exponent;
	}

	public double getRoot(){
		return this.root;
	}

	/**
	 * Multiply two numbers in scientific representation 
	 * @param a
	 * @param b
	 * @return a new SmallNumber
	 */
	public static SmallNumber multiply(SmallNumber a, SmallNumber b){
		if (a.root == 0 || b.root == 0){
			return new SmallNumber(0);
		}
		SmallNumber result = new SmallNumber(a.root * b.root, a.exponent + b.exponent);

		// use update to make sure the result is in the correct representation
		result.update();
		return result;
	}

	/**
	 * Multiply a SmallNumber with a double
	 * @param lambda
	 * @return a new SmallNumber
	 */
	public SmallNumber scalarMultiply(double lambda){
		if (Double.isInfinite(lambda))
			throw new RuntimeException("Unauthorized number (Infinity) used for multiplication with a SmallNumber");
		SmallNumber res = new SmallNumber(this.root * lambda, this.exponent);
		res.update();
		return res;
	}

	/**
	 * Increase the value of a SmallNumber by 'exp' orders of magnitude
	 * @param exp
	 */
	public void addExponent(int exp){
		// TO DO warn against misuses of this method (if infinity for instance)
		this.exponent += exp;
		this.update();
	}

	/**
	 * Return the opposite of a SmallNumber
	 * @return a new SmallNumber
	 */
	public SmallNumber opposite(){
		return new SmallNumber(-this.root, this.exponent);
	}

	/**
	 * Add two numbers in scientific representation
	 * @param a
	 * @param b
	 * @return
	 */
	public static SmallNumber add(SmallNumber a, SmallNumber b){

		if (a.root == 0 || (b.exponent - a.exponent)> approximationThreshold) {
			return b;
		} else if (b.root == 0 || (a.exponent - b.exponent) > approximationThreshold) {
			return a;
		} else {
			SmallNumber c = new SmallNumber(0);
			if (a.exponent > b.exponent) {
				c = new SmallNumber(a.root + b.root * Math.pow(10, b.exponent-a.exponent), a.exponent); 
			} else {
				c = new SmallNumber(b.root + a.root * Math.pow(10, a.exponent-b.exponent), b.exponent);
			}	
			return c;
		}
	}

	/**
	 * If possible, return the number of interest in a 'double'
	 * WARNING: if the number stored in this SmallNumber is too low to be accurately represented in a double, revert will return 0.
	 * @return a double
	 */
	public double revert(){
		if(this.exponent < -308)
			return 0;

		return this.root*Math.pow(10,this.exponent);
	}


	public String toString(){
		String pattern = "0.";
		for (int i=1; i<numbersPrinted; i++) {
			pattern += "#";
		}
		DecimalFormat dF  = new DecimalFormat(pattern);
		return "" + dF.format(this.root) + "E" + this.exponent;
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
		if (this.root <= 0) 
			return Double.NEGATIVE_INFINITY;
		return Math.log(this.root) + this.exponent*Math.log(10);
	}


	public static boolean isInfinite(SmallNumber a){
		return (Double.isInfinite(a.root));
	}

	/**
	 * Average exponent of an array of Small Numbers. Can be useful for testing properties of the likelihood input and results.
	 * @param a
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
/*
		//Simple test 
		
		double aOld = 0.00000000000000000000000000000000000000000000000000000000000000000000000000123645445640000;
		SmallNumber a = new SmallNumber(aOld);
		SmallNumber c = new SmallNumber(0);
		System.out.println("The value of c is " + c.toString());
		SmallNumber b = SmallNumber.multiply(SmallNumber.multiply(SmallNumber.multiply(a, a), SmallNumber.multiply(a, a)),SmallNumber.multiply(SmallNumber.multiply(a, a), SmallNumber.multiply(a, a)));
		double bOld = aOld*aOld*aOld*aOld*aOld*aOld*aOld*aOld;
		System.out.println("The value of b is " + b.toString());
		System.out.println("The value of b using double would have been " + bOld);

		SmallNumber s = new SmallNumber();
		System.out.println(s.toString());
		SmallNumber[] emptyTable = new SmallNumber[4];
		double emptyDouble[] = new double[4];
		for (int i=0; i<4; i++) emptyTable[i] = new SmallNumber();
		System.out.println(SmallNumber.toString(emptyTable));
		System.out.println(emptyDouble[2]);
		 
*/
		// Tests on basic operations
		SmallNumber a = new SmallNumber(4564.3453);
		SmallNumber b = new SmallNumber(89);
		SmallNumber c = SmallNumber.multiply(a, b);
		double lambda = 0.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000005;
		c = c.scalarMultiply(lambda);
		SmallNumber d = SmallNumber.add(a, b);
		System.out.println("With SmallNumber implementation: " + d.toString());

		double aOld = 4564.3453;
		double bOld = 89;
		double cOld = aOld*bOld*lambda;
		double dOld = aOld+bOld;
		System.out.println("With classic double implementation: " + dOld);

		// Test sur scaledNumbers
		SmallNumber[] eq = {new SmallNumber(0), new SmallNumber(0), new SmallNumber(1.5), new SmallNumber(0), new SmallNumber(1., 400), new SmallNumber(1., -200), new SmallNumber(1., -500)};
		double m = SmallNumber.averageExponent(eq);
		
		ScaledNumbers scaeq = SmallNumberScaler.scale(eq, EquationType.EquationOnP);
		System.out.println(SmallNumber.toString(eq) +  "with an average exponent of: " + m + "\t and minimal exponent compared to the set threshold of: " + SmallNumber.compareExponent(eq));
		System.out.println(scaeq.getScalingFactor()[0]);
		System.out.println("\n" + scaeq.getEquation()[0] + " " + scaeq.getEquation()[1] + " " + scaeq.getEquation()[2] + " " + scaeq.getEquation()[3] + " " + scaeq.getEquation()[4] + " " + scaeq.getEquation()[5] + " " + scaeq.getEquation()[6]);

		//		double res = 0*Math.exp(Math.log(10)*(389));
		//		System.out.println(Math.exp(2*389));

		//		SmallNumber ka = new SmallNumber(Double.POSITIVE_INFINITY);
		//		System.out.println(SmallNumber.isInfinite(ka));
		//		System.out.println(ka.toString());
		//		SmallNumber kb = new SmallNumber(0);
		//		System.out.println(SmallNumber.isInfinite(kb));
		 
		


	}

}
