package bdmm.distributions;

import org.junit.Assert;
import org.junit.Test;

public class SmallNumberTest {

    double TOLERANCE = 1e-20;

    @Test
    public void test1() {

        double testedA = 2.2914985084252684E90;
        SmallNumber snA = new SmallNumber(testedA);

        Assert.assertEquals(Math.log(testedA), snA.log(), TOLERANCE);
    }

    @Test
    public void test2() {

        double aOld = 1.2364544564e-75;
        SmallNumber a = new SmallNumber(aOld);
        SmallNumber f = new SmallNumber(Math.exp(-50));
        SmallNumber c = new SmallNumber(0);

        Assert.assertEquals(Math.log(aOld), a.log(), TOLERANCE);

        SmallNumber b = new SmallNumber(1.0)
                .multiplyBy(a)
                .multiplyBy(a)
                .multiplyBy(a)
                .multiplyBy(a)
                .multiplyBy(a)
                .multiplyBy(a)
                .multiplyBy(a)
                .multiplyBy(a);

        double trueLogB = -1379.8530719994023;
        Assert.assertEquals(trueLogB, b.log(), TOLERANCE);


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

    }

    @Test
    public void test3() {

		// Test on scaledNumbers
		double[] eqp = {0, 1, 0.5, 0.8, 0.9, 1.0, 0.6};
		SmallNumber[] eq = {
		        new SmallNumber(0),
                new SmallNumber(0),
                new SmallNumber(1.5),
                new SmallNumber(0),
                new SmallNumber(1., 400),
                new SmallNumber(1., -200),
                new SmallNumber(1., -1000)};
		double m = SmallNumber.averageExponent(eq);

		ScaledNumbers scaeq = SmallNumberScaler.scale(new p0ge_InitialConditions(eqp, eq));
		System.out.println(SmallNumber.toString(eq) +  "with an average exponent of: " + m + "\t and minimal exponent compared to the set threshold of: " + SmallNumber.compareExponent(eq));
		System.out.println(scaeq.getScalingFactor());
		System.out.println("\n"
                + scaeq.getEquation()[0] + " "
                + scaeq.getEquation()[1] + " "
                + scaeq.getEquation()[2] + " "
                + scaeq.getEquation()[3] + " "
                + scaeq.getEquation()[4] + " "
                + scaeq.getEquation()[5] + " "
                + scaeq.getEquation()[6]);

		// TODO: Implement assertions for this test (not sure what it's meant to be doing TBH!)

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
