package bdmmprime.distribution;

import bdmmprime.parameterization.*;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import org.apache.commons.math.special.Gamma;
import org.junit.Test;

import static junit.framework.Assert.assertEquals;


/**
 * Created by Jeremie Scire (jscire) on 26.06.17.
 */

public class BirthDeathMigrationLikelihoodTest {

	double runtime;

	/**
     * The original tests were developed assuming BDSKY/BDMM-like behaviour, i.e. return an oriented
	 * tree probability unless r!=1 in which case return an un-oriented and unlabeled tree probability.
	 * In contrast, BDMM-Prime always returns a labeled tree probability.
     *
	 * This method exists to convert BDSKY/BDMM test probabilities to be labeled tree probabilities,
	 * allowing comparison with BDMM-Prime.
	 *
	 * @param density BDMM-prime probability density object
	 * @return conversion factor
	 */
	private double labeledTreeConversionFactor(BirthDeathMigrationDistribution density) {
		Tree tree = (Tree)density.treeInput.get();
		boolean SAmodel = density.parameterizationInput.get().getRemovalProbs()[0][0] != 1.0;
		double factor = - Gamma.logGamma(tree.getLeafNodeCount() +1);

		if (!SAmodel)
			factor += Math.log(2) * (tree.getLeafNodeCount() - tree.getDirectAncestorNodeCount() - 1);

		return factor;
	}

	/**
	 * Basic test for migration rate change 
	 * Reference from BDMM itself
     * Canonical parameterization
	 */
	@Test 
	public void testLikelihoodMigRateChangeBasicCanonical() {

		// Test for uncoloured tree

        String newick = "(t1[&state=0] : 1.5, t2[&state=1] : 0.5);";

		RealParameter originParam = new RealParameter("2.5");

		Parameterization parameterization = new CanonicalParameterization();
		parameterization.initByName(
		        "typeSet", new TypeSet(2),
                "processLength", originParam,
                "birthRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("2.0"), 2),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0"), 2),
                "birthRateAmongDemes", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.0"), 2),
                "migrationRate", new SkylineMatrixParameter(
                        new RealParameter("1.0"),
                        new RealParameter("0.1 0.2"), 2),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.5"), 2),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0"), 2));


		BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();

		density.initByName(
		        "parameterization", parameterization,
                "frequencies", new RealParameter("0.5 0.5"),
                "conditionOnSurvival", false,
                "tree", new TreeParser(newick,
                        false, false,
                        true,0),
                "typeLabel", "state",
                "parallelize", false);

		double logL = density.calculateLogP();

		System.out.println("Birth-death result: " + logL + "\t- Test LikelihoodMigRateChange 1");

		// Reference BDMM (version 0.2.0) 22/06/2017
		assertEquals(-6.7022069383966025 - labeledTreeConversionFactor(density),
				logL, 1e-5);
	}

	/**
	 * Basic test for migration rate change
	 * Reference from BDMM itself
	 * Canonical parameterization
	 */
	@Test
	public void testLikelihoodMigRateChangeBasicCanonicalRevTime() {

		// Test for uncoloured tree

		String newick = "(t1[&state=0] : 1.5, t2[&state=1] : 0.5);";

		RealParameter originParam = new RealParameter("2.5");

		Parameterization parameterization = new CanonicalParameterization();
		parameterization.initByName(
				"typeSet", new TypeSet(2),
				"processLength", originParam,
				"birthRate", new SkylineVectorParameter(
						null,
						new RealParameter("2.0"), 2),
				"deathRate", new SkylineVectorParameter(
						null,
						new RealParameter("1.0"), 2),
				"birthRateAmongDemes", new SkylineMatrixParameter(
						null,
						new RealParameter("0.0"), 2),
				"migrationRate", new SkylineMatrixParameter(
						new RealParameter("1.5"),
						new RealParameter("0.2 0.1"), 2, originParam),
				"samplingRate", new SkylineVectorParameter(
						null,
						new RealParameter("0.5"), 2),
				"removalProb", new SkylineVectorParameter(
						null,
						new RealParameter("1.0"), 2));


		BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();

		density.initByName(
				"parameterization", parameterization,
				"frequencies", new RealParameter("0.5 0.5"),
				"conditionOnSurvival", false,
				"tree", new TreeParser(newick,
						false, false,
						true,0),
				"typeLabel", "state",
				"parallelize", false);

		double logL = density.calculateLogP();

		System.out.println("Birth-death result: " + logL + "\t- Test LikelihoodMigRateChange 1");

		// Reference BDMM (version 0.2.0) 22/06/2017
		assertEquals(-6.7022069383966025 - labeledTreeConversionFactor(density),
				logL, 1e-5);
	}

    /**
	 * Basic test for migration rate change
	 * Reference from BDMM itself
     * Epi parameterization
	 */
	@Test
	public void testLikelihoodMigRateChangeBasicEpi() {

		// Test for uncoloured tree

        String newick = "(t1[&state=0] : 1.5, t2[&state=1] : 0.5);";

		RealParameter originParam = new RealParameter("2.5");

		Parameterization parameterization = new EpiParameterization();
		parameterization.initByName(
                "typeSet", new TypeSet(2),
                "processLength", originParam,
                "R0", new SkylineVectorParameter(
                        null,
                        new RealParameter(4.0/3.0 + " " + 4.0/3.0)),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.5 1.5")),
                "R0AmongDemes", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.0 0.0")),
                "migrationRate", new SkylineMatrixParameter(
                        new RealParameter("1.0"),
                        new RealParameter("0.1 0.1 0.2 0.2")),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealParameter(1.0/3.0 + " " + 1.0/3.0)),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0 1.0")));


		BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();

		density.initByName(
		        "parameterization", parameterization,
                "frequencies", new RealParameter("0.5 0.5"),
                "conditionOnSurvival", false,
                "tree", new TreeParser(newick,
                        false, false,
                        true,0),
                "typeLabel", "state",
                "parallelize", true);

		double logL = density.calculateLogP();

		System.out.println("Birth-death result: " + logL + "\t- Test LikelihoodMigRateChange 1");

		// Reference BDMM (version 0.2.0) 22/06/2017
		assertEquals(-6.7022069383966025, logL, 1e-5);
	}

	/**
	 * Basic test for removal-probability rate change 
	 * Reference from BDMM itself
	 * @throws Exception
	 */
	@Test 
	public void testLikelihoodRemovalProbChangeBasic() {

		String newick = "((1[&type=0]: 1.5, 2[&type=0]: 0)3[&type=0]: 3.5, 4[&type=0]: 4) ;";

		Parameterization parameterization = new EpiParameterization();
		parameterization.initByName(
                "typeSet", new TypeSet(1),
                "processLength", new RealParameter("6.0"),
                "R0", new SkylineVectorParameter(
                        null,
                        new RealParameter(String.valueOf(4.0/3.0))),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.5")),
                "R0AmongDemes", new SkylineMatrixParameter(null, null),
                "migrationRate", new SkylineMatrixParameter(null, null),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealParameter(String.valueOf(1.0/3.0))),
                "removalProb", new SkylineVectorParameter(
                        new RealParameter("1.0"),
                        new RealParameter("0.3 0.7")));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
		density.initByName(
		        "parameterization", parameterization,
                "frequencies", new RealParameter("1.0"),
                "conditionOnSurvival", false,
                "tree", new TreeParser(newick, false, false, true,0),
                "typeLabel", "type",
                "parallelize", false,
				"useAnalyticalSingleTypeSolution", false);

		double logL = density.calculateLogP();

        // Reference BDMM (version 0.2.0) 22/06/2017
		assertEquals(-21.25413884159791 + labeledTreeConversionFactor(density),
				logL, 1e-5);

		BirthDeathMigrationDistribution densityExact = new BirthDeathMigrationDistribution();
		densityExact.initByName(
				"parameterization", parameterization,
				"frequencies", new RealParameter("1.0"),
				"conditionOnSurvival", false,
				"tree", new TreeParser(newick, false, false, true,0),
				"typeLabel", "type",
				"parallelize", false,
				"useAnalyticalSingleTypeSolution", true);

		double logLExact = densityExact.calculateLogP();

		// Reference BDMM (version 0.2.0) 22/06/2017
		assertEquals(-21.25413884159791 + labeledTreeConversionFactor(density),
				logLExact, 1e-5);
	}

	/**
	 * Direct comparison between numerical and analytical solutions for a tiny example with no rate changes.
	 */
	@Test
	public void tinyAnalyticalTest() {
		String newick = "(1[&type=0]: 1.0, 2[&type=0]: 1.0): 1.0;";

		Parameterization parameterization = new CanonicalParameterization();
		parameterization.initByName(
				"typeSet", new TypeSet(1),
				"processLength", new RealParameter("2.0"),
				"birthRate", new SkylineVectorParameter(
						null,
						new RealParameter("2.0")),
				"deathRate", new SkylineVectorParameter(
						null,
						new RealParameter("0.5")),
				"birthRateAmongDemes", new SkylineMatrixParameter(null, null),
				"migrationRate", new SkylineMatrixParameter(null, null),
				"samplingRate", new SkylineVectorParameter(
						null,
						new RealParameter("0.5")),
				"rhoSampling", new TimedParameter(
						new RealParameter("2.0"),
						new RealParameter("0.5")),
				"removalProb", new SkylineVectorParameter(
						null,
						new RealParameter("0.5")));

		BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
		density.initByName(
				"parameterization", parameterization,
				"frequencies", new RealParameter("1.0"),
				"conditionOnSurvival", false,
				"tree", new TreeParser(newick, false, false, true,0),
				"typeLabel", "type",
				"parallelize", false,
				"useAnalyticalSingleTypeSolution", false);

		double logLnumerical = density.calculateLogP();

		BirthDeathMigrationDistribution densityExact = new BirthDeathMigrationDistribution();
		densityExact.initByName(
				"parameterization", parameterization,
				"frequencies", new RealParameter("1.0"),
				"conditionOnSurvival", false,
				"tree", new TreeParser(newick, false, false, true,0),
				"typeLabel", "type",
				"parallelize", false,
				"useAnalyticalSingleTypeSolution", true);

		double logLanalytical = densityExact.calculateLogP();

		assertEquals(logLnumerical, logLanalytical, 1e-5);
	}

	/**
	 * Two-state test for removal-probability rate change
	 * Reference from BDMM itself
     */
	@Test
	public void testLikelihoodRemovalProbChangeTwoState() {

		String newick = "((1[&type=0]: 1.5, 2[&type=1]: 0)3[&type=0]: 3.5, 4[&type=1]: 4) ;";

        Parameterization parameterization = new EpiParameterization();
		parameterization.initByName(
                "typeSet", new TypeSet(2),
                "processLength", new RealParameter("6.0"),
                "R0", new SkylineVectorParameter(
                        new RealParameter("1.0"),
                        new RealParameter((4.0/3.0) + " 1.1"),
                        2),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        new RealParameter("1.0"),
                        new RealParameter("1.5 1.4"),
                        2),
                "R0AmongDemes", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.0"),
                        2),
                "migrationRate", new SkylineMatrixParameter(
                        new RealParameter("1.0"),
                        new RealParameter("0.2 0.3"),
                        2),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.33"),
                        2),
                "removalProb", new SkylineVectorParameter(
                        new RealParameter("1.0"),
                        new RealParameter("0.3 0.4 0.7 0.6")));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
		density.initByName(
		        "parameterization", parameterization,
                "frequencies", new RealParameter("0.5 0.5"),
                "conditionOnSurvival", false,
                "tree", new TreeParser(newick, false, false, true,0),
                "typeLabel", "type",
                "parallelize", true);

		double logL = density.calculateLogP();

        // Reference BDMM (version 0.2.0) 29/03/2018
		assertEquals(-21.185194919464568 + labeledTreeConversionFactor(density), logL, 1e-5);
	}


	/**
	 * Test simple configuration with one rho-sampling event.
	 * One type, no psi-sampling, no sampled-ancestor
	 * No rate-changes
	 * Reference: R
     */
	@Test
	public void testSingleRho() {

		Tree tree = new TreeParser("((1[&type=0]: 4.5, 2[&type=0]: 4.5):1,3[&type=0]:5.5);",false);

        Parameterization parameterization = new EpiParameterization();
		parameterization.initByName(
		        "typeSet", new TypeSet(1),
                "processLength", tree,
                "R0", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.5")),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.5")),
                "R0AmongDemes", new SkylineMatrixParameter(
                        null,
                        null),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        null),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.0")),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0")),
                "rhoSampling", new TimedParameter(
                        new RealParameter("5.5"),
                        new RealParameter("0.01")));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
		density.initByName(
		        "parameterization", parameterization,
                "frequencies", new RealParameter("1.0"),
                "conditionOnSurvival", false,
				"conditionOnRoot", true,
                "tree", tree,
                "typeLabel", "type",
                "parallelize", false,
				"useAnalyticalSingleTypeSolution", false);

		double logL = density.calculateLogP();

		// this result is from R: LikConstant(2.25,1.5,0.01,c(4.5,5.5),root=1,survival=0)
		assertEquals(-6.761909 + labeledTreeConversionFactor(density), logL, 1e-4);

		density.setInputValue("useAnalyticalSingleTypeSolution", true);
		density.initAndValidate();

		double logLanalytical = density.calculateLogP();
		assertEquals(-6.761909 + labeledTreeConversionFactor(density), logLanalytical, 1e-4);

		// test with conditioned-on-survival tree
		parameterization.setInputValue("processLength", "10");
        parameterization.setInputValue("rhoSampling",
                new TimedParameter(new RealParameter("10"),
                        new RealParameter("0.01")));
        parameterization.initAndValidate();

        density.setInputValue("conditionOnSurvival", true);
		density.setInputValue("conditionOnRoot", false);
		density.setInputValue("useAnalyticalSingleTypeSolution", false);
		density.initAndValidate();

		double logL2 = density.calculateLogP();

        // this result is from R: LikConstant(2.25,1.5,0.01,c(4.5,5.5,5.5+1e-100),root=0,survival=1)
		assertEquals(-7.404227 + labeledTreeConversionFactor(density), logL2, 1e-4);

		density.setInputValue("useAnalyticalSingleTypeSolution", false);
		density.initAndValidate();

		double logL2analytical = density.calculateLogP();

		assertEquals(-7.404227 + labeledTreeConversionFactor(density), logL2analytical, 1e-4);
	}

	/**
	 * Basic test for rho-sampling in the past
	 * One type, 2 tips, one state
	 * No rate changes
	 * @throws Exception
	 */
	@Test
	public void testMultiRho2tips() throws Exception {

		// two tips sampled at the same time
		Tree tree = new TreeParser("(3[&type=0]: 4, 4[&type=0]: 4) ;",false);

		RealParameter originParam = new RealParameter("5.0");

        Parameterization parameterization = new EpiParameterization();
		parameterization.initByName(
                "typeSet", new TypeSet(1),
                "processLength", originParam,
                "R0", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.5")),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.5")),
                "R0AmongDemes", new SkylineMatrixParameter(
                        null,
                        null),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        null),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.0")),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0")),
                "rhoSampling", new TimedParameter(
                        new RealParameter("0.0 2.5"),
                        new RealParameter("1.0 0.2"),
                        originParam));


        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
		density.initByName(
		        "parameterization", parameterization,
                "frequencies", new RealParameter("1.0"),
                "conditionOnSurvival", true,
                "tree", tree,
                "typeLabel", "type",
                "parallelize", false,
				"useAnalyticalSingleTypeSolution", false);

		double logL = density.calculateLogP();

		// this result is from BEAST: BDSKY, not double checked in R
		assertEquals(-10.569863754307026, logL, 1e-4);

		density.setInputValue("useAnalyticalSingleTypeSolution", true);
		density.initAndValidate();

		double logLanalytical = density.calculateLogP();
		assertEquals(-10.569863754307026, logLanalytical, 1e-4);

		// tips sampled at two different times
		tree = new TreeParser("(3[&type=0]: 1.5, 4[&type=0]: 4) ;",false);
		density.setInputValue("tree", tree);
		density.setInputValue("useAnalyticalSingleTypeSolution", false);
		density.initAndValidate();

		double logL2 = density.calculateLogP();

		// this result is from BEAST: BDSKY, not double checked in R
		assertEquals(-8.099631076932816, logL2, 1e-4);

		density.setInputValue("useAnalyticalSingleTypeSolution", false);
		density.initAndValidate();

		double logL2analytical = density.calculateLogP();

		assertEquals(-8.099631076932816, logL2, 1e-4);
	}

	/**
	 * Test with combined multiple-rho-sampling events in the past and psi-sampling
	 * was "TestRhoSasha"
	 * One type, no rate-changes, no sampled-ancestors
	 * 26 tips
	 * @throws Exception
	 */
	@Test
	public void testMultiRhoSampling() throws Exception {
		// Uncoloured tree
		Tree tree = new TreeParser("(((((t1[&type=0]:0.4595008531,t25[&type=0]:0.4595008531)[&type=0]:0.3373053072,t23[&type=0]:0.3567584538)[&type=0]:0.007310819036,t16[&type=0]:0.3489190732)[&type=0]:0.331009529,((t18[&type=0]:0.03315384045,t14[&type=0]:0.03315384045)[&type=0]:0.5063451374,(t10[&type=0]:0.4211543131,t15[&type=0]:0.4211543131)[&type=0]:0.1183446648)[&type=0]:0.5956275305)[&type=0]:0.1158090878,((t19[&type=0]:0.9429393194,((t6[&type=0]:0.363527235,t11[&type=0]:0.4417423167)[&type=0]:0.01881829549,((((t3[&type=0]:0.3071904376,(((t24[&type=0]:0.01065209364,t13[&type=0]:0.01065209364)[&type=0]:0.06076485145,t8[&type=0]:0.07141694509)[&type=0]:0.123620245,(t22[&type=0]:0.1616119808,t2[&type=0]:0.1616119808)[&type=0]:0.03342520927)[&type=0]:0.1121532475)[&type=0]:0.24520579,t9[&type=0]:0.5523962276)[&type=0]:0.3852615426,(((t20[&type=0]:0.2935970782,(t17[&type=0]:0.06569090089,t4[&type=0]:0.06569090089)[&type=0]:0.2279061773)[&type=0]:0.08350780408,(t21[&type=0]:0.05109047139,t5[&type=0]:0.05109047139)[&type=0]:0.3260144109)[&type=0]:0.2298344132,t7[&type=0]:0.6069392955)[&type=0]:0.3307184747)[&type=0]:0.01206284377,t26[&type=0]:0.9497206139)[&type=0]:0.05755333197)[&type=0]:0.03290891884)[&type=0]:0.07263755325,t12[&type=0]:1.112820418)[&type=0]:0.1381151782);",false);

        RealParameter originParam = new RealParameter("2.0");

        Parameterization parameterization = new EpiParameterization();
		parameterization.initByName(
                "typeSet", new TypeSet(1),
                "processLength", originParam,
                "R0", new SkylineVectorParameter(
                        null,
                        new RealParameter(new Double[]{3.0/4.5})),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("4.5")),
                "R0AmongDemes", new SkylineMatrixParameter(
                        null,
                        null),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        null),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealParameter(new Double[]{2.0/4.5})),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0")),
                "rhoSampling", new TimedParameter(
                        new RealParameter("1.0 1.5 2.0"),
                        new RealParameter("0.0 0.05 0.01")));


        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
		density.initByName(
		        "parameterization", parameterization,
                "frequencies", new RealParameter("1.0"),
                "conditionOnSurvival", false,
                "tree", tree,
                "typeLabel", "type",
                "parallelize", true,
				"useAnalyticalSingleTypeSolution", false);

		assertEquals(-124.96086690757612 + labeledTreeConversionFactor(density),
				density.calculateLogP(), 1e-2);     // this result is from BEAST, not double checked in R

		density.setInputValue("useAnalyticalSingleTypeSolution", true);
		density.initAndValidate();
		assertEquals(-124.96086690757612 + labeledTreeConversionFactor(density),
				density.calculateLogP(), 1e-2);     // this result is from BEAST, not double checked in R

        parameterization.setInputValue("rhoSampling",
                new TimedParameter(new RealParameter("0.0 0.5 1.0"),
                        new RealParameter("0.01 0.05 0.0"),
                        originParam));
        parameterization.initAndValidate();
        density.setInputValue("useAnalyticalSingleTypeSolution", false);
        density.initAndValidate();

        assertEquals(-124.96086690757612 + labeledTreeConversionFactor(density),
				density.calculateLogP(), 1e-2);     // this result is from BEAST, not double checked in R

		density.setInputValue("useAnalyticalSingleTypeSolution", true);
		density.initAndValidate();

		assertEquals(-124.96086690757612 + labeledTreeConversionFactor(density),
				density.calculateLogP(), 1e-2);     // this result is from BEAST, not double checked in R
	}

	/**
	 * Test with multiple cases for rho-sampling in the past combined with rate changes
	 * 1 state, no sampled ancestors
	 * 26 tips
     */
	@Test
	public void testMultiRhoWithRateChanges1() {

        Tree tree = new TreeParser("(((((t1[&type=0]:0.4595008531,t25[&type=0]:0.4595008531)[&type=0]:0.3373053072,t23[&type=0]:0.3567584538)[&type=0]:0.007310819036,t16[&type=0]:0.3489190732)[&type=0]:0.331009529,((t18[&type=0]:0.03315384045,t14[&type=0]:0.03315384045)[&type=0]:0.5063451374,(t10[&type=0]:0.4211543131,t15[&type=0]:0.4211543131)[&type=0]:0.1183446648)[&type=0]:0.5956275305)[&type=0]:0.1158090878,((t19[&type=0]:0.9429393194,((t6[&type=0]:0.363527235,t11[&type=0]:0.4417423167)[&type=0]:0.01881829549,((((t3[&type=0]:0.3071904376,(((t24[&type=0]:0.01065209364,t13[&type=0]:0.01065209364)[&type=0]:0.06076485145,t8[&type=0]:0.07141694509)[&type=0]:0.123620245,(t22[&type=0]:0.1616119808,t2[&type=0]:0.1616119808)[&type=0]:0.03342520927)[&type=0]:0.1121532475)[&type=0]:0.24520579,t9[&type=0]:0.5523962276)[&type=0]:0.3852615426,(((t20[&type=0]:0.2935970782,(t17[&type=0]:0.06569090089,t4[&type=0]:0.06569090089)[&type=0]:0.2279061773)[&type=0]:0.08350780408,(t21[&type=0]:0.05109047139,t5[&type=0]:0.05109047139)[&type=0]:0.3260144109)[&type=0]:0.2298344132,t7[&type=0]:0.6069392955)[&type=0]:0.3307184747)[&type=0]:0.01206284377,t26[&type=0]:0.9497206139)[&type=0]:0.05755333197)[&type=0]:0.03290891884)[&type=0]:0.07263755325,t12[&type=0]:1.112820418)[&type=0]:0.1381151782);", false);

        // no rate-change, rho-sampling at present
        RealParameter originParam = new RealParameter("2.0");
        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(1),
                "processLength", originParam,
                "birthRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("2.0")),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.5")),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.5")),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0")),
                "birthRateAmongDemes", new SkylineMatrixParameter(
                        null,
                        null),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        null),
                "rhoSampling", new TimedParameter(
                        new RealParameter("0.0"),
                        new RealParameter("1.0"),
                        originParam));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("parameterization", parameterization,
                "frequencies", new RealParameter("1.0"),
                "conditionOnSurvival", true,
                "tree", tree,
                "typeLabel", "type",
                "parallelize", false,
				"useAnalyticalSingleTypeSolution", false);

        //		System.out.println("\na) Likelihood: " + bdssm.calculateTreeLogLikelihood(tree));
        assertEquals(-21.42666177086957 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-5);

		BirthDeathMigrationDistribution densityExact = new BirthDeathMigrationDistribution();
		densityExact.initByName("parameterization", parameterization,
				"frequencies", new RealParameter("1.0"),
				"conditionOnSurvival", true,
				"tree", tree,
				"typeLabel", "type",
				"parallelize", false,
				"useAnalyticalSingleTypeSolution", true);

		assertEquals(-21.42666177086957 + labeledTreeConversionFactor(density), densityExact.calculateLogP(), 1e-5);
	}


    @Test
    public void testMultiRhoWithRateChanges2() {

        Tree tree = new TreeParser("(((((t1[&type=0]:0.4595008531,t25[&type=0]:0.4595008531)[&type=0]:0.3373053072,t23[&type=0]:0.3567584538)[&type=0]:0.007310819036,t16[&type=0]:0.3489190732)[&type=0]:0.331009529,((t18[&type=0]:0.03315384045,t14[&type=0]:0.03315384045)[&type=0]:0.5063451374,(t10[&type=0]:0.4211543131,t15[&type=0]:0.4211543131)[&type=0]:0.1183446648)[&type=0]:0.5956275305)[&type=0]:0.1158090878,((t19[&type=0]:0.9429393194,((t6[&type=0]:0.363527235,t11[&type=0]:0.4417423167)[&type=0]:0.01881829549,((((t3[&type=0]:0.3071904376,(((t24[&type=0]:0.01065209364,t13[&type=0]:0.01065209364)[&type=0]:0.06076485145,t8[&type=0]:0.07141694509)[&type=0]:0.123620245,(t22[&type=0]:0.1616119808,t2[&type=0]:0.1616119808)[&type=0]:0.03342520927)[&type=0]:0.1121532475)[&type=0]:0.24520579,t9[&type=0]:0.5523962276)[&type=0]:0.3852615426,(((t20[&type=0]:0.2935970782,(t17[&type=0]:0.06569090089,t4[&type=0]:0.06569090089)[&type=0]:0.2279061773)[&type=0]:0.08350780408,(t21[&type=0]:0.05109047139,t5[&type=0]:0.05109047139)[&type=0]:0.3260144109)[&type=0]:0.2298344132,t7[&type=0]:0.6069392955)[&type=0]:0.3307184747)[&type=0]:0.01206284377,t26[&type=0]:0.9497206139)[&type=0]:0.05755333197)[&type=0]:0.03290891884)[&type=0]:0.07263755325,t12[&type=0]:1.112820418)[&type=0]:0.1381151782);", false);

        // rate-changes, rho-sampling in the past
        RealParameter originParam = new RealParameter("2.0");
        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(1),
                "processLength", originParam,
                "R0", new SkylineVectorParameter(
                        new RealParameter("1.0 1.5"),
                        new RealParameter(new Double[]{3.0/4.5, 2.0/1.5, 4.0/1.5})),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        new RealParameter("1.0 1.5"),
                        new RealParameter("4.5 1.5 1.5")),
                "samplingProportion", new SkylineVectorParameter(
                        new RealParameter("1.0 1.5"),
                        new RealParameter(new Double[]{2.0/4.5, 0.5/1.5, 1.0/1.5})),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0")),
                "R0AmongDemes", new SkylineMatrixParameter(
                        null,
                        null),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        null),
                "rhoSampling", new TimedParameter(
                        new RealParameter("2.0"),
                        new RealParameter("0.01")));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("parameterization", parameterization,
                "frequencies", new RealParameter("1.0"),
                "conditionOnSurvival", false,
                "tree", tree,
                "typeLabel", "type",
                "parallelize", false,
				"useAnalyticalSingleTypeSolution", false);

        assertEquals(-87.59718586549747 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-5);

        density.setInputValue("useAnalyticalSingleTypeSolution", true);
        density.initAndValidate();
		assertEquals(-87.59718586549747 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-5);
	}

    @Test
    public void testMultiRhoWithRateChanges3() {

        Tree tree = new TreeParser("(((((t1[&type=0]:0.4595008531,t25[&type=0]:0.4595008531)[&type=0]:0.3373053072,t23[&type=0]:0.3567584538)[&type=0]:0.007310819036,t16[&type=0]:0.3489190732)[&type=0]:0.331009529,((t18[&type=0]:0.03315384045,t14[&type=0]:0.03315384045)[&type=0]:0.5063451374,(t10[&type=0]:0.4211543131,t15[&type=0]:0.4211543131)[&type=0]:0.1183446648)[&type=0]:0.5956275305)[&type=0]:0.1158090878,((t19[&type=0]:0.9429393194,((t6[&type=0]:0.363527235,t11[&type=0]:0.4417423167)[&type=0]:0.01881829549,((((t3[&type=0]:0.3071904376,(((t24[&type=0]:0.01065209364,t13[&type=0]:0.01065209364)[&type=0]:0.06076485145,t8[&type=0]:0.07141694509)[&type=0]:0.123620245,(t22[&type=0]:0.1616119808,t2[&type=0]:0.1616119808)[&type=0]:0.03342520927)[&type=0]:0.1121532475)[&type=0]:0.24520579,t9[&type=0]:0.5523962276)[&type=0]:0.3852615426,(((t20[&type=0]:0.2935970782,(t17[&type=0]:0.06569090089,t4[&type=0]:0.06569090089)[&type=0]:0.2279061773)[&type=0]:0.08350780408,(t21[&type=0]:0.05109047139,t5[&type=0]:0.05109047139)[&type=0]:0.3260144109)[&type=0]:0.2298344132,t7[&type=0]:0.6069392955)[&type=0]:0.3307184747)[&type=0]:0.01206284377,t26[&type=0]:0.9497206139)[&type=0]:0.05755333197)[&type=0]:0.03290891884)[&type=0]:0.07263755325,t12[&type=0]:1.112820418)[&type=0]:0.1381151782);", false);

        // rate-changes, rho-sampling in the past and present
        RealParameter originParam = new RealParameter("2.0");
        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(1),
                "processLength", originParam,
                "R0", new SkylineVectorParameter(
                        new RealParameter("1.0 1.5"),
                        new RealParameter(new Double[]{3.0/4.5, 2.0/1.5, 4.0/1.5})),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        new RealParameter("1.0 1.5"),
                        new RealParameter("4.5 1.5 1.5")),
                "samplingProportion", new SkylineVectorParameter(
                        new RealParameter("1.0 1.5"),
                        new RealParameter(new Double[]{2.0/4.5, 0.5/1.5, 1.0/1.5})),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0")),
                "R0AmongDemes", new SkylineMatrixParameter(
                        null,
                        null),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        null),
                "rhoSampling", new TimedParameter(
                        new RealParameter("1.0 2.0"),
                        new RealParameter("0.05 0.01")));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("parameterization", parameterization,
                "frequencies", new RealParameter("1.0"),
                "conditionOnSurvival", false,
                "tree", tree,
                "typeLabel", "type",
                "parallelize", false,
				"useAnalyticalSingleTypeSolution", false);

        assertEquals(-87.96488 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-1);

        density.setInputValue("useAnalyticalSingleTypeSolution", true);
		assertEquals(-87.96488 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-1);
	}

    @Test
    public void testMultiRhoWithRateChanges4() {

        Tree tree = new TreeParser("(((((t1[&type=0]:0.4595008531,t25[&type=0]:0.4595008531)[&type=0]:0.3373053072,t23[&type=0]:0.3567584538)[&type=0]:0.007310819036,t16[&type=0]:0.3489190732)[&type=0]:0.331009529,((t18[&type=0]:0.03315384045,t14[&type=0]:0.03315384045)[&type=0]:0.5063451374,(t10[&type=0]:0.4211543131,t15[&type=0]:0.4211543131)[&type=0]:0.1183446648)[&type=0]:0.5956275305)[&type=0]:0.1158090878,((t19[&type=0]:0.9429393194,((t6[&type=0]:0.363527235,t11[&type=0]:0.4417423167)[&type=0]:0.01881829549,((((t3[&type=0]:0.3071904376,(((t24[&type=0]:0.01065209364,t13[&type=0]:0.01065209364)[&type=0]:0.06076485145,t8[&type=0]:0.07141694509)[&type=0]:0.123620245,(t22[&type=0]:0.1616119808,t2[&type=0]:0.1616119808)[&type=0]:0.03342520927)[&type=0]:0.1121532475)[&type=0]:0.24520579,t9[&type=0]:0.5523962276)[&type=0]:0.3852615426,(((t20[&type=0]:0.2935970782,(t17[&type=0]:0.06569090089,t4[&type=0]:0.06569090089)[&type=0]:0.2279061773)[&type=0]:0.08350780408,(t21[&type=0]:0.05109047139,t5[&type=0]:0.05109047139)[&type=0]:0.3260144109)[&type=0]:0.2298344132,t7[&type=0]:0.6069392955)[&type=0]:0.3307184747)[&type=0]:0.01206284377,t26[&type=0]:0.9497206139)[&type=0]:0.05755333197)[&type=0]:0.03290891884)[&type=0]:0.07263755325,t12[&type=0]:1.112820418)[&type=0]:0.1381151782);", false);

        // rate-changes, rho-sampling in the past and present, with reversed times
        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "processLength", tree,
                "typeSet", new TypeSet(1),
                "R0", new SkylineVectorParameter(
                        new RealParameter("0.5 1.0 1.1"),
                        new RealParameter(new Double[]{3.0/4.5, 2.0/1.5, 4.0/1.5, 4.0/2.5})),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        new RealParameter("0.5 1.0 1.1"),
                        new RealParameter("4.5 1.5 1.5 2.5")),
                "samplingProportion", new SkylineVectorParameter(
                        new RealParameter("0.5 1.0 1.1"),
                        new RealParameter(new Double[]{2.0/4.5, 0.5/1.5, 1.0/1.5, 2.0/2.5})),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0")),
                "R0AmongDemes", new SkylineMatrixParameter(
                        null,
                        null),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        null),
                "rhoSampling", new TimedParameter(
                        new RealParameter("1.0 " + tree.getRoot().getHeight()),
                        new RealParameter("0.05 0.01")));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("parameterization", parameterization,
                "frequencies", new RealParameter("1.0"),
                "conditionOnSurvival", false,
				"conditionOnRoot", true,
                "tree", tree,
                "typeLabel", "type",
                "parallelize", false,
				"useAnalyticalSingleTypeSolution", false);

		assertEquals(-100.15682190617582 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-1);

		density.setInputValue("useAnalyticalSingleTypeSolution", true);
		assertEquals(-100.15682190617582 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-1);
	}

	/**
	 * Basic 1-dim test
	 * No rate change, 1 state, no rho-sampling
	 * Reference from BDSKY
	 */
	@Test
	public void testLikelihood1dim() {

        Tree tree = new TreeParser( "((3[&state=0] : 1.5, 4[&state=0] : 0.5)[&state=0] : 1 , (1[&state=0] : 2, 2[&state=0] : 1)[&state=0] : 3)[&state=0];",
                false);

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "processLength", new RealParameter("6.0"),
                "typeSet", new TypeSet(1),
                "R0", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.3333333334")),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.5")),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.33333333333")),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0")));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("parameterization", parameterization,
                "frequencies", new RealParameter("1.0"),
                "conditionOnSurvival", false,
                "tree", tree,
                "typeLabel", "state",
                "parallelize", false,
				"useAnalyticalSingleTypeSolution", false);

		assertEquals(-19.019796073623493 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-5);   // Reference BDSKY (version 1.3.3)

		density.setInputValue("useAnalyticalSingleTypeSolution", true);
		density.initAndValidate();
		assertEquals(-19.019796073623493 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-5);   // Reference BDSKY (version 1.3.3)

        density.setInputValue("conditionOnSurvival", true);
		density.setInputValue("useAnalyticalSingleTypeSolution", false);
        density.initAndValidate();

		assertEquals(-18.574104140202046 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-5); // Reference BDSKY (version 1.3.3)

		density.setInputValue("useAnalyticalSingleTypeSolution", true);
		density.initAndValidate();

		assertEquals(-18.574104140202046 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-5); // Reference BDSKY (version 1.3.3)
	}

	/**
	 * 1-dim and 1 rate-change test
	 * reference from BDSKY
     */
	@Test
	public void testLikelihoodRateChange1dim() {

        Tree tree = new TreeParser("((3[&state=0] : 1.5, 4[&state=0] : 0.5)[&state=0] : 1 , (1[&state=0] : 2, 2[&state=0] : 1)[&state=0] : 3)[&state=0];",
                false);

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "processLength", new RealParameter("6.0"),
                "typeSet", new TypeSet(1),
                "R0", new SkylineVectorParameter(
                        new RealParameter("3.0"),
                        new RealParameter("0.6666666667 1.3333333334")),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        new RealParameter("3.0"),
                        new RealParameter("4.5 1.5")),
                "samplingProportion", new SkylineVectorParameter(
                        new RealParameter("3.0"),
                        new RealParameter("0.4444444444 0.33333333333")),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0")));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("parameterization", parameterization,
                "frequencies", new RealParameter("1.0"),
                "conditionOnSurvival", false,
                "tree", tree,
                "typeLabel", "state",
                "parallelize", false,
				"useAnalyticalSingleTypeSolution", false);

        assertEquals(-33.7573 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-4); // Reference BDSKY

		density.setInputValue("useAnalyticalSingleTypeSolution", true);
		density.initAndValidate();
		assertEquals(-33.7573 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-4); // Reference BDSKY
	}

	/**
	 * Compare analytical and numerical for one rate change implementing
	 * truncated sampling (no sampling before specific time).
	 */
	@Test
	public void testLikelihoodRateTrancatedSampling1dim() {

		Tree tree = new TreeParser("((3[&state=0] : 1.5, 4[&state=0] : 0.5)[&state=0] : 1 , (1[&state=0] : 2, 2[&state=0] : 1)[&state=0] : 3)[&state=0];",
				false);

		Parameterization parameterization = new EpiParameterization();
		parameterization.initByName(
				"processLength", new RealParameter("6.0"),
				"typeSet", new TypeSet(1),
				"R0", new SkylineVectorParameter(
						null,
						new RealParameter("1.2")),
				"becomeUninfectiousRate", new SkylineVectorParameter(
						null,
						new RealParameter("1.0")),
				"samplingProportion", new SkylineVectorParameter(
						new RealParameter("2.4"),
						new RealParameter("0 0.33333333333")),
				"removalProb", new SkylineVectorParameter(
						null,
						new RealParameter("1.0")));

		BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
		density.initByName("parameterization", parameterization,
				"frequencies", new RealParameter("1.0"),
				"conditionOnSurvival", false,
				"tree", tree,
				"typeLabel", "state",
				"parallelize", false,
				"useAnalyticalSingleTypeSolution", false);

		double logPNumerical = density.calculateLogP();

		density.setInputValue("useAnalyticalSingleTypeSolution", true);
		density.initAndValidate();
		double logPAnalytical = density.calculateLogP();

		assertEquals(logPNumerical, logPAnalytical, 1e-4);
	}

	/**
	 * 1-dim and 1 rate-change test
	 * reference from BDSKY
	 */
	@Test
	public void testLikelihoodRateChangeCondOnSampling1dim() {

		Tree tree = new TreeParser("((3[&state=0] : 1.5, 4[&state=0] : 0.5)[&state=0] : 1 , (1[&state=0] : 2, 2[&state=0] : 1)[&state=0] : 3)[&state=0];",
				false);

		Parameterization parameterization = new EpiParameterization();
		parameterization.initByName(
				"processLength", new RealParameter("6.0"),
				"typeSet", new TypeSet(1),
				"R0", new SkylineVectorParameter(
						new RealParameter("3.0"),
						new RealParameter("0.6666666667 1.3333333334")),
				"becomeUninfectiousRate", new SkylineVectorParameter(
						new RealParameter("3.0"),
						new RealParameter("4.5 1.5")),
				"samplingProportion", new SkylineVectorParameter(
						new RealParameter("2.4"),
						new RealParameter("0.0 0.33333333333")),
				"removalProb", new SkylineVectorParameter(
						null,
						new RealParameter("1.0")));

		BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
		density.initByName("parameterization", parameterization,
				"frequencies", new RealParameter("1.0"),
				"conditionOnSurvival", true,
				"tree", tree,
				"typeLabel", "state",
				"parallelize", false,
				"useAnalyticalSingleTypeSolution", false);

		double logPnumeric = density.calculateLogP();
		System.out.println("Numerical solution: " + logPnumeric);

//		assertEquals(-33.7573 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-4); // Reference BDSKY

		density.setInputValue("useAnalyticalSingleTypeSolution", true);
		density.initAndValidate();

		double logPanalytic = density.calculateLogP();
		System.out.println("Analytical solution: " + logPnumeric);

		assertEquals(logPnumeric, logPanalytic, 1e-4);

//		assertEquals(-33.7573 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-4); // Reference BDSKY
	}

	/**
	 * Basic tests on 2 types situations with migration or birth among demes
	 * reference from R
	 * @throws Exception
	 */
	@Test 
	public void testLikelihoodCalculationMigTiny() throws Exception {

		// migration and no birth among demes

		Tree tree = new TreeParser("(1[&state=0] : 1.5, 2[&state=1] : 0.5)[&state=0];", false);
        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "processLength", new RealParameter("2.5"),
                "typeSet", new TypeSet(2),
                "R0", new SkylineVectorParameter(
                        null,
                        new RealParameter(Double.toString(4.0/3.0)), 2),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.5"), 2),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealParameter(Double.toString(1.0/3.0)), 2),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.1"), 2),

                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0"), 2));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("parameterization", parameterization,
                "frequencies", new RealParameter("0.5 0.5"),
                "conditionOnSurvival", false,
                "tree", tree,
                "typeLabel", "state",
                "parallelize", false);

		assertEquals(-7.215222 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-6); // result from R

		// no migration, symmetric birth among demes

        parameterization.setInputValue("migrationRate", null);
        parameterization.setInputValue("R0AmongDemes", new SkylineMatrixParameter(
                null,
                new RealParameter("0.0666667"), 2));
        parameterization.initAndValidate();
        density.initAndValidate();

		assertEquals(-7.404888 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-6); // result from R

		// no migration, asymmetric birth among demes

        parameterization.setInputValue("R0AmongDemes", new SkylineMatrixParameter(
                null,
                new RealParameter("0.0666667 0.1"), 2));
        parameterization.initAndValidate();
        density.initAndValidate();

		assertEquals(-7.18723 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-6); // result from R


        // no migration, asymmetric R0, asymmetric birth among demes

        parameterization.setInputValue("R0", new SkylineVectorParameter(
                null,
                new RealParameter("2 1.3333333")));
        parameterization.initAndValidate();
        density.initAndValidate();

		assertEquals(-7.350649 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-6); // result from R

        // no migration, asymmetric R0, birth among demes, BU rate, samp proportion

        parameterization.setInputValue("R0", new SkylineVectorParameter(
                null,
                new RealParameter("2.0 1.5")));
        parameterization.setInputValue("becomeUninfectiousRate", new SkylineVectorParameter(
                null,
                new RealParameter("2.0 1.0")));
        parameterization.setInputValue("samplingProportion", new SkylineVectorParameter(
                null,
                new RealParameter("0.5 0.3")));
        parameterization.setInputValue("R0AmongDemes", new SkylineMatrixParameter(
                null,
                new RealParameter("0.1 0.5")));
        parameterization.initAndValidate();
        density.initAndValidate();

		assertEquals(-6.504139 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-6); // result from R

		// Same params as last test, swapped leaf states

        tree = new TreeParser("(1[&state=1] : 1.5, 2[&state=0] : 0.5)[&state=0];", false);
        density.setInputValue("tree", tree);
        density.initAndValidate();

		assertEquals(-7.700916 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-6); // result from R
	}

	/**
	 * Test migration
	 * 2 types, migration, no birth among demes
	 * Adapted from BDSKY
	 */
	@Test
	public void testLikelihoodCalculationMig() {

        // uncoloured tree, asymmetric types
        Tree tree = new TreeParser(
                "((3[&type=0] : 1.5, 4[&type=1] : 0.5) : 1 , (1[&type=1] : 2, 2[&type=0] : 1) : 3);",
                false);

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "processLength", new RealParameter("6.0"),
                "typeSet", new TypeSet(2),
                "R0", new SkylineVectorParameter(
                        null,
                        new RealParameter((4.0 / 3.0) + " " + 5.0)),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.5 1.25")),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealParameter((1.0 / 3.0) + " " + (1.0/2.0))),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.2 0.1")),

                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0"), 2));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("parameterization", parameterization,
                "frequencies", new RealParameter("0.5 0.5"),
                "conditionOnSurvival", false,
                "tree", tree,
                "typeLabel", "type",
                "parallelize", false);

        assertEquals(-26.53293 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-5);
    }

	/**
	 * Test migration on big tree
	 * 2 types, migration, no birth among demes
	 * Adapted from BDSKY
	 */
    @Test
    public void testLikelihoodCalculationMigBig () {
		// uncoloured tree, 291 tips

        Tree tree = new TreeParser(
                "(((((((t1[&type=0]:0.9803361397,t2[&type=0]:0.9035540882):0.0532383481,t3[&type=0]:0.2637392259):0.6273536528,(t4[&type=0]:0.8624112266,t5[&type=0]:0.3278892266):0.2606245542):0.2941323873,(t6[&type=0]:0.09820114588,t7[&type=0]:0.533115675):0.8625875909):0.7040311908,(((t8[&type=0]:0.8696136218,t9[&type=0]:0.08719484485):0.4204288905,(t10[&type=0]:0.102143287,(t11[&type=0]:0.9850614571,t12[&type=0]:0.7407912319):0.8715072596):0.5182644848):0.524062254,(((((((t13[&type=0]:0.3981794417,(t14[&type=0]:0.03889928572,t15[&type=0]:0.5187105467):0.1127638209):0.3431177251,((t16[&type=0]:0.4239511855,t17[&type=0]:0.001895790454):0.690600364,t18[&type=0]:0.6283850113):0.4073564562):0.6862231812,(((t19[&type=0]:0.9947085041,t20[&type=0]:0.4739363373):0.1873670686,t21[&type=0]:0.151270482):0.803061039,((t22[&type=0]:0.8899249982,((t23[&type=0]:0.1329096023,t24[&type=0]:0.84205155):0.8838408566,(t25[&type=0]:0.7541888549,t26[&type=0]:0.8602364615):0.8912267659):0.771449636):0.1022819551,(((t27[&type=0]:0.3134289116,(t28[&type=0]:0.2446750235,t29[&type=0]:0.8565168788):0.8277210968):0.4307989818,((t30[&type=0]:0.2330717787,t31[&type=0]:0.4438336496):0.6521712865,(t32[&type=0]:0.2534400895,t33[&type=0]:0.7885409284):0.3051449039):0.1196702593):0.4061951274,t34[&type=0]:0.8415271267):0.4365981282):0.753448925):0.1580670979):0.04210642632,(((t35[&type=0]:0.7504386581,t36[&type=0]:0.6328390085):0.9047614154,t37[&type=0]:0.4946133171):0.2264722914,((((t38[&type=0]:0.06683212146,t39[&type=0]:0.479845396):0.9424520086,t40[&type=0]:0.894530142):0.3844042511,(((t41[&type=0]:0.5215392481,t42[&type=0]:0.2366602973):0.8142298241,(t43[&type=0]:0.2968777204,(t44[&type=0]:0.655541793,t45[&type=0]:0.8608812049):0.3564132168):0.04912991729):0.1511388237,t46[&type=0]:0.9031036345):0.1874918914):0.9690212663,(t47[&type=0]:0.07753491728,(t48[&type=0]:0.8349514075,(t49[&type=0]:0.9689748741,t50[&type=0]:0.925813166):0.4534903264):0.3571097804):0.1324767114):0.5515443345):0.3330309158):0.7202291801,((t51[&type=0]:0.6977306763,((t52[&type=0]:0.9157640305,t53[&type=0]:0.4226291834):0.5872618856,t54[&type=0]:0.2063144948):0.1422286083):0.7182746637,t55[&type=0]:0.759545143):0.7437628019):0.2425582204,((t56[&type=0]:0.4614429038,(t57[&type=0]:0.9092229386,((t58[&type=0]:0.1049408391,t59[&type=0]:0.6328130178):0.642241966,((t60[&type=0]:0.264340204,t61[&type=0]:0.5904771155):0.7333205172,(t62[&type=0]:0.9183179205,t63[&type=0]:0.1090340314):0.3010568973):0.3240860389):0.3192155454):0.1835780439):0.5942421539,t64[&type=0]:0.7931551472):0.967891278):0.06263663713,(t65[&type=0]:0.5774453548,((t66[&type=0]:0.07208712469,((t67[&type=0]:0.8918803469,t68[&type=0]:0.5110983853):0.1491188321,t69[&type=0]:0.2471361952):0.9591872343):0.3133718621,(t70[&type=0]:0.944087367,t71[&type=0]:0.7830825299):0.2284035049):0.5492361034):0.1136150162):0.002181729767):0.4548798562):0.4258609388,((((((t72[&type=0]:0.27679418,t73[&type=0]:0.5398862793):0.8871422287,(((((t74[&type=0]:0.2531923286,t75[&type=0]:0.3796772889):0.4489221217,t76[&type=0]:0.2554209188):0.3248268673,t77[&type=0]:0.5372577759):0.5699883625,t78[&type=0]:0.1656995732):0.957750936,(t79[&type=0]:0.1301121258,t80[&type=0]:0.8925942327):0.2838441601):0.5258686764):0.47825964,(t81[&type=0]:0.5749240227,((t82[&type=0]:0.9574132746,(t83[&type=0]:0.00485483068,t84[&type=0]:0.8091488208):0.1985368489):0.3703975577,(((t85[&type=0]:0.3991035291,(t86[&type=0]:0.03201846033,t87[&type=0]:0.8380640063):0.05616304209):0.8414494572,t88[&type=0]:0.6844437125):0.2426782607,((t89[&type=0]:0.7543559887,t90[&type=0]:0.7162597755):0.8230077426,t91[&type=0]:0.08967904118):0.4460245941):0.8679371702):0.51572948):0.4362259945):0.2631344711,(((t92[&type=0]:0.3353162925,((t93[&type=0]:0.4025212794,t94[&type=0]:0.0281926766):0.7965471447,t95[&type=0]:0.1145715592):0.5993301494):0.08854756854,(t96[&type=0]:0.1461353719,((t97[&type=0]:0.3158547124,t98[&type=0]:0.06653800653):0.5634025722,t99[&type=0]:0.9711292514):0.9727503664):0.7684133062):0.4824229684,((t100[&type=0]:0.06834940333,t101[&type=0]:0.7794982188):0.3453287922,(t102[&type=0]:0.627945075,t103[&type=0]:0.1914187325):0.9974814849):0.6312927424):0.04858242651):0.2845227425,((t104[&type=0]:0.6782600286,(t105[&type=0]:0.03190574702,t106[&type=0]:0.5840284519):0.03041352634):0.725893975,(((t107[&type=0]:0.9885271091,t108[&type=0]:0.07126446022):0.8419693699,t109[&type=0]:0.1546431775):0.898004594,t110[&type=0]:0.2500803664):0.1493327522):0.4266726137):0.5946582041,(t111[&type=0]:0.1395377244,(((t112[&type=0]:0.7170655408,(t113[&type=0]:0.976886861,t114[&type=0]:0.9406369971):0.7471234254):0.8065501407,((t115[&type=0]:0.1713845057,(t116[&type=0]:0.7861330248,t117[&type=0]:0.6082276558):0.8413775554):0.3245444677,t118[&type=0]:0.3892389825):0.5992471091):0.7592411407,(((t119[&type=0]:0.535931844,t120[&type=0]:0.09058958571):0.4227561057,(t121[&type=0]:0.5531579193,t122[&type=0]:0.8276180199):0.6653355309):0.0941624688,t123[&type=0]:0.3623022255):0.1494971744):0.3526274569):0.9720881658):0.8149677955):0.6065687414,((((((t124[&type=0]:0.5406888947,t125[&type=0]:0.8892341822):0.06211395678,((t126[&type=0]:0.8203180477,(t127[&type=0]:0.8536844573,t128[&type=0]:0.360511546):0.9030223228):0.9095590916,((t129[&type=0]:0.9110714826,(t130[&type=0]:0.2346256471,t131[&type=0]:0.6523390864):0.1288849309):0.7077432328,(t132[&type=0]:0.4060195235,t133[&type=0]:0.1661393729):0.3910941551):0.205704404):0.8609933471):0.3724007562,((t134[&type=0]:0.1731842053,(t135[&type=0]:0.7232482471,(t136[&type=0]:0.3883952193,((t137[&type=0]:0.6709475764,t138[&type=0]:0.0372075201):0.5473196667,(t139[&type=0]:0.8092764446,t140[&type=0]:0.4123262055):0.2000603897):0.55258787):0.2654263263):0.745555162):0.2956101163,((t141[&type=0]:0.52147611,(t142[&type=0]:0.9462005703,t143[&type=0]:0.5671354234):0.6887917654):0.362258781,t144[&type=0]:0.4798202242):0.8242726682):0.6072624433):0.695287361,((((t145[&type=0]:0.03793937969,t146[&type=0]:0.07275558705):0.3482963489,t147[&type=0]:0.1457363514):0.1479936559,(t148[&type=0]:0.7158309214,((t149[&type=0]:0.2174433649,t150[&type=0]:0.04072828358):0.4112026501,t151[&type=0]:0.6422409331):0.3413406226):0.1693999742):0.6631712937,(((t152[&type=0]:0.2706006162,t153[&type=0]:0.9267972289):0.1387761638,((((t154[&type=0]:0.2563392594,t155[&type=0]:0.3058371837):0.5946117372,t156[&type=0]:0.6161190302):0.6970871226,(t157[&type=0]:0.2388902532,(t158[&type=0]:0.9486316761,t159[&type=0]:0.215360787):0.168830334):0.03888285463):0.1640696453,t160[&type=0]:0.6803096831):0.1418975852):0.4218000816,(((t161[&type=0]:0.8702562298,t162[&type=0]:0.9289729816):0.05807372741,t163[&type=0]:0.3533785399):0.5012762842,(((t164[&type=0]:0.8666574673,t165[&type=0]:0.9603798252):0.7887994377,t166[&type=0]:0.857058729):0.4139410679,(t167[&type=0]:0.5900272813,t168[&type=0]:0.3345388798):0.06017537019):0.9609203783):0.7103463742):0.696603697):0.6451920038):0.1909481271,((((t169[&type=0]:0.9171597108,t170[&type=0]:0.9479122513):0.7170342554,(t171[&type=0]:0.2722596873,((t172[&type=0]:0.1194724559,(t173[&type=0]:0.03922236571,t174[&type=0]:0.6290624789):0.07739861775):0.8598598302,(t175[&type=0]:0.2009421999,(t176[&type=0]:0.06154947914,t177[&type=0]:8.997193072E-4):0.04738179315):0.3235510678):0.3443877005):0.6351028818):0.5525081949,((((t178[&type=0]:0.7599076207,t179[&type=0]:0.2997759853):0.5921433992,t180[&type=0]:0.7098581635):0.3725496214,(t181[&type=0]:0.5053773888,(t182[&type=0]:0.5991492711,(t183[&type=0]:0.5036820578,t184[&type=0]:0.6361607853):0.510631816):0.9604382808):0.2464167587):0.6073093358,(((t185[&type=0]:0.03128415369,(t186[&type=0]:0.5260852403,(t187[&type=0]:0.878767435,t188[&type=0]:0.4992109234):0.5333148066):0.00347468094):0.5590308013,t189[&type=0]:0.3710992143):0.5034162949,(t190[&type=0]:0.778916508,((t191[&type=0]:0.3069154553,(((t192[&type=0]:0.9946115273,t193[&type=0]:0.9138687006):0.5209144899,t194[&type=0]:0.5152770842):0.9462409306,t195[&type=0]:0.7395236609):0.4110851623):0.930918345,(((t196[&type=0]:0.7895439987,((t197[&type=0]:0.4697002599,t198[&type=0]:0.1383787312):0.6911794308,(t199[&type=0]:0.8664436699,t200[&type=0]:0.1959039853):0.8656513852):0.3620497067):0.2839249384,(t201[&type=0]:0.6558795469,t202[&type=0]:0.2103423763):0.969477433):0.9058840063,(t203[&type=0]:0.0856692954,t204[&type=0]:0.4175976661):0.820434629):0.5355881769):0.2263581599):0.4512835185):0.7323478526):0.2479199937):0.1964542414,((t205[&type=0]:0.7537573762,(t206[&type=0]:0.1392466244,(t207[&type=0]:0.5136175761,(t208[&type=0]:0.7852529553,t209[&type=0]:0.07355738804):0.1220811389):0.7572090242):0.1422528555):0.5948274662,(((((t210[&type=0]:0.3068353184,(t211[&type=0]:0.3314456891,((t212[&type=0]:0.5265486804,t213[&type=0]:0.1382007354):0.1814086549,t214[&type=0]:0.9276472756):0.07718444197):0.03486835537):0.1617580003,(t215[&type=0]:0.3328830956,t216[&type=0]:0.8558843595):0.8366736979):0.347376487,t217[&type=0]:0.8222538356):0.2337225529,(t218[&type=0]:0.06199815008,t219[&type=0]:0.45975962):0.179990889):0.0635867205,(t220[&type=0]:0.3214025751,(t221[&type=0]:0.5022090652,t222[&type=0]:0.6454557138):0.6956466341):0.2711792416):0.1847200533):0.1051658324):0.4945860899):0.936143348,(((t223[&type=0]:0.06268779701,((t224[&type=0]:0.3337278806,t225[&type=0]:0.1570303424):0.3089733059,(t226[&type=0]:0.5069784883,t227[&type=0]:0.1434204187):0.2001587199):0.04750720505):0.3600859912,((((t228[&type=0]:0.9994731578,(t229[&type=0]:0.8934116936,t230[&type=0]:0.03698333143):0.8173468311):0.3089058488,((((t231[&type=0]:0.3216121283,t232[&type=0]:0.5232846253):0.8687884973,(t233[&type=0]:0.6280638413,((t234[&type=0]:0.6543256822,t235[&type=0]:0.8677638234):0.8895299246,t236[&type=0]:0.4047793006):0.7147388768):0.3533478715):0.9470084386,t237[&type=0]:0.7769409856):0.4955915695,((t238[&type=0]:0.2772087415,(t239[&type=0]:0.4904922615,(t240[&type=0]:0.05356206303,t241[&type=0]:0.08998329984):0.8154862223):0.5610961432):0.1617916438,(t242[&type=0]:0.5707751412,(t243[&type=0]:0.9836868793,t244[&type=0]:0.1984052949):0.6953297216):0.05552111682):0.9476150468):0.2473166997):0.9623488116,((t245[&type=0]:0.7935025664,t246[&type=0]:0.08509867964):0.3953444003,(t247[&type=0]:0.09163277131,(t248[&type=0]:0.5201428954,t249[&type=0]:0.8055520628):0.7452739514):0.3989078877):0.07581191277):0.9779064963,(((t250[&type=0]:0.943611098,(t251[&type=0]:0.33392801,t252[&type=0]:0.5996331484):0.4291575127):0.4906436009,((((t253[&type=0]:0.7749450852,(t254[&type=0]:0.8616885878,t255[&type=0]:0.585028409):0.06060880423):0.1238881133,((t256[&type=0]:0.7451687793,t257[&type=0]:0.6925335305):0.05338745634,t258[&type=0]:0.3357626374):0.2069296469):0.09644073155,((((t259[&type=0]:0.2258843291,t260[&type=0]:0.2671526412):0.3940743534,(t261[&type=0]:0.5022506947,(t262[&type=0]:0.9498897423,t263[&type=0]:0.1406114365):0.2847759123):0.04320593993):0.6982026948,t264[&type=0]:0.2693712024):0.959781138,(((t265[&type=0]:0.6035173486,t266[&type=0]:0.5529949202):0.9900399651,(t267[&type=0]:0.5455351078,t268[&type=0]:0.3530619899):0.4626278321):0.2735997427,(t269[&type=0]:0.9580646451,(t270[&type=0]:0.3280033092,t271[&type=0]:0.7206294278):0.03739526332):0.4967516926):0.9350089293):0.4371789068):0.1014483059,t272[&type=0]:0.2867298371):0.07522285799):0.06352435821,((t273[&type=0]:0.4001782183,t274[&type=0]:0.7190070178):0.1696753846,(t275[&type=0]:0.5535608665,t276[&type=0]:0.01324651297):0.2691543309):0.8676247413):0.8461736294):0.1769516913):0.344365149,(((t277[&type=0]:0.3245107541,(t278[&type=0]:0.4142541443,t279[&type=0]:0.5857141651):0.819547887):0.0867733527,(t280[&type=0]:0.4938162852,(t281[&type=0]:0.2444119717,t282[&type=0]:0.08141433029):0.05381231918):0.8375963389):0.176160393,((t283[&type=0]:0.4199601968,t284[&type=0]:0.8354801824):0.3150380594,(((t285[&type=0]:0.9818797186,(t286[&type=0]:0.8971825438,((t287[&type=0]:0.5155417006,t288[&type=0]:0.8260786769):0.7060374152,t289[&type=0]:0.6001661876):0.4120474763):0.9949228324):0.8038698458,t290[&type=0]:0.1939124272):0.6380942846,t291[&type=0]:0.3665255161):0.459349304):0.482901911):0.4833473735):0.5903116504):0.9973697898);",
        false);

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "processLength", new RealParameter(Double.toString(tree.getRoot().getHeight()+0.1)),
                "typeSet", new TypeSet(2),
                "R0", new SkylineVectorParameter(
                        null,
                        new RealParameter((4.0 / 3.0) + " " + 5.0)),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.5 1.25")),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealParameter((1.0 / 3.0) + " " + 0.5)),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.2 0.1")),

                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0"), 2));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("parameterization", parameterization,
                "frequencies", new RealParameter("0.5 0.5"),
                "conditionOnSurvival", true,
                "tree", tree,
                "typeLabel", "type",
                "parallelize", true);

		assertEquals(-661.9588648301033 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-5); // result from BEAST, not checked in R

	}

	/**
	 * Test infection among demes
	 * No rate changes
	 * Symmetric configuration
	 * reference from R
     */
	@Test
	public void testLikelihoodCalculationInfAmongDemesSymmetric() {

        // uncoloured, symmetric tree

        Tree tree = new TreeParser("((t3[&type=1]:0.004214277605,t4[&type=1]:0.02157681391):0.229186993,(t2[&type=0]:0.624713651,t1[&type=1]:1.347400211):0.06231047755);",
                false);

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "processLength", new RealParameter(Double.toString(tree.getRoot().getHeight() + 0.02686563367)),
                "typeSet", new TypeSet(2),
                "birthRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("2.0"), 2),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.5"), 2),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.5"), 2),
                "birthRateAmongDemes", new SkylineMatrixParameter(
                        null,
                        new RealParameter("1.0"), 2),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0"), 2));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("parameterization", parameterization,
                "frequencies", new RealParameter("0.5 0.5"),
                "conditionOnSurvival", true,
                "tree", tree,
                "typeLabel", "type",
                "parallelize", false);

        //System.out.println("Log-likelihood " + logL + " - testLikelihoodCalculationInfAmongDemes \t");
        assertEquals(-5.1966118470881 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-3);

    }

    /**
	 * Test infection among demes
	 * No rate changes
	 * Asymmetric configuration
	 * reference from R
     */
	@Test
	public void testLikelihoodCalculationInfAmongDemesAsymmetric() {

        Tree tree = new TreeParser("((3[&type=1]:1.5,4[&type=0]:0.5):1,(1[&type=0]:2,2[&type=1]:1):3);",
                false);

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "processLength", new RealParameter("6.0"),
                "typeSet", new TypeSet(2),
                "birthRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("2.0 6.25"), 2),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.2 0.625"), 2),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.3 0.625"), 2),
                "birthRateAmongDemes", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.2 0.1"), 2),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0"), 2));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("parameterization", parameterization,
                "frequencies", new RealParameter("0.5 0.5"),
                "conditionOnSurvival", true,
                "tree", tree,
                "typeLabel", "type",
                "parallelize", false);

		assertEquals(-26.7939 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-5);  //result from R
	}

	/**
	 * Test of migration and infection among demes with rate changes
	 * 2 types, no SA
	 * Uncoloured tree
	 * Reference from BDMM itself (version 0.2.0 28/06/2017)
	 * @throws Exception
	 */
	@Test
	public void testAmongRateChange() throws Exception {

        Tree tree = new TreeParser("((3[&type=0]:1.5,4[&type=1]:0.5):1,(1[&type=1]:1,2[&type=0]:1):3);",
                false);

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "processLength", new RealParameter("4.1"),
                "typeSet", new TypeSet(2),
                "R0", new SkylineVectorParameter(
                        new RealParameter("1.0"),
                        new RealParameter("6 5 2 2.5"), 2),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        new RealParameter("1.0"),
                        new RealParameter("0.5 0.55 0.45 0.6"), 2),
                "samplingProportion", new SkylineVectorParameter(
                        new RealParameter("1.0"),
                        new RealParameter("0.5 0.45 0.333333 0.35"), 2),
                "R0AmongDemes", new SkylineMatrixParameter(
                        new RealParameter("1.0"),
                        new RealParameter("1.1 1.3 1.2 1.15"), 2),
                "migrationRate", new SkylineMatrixParameter(
                        new RealParameter("1.0"),
                        new RealParameter("0.1 0.15 0.2 0.25")),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0"), 2));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("parameterization", parameterization,
                "frequencies", new RealParameter("0.5 0.5"),
                "conditionOnSurvival", false,
                "tree", tree,
                "typeLabel", "type",
                "parallelize", false);

        assertEquals(-16.466832439520886 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-4); // result from BDMM, 28/06/2017
    }

    /**
	 * Test of migration and infection among demes with rate changes
	 * 2 types, no SA
	 * Uncoloured tree
	 * Reference from BDMM itself (version 0.2.0 28/06/2017)
	 * @throws Exception
	 */
	@Test
	public void testAmongNoRateChange() throws Exception {

		Tree tree = new TreeParser("((3[&type=1]:1.5,4[&type=1]:0.5):1,(1[&type=1]:2,2[&type=1]:1):3);",
                false);

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "processLength", new RealParameter("6.0"),
                "typeSet", new TypeSet(2),
                "R0", new SkylineVectorParameter(
                        null,
                        new RealParameter("0 0"), 2),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0 0.75"), 2),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealParameter("0 0.7"), 2),
                "R0AmongDemes", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0 2"), 2),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.5 0")),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0"), 2));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("parameterization", parameterization,
                "frequencies", new RealParameter("1.0 0.0"),
                "conditionOnSurvival", false,
                "tree", tree,
                "typeLabel", "type",
                "parallelize", false);
		
		assertEquals(-12.1441 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-4); // tanja's result from R
		
	}
	
	/**
	 * Test of migration with 3 types
	 * No rate change, no SA
	 * Reference from BDMM itself 
	 * @throws Exception
	 */
	@Test
	public void testMig3types() throws Exception {

		Tree tree = new TreeParser("((3[&type=2]:1.5,4[&type=1]:0.5):1,(1[&type=1]:1,2[&type=0]:1):3);",
                false);

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "processLength", new RealParameter("4.1"),
                "typeSet", new TypeSet(3),
                "R0", new SkylineVectorParameter(
                        null,
                        new RealParameter("6 2 5")),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.5 0.45 0.55")),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.5 0.333333 0.45")),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.1 0.2 0.15 0.12 0.12 0.15")),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0"), 3));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("parameterization", parameterization,
                "frequencies", new RealParameter((1.0/3.0) + " " + (1.0/3.0) + " " + (1.0/3.0)),
                "conditionOnSurvival", false,
                "tree", tree,
                "typeLabel", "type",
                "parallelize", false);

		assertEquals(-16.88601100061662 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-4); // result from BDMM, version 0.2.0, 06/07/2017
	}

	/**
	 * Test tree with unknown states
	 * With and without migration, with and without rate changes
	 * 2 types, assymetric
     */
	@Test
	public void testUnknownStatesWithoutMigrationOrRateChange() {

        Tree tree = new TreeParser("((3[&type=\"?\"]:1.5,4[&type=\"?\"]:0.5):1,(1[&type=\"?\"]:1,2[&type=\"?\"]:1):3);",
                false);

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "processLength", new RealParameter("4.1"),
                "typeSet", new TypeSet(2),
                "R0", new SkylineVectorParameter(
                        null,
                        new RealParameter("6 2")),
                "R0AmongDemes", new SkylineMatrixParameter(
                        null,
                        new RealParameter("1 1")),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.5 1.0")),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.5 0.333333")),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0"), 2));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("parameterization", parameterization,
                "frequencies", new RealParameter("0.5 0.5"),
                "conditionOnSurvival", false,
                "tree", tree,
                "typeLabel", "type",
                "parallelize", false);

        //	System.out.println("Log-likelihood = " + logL);
        assertEquals(-18.82798 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-4); // tanja's result from R
    }

    /**
	 * Test tree with unknown states
	 * With and without migration, with and without rate changes
	 * 2 types, assymetric
     */
	@Test
	public void testUnknownStatesWithMigration() {

        Tree tree = new TreeParser("((3[&type=\"?\"]:1.5,4[&type=\"?\"]:0.5):1,(1[&type=\"?\"]:1,2[&type=\"?\"]:1):3);",
                false);

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "processLength", new RealParameter("4.1"),
                "typeSet", new TypeSet(2),
                "R0", new SkylineVectorParameter(
                        null,
                        new RealParameter("6 2")),
                "R0AmongDemes", new SkylineMatrixParameter(
                        null,
                        new RealParameter("1 1")),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.5 1.0")),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.5 0.333333")),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0"), 2),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.3 0.4")));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("parameterization", parameterization,
                "frequencies", new RealParameter("0.5 0.5"),
                "conditionOnSurvival", false,
                "tree", tree,
                "typeLabel", "type",
                "parallelize", false);

        assertEquals(-18.986212857895506 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-4); // reference from BDMM - 0.2.0 - 06/07/2017
    }

    /**
	 * Test tree with unknown states
	 * With and without migration, with and without rate changes
	 * 2 types, assymetric
     */
	@Test
	public void testUnknownStatesWithMigrationAndRateChange() {

        Tree tree = new TreeParser("((3[&type=\"?\"]:1.5,4[&type=\"?\"]:0.5):1,(1[&type=\"?\"]:1,2[&type=\"?\"]:1):3);",
                false);

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "processLength", new RealParameter("4.1"),
                "typeSet", new TypeSet(2),
                "R0", new SkylineVectorParameter(
                        new RealParameter("1.5"),
                        new RealParameter("6 5 2 1.5")),
                "R0AmongDemes", new SkylineMatrixParameter(
                        null,
                        new RealParameter("1 1.2")),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        new RealParameter("1.5"),
                        new RealParameter("0.5 1.0 1.0 0.5")),
                "samplingProportion", new SkylineVectorParameter(
                        new RealParameter("1.5"),
                        new RealParameter("0.5 0.45 0.333333 0.4")),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0"), 2),
                "migrationRate", new SkylineMatrixParameter(
                        new RealParameter("1.5"),
                        new RealParameter("0.3 0.35 0.4 0.32")));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("parameterization", parameterization,
                "frequencies", new RealParameter("0.5 0.5"),
                "conditionOnSurvival", false,
                "tree", tree,
                "typeLabel", "type",
                "parallelize", false);

		assertEquals(-17.87099909579358 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-4); // reference from BDMM - 0.2.0 - 06/07/2017
	}
	
	/**
	 * Test on combining migration with rho-sampling
	 * Reference from BDMM
     */
	@Test
	public void testLikelihoodMigrationRhoSampling() {

		Tree tree = new TreeParser("((1[&type=0]: 4.5, 2[&type=1]: 4.5):1,3[&type=0]:5.5);",
                false);

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "processLength", tree,
                "typeSet", new TypeSet(2),
                "R0", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.5 1.4")),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.5 1.3")),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.0"), 2),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0"), 2),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.3 0.4")),
                "rhoSampling", new TimedParameter(
                        new RealParameter("0.0"),
                        new RealParameter("0.01 0.015"),
                        tree));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("parameterization", parameterization,
                "frequencies", new RealParameter("0.6 0.4"),
                "conditionOnSurvival", false,
				"conditionOnRoot", true,
                "tree", tree,
                "typeLabel", "type",
                "parallelize", false);

		assertEquals(-8.906223150087108 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-4);   // Reference from BDMM - version 0.2.0 - 06/07/2017

	}
	
	/**
	 * Basic test on sampled-ancestors lik. calculation.
	 * 2 leaves, 1 SA. 1 type, no rho-sampling, no rate-change
	 * Reference value from BDSKY (23/03/2017)
     */
	@Test
	public void testSALikelihoodMini() {

		Tree tree = new TreeParser("((3[&type=0]: 1.5, 6[&type=0]: 0)5[&type=0]: 3.5, 4[&type=0]: 4) ;",
                false);

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "processLength", new RealParameter("6.0"),
                "typeSet", new TypeSet(1),
                "R0", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.5")),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.5")),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.2")),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.9")));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("parameterization", parameterization,
                "frequencies", new RealParameter("1.0"),
                "conditionOnSurvival", true,
                "tree", tree,
                "typeLabel", "type",
                "parallelize", false,
				"useAnalyticalSingleTypeSolution", false);

		assertEquals(-18.854438107814335 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-4); //Reference value from BDSKY (23/03/2017)

		density.setInputValue("useAnalyticalSingleTypeSolution", true);
		density.initAndValidate();

		assertEquals(-18.854438107814335 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-4); //Reference value from BDSKY (23/03/2017)
	}

	/**
	 * Basic test on sampled-ancestors lik. calculation with multi-rho sampling.
	 * 2 tips, 1 SA. 1 type, no rate-change
	 * Coloured and uncoloured trees
	 * Reference value from BDSKY (06/04/2017)
	 * @throws Exception
	 */
	@Test
	public void testSALikelihoodMultiRho() throws Exception {

		Tree tree = new TreeParser("((3[&type=0]: 1.5, 6[&type=0]: 0)5[&type=0]: 3.5, 4[&type=0]: 4) ;",false);

		RealParameter origin = new RealParameter("6.0");

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "processLength", origin,
                "typeSet", new TypeSet(1),
                "R0", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.5")),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.5")),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.2")),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.9")),
                "rhoSampling", new TimedParameter(
                        new RealParameter("0.0 1.5"),
                        new RealParameter("0.05 0.3"),
                        origin));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("parameterization", parameterization,
                "frequencies", new RealParameter("1.0"),
                "conditionOnSurvival", true,
                "tree", tree,
                "typeLabel", "type",
                "parallelize", false,
				"useAnalyticalSingleTypeSolution", false);

		assertEquals(-22.348462265673483 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-5); //Reference value from BDSKY (06/04/2017)

		density.setInputValue("useAnalyticalSingleTypeSolution", true);
		density.initAndValidate();

		assertEquals(-22.348462265673483 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-5); //Reference value from BDSKY (06/04/2017)
	}

	/**
	 * Test on sampled-ancestors lik. calculation with no sampled ancestor
	 * No rate-change, one state, 4 tips
	 * This state is just there in case something is broken with sampled ancestors,
	 * helps for debugging if combined with testSALikelihoodMini for instance
	 * @throws Exception
	 */
	@Test
	public void testSALikelihoodCalculationWithoutAncestors() throws Exception {

	    Tree tree = new TreeParser("((3[&type=0] : 1.5, 4[&type=0] : 0.5) : 1 , (1[&type=0] : 2, 2[&type=0] : 1) : 3);",
                false);

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "processLength", tree,
                "typeSet", new TypeSet(1),
                "R0", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.5")),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.5")),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.3")),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.9")));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("parameterization", parameterization,
                "frequencies", new RealParameter("1.0"),
                "conditionOnSurvival", true,
				"conditionOnRoot", true,
                "tree", tree,
                "typeLabel", "type",
                "parallelize", false,
				"useAnalyticalSingleTypeSolution", false);

        // Conditioned on root:

		assertEquals(-15.99699690815937 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-4);

		density.setInputValue("useAnalyticalSingleTypeSolution", true);
		density.initAndValidate();

		assertEquals(-15.99699690815937 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-4);

		// Conditioned on origin:

		parameterization.setInputValue("processLength", new RealParameter("10.0"));
		parameterization.initAndValidate();
		density.setInputValue("useAnalyticalSingleTypeSolution", false);
		density.setInputValue("conditionOnRoot", false);
		density.initAndValidate();

		assertEquals(-25.991511346557598 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-4);

		density.setInputValue("useAnalyticalSingleTypeSolution", true);
		density.initAndValidate();

		assertEquals(-25.991511346557598 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-4);
	}
}