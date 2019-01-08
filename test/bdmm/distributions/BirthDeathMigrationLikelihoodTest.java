package bdmm.distributions;

import bdmm.parameterization.*;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;

import static junit.framework.Assert.assertEquals;


/**
 * Created by Jeremie Scire (jscire) on 26.06.17.
 */

public class BirthDeathMigrationLikelihoodTest {

	double runtime;

	/**
	 * Basic test for migration rate change 
	 * Coloured and uncoloured trees 
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
		        "nTypes", 2,
                "origin", originParam,
                "birthRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("2.0 2.0")),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0 1.0")),
                "birthRateAmongDemes", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.0 0.0")),
                "migrationRate", new SkylineMatrixParameter(
                        new RealParameter("1.0"),
                        new RealParameter("0.1 0.1 0.2 0.2")),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.5 0.5")),
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
	 * Basic test for migration rate change
	 * Coloured and uncoloured trees
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
		        "nTypes", 2,
                "origin", originParam,
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
	 * Coloured and uncoloured trees 
	 * Reference from BDMM itself
	 * @throws Exception
	 */
	@Test 
	public void testLikelihoodRemovalProbChangeBasic() {

		String newick = "((1[&type=0]: 1.5, 2[&type=0]: 0)3[&type=0]: 3.5, 4[&type=0]: 4) ;";

		Parameterization parameterization = new EpiParameterization();
		parameterization.initByName(
		        "nTypes", 1,
                "origin", new RealParameter("6.0"),
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
                "parallelize", false);

		double logL = density.calculateLogP();

        // Reference BDMM (version 0.2.0) 22/06/2017
		assertEquals(-21.25413884159791, logL, 1e-5);
	}

	/**
	 * Two-state test for removal-probability rate change
	 * Coloured and uncoloured trees
	 * Reference from BDMM itself
	 * @throws Exception
	 */
	@Test
	public void testLikelihoodRemovalProbChangeTwoState() throws Exception{

		String newick = "((1[&type=0]: 1.5, 2[&type=1]: 0)3[&type=0]: 3.5, 4[&type=1]: 4) ;";

        Parameterization parameterization = new EpiParameterization();
		parameterization.initByName(
		        "nTypes", 2,
                "origin", new RealParameter("6.0"),
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
		assertEquals(-21.185194919464568, logL, 1e-5);
	}


	/**
	 * Test simple configuration with one rho-sampling event for
	 * coloured and uncoloured trees.
	 * One type, no psi-sampling, no sampled-ancestor
	 * No rate-changes
	 * Reference: R
	 * @throws Exception
	 */
	@Test
	public void testSingleRho() throws Exception {

		// Uncoloured-tree test cases
		Tree tree = new TreeParser("((1[&type=0]: 4.5, 2[&type=0]: 4.5):1,3[&type=0]:5.5);",false);

        Parameterization parameterization = new EpiParameterization();
		parameterization.initByName(
		        "nTypes", 1,
                "tree", tree,
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
                "tree", tree,
                "typeLabel", "type",
                "parallelize", false);

		double logL = density.calculateLogP();

		// this result is from R: LikConstant(2.25,1.5,0.01,c(4.5,5.5),root=1,survival=0)
		assertEquals(-6.761909, logL, 1e-4);

		// test with conditioned-on-survival tree
		density.setInputValue("conditionOnSurvival", true);
		density.setInputValue("origin", "10");
		density.initAndValidate();

		double logL2 = density.calculateLogP();

        // this result is from R: LikConstant(2.25,1.5,0.01,c(4.5,5.5,5.5+1e-100),root=0,survival=1)
		assertEquals(-7.404227, logL2, 1e-4);
	}

	/**
	 * Basic test for rho-sampling in the past
	 * Coloured and Uncoloured trees
	 * One type, 2 tips, one state
	 * No rate changes
	 * @throws Exception
	 */
	@Test
	public void testMultiRho2tips() throws Exception {
		// Uncoloured Tree

		// two tips sampled at the same time
		Tree tree = new TreeParser("(3[&type=0]: 4, 4[&type=0]: 4) ;",false);

		BirthDeathMigrationDistribution bdm =  new BirthDeathMigrationDistribution();

		bdm.setInputValue("tree", tree);
		bdm.setInputValue("typeLabel", "type");

		bdm.setInputValue("stateNumber", "1");
		bdm.setInputValue("migrationMatrix", "0.");
		bdm.setInputValue("frequencies", "1");

		bdm.setInputValue("R0", new RealParameter("1.5"));
		bdm.setInputValue("becomeUninfectiousRate", new RealParameter("1.5"));
		bdm.setInputValue("samplingProportion", new RealParameter("0.") );

		bdm.setInputValue("rho", new RealParameter("0.2 1.") );
		bdm.setInputValue("rhoSamplingTimes", new RealParameter("0. 2.5") );
		bdm.setInputValue("reverseTimeArrays", "false false false true");
		bdm.setInputValue("conditionOnSurvival", true);
		bdm.setInputValue("origin", "5.");

		bdm.initAndValidate();

		// this result is from BEAST: BDSKY, not double checked in R
		assertEquals(-10.569863754307026, bdm.calculateLogP(), 1e-4);

		// tips sampled at two different times
		tree = new TreeParser("(3[&type=0]: 1.5, 4[&type=0]: 4) ;",false);
		bdm.setInputValue("tree", tree);
		bdm.initAndValidate();

		// this result is from BEAST: BDSKY, not double checked in R
		assertEquals(-8.099631076932816, bdm.calculateLogP(), 1e-4);
	}

	/**
	 * Test with combined multiple-rho-sampling events in the past and psi-sampling
	 * was "TestRhoSasha"
	 * Coloured and Uncoloured trees
	 * One type, no rate-changes, no sampled-ancestors
	 * 26 tips
	 * @throws Exception
	 */
	@Test
	public void testMultiRhoSampling() throws Exception {
		// Uncoloured tree
		Tree treeU = new TreeParser("(((((t1[&type=0]:0.4595008531,t25[&type=0]:0.4595008531)[&type=0]:0.3373053072,t23[&type=0]:0.3567584538)[&type=0]:0.007310819036,t16[&type=0]:0.3489190732)[&type=0]:0.331009529,((t18[&type=0]:0.03315384045,t14[&type=0]:0.03315384045)[&type=0]:0.5063451374,(t10[&type=0]:0.4211543131,t15[&type=0]:0.4211543131)[&type=0]:0.1183446648)[&type=0]:0.5956275305)[&type=0]:0.1158090878,((t19[&type=0]:0.9429393194,((t6[&type=0]:0.363527235,t11[&type=0]:0.4417423167)[&type=0]:0.01881829549,((((t3[&type=0]:0.3071904376,(((t24[&type=0]:0.01065209364,t13[&type=0]:0.01065209364)[&type=0]:0.06076485145,t8[&type=0]:0.07141694509)[&type=0]:0.123620245,(t22[&type=0]:0.1616119808,t2[&type=0]:0.1616119808)[&type=0]:0.03342520927)[&type=0]:0.1121532475)[&type=0]:0.24520579,t9[&type=0]:0.5523962276)[&type=0]:0.3852615426,(((t20[&type=0]:0.2935970782,(t17[&type=0]:0.06569090089,t4[&type=0]:0.06569090089)[&type=0]:0.2279061773)[&type=0]:0.08350780408,(t21[&type=0]:0.05109047139,t5[&type=0]:0.05109047139)[&type=0]:0.3260144109)[&type=0]:0.2298344132,t7[&type=0]:0.6069392955)[&type=0]:0.3307184747)[&type=0]:0.01206284377,t26[&type=0]:0.9497206139)[&type=0]:0.05755333197)[&type=0]:0.03290891884)[&type=0]:0.07263755325,t12[&type=0]:1.112820418)[&type=0]:0.1381151782);",false);

		BirthDeathMigrationDistribution bdssm =  new BirthDeathMigrationDistribution();
		bdssm.setInputValue("typeLabel", "type");
		bdssm.setInputValue("frequencies", "1");
		bdssm.setInputValue("migrationMatrix", "0.");
		bdssm.setInputValue("stateNumber", 1);

		bdssm.setInputValue("tree", treeU);
		bdssm.setInputValue("origin", new RealParameter("2."));
		bdssm.setInputValue("conditionOnSurvival", false);

		bdssm.setInputValue("R0", new RealParameter(new Double[]{3./4.5}));
		bdssm.setInputValue("becomeUninfectiousRate", new RealParameter("4.5"));
		bdssm.setInputValue("samplingProportion", new RealParameter(new Double[]{2./4.5}));
		bdssm.setInputValue("rho", new RealParameter("0.0 0.05 0.01"));
		bdssm.setInputValue("rhoSamplingTimes","0. 1. 1.5");
		bdssm.setInputValue("reverseTimeArrays","false false false false");
		bdssm.initAndValidate();

		assertEquals(-124.96086690757612, bdssm.calculateTreeLogLikelihood(treeU), 1e-2);     // this result is from BEAST, not double checked in R

		// same with reverse rhoSamplingTimes
		bdssm.setInputValue("rhoSamplingTimes","0. 0.5 1.");
		bdssm.setInputValue("reverseTimeArrays","false false false true");
		bdssm.initAndValidate();

		assertEquals(-124.96086690757612, bdssm.calculateTreeLogLikelihood(treeU), 1e-2);     // this result is from BEAST, not double checked in R
	}

	/**
	 * Test with multiple cases for rho-sampling in the past combined with rate changes
	 * Coloured and Uncoloured trees
	 * 1 state, no sampled ancestors
	 * 26 tips
	 * @throws Exception
	 */
	@Test
	public void testMultiRhoWithRateChanges() throws Exception {

		// uncoloured trees
		Tree tree = new TreeParser("(((((t1[&type=0]:0.4595008531,t25[&type=0]:0.4595008531)[&type=0]:0.3373053072,t23[&type=0]:0.3567584538)[&type=0]:0.007310819036,t16[&type=0]:0.3489190732)[&type=0]:0.331009529,((t18[&type=0]:0.03315384045,t14[&type=0]:0.03315384045)[&type=0]:0.5063451374,(t10[&type=0]:0.4211543131,t15[&type=0]:0.4211543131)[&type=0]:0.1183446648)[&type=0]:0.5956275305)[&type=0]:0.1158090878,((t19[&type=0]:0.9429393194,((t6[&type=0]:0.363527235,t11[&type=0]:0.4417423167)[&type=0]:0.01881829549,((((t3[&type=0]:0.3071904376,(((t24[&type=0]:0.01065209364,t13[&type=0]:0.01065209364)[&type=0]:0.06076485145,t8[&type=0]:0.07141694509)[&type=0]:0.123620245,(t22[&type=0]:0.1616119808,t2[&type=0]:0.1616119808)[&type=0]:0.03342520927)[&type=0]:0.1121532475)[&type=0]:0.24520579,t9[&type=0]:0.5523962276)[&type=0]:0.3852615426,(((t20[&type=0]:0.2935970782,(t17[&type=0]:0.06569090089,t4[&type=0]:0.06569090089)[&type=0]:0.2279061773)[&type=0]:0.08350780408,(t21[&type=0]:0.05109047139,t5[&type=0]:0.05109047139)[&type=0]:0.3260144109)[&type=0]:0.2298344132,t7[&type=0]:0.6069392955)[&type=0]:0.3307184747)[&type=0]:0.01206284377,t26[&type=0]:0.9497206139)[&type=0]:0.05755333197)[&type=0]:0.03290891884)[&type=0]:0.07263755325,t12[&type=0]:1.112820418)[&type=0]:0.1381151782);",false);

		for (int i = 0; i<4; i++){

			switch (i){
			case 0:{ // no rate-change, rho-sampling at present
				BirthDeathMigrationDistribution bdssm =  new BirthDeathMigrationDistribution();
				bdssm.setInputValue("typeLabel", "type");
				bdssm.setInputValue("frequencies", "1");
				bdssm.setInputValue("migrationMatrix", "0.");
				bdssm.setInputValue("stateNumber", 1);

				bdssm.setInputValue("tree", new TreeParser("(((((t1[&type=0]:0.4595008531,t25[&type=0]:0.4595008531)[&type=0]:0.3373053072,t23[&type=0]:0.3567584538)[&type=0]:0.007310819036,t16[&type=0]:0.3489190732)[&type=0]:0.331009529,((t18[&type=0]:0.03315384045,t14[&type=0]:0.03315384045)[&type=0]:0.5063451374,(t10[&type=0]:0.4211543131,t15[&type=0]:0.4211543131)[&type=0]:0.1183446648)[&type=0]:0.5956275305)[&type=0]:0.1158090878,((t19[&type=0]:0.9429393194,((t6[&type=0]:0.363527235,t11[&type=0]:0.4417423167)[&type=0]:0.01881829549,((((t3[&type=0]:0.3071904376,(((t24[&type=0]:0.01065209364,t13[&type=0]:0.01065209364)[&type=0]:0.06076485145,t8[&type=0]:0.07141694509)[&type=0]:0.123620245,(t22[&type=0]:0.1616119808,t2[&type=0]:0.1616119808)[&type=0]:0.03342520927)[&type=0]:0.1121532475)[&type=0]:0.24520579,t9[&type=0]:0.5523962276)[&type=0]:0.3852615426,(((t20[&type=0]:0.2935970782,(t17[&type=0]:0.06569090089,t4[&type=0]:0.06569090089)[&type=0]:0.2279061773)[&type=0]:0.08350780408,(t21[&type=0]:0.05109047139,t5[&type=0]:0.05109047139)[&type=0]:0.3260144109)[&type=0]:0.2298344132,t7[&type=0]:0.6069392955)[&type=0]:0.3307184747)[&type=0]:0.01206284377,t26[&type=0]:0.9497206139)[&type=0]:0.05755333197)[&type=0]:0.03290891884)[&type=0]:0.07263755325,t12[&type=0]:1.112820418)[&type=0]:0.1381151782);", false));
				bdssm.setInputValue("origin", "2.");
				bdssm.setInputValue("conditionOnSurvival", true);
				bdssm.setInputValue("birthRate", new RealParameter("2."));
				bdssm.setInputValue("deathRate", new RealParameter("0.5"));
				bdssm.setInputValue("samplingRate", new RealParameter("0.5"));

				bdssm.setInputValue("rho", new RealParameter("1."));
				bdssm.initAndValidate();


				//		System.out.println("\na) Likelihood: " + bdssm.calculateTreeLogLikelihood(tree));
				assertEquals(-21.42666177086957, bdssm.calculateTreeLogLikelihood(tree), 1e-5);

			}
			case 1:{ // rate-changes, rho-sampling in the past
				BirthDeathMigrationDistribution bdssm =  new BirthDeathMigrationDistribution();
				bdssm.setInputValue("typeLabel", "type");
				bdssm.setInputValue("frequencies", "1");
				bdssm.setInputValue("migrationMatrix", "0.");
				bdssm.setInputValue("stateNumber", 1);

				bdssm.setInputValue("tree", tree);
				bdssm.setInputValue("origin", new RealParameter("2."));
				bdssm.setInputValue("conditionOnSurvival", false);

				bdssm.setInputValue("R0", new RealParameter(new Double[]{3./4.5, 2./1.5, 4./1.5}));   // birthRate = "3. 2. 4."
				bdssm.setInputValue("becomeUninfectiousRate", new RealParameter("4.5 1.5 1.5"));      // deathRate = "2.5 1. .5"
				bdssm.setInputValue("samplingProportion", new RealParameter(new Double[]{2./4.5, .5/1.5, 1./1.5}));                  // samplingRate = "2. 0.5 1."
				bdssm.setInputValue("rho", new RealParameter("0. 0. 0.01"));
				bdssm.setInputValue("birthRateChangeTimes", new RealParameter("0. 1. 1.5"));
				bdssm.setInputValue("deathRateChangeTimes", new RealParameter("0. 1. 1.5"));
				bdssm.setInputValue("samplingRateChangeTimes", new RealParameter("0. 1. 1.5"));
				bdssm.setInputValue("rhoSamplingTimes", new RealParameter("0. 1. 1.5"));

				bdssm.initAndValidate();


				assertEquals(-87.59718586549747, bdssm.calculateTreeLogLikelihood(tree), 1e-1);

			}
			case 2:{ // rate-changes, rho-sampling in the past and present
				BirthDeathMigrationDistribution bdssm =  new BirthDeathMigrationDistribution();
				bdssm.setInputValue("typeLabel", "type");
				bdssm.setInputValue("frequencies", "1");
				bdssm.setInputValue("migrationMatrix", "0.");
				bdssm.setInputValue("stateNumber", 1);

				bdssm.setInputValue("tree", tree);
				bdssm.setInputValue("origin", new RealParameter("2."));
				bdssm.setInputValue("conditionOnSurvival", false);

				bdssm.setInputValue("R0", new RealParameter(new Double[]{3./4.5, 2./1.5, 4./1.5}));   // birthRate = "3. 2. 4."
				bdssm.setInputValue("becomeUninfectiousRate", new RealParameter("4.5 1.5 1.5"));      // deathRate = "2.5 1. .5"
				bdssm.setInputValue("samplingProportion", new RealParameter(new Double[]{2./4.5, .5/1.5, 1./1.5}));                  // samplingRate = "2. 0.5 1."
				bdssm.setInputValue("birthRateChangeTimes", new RealParameter("0. 1. 1.5"));
				bdssm.setInputValue("deathRateChangeTimes", new RealParameter("0. 1. 1.5"));
				bdssm.setInputValue("samplingRateChangeTimes", new RealParameter("0. 1. 1.5"));

				bdssm.setInputValue("rho", new RealParameter("0.05 0.01"));
				bdssm.setInputValue("rhoSamplingTimes","0. 1.");
				bdssm.initAndValidate();


				assertEquals(-87.96488, bdssm.calculateTreeLogLikelihood(tree), 1e-1);
			}

			case 3:{ // rate-changes, rho-sampling in the past and present, with reversed times
				BirthDeathMigrationDistribution bdssm =  new BirthDeathMigrationDistribution();
				bdssm.setInputValue("typeLabel", "type");
				bdssm.setInputValue("frequencies", "1");
				bdssm.setInputValue("migrationMatrix", "0.");
				bdssm.setInputValue("stateNumber", 1);

				bdssm.setInputValue("tree", tree);
				bdssm.setInputValue("conditionOnSurvival", false);

				bdssm.setInputValue("R0", new RealParameter(new Double[]{3./4.5, 2./1.5, 4./1.5, 4./2.5}));   // birthRate = "3. 2. 4. 4."
				bdssm.setInputValue("becomeUninfectiousRate", new RealParameter("4.5 1.5 1.5 2.5"));      // deathRate = "2.5 1. .5 .5"
				bdssm.setInputValue("samplingProportion", new RealParameter(new Double[]{2./4.5, .5/1.5, 1./1.5, 2./2.5}));                  // samplingRate = "2. 0.5 1. 2."
				bdssm.setInputValue("rho", new RealParameter("0.05 0.01"));
				bdssm.setInputValue("rhoSamplingTimes","0. 1.");
				bdssm.setInputValue("intervalTimes", new RealParameter("0. 0.5 1. 1.1"));
				bdssm.initAndValidate();


				assertEquals(-100.15682190617582, bdssm.calculateTreeLogLikelihood(tree), 1e-1);
			}
			}
		}
	}

	/**
	 * Basic 1-dim test
	 * No rate change, 1 state, no rho-sampling
	 * Coloured and Uncoloured trees
	 * Reference from BDSKY
	 * @throws Exception
	 */
	@Test
	public void testLikelihood1dim() throws Exception {

		// Test conditions:
		// 4 tips
		// 1 state
		// 0 migration
		// frequency 1 (unique state)
		// origin at t = -6
		// R0 = 1.33...
		// becomeUninfectiousRate = 1.5
		// samplingProportion (psi-sampling) = 0.333...
		// no rho-sampling
		// no rate changes
		// conditionOnSurvival: false or true

		String tree = "((3[&state=0] : 1.5, 4[&state=0] : 0.5)[&state=0] : 1 , (1[&state=0] : 2, 2[&state=0] : 1)[&state=0] : 3)[&state=0];";

		int nbOfTips = 4;
		String nbOfStates = "1";
		String frequencies = "1.";
		String migration = "0.";
		String origin = "6.";
		String R0  = "1.3333333334";
		String becomeUninfectiousRate = "1.5";
		String samplingProportion = "0.33333333333";

		/////// Uncoloured Tree
		String locations = "1=0,2=0,3=0,4=0" ;
		boolean conditionOnSurvival = false;

		double logL2 = bdm_likelihood( nbOfStates,
				migration,
				frequencies,
				tree, "1.", // origin is defined at 1. instead of 6. because bdm_likelihood adds the height of the root of the tree to that value (here 5.)
				R0,null,
				becomeUninfectiousRate,
				samplingProportion,  null,
				"", locations, nbOfTips, null, conditionOnSurvival);

		System.out.println("Birth-death result 2: " +logL2 + "\t- Test Likelihood1dim");

		assertEquals(-19.019796073623493, logL2, 1e-5);   // Reference BDSKY (version 1.3.3)

		conditionOnSurvival = true;

		double logL2bis = bdm_likelihood( nbOfStates,
				migration,
				frequencies,
				tree, "1.", // origin is defined at 1. instead of 6. because bdm_likelihood adds the height of the root of the tree to that value (here 5.)
				R0,null,
				becomeUninfectiousRate,
				samplingProportion, null,
				"", locations, nbOfTips, null, conditionOnSurvival);

		System.out.println("Birth-death result 2 conditionned on survival: " +logL2bis + "\t- Test Likelihood1dim");

		assertEquals(-18.574104140202046, logL2bis, 1e-5); // Reference BDSKY (version 1.3.3)
	}

	/**
	 * 1-dim and 1 rate-change test
	 * coloured and uncoloured trees
	 * reference from BDSKY
	 * @throws Exception
	 */
	@Test
	public void testLikelihoodRateChange1dim() throws Exception {

		// Test conditions:
		// 4 tips
		// 1 state
		// 0 migration
		// frequency 1 (unique state)
		// origin at t = 1 + rootoftree
		// R0 = 1.33...
		// becomeUninfectiousRate = 1.5
		// samplingProportion (psi-sampling) = 0.44444.. / 0.33333..
		// no rho-sampling
		// 1 rate change
		// conditionOnSurvival = false

		String tree = "((3[&state=0] : 1.5, 4[&state=0] : 0.5)[&state=0] : 1 , (1[&state=0] : 2, 2[&state=0] : 1)[&state=0] : 3)[&state=0];";

		boolean conditionOnSurvival = false;

		int nbOfTips = 4;
		String nbOfStates = "1";
		String frequencies = "1.";
		String migration = "0.";
		String R0  = "0.6666666667 1.3333333334";
		String becomeUninfectiousRate = "4.5 1.5";
		String samplingProportion = "0.4444444444 0.33333333333";
		String locations = "1=0,2=0,3=0,4=0" ;
		String intervalTimes = "0. 3.";

		// uncoloured trees

		double logL = bdm_likelihood( nbOfStates,
				migration,
				frequencies,
				tree, "1.", // origin is defined at 1. instead of 6. because bdm_likelihood adds the height of the root of the tree to that value (here 5.)
				R0,null,
				becomeUninfectiousRate,
				samplingProportion, null,
				"", locations, nbOfTips, intervalTimes, conditionOnSurvival);

		System.out.println("Birth-death result (uncoloured tree): " + logL + "\t- Test LikelihoodRateChange1dim");

		assertEquals(-33.7573, logL, 1e-4); // Reference BDSKY

	}

	/**
	 * Basic tests on 2 types situations with migration or birth among demes
	 * uncoloured trees
	 * reference from R
	 * @throws Exception
	 */
	@Test 
	public void testLikelihoodCalculationMigTiny() throws Exception {

		String tree ="(1[&state=0] : 1.5, 2[&state=1] : 0.5)[&state=0];";

		String orig="1.";
		String stateNumber = "2";
		String migrationMatrix = ".1 .1";
		String frequencies = "0.5 0.5";
		String R0 = Double.toString(4./3.) + " " + Double.toString(4./3.) ;
		String becomeUninfectiousRate = "1.5 1.5";
		String samplingProportion = Double.toString(1./3.) + " " + Double.toString(1./3.);
		String locations = "1=0,2=1" ;
		boolean conditionOnSurvival = false;

		// migration and no birth among demes
		double logL = bdm_likelihood(stateNumber,
				migrationMatrix,
				frequencies,
				tree, orig,
				R0, null,
				becomeUninfectiousRate,
				samplingProportion, null,
				"",
				locations,
				2, null, conditionOnSurvival);

		assertEquals(-7.215222, logL, 1e-6); // result from R

		// no migration, symmetric birth among demes
		String R0AmongDemes = "0.0666667 0.0666667" ;

		logL = bdm_likelihood(stateNumber,
				"0 0", // no migration
				frequencies,
				tree, orig,
				R0, R0AmongDemes,
				becomeUninfectiousRate,
				samplingProportion, null,
				"",
				locations,
				2, null,
				conditionOnSurvival);

		assertEquals(-7.404888, logL, 1e-6); // result from R

		// no migration, asymmetric birth among demes
		R0AmongDemes = "0.0666667 0.1" ;
		logL = bdm_likelihood(stateNumber,
				"0 0",
				frequencies,
				tree, orig,
				R0, R0AmongDemes,
				becomeUninfectiousRate,
				samplingProportion, null,
				"",
				locations,
				2, null,
				conditionOnSurvival);

		//  System.out.println("Log-likelihood (not conditioned on survival) = " + logL);

		assertEquals(-7.18723, logL, 1e-6); // result from R

		R0 = "2 1.3333333" ;
		logL = bdm_likelihood(stateNumber,
				"0 0",
				frequencies,
				tree, orig,
				R0, R0AmongDemes,
				becomeUninfectiousRate,
				samplingProportion, null,
				"",
				locations,
				2, null,
				conditionOnSurvival);

		//         System.out.println("Log-likelihood (not conditioned on survival) = " + logL);

		assertEquals(-7.350649, logL, 1e-6); // result from R

		R0 = "2 1.5" ;
		becomeUninfectiousRate = "2 1" ;
		samplingProportion = "0.5 0.3";
		R0AmongDemes = "0.1 0.5" ;

		logL = bdm_likelihood(stateNumber,
				"0 0", // no migration
				frequencies,
				tree, orig,
				R0, R0AmongDemes,
				becomeUninfectiousRate,
				samplingProportion, null,
				"",
				locations,
				2, null,
				conditionOnSurvival);

		//         System.out.println("Log-likelihood (not conditioned on survival) = " + logL);
		assertEquals(-6.504139, logL, 1e-6); // result from R

		locations = "1=1,2=0";
		logL = bdm_likelihood(stateNumber,
				"0 0", // no migration
				frequencies,
				tree, orig,
				R0, R0AmongDemes,
				becomeUninfectiousRate,
				samplingProportion, null,
				"",
				locations,
				2, null,
				conditionOnSurvival);

		//         System.out.println("Log-likelihood (not conditioned on survival) = " + logL);
		assertEquals(-7.700916, logL, 1e-6); // result from R
	}

	/**
	 * Test migration
	 * 2 types, migration, no birth among demes
	 * small and bigger trees
	 * Adapted from BDSKY
	 * Coloured and Uncoloured trees
	 * @throws Exception
	 */
	@Test
	public void testLikelihoodCalculationMig() throws Exception {

		// uncoloured tree, asymmetric types
		String tree ="((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);"; //
		String orig="1."; //
		String stateNumber = "2";
		String migrationMatrix = ".2 .1";
		String frequencies = "0.5 0.5";

		String R0 = Double.toString(4./3.) + " " + Double.toString(5.) ;
		String becomeUninfectiousRate = "1.5 1.25";
		String samplingProportion = Double.toString(1./3.) + " " + Double.toString(1./2.);

		String locations = "1=1,2=0,3=0,4=1" ;

		boolean conditionOnSurvival = false;

		double logL = bdm_likelihood(stateNumber,
				migrationMatrix,
				frequencies,
				tree, orig,
				R0, null,
				becomeUninfectiousRate,
				samplingProportion, null,
				"", locations, 4, null,
				conditionOnSurvival);

		//		System.out.println("Log-likelihood (NOT conditioned on survival) = " + logL);
		assertEquals(-26.53293, logL, 1e-5);

		// uncoloured tree, 291 tips

		String treeBig = "(((((((t1:0.9803361397,t2:0.9035540882):0.0532383481,t3:0.2637392259):0.6273536528,(t4:0.8624112266,t5:0.3278892266):0.2606245542):0.2941323873,(t6:0.09820114588,t7:0.533115675):0.8625875909):0.7040311908,(((t8:0.8696136218,t9:0.08719484485):0.4204288905,(t10:0.102143287,(t11:0.9850614571,t12:0.7407912319):0.8715072596):0.5182644848):0.524062254,(((((((t13:0.3981794417,(t14:0.03889928572,t15:0.5187105467):0.1127638209):0.3431177251,((t16:0.4239511855,t17:0.001895790454):0.690600364,t18:0.6283850113):0.4073564562):0.6862231812,(((t19:0.9947085041,t20:0.4739363373):0.1873670686,t21:0.151270482):0.803061039,((t22:0.8899249982,((t23:0.1329096023,t24:0.84205155):0.8838408566,(t25:0.7541888549,t26:0.8602364615):0.8912267659):0.771449636):0.1022819551,(((t27:0.3134289116,(t28:0.2446750235,t29:0.8565168788):0.8277210968):0.4307989818,((t30:0.2330717787,t31:0.4438336496):0.6521712865,(t32:0.2534400895,t33:0.7885409284):0.3051449039):0.1196702593):0.4061951274,t34:0.8415271267):0.4365981282):0.753448925):0.1580670979):0.04210642632,(((t35:0.7504386581,t36:0.6328390085):0.9047614154,t37:0.4946133171):0.2264722914,((((t38:0.06683212146,t39:0.479845396):0.9424520086,t40:0.894530142):0.3844042511,(((t41:0.5215392481,t42:0.2366602973):0.8142298241,(t43:0.2968777204,(t44:0.655541793,t45:0.8608812049):0.3564132168):0.04912991729):0.1511388237,t46:0.9031036345):0.1874918914):0.9690212663,(t47:0.07753491728,(t48:0.8349514075,(t49:0.9689748741,t50:0.925813166):0.4534903264):0.3571097804):0.1324767114):0.5515443345):0.3330309158):0.7202291801,((t51:0.6977306763,((t52:0.9157640305,t53:0.4226291834):0.5872618856,t54:0.2063144948):0.1422286083):0.7182746637,t55:0.759545143):0.7437628019):0.2425582204,((t56:0.4614429038,(t57:0.9092229386,((t58:0.1049408391,t59:0.6328130178):0.642241966,((t60:0.264340204,t61:0.5904771155):0.7333205172,(t62:0.9183179205,t63:0.1090340314):0.3010568973):0.3240860389):0.3192155454):0.1835780439):0.5942421539,t64:0.7931551472):0.967891278):0.06263663713,(t65:0.5774453548,((t66:0.07208712469,((t67:0.8918803469,t68:0.5110983853):0.1491188321,t69:0.2471361952):0.9591872343):0.3133718621,(t70:0.944087367,t71:0.7830825299):0.2284035049):0.5492361034):0.1136150162):0.002181729767):0.4548798562):0.4258609388,((((((t72:0.27679418,t73:0.5398862793):0.8871422287,(((((t74:0.2531923286,t75:0.3796772889):0.4489221217,t76:0.2554209188):0.3248268673,t77:0.5372577759):0.5699883625,t78:0.1656995732):0.957750936,(t79:0.1301121258,t80:0.8925942327):0.2838441601):0.5258686764):0.47825964,(t81:0.5749240227,((t82:0.9574132746,(t83:0.00485483068,t84:0.8091488208):0.1985368489):0.3703975577,(((t85:0.3991035291,(t86:0.03201846033,t87:0.8380640063):0.05616304209):0.8414494572,t88:0.6844437125):0.2426782607,((t89:0.7543559887,t90:0.7162597755):0.8230077426,t91:0.08967904118):0.4460245941):0.8679371702):0.51572948):0.4362259945):0.2631344711,(((t92:0.3353162925,((t93:0.4025212794,t94:0.0281926766):0.7965471447,t95:0.1145715592):0.5993301494):0.08854756854,(t96:0.1461353719,((t97:0.3158547124,t98:0.06653800653):0.5634025722,t99:0.9711292514):0.9727503664):0.7684133062):0.4824229684,((t100:0.06834940333,t101:0.7794982188):0.3453287922,(t102:0.627945075,t103:0.1914187325):0.9974814849):0.6312927424):0.04858242651):0.2845227425,((t104:0.6782600286,(t105:0.03190574702,t106:0.5840284519):0.03041352634):0.725893975,(((t107:0.9885271091,t108:0.07126446022):0.8419693699,t109:0.1546431775):0.898004594,t110:0.2500803664):0.1493327522):0.4266726137):0.5946582041,(t111:0.1395377244,(((t112:0.7170655408,(t113:0.976886861,t114:0.9406369971):0.7471234254):0.8065501407,((t115:0.1713845057,(t116:0.7861330248,t117:0.6082276558):0.8413775554):0.3245444677,t118:0.3892389825):0.5992471091):0.7592411407,(((t119:0.535931844,t120:0.09058958571):0.4227561057,(t121:0.5531579193,t122:0.8276180199):0.6653355309):0.0941624688,t123:0.3623022255):0.1494971744):0.3526274569):0.9720881658):0.8149677955):0.6065687414,((((((t124:0.5406888947,t125:0.8892341822):0.06211395678,((t126:0.8203180477,(t127:0.8536844573,t128:0.360511546):0.9030223228):0.9095590916,((t129:0.9110714826,(t130:0.2346256471,t131:0.6523390864):0.1288849309):0.7077432328,(t132:0.4060195235,t133:0.1661393729):0.3910941551):0.205704404):0.8609933471):0.3724007562,((t134:0.1731842053,(t135:0.7232482471,(t136:0.3883952193,((t137:0.6709475764,t138:0.0372075201):0.5473196667,(t139:0.8092764446,t140:0.4123262055):0.2000603897):0.55258787):0.2654263263):0.745555162):0.2956101163,((t141:0.52147611,(t142:0.9462005703,t143:0.5671354234):0.6887917654):0.362258781,t144:0.4798202242):0.8242726682):0.6072624433):0.695287361,((((t145:0.03793937969,t146:0.07275558705):0.3482963489,t147:0.1457363514):0.1479936559,(t148:0.7158309214,((t149:0.2174433649,t150:0.04072828358):0.4112026501,t151:0.6422409331):0.3413406226):0.1693999742):0.6631712937,(((t152:0.2706006162,t153:0.9267972289):0.1387761638,((((t154:0.2563392594,t155:0.3058371837):0.5946117372,t156:0.6161190302):0.6970871226,(t157:0.2388902532,(t158:0.9486316761,t159:0.215360787):0.168830334):0.03888285463):0.1640696453,t160:0.6803096831):0.1418975852):0.4218000816,(((t161:0.8702562298,t162:0.9289729816):0.05807372741,t163:0.3533785399):0.5012762842,(((t164:0.8666574673,t165:0.9603798252):0.7887994377,t166:0.857058729):0.4139410679,(t167:0.5900272813,t168:0.3345388798):0.06017537019):0.9609203783):0.7103463742):0.696603697):0.6451920038):0.1909481271,((((t169:0.9171597108,t170:0.9479122513):0.7170342554,(t171:0.2722596873,((t172:0.1194724559,(t173:0.03922236571,t174:0.6290624789):0.07739861775):0.8598598302,(t175:0.2009421999,(t176:0.06154947914,t177:8.997193072E-4):0.04738179315):0.3235510678):0.3443877005):0.6351028818):0.5525081949,((((t178:0.7599076207,t179:0.2997759853):0.5921433992,t180:0.7098581635):0.3725496214,(t181:0.5053773888,(t182:0.5991492711,(t183:0.5036820578,t184:0.6361607853):0.510631816):0.9604382808):0.2464167587):0.6073093358,(((t185:0.03128415369,(t186:0.5260852403,(t187:0.878767435,t188:0.4992109234):0.5333148066):0.00347468094):0.5590308013,t189:0.3710992143):0.5034162949,(t190:0.778916508,((t191:0.3069154553,(((t192:0.9946115273,t193:0.9138687006):0.5209144899,t194:0.5152770842):0.9462409306,t195:0.7395236609):0.4110851623):0.930918345,(((t196:0.7895439987,((t197:0.4697002599,t198:0.1383787312):0.6911794308,(t199:0.8664436699,t200:0.1959039853):0.8656513852):0.3620497067):0.2839249384,(t201:0.6558795469,t202:0.2103423763):0.969477433):0.9058840063,(t203:0.0856692954,t204:0.4175976661):0.820434629):0.5355881769):0.2263581599):0.4512835185):0.7323478526):0.2479199937):0.1964542414,((t205:0.7537573762,(t206:0.1392466244,(t207:0.5136175761,(t208:0.7852529553,t209:0.07355738804):0.1220811389):0.7572090242):0.1422528555):0.5948274662,(((((t210:0.3068353184,(t211:0.3314456891,((t212:0.5265486804,t213:0.1382007354):0.1814086549,t214:0.9276472756):0.07718444197):0.03486835537):0.1617580003,(t215:0.3328830956,t216:0.8558843595):0.8366736979):0.347376487,t217:0.8222538356):0.2337225529,(t218:0.06199815008,t219:0.45975962):0.179990889):0.0635867205,(t220:0.3214025751,(t221:0.5022090652,t222:0.6454557138):0.6956466341):0.2711792416):0.1847200533):0.1051658324):0.4945860899):0.936143348,(((t223:0.06268779701,((t224:0.3337278806,t225:0.1570303424):0.3089733059,(t226:0.5069784883,t227:0.1434204187):0.2001587199):0.04750720505):0.3600859912,((((t228:0.9994731578,(t229:0.8934116936,t230:0.03698333143):0.8173468311):0.3089058488,((((t231:0.3216121283,t232:0.5232846253):0.8687884973,(t233:0.6280638413,((t234:0.6543256822,t235:0.8677638234):0.8895299246,t236:0.4047793006):0.7147388768):0.3533478715):0.9470084386,t237:0.7769409856):0.4955915695,((t238:0.2772087415,(t239:0.4904922615,(t240:0.05356206303,t241:0.08998329984):0.8154862223):0.5610961432):0.1617916438,(t242:0.5707751412,(t243:0.9836868793,t244:0.1984052949):0.6953297216):0.05552111682):0.9476150468):0.2473166997):0.9623488116,((t245:0.7935025664,t246:0.08509867964):0.3953444003,(t247:0.09163277131,(t248:0.5201428954,t249:0.8055520628):0.7452739514):0.3989078877):0.07581191277):0.9779064963,(((t250:0.943611098,(t251:0.33392801,t252:0.5996331484):0.4291575127):0.4906436009,((((t253:0.7749450852,(t254:0.8616885878,t255:0.585028409):0.06060880423):0.1238881133,((t256:0.7451687793,t257:0.6925335305):0.05338745634,t258:0.3357626374):0.2069296469):0.09644073155,((((t259:0.2258843291,t260:0.2671526412):0.3940743534,(t261:0.5022506947,(t262:0.9498897423,t263:0.1406114365):0.2847759123):0.04320593993):0.6982026948,t264:0.2693712024):0.959781138,(((t265:0.6035173486,t266:0.5529949202):0.9900399651,(t267:0.5455351078,t268:0.3530619899):0.4626278321):0.2735997427,(t269:0.9580646451,(t270:0.3280033092,t271:0.7206294278):0.03739526332):0.4967516926):0.9350089293):0.4371789068):0.1014483059,t272:0.2867298371):0.07522285799):0.06352435821,((t273:0.4001782183,t274:0.7190070178):0.1696753846,(t275:0.5535608665,t276:0.01324651297):0.2691543309):0.8676247413):0.8461736294):0.1769516913):0.344365149,(((t277:0.3245107541,(t278:0.4142541443,t279:0.5857141651):0.819547887):0.0867733527,(t280:0.4938162852,(t281:0.2444119717,t282:0.08141433029):0.05381231918):0.8375963389):0.176160393,((t283:0.4199601968,t284:0.8354801824):0.3150380594,(((t285:0.9818797186,(t286:0.8971825438,((t287:0.5155417006,t288:0.8260786769):0.7060374152,t289:0.6001661876):0.4120474763):0.9949228324):0.8038698458,t290:0.1939124272):0.6380942846,t291:0.3665255161):0.459349304):0.482901911):0.4833473735):0.5903116504):0.9973697898)";
		int nrTaxa = 291;
		int[] states = new int[291];
		Arrays.fill(states, 1);

		orig = ".1";
		locations = "" ;
		for (int i=1; i<=states.length; i++){
			locations = locations + "t" + i + "=" + (states[i-1]-1) + (i<states.length?",":"");
		}

		stateNumber = "2";
		migrationMatrix = "0.2 0.1";
		frequencies = "0.5 0.5";

		R0 = Double.toString(4./3.) + " " + Double.toString(5.) ;
		becomeUninfectiousRate = "1.5 1.25";
		samplingProportion = Double.toString(1./3.) + " " + Double.toString(1./2.);

		conditionOnSurvival = true;

		logL = bdm_likelihood(stateNumber,
				migrationMatrix,
				frequencies,
				treeBig, orig,
				R0,  null,
				becomeUninfectiousRate,
				samplingProportion, null,
				"t", locations, nrTaxa, null,
				conditionOnSurvival);

		//		System.out.println("Log-likelihood ("+(conditionOnSurvival?"":"not ")+"conditioned on survival) " + logL + "\t");
		assertEquals(-661.9588648301033, logL, 1e-5); // result from BEAST, not checked in R

	}

	/**
	 * Test on migration combined with infection among demes
	 * Uncoloured trees
	 * No rate changes
	 * Symmetric and assymetric configurations
	 * reference from R
	 * @throws Exception
	 */
	public void testLikelihoodCalculationInfAmongDemes() throws Exception{

		// uncoloured, symmetric tree

		BirthDeathMigrationDistribution bdm =  new BirthDeathMigrationDistribution();

		String newick = "((t3[&type=1]:0.004214277605,t4[&type=1]:0.02157681391):0.229186993,(t2[&type=0]:0.624713651,t1[&type=1]:1.347400211):0.06231047755);";
		String prefixname = "t";

		int nrTaxa = 4 ;

		String orig = "0.02686563367";
		String stateNumber = "2";
		String frequencies = "0.5 0.5";

		String birth = "2. 2.";
		String birthRateAmongDemes = "1. 1.";
		String deathRate = "0.5 0.5";
		String samplingRate = "0.5 0.5";

		boolean conditionOnSurvival = true;

		Tree tree = new TreeParser();
		tree.setInputValue("adjustTipHeights", "false");
		tree.setInputValue("IsLabelledNewick", "true");
		tree.setInputValue("newick", newick);
		tree.initAndValidate();

		bdm.setInputValue("tree", tree);
		bdm.setInputValue("typeLabel", "type");

		bdm.setInputValue("origin", Double.toString(Double.parseDouble(orig)+tree.getRoot().getHeight()));
		bdm.setInputValue("stateNumber", stateNumber);
		bdm.setInputValue("migrationMatrix", "0 0");
		bdm.setInputValue("frequencies", frequencies);

		bdm.setInputValue("birthRate", birth);
		bdm.setInputValue("birthRateAmongDemes", birthRateAmongDemes);
		bdm.setInputValue("deathRate", deathRate);
		bdm.setInputValue("samplingRate", samplingRate);

		bdm.setInputValue("conditionOnSurvival", conditionOnSurvival);

		bdm.initAndValidate();

		double logL = bdm.calculateLogP();

		//System.out.println("Log-likelihood " + logL + " - testLikelihoodCalculationInfAmongDemes \t");
		assertEquals(-5.1966118470881, logL, 1e-3);

		// uncoloured, assymetric tree

		newick = "((3:1.5,4:0.5):1,(1:2,2:1):3);";
		prefixname = "";

		conditionOnSurvival = true;

		nrTaxa = 4 ;
		int[] states = new int[]{1,2,2,1};

		orig = "1.";

		String locations = "" ;
		for (int i=1; i<=states.length; i++){
			locations = locations + prefixname + i + "=" + (states[i-1]-1) + (i<states.length?",":"");
		}

		stateNumber = "2";
		String migrationMatrix = "0.2 0.1" ; //"0.2 0.1" ; //"1. 1.";
		frequencies = "0.5 0.5";

		bdm =  new BirthDeathMigrationDistribution();

		birth = "2. 6.25";
		birthRateAmongDemes = "0.2 0.1"; // "0. 0.";
		deathRate = "1.2 0.625";
		samplingRate = "0.3 0.625";

		bdm.setInputValue("birthRate", birth);
		bdm.setInputValue("birthRateAmongDemes", birthRateAmongDemes);
		bdm.setInputValue("migrationMatrix", "0 0");

		bdm.setInputValue("deathRate", deathRate);
		bdm.setInputValue("samplingRate", samplingRate);

		tree = new TreeParser();
		tree.setInputValue("adjustTipHeights", "false");
		tree.setInputValue("IsLabelledNewick", "true");
		tree.setInputValue("newick", newick);
		tree.initAndValidate();

		ArrayList<Taxon> taxa = new ArrayList<Taxon>();

		for (int i=1; i<=nrTaxa; i++){
			taxa.add(new Taxon(prefixname+i));
		}

		TraitSet trait = new TraitSet();
		trait.setInputValue("taxa", new TaxonSet(taxa));
		trait.setInputValue("value", locations);
		trait.setInputValue("traitname", "tiptypes");
		trait.initAndValidate();

		bdm.setInputValue("tree", tree);
		bdm.setInputValue("tiptypes", trait);

		bdm.setInputValue("origin", Double.toString(Double.parseDouble(orig)+tree.getRoot().getHeight()));
		bdm.setInputValue("stateNumber", stateNumber);
		bdm.setInputValue("frequencies", frequencies);

		bdm.setInputValue("conditionOnSurvival", conditionOnSurvival);

		bdm.initAndValidate();

		logL = bdm.calculateLogP(); 

		// System.out.println("Log-likelihood" + logL + "\t");
		assertEquals(-26.7939, logL, 1e-5);  //result from R
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

		String tree ="((3:1.5,4:0.5):1,(1:1,2:1):3);";
		String orig=".1";
		String stateNumber = "2";
		String locations = "1=1,2=0,3=0,4=1" ;
		String migrationMatrix = "0.1 0.2 0.15 0.25";
		String frequencies = "0.5 0.5";

		String R0 = "6 2 5 2.5";
		String R0AmongDemes = "1.1 1.2 1.3 1.15 ";
		String becomeUninfectiousRate = "0.5 0.45 0.55 0.6";
		String samplingProportion = "0.5 0.333333 0.45 0.35";
		
		String intervalTimes = "0. 1.";


		boolean conditionOnSurvival = false;

		double logL = bdm_likelihood(stateNumber,
				migrationMatrix,
				frequencies,
				tree, orig,
				R0,R0AmongDemes,
				becomeUninfectiousRate,
				samplingProportion,
				null,
				"", locations, 4,
				intervalTimes, conditionOnSurvival);


		// System.out.println("Log-likelihood = " + logL);
		assertEquals(-16.466832439520886, logL, 1e-4); // result from BDMM, 28/06/2017

		
		//Test without rate change, with asymmetric types
		tree ="((3:1.5,4:0.5):1,(1:2,2:1):3);"; //
		orig="1."; //
		stateNumber = "2";
		migrationMatrix = "0.5 0.";
		frequencies = "1 0";

		// test without rate change
		R0AmongDemes = "0. 2.";
		R0 = "0. 0.";
		becomeUninfectiousRate = "0. 0.75";
		samplingProportion = "0. 0.7";
		locations = "1=1,2=1,3=1,4=1" ;

		conditionOnSurvival = false;

		logL = bdm_likelihood(stateNumber,
				migrationMatrix,
				frequencies,
				tree, orig,
				R0, R0AmongDemes,
				becomeUninfectiousRate,
				samplingProportion, null,
				"", locations, 4, null, conditionOnSurvival);


	//	System.out.println("Log-likelihood = " + logL);
		assertEquals(-12.1441, logL, 1e-4); // tanja's result from R
		
	}
	
	/**
	 * Test of migration with 3 types
	 * Uncoloured tree
	 * No rate change, no SA
	 * Reference from BDMM itself 
	 * @throws Exception
	 */
	public void testMig3types() throws Exception {

		String tree ="((3:1.5,4:0.5):1,(1:1,2:1):3);";
		String orig=".1";
		String stateNumber = "3";
		String locations = "1=1,2=0,3=2,4=1" ;
		String migrationMatrix = "0.1 0.2 0.15 0.12 0.12 0.15";
		String frequencies = Double.toString(1./3.) + " " + Double.toString(1./3.) + " " + Double.toString(1./3.) ;

		String R0 = "6 2 5";

		String becomeUninfectiousRate = "0.5 0.45 0.55";
		String samplingProportion = "0.5 0.333333 0.45";
	
		boolean conditionOnSurvival = false;

		double logL = bdm_likelihood(stateNumber,
				migrationMatrix,
				frequencies,
				tree, orig,
				R0,null,
				becomeUninfectiousRate,
				samplingProportion,
				null,
				"", locations, 4,
				null, conditionOnSurvival);


		 //System.out.println("Log-likelihood = " + logL);
		assertEquals(-16.88601100061662, logL, 1e-4); // result from BDMM, version 0.2.0, 06/07/2017
		
		// coloured tree
		String treeMT ="((3[&state=2]:0.5)[&state=2]:1,4[&state=1]:0.5)[&state=1]:1,(1[&state=1]:1,(2[&state=0]:0.5)[&state=1]:0.5)[&state=1]:3)[&state=1]:0.0;";
//		String orig=".1";
//		String stateNumber = "3";
//		String locations = "1=1,2=0,3=2,4=1" ;
//		String migrationMatrix = "0.1 0.2 0.15 0.12 0.12 0.15";
//		String frequencies = Double.toString(1./3.) + " " + Double.toString(1./3.) + " " + Double.toString(1./3.) ;
//
//		String R0 = "6 2 5";
//
//		String becomeUninfectiousRate = "0.5 0.45 0.55";
//		String samplingProportion = "0.5 0.333333 0.45";
//	
//		boolean conditionOnSurvival = false;
//
//		double logL = bdm_likelihood(stateNumber,
//				migrationMatrix,
//				frequencies,
//				tree, orig,
//				R0,null,
//				becomeUninfectiousRate,
//				samplingProportion,
//				null,
//				"", locations, 4,
//				null, conditionOnSurvival);
//
//
//		 //System.out.println("Log-likelihood = " + logL);
//		assertEquals(-16.88601100061662, logL, 1e-4); // result from BDMM, version 0.2.0, 06/07/2017

	}

	/**
	 * Test uncoloured tree with unknown states
	 * With and without migration, with and without rate changes
	 * 2 types, assymetric
	 * @throws Exception
	 */
	@Test
	public void testUnknownStates() throws Exception {

		String tree ="((3:1.5,4:0.5):1,(1:1,2:1):3);"; //
		String orig=".1"; //
		String stateNumber = "2";
		String migrationMatrix = "0. 0.";
		String frequencies = "0.5 0.5";

		// test without rate change
		String R0 = "6 2";
		String R0AmongDemes = "1. 1.";
		String becomeUninfectiousRate = "0.5 1";
		String samplingProportion = "0.5 0.333333";
		String locations = "1=-1,2=-1,3=-1,4=-1" ;

		boolean conditionOnSurvival = false;

		double logL;

		logL = bdm_likelihood(stateNumber,
				migrationMatrix,
				frequencies,
				tree, orig,
				R0,R0AmongDemes,
				becomeUninfectiousRate,
				samplingProportion,
				null,
				"", locations, 4, null,
				conditionOnSurvival);
		
		//	System.out.println("Log-likelihood = " + logL);
		assertEquals(-18.82798, logL, 1e-4); // tanja's result from R
		
		// with migration
		migrationMatrix = "0.3 0.4";
		
		logL = bdm_likelihood(stateNumber,
				migrationMatrix,
				frequencies,
				tree, orig,
				R0,R0AmongDemes,
				becomeUninfectiousRate,
				samplingProportion,
				null,
				"", locations, 4, null,
				conditionOnSurvival);

		// System.out.println("Log-likelihood = " + logL + " - Unknow states with migration");
		assertEquals(-18.986212857895506, logL, 1e-4); // reference from BDMM - 0.2.0 - 06/07/2017
		
		// with migration and rate changes
		migrationMatrix = "0.3 0.4 0.35 0.32";
		R0 = "6 2 5 1.5";
		R0AmongDemes = "1. 1. 1.2 1.2";
		becomeUninfectiousRate = "0.5 1. 1. 0.5";
		samplingProportion = "0.5 0.333333 0.45 0.4";
		String intervalTimes = "0. 1.5";
		
		logL = bdm_likelihood(stateNumber,
				migrationMatrix,
				frequencies,
				tree, orig,
				R0,R0AmongDemes,
				becomeUninfectiousRate,
				samplingProportion,
				null,
				"", locations, 4,
				intervalTimes,
				conditionOnSurvival);
		
		System.out.println("Log-likelihood = " + logL + " - Unknow states with migration and rate changes");
		assertEquals(-17.87099909579358, logL, 1e-4); // reference from BDMM - 0.2.0 - 06/07/2017
	}
	
	/**
	 * Test on combining migration with rho-sampling
	 * Uncoloured and coloured trees
	 * Reference from BDMM
	 * @throws Exception
	 */
	@Test
	public void testLikelihoodMigrationRhoSampling() throws Exception{
		// Uncoloured tree
		Tree tree = new TreeParser("((1[&type=0]: 4.5, 2[&type=1]: 4.5):1,3[&type=0]:5.5);",false);

		BirthDeathMigrationDistribution bdm =  new BirthDeathMigrationDistribution();

		bdm.setInputValue("tree", tree);

		bdm.setInputValue("typeLabel", "type");
		bdm.setInputValue("stateNumber", "2");
		bdm.setInputValue("migrationMatrix", "0.3 0.4");

		bdm.setInputValue("parallelize", false);
		bdm.setInputValue("migrationMatrixScaleFactor", "1");


		bdm.setInputValue("frequencies", "0.6 0.4");

		bdm.setInputValue("R0", new RealParameter("1.5 1.4"));
		bdm.setInputValue("becomeUninfectiousRate", new RealParameter("1.5 1.3"));
		bdm.setInputValue("samplingProportion", new RealParameter("0. 0.") );

		bdm.setInputValue("rho", new RealParameter("0.01 0.015") );
		bdm.setInputValue("conditionOnSurvival", false);

		bdm.initAndValidate();
		
//		System.out.println("Likelihood: "+  bdm.calculateLogP() + " rho-sampling with migration");

		assertEquals(-8.906223150087108, bdm.calculateLogP(), 1e-4);   // Reference from BDMM - version 0.2.0 - 06/07/2017

	}
	
	/**
	 * Basic test on sampled-ancestors lik. calculation.
	 * 2 leaves, 1 SA. 1 type, no rho-sampling, no rate-change
	 * Coloured and uncoloured trees
	 * Reference value from BDSKY (23/03/2017)
	 * @throws Exception
	 */
	@Test
	public void testSALikelihoodMini() throws Exception {

		// uncoloured tree
		Tree tree = new TreeParser("((3[&type=0]: 1.5, 6[&type=0]: 0)5[&type=0]: 3.5, 4[&type=0]: 4) ;",false);

		BirthDeathMigrationDistribution bdm =  new BirthDeathMigrationDistribution();

		bdm.setInputValue("tree", tree);
		bdm.setInputValue("typeLabel", "type");

		bdm.setInputValue("stateNumber", "1");
		bdm.setInputValue("migrationMatrix", "0.");
		bdm.setInputValue("frequencies", "1");

		bdm.setInputValue("R0", new RealParameter("1.5"));
		bdm.setInputValue("becomeUninfectiousRate", new RealParameter("1.5"));
		bdm.setInputValue("samplingProportion", new RealParameter("0.2") );
		bdm.setInputValue("removalProbability", new RealParameter("0.9") );

		bdm.setInputValue("conditionOnSurvival", true);
		bdm.setInputValue("origin", "6.");

		bdm.initAndValidate();

		double logL = bdm.calculateLogP();
		assertEquals(-18.854438107814335, logL, 1e-4); //Reference value from BDSKY (23/03/2017)
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

		// uncoloured tree
		Tree tree = new TreeParser("((3[&type=0]: 1.5, 6[&type=0]: 0)5[&type=0]: 3.5, 4[&type=0]: 4) ;",false);

		BirthDeathMigrationDistribution bdm =  new BirthDeathMigrationDistribution();

		bdm.setInputValue("tree", tree);
		bdm.setInputValue("typeLabel", "type");

		bdm.setInputValue("stateNumber", "1");
		bdm.setInputValue("migrationMatrix", "0.");
		bdm.setInputValue("frequencies", "1");

		bdm.setInputValue("R0", new RealParameter("1.5"));
		bdm.setInputValue("becomeUninfectiousRate", new RealParameter("1.5"));
		bdm.setInputValue("samplingProportion", new RealParameter("0.2") );
		bdm.setInputValue("removalProbability", new RealParameter("0.9") );

		bdm.setInputValue("conditionOnSurvival", true);
		bdm.setInputValue("origin", "6.");

		bdm.setInputValue("rho", new RealParameter("0.3 0.05"));
		bdm.setInputValue("rhoSamplingTimes", new RealParameter("0 1.5") );

		bdm.setInputValue("reverseTimeArrays", "false false false true");

		bdm.initAndValidate();


		//		double a = bdm.calculateLogP();	
		//		System.out.println(a);

		assertEquals(-22.348462265673483, bdm.calculateLogP(), 1e-5); //Reference value from BDSKY (06/04/2017)

	}

	/**
	 * Test on sampled-ancestors lik. calculation with no sampled ancestor
	 * Coloured and Uncoloured
	 * No rate-change, one state, 4 tips
	 * This state is just there in case something is broken with sampled ancestors,
	 * helps for debugging if combined with testSALikelihoodMini for instance
	 * @throws Exception
	 */
	@Test
	public void testSALikelihoodCalculationWithoutAncestors() throws Exception {

		// uncoloured tree
		BirthDeathMigrationDistribution bdm =  new BirthDeathMigrationDistribution();

		ArrayList<Taxon> taxa = new ArrayList<Taxon>();

		for (int i=1; i<=4; i++){
			taxa.add(new Taxon(""+i));
		}

		Tree tree = new TreeParser();
		tree.setInputValue("taxonset", new TaxonSet(taxa));
		tree.setInputValue("adjustTipHeights", "false");
		tree.setInputValue("IsLabelledNewick", "true");
		tree.setInputValue("newick", "((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);");
		tree.initAndValidate();

		TraitSet trait = new TraitSet();
		trait.setInputValue("taxa", new TaxonSet(taxa));
		trait.setInputValue("value", "1=0,2=0,3=0,4=0");
		trait.setInputValue("traitname", "tiptypes");
		trait.initAndValidate();

		bdm.setInputValue("tree", tree);
		bdm.setInputValue("tiptypes", trait);


		bdm.setInputValue("stateNumber", "1");
		bdm.setInputValue("migrationMatrix", "0.");
		bdm.setInputValue("frequencies", "1");

		bdm.setInputValue("R0", new RealParameter("1.5"));
		bdm.setInputValue("becomeUninfectiousRate", new RealParameter("1.5"));
		bdm.setInputValue("samplingProportion", new RealParameter("0.3") );
		bdm.setInputValue("removalProbability", new RealParameter("0.9") );
		bdm.setInputValue("conditionOnSurvival", true);


		// without a defined origin
		bdm.initAndValidate();

		// likelihood conditioning on at least one sampled individual    - "true" result from BEAST one-deme SA model 09 June 2015 (DK)
		assertEquals(-15.99699690815937, bdm.calculateLogP(), 1e-4);


		// with an origin
		bdm.setInputValue("origin", "10.");
		bdm.initAndValidate();


		// likelihood conditioning on at least one sampled individual    - "true" result from BEAST one-deme SA model 09 June 2015 (DK)
		assertEquals(-25.991511346557598, bdm.calculateLogP(), 1e-4);

	}

	public double bdm_likelihood(String statenumber, String migrationMatrix,
			String frequencies, Tree tree, String typeLabel, String origin,
			String R0, String R0AmongDemes, String becomeUninfectiousRate, String samplingProportion, String removalProbability,
			String intervalTimes, Boolean conditionOnSurvival) throws Exception {

		BirthDeathMigrationDistribution bdm =  new BirthDeathMigrationDistribution();

		bdm.setInputValue("tree", tree);
		bdm.setInputValue("typeLabel", typeLabel);

		bdm.setInputValue("origin", Double.toString(Double.parseDouble(origin)+tree.getRoot().getHeight()));
		bdm.setInputValue("stateNumber", statenumber);
		bdm.setInputValue("migrationMatrix", migrationMatrix);
		bdm.setInputValue("frequencies", frequencies);
		bdm.setInputValue("checkRho", false);

		bdm.setInputValue("R0", R0);

		if (R0AmongDemes != null) bdm.setInputValue("R0AmongDemes", R0AmongDemes);

		if (removalProbability != null) bdm.setInputValue("removalProbability", removalProbability);

		bdm.setInputValue("becomeUninfectiousRate", becomeUninfectiousRate);
		bdm.setInputValue("samplingProportion", samplingProportion);
		bdm.setInputValue("intervalTimes", intervalTimes);

		bdm.setInputValue("conditionOnSurvival", conditionOnSurvival);

		bdm.initAndValidate();

		long startTime = System.currentTimeMillis();
		double logL = bdm.calculateLogP();
		runtime = System.currentTimeMillis() - startTime;

		return logL;

	}

	public double bdm_likelihood(String statenumber, String migrationMatrix,
			String frequencies, Tree tree, TraitSet trait, String origin,
			String R0, String R0AmongDemes, String becomeUninfectiousRate, String samplingProportion, String removalProbability,
			String intervalTimes, Boolean conditionOnSurvival) throws Exception {

		BirthDeathMigrationDistribution bdm =  new BirthDeathMigrationDistribution();

		bdm.setInputValue("tree", tree);
		bdm.setInputValue("tiptypes", trait);

		bdm.setInputValue("origin", Double.toString(Double.parseDouble(origin)+tree.getRoot().getHeight()));
		bdm.setInputValue("stateNumber", statenumber);
		bdm.setInputValue("migrationMatrix", migrationMatrix);
		bdm.setInputValue("frequencies", frequencies);
		bdm.setInputValue("checkRho", false);

		bdm.setInputValue("R0", R0);

		if (R0AmongDemes != null) bdm.setInputValue("R0AmongDemes", R0AmongDemes);

		if (removalProbability != null) bdm.setInputValue("removalProbability", removalProbability);

		bdm.setInputValue("becomeUninfectiousRate", becomeUninfectiousRate);
		bdm.setInputValue("samplingProportion", samplingProportion);
		bdm.setInputValue("intervalTimes", intervalTimes);

		bdm.setInputValue("conditionOnSurvival", conditionOnSurvival);

		bdm.initAndValidate();

		long startTime = System.currentTimeMillis();
		double logL = bdm.calculateLogP();
		runtime = System.currentTimeMillis() - startTime;

		return logL;
	}

	public double bdm_likelihood(String statenumber, String migrationMatrix,
			String frequencies, String newick, String origin,
			String R0, String R0AmongDemes, String becomeUninfectiousRate, String samplingProportion, String removalProbability,
			String prefixname, String locations, int nrTaxa, String intervalTimes, Boolean conditionOnSurvival) throws Exception {

		ArrayList<Taxon> taxa = new ArrayList<Taxon>();

		for (int i=1; i<=nrTaxa; i++){
			taxa.add(new Taxon(prefixname+i));
		}

		Tree tree = new TreeParser();

		tree.setInputValue("taxonset", new TaxonSet(taxa));
		tree.setInputValue("adjustTipHeights", "false");
		tree.setInputValue("IsLabelledNewick", "true");
		tree.setInputValue("newick", newick);
		tree.initAndValidate();

		TraitSet trait = new TraitSet();
		trait.setInputValue("taxa", new TaxonSet(taxa));
		trait.setInputValue("value", locations);
		trait.setInputValue("traitname", "tiptypes");
		trait.initAndValidate();

		return   bdm_likelihood(statenumber, migrationMatrix, frequencies, tree, trait, origin,
				R0, R0AmongDemes, becomeUninfectiousRate, samplingProportion, removalProbability,
				intervalTimes, conditionOnSurvival);
	}

}