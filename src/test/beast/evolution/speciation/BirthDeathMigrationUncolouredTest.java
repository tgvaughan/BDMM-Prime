package test.beast.evolution.speciation;

import beast.core.parameter.RealParameter;
import beast.evolution.speciation.BirthDeathMigrationModelUncoloured;
import beast.evolution.speciation.PiecewiseBirthDeathMigrationDistribution;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TraitSet;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.alignment.Taxon;
import beast.util.TreeParser;
//import beast.util.ZeroBranchSATreeParser;
import junit.framework.TestCase;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * User: Denise
 * Date: Jul 5, 2013
 * Time: 2:11:07 PM
 */
public class BirthDeathMigrationUncolouredTest extends TestCase {


	double runtime;
	int maxEvalsUsed = -1;

	Boolean conditionOnSurvival = false;

	// Initial version of this test from 20th December commit by Denise 
	//    @Test
	//    public void testmultiRho2tipsNoDeath() throws Exception {
	//
	//        Tree tree = new TreeParser("(3[&type=0]: 4, 4[&type=0]: 4) ;", false);
	//
	//        BirthDeathMigrationModelUncoloured bdm = new BirthDeathMigrationModelUncoloured();
	//
	//        bdm.setInputValue("tree", tree);
	//        bdm.setInputValue("typeLabel", "type");
	//        bdm.setInputValue("origin", "5.");
	//
	//        bdm.setInputValue("stateNumber", "1");
	//        bdm.setInputValue("migrationMatrix", "0.");
	//        bdm.setInputValue("frequencies", "1");
	//
	//        bdm.setInputValue("birthRate", new RealParameter("2.25"));
	//        bdm.setInputValue("deathRate", new RealParameter("0."));
	//        bdm.setInputValue("samplingRate", new RealParameter("0."));
	//        bdm.setInputValue("rhoSamplingTimes", new RealParameter("0. 2.5"));
	//        bdm.setInputValue("reverseTimeArrays", "false false false true");
	//
	//        bdm.setInputValue("rho", new RealParameter("0.2 1."));
	//
	////        bdm.setInputValue("conditionOnSurvival", false);
	////        bdm.initAndValidate();
	////        assertEquals(-19.88535688641209, bdm.calculateLogP(), 1e-4);   // this result is from BEAST, not double checked in R
	////
	////        bdm.setInputValue("rho", new RealParameter("0.2 0.6"));
	////        bdm.initAndValidate();
	////        assertEquals(-18.865767180278915, bdm.calculateLogP(), 1e-4);   // this result is from BEAST, not double checked in R
	//
	//        tree = new TreeParser("(3[&type=0]: 1.5, 4[&type=0]: 4) ;",false);
	//        bdm.setInputValue("tree", tree);
	//
	//
	//        //        bdm.setInputValue("rho", new RealParameter("0.2 1."));
	//        //        bdm.initAndValidate();
	//        //        assertEquals(-15.646651247531981, bdm.calculateLogP(), 1e-4);   // this result is from BEAST, not double checked in R
	//
	//        //        bdm.setInputValue("conditionOnSurvival", false);
	//        //        bdm.setInputValue("rho", new RealParameter("0.2 0.6"));
	//        //        bdm.initAndValidate();
	//        //        assertEquals(-15.133091119955177, bdm.calculateLogP(), 1e-4);   // this result is from BEAST, not double checked in R
	//
	//        bdm.setInputValue("conditionOnSurvival", false);
	//        bdm.setInputValue("deathRate", "1.5");
	//        bdm.setInputValue("rho", new RealParameter("0.2 0.6"));
	//        bdm.initAndValidate();
	//        assertEquals(-8.637410990319223, bdm.calculateLogP(), 1e-4);   // this result is from BEAST (BDSKY), not double checked in R
	//
	//    }

	// Version of this test from Denise's email of the 20th December
	@Test
	public void testmultiRho2tipsNoDeath() throws Exception {

		Tree tree = new TreeParser("(3[&type=0]: 4, 4[&type=0]: 4) ;", false);

		BirthDeathMigrationModelUncoloured bdm = new BirthDeathMigrationModelUncoloured();

		bdm.setInputValue("tree", tree);
		bdm.setInputValue("typeLabel", "type");
		bdm.setInputValue("origin", "5.");

		bdm.setInputValue("stateNumber", "1");
		bdm.setInputValue("migrationMatrix", "0.");
		bdm.setInputValue("frequencies", "1");

		bdm.setInputValue("birthRate", new RealParameter("2.25"));
		bdm.setInputValue("deathRate", "1.5");

		bdm.setInputValue("samplingRate", new RealParameter("0."));

		bdm.setInputValue("rhoSamplingTimes", new RealParameter("0. 2.5"));
		bdm.setInputValue("reverseTimeArrays", "false false false true");


		tree = new TreeParser("(3[&type=0]: 1.5, 4[&type=0]: 4) ;",false);
		bdm.setInputValue("tree", tree);

		bdm.setInputValue("conditionOnSurvival", false);
		bdm.setInputValue("rho", new RealParameter("0.2 0.6"));
		bdm.initAndValidate();
		assertEquals(-8.637410990319223, bdm.calculateLogP(), 1e-4);   // this result is from BEAST (BDSKY), not double checked in R

	}

	@Test
	public void testmultiRho2tips() throws Exception {

		Tree tree = new TreeParser("(3[&type=0]: 4, 4[&type=0]: 4) ;",false);

		BirthDeathMigrationModelUncoloured bdm =  new BirthDeathMigrationModelUncoloured();

		bdm.setInputValue("tree", tree);
		bdm.setInputValue("typeLabel", "type");

		bdm.setInputValue("stateNumber", "1");
		bdm.setInputValue("migrationMatrix", "0.");
		bdm.setInputValue("frequencies", "1");

		bdm.setInputValue("R0", new RealParameter("1.5"));
		bdm.setInputValue("becomeUninfectiousRate", new RealParameter("1.5"));
		bdm.setInputValue("samplingProportion", new RealParameter("0.") );
		// TO DO REMOVE COMMENTS
		// changed the rho to get sth non symetrical
		bdm.setInputValue("rho", new RealParameter("0.2 1.") );
		bdm.setInputValue("rhoSamplingTimes", new RealParameter("0. 2.5") );
		bdm.setInputValue("reverseTimeArrays", "false false false true");
		bdm.setInputValue("conditionOnSurvival", true);
		bdm.setInputValue("origin", "5.");

		bdm.initAndValidate();
		// TO DO REMOVE: Jeremie changed the value here (comes from BDSKY)
		assertEquals(-10.569863754307026, bdm.calculateLogP(), 1e-4);   // this result is from BEAST, not double checked in R

		tree = new TreeParser("(3[&type=0]: 1.5, 4[&type=0]: 4) ;",false);
		bdm.setInputValue("tree", tree);
		bdm.initAndValidate();
		// TO DO REMOVE: Jeremie changed the value here (comes from BDSKY)
		assertEquals(-8.099631076932816, bdm.calculateLogP(), 1e-4);   // this result is from BEAST, not double checked in R

	}

	@Test
	public void testmultiRho() throws Exception {

		TreeParser tree = new TreeParser("((3[&type=0]: 1.5, 4[&type=0]: 1.5)[&type=0]: 1 , (1[&type=0]: 2, 2[&type=0]: 2)[&type=0]: 3)[&type=0];", false);

		BirthDeathMigrationModelUncoloured bdm =  new BirthDeathMigrationModelUncoloured();

		bdm.setInputValue("tree", tree);
		bdm.setInputValue("typeLabel", "type");

		bdm.setInputValue("stateNumber", "1");
		bdm.setInputValue("migrationMatrix", "0.");
		bdm.setInputValue("frequencies", "1");

		bdm.setInputValue("origin", "10.");
		bdm.setInputValue("R0", new RealParameter("1.5"));
		bdm.setInputValue("becomeUninfectiousRate", new RealParameter("1.5"));
		bdm.setInputValue("samplingProportion", new RealParameter("0.") );
		bdm.setInputValue("rho", new RealParameter("1. 1.") );
		bdm.setInputValue("rhoSamplingTimes", new RealParameter("0. 2.5") );
		bdm.setInputValue("reverseTimeArrays", "false false false true false false");
		bdm.setInputValue("conditionOnSurvival", false);
		//bdm.setInputValue("useSN", false);

		//         bdm.initAndValidate();
		//         assertEquals(Double.NEGATIVE_INFINITY, bdm.calculateLogP(), 1e-4);
		//
		bdm.setInputValue("rho", new RealParameter("0.5 1.") );
		bdm.initAndValidate();
		assertEquals(-19.29773458054159, bdm.calculateLogP(), 1e-4);   // this result is from BEAST, not double checked in R
		//
		bdm.setInputValue("rho", new RealParameter("0.01 0.01") );
		bdm.initAndValidate();
		assertEquals(-14.676660714677519, bdm.calculateLogP(), 1e-4);   // this result is from BEAST, not double checked in R

		bdm.setInputValue("rho", new RealParameter("0.") );
		bdm.setInputValue("samplingProportion", new RealParameter("0.1") );
		bdm.initAndValidate();
		// TO DO REMOVE COMMENT: I changed the value here (comes from BDSKY) - Jeremie 
		assertEquals(-21.397702433282273, bdm.calculateLogP(), 1e-4);   // this result is from BEAST, not double checked in R 
		// TO DO: remove comment
		// double checked the result in BDSKY and found -21.397702433282273 (consistent with what bdmm gives)


		//         tree.initByName(
		//                  "adjustTipHeights", false,
		//                  "newick", "(3[&type=0]: 4.5, 4[&type=0]: 4.5)[&type=0]: 1 ;",
		//                  "typeLabel", "type");
		//         bdm.setInputValue("tree", tree);
		//         bdm.setInputValue("rho", new RealParameter("0.01") );
		//         bdm.setInputValue("rhoSamplingTimes", new RealParameter("0.") );
		//         bdm.initAndValidate();
		//         assertEquals(-4.48631409442118, bdm.calculateLogP(), 1e-4);   // this result is from BEAST, not double checked in R

	}



	@Test
	public void testSingleRho() throws Exception {

		Tree tree = new TreeParser("((1[&type=0]: 4.5, 2[&type=0]: 4.5):1,3[&type=0]:5.5);",false);

		BirthDeathMigrationModelUncoloured bdm =  new BirthDeathMigrationModelUncoloured();

		bdm.setInputValue("tree", tree);

		bdm.setInputValue("typeLabel", "type");
		bdm.setInputValue("stateNumber", "1");
		bdm.setInputValue("migrationMatrix", "0.");
		bdm.setInputValue("frequencies", "1");

		bdm.setInputValue("R0", new RealParameter("1.5"));
		bdm.setInputValue("becomeUninfectiousRate", new RealParameter("1.5"));
		bdm.setInputValue("samplingProportion", new RealParameter("0.") );
		//bdm.setInputValue("useSN", false);

		bdm.setInputValue("rho", new RealParameter("0.01") );

		bdm.setInputValue("conditionOnSurvival", false);
		bdm.initAndValidate();
		assertEquals(-6.761909, bdm.calculateLogP(), 1e-4);   // this result is from R: LikConstant(2.25,1.5,0.01,c(4.5,5.5),root=1,survival=0)

		//        bdm.setInputValue("conditionOnSurvival", true);
		//        bdm.initAndValidate();
		//        assertEquals(-3.72382, bdm.calculateLogP(), 1e-4);   // this result is from R: LikConstant(2.25,1.5,0.01,c(4.5,5.5),root=1,survival=1)

		bdm.setInputValue("conditionOnSurvival", true);
		bdm.setInputValue("origin", "10");
		//         bdm.setInputValue("originIsRootEdge", false);
		bdm.initAndValidate();
		assertEquals(-7.404227, bdm.calculateLogP(), 1e-4);   // this result is from R: LikConstant(2.25,1.5,0.01,c(4.5,5.5,5.5+1e-100),root=0,survival=1)


	}

	@Test
	public void testRhoSasha() throws Exception {


		Tree tree = new TreeParser("(((((t1[&type=0]:0.4595008531,t25[&type=0]:0.4595008531)[&type=0]:0.3373053072,t23[&type=0]:0.3567584538)[&type=0]:0.007310819036,t16[&type=0]:0.3489190732)[&type=0]:0.331009529,((t18[&type=0]:0.03315384045,t14[&type=0]:0.03315384045)[&type=0]:0.5063451374,(t10[&type=0]:0.4211543131,t15[&type=0]:0.4211543131)[&type=0]:0.1183446648)[&type=0]:0.5956275305)[&type=0]:0.1158090878,((t19[&type=0]:0.9429393194,((t6[&type=0]:0.363527235,t11[&type=0]:0.4417423167)[&type=0]:0.01881829549,((((t3[&type=0]:0.3071904376,(((t24[&type=0]:0.01065209364,t13[&type=0]:0.01065209364)[&type=0]:0.06076485145,t8[&type=0]:0.07141694509)[&type=0]:0.123620245,(t22[&type=0]:0.1616119808,t2[&type=0]:0.1616119808)[&type=0]:0.03342520927)[&type=0]:0.1121532475)[&type=0]:0.24520579,t9[&type=0]:0.5523962276)[&type=0]:0.3852615426,(((t20[&type=0]:0.2935970782,(t17[&type=0]:0.06569090089,t4[&type=0]:0.06569090089)[&type=0]:0.2279061773)[&type=0]:0.08350780408,(t21[&type=0]:0.05109047139,t5[&type=0]:0.05109047139)[&type=0]:0.3260144109)[&type=0]:0.2298344132,t7[&type=0]:0.6069392955)[&type=0]:0.3307184747)[&type=0]:0.01206284377,t26[&type=0]:0.9497206139)[&type=0]:0.05755333197)[&type=0]:0.03290891884)[&type=0]:0.07263755325,t12[&type=0]:1.112820418)[&type=0]:0.1381151782);",false);

		BirthDeathMigrationModelUncoloured bdssm =  new BirthDeathMigrationModelUncoloured();
		bdssm.setInputValue("typeLabel", "type");
		bdssm.setInputValue("frequencies", "1");
		bdssm.setInputValue("migrationMatrix", "0.");
		bdssm.setInputValue("stateNumber", 1);

		bdssm.setInputValue("tree", tree);
		bdssm.setInputValue("origin", new RealParameter("2."));
		bdssm.setInputValue("conditionOnSurvival", false);

		bdssm.setInputValue("R0", new RealParameter(new Double[]{3./4.5}));
		bdssm.setInputValue("becomeUninfectiousRate", new RealParameter("4.5"));
		bdssm.setInputValue("samplingProportion", new RealParameter(new Double[]{2./4.5}));
		bdssm.setInputValue("rho", new RealParameter("0.0 0.05 0.01"));
		bdssm.setInputValue("rhoSamplingTimes","0. 1. 1.5");
		bdssm.setInputValue("reverseTimeArrays","false false false false");
		bdssm.initAndValidate();


		assertEquals(-124.96086690757612, bdssm.calculateTreeLogLikelihood(tree), 1e-2);     // this result is from BEAST, not double checked in R


		// now the same with reverse rhoSamplingTimes
		bdssm =  new BirthDeathMigrationModelUncoloured();
		bdssm.setInputValue("typeLabel", "type");
		bdssm.setInputValue("frequencies", "1");
		bdssm.setInputValue("migrationMatrix", "0.");
		bdssm.setInputValue("stateNumber", 1);

		bdssm.setInputValue("tree", tree);
		bdssm.setInputValue("origin", new RealParameter("2."));
		bdssm.setInputValue("conditionOnSurvival", false);

		bdssm.setInputValue("R0", new RealParameter(new Double[]{3./4.5}));
		bdssm.setInputValue("becomeUninfectiousRate", new RealParameter("4.5"));
		bdssm.setInputValue("samplingProportion", new RealParameter(new Double[]{2./4.5}));
		bdssm.setInputValue("rho", new RealParameter("0.0 0.05 0.01"));
		bdssm.setInputValue("rhoSamplingTimes","0. 0.5 1.");
		bdssm.setInputValue("reverseTimeArrays","false false false true");
		bdssm.initAndValidate();

		assertEquals(-124.96086690757612, bdssm.calculateTreeLogLikelihood(tree), 1e-2);     // this result is from BEAST, not double checked in R

	}


	@Test
	public void test3intsWithRho1dim() throws Exception {


		Tree tree = new TreeParser("(((((t1[&type=0]:0.4595008531,t25[&type=0]:0.4595008531)[&type=0]:0.3373053072,t23[&type=0]:0.3567584538)[&type=0]:0.007310819036,t16[&type=0]:0.3489190732)[&type=0]:0.331009529,((t18[&type=0]:0.03315384045,t14[&type=0]:0.03315384045)[&type=0]:0.5063451374,(t10[&type=0]:0.4211543131,t15[&type=0]:0.4211543131)[&type=0]:0.1183446648)[&type=0]:0.5956275305)[&type=0]:0.1158090878,((t19[&type=0]:0.9429393194,((t6[&type=0]:0.363527235,t11[&type=0]:0.4417423167)[&type=0]:0.01881829549,((((t3[&type=0]:0.3071904376,(((t24[&type=0]:0.01065209364,t13[&type=0]:0.01065209364)[&type=0]:0.06076485145,t8[&type=0]:0.07141694509)[&type=0]:0.123620245,(t22[&type=0]:0.1616119808,t2[&type=0]:0.1616119808)[&type=0]:0.03342520927)[&type=0]:0.1121532475)[&type=0]:0.24520579,t9[&type=0]:0.5523962276)[&type=0]:0.3852615426,(((t20[&type=0]:0.2935970782,(t17[&type=0]:0.06569090089,t4[&type=0]:0.06569090089)[&type=0]:0.2279061773)[&type=0]:0.08350780408,(t21[&type=0]:0.05109047139,t5[&type=0]:0.05109047139)[&type=0]:0.3260144109)[&type=0]:0.2298344132,t7[&type=0]:0.6069392955)[&type=0]:0.3307184747)[&type=0]:0.01206284377,t26[&type=0]:0.9497206139)[&type=0]:0.05755333197)[&type=0]:0.03290891884)[&type=0]:0.07263755325,t12[&type=0]:1.112820418)[&type=0]:0.1381151782);",false);

		for (int i = 3; i<4; i++){


			switch (i){
			case 0:{
				BirthDeathMigrationModelUncoloured bdssm =  new BirthDeathMigrationModelUncoloured();
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


				System.out.println("\na) Likelihood: " + bdssm.calculateTreeLogLikelihood(tree));
				assertEquals(-21.42666177086957, bdssm.calculateTreeLogLikelihood(tree), 1e-7);

			}
			case 1:{
				BirthDeathMigrationModelUncoloured bdssm =  new BirthDeathMigrationModelUncoloured();
				bdssm.setInputValue("typeLabel", "type");
				bdssm.setInputValue("frequencies", "1");
				bdssm.setInputValue("migrationMatrix", "0.");
				bdssm.setInputValue("stateNumber", 1);
				//                    bdssm.setInputValue("tolerance", 1e-16);
				//                    bdssm.setInputValue("maxEvaluations", Integer.MAX_VALUE);

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
			case 2:{
				BirthDeathMigrationModelUncoloured bdssm =  new BirthDeathMigrationModelUncoloured();
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

			case 3:{
				BirthDeathMigrationModelUncoloured bdssm =  new BirthDeathMigrationModelUncoloured();
				bdssm.setInputValue("typeLabel", "type");
				bdssm.setInputValue("frequencies", "1");
				bdssm.setInputValue("migrationMatrix", "0.");
				bdssm.setInputValue("stateNumber", 1);

				bdssm.setInputValue("tree", tree);
				//                    bdssm.setInputValue("origin", new RealParameter("2."));
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


	@Test  //1-dim test from BDSKY with
	public void testLikelihood1dim() throws Exception {

		String tree = "((3[&state=0] : 1.5, 4[&state=0] : 0.5)[&state=0] : 1 , (1[&state=0] : 2, 2[&state=0] : 1)[&state=0] : 3)[&state=0];";

		conditionOnSurvival = false;
		int maxEvals = Integer.MAX_VALUE;
		double tolerance = 1e-15;
		String locations = "1=0,2=0,3=0,4=0" ;

		double logL = bdm_likelihood(tolerance, maxEvals, "1",
				"1",
				"1.",
				tree, "1",
				"0.6666666667",null,
				"4.5",
				"0.4444444444",
				"", locations, 4, null);

		System.out.println("Birth-death result: " +logL);

		assertEquals(-44.28713883581996, logL, 1e-4);   // result from BDSKY in BEAST 9 June 2015
	}


	//    @Test
	//    public void testSALikelihoodCalculation1() throws Exception {
	//
	//        for (int i=0; i<2; i++) {
	//
	//            BirthDeathMigrationModelUncoloured model = new BirthDeathMigrationModelUncoloured();
	//
	//            ZeroBranchSATreeParser tree = (i==0)? new ZeroBranchSATreeParser("((1[&type=0]:1.0)2[&type=0]:1.0)3[&type=0]:0.0", true, false, 1)
	//                                                : new ZeroBranchSATreeParser("((1:1.5,2:0.5):0.5)3:0.0", true, false, 1);
	//
	//            model.setInputValue("tree", tree);
	//            model.setInputValue("origin", new RealParameter("10."));
	//
	//            model.setInputValue("birthRate", new RealParameter("2."));
	//            model.setInputValue("deathRate", new RealParameter("0.99"));
	//            model.setInputValue("samplingRate", new RealParameter("0.5"));
	//
	//            //        model.setInputValue("R0", new RealParameter(new Double[]{2./1.49}));
	//            //        model.setInputValue("becomeUninfectiousRate", new RealParameter("1.49"));
	//            //        model.setInputValue("samplingProportion", new RealParameter(new Double[]{0.5/1.49}) );
	//
	//            model.setInputValue("removalProbability", new RealParameter("0.9"));
	//            model.setInputValue("conditionOnSurvival", false);
	//
	//            model.setInputValue("stateNumber", "1");
	//            model.setInputValue("typeLabel", "type");
	//            model.setInputValue("migrationMatrix", "0.");
	//            model.setInputValue("frequencies", "1");
	//
	//            model.setInputValue("R0", new RealParameter("1.5"));
	//            model.setInputValue("becomeUninfectiousRate", new RealParameter("1.5"));
	//            model.setInputValue("samplingProportion", new RealParameter("0.3"));
	//
	//            model.initAndValidate();
	//
	//            // these values ate calculated with Mathematica
	//            if (i==0) assertEquals(-25.3707, model.calculateTreeLogLikelihood(tree), 1e-5); // likelihood conditioning only on parameters and origin time
	//            else    assertEquals(-22.524157039646802, model.calculateTreeLogLikelihood(tree), 1e-5);
	//        }
	//    }

	@Test
	public void testSALikelihoodCalculationWithoutAncestors() throws Exception {


		BirthDeathMigrationModelUncoloured bdm =  new BirthDeathMigrationModelUncoloured();

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

		bdm.setInputValue("origin", "10.");
		bdm.setInputValue("stateNumber", "1");
		bdm.setInputValue("migrationMatrix", "0.");
		bdm.setInputValue("frequencies", "1");

		bdm.setInputValue("R0", new RealParameter("1.5"));
		bdm.setInputValue("becomeUninfectiousRate", new RealParameter("1.5"));
		bdm.setInputValue("samplingProportion", new RealParameter("0.3") );
		bdm.setInputValue("removalProbability", new RealParameter("0.9") );
		bdm.setInputValue("conditionOnSurvival", true);

		bdm.initAndValidate();


		// likelihood conditioning on at least one sampled individual    - "true" result from BEAST one-deme SA model 09 June 2015 (DK)
		assertEquals(-25.991511346557598, bdm.calculateLogP(), 1e-4);

	}

	@Test
	public void testSALikelihoodCalculationWithoutAncestorsWithoutOrigin() throws Exception {


		BirthDeathMigrationModelUncoloured bdm =  new BirthDeathMigrationModelUncoloured();

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

		bdm.initAndValidate();


		// likelihood conditioning on at least one sampled individual    - "true" result from BEAST one-deme SA model 09 June 2015 (DK)
		assertEquals(-15.99699690815937, bdm.calculateLogP(), 1e-4);

	}

	// @Test
	public void testAmongRateChange() throws Exception {

		String tree ="((3:1.5,4:0.5):1,(1:1,2:1):3);"; //
		String orig=".1"; //
		String stateNumber = "2";
		String migrationMatrix = "0. 0.";
		String frequencies = "0.5 0.5";

		// test without rate change
		String R0 = "6 2";
		String R0AmongDemes = "1. 1.";
		String becomeUninfectiousRate = "0.5 0.5";
		String samplingProportion = "0.5 0.333333";
		String locations = "1=1,2=0,3=0,4=1" ;

		conditionOnSurvival = false;
		int maxEvals = Integer.MAX_VALUE;
		double tolerance = 1e-14;

		double logL;

		logL = bdm_likelihood(tolerance, maxEvals, stateNumber,
				migrationMatrix,
				frequencies,
				tree, orig,
				R0,R0AmongDemes,
				becomeUninfectiousRate,
				samplingProportion,
				"", locations, 4, null);


		System.out.println("Log-likelihood = " + logL);
		assertEquals(-16.683533887631224, logL, 1e-4); // result from BEAST, 8 Jan 2015

	}
	// @Test
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

		conditionOnSurvival = false;
		int maxEvals = Integer.MAX_VALUE;
		double tolerance = 1e-14;

		double logL;

		logL = bdm_likelihood(tolerance, maxEvals, stateNumber,
				migrationMatrix,
				frequencies,
				tree, orig,
				R0,R0AmongDemes,
				becomeUninfectiousRate,
				samplingProportion,
				"", locations, 4, null);


		System.out.println("Log-likelihood = " + logL);
		assertEquals(-18.82798, logL, 1e-4); // tanja's result from R

	}

	// @Test
	public void testEI_tanja() throws Exception {

		String tree ="((3:1.5,4:0.5):1,(1:2,2:1):3);"; //
		String orig="1."; //
		String stateNumber = "2";
		String migrationMatrix = "0.5 0.";
		String frequencies = "1 0";

		// test without rate change
		String R0AmongDemes = "0. 2.";
		String becomeUninfectiousRate = "0. 0.75";
		String samplingProportion = "0. 0.7";
		String locations = "1=1,2=1,3=1,4=1" ;

		conditionOnSurvival = false;
		int maxEvals = Integer.MAX_VALUE;
		double tolerance = 1e-10;

		double logL;

		logL = bdm_likelihood(tolerance, maxEvals, stateNumber,
				migrationMatrix,
				frequencies,
				tree, orig,
				"0. 0.", R0AmongDemes,
				becomeUninfectiousRate,
				samplingProportion,
				"", locations, 4, null);


		System.out.println("Log-likelihood = " + logL);
		assertEquals(-12.1441, logL, 1e-4); // tanja's result from R

	}


	// non-migration example from BDSKY
	//
	// @Test
	public void testLikelihoodCalculationNoMig() throws Exception {

		String tree ="((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);"; //
		String orig="1."; //
		String stateNumber = "1";
		String migrationMatrix = "0";
		String frequencies = "1";

		// test without rate change
		String R0 = Double.toString(4./3.);
		String becomeUninfectiousRate = "1.5";
		String samplingProportion = Double.toString(1./3.);
		String locations = "1=0,2=0,3=0,4=0" ;

		conditionOnSurvival = false;
		int maxEvals = Integer.MAX_VALUE;
		double tolerance = 1e-10;

		double logL;

		logL = bdm_likelihood(tolerance, maxEvals, stateNumber,
				migrationMatrix,
				frequencies,
				tree, orig,
				R0,null,
				becomeUninfectiousRate,
				samplingProportion,
				"", locations, 4, null);


		System.out.println("Log-likelihood = " + logL);
		assertEquals(-19.0198, logL, 1e-4); // -18.5741 if conditionOnSurvival=true

	}


	@Test  //1-dim test from BDSKY with
	public void testLikelihoodRateChange1dim() throws Exception {

		String tree = "((3[&state=0] : 1.5, 4[&state=0] : 0.5)[&state=0] : 1 , (1[&state=0] : 2, 2[&state=0] : 1)[&state=0] : 3)[&state=0];";

		conditionOnSurvival = false;
		int maxEvals = Integer.MAX_VALUE;
		double tolerance = 1e-10;
		String locations = "1=0,2=0,3=0,4=0" ;
		String intervalTimes = "0. 3.";

		double logL = bdm_likelihood(tolerance, maxEvals, "1",
				"1",
				"1.",
				tree, "1",
				"0.6666666667 1.3333333334",null,
				"4.5 1.5",
				"0.4444444444 0.3333333333",
				"", locations, 4, intervalTimes);

		//        bdssm.setInputValue("conditionOnSurvival", false);

		//        bdssm.setInputValue("birthRate", new RealParameter("3. 2."));
		//        bdssm.setInputValue("deathRate", new RealParameter("2.5 1."));
		//        bdssm.setInputValue("samplingRate", new RealParameter("2. 0.5"));
		//        bdssm.setInputValue("intervalTimes", new RealParameter("0. 3."));
		//        bdssm.setInputValue("samplingRateChangeTimes", new RealParameter("3. 0."));
		//
		//        bdssm.setInputValue("reverseTimeArrays", "false false true false");
		//
		//        bdssm.initAndValidate();
		//        bdssm.printTempResults = true;
		System.out.println("Birth-death result: " +logL);

		assertEquals(-33.7573, logL, 1e-2);
	}

	public void testLikelihoodCalculationWithMig() throws Exception{

		String tree = "(((((t1:0.9364152609,t2:0.3325266235):0.6741873075,t3:0.05631965678):0.02027716581,t4:0.7622636168):0.03208635212,((((t5:0.6713889462,(t6:0.1857663356,t7:0.4931635624):0.8884018492):0.6485969226,((t8:0.6032408362,t9:0.8932013339):0.4847603254,t10:0.6412034486):0.1307068868):0.1782681448,((((((t11:0.5446377769,t12:0.5245424444):0.01255171769,t13:0.709872612):0.1737786294,t14:0.5846257175):0.2142320445,(t15:0.7903199508,t16:0.1911238336):0.7025476368):0.8200157876,t17:0.6692489227):0.01014118479,((t18:0.4153864153,t19:0.9988613673):0.5179155183,((t20:0.9572573467,t21:0.1383048228):0.7580099995,t22:0.5800547246):0.4191723403):0.4258894029):0.1695718046):0.8316473491,((((t23:0.6841176334,((t24:0.1905560789,t25:0.9144611149):0.1210783592,t26:0.4941913132):0.642020765):0.9268047926,((t27:0.8871689022,((t28:0.2958402268,t29:0.149307796):0.7057882755,(t30:0.1131704114,t31:0.4348528353):0.4895360142):0.9755446969):0.2343815062,t32:0.8258521033):0.2539390386):0.1383914386,((((t33:0.2146211416,t34:0.8262746611):0.9609895966,((((t35:0.4596964282,t36:0.05147929071):0.5753136226,(t37:0.9932728133,t38:0.5782027193):0.616813526):0.8144772681,((t39:0.2140128147,t40:0.9378008009):0.4376288333,(t41:0.01509191399,(t42:0.7252295294,t43:0.4592927478):0.4014166105):0.9703455286):0.3391084191):0.2390605409,((t44:0.9023445742,(t45:0.7600001141,t46:0.8390259156):0.7530289539):0.8280177859,(((t47:0.7356161084,((t48:0.8990668394,t49:0.476900334):0.401515421,t50:0.8970352295):0.6753761517):0.4681597212,t51:0.9660184374):0.4445622298,(t52:0.5702694552,((t53:0.7867721654,(t54:0.3036163356,t55:0.2289324626):0.101905775):0.9686289642,t56:0.3421533015):0.4520260729):0.6061327022):0.6351576883):0.5171940143):0.8868461328):0.6640388558,(((t57:0.5333570447,t58:0.3410403242):0.5509593328,t59:0.7223718707):0.8887142881,(t60:0.3024317482,(t61:0.7791908968,t62:0.4357111922):0.4007126684):0.7111749286):0.2846378125):0.8696146321,(t63:0.3404058055,(t64:0.5541903777,t65:0.8805724247):0.2429451609):0.7887147672):0.08329026634):0.6959039369,((t66:0.1594602023,(t67:0.3549112917,t68:0.6501219161):0.3202879052):0.7210102496,(((t69:0.8546685255,(t70:0.4183590303,(t71:0.7280786086,t72:0.0501071047):0.2249922389):0.5513544895):0.5412701196,(t73:0.03845808865,t74:0.01786546269):0.2918240293):0.2964613114,t75:0.6997288649):0.5851331023):0.6263443602):0.382944342):0.9287640562):0.7195440128,(((((((((t76:0.908999736,t77:0.7081224946):0.509465781,t78:0.2401761867):0.6602909768,((t79:0.4055437755,t80:0.6022770065):0.444950735,t81:0.3619997192):0.2092859314):0.07271602307,((t82:0.7013580918,t83:0.3649420871):0.7937776847,t84:0.5524289906):0.2912156314):0.1191252382,((t85:0.1511951322,t86:0.6732844713):0.1495890303,t87:0.4700718245):0.1184274761):0.2547142582,(((((t88:0.4415043984,(t89:0.389181149,(t90:0.08476964966,t91:0.002604217967):0.753279258):0.6555457374):0.1691021193,(t92:0.6564390776,t93:0.2826642371):0.9404367092):0.1501046501,(t94:0.9856409656,t95:0.6226283619):0.4212368946):0.782352773,t96:0.7864871565):0.7963304925,(((t97:0.491204008,((t98:0.455413609,t99:0.1336613444):0.6196382644,t100:0.06891984516):0.07895972207):0.4890337239,t101:0.7780050633):0.5802606251,t102:0.5368855898):0.7450205528):0.2824418589):0.6647088744,((((t103:0.6554223788,(t104:0.1252187097,(t105:0.5746273743,t106:0.9021635554):0.2750798387):0.8763033026):0.8557207501,t107:0.5860088668):0.182998168,t108:0.1329866969):0.4255650612,(t109:0.7125762652,((t110:0.05844185059,t111:0.8927688843):0.596926763,t112:0.1202381498):0.4700582502):0.1470693348):0.4912134372):0.4780910299,((((((t113:0.2291344877,((t114:0.07061028155,(t115:0.3895098094,t116:0.8108067587):0.7016420467):0.5853939084,t117:0.801262571):0.8635679013):0.8172595853,t118:0.5230401708):0.2267086026,(t119:0.3480776774,t120:0.3773546459):0.09923658706):0.6337174773,(t121:0.3583609944,(t122:0.9210918313,t123:0.04823979852):0.4626455419):0.3761554039):0.9973847859,(t124:0.9720817744,(t125:0.5189394697,((t126:0.7307705784,t127:0.4378616042):0.7934694779,t128:0.3696364786):0.7001878037):0.6506283809):0.399283418):0.01584472018,t129:0.4401801128):0.8105186713):0.03882670798,(((((t130:0.1942299611,t131:0.1218361468):0.4073948862,(((t132:0.8818691929,t133:0.2246095352):0.4217373263,(t134:0.3414793743,(t135:0.1749610452,(t136:0.4462995038,t137:0.71709141):0.1657548398):0.4825986377):0.5821935786):0.4052431665,t138:0.06300760945):0.5269570656):0.884282737,((t139:0.8558092408,(t140:0.3321465107,t141:0.6650843567):0.5518979621):0.257976112,(t142:0.3446148888,(((((t143:0.5418575266,(t144:0.6116074673,t145:0.006826744182):0.9411114273):0.6997815147,(t146:0.4408475324,(((t147:0.3645979795,t148:0.0404281714):0.2886804543,t149:0.8369950363):0.9826068394,t150:0.3104573903):0.2021026732):0.257088437):0.843516761,(t151:0.7693858896,(t152:0.9971931365,t153:0.146940355):0.305027812):0.6522608057):0.04420206766,(t154:0.5183477278,t155:0.0003721236717):0.03044462833):0.7406412989,(t156:0.7852081677,(((((((t157:0.1363189702,(t158:0.88613495,t159:0.4112141535):0.685479136):0.05440945271,(t160:0.330911414,t161:0.8818563768):0.7312560759):0.1806290515,t162:0.907872692):0.5071503669,t163:0.4420965984):0.2860673463,(t164:0.9970039094,t165:0.09012954589):0.3197214403):0.5374525476,(t166:0.9520867753,t167:0.4960751149):0.6838007318):0.7549668197,(((t168:0.2440482709,t169:0.633632384):0.8411958364,(t170:0.1448847703,(t171:0.4834090637,t172:0.409381151):0.2860563635):0.6396883058):0.7715818523,(t173:0.6423769889,(t174:0.9425992649,t175:0.706059468):0.1856973255):0.09803533414):0.7050393047):0.2558317431):0.1565285209):0.3077181939):0.5945158091):0.1038589557):0.3777994341,((((t176:0.8200896268,t177:0.8051551071):0.4032325067,t178:0.479055285):0.3715889419,(((t179:0.04895586846,t180:0.9665267451):0.3085221592,(t181:0.3709749565,t182:0.3955604453):0.05966877262):0.9807001806,((t183:0.1330345676,t184:0.356874147):0.05452865968,t185:0.85592992):0.6030724624):0.2495698626):0.833886015,(t186:0.581214027,(((t187:0.8901591096,t188:0.9351359161):0.9671047614,(t189:0.2125755912,(t190:0.0701806522,t191:0.07933806628):0.780095432):0.8840237022):0.02146822144,(t192:0.3327144748,(t193:0.3291403668,t194:0.1709528794):0.1115758724):0.3335891911):0.9738473571):0.7402053783):0.8972619711):0.1375731654,((t195:0.1412698505,(t196:0.3687814591,t197:0.337743812):0.37485075):0.8553067574,(t198:0.5955553057,(t199:0.6654327868,t200:0.9634261283):0.977098627):0.7736232406):0.911330268):0.771963502):0.5990205682);";
		int nrTaxa = 200;
		int[] states = new int[]{1, 2, 1, 2, 1, 2, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 2, 1, 2, 2, 1, 2, 1, 2, 2, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 1, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 1, 2, 1, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 2, 2, 1, 2, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 1, 2, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 1, 2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 2, 2, 2, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 1, 2, 2};

		String orig = ".1";

		String locations = "" ;
		for (int i=1; i<=states.length; i++){
			locations = locations + "t" + i + "=" + (states[i-1]-1) + (i<states.length?",":"");
		}

		String stateNumber = "2";
		String migrationMatrix = "0.2 0.1";
		String frequencies = "0.5 0.5";

		String R0 = Double.toString(4./3.) + " " + Double.toString(5.) ;
		String becomeUninfectiousRate = "1.5 1.25";
		String samplingProportion = Double.toString(1./3.) + " " + Double.toString(1./2.);


		int maxEvals = Integer.MAX_VALUE;
		double tolerance = 1e-13;

		double logL;

		conditionOnSurvival = false;

		logL = bdm_likelihood(tolerance, maxEvals, stateNumber,
				migrationMatrix,
				frequencies,
				tree, orig,
				R0,  null,
				becomeUninfectiousRate,
				samplingProportion,
				"t", locations, nrTaxa, null);

		System.out.println("Log-likelihood (conditioned on survival) " + logL + "\t");
		assertEquals(-667.885, logL, 1e-3);

	}

	public void testLikelihoodCalculationInfAmongDemesSymmetric() throws Exception{

		Boolean r0param = false;
		Boolean infectionAmongDemes = true;
		conditionOnSurvival = true;

		String newick = "((3:1.5,4:0.5):1,(1:2,2:1):3);";
		//= "(((((t1:0.9364152609,t2:0.3325266235):0.6741873075,t3:0.05631965678):0.02027716581,t4:0.7622636168):0.03208635212,((((t5:0.6713889462,(t6:0.1857663356,t7:0.4931635624):0.8884018492):0.6485969226,((t8:0.6032408362,t9:0.8932013339):0.4847603254,t10:0.6412034486):0.1307068868):0.1782681448,((((((t11:0.5446377769,t12:0.5245424444):0.01255171769,t13:0.709872612):0.1737786294,t14:0.5846257175):0.2142320445,(t15:0.7903199508,t16:0.1911238336):0.7025476368):0.8200157876,t17:0.6692489227):0.01014118479,((t18:0.4153864153,t19:0.9988613673):0.5179155183,((t20:0.9572573467,t21:0.1383048228):0.7580099995,t22:0.5800547246):0.4191723403):0.4258894029):0.1695718046):0.8316473491,((((t23:0.6841176334,((t24:0.1905560789,t25:0.9144611149):0.1210783592,t26:0.4941913132):0.642020765):0.9268047926,((t27:0.8871689022,((t28:0.2958402268,t29:0.149307796):0.7057882755,(t30:0.1131704114,t31:0.4348528353):0.4895360142):0.9755446969):0.2343815062,t32:0.8258521033):0.2539390386):0.1383914386,((((t33:0.2146211416,t34:0.8262746611):0.9609895966,((((t35:0.4596964282,t36:0.05147929071):0.5753136226,(t37:0.9932728133,t38:0.5782027193):0.616813526):0.8144772681,((t39:0.2140128147,t40:0.9378008009):0.4376288333,(t41:0.01509191399,(t42:0.7252295294,t43:0.4592927478):0.4014166105):0.9703455286):0.3391084191):0.2390605409,((t44:0.9023445742,(t45:0.7600001141,t46:0.8390259156):0.7530289539):0.8280177859,(((t47:0.7356161084,((t48:0.8990668394,t49:0.476900334):0.401515421,t50:0.8970352295):0.6753761517):0.4681597212,t51:0.9660184374):0.4445622298,(t52:0.5702694552,((t53:0.7867721654,(t54:0.3036163356,t55:0.2289324626):0.101905775):0.9686289642,t56:0.3421533015):0.4520260729):0.6061327022):0.6351576883):0.5171940143):0.8868461328):0.6640388558,(((t57:0.5333570447,t58:0.3410403242):0.5509593328,t59:0.7223718707):0.8887142881,(t60:0.3024317482,(t61:0.7791908968,t62:0.4357111922):0.4007126684):0.7111749286):0.2846378125):0.8696146321,(t63:0.3404058055,(t64:0.5541903777,t65:0.8805724247):0.2429451609):0.7887147672):0.08329026634):0.6959039369,((t66:0.1594602023,(t67:0.3549112917,t68:0.6501219161):0.3202879052):0.7210102496,(((t69:0.8546685255,(t70:0.4183590303,(t71:0.7280786086,t72:0.0501071047):0.2249922389):0.5513544895):0.5412701196,(t73:0.03845808865,t74:0.01786546269):0.2918240293):0.2964613114,t75:0.6997288649):0.5851331023):0.6263443602):0.382944342):0.9287640562):0.7195440128,(((((((((t76:0.908999736,t77:0.7081224946):0.509465781,t78:0.2401761867):0.6602909768,((t79:0.4055437755,t80:0.6022770065):0.444950735,t81:0.3619997192):0.2092859314):0.07271602307,((t82:0.7013580918,t83:0.3649420871):0.7937776847,t84:0.5524289906):0.2912156314):0.1191252382,((t85:0.1511951322,t86:0.6732844713):0.1495890303,t87:0.4700718245):0.1184274761):0.2547142582,(((((t88:0.4415043984,(t89:0.389181149,(t90:0.08476964966,t91:0.002604217967):0.753279258):0.6555457374):0.1691021193,(t92:0.6564390776,t93:0.2826642371):0.9404367092):0.1501046501,(t94:0.9856409656,t95:0.6226283619):0.4212368946):0.782352773,t96:0.7864871565):0.7963304925,(((t97:0.491204008,((t98:0.455413609,t99:0.1336613444):0.6196382644,t100:0.06891984516):0.07895972207):0.4890337239,t101:0.7780050633):0.5802606251,t102:0.5368855898):0.7450205528):0.2824418589):0.6647088744,((((t103:0.6554223788,(t104:0.1252187097,(t105:0.5746273743,t106:0.9021635554):0.2750798387):0.8763033026):0.8557207501,t107:0.5860088668):0.182998168,t108:0.1329866969):0.4255650612,(t109:0.7125762652,((t110:0.05844185059,t111:0.8927688843):0.596926763,t112:0.1202381498):0.4700582502):0.1470693348):0.4912134372):0.4780910299,((((((t113:0.2291344877,((t114:0.07061028155,(t115:0.3895098094,t116:0.8108067587):0.7016420467):0.5853939084,t117:0.801262571):0.8635679013):0.8172595853,t118:0.5230401708):0.2267086026,(t119:0.3480776774,t120:0.3773546459):0.09923658706):0.6337174773,(t121:0.3583609944,(t122:0.9210918313,t123:0.04823979852):0.4626455419):0.3761554039):0.9973847859,(t124:0.9720817744,(t125:0.5189394697,((t126:0.7307705784,t127:0.4378616042):0.7934694779,t128:0.3696364786):0.7001878037):0.6506283809):0.399283418):0.01584472018,t129:0.4401801128):0.8105186713):0.03882670798,(((((t130:0.1942299611,t131:0.1218361468):0.4073948862,(((t132:0.8818691929,t133:0.2246095352):0.4217373263,(t134:0.3414793743,(t135:0.1749610452,(t136:0.4462995038,t137:0.71709141):0.1657548398):0.4825986377):0.5821935786):0.4052431665,t138:0.06300760945):0.5269570656):0.884282737,((t139:0.8558092408,(t140:0.3321465107,t141:0.6650843567):0.5518979621):0.257976112,(t142:0.3446148888,(((((t143:0.5418575266,(t144:0.6116074673,t145:0.006826744182):0.9411114273):0.6997815147,(t146:0.4408475324,(((t147:0.3645979795,t148:0.0404281714):0.2886804543,t149:0.8369950363):0.9826068394,t150:0.3104573903):0.2021026732):0.257088437):0.843516761,(t151:0.7693858896,(t152:0.9971931365,t153:0.146940355):0.305027812):0.6522608057):0.04420206766,(t154:0.5183477278,t155:0.0003721236717):0.03044462833):0.7406412989,(t156:0.7852081677,(((((((t157:0.1363189702,(t158:0.88613495,t159:0.4112141535):0.685479136):0.05440945271,(t160:0.330911414,t161:0.8818563768):0.7312560759):0.1806290515,t162:0.907872692):0.5071503669,t163:0.4420965984):0.2860673463,(t164:0.9970039094,t165:0.09012954589):0.3197214403):0.5374525476,(t166:0.9520867753,t167:0.4960751149):0.6838007318):0.7549668197,(((t168:0.2440482709,t169:0.633632384):0.8411958364,(t170:0.1448847703,(t171:0.4834090637,t172:0.409381151):0.2860563635):0.6396883058):0.7715818523,(t173:0.6423769889,(t174:0.9425992649,t175:0.706059468):0.1856973255):0.09803533414):0.7050393047):0.2558317431):0.1565285209):0.3077181939):0.5945158091):0.1038589557):0.3777994341,((((t176:0.8200896268,t177:0.8051551071):0.4032325067,t178:0.479055285):0.3715889419,(((t179:0.04895586846,t180:0.9665267451):0.3085221592,(t181:0.3709749565,t182:0.3955604453):0.05966877262):0.9807001806,((t183:0.1330345676,t184:0.356874147):0.05452865968,t185:0.85592992):0.6030724624):0.2495698626):0.833886015,(t186:0.581214027,(((t187:0.8901591096,t188:0.9351359161):0.9671047614,(t189:0.2125755912,(t190:0.0701806522,t191:0.07933806628):0.780095432):0.8840237022):0.02146822144,(t192:0.3327144748,(t193:0.3291403668,t194:0.1709528794):0.1115758724):0.3335891911):0.9738473571):0.7402053783):0.8972619711):0.1375731654,((t195:0.1412698505,(t196:0.3687814591,t197:0.337743812):0.37485075):0.8553067574,(t198:0.5955553057,(t199:0.6654327868,t200:0.9634261283):0.977098627):0.7736232406):0.911330268):0.771963502):0.5990205682);";
		String prefixname = "";//"t";

		int nrTaxa = 4 ;//200;
		int[] states = new int[]{2,1,1,2};//{2,1,1,2};//
		//{1, 2, 1, 2, 1, 2, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 2, 1, 2, 2, 1, 2, 1, 2, 2, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 1, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 1, 2, 1, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 2, 2, 1, 2, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 1, 2, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 1, 2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 2, 2, 2, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 1, 2, 2};

		String orig = "1.";

		String locations = "" ;
		for (int i=1; i<=states.length; i++){
			locations = locations + prefixname + i + "=" + (states[i-1]-1) + (i<states.length?",":"");
		}

		String stateNumber = "2";
		String migrationMatrix = "0. 0." ; //"1. 1.";
		String frequencies = "0.5 0.5";

		String birth = "2. 2.";
		String birthRateAmongDemes = "1. 1."; // "0. 0.";
		String deathRate = "0.4 0.4";
		String samplingRate = "0.1 0.1";

		String R0 = "4. 4.";
		String R0aD = "2. 2."; //"0. 0.";
		String delta = "0.5 0.5";
		String s = "0.2 0.2";

		int maxEvals = Integer.MAX_VALUE;
		double tolerance = 1e-14;

		double logL;

		BirthDeathMigrationModelUncoloured bdm =  new BirthDeathMigrationModelUncoloured();

		Tree tree = new TreeParser();
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
		bdm.setInputValue("migrationMatrix", migrationMatrix);
		bdm.setInputValue("frequencies", frequencies);

		if (r0param){
			bdm.setInputValue("R0", R0);
			if (infectionAmongDemes) bdm.setInputValue("R0AmongDemes", R0aD);
			bdm.setInputValue("becomeUninfectiousRate", delta);
			bdm.setInputValue("samplingProportion", s);

		}
		else{

			bdm.setInputValue("birthRate", birth);
			if (infectionAmongDemes) bdm.setInputValue("birthRateAmongDemes", birthRateAmongDemes);
			bdm.setInputValue("deathRate", deathRate);
			bdm.setInputValue("samplingRate", samplingRate);
		}
		bdm.setInputValue("maxEvaluations", maxEvals);
		bdm.setInputValue("conditionOnSurvival", conditionOnSurvival);

		bdm.setInputValue("relTolerance", tolerance);

		bdm.initAndValidate();

		long startTime = System.currentTimeMillis();
		logL = bdm.calculateLogP(); //calculateTreeLogLikelihood(coltree.getUncolouredTree());
		runtime = System.currentTimeMillis() - startTime;
		maxEvalsUsed = bdm.maxEvalsUsed;
		//        assertEquals( -19.76, logL,1e-2);

		if (infectionAmongDemes) {
			if (conditionOnSurvival) assertEquals(-28.05224, logL, 1e-3);
			else assertEquals(-28.18969, logL, 1e-3);
		}

		else{
			if (conditionOnSurvival) assertEquals(-19.82658,logL, 1e-3);
			else assertEquals(-20.035714,logL, 1e-3);
		}
		System.out.println("Log-likelihood ("+ (conditionOnSurvival?"":"NOT ") +"conditioned on survival) " + logL + "\t");

	}

	public void testInfAmongDemesNEW() throws Exception{

		conditionOnSurvival = true;

		//sim.bdtypes.stt.taxa(4,l=rbind(c(2,1),c(1,2)),d=c(1,1),s=c(.5,.5))

		String newick = "((t3[&type=1]:0.004214277605,t4[&type=1]:0.02157681391):0.229186993,(t2[&type=0]:0.624713651,t1[&type=1]:1.347400211):0.06231047755);";
		String prefixname = "t";

		int nrTaxa = 4 ;

		String orig = "0.02686563367";


		String stateNumber = "2";
		String frequencies = "0.5 0.5";

		String birth = "2. 2.";
		String birthRateAmongDemes = "1. 1."; // "0. 0.";
		String deathRate = "0.5 0.5";
		String samplingRate = "0.5 0.5";

		double logL;

		BirthDeathMigrationModelUncoloured bdm =  new BirthDeathMigrationModelUncoloured();

		//        ArrayList<Taxon> taxa = new ArrayList<Taxon>();
		//
		//        for (int i=1; i<=nrTaxa; i++){
		//            taxa.add(new Taxon(prefixname+i));
		//        }
		//
		//        TraitSet trait = new TraitSet();
		//        trait.setInputValue("taxa", new TaxonSet(taxa));
		//        trait.setInputValue("value", locations);
		//        trait.setInputValue("traitname", "type");
		//        trait.initAndValidate();

		Tree tree = new TreeParser();
		tree.setInputValue("adjustTipHeights", "false");
		tree.setInputValue("IsLabelledNewick", "true");
		tree.setInputValue("newick", newick);
		//        tree.setInputValue("trait", trait);
		tree.initAndValidate();

		bdm.setInputValue("tree", tree);
		//        bdm.setInputValue("tiptypes", trait);
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

		logL = bdm.calculateLogP();

		//        if (conditionOnSurvival) assertEquals(-28.05224, logL, 1e-3);
		//        else assertEquals(-28.18969, logL, 1e-3);

		System.out.println("Log-likelihood ("+ (conditionOnSurvival?"":"NOT ") +"conditioned on survival) " + logL + "\t");

	}

	public void testLikelihoodCalculationInfAmongDemes() throws Exception{

		Boolean r0param = false;
		Boolean infectionAmongDemes = true;
		conditionOnSurvival = true;

		String newick = "((3:1.5,4:0.5):1,(1:2,2:1):3);";
		//= "(((((t1:0.9364152609,t2:0.3325266235):0.6741873075,t3:0.05631965678):0.02027716581,t4:0.7622636168):0.03208635212,((((t5:0.6713889462,(t6:0.1857663356,t7:0.4931635624):0.8884018492):0.6485969226,((t8:0.6032408362,t9:0.8932013339):0.4847603254,t10:0.6412034486):0.1307068868):0.1782681448,((((((t11:0.5446377769,t12:0.5245424444):0.01255171769,t13:0.709872612):0.1737786294,t14:0.5846257175):0.2142320445,(t15:0.7903199508,t16:0.1911238336):0.7025476368):0.8200157876,t17:0.6692489227):0.01014118479,((t18:0.4153864153,t19:0.9988613673):0.5179155183,((t20:0.9572573467,t21:0.1383048228):0.7580099995,t22:0.5800547246):0.4191723403):0.4258894029):0.1695718046):0.8316473491,((((t23:0.6841176334,((t24:0.1905560789,t25:0.9144611149):0.1210783592,t26:0.4941913132):0.642020765):0.9268047926,((t27:0.8871689022,((t28:0.2958402268,t29:0.149307796):0.7057882755,(t30:0.1131704114,t31:0.4348528353):0.4895360142):0.9755446969):0.2343815062,t32:0.8258521033):0.2539390386):0.1383914386,((((t33:0.2146211416,t34:0.8262746611):0.9609895966,((((t35:0.4596964282,t36:0.05147929071):0.5753136226,(t37:0.9932728133,t38:0.5782027193):0.616813526):0.8144772681,((t39:0.2140128147,t40:0.9378008009):0.4376288333,(t41:0.01509191399,(t42:0.7252295294,t43:0.4592927478):0.4014166105):0.9703455286):0.3391084191):0.2390605409,((t44:0.9023445742,(t45:0.7600001141,t46:0.8390259156):0.7530289539):0.8280177859,(((t47:0.7356161084,((t48:0.8990668394,t49:0.476900334):0.401515421,t50:0.8970352295):0.6753761517):0.4681597212,t51:0.9660184374):0.4445622298,(t52:0.5702694552,((t53:0.7867721654,(t54:0.3036163356,t55:0.2289324626):0.101905775):0.9686289642,t56:0.3421533015):0.4520260729):0.6061327022):0.6351576883):0.5171940143):0.8868461328):0.6640388558,(((t57:0.5333570447,t58:0.3410403242):0.5509593328,t59:0.7223718707):0.8887142881,(t60:0.3024317482,(t61:0.7791908968,t62:0.4357111922):0.4007126684):0.7111749286):0.2846378125):0.8696146321,(t63:0.3404058055,(t64:0.5541903777,t65:0.8805724247):0.2429451609):0.7887147672):0.08329026634):0.6959039369,((t66:0.1594602023,(t67:0.3549112917,t68:0.6501219161):0.3202879052):0.7210102496,(((t69:0.8546685255,(t70:0.4183590303,(t71:0.7280786086,t72:0.0501071047):0.2249922389):0.5513544895):0.5412701196,(t73:0.03845808865,t74:0.01786546269):0.2918240293):0.2964613114,t75:0.6997288649):0.5851331023):0.6263443602):0.382944342):0.9287640562):0.7195440128,(((((((((t76:0.908999736,t77:0.7081224946):0.509465781,t78:0.2401761867):0.6602909768,((t79:0.4055437755,t80:0.6022770065):0.444950735,t81:0.3619997192):0.2092859314):0.07271602307,((t82:0.7013580918,t83:0.3649420871):0.7937776847,t84:0.5524289906):0.2912156314):0.1191252382,((t85:0.1511951322,t86:0.6732844713):0.1495890303,t87:0.4700718245):0.1184274761):0.2547142582,(((((t88:0.4415043984,(t89:0.389181149,(t90:0.08476964966,t91:0.002604217967):0.753279258):0.6555457374):0.1691021193,(t92:0.6564390776,t93:0.2826642371):0.9404367092):0.1501046501,(t94:0.9856409656,t95:0.6226283619):0.4212368946):0.782352773,t96:0.7864871565):0.7963304925,(((t97:0.491204008,((t98:0.455413609,t99:0.1336613444):0.6196382644,t100:0.06891984516):0.07895972207):0.4890337239,t101:0.7780050633):0.5802606251,t102:0.5368855898):0.7450205528):0.2824418589):0.6647088744,((((t103:0.6554223788,(t104:0.1252187097,(t105:0.5746273743,t106:0.9021635554):0.2750798387):0.8763033026):0.8557207501,t107:0.5860088668):0.182998168,t108:0.1329866969):0.4255650612,(t109:0.7125762652,((t110:0.05844185059,t111:0.8927688843):0.596926763,t112:0.1202381498):0.4700582502):0.1470693348):0.4912134372):0.4780910299,((((((t113:0.2291344877,((t114:0.07061028155,(t115:0.3895098094,t116:0.8108067587):0.7016420467):0.5853939084,t117:0.801262571):0.8635679013):0.8172595853,t118:0.5230401708):0.2267086026,(t119:0.3480776774,t120:0.3773546459):0.09923658706):0.6337174773,(t121:0.3583609944,(t122:0.9210918313,t123:0.04823979852):0.4626455419):0.3761554039):0.9973847859,(t124:0.9720817744,(t125:0.5189394697,((t126:0.7307705784,t127:0.4378616042):0.7934694779,t128:0.3696364786):0.7001878037):0.6506283809):0.399283418):0.01584472018,t129:0.4401801128):0.8105186713):0.03882670798,(((((t130:0.1942299611,t131:0.1218361468):0.4073948862,(((t132:0.8818691929,t133:0.2246095352):0.4217373263,(t134:0.3414793743,(t135:0.1749610452,(t136:0.4462995038,t137:0.71709141):0.1657548398):0.4825986377):0.5821935786):0.4052431665,t138:0.06300760945):0.5269570656):0.884282737,((t139:0.8558092408,(t140:0.3321465107,t141:0.6650843567):0.5518979621):0.257976112,(t142:0.3446148888,(((((t143:0.5418575266,(t144:0.6116074673,t145:0.006826744182):0.9411114273):0.6997815147,(t146:0.4408475324,(((t147:0.3645979795,t148:0.0404281714):0.2886804543,t149:0.8369950363):0.9826068394,t150:0.3104573903):0.2021026732):0.257088437):0.843516761,(t151:0.7693858896,(t152:0.9971931365,t153:0.146940355):0.305027812):0.6522608057):0.04420206766,(t154:0.5183477278,t155:0.0003721236717):0.03044462833):0.7406412989,(t156:0.7852081677,(((((((t157:0.1363189702,(t158:0.88613495,t159:0.4112141535):0.685479136):0.05440945271,(t160:0.330911414,t161:0.8818563768):0.7312560759):0.1806290515,t162:0.907872692):0.5071503669,t163:0.4420965984):0.2860673463,(t164:0.9970039094,t165:0.09012954589):0.3197214403):0.5374525476,(t166:0.9520867753,t167:0.4960751149):0.6838007318):0.7549668197,(((t168:0.2440482709,t169:0.633632384):0.8411958364,(t170:0.1448847703,(t171:0.4834090637,t172:0.409381151):0.2860563635):0.6396883058):0.7715818523,(t173:0.6423769889,(t174:0.9425992649,t175:0.706059468):0.1856973255):0.09803533414):0.7050393047):0.2558317431):0.1565285209):0.3077181939):0.5945158091):0.1038589557):0.3777994341,((((t176:0.8200896268,t177:0.8051551071):0.4032325067,t178:0.479055285):0.3715889419,(((t179:0.04895586846,t180:0.9665267451):0.3085221592,(t181:0.3709749565,t182:0.3955604453):0.05966877262):0.9807001806,((t183:0.1330345676,t184:0.356874147):0.05452865968,t185:0.85592992):0.6030724624):0.2495698626):0.833886015,(t186:0.581214027,(((t187:0.8901591096,t188:0.9351359161):0.9671047614,(t189:0.2125755912,(t190:0.0701806522,t191:0.07933806628):0.780095432):0.8840237022):0.02146822144,(t192:0.3327144748,(t193:0.3291403668,t194:0.1709528794):0.1115758724):0.3335891911):0.9738473571):0.7402053783):0.8972619711):0.1375731654,((t195:0.1412698505,(t196:0.3687814591,t197:0.337743812):0.37485075):0.8553067574,(t198:0.5955553057,(t199:0.6654327868,t200:0.9634261283):0.977098627):0.7736232406):0.911330268):0.771963502):0.5990205682);";
		String prefixname = "";//"t";

		int nrTaxa = 4 ;//200;
		int[] states = new int[]{1,2,2,1};//{2,1,1,2};//
		//{1, 2, 1, 2, 1, 2, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 2, 1, 2, 2, 1, 2, 1, 2, 2, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 1, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 1, 2, 1, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 2, 2, 1, 2, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 1, 2, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 1, 2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 2, 2, 2, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 1, 2, 2};

		String orig = "1.";

		String locations = "" ;
		for (int i=1; i<=states.length; i++){
			locations = locations + prefixname + i + "=" + (states[i-1]-1) + (i<states.length?",":"");
		}

		String stateNumber = "2";
		String migrationMatrix = "0.2 0.1" ; //"0.2 0.1" ; //"1. 1.";
		String frequencies = "0.5 0.5";

		BirthDeathMigrationModelUncoloured bdm =  new BirthDeathMigrationModelUncoloured();

		if (!r0param){
			String birth = "2. 6.25";
			String birthRateAmongDemes = "0.2 0.1"; // "0. 0.";
			String deathRate = "1.2 0.625";
			String samplingRate = "0.3 0.625";


			bdm.setInputValue("birthRate", birth);
			if (infectionAmongDemes) {
				bdm.setInputValue("birthRateAmongDemes", birthRateAmongDemes);
				bdm.setInputValue("migrationMatrix", "0 0");
			} else {
				bdm.setInputValue("migrationMatrix", migrationMatrix);
			}
			bdm.setInputValue("deathRate", deathRate);
			bdm.setInputValue("samplingRate", samplingRate);

		}
		else {

			String R0 = "1.33333333 5.";
			String R0aD =  "0.13333333 0.08"; //"0. 0.";//
			String delta = "1.5 1.25";
			String s = "0.2 0.5";

			bdm.setInputValue("R0", R0);
			if (infectionAmongDemes){
				bdm.setInputValue("R0AmongDemes", R0aD);
				bdm.setInputValue("migrationMatrix", "0 0");
			}   else{
				bdm.setInputValue("migrationMatrix", migrationMatrix);

			}
			bdm.setInputValue("becomeUninfectiousRate", delta);
			bdm.setInputValue("samplingProportion", s);

		}
		int maxEvals = Integer.MAX_VALUE;
		double tolerance = 1e-14;

		double logL;


		Tree tree = new TreeParser();
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

		bdm.setInputValue("maxEvaluations", maxEvals);
		bdm.setInputValue("relTolerance", tolerance);

		bdm.initAndValidate();


		long startTime = System.currentTimeMillis();
		logL = bdm.calculateLogP(); //calculateTreeLogLikelihood(coltree.getUncolouredTree());
		runtime = System.currentTimeMillis() - startTime;
		maxEvalsUsed = bdm.maxEvalsUsed;
		//        assertEquals( -19.76, logL,1e-2);

		if (infectionAmongDemes) {
			//            if (bdm.migrationMatrix.get().getValue(0) !=0.) {  // can't check this in R yet
			//                if (conditionOnSurvival) assertEquals(-25.095773220763093, logL, 1e-3);
			//            }
			//            else {
			if (conditionOnSurvival) assertEquals(-26.7939, logL, 1e-5);
			else assertEquals(-27.08992, logL, 1e-5);  //result from R
			//            }
		}

		else{
			if (conditionOnSurvival) assertEquals( -25.34393,logL, 1e-5);  //result from R
			else assertEquals(-25.64794,logL, 1e-5);//result from R
		}
		System.out.println("Log-likelihood ("+ (conditionOnSurvival?"":"NOT ") +"conditioned on survival) " + logL + "\t");

	}

	public void testBigtree() throws Exception{

		String tree = "(((((((t1:0.9803361397,t2:0.9035540882):0.0532383481,t3:0.2637392259):0.6273536528,(t4:0.8624112266,t5:0.3278892266):0.2606245542):0.2941323873,(t6:0.09820114588,t7:0.533115675):0.8625875909):0.7040311908,(((t8:0.8696136218,t9:0.08719484485):0.4204288905,(t10:0.102143287,(t11:0.9850614571,t12:0.7407912319):0.8715072596):0.5182644848):0.524062254,(((((((t13:0.3981794417,(t14:0.03889928572,t15:0.5187105467):0.1127638209):0.3431177251,((t16:0.4239511855,t17:0.001895790454):0.690600364,t18:0.6283850113):0.4073564562):0.6862231812,(((t19:0.9947085041,t20:0.4739363373):0.1873670686,t21:0.151270482):0.803061039,((t22:0.8899249982,((t23:0.1329096023,t24:0.84205155):0.8838408566,(t25:0.7541888549,t26:0.8602364615):0.8912267659):0.771449636):0.1022819551,(((t27:0.3134289116,(t28:0.2446750235,t29:0.8565168788):0.8277210968):0.4307989818,((t30:0.2330717787,t31:0.4438336496):0.6521712865,(t32:0.2534400895,t33:0.7885409284):0.3051449039):0.1196702593):0.4061951274,t34:0.8415271267):0.4365981282):0.753448925):0.1580670979):0.04210642632,(((t35:0.7504386581,t36:0.6328390085):0.9047614154,t37:0.4946133171):0.2264722914,((((t38:0.06683212146,t39:0.479845396):0.9424520086,t40:0.894530142):0.3844042511,(((t41:0.5215392481,t42:0.2366602973):0.8142298241,(t43:0.2968777204,(t44:0.655541793,t45:0.8608812049):0.3564132168):0.04912991729):0.1511388237,t46:0.9031036345):0.1874918914):0.9690212663,(t47:0.07753491728,(t48:0.8349514075,(t49:0.9689748741,t50:0.925813166):0.4534903264):0.3571097804):0.1324767114):0.5515443345):0.3330309158):0.7202291801,((t51:0.6977306763,((t52:0.9157640305,t53:0.4226291834):0.5872618856,t54:0.2063144948):0.1422286083):0.7182746637,t55:0.759545143):0.7437628019):0.2425582204,((t56:0.4614429038,(t57:0.9092229386,((t58:0.1049408391,t59:0.6328130178):0.642241966,((t60:0.264340204,t61:0.5904771155):0.7333205172,(t62:0.9183179205,t63:0.1090340314):0.3010568973):0.3240860389):0.3192155454):0.1835780439):0.5942421539,t64:0.7931551472):0.967891278):0.06263663713,(t65:0.5774453548,((t66:0.07208712469,((t67:0.8918803469,t68:0.5110983853):0.1491188321,t69:0.2471361952):0.9591872343):0.3133718621,(t70:0.944087367,t71:0.7830825299):0.2284035049):0.5492361034):0.1136150162):0.002181729767):0.4548798562):0.4258609388,((((((t72:0.27679418,t73:0.5398862793):0.8871422287,(((((t74:0.2531923286,t75:0.3796772889):0.4489221217,t76:0.2554209188):0.3248268673,t77:0.5372577759):0.5699883625,t78:0.1656995732):0.957750936,(t79:0.1301121258,t80:0.8925942327):0.2838441601):0.5258686764):0.47825964,(t81:0.5749240227,((t82:0.9574132746,(t83:0.00485483068,t84:0.8091488208):0.1985368489):0.3703975577,(((t85:0.3991035291,(t86:0.03201846033,t87:0.8380640063):0.05616304209):0.8414494572,t88:0.6844437125):0.2426782607,((t89:0.7543559887,t90:0.7162597755):0.8230077426,t91:0.08967904118):0.4460245941):0.8679371702):0.51572948):0.4362259945):0.2631344711,(((t92:0.3353162925,((t93:0.4025212794,t94:0.0281926766):0.7965471447,t95:0.1145715592):0.5993301494):0.08854756854,(t96:0.1461353719,((t97:0.3158547124,t98:0.06653800653):0.5634025722,t99:0.9711292514):0.9727503664):0.7684133062):0.4824229684,((t100:0.06834940333,t101:0.7794982188):0.3453287922,(t102:0.627945075,t103:0.1914187325):0.9974814849):0.6312927424):0.04858242651):0.2845227425,((t104:0.6782600286,(t105:0.03190574702,t106:0.5840284519):0.03041352634):0.725893975,(((t107:0.9885271091,t108:0.07126446022):0.8419693699,t109:0.1546431775):0.898004594,t110:0.2500803664):0.1493327522):0.4266726137):0.5946582041,(t111:0.1395377244,(((t112:0.7170655408,(t113:0.976886861,t114:0.9406369971):0.7471234254):0.8065501407,((t115:0.1713845057,(t116:0.7861330248,t117:0.6082276558):0.8413775554):0.3245444677,t118:0.3892389825):0.5992471091):0.7592411407,(((t119:0.535931844,t120:0.09058958571):0.4227561057,(t121:0.5531579193,t122:0.8276180199):0.6653355309):0.0941624688,t123:0.3623022255):0.1494971744):0.3526274569):0.9720881658):0.8149677955):0.6065687414,((((((t124:0.5406888947,t125:0.8892341822):0.06211395678,((t126:0.8203180477,(t127:0.8536844573,t128:0.360511546):0.9030223228):0.9095590916,((t129:0.9110714826,(t130:0.2346256471,t131:0.6523390864):0.1288849309):0.7077432328,(t132:0.4060195235,t133:0.1661393729):0.3910941551):0.205704404):0.8609933471):0.3724007562,((t134:0.1731842053,(t135:0.7232482471,(t136:0.3883952193,((t137:0.6709475764,t138:0.0372075201):0.5473196667,(t139:0.8092764446,t140:0.4123262055):0.2000603897):0.55258787):0.2654263263):0.745555162):0.2956101163,((t141:0.52147611,(t142:0.9462005703,t143:0.5671354234):0.6887917654):0.362258781,t144:0.4798202242):0.8242726682):0.6072624433):0.695287361,((((t145:0.03793937969,t146:0.07275558705):0.3482963489,t147:0.1457363514):0.1479936559,(t148:0.7158309214,((t149:0.2174433649,t150:0.04072828358):0.4112026501,t151:0.6422409331):0.3413406226):0.1693999742):0.6631712937,(((t152:0.2706006162,t153:0.9267972289):0.1387761638,((((t154:0.2563392594,t155:0.3058371837):0.5946117372,t156:0.6161190302):0.6970871226,(t157:0.2388902532,(t158:0.9486316761,t159:0.215360787):0.168830334):0.03888285463):0.1640696453,t160:0.6803096831):0.1418975852):0.4218000816,(((t161:0.8702562298,t162:0.9289729816):0.05807372741,t163:0.3533785399):0.5012762842,(((t164:0.8666574673,t165:0.9603798252):0.7887994377,t166:0.857058729):0.4139410679,(t167:0.5900272813,t168:0.3345388798):0.06017537019):0.9609203783):0.7103463742):0.696603697):0.6451920038):0.1909481271,((((t169:0.9171597108,t170:0.9479122513):0.7170342554,(t171:0.2722596873,((t172:0.1194724559,(t173:0.03922236571,t174:0.6290624789):0.07739861775):0.8598598302,(t175:0.2009421999,(t176:0.06154947914,t177:8.997193072E-4):0.04738179315):0.3235510678):0.3443877005):0.6351028818):0.5525081949,((((t178:0.7599076207,t179:0.2997759853):0.5921433992,t180:0.7098581635):0.3725496214,(t181:0.5053773888,(t182:0.5991492711,(t183:0.5036820578,t184:0.6361607853):0.510631816):0.9604382808):0.2464167587):0.6073093358,(((t185:0.03128415369,(t186:0.5260852403,(t187:0.878767435,t188:0.4992109234):0.5333148066):0.00347468094):0.5590308013,t189:0.3710992143):0.5034162949,(t190:0.778916508,((t191:0.3069154553,(((t192:0.9946115273,t193:0.9138687006):0.5209144899,t194:0.5152770842):0.9462409306,t195:0.7395236609):0.4110851623):0.930918345,(((t196:0.7895439987,((t197:0.4697002599,t198:0.1383787312):0.6911794308,(t199:0.8664436699,t200:0.1959039853):0.8656513852):0.3620497067):0.2839249384,(t201:0.6558795469,t202:0.2103423763):0.969477433):0.9058840063,(t203:0.0856692954,t204:0.4175976661):0.820434629):0.5355881769):0.2263581599):0.4512835185):0.7323478526):0.2479199937):0.1964542414,((t205:0.7537573762,(t206:0.1392466244,(t207:0.5136175761,(t208:0.7852529553,t209:0.07355738804):0.1220811389):0.7572090242):0.1422528555):0.5948274662,(((((t210:0.3068353184,(t211:0.3314456891,((t212:0.5265486804,t213:0.1382007354):0.1814086549,t214:0.9276472756):0.07718444197):0.03486835537):0.1617580003,(t215:0.3328830956,t216:0.8558843595):0.8366736979):0.347376487,t217:0.8222538356):0.2337225529,(t218:0.06199815008,t219:0.45975962):0.179990889):0.0635867205,(t220:0.3214025751,(t221:0.5022090652,t222:0.6454557138):0.6956466341):0.2711792416):0.1847200533):0.1051658324):0.4945860899):0.936143348,(((t223:0.06268779701,((t224:0.3337278806,t225:0.1570303424):0.3089733059,(t226:0.5069784883,t227:0.1434204187):0.2001587199):0.04750720505):0.3600859912,((((t228:0.9994731578,(t229:0.8934116936,t230:0.03698333143):0.8173468311):0.3089058488,((((t231:0.3216121283,t232:0.5232846253):0.8687884973,(t233:0.6280638413,((t234:0.6543256822,t235:0.8677638234):0.8895299246,t236:0.4047793006):0.7147388768):0.3533478715):0.9470084386,t237:0.7769409856):0.4955915695,((t238:0.2772087415,(t239:0.4904922615,(t240:0.05356206303,t241:0.08998329984):0.8154862223):0.5610961432):0.1617916438,(t242:0.5707751412,(t243:0.9836868793,t244:0.1984052949):0.6953297216):0.05552111682):0.9476150468):0.2473166997):0.9623488116,((t245:0.7935025664,t246:0.08509867964):0.3953444003,(t247:0.09163277131,(t248:0.5201428954,t249:0.8055520628):0.7452739514):0.3989078877):0.07581191277):0.9779064963,(((t250:0.943611098,(t251:0.33392801,t252:0.5996331484):0.4291575127):0.4906436009,((((t253:0.7749450852,(t254:0.8616885878,t255:0.585028409):0.06060880423):0.1238881133,((t256:0.7451687793,t257:0.6925335305):0.05338745634,t258:0.3357626374):0.2069296469):0.09644073155,((((t259:0.2258843291,t260:0.2671526412):0.3940743534,(t261:0.5022506947,(t262:0.9498897423,t263:0.1406114365):0.2847759123):0.04320593993):0.6982026948,t264:0.2693712024):0.959781138,(((t265:0.6035173486,t266:0.5529949202):0.9900399651,(t267:0.5455351078,t268:0.3530619899):0.4626278321):0.2735997427,(t269:0.9580646451,(t270:0.3280033092,t271:0.7206294278):0.03739526332):0.4967516926):0.9350089293):0.4371789068):0.1014483059,t272:0.2867298371):0.07522285799):0.06352435821,((t273:0.4001782183,t274:0.7190070178):0.1696753846,(t275:0.5535608665,t276:0.01324651297):0.2691543309):0.8676247413):0.8461736294):0.1769516913):0.344365149,(((t277:0.3245107541,(t278:0.4142541443,t279:0.5857141651):0.819547887):0.0867733527,(t280:0.4938162852,(t281:0.2444119717,t282:0.08141433029):0.05381231918):0.8375963389):0.176160393,((t283:0.4199601968,t284:0.8354801824):0.3150380594,(((t285:0.9818797186,(t286:0.8971825438,((t287:0.5155417006,t288:0.8260786769):0.7060374152,t289:0.6001661876):0.4120474763):0.9949228324):0.8038698458,t290:0.1939124272):0.6380942846,t291:0.3665255161):0.459349304):0.482901911):0.4833473735):0.5903116504):0.9973697898)";
		int nrTaxa = 291;
		int[] states = new int[291];
		Arrays.fill(states, 1);

		String orig = ".1";
		String locations = "" ;   //"1=1,2=0,3=0,4=1";
		for (int i=1; i<=states.length; i++){
			locations = locations + "t" + i + "=" + (states[i-1]-1) + (i<states.length?",":"");
		}


		String stateNumber = "2";
		String migrationMatrix = "0.2 0.1";
		String frequencies = "0.5 0.5";

		String R0 = Double.toString(4./3.) + " " + Double.toString(5.) ;
		String becomeUninfectiousRate = "1.5 1.25";
		String samplingProportion = Double.toString(1./3.) + " " + Double.toString(1./2.);


		int maxEvals = Integer.MAX_VALUE;
		double tolerance = 1e-13;

		double logL;

		conditionOnSurvival = true;

		logL = bdm_likelihood(tolerance, maxEvals, stateNumber,
				migrationMatrix,
				frequencies,
				tree, orig,
				R0,  null,
				becomeUninfectiousRate,
				samplingProportion,
				"t", locations, nrTaxa, null);

		System.out.println("Log-likelihood ("+(conditionOnSurvival?"":"not ")+"conditioned on survival) " + logL + "\t");


		assertEquals(-661.9588648301033, logL, 1e-5); // result from BEAST, not checked in R
	}

	// migration example adapted from BDSKY
	//
	// @Test
	public void testLikelihoodCalculationMig() throws Exception {

		String tree ="((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);"; //
		String orig="1."; //
		String stateNumber = "2";
		String migrationMatrix = ".2 .1";
		String frequencies = "0.5 0.5";

		String R0 = Double.toString(4./3.) + " " + Double.toString(5.) ;
		String becomeUninfectiousRate = "1.5 1.25";
		String samplingProportion = Double.toString(1./3.) + " " + Double.toString(1./2.);

		String locations = "1=1,2=0,3=0,4=1" ;

		int maxEvals = Integer.MAX_VALUE;
		double tolerance = 1e-14;

		double logL;

		//        conditionOnSurvival = true;
		//
		////         for (int i = 0; i<=10; i++){
		////
		////             R0 = Double.toString((2+i*.1)/(3./2.)) + " " + Double.toString((2+i*.1)/(3./2.)) ;
		//
		//        logL = bdm_likelihood(tolerance, maxEvals, stateNumber,
		//                migrationMatrix,
		//                frequencies,
		//                tree, orig,
		//                R0,  null,
		//                becomeUninfectiousRate,
		//                samplingProportion,
		//                "", locations, 4);
		//
		//        System.out.println("Log-likelihood (conditioned on survival) " + logL + "\t");
		//        assertEquals(-27.9720672, logL, 1e-4);

		//         }

		conditionOnSurvival = false;
		logL = bdm_likelihood(tolerance, maxEvals, stateNumber,
				migrationMatrix,
				frequencies,
				tree, orig,
				R0, null,
				becomeUninfectiousRate,
				samplingProportion,
				"", locations, 4, null);

		System.out.println("Log-likelihood (NOT conditioned on survival) = " + logL);
		assertEquals(-26.53293, logL, 1e-5);

	}


	// migration example adapted from BDSKY
	//
	// @Test
	public void testLikelihoodCalculationMigTiny() throws Exception {

		String tree ="(1 : 1.5, 2 : 0.5);"; //
		String orig="1."; //
		String stateNumber = "2";
		String migrationMatrix = ".1 .1";
		String frequencies = "0.5 0.5";

		// test without rate change
		//         String R0 = Double.toString(4./3.) + " 3";
		//         String becomeUninfectiousRate = "1.5 1.";
		//         String samplingProportion = Double.toString(1./3.) + "\t" + Double.toString(2./3.);


		String R0 = Double.toString(4./3.) + " " + Double.toString(4./3.) ;
		String becomeUninfectiousRate = "1.5 1.5";
		String samplingProportion = Double.toString(1./3.) + " " + Double.toString(1./3.);


		String locations = "1=0,2=1" ;

		int maxEvals = Integer.MAX_VALUE;
		double tolerance = 1e-10;

		double logL;

		conditionOnSurvival = false;

		//         for (int i = 0; i<=10; i++){
		//
		//             R0 = Double.toString((2+i*.1)/(3./2.)) + " " + Double.toString((2+i*.1)/(3./2.)) ;

		logL = bdm_likelihood(tolerance, maxEvals, stateNumber,
				migrationMatrix,
				frequencies,
				tree, orig,
				R0, null,
				becomeUninfectiousRate,
				samplingProportion,
				"",
				locations,
				2, null);

		//         System.out.println("Log-likelihood (conditioned on survival) = " + logL);
		System.out.println(logL + "\t");
		assertEquals(-7.215222, logL, 1e-6); // result from R


		String R0AmongDemes = "0.0666667 0.0666667" ;
		logL = bdm_likelihood(tolerance, maxEvals, stateNumber,
				"0 0",
				frequencies,
				tree, orig,
				R0, R0AmongDemes,
				becomeUninfectiousRate,
				samplingProportion,
				"",
				locations,
				2, null);

		//         System.out.println("Log-likelihood (conditioned on survival) = " + logL);
		System.out.println(logL + "\t");
		assertEquals(-7.404888, logL, 1e-6); // result from R

		R0AmongDemes = "0.0666667 0.1" ;
		logL = bdm_likelihood(tolerance, maxEvals, stateNumber,
				"0 0",
				frequencies,
				tree, orig,
				R0, R0AmongDemes,
				becomeUninfectiousRate,
				samplingProportion,
				"",
				locations,
				2, null);

		//         System.out.println("Log-likelihood (conditioned on survival) = " + logL);
		System.out.println(logL + "\t");
		assertEquals(-7.18723, logL, 1e-6); // result from R

		R0 = "2 1.3333333" ;
		logL = bdm_likelihood(tolerance, maxEvals, stateNumber,
				"0 0",
				frequencies,
				tree, orig,
				R0, R0AmongDemes,
				becomeUninfectiousRate,
				samplingProportion,
				"",
				locations,
				2, null);

		//         System.out.println("Log-likelihood (conditioned on survival) = " + logL);
		System.out.println(logL + "\t");
		assertEquals(-7.350649, logL, 1e-6); // result from R

		R0 = "2 1.5" ;
		becomeUninfectiousRate = "2 1" ;
		samplingProportion = "0.5 0.3";
		R0AmongDemes = "0.1 0.5" ;

		logL = bdm_likelihood(tolerance, maxEvals, stateNumber,
				"0 0",
				frequencies,
				tree, orig,
				R0, R0AmongDemes,
				becomeUninfectiousRate,
				samplingProportion,
				"",
				locations,
				2, null);

		//         System.out.println("Log-likelihood (conditioned on survival) = " + logL);
		System.out.println(logL + "\t");
		assertEquals(-6.504139, logL, 1e-6); // result from R

		locations = "1=1,2=0";
		logL = bdm_likelihood(tolerance, maxEvals, stateNumber,
				"0 0",
				frequencies,
				tree, orig,
				R0, R0AmongDemes,
				becomeUninfectiousRate,
				samplingProportion,
				"",
				locations,
				2, null);

		//         System.out.println("Log-likelihood (conditioned on survival) = " + logL);
		System.out.println(logL + "\t");
		assertEquals(-7.700916, logL, 1e-6); // result from R

		//         }

		//         conditionOnSurvival = false;
		//         logL = bdm_likelihood(tolerance, maxEvals, stateNumber,
		//                 migrationMatrix,
		//                 frequencies,
		//                 tree, orig,
		//                 R0,
		//                 becomeUninfectiousRate,
		//                 samplingProportion,
		//                 locations);
		//
		//         System.out.println("Log-likelihood (NOT conditioned on survival) = " + logL);

	}
	//    public void testGeorge() throws Exception {
	//
	////        String tree = "(((\"t4_2\"[&state=2]:0.6132389882,(\"t5_2\"[&state=2]:1.133526122,\"t3_2\"[&state=2]:0.7334865664):0.2213778867):0.00159031537,\"t2_1\"[&state=1]:0.1394206244):0.06661565173,\"t1_1\"[&state=1]:0.7815374183);";
	//
	//        String tree = "(1[&state=1]:0.02439127831,2[&state=2]:0.4577793404);";
	//        String origin = "0.4015719006";
	//
	//        Double birth1 = 1.0776298 ;
	//        Double birth2 = 0.9443746 ;
	//
	//        Double death1 = 1.8010283;
	//        Double death2 = 1.2499457;
	//
	//        Double s1 = 0.5 ;
	//        Double s2 = 0.5 ;
	//
	//        String locations = "1=0,2=0, 3=1, 4=1, 5=1" ;
	//
	//        int maxEvals = Integer.MAX_VALUE;
	//        double tolerance = 1e-10;
	//
	//        double logL;
	//
	//
	////        for (double lambda1=0; lambda1<=2.2; lambda1+=0.1){
	//
	//        logL = bdm_likelihood(tolerance, maxEvals, "2",
	//                "0.4268963 1.1652556", //migmatrix
	//                "0.5 0.5", // frequencies
	//                tree, origin,
	//                Double.toString(birth1 / death1) + " " + Double.toString(birth2/death2),
	//                null,
	//                Double.toString(death1) + " " + Double.toString(death2),
	//                Double.toString(s1*death1) + " " + Double.toString(s2*death2),
	//                "", locations, 2, null);
	//
	////            System.out.println(lambda1 + "\t" + logL);
	//        System.out.println( logL);
	//        //        }
	//
	////   lam11        lam12        lam21        lam22           d1           d2           s1           s2      gamma12      gamma21 probontrans1 probontrans2
	////   1.0776298    0.0000000    0.0000000    0.9443746    1.8010283    1.2499457    0.5000000    0.5000000    0.4268963    1.1652556    0.0000000    0.0000000
	//
	//    }

	public double bdm_likelihood(double tolerance, int maxEvals, String statenumber, String migrationMatrix,
			String frequencies, String newick, String origin,
			String R0, String R0AmongDemes, String becomeUninfectiousRate, String samplingProportion, String prefixname, String locations, int nrTaxa, String intervalTimes) throws Exception {

		return bdm_likelihood(tolerance,maxEvals, statenumber,migrationMatrix,
				frequencies, newick, origin,
				R0, R0AmongDemes, becomeUninfectiousRate, samplingProportion, prefixname, locations, nrTaxa, intervalTimes, false);
	}

	public double bdm_likelihood(double tolerance, int maxEvals, String statenumber, String migrationMatrix,
			String frequencies, String newick, String origin,
			String R0, String R0AmongDemes, String becomeUninfectiousRate, String samplingProportion, String prefixname, String locations, int nrTaxa, String intervalTimes, Boolean SA) throws Exception {



		BirthDeathMigrationModelUncoloured bdm =  new BirthDeathMigrationModelUncoloured();



		ArrayList<Taxon> taxa = new ArrayList<Taxon>();

		for (int i=1; i<=nrTaxa; i++){
			taxa.add(new Taxon(prefixname+i));
		}

		Tree tree = new TreeParser();
		//        tree.setInputValue("singlechild", "true");
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

		bdm.setInputValue("tree", tree);
		bdm.setInputValue("tiptypes", trait);

		bdm.setInputValue("origin", Double.toString(Double.parseDouble(origin)+tree.getRoot().getHeight()));
		bdm.setInputValue("stateNumber", statenumber);
		bdm.setInputValue("migrationMatrix", migrationMatrix);
		bdm.setInputValue("frequencies", frequencies);
		bdm.setInputValue("checkRho", false);

		bdm.setInputValue("R0", R0);

		//        String R0AmongDemes = Double.toString(3./6.) + " " + Double.toString(2.5) ;
		if (R0AmongDemes != null) bdm.setInputValue("R0AmongDemes", R0AmongDemes);

		bdm.setInputValue("becomeUninfectiousRate", becomeUninfectiousRate);
		bdm.setInputValue("samplingProportion", samplingProportion);
		bdm.setInputValue("maxEvaluations", maxEvals);
		bdm.setInputValue("intervalTimes", intervalTimes);

		bdm.setInputValue("conditionOnSurvival", conditionOnSurvival);

		//TO DO uncomment
		//bdm.setInputValue("tolerance", tolerance);

		if (SA) {
			Double[] r = new Double[Integer.parseInt(statenumber)];
			Arrays.fill(r, 1.);
			bdm.setInputValue("removalProbability", new RealParameter(r));
		}

		bdm.initAndValidate();



		long startTime = System.currentTimeMillis();
		double logL = bdm.calculateLogP(); //calculateTreeLogLikelihood(coltree.getUncolouredTree());
		runtime = System.currentTimeMillis() - startTime;
		maxEvalsUsed = bdm.maxEvalsUsed;

		assertEquals(0., 0., 1e-2); // TO DO REMOVE THIS LINE if no reason for it being here is found

		return logL;

	}

	public void testLikelihoodCalculation4() throws Exception {

		BirthDeathMigrationModelUncoloured bdm =  new BirthDeathMigrationModelUncoloured();

		Tree tree = new TreeParser("((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);",false);
		bdm.setInputValue("tree", tree);

		String prefixname = "";
		ArrayList<Taxon> taxa = new ArrayList<Taxon>();
		for (int i=1; i<=4; i++){
			taxa.add(new Taxon(prefixname+i));
		}

		TraitSet trait = new TraitSet();
		trait.setInputValue("taxa", new TaxonSet(taxa));
		trait.setInputValue("value", "1=0,2=0, 3=0, 4=0");
		trait.setInputValue("traitname", "tiptypes");
		trait.initAndValidate();

		bdm.setInputValue("tree", tree);
		bdm.setInputValue("tiptypes", trait);


		bdm.setInputValue("origin", new RealParameter("6.")); // = treeheight + 1
		bdm.setInputValue("stateNumber", "1");
		bdm.setInputValue("migrationMatrix", "0");
		bdm.setInputValue("frequencies", "1");

		//        bdm.setInputValue("intervalNumber", 3);
		bdm.setInputValue("intervalTimes", new RealParameter("0. 3. 4.5"));

		bdm.setInputValue("R0", new RealParameter(new Double[]{2./3., 4./3., 8./3.}));
		bdm.setInputValue("becomeUninfectiousRate", new RealParameter("4.5 1.5 1.5"));
		bdm.setInputValue("samplingProportion", new RealParameter(new Double[]{4./9., 1./3., 2./3.}));

		bdm.initAndValidate();
		bdm.setInputValue("conditionOnSurvival", false);
		double logP = bdm.calculateLogP();

		assertEquals(-37.8056,logP, 1e-2);
		System.out.println("NOT conditioned on survival: " + logP);

		bdm.setInputValue("conditionOnSurvival", true);
		bdm.initAndValidate();


		logP = bdm.calculateLogP();
		assertEquals(-37.30123947394569, logP, 1e-2);
		System.out.println("conditioned on survival: " + logP);

	}

	//    public void testInfAMongDemesTanjasLatvianTestTree() throws Exception{
	//
	//        String tree = "(t142:0.183060207,(t198:0.01180228928,((((((((((t100:0.1634185983,t155:0.176845951):0.06125685224,(t90:0.2391258009,((t14:0.1688367367,t62:0.1382429133):0.02139973524,((t152:0.03832952194,t88:0.04748321911):0.06995234823,t40:0.0986303732):0.08495120652):0.08095905856):0.024090675):0.05502436057,(t9:0.016875208,(t136:0.1599505763,t115:0.1851914767):0.02114991916):0.1157982687):0.02505648067,((t191:0.2238218232,t32:0.2388083634):0.009103003537,t66:0.162597966):0.1403399862):0.01137570248,(t173:0.01281029964,((t106:0.151507325,(t16:0.1367418229,((t181:0.01690190105,t83:0.03139518674):0.1220090972,((t165:0.02348911079,t122:0.05928644408):0.09983758416,t166:0.103003683):0.001320521952):0.05404212282):0.1309743334):0.01165354228,((t144:0.02192190953,((((((t45:0.151904261,t127:0.1228277428):0.01751402187,t58:0.07155311902):0.03119543921,(t153:0.1643371305,(t91:0.1356500251,t6:0.07584656203):0.03523977934):0.002781987005):0.07666044116,(t11:0.03982764634,t71:0.248150788):0.02725485004):0.02437437474,t196:0.2944340121):0.004226398195,((t33:0.06976502685,t30:0.07359772989):0.04680491651,t180:0.03113827907):0.1504993392):0.02464072669):0.02582509488,(t59:0.04786010356,(t46:0.088783985,(t111:6.758411617E-4,t51:0.008484605266):0.03406107085):0.08959686263):0.08964185833):0.005347847047):7.839154956E-4):0.04246678371):0.006323412459,(((t148:0.2616973451,(t25:0.1485161003,t171:0.07110687729):0.1064329955):0.0207787381,t55:0.03068826073):0.0307313671,((((t84:0.05622204425,((t120:0.006585266951,t50:0.01291705632):0.0356409028,t26:0.005108493439):0.0142455076):0.1461389222,(t12:0.1882047365,t151:0.1693527171):0.02716535138):0.08346882642,t156:0.01134880599):0.05232911647,((t108:0.01025212681,t131:0.1257695324):0.1437260925,((t15:0.08139845545,t60:0.01716768662):0.05235158123,t21:0.03437918285):0.1779913884):0.04800565135):0.01589697111):0.0295971305):0.06365059958,(t121:0.0570598473,(((t101:0.003757120065,(t73:0.1331208962,t104:0.1390293964):0.03226618599):0.07895021643,t70:0.2029327018):0.03232810914,(t42:0.0912770908,t170:0.1009290436):0.1576475743):0.1702905689):1.49960922E-4):0.02622749623,(t190:0.3863223765,((((t20:0.006144675272,(t17:0.2076446847,t7:0.186814124):0.04757474124):0.05763698983,((((t195:0.05793824162,t184:0.06303464549):0.08581589396,t56:0.02170833407):0.1289768065,(t49:0.1520387747,(t116:0.01868664063,(t75:0.1386339963,(t102:0.09163444097,t178:0.01597096436):0.06122327271):0.0142621277):0.04638068884):0.05837448127):0.006005043284,(t23:0.1356035764,(t124:0.1095003806,t94:0.01166115622):0.005174865347):0.1123754984):0.08549847623):0.03848203141,(t133:0.064454232,(t80:0.1426635889,(t29:0.003532335453,((t36:0.08427990558,t197:0.02731956008):0.02526265024,t93:0.06034347072):0.05067758546):0.0842293226):0.1037653608):0.04149527346):0.001087278236,((((((((t105:0.05418593228,t5:0.05347499934):0.1092575268,(t87:0.05126456021,((t172:0.08901087761,(t193:0.1104958027,t63:0.04173081124):0.01167865165):0.04620158527,t76:0.05538311904):0.0260737547):0.03139896852):0.04302391868,((t168:0.01747855094,t182:0.004321316462):0.1154562131,(t176:0.2036736767,t22:0.03936042753):0.01222087371):0.03189806639):0.03843343927,(t98:0.3013555011,t154:0.2070581369):0.009036684462):0.02478079121,t24:0.02846779585):0.01158443985,t47:0.1386745288):0.04022715139,(((t126:0.06431536312,t140:0.03118945665):0.05312353763,t103:0.08847651209):0.07393791882,(t174:0.1347066037,(((t129:0.04628905913,t199:0.01558097455):0.01286585011,t187:0.03053488837):0.02995162052,t157:0.06532032986):0.04301296567):0.05240200038):0.1838428476):0.01762093497,(t77:0.05231337331,(t138:0.1683801303,((t86:0.0728897394,t119:0.05609098346):0.1437950228,t125:0.1924713237):0.1131018254):0.04496657346):0.01944226675):0.004583006969):0.05147142928):0.03956163653):8.491132634E-4,(((((((t163:0.1201147561,t44:0.08285147313):0.1228651137,((t194:0.1700679131,(((t68:0.1751727346,((t78:0.1008519948,(t82:0.06423817278,t145:0.03395378342):0.006710614292):0.03493791633,t54:0.08951110861):0.0298547759):0.02673622034,t169:0.05622630238):2.269156301E-4,((t43:0.09187936242,(t123:0.07939197877,t159:0.00885262305):0.09097112099):0.01686395903,t200:0.112520717):0.00410229035):0.0392765129):0.01044584409,(((t69:0.0993924397,t99:0.03716855707):0.1057649837,t72:0.2466977489):6.537301731E-4,(((t183:0.1221756465,(t8:0.04600100446,t185:0.1598975905):0.01236092288):0.03134936472,t64:0.2331719844):0.005748850337,(t35:0.09956972948,t95:0.1455439657):0.08297895819):0.01833473632):0.008558418572):0.004363695617):0.04776935424,(t48:0.1838206114,t162:0.04042607718):0.06426449193):0.01183601753,t107:0.05148569515):0.1154398066,t1:0.0528720883):0.01094055235,((t164:0.08798583573,t192:0.1655549963):0.2643199004,(((t188:0.07537181053,(t10:0.09000280238,t118:0.06889676539):0.0434599342):0.01906353948,t74:0.1521029803):0.116143132,(((((t177:0.1221571447,t158:0.182846772):0.03050588246,t161:0.174193549):0.02770764754,t113:0.1746231004):0.00467812037,t117:0.1940242414):0.0877704612,((t160:0.0195669223,((t141:0.02093340124,t167:0.1178564624):0.08021994961,(t189:0.2029063144,t79:0.04413767462):0.01343368462):0.1080166967):0.006278860901,((t67:0.02758579073,(t135:0.211053921,t175:0.1493736704):0.03223408488):0.02812430545,(t41:0.1852120201,(t65:0.10096633,t39:0.06974880375):0.0962984252):0.06948871982):0.0500274026):0.006290958099):0.002201313325):0.08829518442):0.02973761663):0.03205348426,((((t96:0.02855231257,t57:0.2366603091):0.08703694777,t128:0.2909812592):0.1068902459,((((t114:0.03874104846,t147:0.07489950667):0.01264666222,t38:0.1080981363):0.1915099756,t2:0.07914537186):0.1364513481,(t186:0.2896596271,((t27:0.2661713312,((t61:0.02629957646,t85:0.067533441):0.06356071446,(t150:0.1117478683,t97:0.1494987117):0.002895683023):0.09962172567):0.02035645508,((t19:0.04944628219,(t37:0.04422596834,t110:0.1635019172):0.006321885964):0.0986250419,(t139:0.09913687139,((((t28:0.05041488701,(t149:0.07203828491,t31:0.04965796547):0.0348763622):0.0174977869,t179:0.0432033841):0.02049505061,t146:0.1145895663):0.002492846607,(t52:0.06280646817,t89:0.01153453577):0.0872711871):0.1243045497):0.01218030615):0.001171646263):0.01605661203):0.1364112954):0.01345028016):0.03194573119,((t134:0.2725510659,(t143:0.06613600244,t13:0.05002345453):0.08925245306):0.06279893606,(((t130:0.2396989734,t137:0.02336619634):0.04278756619,(t18:0.08261003708,t81:0.1645683621):0.1034035373):0.01146696761,t132:0.08196130812):0.03957785721):0.1156576885):0.003568506952):0.01052412183):0.04825274753,((t3:0.2170329863,t4:0.04295718708):0.1287823831,(((t53:0.08968123545,t34:0.07908222397):0.005881300452,t92:0.04965904536):0.08001562771,(t109:0.02032432217,t112:0.05469748807):0.05596055895):0.1690398883):0.2018440042):0.05643582614):0.02748407343);";
	//        int nrTaxa = 200;
	//
	//        String orig = "0.136171103";
	//
	//        String locations = "t198=0, t142=1, t121=0, t1=0, t77=0, t173=0, t55=0, t47=0, t133=0, t24=0, t107=0, t144=0, t160=0, t156=0, t20=0, t2=0, t132=0, t137=0, t139=0, t59=0, t17=0, t162=0, t67=0, t9=0, t11=0, t96=1, t143=0, t13=0, t4=0, t169=0, t106=0, t101=0, t22=0, t154=0, t8=0, t79=0, t188=0, t37=0, t29=1, t138=1, t70=0, t116=0, t108=0, t76=0, t87=0, t6=0, t94=0, t19=0, t171=1, t80=0, t99=1, t194=1, t43=0, t44=1, t56=1, t111=0, t51=0, t141=0, t18=0, t14=1, t103=1, t93=0, t180=0, t58=0, t115=1, t168=0, t182=0, t179=0, t161=0, t200=0, t64=1, t21=0, t63=0, t155=0, t113=1, t109=1, t164=1, t16=0, t131=0, t100=0, t157=0, t48=0, t159=0, t61=0, t177=0, t147=0, t105=0, t118=0, t46=0, t178=0, t114=0, t185=0, t85=1, t69=0, t40=0, t128=1, t145=0, t66=0, t175=0, t28=0, t190=0, t74=1, t192=0, t60=0, t35=0, t49=0, t183=0, t7=0, t119=0, t112=1, t25=0, t197=0, t10=0, t187=0, t166=1, t82=1, t126=0, t148=0, t39=0, t54=0, t136=1, t123=0, t5=0, t172=1, t149=0, t26=0, t65=0, t153=0, t117=0, t15=0, t199=0, t89=0, t81=0, t140=0, t124=0, t102=0, t125=0, t42=0, t34=0, t176=0, t104=0, t150=0, t33=0, t146=0, t97=0, t195=0, t68=0, t90=0, t78=0, t174=0, t196=0, t72=0, t186=0, t88=0, t84=0, t92=0, t189=0, t62=0, t152=0, t170=0, t41=0, t38=0, t12=1, t134=0, t130=0, t36=0, t30=0, t167=0, t181=0, t83=0, t23=0, t3=0, t165=0, t73=0, t91=0, t135=0, t31=0, t50=0, t151=0, t163=0, t129=0, t127=0, t75=0, t120=0, t86=0, t53=0, t57=0, t122=0, t110=0, t191=1, t95=1, t32=0, t52=0, t27=0, t193=0, t98=0, t184=1, t158=0, t71=1, t45=0" ;
	//
	//        String stateNumber = "2";
	//        String migrationMatrix = "0. 0.";
	//        String frequencies = "0.5 0.5";
	//
	//        String R0 = Double.toString(15./6.) + " 3.5" ;
	//        String becomeUninfectiousRate = "6. 2.";
	//        String samplingProportion = "0.25 0.25";
	//
	//
	//        int maxEvals = 10000; //Integer.MAX_VALUE;
	//        double tolerance = 1e-12;
	//
	//        double logL;
	//
	//        conditionOnSurvival = false;
	//
	//        logL = bdm_likelihood(tolerance, maxEvals, stateNumber,
	//                migrationMatrix,
	//                frequencies,
	//                tree, orig,
	//                R0,
	//                Double.toString(3./6.) + " " + Double.toString(2.5),
	//                becomeUninfectiousRate,
	//                samplingProportion,
	//                "t", locations, nrTaxa, null);
	//
	//        System.out.println("Log-likelihood ("+(conditionOnSurvival?"":"not")+" conditioned on survival) " + logL + "\t");
	//        assertEquals(-589.232, logL, 1e-3);
	//
	//    }
	
	// TO DO REMOVE COMMENT
	// THIS TEST WAS COMMENTED OUT UNTIL NOW 23/02/17
//	    public void testInfAMongDemesTanjasLatvianTestTree_bds() throws Exception{
//	
//	        String tree = "(t142:0.183060207,(t198:0.01180228928,((((((((((t100:0.1634185983,t155:0.176845951):0.06125685224,(t90:0.2391258009,((t14:0.1688367367,t62:0.1382429133):0.02139973524,((t152:0.03832952194,t88:0.04748321911):0.06995234823,t40:0.0986303732):0.08495120652):0.08095905856):0.024090675):0.05502436057,(t9:0.016875208,(t136:0.1599505763,t115:0.1851914767):0.02114991916):0.1157982687):0.02505648067,((t191:0.2238218232,t32:0.2388083634):0.009103003537,t66:0.162597966):0.1403399862):0.01137570248,(t173:0.01281029964,((t106:0.151507325,(t16:0.1367418229,((t181:0.01690190105,t83:0.03139518674):0.1220090972,((t165:0.02348911079,t122:0.05928644408):0.09983758416,t166:0.103003683):0.001320521952):0.05404212282):0.1309743334):0.01165354228,((t144:0.02192190953,((((((t45:0.151904261,t127:0.1228277428):0.01751402187,t58:0.07155311902):0.03119543921,(t153:0.1643371305,(t91:0.1356500251,t6:0.07584656203):0.03523977934):0.002781987005):0.07666044116,(t11:0.03982764634,t71:0.248150788):0.02725485004):0.02437437474,t196:0.2944340121):0.004226398195,((t33:0.06976502685,t30:0.07359772989):0.04680491651,t180:0.03113827907):0.1504993392):0.02464072669):0.02582509488,(t59:0.04786010356,(t46:0.088783985,(t111:6.758411617E-4,t51:0.008484605266):0.03406107085):0.08959686263):0.08964185833):0.005347847047):7.839154956E-4):0.04246678371):0.006323412459,(((t148:0.2616973451,(t25:0.1485161003,t171:0.07110687729):0.1064329955):0.0207787381,t55:0.03068826073):0.0307313671,((((t84:0.05622204425,((t120:0.006585266951,t50:0.01291705632):0.0356409028,t26:0.005108493439):0.0142455076):0.1461389222,(t12:0.1882047365,t151:0.1693527171):0.02716535138):0.08346882642,t156:0.01134880599):0.05232911647,((t108:0.01025212681,t131:0.1257695324):0.1437260925,((t15:0.08139845545,t60:0.01716768662):0.05235158123,t21:0.03437918285):0.1779913884):0.04800565135):0.01589697111):0.0295971305):0.06365059958,(t121:0.0570598473,(((t101:0.003757120065,(t73:0.1331208962,t104:0.1390293964):0.03226618599):0.07895021643,t70:0.2029327018):0.03232810914,(t42:0.0912770908,t170:0.1009290436):0.1576475743):0.1702905689):1.49960922E-4):0.02622749623,(t190:0.3863223765,((((t20:0.006144675272,(t17:0.2076446847,t7:0.186814124):0.04757474124):0.05763698983,((((t195:0.05793824162,t184:0.06303464549):0.08581589396,t56:0.02170833407):0.1289768065,(t49:0.1520387747,(t116:0.01868664063,(t75:0.1386339963,(t102:0.09163444097,t178:0.01597096436):0.06122327271):0.0142621277):0.04638068884):0.05837448127):0.006005043284,(t23:0.1356035764,(t124:0.1095003806,t94:0.01166115622):0.005174865347):0.1123754984):0.08549847623):0.03848203141,(t133:0.064454232,(t80:0.1426635889,(t29:0.003532335453,((t36:0.08427990558,t197:0.02731956008):0.02526265024,t93:0.06034347072):0.05067758546):0.0842293226):0.1037653608):0.04149527346):0.001087278236,((((((((t105:0.05418593228,t5:0.05347499934):0.1092575268,(t87:0.05126456021,((t172:0.08901087761,(t193:0.1104958027,t63:0.04173081124):0.01167865165):0.04620158527,t76:0.05538311904):0.0260737547):0.03139896852):0.04302391868,((t168:0.01747855094,t182:0.004321316462):0.1154562131,(t176:0.2036736767,t22:0.03936042753):0.01222087371):0.03189806639):0.03843343927,(t98:0.3013555011,t154:0.2070581369):0.009036684462):0.02478079121,t24:0.02846779585):0.01158443985,t47:0.1386745288):0.04022715139,(((t126:0.06431536312,t140:0.03118945665):0.05312353763,t103:0.08847651209):0.07393791882,(t174:0.1347066037,(((t129:0.04628905913,t199:0.01558097455):0.01286585011,t187:0.03053488837):0.02995162052,t157:0.06532032986):0.04301296567):0.05240200038):0.1838428476):0.01762093497,(t77:0.05231337331,(t138:0.1683801303,((t86:0.0728897394,t119:0.05609098346):0.1437950228,t125:0.1924713237):0.1131018254):0.04496657346):0.01944226675):0.004583006969):0.05147142928):0.03956163653):8.491132634E-4,(((((((t163:0.1201147561,t44:0.08285147313):0.1228651137,((t194:0.1700679131,(((t68:0.1751727346,((t78:0.1008519948,(t82:0.06423817278,t145:0.03395378342):0.006710614292):0.03493791633,t54:0.08951110861):0.0298547759):0.02673622034,t169:0.05622630238):2.269156301E-4,((t43:0.09187936242,(t123:0.07939197877,t159:0.00885262305):0.09097112099):0.01686395903,t200:0.112520717):0.00410229035):0.0392765129):0.01044584409,(((t69:0.0993924397,t99:0.03716855707):0.1057649837,t72:0.2466977489):6.537301731E-4,(((t183:0.1221756465,(t8:0.04600100446,t185:0.1598975905):0.01236092288):0.03134936472,t64:0.2331719844):0.005748850337,(t35:0.09956972948,t95:0.1455439657):0.08297895819):0.01833473632):0.008558418572):0.004363695617):0.04776935424,(t48:0.1838206114,t162:0.04042607718):0.06426449193):0.01183601753,t107:0.05148569515):0.1154398066,t1:0.0528720883):0.01094055235,((t164:0.08798583573,t192:0.1655549963):0.2643199004,(((t188:0.07537181053,(t10:0.09000280238,t118:0.06889676539):0.0434599342):0.01906353948,t74:0.1521029803):0.116143132,(((((t177:0.1221571447,t158:0.182846772):0.03050588246,t161:0.174193549):0.02770764754,t113:0.1746231004):0.00467812037,t117:0.1940242414):0.0877704612,((t160:0.0195669223,((t141:0.02093340124,t167:0.1178564624):0.08021994961,(t189:0.2029063144,t79:0.04413767462):0.01343368462):0.1080166967):0.006278860901,((t67:0.02758579073,(t135:0.211053921,t175:0.1493736704):0.03223408488):0.02812430545,(t41:0.1852120201,(t65:0.10096633,t39:0.06974880375):0.0962984252):0.06948871982):0.0500274026):0.006290958099):0.002201313325):0.08829518442):0.02973761663):0.03205348426,((((t96:0.02855231257,t57:0.2366603091):0.08703694777,t128:0.2909812592):0.1068902459,((((t114:0.03874104846,t147:0.07489950667):0.01264666222,t38:0.1080981363):0.1915099756,t2:0.07914537186):0.1364513481,(t186:0.2896596271,((t27:0.2661713312,((t61:0.02629957646,t85:0.067533441):0.06356071446,(t150:0.1117478683,t97:0.1494987117):0.002895683023):0.09962172567):0.02035645508,((t19:0.04944628219,(t37:0.04422596834,t110:0.1635019172):0.006321885964):0.0986250419,(t139:0.09913687139,((((t28:0.05041488701,(t149:0.07203828491,t31:0.04965796547):0.0348763622):0.0174977869,t179:0.0432033841):0.02049505061,t146:0.1145895663):0.002492846607,(t52:0.06280646817,t89:0.01153453577):0.0872711871):0.1243045497):0.01218030615):0.001171646263):0.01605661203):0.1364112954):0.01345028016):0.03194573119,((t134:0.2725510659,(t143:0.06613600244,t13:0.05002345453):0.08925245306):0.06279893606,(((t130:0.2396989734,t137:0.02336619634):0.04278756619,(t18:0.08261003708,t81:0.1645683621):0.1034035373):0.01146696761,t132:0.08196130812):0.03957785721):0.1156576885):0.003568506952):0.01052412183):0.04825274753,((t3:0.2170329863,t4:0.04295718708):0.1287823831,(((t53:0.08968123545,t34:0.07908222397):0.005881300452,t92:0.04965904536):0.08001562771,(t109:0.02032432217,t112:0.05469748807):0.05596055895):0.1690398883):0.2018440042):0.05643582614):0.02748407343);";
//	        int nrTaxa = 200;
//	
//	        String orig = "0.136171103";
//	
//	        String locations = "t198=0, t142=1, t121=0, t1=0, t77=0, t173=0, t55=0, t47=0, t133=0, t24=0, t107=0, t144=0, t160=0, t156=0, t20=0, t2=0, t132=0, t137=0, t139=0, t59=0, t17=0, t162=0, t67=0, t9=0, t11=0, t96=1, t143=0, t13=0, t4=0, t169=0, t106=0, t101=0, t22=0, t154=0, t8=0, t79=0, t188=0, t37=0, t29=1, t138=1, t70=0, t116=0, t108=0, t76=0, t87=0, t6=0, t94=0, t19=0, t171=1, t80=0, t99=1, t194=1, t43=0, t44=1, t56=1, t111=0, t51=0, t141=0, t18=0, t14=1, t103=1, t93=0, t180=0, t58=0, t115=1, t168=0, t182=0, t179=0, t161=0, t200=0, t64=1, t21=0, t63=0, t155=0, t113=1, t109=1, t164=1, t16=0, t131=0, t100=0, t157=0, t48=0, t159=0, t61=0, t177=0, t147=0, t105=0, t118=0, t46=0, t178=0, t114=0, t185=0, t85=1, t69=0, t40=0, t128=1, t145=0, t66=0, t175=0, t28=0, t190=0, t74=1, t192=0, t60=0, t35=0, t49=0, t183=0, t7=0, t119=0, t112=1, t25=0, t197=0, t10=0, t187=0, t166=1, t82=1, t126=0, t148=0, t39=0, t54=0, t136=1, t123=0, t5=0, t172=1, t149=0, t26=0, t65=0, t153=0, t117=0, t15=0, t199=0, t89=0, t81=0, t140=0, t124=0, t102=0, t125=0, t42=0, t34=0, t176=0, t104=0, t150=0, t33=0, t146=0, t97=0, t195=0, t68=0, t90=0, t78=0, t174=0, t196=0, t72=0, t186=0, t88=0, t84=0, t92=0, t189=0, t62=0, t152=0, t170=0, t41=0, t38=0, t12=1, t134=0, t130=0, t36=0, t30=0, t167=0, t181=0, t83=0, t23=0, t3=0, t165=0, t73=0, t91=0, t135=0, t31=0, t50=0, t151=0, t163=0, t129=0, t127=0, t75=0, t120=0, t86=0, t53=0, t57=0, t122=0, t110=0, t191=1, t95=1, t32=0, t52=0, t27=0, t193=0, t98=0, t184=1, t158=0, t71=1, t45=0" ;
//	
//	        String stateNumber = "2";
//	        String migrationMatrix = "0. 0.";
//	        String frequencies = "0.5 0.5";
//	
//	        String b = "15. 7." ;
//	        String b_ij = "3. 5." ;
//	        String d = "4.5 1.5";
//	        String s = "1.5 0.5";
//	
//	
//	        int maxEvals = 10000; //Integer.MAX_VALUE;
//	        double tolerance = 1e-12;
//	
//	        double logL;
//	
//	        conditionOnSurvival = false;
//	
//	        logL = bdm_likelihood_bds(tolerance, maxEvals, stateNumber,
//	                migrationMatrix,
//	                frequencies,
//	                tree, orig,
//	                b,
//	                b_ij,
//	                d,
//	                s,
//	                "t", locations, nrTaxa);
//	        
//	
//	        System.out.println("Log-likelihood ("+(conditionOnSurvival?"":"not")+" conditioned on survival) " + logL + "\t");
//	        assertEquals(-589.232, logL, 1e-3);
//	
//	    }
//

	//    public void testCorrelations() throws Exception{
	//
	//
	//        String tree = "(t142:0.183060207,(t198:0.01180228928,((((((((((t100:0.1634185983,t155:0.176845951):0.06125685224,(t90:0.2391258009,((t14:0.1688367367,t62:0.1382429133):0.02139973524,((t152:0.03832952194,t88:0.04748321911):0.06995234823,t40:0.0986303732):0.08495120652):0.08095905856):0.024090675):0.05502436057,(t9:0.016875208,(t136:0.1599505763,t115:0.1851914767):0.02114991916):0.1157982687):0.02505648067,((t191:0.2238218232,t32:0.2388083634):0.009103003537,t66:0.162597966):0.1403399862):0.01137570248,(t173:0.01281029964,((t106:0.151507325,(t16:0.1367418229,((t181:0.01690190105,t83:0.03139518674):0.1220090972,((t165:0.02348911079,t122:0.05928644408):0.09983758416,t166:0.103003683):0.001320521952):0.05404212282):0.1309743334):0.01165354228,((t144:0.02192190953,((((((t45:0.151904261,t127:0.1228277428):0.01751402187,t58:0.07155311902):0.03119543921,(t153:0.1643371305,(t91:0.1356500251,t6:0.07584656203):0.03523977934):0.002781987005):0.07666044116,(t11:0.03982764634,t71:0.248150788):0.02725485004):0.02437437474,t196:0.2944340121):0.004226398195,((t33:0.06976502685,t30:0.07359772989):0.04680491651,t180:0.03113827907):0.1504993392):0.02464072669):0.02582509488,(t59:0.04786010356,(t46:0.088783985,(t111:6.758411617E-4,t51:0.008484605266):0.03406107085):0.08959686263):0.08964185833):0.005347847047):7.839154956E-4):0.04246678371):0.006323412459,(((t148:0.2616973451,(t25:0.1485161003,t171:0.07110687729):0.1064329955):0.0207787381,t55:0.03068826073):0.0307313671,((((t84:0.05622204425,((t120:0.006585266951,t50:0.01291705632):0.0356409028,t26:0.005108493439):0.0142455076):0.1461389222,(t12:0.1882047365,t151:0.1693527171):0.02716535138):0.08346882642,t156:0.01134880599):0.05232911647,((t108:0.01025212681,t131:0.1257695324):0.1437260925,((t15:0.08139845545,t60:0.01716768662):0.05235158123,t21:0.03437918285):0.1779913884):0.04800565135):0.01589697111):0.0295971305):0.06365059958,(t121:0.0570598473,(((t101:0.003757120065,(t73:0.1331208962,t104:0.1390293964):0.03226618599):0.07895021643,t70:0.2029327018):0.03232810914,(t42:0.0912770908,t170:0.1009290436):0.1576475743):0.1702905689):1.49960922E-4):0.02622749623,(t190:0.3863223765,((((t20:0.006144675272,(t17:0.2076446847,t7:0.186814124):0.04757474124):0.05763698983,((((t195:0.05793824162,t184:0.06303464549):0.08581589396,t56:0.02170833407):0.1289768065,(t49:0.1520387747,(t116:0.01868664063,(t75:0.1386339963,(t102:0.09163444097,t178:0.01597096436):0.06122327271):0.0142621277):0.04638068884):0.05837448127):0.006005043284,(t23:0.1356035764,(t124:0.1095003806,t94:0.01166115622):0.005174865347):0.1123754984):0.08549847623):0.03848203141,(t133:0.064454232,(t80:0.1426635889,(t29:0.003532335453,((t36:0.08427990558,t197:0.02731956008):0.02526265024,t93:0.06034347072):0.05067758546):0.0842293226):0.1037653608):0.04149527346):0.001087278236,((((((((t105:0.05418593228,t5:0.05347499934):0.1092575268,(t87:0.05126456021,((t172:0.08901087761,(t193:0.1104958027,t63:0.04173081124):0.01167865165):0.04620158527,t76:0.05538311904):0.0260737547):0.03139896852):0.04302391868,((t168:0.01747855094,t182:0.004321316462):0.1154562131,(t176:0.2036736767,t22:0.03936042753):0.01222087371):0.03189806639):0.03843343927,(t98:0.3013555011,t154:0.2070581369):0.009036684462):0.02478079121,t24:0.02846779585):0.01158443985,t47:0.1386745288):0.04022715139,(((t126:0.06431536312,t140:0.03118945665):0.05312353763,t103:0.08847651209):0.07393791882,(t174:0.1347066037,(((t129:0.04628905913,t199:0.01558097455):0.01286585011,t187:0.03053488837):0.02995162052,t157:0.06532032986):0.04301296567):0.05240200038):0.1838428476):0.01762093497,(t77:0.05231337331,(t138:0.1683801303,((t86:0.0728897394,t119:0.05609098346):0.1437950228,t125:0.1924713237):0.1131018254):0.04496657346):0.01944226675):0.004583006969):0.05147142928):0.03956163653):8.491132634E-4,(((((((t163:0.1201147561,t44:0.08285147313):0.1228651137,((t194:0.1700679131,(((t68:0.1751727346,((t78:0.1008519948,(t82:0.06423817278,t145:0.03395378342):0.006710614292):0.03493791633,t54:0.08951110861):0.0298547759):0.02673622034,t169:0.05622630238):2.269156301E-4,((t43:0.09187936242,(t123:0.07939197877,t159:0.00885262305):0.09097112099):0.01686395903,t200:0.112520717):0.00410229035):0.0392765129):0.01044584409,(((t69:0.0993924397,t99:0.03716855707):0.1057649837,t72:0.2466977489):6.537301731E-4,(((t183:0.1221756465,(t8:0.04600100446,t185:0.1598975905):0.01236092288):0.03134936472,t64:0.2331719844):0.005748850337,(t35:0.09956972948,t95:0.1455439657):0.08297895819):0.01833473632):0.008558418572):0.004363695617):0.04776935424,(t48:0.1838206114,t162:0.04042607718):0.06426449193):0.01183601753,t107:0.05148569515):0.1154398066,t1:0.0528720883):0.01094055235,((t164:0.08798583573,t192:0.1655549963):0.2643199004,(((t188:0.07537181053,(t10:0.09000280238,t118:0.06889676539):0.0434599342):0.01906353948,t74:0.1521029803):0.116143132,(((((t177:0.1221571447,t158:0.182846772):0.03050588246,t161:0.174193549):0.02770764754,t113:0.1746231004):0.00467812037,t117:0.1940242414):0.0877704612,((t160:0.0195669223,((t141:0.02093340124,t167:0.1178564624):0.08021994961,(t189:0.2029063144,t79:0.04413767462):0.01343368462):0.1080166967):0.006278860901,((t67:0.02758579073,(t135:0.211053921,t175:0.1493736704):0.03223408488):0.02812430545,(t41:0.1852120201,(t65:0.10096633,t39:0.06974880375):0.0962984252):0.06948871982):0.0500274026):0.006290958099):0.002201313325):0.08829518442):0.02973761663):0.03205348426,((((t96:0.02855231257,t57:0.2366603091):0.08703694777,t128:0.2909812592):0.1068902459,((((t114:0.03874104846,t147:0.07489950667):0.01264666222,t38:0.1080981363):0.1915099756,t2:0.07914537186):0.1364513481,(t186:0.2896596271,((t27:0.2661713312,((t61:0.02629957646,t85:0.067533441):0.06356071446,(t150:0.1117478683,t97:0.1494987117):0.002895683023):0.09962172567):0.02035645508,((t19:0.04944628219,(t37:0.04422596834,t110:0.1635019172):0.006321885964):0.0986250419,(t139:0.09913687139,((((t28:0.05041488701,(t149:0.07203828491,t31:0.04965796547):0.0348763622):0.0174977869,t179:0.0432033841):0.02049505061,t146:0.1145895663):0.002492846607,(t52:0.06280646817,t89:0.01153453577):0.0872711871):0.1243045497):0.01218030615):0.001171646263):0.01605661203):0.1364112954):0.01345028016):0.03194573119,((t134:0.2725510659,(t143:0.06613600244,t13:0.05002345453):0.08925245306):0.06279893606,(((t130:0.2396989734,t137:0.02336619634):0.04278756619,(t18:0.08261003708,t81:0.1645683621):0.1034035373):0.01146696761,t132:0.08196130812):0.03957785721):0.1156576885):0.003568506952):0.01052412183):0.04825274753,((t3:0.2170329863,t4:0.04295718708):0.1287823831,(((t53:0.08968123545,t34:0.07908222397):0.005881300452,t92:0.04965904536):0.08001562771,(t109:0.02032432217,t112:0.05469748807):0.05596055895):0.1690398883):0.2018440042):0.05643582614):0.02748407343);";
	//        int nrTaxa = 200;
	//
	//        String orig = "0.136171103";
	//
	//        String locations = "t198=0, t142=1, t121=0, t1=0, t77=0, t173=0, t55=0, t47=0, t133=0, t24=0, t107=0, t144=0, t160=0, t156=0, t20=0, t2=0, t132=0, t137=0, t139=0, t59=0, t17=0, t162=0, t67=0, t9=0, t11=0, t96=1, t143=0, t13=0, t4=0, t169=0, t106=0, t101=0, t22=0, t154=0, t8=0, t79=0, t188=0, t37=0, t29=1, t138=1, t70=0, t116=0, t108=0, t76=0, t87=0, t6=0, t94=0, t19=0, t171=1, t80=0, t99=1, t194=1, t43=0, t44=1, t56=1, t111=0, t51=0, t141=0, t18=0, t14=1, t103=1, t93=0, t180=0, t58=0, t115=1, t168=0, t182=0, t179=0, t161=0, t200=0, t64=1, t21=0, t63=0, t155=0, t113=1, t109=1, t164=1, t16=0, t131=0, t100=0, t157=0, t48=0, t159=0, t61=0, t177=0, t147=0, t105=0, t118=0, t46=0, t178=0, t114=0, t185=0, t85=1, t69=0, t40=0, t128=1, t145=0, t66=0, t175=0, t28=0, t190=0, t74=1, t192=0, t60=0, t35=0, t49=0, t183=0, t7=0, t119=0, t112=1, t25=0, t197=0, t10=0, t187=0, t166=1, t82=1, t126=0, t148=0, t39=0, t54=0, t136=1, t123=0, t5=0, t172=1, t149=0, t26=0, t65=0, t153=0, t117=0, t15=0, t199=0, t89=0, t81=0, t140=0, t124=0, t102=0, t125=0, t42=0, t34=0, t176=0, t104=0, t150=0, t33=0, t146=0, t97=0, t195=0, t68=0, t90=0, t78=0, t174=0, t196=0, t72=0, t186=0, t88=0, t84=0, t92=0, t189=0, t62=0, t152=0, t170=0, t41=0, t38=0, t12=1, t134=0, t130=0, t36=0, t30=0, t167=0, t181=0, t83=0, t23=0, t3=0, t165=0, t73=0, t91=0, t135=0, t31=0, t50=0, t151=0, t163=0, t129=0, t127=0, t75=0, t120=0, t86=0, t53=0, t57=0, t122=0, t110=0, t191=1, t95=1, t32=0, t52=0, t27=0, t193=0, t98=0, t184=1, t158=0, t71=1, t45=0" ;
	//
	//        String stateNumber = "2";
	//        Double[] M = {1., 1.};
	//        String frequencies = "0.5 0.5";
	//
	//        Double[] b;
	//        Double[] d = {1.,1.};
	//        Double[] s;
	//        Double[] bij;
	//        Double[]  delta;
	//        Double[] R0;
	//        Double[]  samplingProportion;
	//        Double[]  R0AmongDemes;
	//        Double psi;
	//
	//        int maxEvals = 10000; //Integer.MAX_VALUE;
	//        double tolerance = 1e-12;
	//
	//        conditionOnSurvival = false;
	//
	//        double logL;
	//        System.out.println("b\t\td\ts\tm\t\tb-d-s\tb-d-s-m\tm*s\tb*s\t\t\tlogL");
	//
	//        for (double i = 1; i<=5; i+=0.25){
	//
	//            b = new Double[]{i, i};
	//
	//            psi = 0.5 * ((i - d[0]) + Math.sqrt((d[0] - i) * (d[0] - i) - .04));
	//
	//            s = new Double[] {psi,psi};
	//
	////            d = new Double[]{b[0]-s[0]-0.5, b[0]-s[0]-0.5};     // assume b[0]-d[0]-s[0]=0.5
	//            //            d = new Double[]{b[0]-s[0], b[0]-s[0]};     // assume b-d-s=0
	//            bij = new Double[]{3.,3.};
	//            M = new Double[]{b[0]-s[0]-d[0],b[0]-s[0]-d[0]};     // assume b-d-s=M
	////            bij = new Double[]{3.,3.};
	//
	//
	//            delta = new Double[]{b[0]+s[0],b[1]+d[1]};
	//            R0 = new Double[]{b[0]/delta[0], b[1]/delta[1]} ;
	//            samplingProportion = new Double[]{s[0]/delta[0], s[1]/delta[1]};
	//            R0AmongDemes = new Double[]{bij[0]/delta[0], bij[1]/delta[1]} ;
	//
	//            logL = bdm_likelihood(tolerance, maxEvals, stateNumber,
	//                    Double.toString(M[0])+" "+ Double.toString(M[1]),
	//                    frequencies,
	//                    tree, orig,
	//                    Double.toString(R0[0])+" "+ Double.toString(R0[1]),
	//                    null,//Double.toString(R0AmongDemes[0])+" "+ Double.toString(R0AmongDemes[1]),
	//                    Double.toString(delta[0])+" "+ Double.toString(delta[1]),
	//                    Double.toString(samplingProportion[0])+" "+ Double.toString(samplingProportion[1]),
	//                    "t", locations, nrTaxa, null);
	//
	////            System.out.print("R0 = "+Double.toString(R0[0])+" "+ Double.toString(R0[1])+ ", delta = " + Double.toString(becomeUninfectiousRate[0])+" "+ Double.toString(becomeUninfectiousRate[1]) + "\t\t");
	//            System.out.print(i + "\t"+ Math.round(d[0]) + "\t"+ Math.round(1000.*s[0])/1000.+ "\t"+ Math.round(1000.*M[0])/1000.+ "\t"+
	//                    Math.round(1000.*(b[0]-d[0]-s[0]))/1000.+ "\t"+ Math.round(1000.*(b[0]-d[0]-s[0]-M[0]))/1000. + "\t"+ Math.round(1000.*(M[0]*s[0]))/1000.+ "\t"+ Math.round(1000.*(b[0]*s[0]))/1000. + "\t\t\t");
	//            System.out.println(logL);
	//        }
	//
	//
	//
	//    }


	//    public void testCorrelationsAmongdemeInfections() throws Exception{
	//
	//        String stateNumber = "2";
	//        String frequencies = ".5 .5";
	//
	//        String newickStr = "";
	//        String locations = "";
	//        String origin = "";
	//
	//        String filename;
	//        File folder = new File("/Users/Denise/Documents/Projects/BDM/SIMS/tanjasims/equalRatesperLoc");
	//        File[] listOfFiles = folder.listFiles();
	//
	//        String files= "";
	//
	//        BufferedReader bf;
	//        BufferedWriter bw = new BufferedWriter(new PrintWriter("temp"));
	//
	//        for (int i = 0; i < listOfFiles.length; i++)
	//        {
	//
	//            if (listOfFiles[i].isFile())
	//            {
	//                filename = listOfFiles[i].getName();
	//                if (filename.endsWith(".xml") )
	//                {
	//                    files += filename + "\t";
	//                }
	//            }
	//        }
	//
	//        String[] filenames = files.split("\t");
	//
	//        for (int t=0;t<filenames.length;t++){
	//
	//            bf = new BufferedReader(new FileReader(filenames[t]));
	//
	//            String line;
	//            while ((line = bf.readLine()) != null) {
	//                if (line.contains("<input name='newick'>")){
	//                    line = bf.readLine();
	//                    newickStr = line.trim();
	//                    break;
	//                }
	//            }
	//
	//            if (newickStr.contains("[tree]")) break;
	//
	//            Tree tree = new TreeParser(newickStr);
	//
	//            while ((line = bf.readLine()) != null) {
	//                if (line.contains("id='tipTypes'")){
	//                    locations = (line.split("value",2)[1]).split("\"",3)[1];
	//                    break;
	//                }
	//            }
	//
	//            while ((line = bf.readLine()) != null) {
	//                if (line.contains("id=\"origin\"")) {
	//                    origin = (line.split("value",2)[1]).split("\"",3)[1];
	////                    origin = Double.toString(Double.parseDouble(origin) + tree.getRoot().getHeight());
	//                    break;
	//                }
	//            }
	//            bf.close();
	//
	//
	//            for (int par= 0; par<1;par++){
	//
	//                Double birth = 11.;
	//                Double b = birth;
	//                Double birth_ij = 4.;
	//                Double b_ij = birth;
	//                Double death = 3.;
	//                Double d = death;
	//                Double sampling = 1.;
	//                Double s = sampling;
	//
	//
	//                switch (par) {
	//                    case 0: bw = new BufferedWriter(new PrintWriter(filenames[t]+"likelihood_b.txt"));
	//                        bw.write("b\tlogL\n");
	//                        break;
	//                    case 1: bw = new BufferedWriter(new PrintWriter(filenames[t]+"likelihood_d.txt"));
	//                        bw.write("d\tlogL\n");
	//                        break;
	//                    case 2: bw = new BufferedWriter(new PrintWriter(filenames[t]+"likelihood_s.txt"));
	//                        bw.write("s\tlogL\n");
	//                        break;
	//                    case 3: bw = new BufferedWriter(new PrintWriter(filenames[t]+"likelihood_bij.txt"));
	//                        bw.write("b_ij\tlogL\n");
	//                        break;
	//                }
	//
	//                for (double i = 1.; i<2; i+=1){
	////                    for (double i = 0.; i<2.01; i+=0.1){
	//
	//                    switch (par) {
	//                        case 0: b = i*birth; break;
	//                        case 1: d = i*death;  break;
	//                        case 2: s = i*sampling;   break;
	//                        case 3: b_ij = i*birth_ij;  break;
	//                    }
	//
	//                    int maxEvals = Integer.MAX_VALUE;
	//                    double tolerance = 1e-14;
	//
	//                    double logL;
	//
	////                    logL = bdm_likelihood(tolerance, maxEvals, stateNumber,
	////                            "0 0",
	////                            frequencies,
	////                            newickStr, origin,
	////                            Double.toString(b/(d+s))+" "+ Double.toString(b/(d+s)),
	////                            Double.toString(b_ij/(d+s))+" "+ Double.toString(b_ij/(d+s)),
	////                            Double.toString(d+s)+" "+ Double.toString(d+s),
	////                            Double.toString(s/(s+d))+" "+ Double.toString(s/(s+d)),
	////                            "t",locations, 100);
	//
	//                    logL = bdm_likelihood_bds(tolerance, maxEvals, stateNumber,
	//                                                "0 0",
	//                                                frequencies,
	//                                                newickStr, origin,
	//                                                Double.toString(b)+" "+ Double.toString(b),
	//                                                Double.toString(b_ij)+" "+ Double.toString(b_ij),
	//                                                Double.toString(d)+" "+ Double.toString(d),
	//                                                Double.toString(s)+" "+ Double.toString(s),
	//                                                "t",locations, 100);
	////                    System.out.println(logL);
	//
	//                    switch (par) {
	//                        case 0: if (!Double.isInfinite(logL)) bw.write(b + "\t" + logL + "\n");  break;
	//                        case 1: if (!Double.isInfinite(logL)) bw.write(d + "\t" + logL + "\n");  break;
	//                        case 2: if (!Double.isInfinite(logL)) bw.write(s + "\t" + logL + "\n");  break;
	//                        case 3: if (!Double.isInfinite(logL)) bw.write(b_ij + "\t" + logL + "\n");  break;
	//                    }
	//                }
	//
	//                bw.close();
	//
	//            }
	//        }
	//    }

	public double bdm_likelihood_bds(double tolerance, int maxEvals, String statenumber, String migrationMatrix,
			String frequencies, String newick, String origin,
			String b, String b_ij, String d, String s, String prefixname, String locations, int nrTaxa) throws Exception {

		BirthDeathMigrationModelUncoloured bdm =  new BirthDeathMigrationModelUncoloured();

		ArrayList<Taxon> taxa = new ArrayList<Taxon>();

		for (int i=1; i<=nrTaxa; i++){
			taxa.add(new Taxon(prefixname+i));
		}

		Tree tree = new TreeParser();
		//        tree.setInputValue("singlechild", "true");
		tree.setInputValue("taxonset", new TaxonSet(taxa));
		tree.setInputValue("adjustTipHeights", "false");
		tree.setInputValue("IsLabelledNewick", "true");
		tree.setInputValue("newick", newick);
		//        tree.setInputValue("newick", "(((((t44:0.182968629,t62:0.3019052714):0.09847167097,(((((((t23:0.040758673,(t74:0.1307127987,t35:0.1154022735):0.07030908932):0.01715989135,(t82:0.1062952955,((t15:0.007155943854,t84:0.1719674607):0.004667546686,t94:0.07572352265):0.006718379801):0.02787255584):0.03983087577,((t18:0.1503687104,t41:0.00517218553):0.05055719264,t39:0.06084914919):0.05056031716):0.07357539995,(((t8:0.08983076687,t70:0.0602673479):0.05197103095,t6:0.2196846108):0.0901707549,(((t34:0.0586397845,t40:0.1023385556):0.01490358549,(t9:0.06963627153,t98:0.183394327):0.02337527028):0.02018601816,(t32:0.02809476828,t53:0.1707374525):0.001785851111):0.08559722349):0.01779689903):0.01418885378,(t92:0.1762833773,(t38:0.1539304027,(t27:0.1170180412,(t61:0.1227178702,t2:0.144810821):0.003192683375):0.001139845915):0.1843146936):0.007920800243):0.04112424996,((t16:0.01715790293,t43:0.02674833642):0.1545677883,((t86:0.201254029,((t63:0.1995157851,t13:0.2257755904):0.0674771468,(t14:0.132152036,t72:0.07387237694):0.0238011846):0.01548476507):0.008522214494,((t99:0.2272449936,((t58:0.06348799963,t65:0.03208196517):0.03970254066,(t5:0.02129213485,t4:0.04467775793):0.04580749779):0.1174956489):0.06393384565,(t11:0.2242223769,t19:0.08351957165):0.04015431429):0.04286224778):0.02560173234):0.01685236952):0.01204236248,((((t83:0.2256751468,t10:0.1056156082):0.05991535337,t81:0.285636386):0.07570776409,t1:0.1975438576):0.005568286208,(t52:0.05442396586,(t3:0.007200988319,t29:0.02201763651):0.1364240987):0.2013785717):0.01555024027):0.002032919483):0.0778253527,t60:0.1456192084):0.05286503306,(t77:0.2994429614,((t46:0.2040567792,(t89:0.09762328879,(t71:0.2126511715,(t69:0.04612659646,t75:0.05852091065):0.03041682186):0.05470761742):0.004322090681):0.01882752463,((t26:0.1767516723,t51:0.1681938365):0.01904119246,t73:0.1268723928):0.1031070826):0.0909379871):0.1430099209):0.02992172794,((((t80:0.120713748,t64:0.2067947338):0.1052132031,(t42:0.2389715643,t66:0.139852688):0.07058070354):0.1500738164,(((((t45:0.03270753351,t36:0.1068647215):0.05785459597,(t93:0.1752902944,((t90:0.0869827058,t17:0.01426846283):0.02048795492,t47:0.1216384825):0.02328058265):0.1358059751):0.008304575197,((((t87:0.2313098699,(t7:0.08291833377,t49:0.1970214693):0.004361280936):0.03206062526,(t59:0.01513464829,t68:0.131590872):0.1026827717):0.0460329278,(((t24:0.1785957156,t67:0.1356484613):0.1086621308,(t37:0.1409938951,t30:0.09531257371):0.1263314739):0.004909886847,(t28:0.1430192535,t91:0.09670155966):0.08685251483):0.02460117621):0.001695946365,t78:0.2596209727):0.004636910093):0.07660591504,t95:5.28382819E-4):0.008793314598,(((t21:0.1780868575,t55:0.2175705194):0.05985500293,t96:0.008640154697):0.05918641693,((t48:0.02505668795,t100:0.00625710341):0.1009423101,(t33:0.2570926519,(t88:0.2099263996,(t31:0.1174383217,t20:0.02089920882):0.05057703225):0.06284594847):0.04055708563):0.02674933951):0.06568593097):0.04955191208):0.02003232425,(t50:0.3174130563,(((t57:0.01339276685,((t97:0.05614760524,t22:0.05673482988):0.05596333356,t25:0.1436950487):0.02793840556):0.1306190568,t56:0.2900430941):0.01956284148,(t79:0.1893550417,(t12:0.09617127232,((t76:0.0525309926,t54:0.07045731661):0.01430433591,t85:0.04474821384):7.137115834E-4):0.105589246):0.127351523):0.08510348151):0.06736886999):0.08222055288);");
		tree.initAndValidate();



		TraitSet trait = new TraitSet();
		trait.setInputValue("taxa", new TaxonSet(taxa));
		trait.setInputValue("value", locations);
		trait.setInputValue("traitname", "tiptypes");
		trait.initAndValidate();

		bdm.setInputValue("tree", tree);
		bdm.setInputValue("tiptypes", trait);

		bdm.setInputValue("origin", Double.toString(Double.parseDouble(origin)+tree.getRoot().getHeight()));
		bdm.setInputValue("stateNumber", statenumber);
		bdm.setInputValue("migrationMatrix", migrationMatrix);
		bdm.setInputValue("frequencies", frequencies);

		bdm.setInputValue("birthRate", b);

		//        String R0AmongDemes = Double.toString(3./6.) + " " + Double.toString(2.5) ;
		if (b_ij != null) bdm.setInputValue("birthRateAmongDemes", b_ij);

		bdm.setInputValue("deathRate", d);
		bdm.setInputValue("samplingRate", s);
		bdm.setInputValue("maxEvaluations", maxEvals);
		bdm.setInputValue("conditionOnSurvival", "false");

		// TO DO uncomment or remove
		//bdm.setInputValue("tolerance", tolerance);

		bdm.initAndValidate();


		long startTime = System.currentTimeMillis();
		double logL = bdm.calculateLogP(); //calculateTreeLogLikelihood(coltree.getUncolouredTree());
		runtime = System.currentTimeMillis() - startTime;
		maxEvalsUsed = bdm.maxEvalsUsed;

		// TO DO REMOVE if not necessary
		assertEquals(0., 0., 1e-2);

		return logL;

	}

}
