package test.beast.evolution.speciation;

import beast.evolution.speciation.BirthDeathMigrationModel;
import beast.evolution.tree.MultiTypeTreeFromNewick;
import junit.framework.TestCase;
import org.junit.Test;

/**
 * Created by Julija Pecerska on 13/07/16.
 */
public class PiecewiseBirthDeathMigrationDistributionTest extends TestCase {
    private BirthDeathMigrationModel inputTestSetup() {
        String treeString = "(((((t1[&type=0]:0.4595008531,t25[&type=0]:0.4595008531)[&type=0]:0.3373053072,t23[&type=0]:0.3567584538)[&type=0]:0.007310819036,t16[&type=0]:0.3489190732)[&type=0]:0.331009529,((t18[&type=0]:0.03315384045,t14[&type=0]:0.03315384045)[&type=0]:0.5063451374,(t10[&type=0]:0.4211543131,t15[&type=0]:0.4211543131)[&type=0]:0.1183446648)[&type=0]:0.5956275305)[&type=0]:0.1158090878,((t19[&type=0]:0.9429393194,((t6[&type=0]:0.363527235,t11[&type=0]:0.4417423167)[&type=0]:0.01881829549,((((t3[&type=0]:0.3071904376,(((t24[&type=0]:0.01065209364,t13[&type=0]:0.01065209364)[&type=0]:0.06076485145,t8[&type=0]:0.07141694509)[&type=0]:0.123620245,(t22[&type=0]:0.1616119808,t2[&type=0]:0.1616119808)[&type=0]:0.03342520927)[&type=0]:0.1121532475)[&type=0]:0.24520579,t9[&type=0]:0.5523962276)[&type=0]:0.3852615426,(((t20[&type=0]:0.2935970782,(t17[&type=0]:0.06569090089,t4[&type=0]:0.06569090089)[&type=0]:0.2279061773)[&type=0]:0.08350780408,(t21[&type=0]:0.05109047139,t5[&type=0]:0.05109047139)[&type=0]:0.3260144109)[&type=0]:0.2298344132,t7[&type=0]:0.6069392955)[&type=0]:0.3307184747)[&type=0]:0.01206284377,t26[&type=0]:0.9497206139)[&type=0]:0.05755333197)[&type=0]:0.03290891884)[&type=0]:0.07263755325,t12[&type=0]:1.112820418)[&type=0]:0.1381151782)[&type=0];";

        MultiTypeTreeFromNewick tree = new MultiTypeTreeFromNewick();
        tree.initByName(
                "adjustTipHeights", false,
                "newick", treeString,
                "typeLabel", "type");
        BirthDeathMigrationModel bdssm =  new BirthDeathMigrationModel();
        bdssm.setInputValue("frequencies", "1");
        bdssm.setInputValue("migrationMatrix", "0.");
        bdssm.setInputValue("stateNumber", 1);
        bdssm.setInputValue("tree", tree);
        bdssm.setInputValue("conditionOnSurvival", false);

        return bdssm;
    }

    private void validateShouldFail(BirthDeathMigrationModel bdmm, Boolean fail) {
        Boolean result = fail;
        try {
            bdmm.initAndValidate();
        } catch (Exception e) {
            result = !result;
        } finally {
            if (result) {
                fail();
            }
        }
    }

    @Test
    public void testModelInputs() throws Exception {
        BirthDeathMigrationModel bdmm1 = inputTestSetup();
        bdmm1.setInputValue("birthRate", "1");
        bdmm1.setInputValue("deathRate", "1");
        bdmm1.setInputValue("samplingRate", "1");
        validateShouldFail(bdmm1, false);

        BirthDeathMigrationModel bdmm2 = inputTestSetup();
        bdmm2.setInputValue("birthRate", "1");
        bdmm2.setInputValue("deathRate", "1");
        validateShouldFail(bdmm2, true);

        BirthDeathMigrationModel bdmm3 = inputTestSetup();
        bdmm3.setInputValue("birthRate", "1");
        bdmm3.setInputValue("samplingRate", "1");
        validateShouldFail(bdmm3, true);

        BirthDeathMigrationModel bdmm4 = inputTestSetup();
        bdmm4.setInputValue("birthRate", "1");
        validateShouldFail(bdmm4, true);

        BirthDeathMigrationModel bdmm5 = inputTestSetup();
        bdmm5.setInputValue("R0", "1");
        bdmm5.setInputValue("becomeUninfectiousRate", "1.5");
        bdmm5.setInputValue("samplingProportion", "0.3");
        validateShouldFail(bdmm5, false);

        BirthDeathMigrationModel bdmm6 = inputTestSetup();
        bdmm6.setInputValue("R0", "1");
        validateShouldFail(bdmm6, true);

        BirthDeathMigrationModel bdmm7 = inputTestSetup();
        bdmm7.setInputValue("birthRate", "1");
        bdmm7.setInputValue("R0", "1");
        validateShouldFail(bdmm7, true);

        BirthDeathMigrationModel bdmm8 = inputTestSetup();
        bdmm8.setInputValue("R0", "1");
        bdmm8.setInputValue("becomeUninfectiousRate", "1.5");
        validateShouldFail(bdmm8, true);

        BirthDeathMigrationModel bdmm9 = inputTestSetup();
        bdmm9.setInputValue("R0", "1");
        bdmm9.setInputValue("samplingProportion", "0.3");
        validateShouldFail(bdmm9, true);

        BirthDeathMigrationModel bdmm10 = inputTestSetup();
        bdmm10.setInputValue("R0_base", "1");
        bdmm10.setInputValue("R0_ratio", "0.8");
        bdmm10.setInputValue("becomeUninfectiousRate", "1.5");
        bdmm10.setInputValue("samplingProportion", "0.3");
        validateShouldFail(bdmm10, false);

        BirthDeathMigrationModel bdmm11 = inputTestSetup();
        bdmm11.setInputValue("R0_base", "1");
        bdmm11.setInputValue("becomeUninfectiousRate", "1.5");
        bdmm11.setInputValue("samplingProportion", "0.3");
        validateShouldFail(bdmm11, true);

        BirthDeathMigrationModel bdmm12 = inputTestSetup();
        bdmm12.setInputValue("R0_ratio", "0.8");
        bdmm12.setInputValue("becomeUninfectiousRate", "1.5");
        bdmm12.setInputValue("samplingProportion", "0.3");
        validateShouldFail(bdmm12, true);

        BirthDeathMigrationModel bdmm13 = inputTestSetup();
        bdmm13.setInputValue("R0_base", "1");
        bdmm13.setInputValue("R0_ratio", "0.8");
        validateShouldFail(bdmm13, true);

        BirthDeathMigrationModel bdmm14 = inputTestSetup();
        bdmm14.setInputValue("R0", "1");
        bdmm14.setInputValue("R0_base", "1");
        bdmm14.setInputValue("R0_ratio", "0.8");
        bdmm14.setInputValue("becomeUninfectiousRate", "1.5");
        bdmm14.setInputValue("samplingProportion", "0.3");
        validateShouldFail(bdmm14, true);

        BirthDeathMigrationModel bdmm15 = inputTestSetup();
        bdmm15.setInputValue("R0", "1");
        bdmm15.setInputValue("R0_base", "1");
        bdmm15.setInputValue("becomeUninfectiousRate", "1.5");
        bdmm15.setInputValue("samplingProportion", "0.3");
        validateShouldFail(bdmm15, true);

        BirthDeathMigrationModel bdmm16 = inputTestSetup();
        bdmm16.setInputValue("R0", "1");
        bdmm16.setInputValue("R0_ratio", "0.8");
        bdmm16.setInputValue("becomeUninfectiousRate", "1.5");
        bdmm16.setInputValue("samplingProportion", "0.3");
        validateShouldFail(bdmm16, true);

        BirthDeathMigrationModel bdmm17 = inputTestSetup();
        bdmm17.setInputValue("R0_base", "1");
        bdmm17.setInputValue("R0_ratio", "0.8");
        bdmm17.setInputValue("deathRate", "1");
        bdmm17.setInputValue("samplingRate", "1");
        validateShouldFail(bdmm17, true);
    }

    @Test
    public void testRatios() throws Exception {
        String treeString = "((3[&state=0] : 1.5, 4[&state=1] : 0.5)[&state=1] : 1 , (1[&state=0] : 2, 2[&state=0] : 1)[&state=1] : 3)[&state=0];";

        MultiTypeTreeFromNewick tree = new MultiTypeTreeFromNewick();
        tree.initByName(
                "adjustTipHeights", false,
                "newick", treeString,
                "typeLabel", "type");

        BirthDeathMigrationModel bdmm = new BirthDeathMigrationModel();

        bdmm.setInputValue("frequencies", "0.5 0.5");
        bdmm.setInputValue("migrationMatrix", "0.2 0.2");
        bdmm.setInputValue("stateNumber", 2);
        bdmm.setInputValue("tree", tree);
        bdmm.setInputValue("conditionOnSurvival", false);
        bdmm.setInputValue("intervalTimes", "0.0 3.0");
        bdmm.setInputValue("becomeUninfectiousRate", "4.5 1.5 4.5 1.5");
        bdmm.setInputValue("samplingProportion", "0.4444444444 0.3333333333 0.4444444444 0.3333333333");

        bdmm.setInputValue("R0", "4.5 1.5 3.6 1.2");

        bdmm.initAndValidate();
        double likelihood = bdmm.calculateTreeLogLikelihood(tree);

        BirthDeathMigrationModel bdmm2 = new BirthDeathMigrationModel();

        bdmm2.setInputValue("frequencies", "0.5 0.5");
        bdmm2.setInputValue("migrationMatrix", "0.2 0.2");
        bdmm2.setInputValue("stateNumber", 2);
        bdmm2.setInputValue("tree", tree);
        bdmm2.setInputValue("conditionOnSurvival", false);
        bdmm2.setInputValue("intervalTimes", "0.0 3.0");
        bdmm2.setInputValue("becomeUninfectiousRate", "4.5 1.5 4.5 1.5");
        bdmm2.setInputValue("samplingProportion", "0.4444444444 0.3333333333 0.4444444444 0.3333333333");

        bdmm2.setInputValue("R0_base", "4.5 1.5");
        bdmm2.setInputValue("R0_ratio", "0.8");

        bdmm2.initAndValidate();
        double likelihood2 = bdmm2.calculateTreeLogLikelihood(tree);

        assertEquals(likelihood, likelihood2, 1e-10);
    }
}
