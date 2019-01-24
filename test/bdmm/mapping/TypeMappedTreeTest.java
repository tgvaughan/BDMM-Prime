package bdmm.mapping;

import bdmm.distributions.BirthDeathMigrationDistribution;
import bdmm.parameterization.EpiParameterization;
import bdmm.parameterization.Parameterization;
import bdmm.parameterization.SkylineMatrixParameter;
import bdmm.parameterization.SkylineVectorParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import org.junit.Test;

import static junit.framework.Assert.assertEquals;

public class TypeMappedTreeTest {

    @Test
    public void testBackwardIntegration() {

        Tree tree = new TreeParser(
                "((3[&type=0] : 1.5, 4[&type=1] : 0.5) : 1 , (1[&type=1] : 2, 2[&type=0] : 1) : 3);",
                false);

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "origin", new RealParameter("6.0"),
                "nTypes", 2,
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

        RealParameter frequencies = new RealParameter("0.5 0.5");

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("parameterization", parameterization,
                "frequencies", frequencies,
                "conditionOnSurvival", false,
                "tree", tree,
                "typeLabel", "type",
                "parallelize", false);

        double logProbTrue = density.calculateLogP();

        System.out.println(logProbTrue);

        TypeMappedTree typeMappedTree = new TypeMappedTree();
        typeMappedTree.initByName(
                "parameterization", parameterization,
                "frequencies", frequencies,
                "untypedTree", tree,
                "typeLabel", "type");

        double[] y = typeMappedTree.backwardsIntegrateSubtree(tree.getRoot(), 0.0);
        double logScaleFactor = typeMappedTree.geScaleFactors[tree.getRoot().getNr()];

        double logProb = 0.0;
        for (int type=0; type<parameterization.getNTypes(); type++) {
            logProb += y[type+parameterization.getNTypes()]*Math.exp(logScaleFactor)*frequencies.getValue(type);
        }

        logProb = Math.log(logProb);

        System.out.println(logProb);


        assertEquals(logProbTrue, logProb, 1e-5);

    }
}
