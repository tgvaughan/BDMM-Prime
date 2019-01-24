package bdmm.mapping;

import bdmm.parameterization.EpiParameterization;
import bdmm.parameterization.Parameterization;
import bdmm.parameterization.SkylineMatrixParameter;
import bdmm.parameterization.SkylineVectorParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import org.junit.Test;

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

        TypeMappedTree typeMappedTree = new TypeMappedTree();
        typeMappedTree.initByName(
                "parameterization", parameterization,
                "frequencies", new RealParameter("0.5 0.5"),
                "untypedTree", tree,
                "typeLabel", "type");


    }
}
