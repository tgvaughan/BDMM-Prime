package bdmmprime.trajectories;

import bdmmprime.distribution.BirthDeathMigrationDistribution;
import bdmmprime.parameterization.*;
import bdmmprime.trajectories.simulation.SimulatedTree;
import beast.core.parameter.RealParameter;
import beast.util.Randomizer;
import beast.util.TreeParser;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class TrajectorySamplerTest {

    @Test
    public void untypedSimpleLikelihoodTest() {
        Randomizer.setSeed(53);
//        Randomizer.setSeed(42);

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(1),
                "origin", new RealParameter("5.0"),
                "birthRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("2.0")),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0")),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.5")),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0")));

        RealParameter finalSampleOffset = new RealParameter("0.0");

        SimulatedTree simulatedTree = new SimulatedTree();
        simulatedTree.initByName(
                "parameterization", parameterization,
                "finalSampleOffset", finalSampleOffset,
                "frequencies", new RealParameter("1.0"),
                "minSamples", 2);

//        System.out.println(simulatedTree);
        System.out.println("Final sample offset: " + finalSampleOffset.getValue());

        SampledTrajectory sampledTrajectory = new SampledTrajectory();
        sampledTrajectory.initByName("typeMappedTree", simulatedTree,
                "parameterization", parameterization,
                "nParticles", 100000);

        double logProbEst = sampledTrajectory.getLogTreeProbEstimate();
        System.out.println("Log probability estimate: " + logProbEst);

        BirthDeathMigrationDistribution bdmm = new BirthDeathMigrationDistribution();
        bdmm.initByName("parameterization", parameterization,
                "frequencies", new RealParameter("1.0"),
                "typeLabel", "type",
                "conditionOnSurvival", false,
                "tree", simulatedTree);

        double logProbTrue = bdmm.calculateLogP();

        System.out.println("Log probability true: " + logProbTrue);

        assertEquals(logProbEst, logProbTrue, 1e-1);
    }

    @Test
    public void untypedSimpleTLLikelihoodTest() {
//        Randomizer.setSeed(53);
        Randomizer.setSeed(42);

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(1),
                "origin", new RealParameter("5.0"),
                "birthRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("2.0")),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0")),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.5")),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0")));

        RealParameter finalSampleOffset = new RealParameter("0.0");

        SimulatedTree simulatedTree = new SimulatedTree();
        simulatedTree.initByName(
                "parameterization", parameterization,
                "finalSampleOffset", finalSampleOffset,
                "frequencies", new RealParameter("1.0"),
                "minSamples", 2);

//        System.out.println(simulatedTree);
        System.out.println("Final sample offset: " + finalSampleOffset.getValue());

        SampledTrajectory sampledTrajectory = new SampledTrajectory();
        sampledTrajectory.initByName("typeMappedTree", simulatedTree,
                "parameterization", parameterization,
                "nParticles", 10000,
                "useTauLeaping", true,
                "minLeapCount", 100,
                "epsilon", 0.01);

        double logProbEst = sampledTrajectory.getLogTreeProbEstimate();
        System.out.println("Log probability estimate: " + logProbEst);

        BirthDeathMigrationDistribution bdmm = new BirthDeathMigrationDistribution();
        bdmm.initByName("parameterization", parameterization,
                "frequencies", new RealParameter("1.0"),
                "typeLabel", "type",
                "conditionOnSurvival", false,
                "tree", simulatedTree);

        double logProbTrue = bdmm.calculateLogP();

        System.out.println("Log probability true: " + logProbTrue);

        assertEquals(logProbTrue, logProbEst, 1e-1);
    }

    @Test
    public void untypedRateShiftLikelihoodTest() {
        Randomizer.setSeed(42);

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(1),
                "origin", new RealParameter("5.0"),
                "birthRate", new SkylineVectorParameter(
                        new RealParameter("2.5"),
                        new RealParameter("2.0 1.0"), 1),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0")),
                "samplingRate", new SkylineVectorParameter(
                        new RealParameter("2"),
                        new RealParameter("0.0 0.5"), 1),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0")));


        RealParameter finalSampleOffset = new RealParameter("0.0");

        SimulatedTree simulatedTree = new SimulatedTree();
        simulatedTree.initByName(
                "parameterization", parameterization,
                "finalSampleOffset", finalSampleOffset,
                "frequencies", new RealParameter("1.0"),
                "minSamples", 2);

//        System.out.println(simulatedTree);
        System.out.println("Final sample offset: " + finalSampleOffset.getValue());

        SampledTrajectory sampledTrajectory = new SampledTrajectory();
        sampledTrajectory.initByName("typeMappedTree", simulatedTree,
                "parameterization", parameterization,
                "nParticles", 100000);

        double logProbEst = sampledTrajectory.getLogTreeProbEstimate();
        System.out.println("Log probability estimate: " + logProbEst);

        BirthDeathMigrationDistribution bdmm = new BirthDeathMigrationDistribution();
        bdmm.initByName("parameterization", parameterization,
                "frequencies", new RealParameter("1.0"),
                "typeLabel", "type",
                "conditionOnSurvival", false,
                "tree", simulatedTree);

        double logProbTrue = bdmm.calculateLogP();

        System.out.println("Log probability true: " + logProbTrue);

        assertEquals(logProbEst, logProbTrue, 1e-1);
    }


    @Test
    public void untypedRateShiftTLLikelihoodTest() {
        Randomizer.setSeed(42);

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(1),
                "origin", new RealParameter("5.0"),
                "birthRate", new SkylineVectorParameter(
                        new RealParameter("2.5"),
                        new RealParameter("2.0 1.0"), 1),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0")),
                "samplingRate", new SkylineVectorParameter(
                        new RealParameter("2"),
                        new RealParameter("0.0 0.5"), 1),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0")));


        RealParameter finalSampleOffset = new RealParameter("0.0");

        SimulatedTree simulatedTree = new SimulatedTree();
        simulatedTree.initByName(
                "parameterization", parameterization,
                "finalSampleOffset", finalSampleOffset,
                "frequencies", new RealParameter("1.0"),
                "minSamples", 2);

//        System.out.println(simulatedTree);
        System.out.println("Final sample offset: " + finalSampleOffset.getValue());

        SampledTrajectory sampledTrajectory = new SampledTrajectory();
        sampledTrajectory.initByName("typeMappedTree", simulatedTree,
                "parameterization", parameterization,
                "nParticles", 10000,
                "useTauLeaping", true,
                "epsilon", 0.01,
                "minLeapCount", 100);

        double logProbEst = sampledTrajectory.getLogTreeProbEstimate();
        System.out.println("Log probability estimate: " + logProbEst);

        BirthDeathMigrationDistribution bdmm = new BirthDeathMigrationDistribution();
        bdmm.initByName("parameterization", parameterization,
                "frequencies", new RealParameter("1.0"),
                "typeLabel", "type",
                "conditionOnSurvival", false,
                "tree", simulatedTree);

        double logProbTrue = bdmm.calculateLogP();

        System.out.println("Log probability true: " + logProbTrue);

        assertEquals(logProbTrue, logProbEst, 1e-1);
    }

    @Test
    public void untypedRhoSamplingLikelihoodTest() {
        Randomizer.setSeed(53);

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(1),
                "origin", new RealParameter("5.0"),
                "birthRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("2.0")),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0")),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.5")),
                "rhoSampling", new TimedParameter(
                        new RealParameter("4.0"),
                        new RealParameter("0.5")),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0")));

        RealParameter finalSampleOffset = new RealParameter("0.0");

        SimulatedTree simulatedTree = new SimulatedTree();
        simulatedTree.initByName(
                "parameterization", parameterization,
                "finalSampleOffset", finalSampleOffset,
                "frequencies", new RealParameter("1.0"),
                "minSamples", 2);

        System.out.println(simulatedTree);
        System.out.println("Final sample offset: " + finalSampleOffset.getValue());

        SampledTrajectory sampledTrajectory = new SampledTrajectory();
        sampledTrajectory.initByName("typeMappedTree", simulatedTree,
                "parameterization", parameterization,
                "nParticles", 100000);

        double logProbEst = sampledTrajectory.getLogTreeProbEstimate();
        System.out.println("Log probability estimate: " + logProbEst);

        BirthDeathMigrationDistribution bdmm = new BirthDeathMigrationDistribution();
        bdmm.initByName("parameterization", parameterization,
                "frequencies", new RealParameter("1.0"),
                "typeLabel", "type",
                "conditionOnSurvival", false,
                "tree", simulatedTree);

        double logProbTrue = bdmm.calculateLogP();

        System.out.println("Log probability true: " + logProbTrue);

        assertEquals(logProbEst, logProbTrue, 1e-1);
    }

    @Test
    public void tinyTypedTreeLikelihoodTest() {

        RealParameter frequencies = new RealParameter("0.5 0.5");
        int nTypes = 2;

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(nTypes),
                "origin", new RealParameter("1.2"),
                "birthRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("2.0"), nTypes),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0"), nTypes),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.5"), nTypes),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0"), nTypes),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        new RealParameter("1.0"), nTypes));

        TreeParser typedTree = new TreeParser("((0[&type=\"0\"]:0.75)2[&type=\"1\"]:0.25,1[&type=\"1\"]:0.5)3[&type=\"1\"]:0.2;",
                false, true, true, 0);

        System.out.println(typedTree);

        RealParameter finalSampleOffset =  new RealParameter("0.0");

        SampledTrajectory sampledTrajectory = new SampledTrajectory();
        sampledTrajectory.initByName("typeMappedTree", typedTree,
                "parameterization", parameterization,
                "finalSampleOffset", finalSampleOffset,
                "frequencies", frequencies,
                "nParticles", 100000,
                "resampThresh", 0.0,
                "useTauLeaping", false);

        double logProbEst = sampledTrajectory.getLogTreeProbEstimate();
        System.out.println("Log probability estimate: " + logProbEst);

        double logProbTrue = -4.365784; // From R script validation/trajectories/multi_type/mttree_prob.R
        System.out.println("Log probability true: " + logProbTrue);

        assertEquals(logProbTrue, logProbEst, 1e-1);
    }

}
