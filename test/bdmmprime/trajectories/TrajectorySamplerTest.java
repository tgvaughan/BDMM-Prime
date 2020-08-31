package bdmmprime.trajectories;

import bdmmprime.distributions.BirthDeathMigrationDistribution;
import bdmmprime.mapping.TypeMappedTree;
import bdmmprime.parameterization.*;
import bdmmprime.trajectories.simulation.SimulatedTree;
import bdmmprime.trajectories.simulation.UntypedTreeFromTypedTree;
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
                "finalSampleOffset", new RealParameter("0.0"),
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

        SimulatedTree simulatedTree = new SimulatedTree();
        simulatedTree.initByName(
                "parameterization", parameterization,
                "frequencies", new RealParameter("1.0"),
                "minSamples", 2);

//        System.out.println(simulatedTree);
        System.out.println("Final sample offset: " + parameterization.getFinalSampleOffset());

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
                "finalSampleOffset", new RealParameter("0.0"),
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

        SimulatedTree simulatedTree = new SimulatedTree();
        simulatedTree.initByName(
                "parameterization", parameterization,
                "frequencies", new RealParameter("1.0"),
                "minSamples", 2);

//        System.out.println(simulatedTree);
        System.out.println("Final sample offset: " + parameterization.getFinalSampleOffset());

        SampledTrajectory sampledTrajectory = new SampledTrajectory();
        sampledTrajectory.initByName("typeMappedTree", simulatedTree,
                "parameterization", parameterization,
                "nParticles", 10000,
                "useTauLeaping", true,
                "stepsPerInterval", 100);

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
    public void untypedRateShiftLikelihoodTest() {
        Randomizer.setSeed(42);

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(1),
                "origin", new RealParameter("5.0"),
                "finalSampleOffset", new RealParameter("0.0"),
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

        SimulatedTree simulatedTree = new SimulatedTree();
        simulatedTree.initByName(
                "parameterization", parameterization,
                "frequencies", new RealParameter("1.0"),
                "minSamples", 2);

//        System.out.println(simulatedTree);
        System.out.println("Final sample offset: " + parameterization.getFinalSampleOffset());

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
                "finalSampleOffset", new RealParameter("0.0"),
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

        SimulatedTree simulatedTree = new SimulatedTree();
        simulatedTree.initByName(
                "parameterization", parameterization,
                "frequencies", new RealParameter("1.0"),
                "minSamples", 2);

//        System.out.println(simulatedTree);
        System.out.println("Final sample offset: " + parameterization.getFinalSampleOffset());

        SampledTrajectory sampledTrajectory = new SampledTrajectory();
        sampledTrajectory.initByName("typeMappedTree", simulatedTree,
                "parameterization", parameterization,
                "nParticles", 10000,
                "useTauLeaping", true,
                "stepsPerInterval", 100);

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
    public void untypedRhoSamplingLikelihoodTest() {
        Randomizer.setSeed(53);

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(1),
                "origin", new RealParameter("5.0"),
                "finalSampleOffset", new RealParameter("0.0"),
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

        SimulatedTree simulatedTree = new SimulatedTree();
        simulatedTree.initByName(
                "parameterization", parameterization,
                "frequencies", new RealParameter("1.0"),
                "minSamples", 2);

        System.out.println(simulatedTree);
        System.out.println("Final sample offset: " + parameterization.getFinalSampleOffset());

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
    public void typedSimpleLikelihoodTest() {
//        Randomizer.setSeed(1);
        System.out.println(Randomizer.getSeed());

        RealParameter frequencies = new RealParameter("0.5 0.5");
        int nTypes = 2;

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(nTypes),
                "origin", new RealParameter("5.0"),
                "finalSampleOffset", new RealParameter("0.0"), // To be set by simulation
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

//        SimulatedTree simulatedTree = new SimulatedTree();
//        simulatedTree.initByName(
//                "parameterization", parameterization,
//                "frequencies", frequencies,
//                "minSamples", 2,
//                "simulateUntypedTree", false);

        TreeParser simulatedTree = new TreeParser("(0[&type=0]:1.0,1[&type=0]:0.5)2[&type=0]:0.0;",
                false, true, true, 0);


        System.out.println(simulatedTree);
        System.out.println("Final sample offset: " + parameterization.getFinalSampleOffset());

        UntypedTreeFromTypedTree untypedSimulatedTree = new UntypedTreeFromTypedTree();
        untypedSimulatedTree.initByName(
                "typedTree", simulatedTree,
                "typeLabel", "type");
        System.out.println(untypedSimulatedTree);

//        SampledTrajectory sampledTrajTrueTree = new SampledTrajectory();
//        sampledTrajTrueTree.initByName(
//                "typeMappedTree", simulatedTree,
//                "parameterization", parameterization,
//                "frequencies", frequencies,
//                "nParticles", 1000,
//                "useTauLeaping", false,
//                "stepsPerInterval", 5);
//
//        System.out.println("Estimate from true typed tree:");
//        System.out.println(sampledTrajTrueTree.getLogTreeProbEstimate());

        TypeMappedTree typeMappedTree = new TypeMappedTree();
        typeMappedTree.initByName(
                "parameterization", parameterization,
                "frequencies", frequencies,
                "typeLabel", "type",
                "untypedTree", untypedSimulatedTree);

        SampledTrajectory sampledTrajectory = new SampledTrajectory();
        sampledTrajectory.initByName("typeMappedTree", typeMappedTree,
                "parameterization", parameterization,
                "frequencies", frequencies,
                "nParticles", 1000,
                "resampThresh", 0.5,
                "useTauLeaping", false,
                "stepsPerInterval", 5);

        double [] logProbEsts = new double[5];
        double maxLogProbEst = Double.NEGATIVE_INFINITY;
        for (int i=0; i<logProbEsts.length; i++) {
            typeMappedTree.initAndValidate();
            double thisEst = sampledTrajectory.getLogTreeProbEstimate();
            if (thisEst > maxLogProbEst)
                maxLogProbEst = thisEst;
            logProbEsts[i] = thisEst;
            System.out.println(thisEst);
        }

        double logProbEst = 0.0;
        for (int i=0; i<logProbEsts.length; i++)
            logProbEst += Math.exp(logProbEsts[i]-maxLogProbEst);
        logProbEst = Math.log(logProbEst/logProbEsts.length) + maxLogProbEst;

        System.out.println("Log probability estimate: " + logProbEst);

        BirthDeathMigrationDistribution bdmm = new BirthDeathMigrationDistribution();
        bdmm.initByName("parameterization", parameterization,
                "frequencies", frequencies,
                "typeLabel", "type",
                "conditionOnSurvival", false,
                "tree", untypedSimulatedTree);

        double logProbTrue = bdmm.calculateLogP();

        System.out.println("BDMM: Log probability true: " + logProbTrue);

        assertEquals(logProbEst, logProbTrue, 1e-1);
    }

    @Test
    public void tinyTypedTreeLikelihoodTest() {

        RealParameter frequencies = new RealParameter("0.5 0.5");
        int nTypes = 2;

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(nTypes),
                "origin", new RealParameter("1.0"),
                "finalSampleOffset", new RealParameter("0.0"),
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

//        TreeParser typedTree = new TreeParser("((0[&type=\"0\"]:0.75)2[&type=\"1\"]:0.25,1[&type=\"1\"]:0.5)3[&type=\"1\"]:0.2;",
//                false, true, true, 0);
        TreeParser typedTree = new TreeParser("(0[&type=\"0\"]:0.75)1[&type=\"1\"]:0.25",
                false, true, true, 0);

        System.out.println(typedTree);

        SampledTrajectory sampledTrajectory = new SampledTrajectory();
        sampledTrajectory.initByName("typeMappedTree", typedTree,
                "parameterization", parameterization,
                "frequencies", frequencies,
                "nParticles", 5,
                "resampThresh", 0.0,
                "useTauLeaping", false,
                "stepsPerInterval", 5);

        double logProbEst = sampledTrajectory.getLogTreeProbEstimate();
        System.out.println("Log probability estimate: " + logProbEst);

        double logProbTrue = -8.847311; // From R
        System.out.println("Log probability true: " + logProbTrue);

        assertEquals(logProbEst, logProbTrue, 1e-1);
    }

}
