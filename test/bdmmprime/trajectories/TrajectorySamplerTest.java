package bdmmprime.trajectories;

import bdmmprime.distributions.BirthDeathMigrationDistribution;
import bdmmprime.parameterization.*;
import bdmmprime.trajectories.simulation.SimulatedTree;
import beast.core.parameter.RealParameter;
import beast.util.Randomizer;
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
}
