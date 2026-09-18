package bdmmprime.distribution.flow.benchmark;

import bdmmprime.distribution.flow.BirthDeathMigrationDistribution;
import bdmmprime.distribution.flow.flowSystems.InitialMatrixStrategy;
import bdmmprime.parameterization.*;
import bdmmprime.trajectories.simulation.SimulatedTree;
import beast.base.evolution.tree.Tree;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.type.Simplex;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.LinkedList;
import java.util.List;

public class Benchmark {

    static void main(String[] args) {
        int NUM_TRIALS = 30_000;

        ParameterizationSampler sampler = new ParameterizationSampler();

        List<BenchmarkResult> results = runBenchmarks(NUM_TRIALS, sampler);
        writeResults(results, "results.csv");
    }

    static List<BenchmarkResult> runBenchmarks(int numTrials, ParameterizationSampler sampler) {
        List<BenchmarkResult> results = new LinkedList<>();

        long start = System.currentTimeMillis();

        for (int i = 0; i < numTrials; i++) {
            InitialMatrixStrategy[] initialStateStrategies = new InitialMatrixStrategy[]{
                    InitialMatrixStrategy.identity,
                    InitialMatrixStrategy.random,
                    InitialMatrixStrategy.average_inverse,
            };
            Boolean[] choices = new Boolean[]{
                    false, true
            };

            Parameterization parameterization = sampler.sampleParameterization();
            Simplex startTypePriorProbs = sampler.sampleStartTypePriorProbs(parameterization);
            Tree tree = simulateTree(parameterization, startTypePriorProbs);
            int minNumIntervals = sampler.sampleMinIntervals();
            boolean parallelized = false;
            boolean useSplitting = false;

            for (Boolean useInverseFlow : choices) {
                for (InitialMatrixStrategy strategy : initialStateStrategies) {
                    BenchmarkRun bdmmRun = runBDMMBenchmark(tree, parameterization, startTypePriorProbs, parallelized);
                    BenchmarkRun flowRun = runFlowBenchmark(tree, parameterization, startTypePriorProbs, useInverseFlow, strategy, parallelized);
                    BenchmarkResult result = new BenchmarkResult(
                            start + i, parameterization, tree, flowRun, bdmmRun, useInverseFlow, useSplitting, strategy, minNumIntervals, parallelized
                    );
                    results.add(result);
                }
            }

            if (i % 100 == 0) {
                System.out.println(i);
                writeResults(results, "results.csv");
            }
        }

        return results;
    }

    static Tree simulateTree(Parameterization parameterization, Simplex startTypePriorProbs) throws IllegalStateException {
        SimulatedTree simulatedTree = new SimulatedTree();
        simulatedTree.parameterizationInput.setTypedValue(parameterization, simulatedTree);
        simulatedTree.finalSampleOffsetInput.setTypedValue(new RealScalarParam<>(0.0, NonNegativeReal.INSTANCE), simulatedTree);
        simulatedTree.startTypePriorProbsInput.setTypedValue(startTypePriorProbs, simulatedTree);
        simulatedTree.minSamplesInput.setTypedValue(2, simulatedTree);
        simulatedTree.simulateUntypedTreeInput.setTypedValue(true, simulatedTree);
        simulatedTree.initAndValidate();
        return simulatedTree;
    }

    static BenchmarkRun runFlowBenchmark(
            Tree tree,
            Parameterization parameterization,
            Simplex startTypePriorProbs,
            boolean useInverseFlow,
            InitialMatrixStrategy initialStateStrategy,
            boolean parallelized
    ) {
        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.parameterizationInput.setTypedValue(parameterization, density);
        density.treeInput.setTypedValue(tree, density);
        density.startTypePriorProbsInput.setTypedValue(startTypePriorProbs, density);
        density.typeLabelInput.setTypedValue("type", density);
        density.initialMatrixStrategyInput.setTypedValue(initialStateStrategy, density);
        density.useInverseFlowInput.setTypedValue(useInverseFlow, density);
        density.parallelizeInput.setTypedValue(parallelized, density);
        density.initAndValidate();

        long start = System.nanoTime();
        double likelihood = density.calculateLogP();
        long duration = System.nanoTime() - start;

        return new BenchmarkRun(duration, likelihood);
    }

    static BenchmarkRun runBDMMBenchmark(Tree tree, Parameterization parameterization, Simplex startTypePriorProbs, boolean parallelized) {
        bdmmprime.distribution.classic.BirthDeathMigrationDistribution density = new bdmmprime.distribution.classic.BirthDeathMigrationDistribution();
        density.initByName(
                "parameterization", parameterization,
                "tree", tree,
                "startTypePriorProbs", startTypePriorProbs,
                "typeLabel", "type",
                "parallelize", false
        );
        density.parameterizationInput.setTypedValue(parameterization, density);
        density.treeInput.setTypedValue(tree, density);
        density.startTypePriorProbsInput.setTypedValue(startTypePriorProbs, density);
        density.typeLabelInput.setTypedValue("type", density);
        density.parallelizeInput.setTypedValue(parallelized, density);
        density.initAndValidate();

        long start = System.nanoTime();
        double likelihood = density.calculateLogP();
        long duration = System.nanoTime() - start;

        return new BenchmarkRun(duration, likelihood);
    }

    static void writeResults(List<BenchmarkResult> results, String fileName) {
        File file = new File(fileName);
        boolean writeHeader = !file.exists() || file.length() == 0;

        try (FileWriter fileWriter = new FileWriter(fileName, true)) {
            if (writeHeader) {
                fileWriter.write(results.get(0).getHeaders());
                fileWriter.write("\n");
            }

            for (BenchmarkResult result : results) {
                fileWriter.write(result.toString());
                fileWriter.write("\n");
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        LoggedMetric.storeMetrics("metrics.csv");
    }

}
