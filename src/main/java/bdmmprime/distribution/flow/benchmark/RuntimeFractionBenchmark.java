/*
 * Copyright (c) 2017-2026 ETH Zürich
 *
 * This file is part of bdmm-prime.
 *
 * bdmm-prime is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at
 * your option) any later version.
 *
 * bdmm-prime is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with bdmm-prime. If not, see <https://www.gnu.org/licenses/>.
 */

package bdmmprime.distribution.flow.benchmark;

import beast.base.inference.MCMC;
import jdk.jfr.Recording;
import jdk.jfr.consumer.RecordedEvent;
import jdk.jfr.consumer.RecordedFrame;
import jdk.jfr.consumer.RecordedStackTrace;
import jdk.jfr.consumer.RecordingFile;

import java.nio.file.Files;
import java.nio.file.Path;
import java.time.Duration;
import java.util.EnumMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

/**
 * Standalone entry point that profiles short MCMC chains of the BEAST analyses and attributes the
 * computational cost to the tree prior, the tree likelihood and everything else, writing the
 * number of profiler samples per operation to CSV.
 * <p>
 * The profiling uses the sampling profiler of the Java Flight Recorder: at a fixed interval, it
 * records the call stack of every thread executing Java code, including worker threads. A sample
 * is attributed to the operation of the innermost stack frame that belongs to the tree prior or
 * the tree likelihood, so no instrumentation of the model code is needed.
 * <p>
 * The optional first argument is the {@code beast_benchmarks} directory of the
 * BDMM-Flow-Supplementary repository.
 */
public class RuntimeFractionBenchmark {

    static final String RESULTS_FILE = "runtime_fraction_results.csv";

    static final long SEED = 1;
    static final Duration SAMPLING_PERIOD = Duration.ofMillis(1);

    static final String EXECUTION_SAMPLE_EVENT = "jdk.ExecutionSample";

    /**
     * A tree prior implementation to profile.
     */
    record Configuration(String method, Map<String, String> treePriorInputs) {
    }

    static final List<Configuration> CONFIGURATIONS = List.of(
            new Configuration("IF AIP", Map.of()),
            new Configuration("BDMM-Prime", Map.of("engine", "classic"))
    );

    /**
     * The parts of an analysis to which the computational cost is attributed.
     */
    enum Operation {
        TREE_PRIOR("Tree Prior"),
        TREE_LIKELIHOOD("Tree Likelihood"),
        OTHER("Other");

        final String label;

        Operation(String label) {
            this.label = label;
        }
    }

    static void main(String[] args) throws Exception {
        Path benchmarkDirectory = BeastChainRunner.getBenchmarkDirectory(args);
        Path resultsFile = Path.of(RESULTS_FILE);

        BeastChainRunner.setUpBeast();
        BeastChainRunner.createResultsFile(resultsFile, getHeader());

        for (BeastChainRunner.Dataset dataset : BeastChainRunner.DATASETS) {
            Path xmlFile = benchmarkDirectory.resolve(dataset.xmlPath());

            for (Configuration configuration : CONFIGURATIONS) {
                BeastChainRunner.runChain(loadChain(xmlFile, configuration, dataset.getWarmupChainLength()));

                Map<Operation, Long> sampleCounts = profileChain(loadChain(xmlFile, configuration, dataset.chainLength()));
                for (Operation operation : Operation.values()) {
                    BeastChainRunner.appendResult(
                            resultsFile,
                            getLine(dataset, configuration, operation, sampleCounts.get(operation))
                    );
                }
            }
        }

        BeastChainRunner.tearDownBeast();
    }

    static String getHeader() {
        StringJoiner joiner = new StringJoiner(",");

        joiner.add("dataset");
        joiner.add("method");
        joiner.add("operation");
        joiner.add("chain_length");
        joiner.add("sample_count");

        return joiner.toString();
    }

    static MCMC loadChain(Path xmlFile, Configuration configuration, long chainLength) throws Exception {
        Map<String, String> treePriorInputs = new LinkedHashMap<>(BeastChainRunner.DEFAULT_TREE_PRIOR_INPUTS);
        treePriorInputs.putAll(configuration.treePriorInputs());

        return BeastChainRunner.loadChain(xmlFile, treePriorInputs, chainLength, SEED);
    }

    /**
     * Runs the chain under the flight recorder and counts the execution samples per operation.
     */
    static Map<Operation, Long> profileChain(MCMC mcmc) throws Exception {
        Path recordingFile = Files.createTempFile("bdmm-benchmark", ".jfr");
        recordingFile.toFile().deleteOnExit();

        try (Recording recording = new Recording()) {
            recording.enable(EXECUTION_SAMPLE_EVENT).withPeriod(SAMPLING_PERIOD);
            recording.start();
            BeastChainRunner.runChain(mcmc);
            recording.stop();
            recording.dump(recordingFile);
        }

        Map<Operation, Long> sampleCounts = new EnumMap<>(Operation.class);
        for (Operation operation : Operation.values()) {
            sampleCounts.put(operation, 0L);
        }

        for (RecordedEvent event : RecordingFile.readAllEvents(recordingFile)) {
            if (event.getEventType().getName().equals(EXECUTION_SAMPLE_EVENT)) {
                sampleCounts.merge(classifySample(event.getStackTrace()), 1L, Long::sum);
            }
        }

        return sampleCounts;
    }

    /**
     * Attributes a sample to the operation of its innermost recognized frame. The frames are
     * ordered from the innermost outwards, so library code called by the tree prior or the tree
     * likelihood (e.g. ODE integrators or substitution models) is attributed to its caller. This
     * also holds for stacks truncated by the flight recorder, as truncation drops outer frames.
     */
    static Operation classifySample(RecordedStackTrace stackTrace) {
        if (stackTrace == null)
            return Operation.OTHER;

        for (RecordedFrame frame : stackTrace.getFrames()) {
            String className = frame.getMethod().getType().getName();

            if (className.startsWith("bdmmprime.distribution.") || className.startsWith("bdmmprime.parameterization."))
                return Operation.TREE_PRIOR;

            if (className.contains(".evolution.likelihood."))
                return Operation.TREE_LIKELIHOOD;
        }

        return Operation.OTHER;
    }

    static String getLine(BeastChainRunner.Dataset dataset, Configuration configuration, Operation operation, long sampleCount) {
        StringJoiner joiner = new StringJoiner(",");

        joiner.add(dataset.name());
        joiner.add(configuration.method());
        joiner.add(operation.label);
        joiner.add(Long.toString(dataset.chainLength()));
        joiner.add(Long.toString(sampleCount));

        return joiner.toString();
    }

}
