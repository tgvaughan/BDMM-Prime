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
import jdk.jfr.consumer.RecordedStackTrace;
import jdk.jfr.consumer.RecordingFile;

import java.nio.file.Files;
import java.nio.file.Path;
import java.time.Duration;
import java.util.Arrays;
import java.util.EnumMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.concurrent.locks.LockSupport;

/**
 * Standalone entry point that profiles short MCMC chains of the BEAST analyses and attributes the
 * cost to the tree prior, the tree likelihood and everything else, writing the number of profiler
 * samples per operation to CSV. Two measures of cost are recorded in the same run:
 * <ul>
 *     <li>CPU time: the sampling profiler of the Java Flight Recorder records, at a fixed interval,
 *     the call stack of every thread executing Java code, including worker threads. Parallel work
 *     therefore counts once per busy thread, while time spent outside of Java code (e.g. garbage
 *     collection or waiting for other threads) is not counted.</li>
 *     <li>Wall-clock time: a sampler thread records, at a fixed interval, the call stack of the
 *     thread running the MCMC loop, whether it is computing or waiting for worker threads. As this
 *     thread drives every step of the chain, the samples are spread over the operations in
 *     proportion to the elapsed time spent in them.</li>
 * </ul>
 * A sample is attributed to the operation of the innermost stack frame that belongs to the tree
 * prior or the tree likelihood, so no instrumentation of the model code is needed.
 * <p>
 * The optional first argument is the {@code beast_benchmarks} directory of the
 * BDMM-Flow-Supplementary repository.
 */
public class RuntimeFractionBenchmark {

    static final String RESULTS_FILE = "runtime_fraction_results.csv";

    static final long SEED = 1;
    static final Duration SAMPLING_PERIOD = Duration.ofMillis(1);

    static final String EXECUTION_SAMPLE_EVENT = "jdk.ExecutionSample";
    static final String WALL_CLOCK_SAMPLER_NAME = "wall-clock-sampler";

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
     * The parts of an analysis to which the cost is attributed.
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

    /**
     * The measures of cost recorded by the profiler.
     */
    enum Measure {
        CPU_TIME("cpu_time"),
        WALL_CLOCK_TIME("wall_clock_time");

        final String label;

        Measure(String label) {
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

                Map<Measure, Map<Operation, Long>> sampleCounts = profileChain(loadChain(xmlFile, configuration, dataset.chainLength()));
                for (Measure measure : Measure.values()) {
                    for (Operation operation : Operation.values()) {
                        BeastChainRunner.appendResult(
                                resultsFile,
                                getLine(dataset, configuration, measure, operation, sampleCounts.get(measure).get(operation))
                        );
                    }
                }
            }
        }

        BeastChainRunner.tearDownBeast();
    }

    static String getHeader() {
        StringJoiner joiner = new StringJoiner(",");

        joiner.add("dataset");
        joiner.add("method");
        joiner.add("measure");
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
     * Runs the chain under both profilers and counts the samples per measure and operation.
     */
    static Map<Measure, Map<Operation, Long>> profileChain(MCMC mcmc) throws Exception {
        Path recordingFile = Files.createTempFile("bdmm-benchmark", ".jfr");
        recordingFile.toFile().deleteOnExit();

        WallClockSampler wallClockSampler = new WallClockSampler(Thread.currentThread());

        try (Recording recording = new Recording()) {
            recording.enable(EXECUTION_SAMPLE_EVENT).withPeriod(SAMPLING_PERIOD);
            recording.start();
            wallClockSampler.start();

            BeastChainRunner.runChain(mcmc);

            wallClockSampler.stopSampling();
            recording.stop();
            recording.dump(recordingFile);
        }

        Map<Measure, Map<Operation, Long>> sampleCounts = new EnumMap<>(Measure.class);
        sampleCounts.put(Measure.CPU_TIME, countCpuSamples(recordingFile));
        sampleCounts.put(Measure.WALL_CLOCK_TIME, wallClockSampler.getSampleCounts());
        return sampleCounts;
    }

    /**
     * Counts the execution samples of the flight recording per operation, leaving out the
     * samples of the wall-clock sampler itself.
     */
    static Map<Operation, Long> countCpuSamples(Path recordingFile) throws Exception {
        Map<Operation, Long> sampleCounts = createEmptySampleCounts();

        for (RecordedEvent event : RecordingFile.readAllEvents(recordingFile)) {
            if (!event.getEventType().getName().equals(EXECUTION_SAMPLE_EVENT))
                continue;

            if (WALL_CLOCK_SAMPLER_NAME.equals(event.getThread("sampledThread").getJavaName()))
                continue;

            sampleCounts.merge(classifyStackTrace(event.getStackTrace()), 1L, Long::sum);
        }

        return sampleCounts;
    }

    static Map<Operation, Long> createEmptySampleCounts() {
        Map<Operation, Long> sampleCounts = new EnumMap<>(Operation.class);
        for (Operation operation : Operation.values()) {
            sampleCounts.put(operation, 0L);
        }
        return sampleCounts;
    }

    static Operation classifyStackTrace(RecordedStackTrace stackTrace) {
        if (stackTrace == null)
            return Operation.OTHER;

        return classifyFrames(stackTrace.getFrames().stream().map(frame -> frame.getMethod().getType().getName()).toList());
    }

    /**
     * Attributes a sample to the operation of its innermost recognized frame, given the class
     * names of its frames from the innermost outwards. Library code called by the tree prior or
     * the tree likelihood (e.g. ODE integrators, substitution models or waiting for worker threads)
     * is thereby attributed to its caller. This also holds for stacks truncated by the flight
     * recorder, as truncation drops outer frames.
     */
    static Operation classifyFrames(List<String> classNames) {
        for (String className : classNames) {
            if (className.startsWith("bdmmprime.distribution.") || className.startsWith("bdmmprime.parameterization."))
                return Operation.TREE_PRIOR;

            if (className.contains(".evolution.likelihood."))
                return Operation.TREE_LIKELIHOOD;
        }

        return Operation.OTHER;
    }

    static String getLine(
            BeastChainRunner.Dataset dataset,
            Configuration configuration,
            Measure measure,
            Operation operation,
            long sampleCount
    ) {
        StringJoiner joiner = new StringJoiner(",");

        joiner.add(dataset.name());
        joiner.add(configuration.method());
        joiner.add(measure.label);
        joiner.add(operation.label);
        joiner.add(Long.toString(dataset.chainLength()));
        joiner.add(Long.toString(sampleCount));

        return joiner.toString();
    }

    /**
     * Samples the call stack of one thread at a fixed interval, regardless of whether the thread
     * is running, waiting or blocked.
     */
    static class WallClockSampler extends Thread {

        private final Thread target;
        private final Map<Operation, Long> sampleCounts = createEmptySampleCounts();
        private volatile boolean isSampling = true;

        WallClockSampler(Thread target) {
            super(WALL_CLOCK_SAMPLER_NAME);
            this.target = target;
            this.setDaemon(true);
        }

        @Override
        public void run() {
            while (this.isSampling) {
                List<String> classNames = Arrays.stream(this.target.getStackTrace())
                        .map(StackTraceElement::getClassName)
                        .toList();
                this.sampleCounts.merge(classifyFrames(classNames), 1L, Long::sum);

                LockSupport.parkNanos(SAMPLING_PERIOD.toNanos());
            }
        }

        void stopSampling() throws InterruptedException {
            this.isSampling = false;
            this.join();
        }

        Map<Operation, Long> getSampleCounts() {
            return this.sampleCounts;
        }
    }

}
