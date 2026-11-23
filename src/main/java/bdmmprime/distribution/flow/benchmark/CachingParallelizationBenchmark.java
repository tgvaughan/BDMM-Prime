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

import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

/**
 * Standalone entry point that times short MCMC chains of the BEAST analyses with the flow tree
 * prior run with and without parallelization and flow caching, writing the runtimes to CSV.
 * Only the tree prior settings change between configurations; the tree likelihood always uses one
 * thread per core.
 * <p>
 * The optional first argument is the {@code beast_benchmarks} directory of the
 * BDMM-Flow-Supplementary repository.
 */
public class CachingParallelizationBenchmark {

    static final String RESULTS_FILE = "caching_parallelization_results.csv";

    static final int NUM_REPLICATES = 3;

    /**
     * A combination of tree prior settings to time.
     */
    record Configuration(String method, boolean parallelize, boolean useFlowCaching) {
    }

    static final List<Configuration> CONFIGURATIONS = List.of(
            new Configuration("No parallelization, no caching", false, false),
            new Configuration("Parallelization, no caching", true, false),
            new Configuration("Parallelization and caching", true, true)
    );

    static void main(String[] args) throws Exception {
        Path benchmarkDirectory = BeastChainRunner.getBenchmarkDirectory(args);
        Path resultsFile = Path.of(RESULTS_FILE);

        BeastChainRunner.setUpBeast();
        BeastChainRunner.createResultsFile(resultsFile, getHeader());

        for (BeastChainRunner.Dataset dataset : BeastChainRunner.DATASETS) {
            Path xmlFile = benchmarkDirectory.resolve(dataset.xmlPath());

            for (Configuration configuration : CONFIGURATIONS) {
                runConfiguration(xmlFile, configuration, dataset.getWarmupChainLength(), 1);
            }

            // the configurations alternate within each replicate so that drifts in machine load affect all equally
            for (int replicate = 1; replicate <= NUM_REPLICATES; replicate++) {
                for (Configuration configuration : CONFIGURATIONS) {
                    long duration = runConfiguration(xmlFile, configuration, dataset.chainLength(), replicate);
                    BeastChainRunner.appendResult(resultsFile, getLine(dataset, configuration, replicate, duration));
                }
            }
        }

        BeastChainRunner.tearDownBeast();
    }

    static String getHeader() {
        StringJoiner joiner = new StringJoiner(",");

        joiner.add("dataset");
        joiner.add("method");
        joiner.add("parallelized");
        joiner.add("flow_caching");
        joiner.add("replicate");
        joiner.add("chain_length");
        joiner.add("runtime");

        return joiner.toString();
    }

    /**
     * Runs one chain with the given configuration. Within a replicate, all configurations use the
     * replicate as seed and therefore propose the same moves.
     */
    static long runConfiguration(Path xmlFile, Configuration configuration, long chainLength, long seed) throws Exception {
        Map<String, String> treePriorInputs = new LinkedHashMap<>(BeastChainRunner.DEFAULT_TREE_PRIOR_INPUTS);
        treePriorInputs.put("parallelize", Boolean.toString(configuration.parallelize()));
        treePriorInputs.put("useFlowCaching", Boolean.toString(configuration.useFlowCaching()));

        MCMC mcmc = BeastChainRunner.loadChain(xmlFile, treePriorInputs, chainLength, seed);
        return BeastChainRunner.runChain(mcmc);
    }

    /**
     * Formats one result row, with the runtime in seconds.
     */
    static String getLine(BeastChainRunner.Dataset dataset, Configuration configuration, int replicate, long duration) {
        StringJoiner joiner = new StringJoiner(",");

        joiner.add(dataset.name());
        joiner.add("\"%s\"".formatted(configuration.method()));
        joiner.add(Boolean.toString(configuration.parallelize()));
        joiner.add(Boolean.toString(configuration.useFlowCaching()));
        joiner.add(Integer.toString(replicate));
        joiner.add(Long.toString(dataset.chainLength()));
        joiner.add(Double.toString(duration / 1e9));

        return joiner.toString();
    }

}
