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

import bdmmprime.distribution.flow.flowSystems.InitialMatrixStrategy;
import bdmmprime.parameterization.Parameterization;
import beast.base.evolution.tree.Tree;

import java.util.*;


/**
 * One row of the benchmark output: the flow and classic runs over the same tree and
 * parameterization, together with the settings under which the flow run was performed.
 */
public class BenchmarkResult {

    long trial;
    Parameterization parameterization;
    Tree tree;
    BenchmarkRun flowRun;
    BenchmarkRun bdmmRun;
    boolean useInverseFlow;
    boolean useSplitting;
    InitialMatrixStrategy initialStateStrategy;
    int minNumInterval;
    boolean parallelized;

    public BenchmarkResult(
            long trial,
            Parameterization parameterization,
            Tree tree,
            BenchmarkRun flowRun,
            BenchmarkRun bdmmRun,
            boolean useInverseFlow,
            boolean useSplitting,
            InitialMatrixStrategy initialStateStrategy,
            int minNumInterval,
            boolean parallelized
    ) {
        this.trial = trial;
        this.parameterization = parameterization;
        this.tree = tree;
        this.flowRun = flowRun;
        this.bdmmRun = bdmmRun;
        this.useInverseFlow = useInverseFlow;
        this.useSplitting = useSplitting;
        this.initialStateStrategy = initialStateStrategy;
        this.minNumInterval = minNumInterval;
        this.parallelized = parallelized;
    }

    @Override
    public String toString() {
        StringJoiner joiner = new StringJoiner(",");

        joiner.add(Long.toString(this.trial));

        joiner.add(Integer.toString(this.tree.getNodeCount()));
        joiner.add(Integer.toString(this.tree.getLeafNodeCount()));

        joiner.add(Integer.toString(this.parameterization.getNTypes()));
        joiner.add(Double.toString(this.parameterization.getTotalProcessLength()));

        joiner.add(Double.toString(this.flowRun.likelihood));
        joiner.add(Long.toString(this.flowRun.duration));

        joiner.add(Double.toString(this.bdmmRun.likelihood));
        joiner.add(Long.toString(this.bdmmRun.duration));

        joiner.add(Boolean.toString(this.useInverseFlow));
        joiner.add(Boolean.toString(this.useSplitting));
        joiner.add(this.initialStateStrategy.name());
        joiner.add(Integer.toString(this.minNumInterval));
        joiner.add(Boolean.toString(this.parallelized));

        return joiner.toString();
    }

    public String getHeaders() {
        StringJoiner joiner = new StringJoiner(",");

        joiner.add("trial");

        joiner.add("node_count");
        joiner.add("leaf_count");

        joiner.add("types_count");
        joiner.add("process_length");

        joiner.add("flow_likelihood");
        joiner.add("flow_duration");

        joiner.add("bdmm_likelihood");
        joiner.add("bdmm_duration");

        joiner.add("use_inverse_flow");
        joiner.add("use_splitting");
        joiner.add("initial_state_strategy");
        joiner.add("min_num_intervals");
        joiner.add("parallelized");

        return joiner.toString();
    }
}
