package bdmmflow.benchmark;

import bdmmprime.parameterization.Parameterization;
import beast.base.evolution.tree.Tree;

import java.util.*;


public class BenchmarkResult {

    long trial;
    Parameterization parameterization;
    Tree tree;
    BenchmarkRun flowRun;
    BenchmarkRun bdmmRun;
    boolean useInverseFlow;
    boolean useSplitting;
    String initialStateStrategy;
    int minNumInterval;
    boolean parallelized;

    List<String> flowMetricNames;
    List<String> bdmmMetricNames;

    public BenchmarkResult(
            long trial,
            Parameterization parameterization,
            Tree tree,
            BenchmarkRun flowRun,
            BenchmarkRun bdmmRun,
            boolean useInverseFlow,
            boolean useSplitting,
            String initialStateStrategy,
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

        this.flowMetricNames = new ArrayList<>(this.flowRun.loggedMetrics.keySet());
        this.bdmmMetricNames = new ArrayList<>(this.bdmmRun.loggedMetrics.keySet());
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
        joiner.add(this.initialStateStrategy);
        joiner.add(Integer.toString(this.minNumInterval));
        joiner.add(Boolean.toString(this.parallelized));

        for (String metricName : this.flowMetricNames) {
            joiner.add(this.flowRun.loggedMetrics.get(metricName));
        }
        for (String metricName : this.bdmmMetricNames) {
            joiner.add(this.bdmmRun.loggedMetrics.get(metricName));
        }

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

        for (String metricName : this.flowMetricNames) {
            joiner.add("flow_" + metricName);
        }
        for (String metricName : this.bdmmMetricNames) {
            joiner.add("bdmm_" + metricName);
        }

        return joiner.toString();
    }
}
