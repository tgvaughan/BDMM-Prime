package bdmmflow.benchmark;

import java.util.HashMap;
import java.util.Map;

public class BenchmarkRun {
    public static Map<String, String> lastLoggedMetrics = new HashMap<>();

    public static void logMetric(String key, String value) {
        lastLoggedMetrics.put(key, value);
    }

    public static void addToMetric(String key, Double value) {
        lastLoggedMetrics.merge(
                key,
                String.valueOf(value),
                (x, y) -> String.valueOf((Double.parseDouble(x) + Double.parseDouble(y)))
        );
    }

    long duration;
    double likelihood;
    Map<String, String> loggedMetrics;

    public BenchmarkRun(long duration, double likelihood) {
        this.duration = duration;
        this.likelihood = likelihood;

        this.loggedMetrics = this.lastLoggedMetrics;
        this.lastLoggedMetrics = new HashMap<>();
    }
}
