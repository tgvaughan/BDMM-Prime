package bdmmflow.benchmark;

import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class LoggedMetric {
    public static List<List<String>> loggedMetrics = new ArrayList<>();

    public static void logMetric(String algorithm, String metricName, double time, double value) {
        loggedMetrics.add(List.of(algorithm, metricName, String.valueOf(time), String.valueOf(value)));
    }

    public static void storeMetrics(String fileName) {
        try (FileWriter fileWriter = new FileWriter(fileName)) {
            fileWriter.write("algorithm,metric,time,value\n");

            for (List<String> result : loggedMetrics) {
                fileWriter.write(String.join(",", result));
                fileWriter.write("\n");
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
}
