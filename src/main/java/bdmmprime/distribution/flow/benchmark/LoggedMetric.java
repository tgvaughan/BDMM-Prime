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
