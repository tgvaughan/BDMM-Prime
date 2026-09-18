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
