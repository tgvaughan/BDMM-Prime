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

package bdmmprime.distribution.flow.intervals;

import bdmmprime.parameterization.Parameterization;
import bdmmprime.util.Utils;

import java.util.ArrayList;
import java.util.List;

public class IntervalUtils {
    /**
     * Returns the list of intervals for the given parameterization intervals.
     */
    public static List<Interval> getIntervals(Parameterization parameterization) {
        List<Interval> intervals = new ArrayList<>();

        int currentParameterizationInterval = 0;
        double currentStartTime = Math.min(0.0, parameterization.getIntervalEndTimes()[0]);

        while (currentParameterizationInterval < parameterization.getTotalIntervalCount()) {
            double currentParameterizationIntervalEndTime = parameterization.getIntervalEndTimes()[currentParameterizationInterval];

            if (Utils.equalWithPrecision(currentParameterizationIntervalEndTime, currentStartTime) || currentParameterizationIntervalEndTime < currentStartTime) {
                // the current interval is empty or ends before it starts. this can happen when the interval has a negative end time.
                currentParameterizationInterval += 1;
                continue;
            }

            intervals.add(
                    new Interval(
                            intervals.size(),
                            currentParameterizationInterval,
                            currentStartTime,
                            currentParameterizationIntervalEndTime
                    )
            );

            currentStartTime = currentParameterizationIntervalEndTime;
            currentParameterizationInterval += 1;
        }

        return intervals;
    }

}
