package bdmmflow.intervals;

import bdmmprime.parameterization.Parameterization;
import bdmmprime.util.Utils;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

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
