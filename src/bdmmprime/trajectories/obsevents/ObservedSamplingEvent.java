package bdmmprime.trajectories.obsevents;

import bdmmprime.parameterization.Parameterization;
import bdmmprime.trajectories.Trajectory;
import bdmmprime.trajectories.trajevents.DeathEvent;
import beast.util.Randomizer;

public class ObservedSamplingEvent extends ObservedEvent {
    public int nLeaves, nSampledAncestors;

    public ObservedSamplingEvent(double time, int type, int nLeaves, int nSampledAncestors) {
        this.time = time;
        this.type = type;
        this.nLeaves = nLeaves;
        this.nSampledAncestors = nSampledAncestors;
        this.multiplicity = nLeaves + nSampledAncestors;
    }

    @Override
    public int[] getNextLineageCounts() {
        int[] nextLineageCounts = lineages.clone();
        nextLineageCounts[type] -= nLeaves;

        return nextLineageCounts;
    }

    @Override
    public double applyToTrajectory(Parameterization param, int interval, Trajectory trajectory) {

        // TODO: Add in treatment of rho sampling.

        double logWeightContrib = 0;

        int s = type;


        for (int i=0; i<nLeaves; i++) {
            double sampling_prop = trajectory.currentState[s]*param.getSamplingRates()[interval][s];
            logWeightContrib += Math.log(sampling_prop);

            boolean isRemoval = Randomizer.nextDouble() < param.getRemovalProbs()[interval][s];
            if (isRemoval) {
                trajectory.addEvent(new DeathEvent(time, s));
            } else {
                logWeightContrib += Math.log(1.0 - (lineages[s]-1.0)/trajectory.currentState[s]);
            }
        }

        for (int i=0; i<nSampledAncestors; i++) {
            double sampling_prop = trajectory.currentState[s]*param.getSamplingRates()[interval][s];
            logWeightContrib += Math.log((1.0-param.getRemovalProbs()[interval][s])*sampling_prop);
        }

        return logWeightContrib;
    }
}
