package bdmmprime.trajectories.obsevents;

import bdmmprime.parameterization.Parameterization;
import bdmmprime.trajectories.Trajectory;
import bdmmprime.trajectories.trajevents.BirthEvent;

public class CoalescenceEvent extends ObservedEvent {

    public int childType1, childType2;

    public CoalescenceEvent(double time, int parentType, int childType1, int childType2) {
        this.time = time;
        this.type = parentType;
        this.childType1 = childType1;
        this.childType2 = childType2;
    }

    @Override
    public int[] getNextLineageCounts() {
        int[] nextLineageCounts = lineages.clone();

        nextLineageCounts[type] -= 1;
        nextLineageCounts[childType1] += 1;
        nextLineageCounts[childType2] += 1;

        return nextLineageCounts;
    }

    @Override
    public double applyToTrajectory(Parameterization param, int interval, Trajectory trajectory) {

        int s = type;
        int sp = Math.max(childType1, childType2);
        int spp = Math.min(childType1, childType2);

        double logWeightContrib = trajectory.currentState[s]*param.getBirthRates()[interval][s][sp][spp];

        trajectory.addEvent(new BirthEvent(time, s, sp, spp));

        logWeightContrib -= sp == spp
                ? Math.log(0.5*(trajectory.currentState[sp]*(trajectory.currentState[sp] - 1.0)))
                : Math.log(trajectory.currentState[sp]*trajectory.currentState[spp]);

        return logWeightContrib;
    }
}
