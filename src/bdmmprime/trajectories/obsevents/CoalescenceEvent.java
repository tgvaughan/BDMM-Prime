package bdmmprime.trajectories.obsevents;

import bdmmprime.parameterization.Parameterization;
import bdmmprime.trajectories.Trajectory;
import bdmmprime.trajectories.trajevents.BirthEvent;
import bdmmprime.trajectories.trajevents.CrossBirthEvent;

public class CoalescenceEvent extends ObservedEvent {

    public int childType1, childType2;

    public CoalescenceEvent(double time, int parentType, int childType1, int childType2, int multiplicity) {
        this.time = time;
        this.type = parentType;
        this.childType1 = childType1;
        this.childType2 = childType2;
        this.multiplicity = multiplicity;
    }

    @Override
    public int[] getNextLineageCounts() {
        int[] nextLineageCounts = lineages.clone();

        nextLineageCounts[type] -= multiplicity;
        nextLineageCounts[childType1] += multiplicity;
        nextLineageCounts[childType2] += multiplicity;

        return nextLineageCounts;
    }

    @Override
    public double applyToTrajectory(Parameterization param, int interval, Trajectory trajectory) {
        double logWeightContrib = 0.0;

        int s = type;
        int sc1 = childType1;
        int sc2 = childType2;

        int sp = sc1 != s ? sc1 : sc2;

        if (sp == s) {
            // Birth

            double birth_prop = trajectory.currentState[s]*param.getBirthRates()[interval][s];
            logWeightContrib += Math.log(birth_prop)
                    - Math.log(0.5*(trajectory.currentState[s]*(trajectory.currentState[s] + 1.0)));

            trajectory.addEvent(new BirthEvent(time, s));

        } else {
            // Cross-birth

            double crossbirth_prop = trajectory.currentState[s]*param.getCrossBirthRates()[interval][s][sp];
            logWeightContrib += Math.log(crossbirth_prop)
                    - Math.log(trajectory.currentState[s]*(trajectory.currentState[sp] + 1.0));

            trajectory.addEvent(new CrossBirthEvent(time, s, sp));
        }

        return logWeightContrib;
    }
}
