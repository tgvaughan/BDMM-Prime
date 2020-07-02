package bdmmprime.trajectories.obsevents;

import bdmmprime.parameterization.Parameterization;
import bdmmprime.trajectories.Trajectory;
import bdmmprime.trajectories.trajevents.CrossBirthEvent;
import bdmmprime.trajectories.trajevents.MigrationEvent;
import beast.util.Randomizer;

public class TypeChangeEvent extends ObservedEvent {

    public int childType;

    public TypeChangeEvent(double time, int parentType, int childType, int multiplicity) {
        this.time = time;
        this.type = parentType;
        this.childType = childType;
        this.multiplicity = multiplicity;
    }

    @Override
    public int[] getNextLineageCounts() {
        int[] nextLineageCounts = lineages.clone();
        nextLineageCounts[type] -= multiplicity;
        nextLineageCounts[childType] += multiplicity;

        return nextLineageCounts;
    }

    @Override
    public double applyToTrajectory(Parameterization param, int interval, Trajectory trajectory) {

        double logWeightContrib = 0;

        int s = type;
        int sp = childType;

        for (int i=0; i<multiplicity; i++) {
            double migration_prop = trajectory.currentState[s] * param.getMigRates()[interval][s][sp];
            double crossbirth_prop = trajectory.currentState[s] * param.getCrossBirthRates()[interval][s][sp];

            logWeightContrib += Math.log(migration_prop + crossbirth_prop);

            boolean isMigration;
            if (migration_prop == 0.0)
                isMigration = false;
            else if (crossbirth_prop == 0.0)
                isMigration = true;
            else {
                isMigration = Randomizer.nextDouble() * (migration_prop + crossbirth_prop) < migration_prop;
            }

            if (isMigration) {
                logWeightContrib += -Math.log(trajectory.currentState[sp] + 1);
                trajectory.addEvent(new MigrationEvent(time, s, sp));

            } else {
                logWeightContrib += Math.log(crossbirth_prop * (1.0 - lineages[s] / trajectory.currentState[s]));

                trajectory.addEvent(new CrossBirthEvent(time, s, sp));
            }
        }

        return logWeightContrib;
    }
}
