package bdmmprime.trajectories.obsevents;

import bdmmprime.parameterization.Parameterization;
import bdmmprime.trajectories.Trajectory;
import bdmmprime.trajectories.trajevents.BirthEvent;
import bdmmprime.trajectories.trajevents.MigrationEvent;
import beast.util.Randomizer;

public class TypeChangeEvent extends ObservedEvent {

    public int childType;

    public TypeChangeEvent(double time, int parentType, int childType, int multiplicity) {
        this.time = time;
        this.type = parentType;
        this.childType = childType;
    }

    @Override
    public int[] getNextLineageCounts() {
        int[] nextLineageCounts = lineages.clone();
        nextLineageCounts[type] -= 1;
        nextLineageCounts[childType] += 1;

        return nextLineageCounts;
    }

    double[][] birth_props;
    @Override
    public double applyToTrajectory(Parameterization param, int interval, Trajectory trajectory) {

        if (birth_props == null)
            birth_props = new double[param.getNTypes()][param.getNTypes()];

        double logWeightContrib = 0;

        int s = type;
        int sc = childType;

        double migration_prop = trajectory.currentState[s] * param.getMigRates()[interval][s][sc];
        double birth_prop_tot = 0.0;
        for (int sc_other=0; sc_other < param.getNTypes(); sc_other++) {
            int sp, spp;
            if (sc_other <= sc) {
                sp = sc;
                spp = sc_other;
            } else {
                sp = sc_other;
                spp = sc;
            }
            birth_props[sp][spp] = trajectory.currentState[s] * param.getBirthRates()[interval][s][sp][spp];
            birth_prop_tot += birth_props[sp][spp];
        }

        logWeightContrib += Math.log(migration_prop + birth_prop_tot);

        double u;

        boolean isMigration;
        if (migration_prop == 0.0)
            isMigration = false;
        else if (birth_prop_tot == 0.0)
            isMigration = true;
        else {
            u = Randomizer.nextDouble() * (migration_prop + birth_prop_tot);
            isMigration = u < migration_prop;
        }

        if (isMigration) {
            logWeightContrib += -Math.log(trajectory.currentState[sc] + 1);
            trajectory.addEvent(new MigrationEvent(time, s, sc));
            return logWeightContrib;

        } else {
            u = -migration_prop;

            for (int sc_other = 0; sc_other < param.getNTypes(); sc_other++) {
                int sp, spp;
                if (sc_other <= sc) {
                    sp = sc;
                    spp = sc_other;
                } else {
                    sp = sc_other;
                    spp = sc;
                }

                if (sc_other == param.getNTypes() || u < birth_props[sp][spp]) {
                    trajectory.addEvent(new BirthEvent(time, s, sp, spp));
                    logWeightContrib +=
                            Math.log(1.0 - lineages[sc_other] / trajectory.currentState[sc_other])
                            - Math.log(trajectory.currentState[sc]);
                    return logWeightContrib;
                }
                u -= birth_props[sp][spp];
            }

            throw new IllegalStateException("Type change event selection reaction fell through.");
        }
    }
}
