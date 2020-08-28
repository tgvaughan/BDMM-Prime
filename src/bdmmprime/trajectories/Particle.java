package bdmmprime.trajectories;

import bdmmprime.parameterization.Parameterization;
import bdmmprime.trajectories.obsevents.ObservedEvent;
import bdmmprime.trajectories.trajevents.*;
import beast.util.Randomizer;

import java.util.Random;

/**
 * This class wraps up the particle trajectories with the state required to simulate their evolution.
 *
 * It's pretty ugly, honestly.  The simulation algorithm itself is composed
 * of a bunch of side-effect-only methods that (ab)use fields to hold the shared state.
 *
 * However, this is the structure I find easiest to debug and maintain.
 */
public class Particle {

    Parameterization param;
    int nTypes;

    Trajectory trajectory;
    double logWeight;

    /**
     * Propensities
     */
    double a_tot, a_illegal_tot;
    double[] a_birth, a_death;
    double[][] a_migration, a_crossbirth;

    boolean useTauLeaping;
    int stepsPerInterval;

    public Particle(Parameterization param, double[] initialState, boolean useTauLeaping, int stepsPerInterval) {
        this.param = param;
        nTypes = param.getNTypes();
        this.useTauLeaping = useTauLeaping;
        this.stepsPerInterval = stepsPerInterval;

        a_birth = new double[nTypes];
        a_death = new double[nTypes];
        a_migration = new double[nTypes][nTypes];
        a_crossbirth = new double[nTypes][nTypes];

        trajectory = new Trajectory(initialState);
        logWeight = 0.0;
    }

    /**
     * Assign trajectory from another particle and set weight to zero.
     * Used during particle ensemble resampling.
     *
     * @param other Particle to assign from
     */
    public void assignFrom(Particle other) {
        trajectory.assignFrom(other.trajectory);
        logWeight = 0.0;
    }

    double t;
    int interval;

    public void propagateParticle(double tStart, int intervalStart, ObservedEvent observedEvent) {
        t = tStart;
        interval = intervalStart;

        if (!trajectory.currentStateValid(observedEvent.lineages))
            throw new IllegalStateException("Particle state incompatible with next observation.");


        if (logWeight == Double.NEGATIVE_INFINITY)
           return;

        while (true) {

            // Step particle
            double tmax = Math.min(param.getIntervalEndTimes()[interval], observedEvent.time);

            if (useTauLeaping)
                stepParticleTauLeaping(observedEvent, tmax);
            else
                stepParticleGillespie(observedEvent, tmax);

            if (logWeight == Double.NEGATIVE_INFINITY)
                return;

            if (!trajectory.currentStateValid(observedEvent.lineages))
                throw new IllegalStateException("Unobserved event produced illegal state.");

            // Test for end of interval or simulation
            if (param.getIntervalEndTimes()[interval] < observedEvent.time) {
                // Include probability of seeing no rho-samples:
                for (int s = 0; s < nTypes; s++) {
                    if (param.getRhoValues()[interval][s] > 0.0)
                        logWeight += trajectory.currentState[s]
                                * Math.log(1.0 - param.getRhoValues()[interval][s]);
                }

                interval += 1;
            } else {
                break;
            }
        }

        // Compute tree event contribution
        if (!observedEvent.isFinalEvent())
            logWeight += observedEvent.applyToTrajectory(param, interval, trajectory);

        if (!trajectory.currentStateValid())
            throw new IllegalStateException("Observed event produced illegal state.");
    }


    public void computePropensities (ObservedEvent observedEvent) {

        double a_temp, p_obs;
        a_tot = 0.0;
        a_illegal_tot = 0.0;

        for (int s=0; s<nTypes; s++) {
            a_temp = trajectory.currentState[s]*param.getBirthRates()[interval][s];
            if (a_temp > 0) {
                p_obs = observedEvent.lineages[s] * (observedEvent.lineages[s] - 1.0)
                        / (trajectory.currentState[s] * (trajectory.currentState[s] + 1.0));
                a_birth[s] = a_temp * (1 - p_obs);
                a_tot += a_birth[s];
                a_illegal_tot += a_temp * p_obs;
            } else {
                a_birth[s] = 0.0;
            }

            a_temp = trajectory.currentState[s] * param.getDeathRates()[interval][s];
            if (trajectory.currentState[s] > observedEvent.lineages[s]) {
                a_death[s] = a_temp;
                a_tot += a_death[s];
            } else {
                a_death[s] = 0.0;
                a_illegal_tot += a_temp;
            }

            a_illegal_tot += trajectory.currentState[s]*param.getSamplingRates()[interval][s];

            for (int sp=0; sp<nTypes; sp++) {
                if (sp == s)
                    continue;

                // Migration
                a_temp = trajectory.currentState[s] * param.getMigRates()[interval][s][sp];
                if (trajectory.currentState[s]>observedEvent.lineages[s]) {
                    p_obs = observedEvent.lineages[sp] / (trajectory.currentState[sp] + 1.0);
                    a_migration[s][sp] = a_temp * (1.0 - p_obs);
                    a_tot += a_migration[s][sp];
                    a_illegal_tot += a_temp * p_obs;
                } else {
                    a_migration[s][sp] = 0.0;
                    a_illegal_tot += a_temp;
                }

                // Cross birth
                a_temp = trajectory.currentState[s]*param.getCrossBirthRates()[interval][s][sp];
                if (a_temp > 0.0) {
                    // The following probability is for _any_ observable event produced as a result
                    // of a cross-birth, either a type change or a coalescence.
                    // This is why we don't include a factor lineages[s]/currentState[s].
                    p_obs = observedEvent.lineages[sp] / (trajectory.currentState[sp] + 1.0);
                    a_crossbirth[s][sp] = a_temp * (1.0 - p_obs);
                    a_tot += a_crossbirth[s][sp];
                    a_illegal_tot += a_temp * p_obs;
                } else {
                    a_crossbirth[s][sp] = 0.0;
                }
            }
        }
    }

    public void stepParticleGillespie(ObservedEvent observedEvent, double tmax) {

        while (true) {
            // Compute rates
            computePropensities(observedEvent);

            double tprime;
            if (a_tot > 0.0)
                tprime = t + Randomizer.nextExponential(a_tot);
            else
                tprime = Double.POSITIVE_INFINITY;

            double tnew = Math.min(tprime, tmax);

            // Update weight and time

            logWeight += -a_illegal_tot * (tnew - t);
            t = tnew;

            if (tprime > tmax)
                return;

            // Implement event

            double u = Randomizer.nextDouble() * a_tot;

            TrajectoryEvent event = null;
            for (int s = 0; s < nTypes; s++) {
                if (u < a_birth[s]) {
                    event = new BirthEvent(t, s);
                    break;
                }
                u -= a_birth[s];

                if (u < a_death[s]) {
                    event = new DeathEvent(t, s);
                    break;
                }
                u -= a_death[s];

                for (int sp = 0; sp < nTypes; sp++) {
                    if (sp == s)
                        continue;

                    if (u < a_migration[s][sp]) {
                        event = new MigrationEvent(t, s, sp);
                        break;
                    }
                    u -= a_migration[s][sp];

                    if (u < a_crossbirth[s][sp]) {
                        event = new CrossBirthEvent(t, s, sp);
                        break;
                    }
                    u -= a_crossbirth[s][sp];
                }
                if (event != null)
                    break;
            }

            if (event == null) {
                throw new IllegalStateException("Event selection loop fell through.");
            }

            trajectory.addEvent(event);
        }

    }

    public void stepParticleTauLeaping(ObservedEvent observedEvent, double tmax) {

        double t0 = t;
        double dt = (tmax - t0) / stepsPerInterval;
        for (int step = 1; step <= stepsPerInterval; step += 1) {
            t = t0 + step*dt;

            computePropensities(observedEvent);

            logWeight += -dt*a_illegal_tot;

            for (int s = 0; s < nTypes; s++) {
                if (a_birth[s] > 0) {
                    int nBirths = (int)Randomizer.nextPoisson(a_birth[s]*dt);
                    if (nBirths > 0)
                        trajectory.addEvent(new BirthEvent(t, s, nBirths));
                }

                if (a_death[s] > 0) {
                    int nDeaths = (int)Randomizer.nextPoisson(a_death[s]*dt);
                    if (nDeaths > 0)
                        trajectory.addEvent(new DeathEvent(t, s, nDeaths));
                }

                for (int sp = 0; sp < nTypes; sp++) {
                    if (sp == s)
                        continue;

                    if (a_migration[s][sp] > 0) {
                        int nMigs = (int)Randomizer.nextPoisson(a_migration[s][sp]*dt);
                        if (nMigs > 0)
                            trajectory.addEvent(new MigrationEvent(t, s, sp, nMigs));
                    }

                    if (a_crossbirth[s][sp] > 0) {
                        int nCrossBirths = (int)Randomizer.nextPoisson(a_crossbirth[s][sp]*dt);
                        if (nCrossBirths > 0)
                            trajectory.addEvent(new CrossBirthEvent(t, s, sp, nCrossBirths));
                    }
                }
            }

            if (!trajectory.currentStateValid(observedEvent.lineages)) {
                logWeight = Double.NEGATIVE_INFINITY;
                return;
            }

        }
    }
}
