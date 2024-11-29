/*
 * Copyright (C) 2019-2024 Tim Vaughan
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package bdmmprime.trajectories;

import bdmmprime.parameterization.Parameterization;
import bdmmprime.trajectories.obsevents.ObservedEvent;
import bdmmprime.trajectories.trajevents.*;
import beast.base.util.Randomizer;

import static bdmmprime.util.Utils.*;

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
    double[][] a_migration, a_crossbirth2;
    double[][][] a_crossbirth3;

    boolean useTauLeaping;
    int minLeapCount;
    double epsilon;

    public Particle(Parameterization param, double[] initialState, boolean useTauLeaping, int minLeapCount,
                    double epsilon) {
        this.param = param;
        nTypes = param.getNTypes();
        this.useTauLeaping = useTauLeaping;
        this.minLeapCount = minLeapCount;
        this.epsilon = epsilon;

        a_birth = new double[nTypes];
        a_death = new double[nTypes];
        a_migration = new double[nTypes][nTypes];
        a_crossbirth2 = new double[nTypes][nTypes];
        a_crossbirth3 = new double[nTypes][nTypes][nTypes];

        trajectory = new Trajectory(initialState);
        logWeight = 0.0;
    }

    /**
     * Assign trajectory from another particle and set weight to zero.
     * Used during particle ensemble resampling.
     *
     * @param other Particle to assign from
     */
    public void assignTrajAndZeroWeight(Particle other) {
        trajectory.assignFrom(other.trajectory);
        logWeight = 0.0;
    }

    double t;
    int interval;

    public void propagateParticle(double tStart, int intervalStart, ObservedEvent observedEvent) {
        t = tStart;
        interval = intervalStart;

        if (logWeight == Double.NEGATIVE_INFINITY)
            return;

        if (!trajectory.currentStateValid(observedEvent.lineages))
            throw new IllegalStateException("Particle state incompatible with next observation.");

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
            if (lessThanWithPrecision(param.getIntervalEndTimes()[interval], observedEvent.time)) {
                // Include probability of seeing no rho-samples:

                if (greaterThanWithPrecision(param.getIntervalEndTimes()[interval], tStart)) {
                    for (int s = 0; s < nTypes; s++) {
                        if (param.getRhoValues()[interval][s] > 0.0)
                            logWeight += trajectory.currentState[s]
                                    * Math.log(1.0 - param.getRhoValues()[interval][s]);
                    }
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

                // Cross birth 2
                a_temp = trajectory.currentState[s]*param.getCrossBirthRates2()[interval][s][sp];
                if (a_temp > 0.0) {
                    // The following probability is for _any_ observable event produced as a result
                    // of a cross-birth, either a type change or a coalescence.
                    // This is why we don't include a factor lineages[s]/currentState[s].
                    p_obs = observedEvent.lineages[sp] / (trajectory.currentState[sp] + 1.0);
                    a_crossbirth2[s][sp] = a_temp * (1.0 - p_obs);
                    a_tot += a_crossbirth2[s][sp];
                    a_illegal_tot += a_temp * p_obs;
                } else {
                    a_crossbirth2[s][sp] = 0.0;
                }

                if (param.hasCrossBirthRates3()) {
                    for (int spp=0; spp<=sp; spp++) {
                        if (spp == s)
                            continue;

                        // Cross birth 3
                        a_temp = trajectory.currentState[s]*param.getCrossBirthRates3()[interval][s][sp][spp];
                        if (a_temp > 0.0) {
                            // Like for crossbirth2, the following probability is for _any_ observable
                            // event produced as a result of a  cross-birth.
                            if (sp == spp) {
                                p_obs = observedEvent.lineages[sp] / (trajectory.currentState[sp] + 2.0);
                            } else {
                                p_obs = observedEvent.lineages[sp] / (trajectory.currentState[sp] + 1.0)
                                        + observedEvent.lineages[spp] / (trajectory.currentState[spp] + 1.0)
                                        - (observedEvent.lineages[sp] / (trajectory.currentState[sp] + 1.0))
                                        * (observedEvent.lineages[spp] / (trajectory.currentState[spp] + 1.0));
                            }
                            a_crossbirth3[s][sp][spp] = a_temp * (1.0 - p_obs);
                            a_tot += a_crossbirth3[s][sp][spp];
                            a_illegal_tot += a_temp * p_obs;
                        } else {
                            a_crossbirth3[s][sp][spp] = 0.0;
                        }

                    }
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

                    if (u < a_crossbirth2[s][sp]) {
                        event = new CrossBirthEvent2(t, s, sp);
                        break;
                    }
                    u -= a_crossbirth2[s][sp];

                    if (param.hasCrossBirthRates3()) {
                        for (int spp=0; spp<=sp; spp++) {
                            if (spp == s)
                                continue;

                            if (u < a_crossbirth3[s][sp][spp]) {
                                event = new CrossBirthEvent3(t, s, sp, spp);
                                break;
                            }
                            u -= a_crossbirth3[s][sp][spp];
                        }
                        if (event != null)
                            break;
                    }
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

    double[] mu, sigma2;
    /**
     * Estimate tau-leaping step size for a given epsilon using the approach
     * from Cao et al., JCP 124, 044109 (2006), doi:10.1063/1.2159468
     *
     * @return tau
     */
    public double getTau() {

        double tau = Double.POSITIVE_INFINITY;

       if (mu == null)
           mu = new double[nTypes];

       if (sigma2 == null)
           sigma2 = new double[nTypes];

       for (int i=0; i<nTypes; i++) {
           mu[i] = a_birth[i] - a_death[i];
           sigma2[i] = a_birth[i] + a_death[i];

           for (int j = 0; j < nTypes; j++) {
               if (j == i)
                   continue;

               mu[i] -= a_migration[i][j];
               mu[j] += a_migration[i][j] + a_crossbirth2[i][j];

               sigma2[i] += a_migration[i][j];
               sigma2[j] += a_migration[i][j] + a_crossbirth2[i][j];

               if (param.hasCrossBirthRates3()) {
                   for (int k = 0; k <= j; k++) {
                       if (k == i)
                           continue;

                       if (k==j) {
                           mu[i] -= a_crossbirth3[i][j][j];
                           mu[j] += 2*a_crossbirth3[i][j][j];

                           sigma2[i] += a_crossbirth3[i][j][j];
                           sigma2[j] += 4*a_crossbirth3[i][j][j]; // (nu=2, so nu^2 =4)

                       } else {
                           mu[i] -= a_crossbirth3[i][j][k];
                           mu[j] += a_crossbirth3[i][j][k];
                           mu[k] += a_crossbirth3[i][j][k];

                           sigma2[i] += a_crossbirth3[i][j][k];
                           sigma2[j] += a_crossbirth3[i][j][k];
                           sigma2[k] += a_crossbirth3[i][j][k];
                       }
                   }
               }
           }
       }

       for (int i=0; i<nTypes; i++) {
           if (mu[i] != 0)
               tau = Math.min(tau, Math.max(epsilon*trajectory.currentState[i], 1)/Math.abs(mu[i]));

           if (sigma2[i] >= 0)
               tau = Math.min(tau, Math.pow(Math.max(epsilon*trajectory.currentState[i], 1),2)/sigma2[i]);
       }

       return Math.min(tau, param.getTotalProcessLength()/minLeapCount);
    }

    public void stepParticleTauLeaping(ObservedEvent observedEvent, double tmax) {

        while (t < tmax) {

            computePropensities(observedEvent);

            double tnext = Math.min(tmax, t + getTau());
            double dt = tnext - t;
            t = tnext;

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

                    if (a_crossbirth2[s][sp] > 0) {
                        int nCrossBirths = (int)Randomizer.nextPoisson(a_crossbirth2[s][sp]*dt);
                        if (nCrossBirths > 0)
                            trajectory.addEvent(new CrossBirthEvent2(t, s, sp, nCrossBirths));
                    }

                    if (param.hasCrossBirthRates3()) {
                        for (int spp=0; spp<=sp; spp++) {
                            if (spp == s)
                                continue;

                            if (a_crossbirth3[s][sp][spp] > 0) {
                                int nCrossBirths = (int)Randomizer.nextPoisson(a_crossbirth3[s][sp][spp]*dt);
                                if (nCrossBirths>0)
                                    trajectory.addEvent((new CrossBirthEvent3(t, s, sp, spp, nCrossBirths)));
                            }
                        }
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
