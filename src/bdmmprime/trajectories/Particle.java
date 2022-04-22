package bdmmprime.trajectories;

import bdmmprime.parameterization.Parameterization;
import bdmmprime.trajectories.obsevents.ObservedEvent;
import bdmmprime.trajectories.trajevents.*;
import beast.util.Randomizer;

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
    double[] a_death;
    double[][] a_migration;
    double[][][] a_birth;

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

        a_birth = new double[nTypes][nTypes][nTypes];
        a_death = new double[nTypes];
        a_migration = new double[nTypes][nTypes];

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

    /**
     *  Compute propensities for particle simulation.
     *
     * @param observedEvent Next observed event object (contains details of
     *                      extant lineages in this interval)
     */
    public void computePropensities (ObservedEvent observedEvent) {

        double a_temp, p_obs;
        a_tot = 0.0;
        a_illegal_tot = 0.0;

        for (int s=0; s<nTypes; s++) {

            int ks = observedEvent.lineages[s];
            double Ns = trajectory.currentState[s];

            a_temp = Ns * param.getDeathRates()[interval][s];
            if (Ns > ks) {
                a_death[s] = a_temp;
                a_tot += a_death[s];
            } else {
                a_death[s] = 0.0;
                a_illegal_tot += a_temp;
            }

            a_illegal_tot += Ns*param.getSamplingRates()[interval][s];

            for (int sp=0; sp<nTypes; sp++) {

                int ksp = observedEvent.lineages[sp];
                double Nsp = trajectory.currentState[sp];

                if (sp != s) {
                    // Migration
                    a_temp = Ns * param.getMigRates()[interval][s][sp];
                    if (Ns > ks) {
                        p_obs = ksp / (Nsp + 1.0);
                        a_migration[s][sp] = a_temp * (1.0 - p_obs);
                        a_tot += a_migration[s][sp];
                        a_illegal_tot += a_temp * p_obs;
                    } else {
                        a_migration[s][sp] = 0.0;
                        a_illegal_tot += a_temp;
                    }
                }

                for (int spp=0; spp<=sp; spp++) {
                    // Birth
                    a_temp = Ns * param.getBirthRates()[interval][s][sp][spp];
                    if (a_temp > 0.0) {

                        int kspp = observedEvent.lineages[spp];
                        double Nspp = trajectory.currentState[spp];

                        if (sp == s && spp == s) {

                            // Only observable event is coalescence
                            p_obs = ks*(ks-1.0)/((Ns+1.0)*Ns);
                        } else if (sp == s) {
                            // sp == s and spp != s

                            // Need lineage in spp to be involved for both type change and
                            // for coalescence
                            p_obs = kspp/(Nspp+1.0);

                        } else if (spp == s) {
                            // sp != s and spp == s

                            // Need lineage in sp to be involved for both type change and
                            // for coalescence
                            p_obs = ksp/(Nsp+1.0);

                        } else if (sp == spp) {
                            // sp != s and spp !=s but sp == spp

                            // Need at least one lineage in either sp or spp to
                            // be involved.  (I.e. 1-p_neither.)  Note that since
                            // sp != s, Nsp is incremented by 2 by this event.

                            p_obs = ksp/(Nsp+2.0) + ksp/(Nsp+1.0) - ksp*ksp/((Nsp+2.0)*(Nsp+1.0));

                        } else {
                            // sp != s and spp != s and sp != spp

                            // Need at least one lineage in either sp or spp to
                            // be involved. (I.e. 1-p_neither.)


                            p_obs = ksp/(Nsp+1.0) + kspp/(Nspp+1.0) - ksp*kspp/((Nsp+1.0)*(Nspp+1.0));
                        }

                        a_birth[s][sp][spp] = a_temp * (1.0 - p_obs);
                        a_tot += a_birth[s][sp][spp];
                        a_illegal_tot += a_temp * p_obs;
                    } else {
                        a_birth[s][sp][spp] = 0.0;
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

                if (u < a_death[s]) {
                    event = new DeathEvent(t, s);
                    break;
                }
                u -= a_death[s];

                for (int sp = 0; sp < nTypes; sp++) {

                    if (sp != s) {
                        if (u < a_migration[s][sp]) {
                            event = new MigrationEvent(t, s, sp);
                            break;
                        }
                        u -= a_migration[s][sp];
                    }

                    for (int spp=0; spp<=sp; spp++) {
                        if (u < a_birth[s][sp][spp]) {
                            event = new BirthEvent(t, s, sp, spp);
                            break;
                        }
                        u -= a_birth[s][sp][spp];
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
     * @return suggested tau
     */
    public double getTau() {

        if (mu == null)
            mu = new double[nTypes];

        if (sigma2 == null)
            sigma2 = new double[nTypes];

        for (int i=0; i<nTypes; i++) {
            mu[i] = -a_death[i];
            sigma2[i] = a_death[i];

            for (int j = 0; j < nTypes; j++) {
                if (i != j) {
                    mu[i] += a_migration[i][j];
                    sigma2[i] += a_migration[i][j];

                    mu[j] -= a_migration[i][j];
                    sigma2[j] += a_migration[i][j];
                }

                for (int k = 0; k <= j; k++) {
                    if (j == i && k == i) {
                        mu[i] += a_birth[i][i][i];
                        sigma2[i] += a_birth[i][i][i];
                    } else if (j == i) {
                        mu[k] += a_birth[i][i][k];
                        sigma2[k] += a_birth[i][i][k];
                    } else if (k == i) {
                        mu[j] += a_birth[i][j][i];
                        sigma2[j] += a_birth[i][j][i];
                    } else if (j == k) {
                        mu[i] -= a_birth[i][j][j];
                        sigma2[i] += a_birth[i][j][j];
                        mu[j] += 2 * a_birth[i][j][j];
                        sigma2[j] += 4 * a_birth[i][j][j];
                    } else {
                        mu[i] -= a_birth[i][j][k];
                        sigma2[i] += a_birth[i][j][k];
                        mu[j] += a_birth[i][j][k];
                        sigma2[j] += a_birth[i][j][k];
                        mu[k] += a_birth[i][j][k];
                        sigma2[k] += a_birth[i][j][k];
                    }
                }
            }
        }

        double tau = Double.POSITIVE_INFINITY;

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
                if (a_death[s] > 0) {
                    int nDeaths = (int)Randomizer.nextPoisson(a_death[s]*dt);
                    if (nDeaths > 0)
                        trajectory.addEvent(new DeathEvent(t, s, nDeaths));
                }

                for (int sp = 0; sp < nTypes; sp++) {
                    if (sp != s) {
                        if (a_migration[s][sp] > 0) {
                            int nMigs = (int) Randomizer.nextPoisson(a_migration[s][sp] * dt);
                            if (nMigs > 0)
                                trajectory.addEvent(new MigrationEvent(t, s, sp, nMigs));
                        }
                    }

                    for (int spp=0; spp<=sp; spp++) {
                        if (a_birth[s][sp][spp] > 0) {
                            int nBirths = (int) Randomizer.nextPoisson(a_birth[s][sp][spp] * dt);
                            if (nBirths > 0)
                                trajectory.addEvent(new BirthEvent(t, s, sp, spp, nBirths));
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
