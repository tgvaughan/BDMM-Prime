package bdmmprime.trajectories.simulation;

import bdmmprime.parameterization.*;
import bdmmprime.trajectories.*;
import bdmmprime.trajectories.trajevents.*;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;
import beast.math.Binomial;
import beast.util.Randomizer;
import cern.jet.random.engine.RandomEngine;

import static org.apache.commons.math3.stat.StatUtils.sum;

/**
 * Simulates a tree from a multi-type birth-death skyline process.
 *
 * Note that the time of origin is also sampled, as this is a random
 * variable that depends on the time of the most recent sample in the tree.
 */
public class SimulatedTree extends Tree {

    public Input<Parameterization> parameterizationInput = new Input<>("parameterization",
            "BDMM parameterization",
            Input.Validate.REQUIRED);

    public Input<RealParameter> frequenciesInput = new Input<>("frequencies",
            "The equilibrium frequencies for each type",
            Input.Validate.REQUIRED);

    public Input<RealParameter> simulationTimeInput = new Input<>("simulationTime",
            "Time to run simulation for.",
            Input.Validate.REQUIRED);

    Parameterization param;
    RealParameter frequencies;
    RealParameter simulationTime;

    double[] a_birth, a_death, a_sampling;
    double[][] a_migration, a_crossbirth;

    int nTypes;

    @Override
    public void initAndValidate() {
        param = parameterizationInput.get();
        frequencies = frequenciesInput.get();
        simulationTime = simulationTimeInput.get();

        nTypes = param.getNTypes();

        a_birth = new double[nTypes];
        a_death = new double[nTypes];
        a_sampling = new double[nTypes];
        a_migration = new double[nTypes][nTypes];
        a_crossbirth = new double[nTypes][nTypes];

        super.initAndValidate();
    }

    Trajectory simulateTrajectory() {
        double t = 0;
        int interval = 0;

        double[] initialState = new double[nTypes];
        int startType;
        double u = Randomizer.nextDouble()*sum(frequencies.getDoubleValues());
        for (startType=0; startType<nTypes-1; startType++) {
            if (u < frequencies.getValue(startType))
                break;
            u -= frequencies.getValue(startType);
        }
        initialState[startType] = 1.0;

        Trajectory traj = new Trajectory(initialState);
        while (true) {

            double a_tot = 0.0;
            for (int s=0; s<nTypes; s++) {
                a_birth[s] = traj.currentState[s] * param.getBirthRates()[interval][s];
                a_death[s] = traj.currentState[s] * param.getDeathRates()[interval][s];
                a_sampling[s] = traj.currentState[s] * param.getSamplingRates()[interval][s];
                a_tot += a_birth[s] + a_death[s] + a_sampling[s];

                for (int sp = 0; sp < nTypes; sp++) {
                    if (sp == s)
                        continue;

                    a_migration[s][sp] = traj.currentState[s] * param.getMigRates()[interval][s][sp];
                    a_crossbirth[s][sp] = traj.currentState[s] * param.getCrossBirthRates()[interval][s][sp];
                    a_tot += a_migration[s][sp] + a_crossbirth[s][sp];
                }
            }

            if (a_tot > 0)
                t += Randomizer.nextExponential(a_tot);
            else
                t += Double.POSITIVE_INFINITY;

            if (param.getIntervalEndTimes()[interval] < simulationTime.getValue()
                    && t > param.getIntervalEndTimes()[interval]) {
                t = param.getIntervalEndTimes()[interval];

                for (int s=0; s<nTypes; s++) {
                    double rho = param.getRhoValues()[interval][s];
                    if (rho > 0) {
                        int nRemoveSamp = nextBinomial((int)Math.round(traj.currentState[s]),
                                rho*param.getRemovalProbs()[interval][s]);

                        if (nRemoveSamp>0)
                            traj.addEvent(new SamplingWithRemoval(t, s, nRemoveSamp));

                        int nNoRemoveSamp = nextBinomial((int)Math.round(traj.currentState[s]),
                                rho*(1.0 - param.getRemovalProbs()[interval][s]));

                        if (nNoRemoveSamp>0)
                            traj.addEvent(new SamplingWithoutRemoval(t, s, nNoRemoveSamp));
                    }
                }

                interval += 1;
                continue;
            }

            if (t > simulationTime.getValue()) {
                t = simulationTime.getValue();
                break;
            }

            u = Randomizer.nextDouble()*a_tot;

            TrajectoryEvent event = null;
            for (int s=0; s<nTypes; s++) {

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

                if (u < a_sampling[s]) {
                    if (u < param.getRemovalProbs()[interval][s]*a_sampling[s])
                        event = new SamplingWithRemoval(t, s);
                    else
                        event = new SamplingWithoutRemoval(t, s);
                    break;
                }
                u -= a_sampling[s];

                for (int sp=0; sp<nTypes; sp++) {
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

            if (event == null)
                throw new IllegalStateException("Event selection loop fell through.");

            traj.addEvent(event);
        }

        return traj;
    }

    public int nextBinomial(int n, double p) {
        double u = Randomizer.nextDouble();

        double acc = Math.pow(1-p, n);
        int m = 0;

        while (u > acc && m <= n) {
            m += 1;
            acc += Math.exp(Binomial.logChoose(n, m) + (n - m) * Math.log(1 - p) + m * Math.log(p));
        }

        return m;
    }

    /**
     * Debug method for testing
     * @param args unused
     */
    public static void main(String[] args) {

        Parameterization param = new CanonicalParameterization();

        param.initByName(
                "typeSet", new TypeSet(2),
                "origin", "3.0",
                "birthRate", new SkylineVectorParameter(
                        new RealParameter("4.0"),
                        new RealParameter("2.0 0.0"), 2),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0"), 2),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.1"), 2),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0"), 2),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.5"), 2),
                "rhoSampling", new TimedParameter(
                        new RealParameter("4.0"),
                        new RealParameter("0.5 0.0"), 2));

        SimulatedTree sim = new SimulatedTree();
        sim.initByName("parameterization", param,
                "frequencies", new RealParameter("0.5 0.5"),
                "simulationTime", "5.0");

        Trajectory traj = sim.simulateTrajectory();
        traj.dump();

    }

}
