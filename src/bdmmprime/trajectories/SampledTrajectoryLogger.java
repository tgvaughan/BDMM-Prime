package bdmmprime.trajectories;

import bdmmprime.parameterization.Parameterization;
import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Loggable;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

public class SampledTrajectoryLogger extends BEASTObject implements Loggable {

    public Input<Tree> mappedTreeInput = new Input<>("typeMappedTree",
            "Tree with stochastically mapped types.", Input.Validate.REQUIRED);

    public Input<Parameterization> parameterizationInput = new Input<>("parametermization",
            "Multi-type birth-death parameterization.", Input.Validate.REQUIRED);

    public Input<Integer> nParticlesInput = new Input<>("nParticles",
            "Number of particles to use in filtering calculation.", 1000);

    public Input<String> typeLabelInput = new Input<>("typeLabel",
            "Type label used for traits in generated metadata.",
            Input.Validate.REQUIRED);

    Tree mappedTree;
    String typeLabel;
    Parameterization param;
    int nTypes, nParticles;

    double[] a_birth, a_death, a_sampling;
    double[][] a_migration, a_crossbirth;
    @Override
    public void initAndValidate() {

        mappedTree = mappedTreeInput.get();
        typeLabel = typeLabelInput.get();
        param = parameterizationInput.get();
        nTypes = param.getNTypes();
        nParticles = nParticlesInput.get();

        a_birth = new double[nTypes];
        a_death = new double[nTypes];
        a_sampling = new double[nTypes];
        a_migration = new double[nTypes][nTypes];
        a_crossbirth = new double[nTypes][nTypes];
    }

    public Trajectory simulateTrajectory() {
        Node[] treeEvents = new Node[mappedTree.getNodeCount()];
        System.arraycopy(mappedTree.getNodesAsArray(), 0, treeEvents, 0, treeEvents.length);
        Arrays.sort(treeEvents, Comparator.comparingDouble(Node::getHeight).reversed());

        int rootType = (int)treeEvents[0].getMetaData(typeLabel);
        int[][] lineageCounts = new int[treeEvents.length][nTypes];
        for (int s=0; s<nTypes; s++)
            lineageCounts[0][s] = s == rootType ? 1 : 0;

        for (int i=1; i<treeEvents.length; i++) {
            System.arraycopy(lineageCounts[i-1], 0, lineageCounts[i], 0, nTypes);
            lineageCounts[i][(int)treeEvents[i-1].getMetaData(typeLabel)] -= 1;
            for (Node child : treeEvents[i-1].getChildren())
                lineageCounts[i][(int)child.getMetaData(typeLabel)] += 1;
        }

        // Initialize particles
        double[] initialState = new double[param.getNTypes()];
        initialState[rootType] = 1.0;

        Trajectory[] particleTrajectories = new Trajectory[nParticles];
        double[] logParticleWeights = new double[nParticles];
        double[] particleWeights = new double[nParticles];
        for (int p=0; p<nParticles; p++) {
            particleTrajectories[p] = new Trajectory(initialState);
            logParticleWeights[p] = 0.0;
        }

        Trajectory[] particleTrajectoriesPrime = new Trajectory[nParticles];

        // Iterate over tree events:

        double t = 0.0;

        for (int i=0; i<treeEvents.length; i++) {

            // Propagate particles to next event

            double maxLogWeight = Double.NEGATIVE_INFINITY;
            for (int p=0; p<nParticles; p++) {
                logParticleWeights[p] = propagateParticle(particleTrajectories[p], lineageCounts[i], t, treeEvents[i]);

                if (logParticleWeights[p] > maxLogWeight)
                    maxLogWeight = logParticleWeights[p];
            }

            if (maxLogWeight == Double.NEGATIVE_INFINITY) {
                throw new IllegalStateException("Particle ensemble depleted.");
            }

            // Compute sum of scaled weights:
            double sumOfScaledWeights = 0.0;
            for (int p=0; p<nParticles; p++) {
                particleWeights[p] = Math.exp(logParticleWeights[p]-maxLogWeight);
                sumOfScaledWeights += particleWeights[p];
            }

            // Normalize weights:
            for (int p=0; p<nParticles; p++)
                particleWeights[p] = particleWeights[p]/sumOfScaledWeights;

            // Resample particle ensemble

            ReplacementSampler sampler = new ReplacementSampler(particleWeights);
            for (int p=0; p<nParticles; p++)
                particleTrajectoriesPrime[p] = particleTrajectories[sampler.next()].copy();

            Trajectory[] tmp = particleTrajectories;
            particleTrajectories = particleTrajectoriesPrime;
            particleTrajectoriesPrime = tmp;
        }

        // WLOG choose 0th particle as the particle to return:

        return particleTrajectories[0];
    }

    double propagateParticle(Trajectory trajectory, int[] lineageCounts, double t, Node nextTreeEvent) {
        double logWeight = 0.0;

        int interval = param.getIntervalIndex(t);

        while (true) {
            // Compute rates

            double a_temp;
            double a_tot = 0.0;
            double a_illegal_tot = 0.0;

            for (int s=0; s<nTypes; s++) {
                a_temp = trajectory.currentState[s]*param.getBirthRates()[interval][s];
                double p_obs = lineageCounts[s]*(lineageCounts[s]-1.0)/(trajectory.currentState[s]*(trajectory.currentState[s]+1));

                a_birth[s] = a_temp*(1-p_obs);
                a_tot += a_birth[s];
                a_illegal_tot += a_temp*p_obs;

                a_temp = trajectory.currentState[s] * param.getDeathRates()[interval][s];
                if (trajectory.currentState[s] > lineageCounts[s]) {
                    a_death[s] = a_temp;
                    a_tot += a_temp;
                } else {
                    a_death[s] = 0.0;
                    a_illegal_tot += a_temp;
                }

                a_illegal_tot += a_sampling[s];

                for (int sp=0; sp<nTypes; sp++) {
                    if (sp == s)
                        continue;

                    a_temp = trajectory.currentState[s]*param.getMigRates()[interval][s][sp];
                    if (trajectory.currentState[s] > lineageCounts[s]) {
                        a_migration[s][sp] = a_temp;
                        a_tot += a_temp;
                    } else {
                        a_migration[s][sp] = 0.0;
                        a_illegal_tot += a_temp;
                    }

                    a_crossbirth[s][sp] = trajectory.currentState[s]*param.getCrossBirthRates()[interval][s][sp];
                    a_tot += a_crossbirth[s][sp];
                }
            }

            // Sample time

            double tprime;
            if (a_tot > 0.0)
                tprime = t + Randomizer.nextExponential(a_tot);
            else
                tprime = Double.POSITIVE_INFINITY;

            // Test time

            if (param.getIntervalEndTimes()[interval] < param.getNodeTime(nextTreeEvent)) {
                if (tprime > param.getIntervalEndTimes()[interval]) {
                    logWeight += -a_illegal_tot*(tprime - t);
                    t = param.getIntervalEndTimes()[interval];
                    continue;
                }
            } else {
                if (t > param.getNodeTime(nextTreeEvent)) {
                    logWeight += -a_illegal_tot*(tprime - t);
                    t = param.getNodeTime(nextTreeEvent);
                    break;
                }
            }

            // Update weight
            logWeight += -a_illegal_tot*(tprime - t);

            // Implement event

            double u = Randomizer.nextDouble()*a_tot;

            TrajectoryEvent event = null;
            for (int s = 0; s<nTypes; s++) {
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

            if (event == null) {
                throw new IllegalStateException("Event selection loop fell through.");
            }

            trajectory.addEvent(event);
        }

        // Compute tree event probability



        return logWeight;
    }

    @Override
    public void init(PrintStream out) {

    }

    @Override
    public void log(long sample, PrintStream out) {

    }

    @Override
    public void close(PrintStream out) {

    }
}
