package bdmmprime.trajectories;

import bdmmprime.parameterization.Parameterization;
import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Loggable;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

public class SampledTrajectoryLogger extends BEASTObject implements Loggable {

    public Input<Tree> mappedTreeInput = new Input<>("typeMappedTree",
            "Tree with stochastically mapped types.", Input.Validate.REQUIRED);

    public Input<Parameterization> parameterizationInput = new Input<>("parameterization",
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
            if (treeEvents[i].isFake())
                continue;

            // Propagate particles to next event

            int interval = param.getIntervalIndex(t);

            double maxLogWeight = Double.NEGATIVE_INFINITY;
            for (int p=0; p<nParticles; p++) {
                logParticleWeights[p] = propagateParticle(particleTrajectories[p], lineageCounts[i], t, interval, treeEvents[i]);

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

            t = param.getNodeTime(treeEvents[i]);
        }

        // WLOG choose 0th particle as the particle to return:

        return particleTrajectories[0];
    }

    double propagateParticle(Trajectory trajectory, int[] lineageCounts, double t, int interval, Node nextTreeEvent) {
        double logWeight = 0.0;

        System.out.println("Event " + nextTreeEvent.getNr());

        while (true) {
            // Compute rates

            double a_temp, p_obs;
            double a_tot = 0.0;
            double a_illegal_tot = 0.0;

            for (int s=0; s<nTypes; s++) {
                a_temp = trajectory.currentState[s]*param.getBirthRates()[interval][s];
                if (a_temp > 0) {
                    p_obs = lineageCounts[s] * (lineageCounts[s] - 1.0) / (trajectory.currentState[s] * (trajectory.currentState[s] + 1));
                    a_birth[s] = a_temp * (1 - p_obs);
                    a_tot += a_birth[s];
                    a_illegal_tot += a_temp * p_obs;
                } else
                    a_birth[s] = 0.0;

                a_temp = trajectory.currentState[s] * param.getDeathRates()[interval][s];
                if (trajectory.currentState[s] > lineageCounts[s]) {
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

                    a_temp = trajectory.currentState[s]*param.getMigRates()[interval][s][sp];
                    p_obs = lineageCounts[sp]/(trajectory.currentState[sp]+1);
                    a_migration[s][sp] = a_temp*(1-p_obs);
                    a_tot += a_migration[s][sp];
                    a_illegal_tot += a_temp*p_obs;

                    a_temp = trajectory.currentState[s]*param.getCrossBirthRates()[interval][s][sp];
                    if (a_temp > 0.0) {
                        p_obs = lineageCounts[s] * lineageCounts[sp] / (trajectory.currentState[s] * (trajectory.currentState[sp] + 1));
                        a_crossbirth[s][sp] = a_temp * (1 - p_obs);
                        a_tot += a_crossbirth[s][sp];
                        a_illegal_tot += a_temp * p_obs;
                    } else {
                        a_crossbirth[s][sp] = 0.0;
                    }
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
                    // TODO: Add in treatment of rho sampling

                    logWeight += -a_illegal_tot*(tprime - t);
                    t = param.getIntervalEndTimes()[interval];
                    interval += 1;
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

            // Update time
            t = tprime;

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

        if (nextTreeEvent.getChildCount() == 0) {
            // Sample

            int s = (int)nextTreeEvent.getMetaData(typeLabel);

            double sampling_prop = trajectory.currentState[s]*param.getSamplingRates()[interval][s];
            if (nextTreeEvent.isDirectAncestor()) {
                logWeight += Math.log((1.0-param.getRemovalProbs()[interval][s])*sampling_prop);
            } else {
                logWeight += Math.log(sampling_prop);

                boolean isRemoval = Randomizer.nextDouble() < param.getRemovalProbs()[interval][s];
                if (isRemoval) {
                    trajectory.addEvent(new DeathEvent(t, s));
                } else {
                    logWeight += Math.log(1.0 - (lineageCounts[s]-1.0)/trajectory.currentState[s]);
                }
            }

            // TODO: Add in treatment of rho sampling.

        } else if (nextTreeEvent.getChildCount() == 1) {
            // Migration or cross-birth

            int s = (int)nextTreeEvent.getMetaData(typeLabel);
            int sp = (int)nextTreeEvent.getChild(0).getMetaData(typeLabel);

            double migration_prop = trajectory.currentState[s]*param.getMigRates()[interval][s][sp] ;
            double crossbirth_prop = trajectory.currentState[s]*param.getCrossBirthRates()[interval][s][sp];

            logWeight += Math.log(migration_prop + crossbirth_prop);

            boolean isMigration;
            if (migration_prop == 0.0)
                isMigration = false;
            else if (crossbirth_prop == 0.0)
                isMigration = true;
            else {
                isMigration = Randomizer.nextDouble()*(migration_prop + crossbirth_prop) < migration_prop;
            }

            if (isMigration) {
                logWeight += -Math.log(trajectory.currentState[sp]+1);
                trajectory.addEvent(new MigrationEvent(t, s, sp));

            } else {
                logWeight += Math.log(crossbirth_prop*(1.0 - lineageCounts[s]/trajectory.currentState[s]));

                trajectory.addEvent(new CrossBirthEvent(t, s, sp));
            }

        } else if (nextTreeEvent.getChildCount() == 2) {
            // Birth or cross-birth

            int s = (int)nextTreeEvent.getMetaData(typeLabel);
            int sc1 = (int)nextTreeEvent.getChild(0).getMetaData(typeLabel);
            int sc2 = (int)nextTreeEvent.getChild(1).getMetaData(typeLabel);

            int sp = sc1 != s ? sc1 : sc2;

            if (sp == s) {
                // Birth

                double birth_prop = trajectory.currentState[s]*param.getBirthRates()[interval][s];
                logWeight += Math.log(birth_prop) - Math.log(0.5*(trajectory.currentState[s]*(trajectory.currentState[s]+1)));

            } else {
                // Cross-birth

                double crossbirth_prop = trajectory.currentState[s]*param.getCrossBirthRates()[interval][s][sp];
                logWeight += Math.log(crossbirth_prop) - Math.log(trajectory.currentState[s]*trajectory.currentState[sp]);

            }

        }

        return logWeight;
    }

    @Override
    public void init(PrintStream out) {
        if (getID() == null)
            out.print("trajectory\t");
        else
            out.print(getID() + "\t");
    }

    @Override
    public void log(long sample, PrintStream out) {
        Trajectory traj = simulateTrajectory();

        if (traj == null) {
            out.print("NA\t");
            return;
        }

        List<double[]> states = traj.getStateList();
        List<Double> eventimes = traj.getEventTimes();

        for (int i=0; i<states.size(); i++ ) {
            out.print(",");

            if (i==0)
                out.print(0.0);
            else
                out.print("," + eventimes.get(i-1));

            for (int s=0; s<nTypes; s++) {
                if (s > 0)
                    out.print(":");
                out.print(states.get(0)[s]);
            }
        }

        out.print("\t");
    }

    @Override
    public void close(PrintStream out) {

    }
}
