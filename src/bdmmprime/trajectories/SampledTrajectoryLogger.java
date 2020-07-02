package bdmmprime.trajectories;

import bdmmprime.parameterization.Parameterization;
import bdmmprime.trajectories.obsevents.CoalescenceEvent;
import bdmmprime.trajectories.obsevents.ObservedEvent;
import bdmmprime.trajectories.obsevents.ObservedSamplingEvent;
import bdmmprime.trajectories.obsevents.TypeChangeEvent;
import bdmmprime.trajectories.trajevents.*;
import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Loggable;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

import java.io.PrintStream;
import java.util.*;
import java.util.stream.Collectors;

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

    double[] a_birth, a_death;
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
        a_migration = new double[nTypes][nTypes];
        a_crossbirth = new double[nTypes][nTypes];
    }

    public Trajectory simulateTrajectory() {
        List<ObservedEvent> observedEvents = getObservedEventList(mappedTree);

        int rootType = observedEvents.get(0).type;

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

        for (ObservedEvent observedEvent : observedEvents) {

            // Propagate particles to next event

            int interval = param.getIntervalIndex(t);

            double maxLogWeight = Double.NEGATIVE_INFINITY;
            for (int p=0; p<nParticles; p++) {
                logParticleWeights[p] = propagateParticle(particleTrajectories[p], t, interval, observedEvent);

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

            t = observedEvent.time;
        }

        // WLOG choose 0th particle as the particle to return:

        return particleTrajectories[0];
    }

    double propagateParticle(Trajectory trajectory, double t, int interval, ObservedEvent observedEvent) {
        double logWeight = 0.0;

        while (true) {
            // Compute rates

            double a_temp, p_obs;
            double a_tot = 0.0;
            double a_illegal_tot = 0.0;

            for (int s=0; s<nTypes; s++) {
                a_temp = trajectory.currentState[s]*param.getBirthRates()[interval][s];
                if (a_temp > 0) {
                    p_obs = observedEvent.lineages[s] * (observedEvent.lineages[s] - 1.0) / (trajectory.currentState[s] * (trajectory.currentState[s] + 1));
                    a_birth[s] = a_temp * (1 - p_obs);
                    a_tot += a_birth[s];
                    a_illegal_tot += a_temp * p_obs;
                } else
                    a_birth[s] = 0.0;

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

                    a_temp = trajectory.currentState[s] * param.getMigRates()[interval][s][sp];
                    if (trajectory.currentState[s]>observedEvent.lineages[s]) {
                        p_obs = observedEvent.lineages[sp] / (trajectory.currentState[sp] + 1);
                        a_migration[s][sp] = a_temp * (1 - p_obs);
                        a_tot += a_migration[s][sp];
                        a_illegal_tot += a_temp * p_obs;
                    } else {
                        a_migration[s][sp] = 0.0;
                        a_illegal_tot += a_temp;
                    }

                    a_temp = trajectory.currentState[s]*param.getCrossBirthRates()[interval][s][sp];
                    if (a_temp > 0.0) {
                        p_obs = observedEvent.lineages[s] * observedEvent.lineages[sp] / (trajectory.currentState[s] * (trajectory.currentState[sp] + 1));
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

            if (param.getIntervalEndTimes()[interval] < observedEvent.time) {
                if (tprime > param.getIntervalEndTimes()[interval]) {
                    // TODO: Add in treatment of rho sampling

                    logWeight += -a_illegal_tot*(tprime - t);
                    t = param.getIntervalEndTimes()[interval];
                    interval += 1;
                    continue;
                }
            } else {
                if (t > observedEvent.time) {
                    logWeight += -a_illegal_tot*(tprime - t);
                    t = observedEvent.time;
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

        // Compute tree event contribution
        logWeight += observedEvent.applyToTrajectory(param, interval, trajectory);

        return logWeight;
    }

    List<ObservedEvent> getObservedEventList(Tree tree) {

        // Extract sample times first:
        List<Node> sampleNodes = new ArrayList<>(tree.getExternalNodes());
        sampleNodes.sort(Comparator.comparingDouble(Node::getHeight));
        Collections.reverse(sampleNodes);

        List<ObservedEvent> eventList = new ArrayList<>();
        ObservedSamplingEvent[] thisSamplingEvent = new ObservedSamplingEvent[param.getNTypes()];
        for (Node node : sampleNodes) {
            double t = param.getNodeTime(node);
            int type = getNodeType(node, typeLabel);

            if (thisSamplingEvent[type] == null || Math.abs(t-thisSamplingEvent[type].time) > 1e-10) {
                thisSamplingEvent[type] = new ObservedSamplingEvent(t, type, 0, 0);
                eventList.add(thisSamplingEvent[type]);
            }

            if (node.isDirectAncestor())
                thisSamplingEvent[type].nSampledAncestors += 1;
            else
                thisSamplingEvent[type].nLeaves += 1;
        }

        List<Node> internalNodes = tree.getInternalNodes().stream()
                .filter(n -> !n.isFake())
                .collect(Collectors.toList());

        for (Node node : internalNodes) {
            if (node.getChildCount() == 1) {
                // Observed type change

                eventList.add(new TypeChangeEvent(param.getNodeTime(node),
                        getNodeType(node, typeLabel),
                        getNodeType(node.getChild(0), typeLabel),1));
            } else {
                // Coalescence

                eventList.add(new CoalescenceEvent(param.getNodeTime(node),
                        getNodeType(node, typeLabel),
                        getNodeType(node.getChild(0), typeLabel),
                        getNodeType(node.getChild(1), typeLabel), 1));
            }
        }

        eventList.sort(Comparator.comparingDouble(e -> e.time));

        // Compute lineage counts

        ObservedEvent rootEvent = eventList.get(0);
        int rootType = rootEvent.type;
        rootEvent.lineages = new int[nTypes];
        for (int s=0; s<nTypes; s++)
            rootEvent.lineages[s] = s == rootType ? 1 : 0;

        ObservedEvent prevEvent = null;
        for (ObservedEvent event : eventList) {
            if (prevEvent == null) {
                event.lineages = new int[nTypes];
                for (int s=0; s<nTypes; s++)
                    event.lineages[s] = s == event.type ? 1 : 0;
            } else {
                event.lineages = prevEvent.getNextLineageCounts();
            }
        }

        return eventList;
    }

    private int getNodeType(Node node, String typeLabel) {
        return Integer.parseInt((String)node.getMetaData(typeLabel));
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
        List<Double> eventTimes = traj.getEventTimes();

        for (int i=0; i<states.size(); i++ ) {
            if (i==0)
                out.print(0.0);
            else
                out.print("," + eventTimes.get(i-1));

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
