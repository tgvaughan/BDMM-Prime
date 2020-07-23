package bdmmprime.trajectories;

import bdmmprime.parameterization.Parameterization;
import bdmmprime.trajectories.obsevents.*;
import bdmmprime.trajectories.trajevents.*;
import beast.core.CalculationNode;
import beast.core.Function;
import beast.core.Input;
import beast.core.Loggable;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import static org.apache.commons.math.special.Gamma.logGamma;

public class SampledTrajectory extends CalculationNode implements Loggable {

    public Input<Tree> mappedTreeInput = new Input<>("typeMappedTree",
            "Tree with stochastically mapped types.", Input.Validate.REQUIRED);

    public Input<Parameterization> parameterizationInput = new Input<>("parameterization",
            "Multi-type birth-death parameterization.", Input.Validate.REQUIRED);

    public Input<Integer> nParticlesInput = new Input<>("nParticles",
            "Number of particles to use in filtering calculation.", 1000);

    public Input<Double> resampThreshInput = new Input<>("resampThresh",
            "Particle distribution resampling occurs when the relative effective particle count " +
                    "drops below this values.", 0.5);

    public Input<String> typeLabelInput = new Input<>("typeLabel",
            "Type label used for traits in generated metadata.",
            "type");

    public Input<Boolean> resampleOnLogInput = new Input<>("resampleOnLog",
            "If true, trajectory simulations will be performed at the logging stage.",
            true);

    Tree mappedTree;
    String typeLabel;
    Parameterization param;
    int nTypes, nParticles;
    double resampThresh;

    boolean resampleOnLog;

    double[] a_birth, a_death;
    double[][] a_migration, a_crossbirth;

    @Override
    public void initAndValidate() {

        mappedTree = mappedTreeInput.get();
        typeLabel = typeLabelInput.get();
        param = parameterizationInput.get();
        nTypes = param.getNTypes();
        nParticles = nParticlesInput.get();
        resampThresh = resampThreshInput.get();

        resampleOnLog = resampleOnLogInput.get();

        a_birth = new double[nTypes];
        a_death = new double[nTypes];
        a_migration = new double[nTypes][nTypes];
        a_crossbirth = new double[nTypes][nTypes];
    }

    public double logTreeProbEstimate;

    public Trajectory sampleTrajectory() {
        logTreeProbEstimate = 0.0;

        List<ObservedEvent> observedEvents = getObservedEventList(mappedTree);

        int rootType = observedEvents.get(0).type;

        // Initialize particles

        double[] initialState = new double[param.getNTypes()];
        initialState[rootType] = 1.0;

        Trajectory[] particleTrajectories = new Trajectory[nParticles];
        Trajectory[] particleTrajectoriesPrime = new Trajectory[nParticles];

        double[] logParticleWeights = new double[nParticles];
        double[] particleWeights = new double[nParticles];

        for (int p=0; p<nParticles; p++) {
            particleTrajectories[p] = new Trajectory(initialState);
            logParticleWeights[p] = 0.0;
        }

        // Iterate over tree events:

        double t = 0.0;

        for (ObservedEvent observedEvent : observedEvents) {

            // Propagate particles to next event

            int interval = param.getIntervalIndex(t);

            double maxLogWeight = Double.NEGATIVE_INFINITY;
            for (int p=0; p<nParticles; p++) {
                if (!particleTrajectories[p].currentStateValid(observedEvent.lineages))
                    throw new IllegalStateException("Particle state incompatible with next observation.");

                if (logParticleWeights[p] > Double.NEGATIVE_INFINITY)
                    logParticleWeights[p] += propagateParticle(particleTrajectories[p], t, interval, observedEvent);

                if (logParticleWeights[p] > maxLogWeight)
                    maxLogWeight = logParticleWeights[p];
            }

            if (maxLogWeight == Double.NEGATIVE_INFINITY)
                throw new IllegalStateException("Particle ensemble depleted.");

            // Compute sum of scaled weights:
            double sumOfScaledWeights = 0.0, sumOfSquaredScaledWeights = 0.0;
            for (int p = 0; p < nParticles; p++) {
                particleWeights[p] = Math.exp(logParticleWeights[p] - maxLogWeight);
                sumOfScaledWeights += particleWeights[p];
                sumOfSquaredScaledWeights += particleWeights[p]*particleWeights[p];
            }

            double Neff = sumOfScaledWeights*sumOfScaledWeights/sumOfSquaredScaledWeights;

            if (observedEvent.isFinalEvent() || Neff < resampThresh*nParticles) {

                logTreeProbEstimate += Math.log(sumOfScaledWeights/nParticles) + maxLogWeight;

                // Normalize weights:
                for (int p = 0; p < nParticles; p++)
                    particleWeights[p] = particleWeights[p] / sumOfScaledWeights;

                // Resample particle ensemble

                ReplacementSampler sampler = new ReplacementSampler(particleWeights);
                for (int p = 0; p < nParticles; p++) {
                    particleTrajectoriesPrime[p] = particleTrajectories[sampler.next()].copy();
                    logParticleWeights[p] = 0.0;
                }

                Trajectory[] tmp = particleTrajectories;
                particleTrajectories = particleTrajectoriesPrime;
                particleTrajectoriesPrime = tmp;

            }

            t = observedEvent.time;
        }

        // WLOG choose 0th particle as the particle to return:

        return particleTrajectories[0];
    }

    /**
     * Use the particle filter to estimate the marginal probability of the coloured tree.
     * Used for testing/validation.
     *
     * @return log tree probability
     */
    public double getLogTreeProbEstimate() {
        sampleTrajectory();
        return logTreeProbEstimate - logGamma(mappedTree.getLeafNodeCount() + 1);
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

            // Sample time

            double tprime;
            if (a_tot > 0.0)
                tprime = t + Randomizer.nextExponential(a_tot);
            else
                tprime = Double.POSITIVE_INFINITY;

            double tnew = Math.min(tprime, Math.min(param.getIntervalEndTimes()[interval], observedEvent.time));

            // Update weight and time

            logWeight += -a_illegal_tot*(tnew - t);
            t = tnew;

            // Test for end of interval or simulation

            if (param.getIntervalEndTimes()[interval] < observedEvent.time) {
                if (tprime > param.getIntervalEndTimes()[interval]) {

                    // Include probability of seeing no rho-samples:
                    for (int s = 0; s<nTypes; s++) {
                        if(param.getRhoValues()[interval][s] > 0.0)
                            logWeight += trajectory.currentState[s]*Math.log(1.0 - param.getRhoValues()[interval][s]);
                    }

                    interval += 1;
                    continue;
                }
            } else {
                if (tprime > observedEvent.time) {
                    break;
                }
            }

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
            if (!trajectory.currentStateValid(observedEvent.lineages))
                throw new IllegalStateException("Unobserved event produced illegal state.");
        }

        // Compute tree event contribution
        if (!observedEvent.isFinalEvent())
            logWeight += observedEvent.applyToTrajectory(param, interval, trajectory);

        if (!trajectory.currentStateValid())
            throw new IllegalStateException("Observed event produced illegal state.");

        return logWeight;
    }

    /**
     * Condense tree down into a list of observed events for use by the particle filter.
     *
     * @param tree Input tree.
     * @return List of observed event objects.
     *
     */
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

            prevEvent = event;
        }

        // Add event marking end of observation period:
        ObservationEndEvent observationEndEvent = new ObservationEndEvent(param.getTotalProcessLength());
        observationEndEvent.lineages = new int[nTypes];
        eventList.add(observationEndEvent);

        return eventList;
    }

    private int getNodeType(Node node, String typeLabel) {
        if (param.getNTypes()>1)
            return Integer.parseInt((String)node.getMetaData(typeLabel));
        else
            return 0;
    }

    @Override
    public void init(PrintStream out) {
        if (getID() == null)
            out.print("trajectory\t");
        else
            out.print(getID() + "\t");
    }

    Trajectory traj = null;
    long prevSimulationSample = -1;

    @Override
    public void log(long sample, PrintStream out) {
        if (resampleOnLog && prevSimulationSample != sample) {
//            System.out.println("Sampling traj with origin=" + param.originInput.get().getValue() +
//                    " and a tree with " + mappedTree.getLeafNodeCount() + " leaves");
            traj = sampleTrajectory();
            prevSimulationSample = sample;
        }

        if (traj == null)
            out.print("NA");
        else
            out.print(traj);

        out.print("\t");
    }

    @Override
    public void close(PrintStream out) {

    }
}
