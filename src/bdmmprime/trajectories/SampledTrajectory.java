package bdmmprime.trajectories;

import bdmmprime.parameterization.Parameterization;
import bdmmprime.trajectories.obsevents.*;
import bdmmprime.util.Utils;
import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Loggable;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

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

    public Input<Boolean> useTauLeapingInput = new Input<>("useTauLeaping",
            "If true, use tau-leaping to speed up trajectory simulation.", false);

    public Input<Integer> stepsPerIntervalInput = new Input<>("stepsPerInterval",
            "Number of tau-leaping steps per inter-observation interval.", 10);

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

    boolean useTauLeaping;
    int stepsPerInterval;

    /**
     * Arrays which hold propensities during BD simulations.
     * Fields because we don't want to reallocate these for each particle step.
     */
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

        useTauLeaping = useTauLeapingInput.get();
        stepsPerInterval = stepsPerIntervalInput.get();

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

        Particle[] particles = new Particle[nParticles];
        Particle[] particlesPrime = new Particle[nParticles];

        double[] particleWeights = new double[nParticles];

        for (int p=0; p<nParticles; p++) {
            particles[p] = new Particle(param, initialState, useTauLeaping, stepsPerInterval);
            particlesPrime[p] = new Particle(param, initialState, useTauLeaping, stepsPerInterval);
        }

        // Iterate over tree events:

        double t = 0.0;

        for (ObservedEvent observedEvent : observedEvents) {

            // Propagate particles to next event

            int interval = param.getIntervalIndex(t);

            double maxLogWeight = Double.NEGATIVE_INFINITY;
            for (int p=0; p<nParticles; p++) {
                particles[p].propagateParticle(t, interval, observedEvent);

                if (particles[p].logWeight > maxLogWeight)
                    maxLogWeight = particles[p].logWeight;
            }

            if (maxLogWeight == Double.NEGATIVE_INFINITY)
                throw new IllegalStateException("Particle ensemble depleted.");

            // Compute sum of scaled weights:
            double sumOfScaledWeights = 0.0, sumOfSquaredScaledWeights = 0.0;
            for (int p = 0; p < nParticles; p++) {
                particleWeights[p] = Math.exp(particles[p].logWeight - maxLogWeight);
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
                for (int p = 0; p < nParticles; p++)
                    particlesPrime[p].assignFrom(particles[sampler.next()]);

                Particle[] tmp = particles;
                particles = particlesPrime;
                particlesPrime = tmp;

            }

            t = observedEvent.time;
        }

        // WLOG choose 0th particle as the particle to return:

        return particles[0].trajectory;
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

            if (thisSamplingEvent[type] == null || !Utils.equalWithPrecision(t,thisSamplingEvent[type].time)) {
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