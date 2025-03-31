/*
 * Copyright (C) 2019-2024 ETH Zurich
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

import bdmmprime.distribution.BirthDeathMigrationDistribution;
import bdmmprime.mapping.TypeMappedTree;
import bdmmprime.parameterization.Parameterization;
import bdmmprime.trajectories.obsevents.*;
import bdmmprime.util.Utils;
import beast.base.inference.CalculationNode;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.inference.parameter.RealParameter;
import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import static org.apache.commons.math.special.Gamma.logGamma;

public class SampledTrajectory extends CalculationNode implements Loggable {

    public Input<Tree> mappedTreeInput = new Input<>("typeMappedTree",
            "Tree with stochastically mapped types.", Input.Validate.REQUIRED);

    public Input<Parameterization> parameterizationInput = new Input<>("parameterization",
            "Multi-type birth-death parameterization.", Input.Validate.REQUIRED);

    public Input<Function> finalSampleOffsetInput = new Input<>("finalSampleOffset",
            "If provided, the difference in time between the final sample and the end of the BD process.",
            new RealParameter("0.0"));

    public Input<Function> startTypePriorProbsInput = new Input<>("startTypePriorProbs",
            "The prior probabilities for the starting type. Only needed for testing tree prob estimates.");

    public Input<Integer> nParticlesInput = new Input<>("nParticles",
            "Number of particles to use in filtering calculation.", 1000);

    public Input<Double> resampThreshInput = new Input<>("resampThresh",
            "Particle distribution resampling occurs when the relative effective particle count " +
                    "drops below this values.", 0.5);

    public Input<Boolean> useTauLeapingInput = new Input<>("useTauLeaping",
            "If true, use tau-leaping to speed up trajectory simulation.", false);

    public Input<Integer> minLeapCountInput = new Input<>("minLeapCount",
            "Minimum number of of tau-leaping steps executed over the tree.", 100);

    public Input<Double> epsilonInput = new Input<>("epsilon",
            "Tolerance parameter for selecting tau leap length.", 0.03);

    public Input<String> typeLabelInput = new Input<>("typeLabel",
            "Type label used for traits in generated metadata.",
            "type");

    public Input<Boolean> resampleOnLogInput = new Input<>("resampleOnLog",
            "If true, trajectory simulations will be performed at the logging stage.",
            true);

    public Input<BirthDeathMigrationDistribution> bdmmDistribInput = new Input<>("bdmmDistrib",
            "If provided, extract the parameterization from here.",
            Input.Validate.XOR, parameterizationInput);

    Tree mappedTree;
    String typeLabel;
    Parameterization param;
    Function finalSampleOffset;

    int nTypes, nParticles;
    double resampThresh;

    boolean resampleOnLog;

    boolean useTauLeaping;
    int minLeapCount;
    double epsilon;

    Particle[] particles, particlesPrime;
    double[] particleWeights;

    @Override
    public void initAndValidate() {

        if (parameterizationInput.get() != null)
            param = parameterizationInput.get();
        else
            param = bdmmDistribInput.get().parameterizationInput.get();

        finalSampleOffset = finalSampleOffsetInput.get();

        mappedTree = mappedTreeInput.get();
        typeLabel = typeLabelInput.get();
        nTypes = param.getNTypes();
        nParticles = nParticlesInput.get();
        resampThresh = resampThreshInput.get();

        resampleOnLog = resampleOnLogInput.get();

        useTauLeaping = useTauLeapingInput.get();
        minLeapCount = minLeapCountInput.get();
        epsilon = epsilonInput.get();

        particles = new Particle[nParticles];
        particlesPrime = new Particle[nParticles];
        particleWeights = new double[nParticles];
    }

    public double logTreeProbEstimate;

    private long calculationTime;

    public Trajectory sampleTrajectory() {
        // Record start time of calculation for performance monitoring
        calculationTime = System.currentTimeMillis();

        logTreeProbEstimate = 0.0;

        List<ObservedEvent> observedEvents = getObservedEventList(mappedTree);

        int rootType = observedEvents.get(0).type;

        // Initialize particles

        double[] initialState = new double[param.getNTypes()];
        initialState[rootType] = 1.0;

        for (int p=0; p<nParticles; p++) {
            particles[p] = new Particle(param, initialState, useTauLeaping, minLeapCount, epsilon);
            particlesPrime[p] = new Particle(param, initialState, useTauLeaping, minLeapCount, epsilon);
            particleWeights[p] = 0.0;
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

            if (maxLogWeight == Double.NEGATIVE_INFINITY) {
                Log.warning.println("Particle ensemble depleted. Consider re-running with a larger number of particles.");
                return null;
            }

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
                    particlesPrime[p].assignTrajAndZeroWeight(particles[sampler.next()]);

                Particle[] tmp = particles;
                particles = particlesPrime;
                particlesPrime = tmp;

            }

            t = observedEvent.time;
        }


        // Store duration of calculation for performance monitoring
        calculationTime = System.currentTimeMillis() - calculationTime;

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
//        System.out.println(mappedTree + ";");


        if (sampleTrajectory() == null)
            return Double.NaN;

        int rootType = getNodeType(mappedTree.getRoot(), typeLabel);

        // Contribution of root type to log likelihood:
        if (nTypes > 1) {
            if (startTypePriorProbsInput.get() == null)
                throw new IllegalArgumentException("Must provide frequency argument to calculate multi-type tree prob.");

            // This is actually the initial weight of the particles accounting for
            // the probability of a trajectory initialized with the chosen origin state
            // distribution having the observed state.  However, because _all_ of the
            // particles have this initial value, it doesn't affect the weight distribution.
            // It affects the tree prob estimate though, which is important for testing.
            logTreeProbEstimate += Math.log(startTypePriorProbsInput.get().getArrayValue(rootType));
        }

        return logTreeProbEstimate - logGamma(mappedTree.getLeafNodeCount() + 1);
    }

    /**
     * Return the duration in milliseconds of the previous trajectory sampling
     * calculation.  This is used for performance monitoring.
     *
     * @return duration of previous trajectory sampling calculation
     */
    public long getPrevCalculationTime() {
        return calculationTime;
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
            double t = param.getNodeTime(node, finalSampleOffset.getArrayValue());
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
                .toList();

        for (Node node : internalNodes) {
            if (node.getChildCount() == 1) {
                // Observed type change

                eventList.add(new TypeChangeEvent(param.getNodeTime(node, finalSampleOffset.getArrayValue()),
                        getNodeType(node, typeLabel),
                        getNodeType(node.getChild(0), typeLabel),1));
            } else {
                // Coalescence

                eventList.add(new CoalescenceEvent(param.getNodeTime(node, finalSampleOffset.getArrayValue()),
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
        if (param.getNTypes()>1) {
            Object metaData = node.getMetaData(typeLabel);
            if (metaData instanceof Integer)
                return (Integer) metaData;
            else if (metaData instanceof String)
                return Integer.parseInt((String) metaData);
            else
                throw new RuntimeException("Error parsing node type in type-mapped tree.");
        } else
            return 0;
    }

    @Override
    public void init(PrintStream out) {
        Trajectory.init(out);
    }

    Trajectory traj = null;
    long prevSimulationSample = -1;

    @Override
    public void log(long sample, PrintStream out) {
        if (sample==0) {
            Trajectory.logEmpty(out);
            return;
        }

        if (mappedTree instanceof TypeMappedTree)
            ((TypeMappedTree) mappedTree).remapForLog(sample);

        if (resampleOnLog && prevSimulationSample != sample) {
            traj = sampleTrajectory();
            prevSimulationSample = sample;
        }

        if (traj==null)
            Trajectory.logEmpty(out);
        else
            Trajectory.log(sample, traj.getStateList(), traj.events,
                    param.getTotalProcessLength(), out);
    }

    @Override
    public void close(PrintStream out) {

    }
}
