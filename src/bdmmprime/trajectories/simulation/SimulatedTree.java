/*
 * Copyright (C) 2019-2025 ETH Zurich
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

package bdmmprime.trajectories.simulation;

import bdmmprime.parameterization.*;
import bdmmprime.trajectories.Trajectory;
import bdmmprime.trajectories.trajevents.*;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import static bdmmprime.util.Utils.nextBinomial;
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

    public Input<Function> finalSampleOffsetInput = new Input<>("finalSampleOffset",
            "The difference in time between the final sample and the end of the BD process. " +
                    "Will be set by the simulator.", Input.Validate.REQUIRED);

    public Input<RealParameter> startTypePriorProbsInput = new Input<>("startTypePriorProbs",
            "The type probabilities for the first individual.",
            Input.Validate.REQUIRED);

    public Input<Integer> minSamplesInput = new Input<>("minSamples",
            "Minimum number of samples to accept in simulated trajectory.", 1);

    public Input<String> typeLabelInput = new Input<>("typeLabel",
            "Used to label tree nodes with corresponding types.",
            "type");

    public Input<String> trajFileNameInput = new Input<>("trajFileName",
            "Name of file to write simulated trajectory to.");

    public Input<String> treeFileNameInput = new Input<>("treeFileName",
            "Name of file to write simulated tree to.");

    public Input<Boolean> simulateUntypedTreeInput = new Input<>("simulateUntypedTree",
            "If true, an untyped tree will be simulated (i.e. migration events will be removed).",
            false);

    int nTypes;
    String typeLabel;


    boolean simulateUntypedTree;

    public Trajectory traj;
    public TrajectorySimulator trajectorySimulator;

    @Override
    public void initAndValidate() {
        Parameterization param = parameterizationInput.get();

        int minSamples = minSamplesInput.get();
        typeLabel = typeLabelInput.get();
        simulateUntypedTree = simulateUntypedTreeInput.get();
        trajectorySimulator = new TrajectorySimulator(parameterizationInput.get(),
                startTypePriorProbsInput.get());

        traj = null;
        do {
            traj = trajectorySimulator.simulateTrajectory();
        } while (traj.getSampleCount() < Math.max(minSamples,1));

        RealParameter fso = (RealParameter) finalSampleOffsetInput.get();
        fso.setValue(param.processLengthInput.get().getArrayValue() - traj.getFinalSampleTime());

        if (trajFileNameInput.get() != null) {
            try (PrintStream out = new PrintStream(trajFileNameInput.get())) {

                traj.dump(out, false);

            } catch (FileNotFoundException e) {
                System.err.println("Error writing trajectory to file.");
            }
        }

        Tree tree = simulateTree();
        assignFromWithoutID(tree);

        if (treeFileNameInput.get() != null) {
            try (PrintStream out = new PrintStream(treeFileNameInput.get())) {

                out.println(tree);

            } catch (FileNotFoundException e) {
                System.err.println("Error writing tree to file.");
            }
        }

        super.initAndValidate();
    }

    /**
     * Simulate tree by iterating over simulated trajectory events
     * in reverse.
     *
     * @return simulated tree
     */
    public Tree simulateTree() {

        List<TrajectoryEvent> events = new ArrayList<>(traj.events);
        Collections.reverse(events);

        double[] state = traj.currentState.clone();

        List<List<Node>> activeLineages = new ArrayList<>();
        for (int s=0; s<nTypes; s++)
            activeLineages.add(new ArrayList<>());

        NodeFactory nodeFactory = new NodeFactory(traj.getFinalSampleTime(), traj.getSampleCount(),
                typeLabel, parameterizationInput.get().getTypeSet());

        for (TrajectoryEvent event : events) {
                event.simulateTreeEvent(state, activeLineages, nodeFactory, simulateUntypedTree);
                event.reverseUpdateState(state);
        }

        int nRemainingLineages = 0;
        for (int s=0; s<nTypes; s++)
            nRemainingLineages += activeLineages.get(s).size();

        if (nRemainingLineages != 1)
            throw new IllegalStateException("Number of remaining lineages not equal to 1.");

        Node root = null;
        for (int s=0; s<nTypes; s++) {
            if (!activeLineages.get(s).isEmpty()) {
                root = activeLineages.get(s).get(0);
                break;
            }
        }

        if (root == null)
            throw new IllegalStateException("Tree simulation failed.");

        return new Tree(root);
    }

    @Override
    public void log(long sample, PrintStream out) {
        Tree tree = (Tree) getCurrent();
        out.print("tree STATE_" + sample + " = ");
        final int[] dummy = new int[1];
        final String newick = tree.getRoot().toSortedNewick(new int[1], true);
        out.print(newick);
        out.print(";");
    }

    /**
     * Debug method for testing
     * @param args unused
     */
    public static void main(String[] args) {

        Randomizer.setSeed(1);

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
                        new RealParameter("0.1"), 2),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.5"), 2),
                "rhoSampling", new TimedParameter(
                        new RealParameter("2.0 3.0 4.0"),
                        new RealParameter("0.5 0.5 0.5"), 2));

        SimulatedTree sim = new SimulatedTree();
        sim.initByName("parameterization", param,
                "frequencies", new RealParameter("0.5 0.5"),
                "simulationTime", "5.0");

        System.out.println(sim.getRoot().toSortedNewick(new int[1], true));
    }
}
