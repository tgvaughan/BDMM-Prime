package bdmmprime.trajectories.trajevents;

import bdmmprime.trajectories.simulation.NodeFactory;
import beast.evolution.tree.Node;
import beast.util.Randomizer;

import java.util.List;

public class CrossBirthEvent extends TrajectoryEvent {

    int srcType, destType;

    public CrossBirthEvent(double time, int srcType, int destType, int multiplicity) {
        this.time = time;
        this.srcType = srcType;
        this.destType = destType;
        this.multiplicity = multiplicity;
    }

    public CrossBirthEvent(double time, int srcType, int destType) {
        this(time, srcType, destType, 1);
    }

    @Override
    public void updateState(double[] state) {
        state[destType] += 1;
    }

    @Override
    public void reverseUpdateState(double[] state) {
        state[destType] -= 1;
    }

    @Override
    public void simulateSingleTreeEvent(double[] state, List<List<Node>> activeLineages, NodeFactory nodeFactory) {
        if (activeLineages.get(destType).isEmpty() || activeLineages.get(srcType).isEmpty())
            return;

        double pCoal = activeLineages.get(destType).size()*activeLineages.get(srcType).size()
                / (state[destType]*state[srcType]);

        if (Randomizer.nextDouble() >= pCoal)
            return;

        Node child1 = activeLineages.get(srcType).remove(Randomizer.nextInt(activeLineages.get(srcType).size()));
        Node child2 = activeLineages.get(destType).remove(Randomizer.nextInt(activeLineages.get(destType).size()));

        Node parent = nodeFactory.newIntNode(srcType, time);
        parent.addChild(child1);
        parent.addChild(child2);

        activeLineages.get(srcType).add(parent);
    }

    @Override
    public String toString() {
        return "CrossBirthEvent{" +
                "srcType=" + srcType +
                ", destType=" + destType +
                ", multiplicity=" + multiplicity +
                '}';
    }
}
