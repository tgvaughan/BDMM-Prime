package bdmmprime.trajectories.trajevents;

import bdmmprime.trajectories.simulation.NodeFactory;
import beast.evolution.tree.Node;
import beast.util.Randomizer;

import java.util.List;

public class MigrationEvent extends TrajectoryEvent {

    int srcType, destType;

    public MigrationEvent(double time, int srcType, int destType, int multiplicity) {
        this.time = time;
        this.srcType = srcType;
        this.destType = destType;
        this.multiplicity = multiplicity;
    }

    public MigrationEvent(double time, int srcType, int destType) {
        this(time, srcType, destType, 1);
    }

    @Override
    public void updateState(double[] state) {
        state[srcType] -= 1;
        state[destType] += 1;
    }

    @Override
    public void reverseUpdateState(double[] state) {
        state[srcType] += 1;
        state[destType] -= 1;
    }

    @Override
    public void simulateSingleTreeEvent(double[] state, List<List<Node>> activeLineages, NodeFactory nodeFactory) {
        if (activeLineages.get(destType).isEmpty())
            return;

        double pMig = activeLineages.get(destType).size()/state[destType];

        if (Randomizer.nextDouble() >= pMig)
            return;

        Node child = activeLineages.get(destType).remove(Randomizer.nextInt(activeLineages.get(destType).size()));

        Node parent = nodeFactory.newIntNode(srcType, time);
        parent.addChild(child);

        activeLineages.get(srcType).add(parent);
    }

    @Override
    public String toString() {
        return "MigrationEvent{" +
                "srcType=" + srcType +
                ", destType=" + destType +
                ", multiplicity=" + multiplicity +
                '}';
    }
}
