package bdmmprime.trajectories.trajevents;

import bdmmprime.trajectories.simulation.NodeFactory;
import beast.base.evolution.tree.Node;
import beast.base.util.Randomizer;

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
        state[srcType] -= multiplicity;
        state[destType] += multiplicity;
    }

    @Override
    public void reverseUpdateState(double[] state) {
        state[srcType] += multiplicity;
        state[destType] -= multiplicity;
    }

    @Override
    public void simulateSingleTreeEvent(double[] state, List<List<Node>> activeLineages, NodeFactory nodeFactory,
                                        Boolean untypedTree) {
        if (activeLineages.get(destType).isEmpty())
            return;

        double pMig = activeLineages.get(destType).size()/state[destType];

        if (Randomizer.nextDouble() >= pMig)
            return;

        Node child = activeLineages.get(destType).remove(Randomizer.nextInt(activeLineages.get(destType).size()));

        if (untypedTree) {
            activeLineages.get(srcType).add(child);

        } else {
            Node parent = nodeFactory.newIntNode(srcType, time);
            parent.addChild(child);

            activeLineages.get(srcType).add(parent);
        }
    }

    @Override
    public String toString() {
        return "MigrationEvent{" +
                "srcType=" + srcType +
                ", destType=" + destType +
                ", multiplicity=" + multiplicity +
                '}';
    }

    @Override
    public String getEventCode() {
        return "M\t" + srcType + "\t" + destType + "\t" + multiplicity;
    }
}
