package bdmmprime.trajectories.trajevents;

import bdmmprime.trajectories.simulation.NodeFactory;
import beast.evolution.tree.Node;

import java.util.List;

public class DeathEvent extends TrajectoryEvent {

    int type;

    public DeathEvent(double time, int type, int multiplicity) {
        this.time = time;
        this.type = type;
        this.multiplicity = multiplicity;
    }

    public DeathEvent(double time, int type) {
        this(time, type, 1);
    }

    @Override
    public void updateState(double[] state) {
        state[type] -= multiplicity;
    }

    public void reverseUpdateState(double[] state) {
        state[type] += multiplicity;
    }

    @Override
    public void simulateSingleTreeEvent(double[] state, List<List<Node>> activeLineages, NodeFactory nodeFactory,
                                        Boolean untypedTree) {
        // Death events don't affect the tree.
    }

    @Override
    public String toString() {
        return "DeathEvent{" +
                "type=" + type +
                ", multiplicity=" + multiplicity +
                '}';
    }

    @Override
    public String getEventCode() {
        return "D" + ":" + type + ":::" + multiplicity;
    }
}
