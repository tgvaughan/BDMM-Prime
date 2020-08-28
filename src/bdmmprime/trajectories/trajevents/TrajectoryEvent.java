package bdmmprime.trajectories.trajevents;

import bdmmprime.trajectories.simulation.NodeFactory;
import beast.evolution.tree.Node;

import java.util.List;

public abstract class TrajectoryEvent {
    public double time;
    public int multiplicity;

    public abstract void updateState(double[] state);

    public abstract void reverseUpdateState(double[] state);

    public boolean isSamplingEvent() {
        return false;
    }

    public void simulateSingleTreeEvent(double[] state, List<List<Node>> activeLineages, NodeFactory factory,
                                        Boolean untypedTree) {
        throw new UnsupportedOperationException("Tree event simulation unsupported for this event type.");
    }

    public void simulateTreeEvent(double[] state, List<List<Node>> activeLineages, NodeFactory factory,
                                  Boolean untypedTree) {
        for (int i=0; i<multiplicity; i++)
            simulateSingleTreeEvent(state, activeLineages, factory, untypedTree);
    }
}
