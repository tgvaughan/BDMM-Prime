package bdmmprime.trajectories.trajevents;

import bdmmprime.trajectories.simulation.NodeFactory;
import beast.base.evolution.tree.Node;

import java.util.List;

public abstract class TrajectoryEvent {
    public double time;
    public int multiplicity;

    public abstract void updateState(double[] state);

    public abstract void reverseUpdateState(double[] state);

    public boolean isSamplingEvent() {
        return false;
    }

    /**
     * Simulate a single tree event, paying no attention to the multiplicity field.
     *
     * @param state Number of individuals immediately after the event.
     * @param activeLineages Lineages active immediately after the event.
     * @param factory Factory object for creating new tree nodes.
     * @param untypedTree If true, don't record types/type changes on resulting tree.
     */
    public void simulateSingleTreeEvent(double[] state, List<List<Node>> activeLineages, NodeFactory factory,
                                        Boolean untypedTree) {
        throw new UnsupportedOperationException("Tree event simulation unsupported for this event type.");
    }

    /**
     * Simulate an event on a sampled tree. (Reverse time.)
     *
     * @param state Number of individuals immediately after the event.
     * @param activeLineages Lineages active immediately after the event.
     * @param factory Factory object for generating new tree nodes.
     * @param untypedTree If true, don't record types/type changes on resulting tree.
     */
    public void simulateTreeEvent(double[] state, List<List<Node>> activeLineages, NodeFactory factory,
                                  Boolean untypedTree) {
        for (int i=0; i<multiplicity; i++)
            simulateSingleTreeEvent(state, activeLineages, factory, untypedTree);
    }

    public abstract String getEventCode();
}
