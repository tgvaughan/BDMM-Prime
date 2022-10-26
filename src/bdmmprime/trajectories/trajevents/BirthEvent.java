package bdmmprime.trajectories.trajevents;

import bdmmprime.trajectories.simulation.NodeFactory;
import beast.base.evolution.tree.Node;
import beast.base.util.Randomizer;

import java.util.List;

public class BirthEvent extends TrajectoryEvent {

    int type;

    public BirthEvent(double time, int type, int multiplicity) {
        this.time = time;
        this.type = type;
        this.multiplicity = multiplicity;
    }

    public BirthEvent(double time, int type) {
        this(time, type, 1);
    }

    @Override
    public void updateState(double[] state) {
        state[type] += multiplicity;
    }

    @Override
    public void reverseUpdateState(double[] state) {
        state[type] -= multiplicity;
    }

    @Override
    public void simulateSingleTreeEvent(double[] state, List<List<Node>> activeLineages, NodeFactory nodeFactory,
                                        Boolean untypedTree) {
        double probCoal = activeLineages.get(type).size()*(activeLineages.get(type).size()-1)
                /(state[type]*(state[type]-1));

        if (Randomizer.nextDouble() >= probCoal)
            return;

        Node child1 = activeLineages.get(type).remove(Randomizer.nextInt(activeLineages.get(type).size()));
        Node child2 = activeLineages.get(type).remove(Randomizer.nextInt(activeLineages.get(type).size()));

        Node parent = nodeFactory.newIntNode(untypedTree ? -1 : type, time);
        parent.addChild(child1);
        parent.addChild(child2);

        activeLineages.get(type).add(parent);
    }

    @Override
    public String toString() {
        return "BirthEvent{" +
                "type=" + type +
                ", multiplicity=" + multiplicity +
                '}';
    }

    @Override
    public String getEventCode() {
        return "B" + ":" +  type + "::" + multiplicity;
    }
}
