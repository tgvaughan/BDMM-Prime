package bdmmprime.trajectories.trajevents;

import bdmmprime.trajectories.simulation.NodeFactory;
import beast.evolution.tree.Node;
import beast.util.Randomizer;

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
        state[type] += 1;
    }

    @Override
    public void reverseUpdateState(double[] state) {
        state[type] -= 1;
    }

    @Override
    public void simulateSingleTreeEvent(double[] state, List<List<Node>> activeLineages, NodeFactory nodeFactory) {
        if (activeLineages.get(type).isEmpty())
            return;

        double probCoal = activeLineages.get(type).size()*(activeLineages.get(type).size()-1)
                /(state[type]*(state[type]-1));

        if (Randomizer.nextDouble() >= probCoal)
            return;

        Node child1 = activeLineages.get(type).remove(Randomizer.nextInt(activeLineages.get(type).size()));
        Node child2 = activeLineages.get(type).remove(Randomizer.nextInt(activeLineages.get(type).size()));

        Node parent = nodeFactory.newIntNode(type, time);
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
}
