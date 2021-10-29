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
        state[destType] += multiplicity;
    }

    @Override
    public void reverseUpdateState(double[] state) {
        state[destType] -= multiplicity;
    }

    @Override
    public void simulateSingleTreeEvent(double[] state, List<List<Node>> activeLineages, NodeFactory nodeFactory,
                                        Boolean untypedTree) {
        if (activeLineages.get(destType).isEmpty())
            return;

        double pObsStateChange = activeLineages.get(destType).size()/state[destType];
        double pCoal = pObsStateChange*activeLineages.get(srcType).size()/state[srcType];

        double u = Randomizer.nextDouble();

        if (u < pCoal) {
            // Coalescence

            Node child1 = activeLineages.get(srcType).remove(Randomizer.nextInt(activeLineages.get(srcType).size()));
            Node child2 = activeLineages.get(destType).remove(Randomizer.nextInt(activeLineages.get(destType).size()));

            Node parent = nodeFactory.newIntNode(untypedTree ? -1 : srcType, time);
            parent.addChild(child1);
            parent.addChild(child2);

            activeLineages.get(srcType).add(parent);

        } else if (u < pObsStateChange) {
            // Lineage state change

            Node child = activeLineages.get(destType).remove(Randomizer.nextInt(activeLineages.get(destType).size()));
            if (untypedTree) {
                activeLineages.get(srcType).add(child);
            } else {
                Node parent = nodeFactory.newIntNode(srcType, time);
                parent.addChild(child);
                activeLineages.get(srcType).add(parent);
            }
        }
    }

    @Override
    public String toString() {
        return "CrossBirthEvent{" +
                "srcType=" + srcType +
                ", destType=" + destType +
                ", multiplicity=" + multiplicity +
                '}';
    }

    @Override
    public String getEventCode() {
        return "C" + ":" + srcType + ":" + destType + ":" + multiplicity;
    }
}
