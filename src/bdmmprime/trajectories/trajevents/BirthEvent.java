package bdmmprime.trajectories.trajevents;

import bdmmprime.trajectories.simulation.NodeFactory;
import beast.evolution.tree.Node;
import beast.util.Randomizer;

import java.util.List;

public class BirthEvent extends TrajectoryEvent {

    int srcType, destType1, destType2;

    public BirthEvent(double time, int srcType, int destType1, int destType2, int multiplicity) {
        this.time = time;
        this.srcType = srcType;
        this.destType1 = destType1;
        this.destType2 = destType2;
        this.multiplicity = multiplicity;
    }

    public BirthEvent(double time, int srcType, int destType1, int destType2) {
        this(time, srcType, destType1, destType2, 1);
    }

    @Override
    public void updateState(double[] state) {
        state[srcType] -= multiplicity;
        state[destType1] += multiplicity;
        state[destType2] += multiplicity;
    }

    @Override
    public void reverseUpdateState(double[] state) {
        state[srcType] += multiplicity;
        state[destType1] -= multiplicity;
        state[destType2] -= multiplicity;
    }

    @Override
    public void simulateSingleTreeEvent(double[] state, List<List<Node>> activeLineages, NodeFactory nodeFactory,
                                        Boolean untypedTree) {

        if (activeLineages.get(destType1).isEmpty() && activeLineages.get(destType2).isEmpty())
            return;

        boolean included1, included2;

        included1 = Randomizer.nextDouble() < activeLineages.get(destType1).size()/state[destType1];

        double p_included2;
        if (destType2 == destType1) {
            if (included1)
                p_included2 = (activeLineages.get(destType2).size()-1)/(state[destType2]-1);
            else
                p_included2 = activeLineages.get(destType2).size()/(state[destType2]-1);
        } else {
            p_included2 = activeLineages.get(destType2).size()/state[destType2];
        }

        included2 = Randomizer.nextDouble() < p_included2;

        if (!included1 && !included2)
            return;

        if (included1 && included2) {
            // Coalescence
            Node child1 = activeLineages.get(destType1).remove(Randomizer.nextInt(activeLineages.get(destType1).size()));
            Node child2 = activeLineages.get(destType2).remove(Randomizer.nextInt(activeLineages.get(destType2).size()));

            Node parent = nodeFactory.newIntNode(untypedTree ? -1 : srcType, time);
            parent.addChild(child1);
            parent.addChild(child2);

            activeLineages.get(srcType).add(parent);

        } else {
            Node child = null;

            if (included1 && destType1 != srcType)
                child = activeLineages.get(destType1).remove(Randomizer.nextInt(activeLineages.get(destType1).size()));
            else if (included2 && destType2 != srcType)
                child = activeLineages.get(destType2).remove(Randomizer.nextInt(activeLineages.get(destType2).size()));

            if (child != null) {
                Node parent = nodeFactory.newIntNode(untypedTree ? -1 : srcType, time);
                parent.addChild(child);
                activeLineages.get(srcType).add(parent);
            }
        }
    }

    @Override
    public String toString() {
        return "BirthEvent{" +
                "srcType=" + srcType +
                ", destType1=" + destType1 +
                ", destType2=" + destType2 +
                ", multiplicity=" + multiplicity +
                '}';
    }

    @Override
    public String getEventCode() {
        return "B" + ":" +  srcType + ":" + destType1 + ":" + destType2 + ":" + multiplicity;
    }
}
