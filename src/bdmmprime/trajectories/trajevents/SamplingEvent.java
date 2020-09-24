package bdmmprime.trajectories.trajevents;

import bdmmprime.trajectories.simulation.NodeFactory;
import beast.evolution.tree.Node;
import beast.util.Randomizer;

import java.util.ArrayList;
import java.util.List;

public class SamplingEvent extends TrajectoryEvent {

    int type;
    int nRemoveSamp, nNoRemoveSamp;

    public SamplingEvent(double time, int type, int nRemoveSamp, int nNoRemoveSamp) {
        this.time = time;
        this.type = type;
        this.nRemoveSamp = nRemoveSamp;
        this.nNoRemoveSamp = nNoRemoveSamp;
        this.multiplicity = nRemoveSamp + nNoRemoveSamp;
    }

    @Override
    public void updateState(double[] state) {
        state[type] -= nRemoveSamp;
    }

    @Override
    public void reverseUpdateState(double[] state) {
        state[type] += nRemoveSamp;
    }

    @Override
    public void simulateTreeEvent(double[] state, List<List<Node>> activeLineages, NodeFactory factory,
                                  Boolean untypedTree) {

        // Add nodes corresponding to sampling WITHOUT removal

        List<Node> unsampledLineages = new ArrayList<>(activeLineages.get(type));
        double unusedLineages = state[type];

        for (int i=0; i<nNoRemoveSamp; i++) {
            Node sampledNode = factory.newLeafNode(type, time);
            double pSampledAncestor = unsampledLineages.size() / unusedLineages;

            if (pSampledAncestor == 1.0 || (pSampledAncestor > 0.0 && Randomizer.nextDouble() < pSampledAncestor)) {
                Node child = unsampledLineages.remove(Randomizer.nextInt(unsampledLineages.size()));
                activeLineages.get(type).remove(child);

                Node fake = factory.newIntNode(untypedTree ? -1 : type, time);
                fake.addChild(child);
                fake.addChild(sampledNode);

                activeLineages.get(type).add(fake);
            } else {
                activeLineages.get(type).add(sampledNode);
            }

            unusedLineages -= 1;
        }

        // Add nodes corresponding to sampling WITH removal

        for (int i=0; i<nRemoveSamp; i++)
            activeLineages.get(type).add(factory.newLeafNode(type, time));
    }

    @Override
    public boolean isSamplingEvent() {
        return true;
    }

    @Override
    public String toString() {
        return "SamplingEvent{" +
                "type=" + type +
                ", nRemoveSamp=" + nRemoveSamp +
                ", nNoRemoveSamp=" + nNoRemoveSamp +
                '}';
    }

    @Override
    public String getEventCode() {
        return "S" + ":" + type + "::" + (nRemoveSamp + nNoRemoveSamp);
    }
}
