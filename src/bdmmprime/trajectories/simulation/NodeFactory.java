package bdmmprime.trajectories.simulation;

import bdmmprime.parameterization.TypeSet;
import beast.evolution.tree.Node;

public class NodeFactory {

    double origin;
    int nextLeafNr, nextIntNr;

    String typeLabel;
    TypeSet typeSet;

    public NodeFactory(double origin, int nSamples, String typeLabel, TypeSet typeSet) {
        this.origin = origin;
        this.nextIntNr = nSamples;
        this.nextLeafNr = 0;
        this.typeLabel = typeLabel;
        this.typeSet = typeSet;
    }

    private Node newNode(int type, double time, int nextNodeNr) {
        Node node = new Node(String.valueOf(nextLeafNr));
        node.setNr(nextNodeNr);

        node.setMetaData(typeLabel, type);
        node.setHeight(origin-time);
        node.metaDataString = String.format("%s=\"%s\"", typeLabel, typeSet.getTypeName(type));

        return node;
    }

    public Node newLeafNode(int type, double time) {
        Node node = newNode(type, time, nextLeafNr);
        nextLeafNr += 1;
        return node;
    }

    public Node newIntNode(int type, double time) {
        Node node = newNode(type, time, nextIntNr);
        nextIntNr += 1;
        return node;
    }
}
