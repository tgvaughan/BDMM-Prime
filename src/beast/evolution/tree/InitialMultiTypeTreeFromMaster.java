package beast.evolution.tree;

import beast.core.Input;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.core.parameter.RealParameter;
import beast.util.Randomizer;
import master.BeastTreeFromMaster;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * User: Denise
 * Date: 11.06.14
 * Time: 14:35
 */
public class InitialMultiTypeTreeFromMaster extends MultiTypeTree implements StateNodeInitialiser {

    public Input<BeastTreeFromMaster> masterTreeInput = new Input<BeastTreeFromMaster>(
            "masterTree",
            "The tree from which traits should be inherited", Input.Validate.REQUIRED);

    public Input<Double> muInput = new Input<Double>("mu",
            "Migration rate for tree proposal", Input.Validate.REQUIRED);

    public Input<Double> popSizeInput = new Input<Double>("popSize",
            "Population size for tree proposal", Input.Validate.REQUIRED);

    public Input<Integer> nTypes = new Input<>("nTypes", "Number of types", Input.Validate.REQUIRED);

    public Input<Boolean> random = new Input<>("random", "Built random tree with traits from BeastTreeFromMaster? If false, tree will copy BeastTreeFromMaster (default true)", true);

    StateNode tree;

    @Override
    public void initAndValidate() throws Exception {

        initStateNodes();
        super.initAndValidate();
    }

    @Override
    public void initStateNodes() throws Exception {

        typeLabel = typeLabelInput.get();
//        nTypes = nTypesInput.get();

        BeastTreeFromMaster masterTree = masterTreeInput.get();

        TraitSet typeTrait = new TraitSet();
        TraitSet dateTrait = new TraitSet();

        String types = "";
        String dates = "";

        for (Node beastNode : masterTree.getExternalNodes()){

            dates += beastNode.getID() + "=" + beastNode.getHeight() +",";
            types += beastNode.getID() + "=" + ((int[])beastNode.getMetaData("location"))[0] +",";

        }

        dates = dates.substring(0,dates.length()-1);
        types = types.substring(0,types.length()-1);

        typeTrait.initByName("value", types, "taxa", m_taxonset.get(), "traitname", "type");
        dateTrait.initByName("value", dates, "taxa", m_taxonset.get(), "traitname", "date-backward");

        MigrationModel migModel = new MigrationModel();

        Double[] temp = new Double[nTypes.get()];
        Arrays.fill(temp, muInput.get());
        migModel.setInputValue("rateMatrix", new RealParameter(temp));
        Arrays.fill(temp, popSizeInput.get());
        migModel.setInputValue("popSizes", new RealParameter(temp));
        migModel.initAndValidate();

        if (random.get()) {
            tree = new StructuredCoalescentMultiTypeTree();

            tree.setInputValue("migrationModel", migModel);
        }
        else{
            Node oldRoot = masterTree.getRoot();
            MultiTypeNode newRoot = new MultiTypeNode();

            newRoot.height = oldRoot.height;
            newRoot.nTypeChanges = 0;
            newRoot.changeTimes.addAll(new ArrayList<Double>());
            newRoot.changeTypes.addAll(new ArrayList<Integer>());
            newRoot.nodeType = 0;

            newRoot.labelNr = oldRoot.labelNr;

            newRoot.addChild(copyFromFlatNode(oldRoot.getLeft()));
            newRoot.addChild(copyFromFlatNode(oldRoot.getRight()));

            tree = new MultiTypeTree(newRoot);
        }

        tree.setInputValue("trait",typeTrait);
        tree.setInputValue("trait",dateTrait);
        tree.initAndValidate();

        setInputValue("trait",dateTrait);
        setInputValue("trait",typeTrait);

        assignFromWithoutID(tree);
    }

    @Override
    public void getInitialisedStateNodes(List<StateNode> stateNodeList) {
        stateNodeList.add(this);
    }

    MultiTypeNode copyFromFlatNode(Node node){

        MultiTypeNode mNode  = new MultiTypeNode();

        mNode.height = node.height;
        mNode.parent = node.parent;

        mNode.nTypeChanges = 0;
        mNode.changeTimes.addAll(new ArrayList<Double>());
        mNode.changeTypes.addAll(new ArrayList<Integer>());
        mNode.nodeType = 0;

        mNode.labelNr = node.labelNr;

        if (node.isLeaf()){

            int type = ((int[]) node.getMetaData("location"))[0];

            if (type!=0) {

                mNode.setNodeType(type);
                mNode.addChange(0, (node.getHeight() + (node.getParent().getHeight() -node.getHeight()) * Randomizer.nextDouble()));
            }

        } else {

            mNode.addChild(copyFromFlatNode(node.getLeft()));
            mNode.addChild(copyFromFlatNode(node.getRight()));
        }

         return mNode;
    }

}
