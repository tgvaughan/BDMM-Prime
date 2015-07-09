package beast.evolution.tree;

import beast.core.Input;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.core.parameter.RealParameter;
import master.BeastTreeFromMaster;

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

        tree = new StructuredCoalescentMultiTypeTree();

        tree.setInputValue("migrationModel",migModel);
//        tree.setInputValue("nTypes",migModel.getNTypes());
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

}
