package beast.core.util;

import beast.core.Description;
import beast.core.Input;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.*;
import master.BeastTreeFromMaster;

import java.util.List;

/**
 * User: Denise.kuehnert@gmail.com
 * Date: 07.03.17
 */
@Description("Make a random tree with tip dates and states obtained from MASTER simulation")
public class RandomCoalescentTreeFromMaster extends RandomTree implements StateNodeInitialiser {

    public Input<BeastTreeFromMaster> masterTreeInput = new Input<BeastTreeFromMaster>(
            "masterTree",
            "The tree from which traits should be inherited", Input.Validate.REQUIRED);

    public Input<Double> popSize = new Input<Double>("popSize", "popSize for tree initialization",1.);

    public Input<Boolean> typesKnown = new Input<Boolean>("typesKnown", "Set to false if the tree should not known the leave types. Default true.", true);


    StateNode tree;

    @Override
    public void initAndValidate() {

        init();

        super.initAndValidate();

    }

    public void init() {

        BeastTreeFromMaster masterTree = masterTreeInput.get();

        int taxonCount = masterTree.getLeafNodeCount();
        System.out.println("taxonCount: " + taxonCount);
        nodeCount = masterTree.getNodeCount();

        tree  = new RandomTree();

        TraitSet typeTrait = new TraitSet();
        TraitSet dateTrait = new TraitSet();

        Alignment taxa = taxaInput.get();
        TaxonSet taxonset = new TaxonSet();

        Alignment actualTaxa = new Alignment();
        for (int i=0; i<taxonCount; i++)
            actualTaxa.setInputValue("sequence", taxa.sequenceInput.get().get(i));
        actualTaxa.initAndValidate();

        taxonset.initByName("alignment", actualTaxa);
        setInputValue("taxa",actualTaxa);

        String types = "";
        String dates = "";

        for (Node beastNode : masterTree.getExternalNodes()){

            dates += beastNode.getID() + "=" + beastNode.getHeight() +",";
            if (typesKnown.get())
                types += beastNode.getID() + "=" + (beastNode.getMetaData("location")) +",";
            else
                types += beastNode.getID() + "= -1,"; // -1 refers to unknown type in bdmm

        }

        dates = dates.substring(0,dates.length()-1);
        types = types.substring(0,types.length()-1);

        typeTrait.initByName("value", types, "taxa", taxonset, "traitname", "type");
        dateTrait.initByName("value", dates, "taxa", taxonset, "traitname", "date-backward");

        setInputValue("trait",dateTrait);
        setInputValue("trait",typeTrait);

    }

    public void getInitialisedStateNodes(List<StateNode> stateNodes){

        stateNodes.add(this);
    }

}



