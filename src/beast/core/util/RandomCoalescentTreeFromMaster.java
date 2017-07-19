package beast.core.util;

import beast.core.Description;
import beast.core.Input;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.core.parameter.IntegerParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.*;
import master.BeastTreeFromMaster;

import java.util.ArrayList;
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
    public Input<Boolean> printTypes = new Input<Boolean>("printTypes", "print the tip types of simulated tree (default false)",false);
    public Input<List<IntegerParameter>> changeType = new Input<>("changeType", "Change the type of tip tip with type changeType at index 0 to the changeType at index 1.", new ArrayList<>(), Input.Validate.OPTIONAL);


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
        for (int i=0; i<taxonCount; i++) {
            actualTaxa.setInputValue("sequence", taxa.sequenceInput.get().get(i));
        }
        actualTaxa.initAndValidate();

        taxonset.initByName("alignment", actualTaxa);
        setInputValue("taxa",actualTaxa);

        String types = "";
        String dates = "";
        Boolean changeTypes = (changeType.get()!=null);
        int[] changeTypeArray = new int[2*changeType.get().size()];
        boolean changed;

        if (changeTypes){
            for (int i=0; i<changeType.get().size(); i++) {
                changeTypeArray[i*2] = changeType.get().get(i).getValue(0);
                changeTypeArray[i*2+1] = changeType.get().get(i).getValue(1);
            }
        }

        if (printTypes.get())
            System.out.println("Tip types of simulated tree:");

        for (Node beastNode : masterTree.getExternalNodes()){

            dates += beastNode.getID() + "=" + beastNode.getHeight() +",";
            if (typesKnown.get()) {
                if (changeTypes){
                    changed = false;
                    for (int i=0; i<changeType.get().size() && !changed; i++) {
                        if (changeTypeArray[i*2] == (Integer) beastNode.getMetaData("location")) {
                            beastNode.setMetaData("location", changeTypeArray[i*2+1]);
                            changed=true;
                }}}
                types += beastNode.getID() + "=" + (beastNode.getMetaData("location")) + ",";
            }
            else
                types += beastNode.getID() + "= -1,"; // -1 refers to unknown type in bdmm

            if (printTypes.get())
                System.out.println(beastNode.getID() + "=" + (beastNode.getMetaData("location")) +",");
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



