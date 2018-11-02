package bdmm.util;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Traitset for tip dates obtained from input tree")
public class TipDatesFromTree extends TraitSet{

    public Input<Tree> treeInput = new Input<>("tree", "Tree from which to " +
            "extract tip dates.", Input.Validate.REQUIRED);

    public TipDatesFromTree() {
        traitsInput.setRule(Input.Validate.OPTIONAL);
        traitNameInput.setRule(Input.Validate.OPTIONAL);
    }

    @Override
    public void initAndValidate() {

        traitNameInput.setValue("date-backward", this);

        StringBuilder valueBuilder = new StringBuilder();

        boolean isFirst = true;
        for (Node leaf : treeInput.get().getExternalNodes()) {
            if (isFirst)
                isFirst = false;
            else
                valueBuilder.append(",");
            valueBuilder.append(leaf.getID() + "=" + leaf.getHeight());
        }

        traitsInput.setValue(valueBuilder.toString(), this);

        super.initAndValidate();
    }
}
