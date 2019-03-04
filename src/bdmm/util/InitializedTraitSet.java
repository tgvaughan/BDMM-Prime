package bdmm.util;

import beast.core.Input;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;

import java.util.Arrays;
import java.util.stream.Collectors;

public class InitializedTraitSet extends TraitSet {

    public Input<Tree> treeInput = new Input<>("tree",
    "Tree from which to extract taxa.",
            Input.Validate.REQUIRED);

    public InitializedTraitSet() {
        traitsInput.setRule(Input.Validate.OPTIONAL);
    }

    @Override
    public void initAndValidate() {

        if (traitsInput.get() == null) {
            String value = Arrays.stream(treeInput.get().getTaxaNames())
                    .map(n -> n + "=NOT_SET")
                    .collect(Collectors.joining(","));

            traitsInput.setValue(value, this);
        }

        super.initAndValidate();
    }
}
