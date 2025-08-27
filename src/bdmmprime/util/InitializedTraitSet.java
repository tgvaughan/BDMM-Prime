package bdmmprime.util;

import bdmmprime.beauti.EpochVisualizerPane;
import bdmmprime.parameterization.SkylineParameter;
import beast.base.core.Input;
import beast.base.evolution.tree.TraitSet;
import beast.base.evolution.tree.Tree;

import java.util.stream.Collectors;

/**
 * TraitSet in which value input is optional and the values
 * are initialized to a stand-in value.  Used by the BDMM-Prime BEAUti template,
 * where the trait set must be specified before any value (or indeed
 * the taxa themselves) can be known.
 *
 * Use of this trait set class in the BEAUti template is also important as
 * it causes BDMM-Prime's own traitset input editor to be used rather than
 * some other input editor (such as the one that MTT provdes).
 */
public class InitializedTraitSet extends TraitSet {

    public InitializedTraitSet() {
        traitsInput.setRule(Input.Validate.OPTIONAL);
    }

    @Override
    public void initAndValidate() {

        if (traitsInput.get() == null) {
            String value = taxaInput.get().getTaxaNames().stream()
                    .map(n -> n + "=NOT_SET")
                    .collect(Collectors.joining(","));

            traitsInput.setValue(value, this);
        }

        super.initAndValidate();
    }

    /**
     * Create new EpochVisualizerPane object.  Allows packages derived
     * from BDMM-Prime to specify different visualizers.
     *
     * @param tree Tree whose tips to visualize
     * @param param Skyline parameter whose epochs to visualize
     * @return epoch vizualizer pane to be included in the input editor.
     */
    public EpochVisualizerPane getNewVisualizer(Tree tree, SkylineParameter param) {
        return new EpochVisualizerPane(tree, this, param);
    }
}
