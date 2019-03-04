package beast.app.bdmm.beauti;

import bdmm.parameterization.Parameterization;
import beast.app.beauti.BeautiDoc;
import beast.app.draw.BEASTObjectInputEditor;
import beast.core.BEASTInterface;
import beast.core.Input;

public class ParameterizationInputEditor extends BEASTObjectInputEditor {

    public ParameterizationInputEditor(BeautiDoc doc) {
        super(doc);
    }

    @Override
    public Class<?> type() {
        return Parameterization.class;
    }

    @Override
    public void init(Input<?> input, BEASTInterface beastObject, int itemNr,
                     ExpandOption isExpandOption, boolean addButtons) {
        super.init(input, beastObject, itemNr, isExpandOption, addButtons);
    }
}
