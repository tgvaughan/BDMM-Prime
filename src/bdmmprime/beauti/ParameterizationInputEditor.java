package bdmmprime.beauti;

import bdmmprime.parameterization.Parameterization;
import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beastfx.app.inputeditor.BEASTObjectInputEditor;
import beastfx.app.inputeditor.BeautiDoc;

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

        m_bAddButtons = addButtons;
        m_input = input;
        m_beastObject = beastObject;
        this.itemNr = itemNr;

        super.init(input, beastObject, itemNr, isExpandOption, addButtons);
    }
}
