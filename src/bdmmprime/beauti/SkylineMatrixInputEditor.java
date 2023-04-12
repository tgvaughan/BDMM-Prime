package bdmmprime.beauti;

import bdmmprime.parameterization.SkylineMatrixParameter;
import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.util.FXUtils;
import javafx.scene.control.Label;
import javafx.scene.layout.HBox;

public class SkylineMatrixInputEditor extends SkylineInputEditor {

    SkylineMatrixParameter skylineMatrix;

    public SkylineMatrixInputEditor(BeautiDoc doc) {
        super(doc);
    }

    @Override
    public Class<?> type() {
        return SkylineMatrixParameter.class;
    }

    @Override
    public void init(Input<?> input, BEASTInterface beastObject, int itemNr,
                     ExpandOption isExpandOption, boolean addButtons) {

        skylineMatrix = (SkylineMatrixParameter) input.get();

        super.init(input, beastObject, itemNr, isExpandOption, addButtons);

        if (skylineMatrix.getNTypes()==1) {
            HBox box = FXUtils.newHBox();
            box.getChildren().add(new Label("Insufficient types in model."));
            pane.getChildren().add(box);
        }

        // Custom cell renderer to gray out diagonal entries in values table:

    }

    @Override
    void ensureValuesConsistency(boolean scalar) {
        // TODO
    }

    @Override
    void updateValuesUI() {
        // TODO
    }

    @Override
    String getChangeTimesParameterID() {
        int idx = skylineMatrix.getID().indexOf("SM");
        String prefix = skylineMatrix.getID().substring(0, idx);
        String suffix = skylineMatrix.getID().substring(idx+2);

        return prefix + "ChangeTimes" + suffix;
    }
}
