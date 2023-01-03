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

        // Custom cell renderer to gray out diagonal entries in values table:

    }

    @Override
    void ensureParamsConsistent() {

        skylineMatrix.typeSetInput.get().initAndValidate();

        int nTypes = skylineMatrix.typeSetInput.get().getNTypes();
        int nIntervals = skylineMatrix.getChangeCount() + 1;

        RealParameter valuesParam = (RealParameter)skylineMatrix.skylineValuesInput.get();
        int valuesPerInterval = valuesParam.getDimension() / nIntervals;

        if (valuesPerInterval == 1 || valuesPerInterval == nTypes*(nTypes-1)) {

            if (nTypes==1)
                valuesParam.isEstimatedInput.setValue(false, valuesParam);

            skylineMatrix.initAndValidate();
            return;
        }

        StringBuilder valueBuilder = new StringBuilder();

        for (int interval=0; interval<nIntervals; interval++) {
            for (int type1Idx = 0; type1Idx < nTypes; type1Idx++) {
                for (int type2Idx = 0; type2Idx < nTypes; type2Idx++) {
                    if (type1Idx == type2Idx)
                        continue;

                    valueBuilder.append(" 0.0");
                }
            }
        }

        valuesParam.valuesInput.setValue(valueBuilder.toString(), valuesParam);

        valuesParam.initAndValidate();

        skylineMatrix.initAndValidate();
    }

    @Override
    void loadFromModel() {
        if (skylineMatrix.getNTypes()==1) {
            mainInputBox.setVisible(false);
            HBox box = FXUtils.newHBox();
            box.getChildren().add(new Label("Insufficient types in model."));
            getChildren().add(box);
            return;
        }

        super.loadFromModel();
    }

    @Override
    void saveToModel() {
        if (skylineMatrix.getNTypes()==1) {
            RealParameter skylineValues = (RealParameter)skylineMatrix.skylineValuesInput.get();
            skylineValues.isEstimatedInput.setValue(
                    false, skylineValues);
        }

        super.saveToModel();
    }

    @Override
    String getChangeTimesParameterID() {
        int idx = skylineMatrix.getID().indexOf("SM");
        String prefix = skylineMatrix.getID().substring(0, idx);
        String suffix = skylineMatrix.getID().substring(idx+2);

        return prefix + "ChangeTimes" + suffix;
    }
}
