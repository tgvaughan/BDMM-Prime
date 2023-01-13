package bdmmprime.beauti;

import bdmmprime.parameterization.SkylineVectorParameter;
import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.inputeditor.InputEditor;
import javafx.beans.binding.Bindings;
import javafx.beans.value.ObservableValueBase;
import javafx.scene.control.TableColumn;

public class SkylineVectorInputEditor extends SkylineInputEditor {

    SkylineVectorParameter skylineVector;

    public SkylineVectorInputEditor(BeautiDoc doc) {
        super(doc);
    }

    @Override
    public Class<?> type() {
        return SkylineVectorParameter.class;
    }

    @Override
    public void init(Input<?> input, BEASTInterface beastObject, int itemNr,
                     InputEditor.ExpandOption isExpandOption, boolean addButtons) {

        super.init(input, beastObject, itemNr, isExpandOption, addButtons);

        skylineVector = (SkylineVectorParameter) input.get();
        skylineVector.initAndValidate();
        int nChanges = skylineVector.getChangeCount();
        int nTypes = skylineVector.getNTypes();

        RealParameter valuesParameter = (RealParameter) skylineVector.skylineValuesInput.get();


        TableColumn<ValuesTableEntry, String> typeCol = new TableColumn<>("Type");
        typeCol.setCellValueFactory(p -> new ObservableValueBase<>() {
            @Override
            public String getValue() {
                int type = ((VectorValuesEntry) p.getValue()).type;
                return type < 0
                        ? "ALL"
                        : skylineVector.typeSetInput.get().getTypeName(type);
            }
        });
        valuesTable.getColumns().add(typeCol);
        for (int i=0; i<nChanges+1; i++) {
            TableColumn<ValuesTableEntry, String> col = new TableColumn<>("Epoch " + (i+1));
            int epochIdx = i;
            col.setCellValueFactory(p -> new ObservableValueBase<>() {
                @Override
                public String getValue() {
                    int type = ((VectorValuesEntry)p.getValue()).type;
                    return String.valueOf(type<0
                            ? valuesParameter.getValue(epochIdx)
                            : valuesParameter.getValue(epochIdx*nTypes + type));
                }
            });
            valuesTable.getColumns().add(col);
        }

        if (valuesParameter.getDimension() / (nChanges+1) > 1) {
            for (int type=0; type<nTypes; type++)
                valuesTable.getItems().add(new VectorValuesEntry(type));
        } else {
            valuesTable.getItems().add(new VectorValuesEntry(-1));
        }

        valuesTable.setFixedCellSize(25);
        valuesTable.prefHeightProperty().bind(valuesTable.fixedCellSizeProperty()
                .multiply(Bindings.size(valuesTable.getItems()).add(1.1)));

    }

    @Override
    String getChangeTimesParameterID() {
        int idx = skylineVector.getID().indexOf("SV");
        String prefix = skylineVector.getID().substring(0, idx);
        String suffix = skylineVector.getID().substring(idx+2);

        return prefix + "ChangeTimes" + suffix;
    }

    @Override
    void ensureValuesConsistency() {
        int nTypes = skylineParameter.typeSetInput.get().getNTypes();
        int nEpochs = skylineParameter.changeTimesInput.get() == null
                ? 1
                : skylineParameter.changeTimesInput.get().getDimension() + 1;
        RealParameter valuesParam = (RealParameter) skylineParameter.skylineValuesInput.get();
        int valuesPerEpoch = valuesParam.getDimension() / nEpochs;

        if (valuesParam.getDimension() % nEpochs != 0
                || (valuesPerEpoch != 1 && valuesPerEpoch != nTypes)) {
            valuesParam.setDimension(nTypes*nEpochs);
        }

        if (skylineParameter.changeTimesInput.get() != null)
            ((RealParameter)skylineParameter.changeTimesInput.get()).initAndValidate();
        valuesParam.initAndValidate();
        skylineParameter.initAndValidate();
    }

    public static class VectorValuesEntry extends ValuesTableEntry {
        public int type;

        public VectorValuesEntry(int type) {
            this.type = type;
        }
    }
}
