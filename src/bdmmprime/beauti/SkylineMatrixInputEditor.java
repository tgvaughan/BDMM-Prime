package bdmmprime.beauti;

import bdmmprime.parameterization.SkylineMatrixParameter;
import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.inputeditor.InputEditor;
import javafx.beans.binding.Bindings;
import javafx.beans.value.ObservableValueBase;
import javafx.scene.control.Label;
import javafx.scene.control.TableColumn;
import javafx.scene.control.cell.TextFieldTableCell;
import javafx.util.converter.DoubleStringConverter;

import java.util.ArrayList;
import java.util.List;

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
                     InputEditor.ExpandOption isExpandOption, boolean addButtons) {

        super.init(input, beastObject, itemNr, isExpandOption, addButtons);

        skylineMatrix = (SkylineMatrixParameter) input.get();
        skylineMatrix.initAndValidate();

        if (skylineMatrix.getNTypes() == 1) {
            mainInputBox.setVisible(false);
            mainInputBox.setManaged(false);
            pane.getChildren().add(new Label("Insufficient types for this parameter."));
        } else {
            updateValuesUI();
            valuesTable.setFixedCellSize(25);
            valuesTable.prefHeightProperty().bind(valuesTable.fixedCellSizeProperty()
                    .multiply(Bindings.size(valuesTable.getItems()).add(2.2)));
        }

    }

    @Override
    String getChangeTimesParameterID() {
        int idx = skylineMatrix.getID().indexOf("SM");
        String prefix = skylineMatrix.getID().substring(0, idx);
        String suffix = skylineMatrix.getID().substring(idx+2);

        return prefix + "ChangeTimes" + suffix;
    }

    @Override
    void ensureValuesConsistency(boolean scalar) {
        int nTypes = skylineParameter.typeSetInput.get().getNTypes();
        int nEpochs = skylineParameter.changeTimesInput.get() == null
                ? 1
                : skylineParameter.changeTimesInput.get().getDimension() + 1;
        RealParameter valuesParam = (RealParameter) skylineParameter.skylineValuesInput.get();

        System.out.println("Number of epochs: " + nEpochs);

        if (scalar)
            valuesParam.setDimension(nEpochs);
        else
            valuesParam.setDimension(nTypes*nTypes*nEpochs);

        if (skylineParameter.changeTimesInput.get() != null)
            ((RealParameter)skylineParameter.changeTimesInput.get()).initAndValidate();
        sanitiseRealParameter(valuesParam);
        skylineParameter.initAndValidate();
    }

    @Override
    void updateValuesUI() {
        valuesTable.getColumns().clear();
        valuesTable.getItems().clear();

        int nChanges = skylineMatrix.getChangeCount();
        int nTypes = skylineMatrix.getNTypes();

        RealParameter valuesParameter = (RealParameter) skylineMatrix.skylineValuesInput.get();
        boolean scalar = valuesParameter.getDimension() / (nChanges+1) == 1;

        TableColumn<ValuesTableEntry, String> typeCol = new TableColumn<>("From Type");
        typeCol.setCellValueFactory(p -> new ObservableValueBase<>() {
            @Override
            public String getValue() {
                int type = ((MatrixValuesEntry) p.getValue()).fromType;
                return type < 0
                        ? "ALL"
                        : skylineMatrix.typeSetInput.get().getTypeName(type);
            }
        });
        valuesTable.getColumns().add(typeCol);

        List<Integer> toTypes = new ArrayList<>();
        if (scalar)
            toTypes.add(-1);
        else {
            for (int t=0; t<nTypes; t++) {
                toTypes.add(t);
            }
        }

        for (int i=0; i<nChanges+1; i++) {
            TableColumn<ValuesTableEntry, Object> epochCol = new TableColumn<>("Epoch " + (i + 1));
            valuesTable.getColumns().add(epochCol);

            int epochIdx = i;

            for (int toType : toTypes) {
                TableColumn<ValuesTableEntry, Double> col = toType == -1
                        ? new TableColumn<>("to ALL")
                        : new TableColumn<>("to " + skylineMatrix.typeSetInput.get().getTypeName(toType));

                col.setCellValueFactory(p -> new ObservableValueBase<>() {
                    @Override
                    public Double getValue() {
                        int from = ((MatrixValuesEntry) p.getValue()).fromType;
                        return from < 0
                                ? valuesParameter.getValue(epochIdx)
                                : valuesParameter.getValue(epochIdx * nTypes * nTypes + nTypes * from + toType);
                    }
                });
                col.setCellFactory(TextFieldTableCell.forTableColumn(new DoubleStringConverter()));
                col.setOnEditCommit(e -> {
                    int from = ((MatrixValuesEntry) e.getTableView()
                            .getItems().get(e.getTablePosition().getRow())).fromType;
                    if (from < 0) {
                        valuesParameter.setValue(epochIdx, e.getNewValue());
                    } else {
                        valuesParameter.setValue(epochIdx * nTypes * nTypes + nTypes * from + toType, e.getNewValue());
                    }
                    sanitiseRealParameter(valuesParameter);
                    System.out.println(skylineMatrix);
                });
                epochCol.getColumns().add(col);
            }
        }

        if (!scalar) {
            for (int from=0; from<nTypes; from++) {
                valuesTable.getItems().add(new MatrixValuesEntry(from));
            }
        } else {
            valuesTable.getItems().add(new MatrixValuesEntry(-1));
        }
    }

    public static class MatrixValuesEntry extends ValuesTableEntry {
        public int fromType;

        public MatrixValuesEntry(int fromType) {
            this.fromType = fromType;
        }
    }
}
