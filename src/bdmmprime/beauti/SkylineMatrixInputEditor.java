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
import javafx.scene.control.TableCell;
import javafx.scene.control.TableColumn;
import javafx.scene.control.cell.TextFieldTableCell;
import javafx.util.Callback;
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
            System.out.println(valuesTable.getLayoutBounds());
            valuesTable.setFixedCellSize(25);
            valuesTable.prefHeightProperty().bind(valuesTable.fixedCellSizeProperty()
                    .multiply(Bindings.size(valuesTable.getItems()).add(3.0)));
        }

    }

    @Override
    void ensureValuesConsistency() {
        int nTypes = skylineParameter.typeSetInput.get().getNTypes();
        int nEpochs = skylineParameter.changeTimesInput.get() == null
                ? 1
                : skylineParameter.changeTimesInput.get().getDimension() + 1;
        RealParameter valuesParam = (RealParameter) skylineParameter.skylineValuesInput.get();

        System.out.println("Number of epochs: " + nEpochs);

        if (skylineParameter.isScalarInput.get())
            valuesParam.setDimension(nEpochs);
        else
            valuesParam.setDimension(nTypes*(nTypes-1)*nEpochs);

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
        if (skylineParameter.isScalarInput.get())
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
                        if (from < 0)
                            return valuesParameter.getValue(epochIdx);

                        if (from == toType)
                            return Double.NaN;

                        return valuesParameter.getValue(
                                epochIdx*nTypes*(nTypes-1)
                                + (nTypes-1)*from
                                + (toType<from ? toType : toType-1));
                    }
                });
//                col.setCellFactory(TextFieldTableCell.forTableColumn(new DoubleStringConverter()));
                col.setCellFactory(new Callback<>() {
                    @Override
                    public TableCell<ValuesTableEntry, Double> call(TableColumn<ValuesTableEntry, Double> param) {
                        return new TextFieldTableCell<>(new DoubleStringConverter()) {
                            @Override
                            public void updateItem(Double item, boolean empty) {
                                super.updateItem(item, empty);

                                if (item == null || empty) {
                                    setText("");
                                    setStyle("");
                                    setEditable(false);
                                } else if (Double.isNaN(item)) {
                                    setText("");
                                    setStyle("-fx-background-color: black");
                                    setEditable(false);
                                } else {
                                    setText(Double.toString(item));
                                    setEditable(true);
                                }
                            }
                        };
                    }
                });
                col.setOnEditCommit(e -> {
                    int from = ((MatrixValuesEntry) e.getTableView()
                            .getItems().get(e.getTablePosition().getRow())).fromType;
                    if (from < 0) {
                        valuesParameter.setValue(epochIdx, e.getNewValue());
                    } else {
                        valuesParameter.setValue(epochIdx*nTypes*(nTypes-1)
                                + (nTypes-1)*from
                                + (toType<from ? toType : toType-1), e.getNewValue());
                    }
                    ensureValuesConsistency();
                    sanitiseRealParameter(valuesParameter);
                    System.out.println(skylineMatrix);
                });
                epochCol.getColumns().add(col);
            }
        }

        if (!skylineParameter.isScalarInput.get()) {
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
