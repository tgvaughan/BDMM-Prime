/*
 * Copyright (C) 2019-2025 ETH Zurich
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package bdmmprime.beauti;

import bdmmprime.parameterization.SkylineVectorParameter;
import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;
import beast.base.math.matrixalgebra.Vector;
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.inputeditor.InputEditor;
import javafx.beans.binding.Bindings;
import javafx.beans.value.ObservableValueBase;
import javafx.scene.control.TableColumn;
import javafx.scene.control.cell.TextFieldTableCell;
import javafx.util.converter.DoubleStringConverter;

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

        updateValuesUI();

        valuesTable.setFixedCellSize(25);
        valuesTable.prefHeightProperty().bind(valuesTable.fixedCellSizeProperty()
                .multiply(Bindings.size(valuesTable.getItems()).add(1.1)));

    }

    @Override
    void ensureValuesConsistency() {
        int nTypes = skylineParameter.typeSetInput.get().getNTypes();
        int nEpochs = skylineParameter.changeTimesInput.get() == null
                ? 1
                : skylineParameter.changeTimesInput.get().getDimension() + 1;
        RealParameter valuesParam = (RealParameter) skylineParameter.skylineValuesInput.get();

//        System.out.println("Number of epochs: " + nEpochs);

        if (skylineParameter.isScalarInput.get())
            valuesParam.setDimension(nEpochs);
        else
            valuesParam.setDimension(nTypes*nEpochs);

        if (skylineParameter.changeTimesInput.get() != null)
            ((RealParameter)skylineParameter.changeTimesInput.get()).initAndValidate();
        sanitiseRealParameter(valuesParam);
        skylineParameter.initAndValidate();
    }

    @Override
    void updateValuesUI() {
        valuesTable.getColumns().clear();
        valuesTable.getItems().clear();

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
            TableColumn<ValuesTableEntry, Double> col = new TableColumn<>("Epoch " + (i+1));
            int epochIdx = i;
            col.setCellValueFactory(p -> new ObservableValueBase<>() {
                @Override
                public Double getValue() {
                    int type = ((VectorValuesEntry)p.getValue()).type;
                    return type<0
                            ? valuesParameter.getValue(epochIdx)
                            : valuesParameter.getValue(epochIdx*nTypes + type);
                }
            });
            col.setCellFactory(TextFieldTableCell.forTableColumn(new DoubleStringConverter()));
            col.setOnEditCommit(e -> {
                int type = ((VectorValuesEntry) e.getTableView()
                        .getItems().get(e.getTablePosition().getRow())).type;
                if (type < 0) {
                    valuesParameter.setValue(epochIdx, e.getNewValue());
                } else {
                    valuesParameter.setValue(epochIdx*nTypes + type, e.getNewValue());
                }
                sanitiseRealParameter(valuesParameter);
            });
            valuesTable.getColumns().add(col);
        }

        if (valuesParameter.getDimension() / (nChanges+1) > 1) {
            for (int type=0; type<nTypes; type++)
                valuesTable.getItems().add(new VectorValuesEntry(type));
        } else {
            valuesTable.getItems().add(new VectorValuesEntry(-1));
        }
    }

    public static class VectorValuesEntry extends ValuesTableEntry {
        public int type;

        public VectorValuesEntry(int type) {
            this.type = type;
        }
    }
}
