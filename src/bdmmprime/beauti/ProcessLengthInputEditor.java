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

import bdmmprime.distribution.BirthDeathMigrationDistribution;
import bdmmprime.parameterization.Parameterization;
import bdmmprime.util.ProcessLength;
import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beast.base.evolution.tree.TraitSet;
import beast.base.inference.parameter.RealParameter;
import beastfx.app.inputeditor.BEASTObjectInputEditor;
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.inputeditor.InputEditor;
import beastfx.app.inputeditor.InputEditorFactory;
import beastfx.app.util.FXUtils;
import javafx.geometry.Orientation;
import javafx.scene.control.*;
import javafx.scene.layout.*;
import javafx.scene.paint.Color;

public class ProcessLengthInputEditor extends InputEditor.Base {

    public ProcessLengthInputEditor(BeautiDoc doc) {
        super(doc);
    }

    @Override
    public Class<?> type() {
        return ProcessLength.class;
    }

    @Override
    public void init(Input<?> input, BEASTInterface beastObject, int itemNr, ExpandOption isExpandOption, boolean addButtons) {

        m_bAddButtons = addButtons;
        m_input = input;
        m_beastObject = beastObject;
        this.itemNr = itemNr;

        pane = FXUtils.newHBox();
        VBox mainInputBox = FXUtils.newVBox();

        ProcessLength processLength = (ProcessLength)input.get();
        BirthDeathMigrationDistribution bdmm = getDistrib(processLength);

        if (m_bAddButtons) {
            addInputLabel("Process Conditioning",
                    "Specify details of the length and conditioning of the birth-death process.");
        }


        if (processLength.originInput.get() != null
                && processLength.originInput.get() instanceof RealParameter originParam) {
            HBox boxHoriz = FXUtils.newHBox();
            Label label = new Label("Process length:");
            boxHoriz.getChildren().add(label);
            TextField originTextField = new TextField(originParam.getValue().toString());
            originTextField.setOnAction(e -> {
                originParam.valuesInput.setValue(originTextField.getText(), originParam);
                originParam.initAndValidate();
                refreshPanel();
                sync();
            });
            boxHoriz.getChildren().add(originTextField);
            mainInputBox.getChildren().add(boxHoriz);
        }


        CheckBox survivalConditionedCheckBox = new CheckBox("Condition on Survival");
        if (bdmm.conditionOnRootInput.get()) {
            bdmm.conditionOnSurvivalInput.setValue(true, bdmm);
            survivalConditionedCheckBox.setDisable(true);
        }
        mainInputBox.getChildren().add(survivalConditionedCheckBox);
        survivalConditionedCheckBox.setSelected(bdmm.conditionOnSurvivalInput.get());
        survivalConditionedCheckBox.setOnAction(e -> {
            bdmm.conditionOnSurvivalInput.setValue(survivalConditionedCheckBox.isSelected(), bdmm);
        });

        CheckBox rootConditionedCheckBox = new CheckBox("Condition on Root");
        rootConditionedCheckBox.setSelected(bdmm.conditionOnRootInput.get());
        mainInputBox.getChildren().add(rootConditionedCheckBox);
        rootConditionedCheckBox.setOnAction(e -> {
            bdmm.conditionOnRootInput.setValue(rootConditionedCheckBox.isSelected(), bdmm);
            refreshPanel();
            sync();
        });

        Separator sep = new Separator(Orientation.HORIZONTAL);
        sep.prefWidthProperty().bind(widthProperty());
        getChildren().add(sep);

        mainInputBox.setBorder(new Border(new BorderStroke(Color.LIGHTGRAY,
                BorderStrokeStyle.SOLID, null, null)));

        // Ensure input editor fills space allocated to it
        HBox.setHgrow(mainInputBox, Priority.ALWAYS);
        pane.prefWidthProperty().bind(widthProperty());

        pane.getChildren().add(mainInputBox);
        getChildren().add(pane);
    }

    protected BirthDeathMigrationDistribution getDistrib(ProcessLength processLength) {

        for (BEASTInterface out : processLength.getOutputs()) {
            if (out instanceof Parameterization param) {
                for (BEASTInterface out2 : param.getOutputs()) {
                    if (out2 instanceof  BirthDeathMigrationDistribution bdmmPrimeDistrib) {
                        return bdmmPrimeDistrib;
                    }
                }
            }
        }

        throw new IllegalStateException("ProcessLength not connected " +
                "to a BirthDeathMigrationModel.");
    }
}
