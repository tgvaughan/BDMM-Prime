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
import beast.base.inference.parameter.RealParameter;
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.inputeditor.InputEditor;
import beastfx.app.util.FXUtils;
import javafx.geometry.Insets;
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

        ToggleGroup toggleGroup = new ToggleGroup();

        HBox boxHoriz = FXUtils.newHBox();
        RadioButton rootButton = new RadioButton("Condition on Root");
        rootButton.setTooltip(new Tooltip(
                "Condition analysis (and place priors) on the age of the tree root."));
        rootButton.setToggleGroup(toggleGroup);
        boxHoriz.getChildren().add(rootButton); // For alignment with other elements
        mainInputBox.getChildren().add(boxHoriz);

        boxHoriz = FXUtils.newHBox();
        RadioButton originButton = new RadioButton("Condition on Origin");
        originButton.setTooltip(new Tooltip(
                "Condition analysis (and place priors) on the age " +
                        "of the birth-death process origin."));
        originButton.setToggleGroup(toggleGroup);
        boxHoriz.getChildren().add(originButton);
        CheckBox originEstimate = new CheckBox("Estimate");
        boxHoriz.getChildren().add(originEstimate);
        mainInputBox.getChildren().add(boxHoriz);

        boxHoriz = FXUtils.newHBox();
        Label originLabel = new Label("Initial age of Origin:");
        boxHoriz.getChildren().add(originLabel);
        TextField originTextField = new TextField();
        boxHoriz.getChildren().add(originTextField);
        mainInputBox.getChildren().add(boxHoriz);

        boxHoriz = FXUtils.newHBox();
        CheckBox survivalConditionedCheckBox = new CheckBox("Condition on Survival");
        survivalConditionedCheckBox.setTooltip(new Tooltip(
                "Condition analysis on having sampled at least one individual."));
        survivalConditionedCheckBox.setSelected(bdmm.conditionOnSurvivalInput.get());
        boxHoriz.getChildren().add(survivalConditionedCheckBox); // For alignment with other elements
        mainInputBox.getChildren().add(boxHoriz);

        if (bdmm.conditionOnRootInput.get()) {

            rootButton.setSelected(true);
            bdmm.conditionOnSurvivalInput.setValue(true, bdmm);

            survivalConditionedCheckBox.setDisable(true);
            originLabel.setDisable(true);
            originTextField.setDisable(true);
            originEstimate.setSelected(false);
            originEstimate.setDisable(true);
            processLength.isEstimatedInput.setValue(true, processLength);

        } else if (processLength.originInput.get() != null
                && processLength.originInput.get() instanceof RealParameter originParam) {

            originButton.setSelected(true);
            originTextField.setText(originParam.getValue().toString());
            originEstimate.setSelected(originParam.isEstimated());
            processLength.isEstimatedInput.setValue(originParam.isEstimated(), processLength);

            originTextField.setOnAction(e -> {
                originParam.valuesInput.setValue(originTextField.getText(), originParam);
                originParam.initAndValidate();
                refreshPanel();
                sync();
            });

            originEstimate.setOnAction(e -> {
                originParam.isEstimatedInput.setValue(originEstimate.isSelected(), originParam);
                processLength.isEstimatedInput.setValue(originEstimate.isSelected(), originParam);
                refreshPanel();
                sync();
            });
        }

        originButton.setOnAction(e -> {
            bdmm.conditionOnRootInput.setValue(false, bdmm);
            refreshPanel();
            sync();
        });

        rootButton.setOnAction(e -> {
            bdmm.conditionOnRootInput.setValue(true, bdmm);
            bdmm.conditionOnSurvivalInput.setValue(true, bdmm);
            refreshPanel();
            sync();
        });

        survivalConditionedCheckBox.setOnAction(e -> {
            bdmm.conditionOnSurvivalInput.setValue(survivalConditionedCheckBox.isSelected(), bdmm);
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
