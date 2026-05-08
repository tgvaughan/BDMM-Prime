/*
 * Copyright (c) 2017-2026 ETH Zürich
 *
 * This file is part of bdmm-prime.
 *
 * bdmm-prime is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at
 * your option) any later version.
 *
 * bdmm-prime is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with bdmm-prime. If not, see <https://www.gnu.org/licenses/>.
 */

package bdmmprime.beauti;

import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beast.base.spec.inference.parameter.SimplexParam;
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.inputeditor.InputEditor;
import beastfx.app.inputeditor.WrappedOptionPane;
import beastfx.app.util.FXUtils;
import javafx.geometry.Pos;
import javafx.scene.control.*;
import javafx.scene.layout.*;
import javafx.scene.paint.Color;

import java.util.stream.Collectors;

public class SimplexParamInputEditor extends InputEditor.Base {

    public SimplexParamInputEditor(BeautiDoc doc) {
        super(doc);
    }

    public SimplexParamInputEditor() { super(); }

    @Override
    public Class<?> type() {
        return SimplexParam.class;
    }

    @Override
    public void init(Input<?> input, BEASTInterface beastObject, int itemNr, ExpandOption isExpandOption, boolean addButtons) {
        SimplexParam simplexParam = (SimplexParam) input.get();
        int demes = simplexParam.size();

        if (demes < 2)
            return;

        m_bAddButtons = addButtons;
        m_input = input;
        m_beastObject = beastObject;
        this.itemNr = itemNr;
        pane = FXUtils.newHBox();

        addInputLabel();

        HBox mainInputBox = FXUtils.newHBox();

        Label valuesLabel = new Label("Must sum to one:");
        mainInputBox.getChildren().add(valuesLabel);

        TextField textField = new TextField(getParamString(simplexParam));
        textField.setOnAction(_ -> {
            String oldString = getParamString(simplexParam);

            if (textField.getText().split(" ").length != demes) {
                WrappedOptionPane.showWrappedMessageDialog("Number of probability " +
                        "values must match number of types (" + demes + ").");
                textField.setText(oldString);
            } else {
                try {
                    simplexParam.valuesInput.setValue(textField.getText(), simplexParam);
                    simplexParam.initAndValidate();
                } catch (IllegalArgumentException _) {
                    WrappedOptionPane.showWrappedMessageDialog("Invalid values. " +
                            "Probabililities must be positive and add up to 1.");
                    textField.setText(oldString);
                    simplexParam.valuesInput.setValue(oldString, simplexParam);
                    simplexParam.initAndValidate();
                }
            }
        });
        mainInputBox.getChildren().add(textField);

        pane.getChildren().add(mainInputBox);
        getChildren().add(pane);

        // Layout
        textField.prefWidthProperty().bind(mainInputBox.widthProperty()
                .subtract(valuesLabel.widthProperty()).subtract(20));
        HBox.setHgrow(mainInputBox, Priority.ALWAYS);
        mainInputBox.setAlignment(Pos.CENTER);
        pane.prefWidthProperty().bind(widthProperty());

        mainInputBox.setBorder(new Border(new BorderStroke(Color.LIGHTGRAY,
                BorderStrokeStyle.SOLID, null, null)));
    }

    private String getParamString(SimplexParam simplexParam) {
        return simplexParam.getElements().stream()
                .map(Object::toString)
                .collect(Collectors.joining(" "));
    }
}
