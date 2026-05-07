package bdmmprime.beauti;

import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beast.base.spec.inference.parameter.SimplexParam;
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.inputeditor.InputEditor;
import beastfx.app.inputeditor.WrappedOptionPane;
import beastfx.app.util.FXUtils;
import javafx.scene.control.CheckBox;
import javafx.scene.control.Label;
import javafx.scene.control.Separator;
import javafx.scene.control.TextField;
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

    SimplexParam simplexParam;
    int demes;

    @Override
    public void init(Input<?> input, BEASTInterface beastObject, int itemNr, ExpandOption isExpandOption, boolean addButtons) {
        simplexParam = (SimplexParam) input.get();
        demes = simplexParam.size();

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

        TextField valuesTextField = new TextField(getParamString());
        valuesTextField.setOnAction(_ -> {
            String oldString = getParamString();
            if (valuesTextField.getText().split(" ").length != demes) {
                WrappedOptionPane.showWrappedMessageDialog(
                        "Number of probability values must match number of " +
                                "types (currently " + demes + ").");
                valuesTextField.setText(oldString);
                return;
            }

            try {
                simplexParam.valuesInput.setValue(valuesTextField.getText(), simplexParam);
                simplexParam.initAndValidate();
            } catch (IllegalArgumentException _) {
                WrappedOptionPane.showWrappedMessageDialog("Invalid values.  " +
                        "Probabilities must be positive and must sum to 1.");
                simplexParam.valuesInput.setValue(oldString, simplexParam);
                simplexParam.initAndValidate();
                valuesTextField.setText(oldString);
            }
        });
        mainInputBox.getChildren().add(valuesTextField);

        mainInputBox.setBorder(new Border(new BorderStroke(Color.LIGHTGRAY,
                BorderStrokeStyle.SOLID, null, null)));

        HBox.setHgrow(mainInputBox, Priority.ALWAYS);

        pane.getChildren().add(mainInputBox);

        getChildren().add(pane);
    }

    private String getParamString() {
        return simplexParam.getElements().stream()
                .map(Object::toString)
                .collect(Collectors.joining(" "));
    }
}
