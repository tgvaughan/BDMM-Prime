package bdmmprime.beauti;

import bdmmprime.distribution.BirthDeathMigrationDistribution;
import bdmmprime.parameterization.SkylineParameter;
import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beast.base.evolution.tree.TraitSet;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.RealParameter;
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.inputeditor.InputEditor;
import beastfx.app.util.FXUtils;
import javafx.beans.property.SimpleDoubleProperty;
import javafx.beans.value.ObservableDoubleValue;
import javafx.beans.value.ObservableValue;
import javafx.beans.value.ObservableValueBase;
import javafx.collections.FXCollections;
import javafx.event.EventHandler;
import javafx.geometry.Insets;
import javafx.scene.control.*;
import javafx.scene.control.cell.PropertyValueFactory;
import javafx.scene.control.cell.TextFieldTableCell;
import javafx.scene.layout.HBox;
import javafx.scene.layout.VBox;
import javafx.util.Callback;

import java.util.List;

public abstract class SkylineInputEditor extends InputEditor.Base {

    SkylineParameter skylineParameter;

    Spinner<Integer> changeCountSpinner;
    CheckBox scalarRatesCheckBox, timesAreAgesCheckBox;

    VBox changeTimesBox;
    HBox changeTimesEntryRow;

    TableView<RealParameter> valuesTable;

    CheckBox estimateValuesCheckBox, estimateTimesCheckBox;

    VBox mainInputBox;

    CheckBox visualizerCheckBox;
    EpochVisualizerPanel epochVisualizer;

    private boolean modelSaveInProcess = false;

    public SkylineInputEditor(BeautiDoc doc) {
        super(doc);
    }

    @Override
    public void init(Input<?> input, BEASTInterface beastObject, int itemNr,
                     ExpandOption isExpandOption, boolean addButtons) {

        m_bAddButtons = addButtons;
        m_input = input;
        m_beastObject = beastObject;
        this.itemNr = itemNr;
        pane = FXUtils.newVBox();

        skylineParameter = (SkylineParameter)input.get();

        addInputLabel();

        HBox boxHoriz;

        // Add elements specific to change times

        mainInputBox = FXUtils.newVBox();

        boxHoriz = FXUtils.newHBox();
        Label changePointLabel = new Label("Number of change times:");
        changeCountSpinner = new Spinner<>(0, 0, Integer.MAX_VALUE);
        boxHoriz.getChildren().add(changePointLabel);
        boxHoriz.getChildren().add(changeCountSpinner);

        mainInputBox.getChildren().add(boxHoriz);

        changeTimesBox = FXUtils.newVBox();
        changeTimesEntryRow = FXUtils.newHBox();
        changeTimesEntryRow.getChildren().add(new Label("Change times:"));
        changeTimesBox.getChildren().add(changeTimesEntryRow);
        HBox changeTimesBoxRow = FXUtils.newHBox();
        timesAreAgesCheckBox = new CheckBox("Times specified as ages");
        changeTimesBoxRow.getChildren().add(timesAreAgesCheckBox);
        estimateTimesCheckBox = new CheckBox("Estimate change times");
        changeTimesBoxRow.getChildren().add(estimateTimesCheckBox);
        changeTimesBox.getChildren().add(changeTimesBoxRow);

        mainInputBox.getChildren().add(changeTimesBox);

        // Add elements specific to values

        boxHoriz = FXUtils.newHBox();
        boxHoriz.getChildren().add(new Label("Values:"));
        valuesTable = new TableView<>();
        VBox valuesTableBoxCol = FXUtils.newVBox();
        valuesTableBoxCol.getChildren().add(valuesTable);
        boxHoriz.getChildren().add(valuesTableBoxCol);

        mainInputBox.getChildren().add(boxHoriz);

        boxHoriz = FXUtils.newHBox();
        scalarRatesCheckBox = new CheckBox("Scalar values");
        boxHoriz.getChildren().add(scalarRatesCheckBox);
        estimateValuesCheckBox = new CheckBox("Estimate values");
        boxHoriz.getChildren().add(estimateValuesCheckBox);

        mainInputBox.getChildren().add(boxHoriz);

        boxHoriz = FXUtils.newHBox();
        visualizerCheckBox = new CheckBox("Display visualization");
        boxHoriz.getChildren().add(visualizerCheckBox);
        mainInputBox.getChildren().add(boxHoriz);

//        epochVisualizer = new EpochVisualizerPanel(getTree(), getTypeTraitSet(), skylineParameter);
//        mainInputBox.getChildren().add(epochVisualizer);

        getChildren().add(mainInputBox);

        ensureParamsConsistent();

        loadFromModel();

        // Add event listeners:
        changeCountSpinner.valueProperty().addListener(e -> saveToModel());
        timesAreAgesCheckBox.selectedProperty().addListener(e -> saveToModel());
        estimateTimesCheckBox.selectedProperty().addListener(e -> saveToModel());

        scalarRatesCheckBox.selectedProperty().addListener(e -> saveToModel());
        estimateValuesCheckBox.selectedProperty().addListener(e -> saveToModel());

        visualizerCheckBox.selectedProperty().addListener(e -> saveToModel());
    }


    /**
     * Ensures that the dimensions of the values RealParameter is consistent
     * with the number of types specified by the TypeSet.
     *
     * This is necessary because the TypeSet can be modified by other input
     * editors, such as the type trait set input editor.
     *
     * Once returning from this method, the number of types reported by the
     * skyline parameter itself should match the number reported by the typeset.
     */
    abstract void ensureParamsConsistent();

    /**
     * Configure inputs for change times.
     *
     * @param nChanges number of change times (number of columns to configure)
     * @param parameter change times parameter
     */
    void addChangeTimesColumns(int nChanges, RealParameter parameter) {
        for (int i=0; i<nChanges; i++) {
            changeTimesEntryRow.getChildren().add(new Label("Epoch " + (i+1) + "->" + (i+2) + ": "));
            TextField textField = new TextField(parameter.getValue(i).toString());

            textField.setPrefWidth(50);
            textField.setPadding(new Insets(0));
            HBox.setMargin(textField, new Insets(0, 10, 0, 0));

            int index = i;
            textField.setOnAction(event -> {
                parameter.setValue(index, Double.valueOf(textField.getText()));
            });

            changeTimesEntryRow.getChildren().add(textField);
        }
    }

    /**
     * Populate GUI elements with values/dimensions from current BEASTObject model.
     * Called immediately after init(), and thus after every refreshPanel().
     */
    void loadFromModel() {

        int nChanges = skylineParameter.getChangeCount();
        int nTypes = skylineParameter.getNTypes();

        changeCountSpinner.getValueFactory().setValue(nChanges);

        // Load change times:

        if (nChanges > 0) {
            addChangeTimesColumns(nChanges, (RealParameter)skylineParameter.changeTimesInput.get());

            changeTimesBox.setVisible(true);

            estimateTimesCheckBox.setSelected(((RealParameter)skylineParameter.changeTimesInput.get()).isEstimatedInput.get());
        } else {
            changeTimesBox.setVisible(false);
        }

        timesAreAgesCheckBox.setSelected(skylineParameter.timesAreAgesInput.get());

        // Load values

        RealParameter valuesParameter = (RealParameter)skylineParameter.skylineValuesInput.get();

        if (valuesParameter.getDimension()==(nChanges+1)) {
            if (nTypes>1) {
                scalarRatesCheckBox.setSelected(true);
                scalarRatesCheckBox.disableProperty().set(false);
//                epochVisualizer.setScalar(true);
            } else {
                scalarRatesCheckBox.setSelected(false);
                scalarRatesCheckBox.disableProperty().set(true);
//                epochVisualizer.setScalar(false);
            }
        } else {
            scalarRatesCheckBox.setSelected(false);
//            epochVisualizer.setScalar(false);
        }

        estimateValuesCheckBox.setSelected(valuesParameter.isEstimatedInput.get());

        visualizerCheckBox.setSelected(skylineParameter.epochVisualizerDisplayed);
//        epochVisualizer.setVisible(skylineParameter.epochVisualizerDisplayed);
    }

    /**
     * Write anything in the GUI elements to the BEASTObject model.  Called
     * immediately when any of the change event handlers for the GUI objects
     * fire, i.e. when any of the values in the GUI are changed.
     *
     * Ends with a call to refreshPanel().
     */
    void saveToModel() {

        if (modelSaveInProcess)
            return;

        modelSaveInProcess = true;

        int nChanges = changeCountSpinner.getValue();

        // Save values

        RealParameter skylineValuesParam = (RealParameter)skylineParameter.skylineValuesInput.get();
        skylineValuesParam.isEstimatedInput.setValue(estimateValuesCheckBox.isSelected(), skylineValuesParam);
        skylineValuesParam.initAndValidate();

        // Save change times

        RealParameter changeTimesParam = (RealParameter)skylineParameter.changeTimesInput.get();
        if (nChanges>0) {
            if (changeTimesParam == null) {
                String changeTimeId = getChangeTimesParameterID();
                changeTimesParam = (RealParameter) doc.pluginmap.get(changeTimeId);
                if (changeTimesParam == null) {
                    changeTimesParam = new RealParameter("0.0");
                    changeTimesParam.setID(changeTimeId);
                }
            }
            changeTimesParam.setDimension(nChanges);
            changeTimesParam.isEstimatedInput.setValue(estimateTimesCheckBox.isSelected(), changeTimesParam);

        } else {
            if (changeTimesParam != null)
                changeTimesParam.isEstimatedInput.setValue(false, changeTimesParam);
        }

        if (nChanges>0) {
            skylineParameter.setInputValue("changeTimes", changeTimesParam);
            changeTimesParam.initAndValidate();
        } else {
            skylineParameter.setInputValue("changeTimes", null);
        }

        skylineParameter.timesAreAgesInput.setValue(timesAreAgesCheckBox.isSelected(), skylineParameter);

//        skylineParameter.epochVisualizerDisplayed = visualizerCheckBox.isSelected();

        skylineParameter.initAndValidate();

        modelSaveInProcess = false;

        refreshPanel();
    }


    abstract String getChangeTimesParameterID();

    private String getPartitionID() {
        return skylineParameter.getID().split("\\.t:")[1];
    }

}
