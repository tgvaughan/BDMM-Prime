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

    TableView<RealParameter> changeTimesTable;
    VBox changeTimesBox;

    TableView valuesTable;
    List<TableColumn<RealParameter,String>> valuesColumns;

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
        HBox changeTimesBoxRow = FXUtils.newHBox();
        changeTimesBoxRow.getChildren().add(new Label("Change times:"));
        changeTimesTable = new TableView<>();
        VBox changeTimesTableBoxCol = FXUtils.newVBox();
        changeTimesTableBoxCol.getChildren().add(changeTimesTable);
        changeTimesBoxRow.getChildren().add(changeTimesTableBoxCol);
        changeTimesBox.getChildren().add(changeTimesBoxRow);
        changeTimesBoxRow = FXUtils.newHBox();
        timesAreAgesCheckBox = new CheckBox("Times specified as ages");
        changeTimesBoxRow.getChildren().add(timesAreAgesCheckBox);
        estimateTimesCheckBox = new CheckBox("Estimate change times");
        changeTimesBoxRow.getChildren().add(estimateTimesCheckBox);
        changeTimesBox.getChildren().add(changeTimesBoxRow);

        mainInputBox.getChildren().add(changeTimesBox);

        // Add elements specific to values

        boxHoriz = FXUtils.newHBox();
        boxHoriz.getChildren().add(new Label("Values:"));
        valuesTable = new TableView<RealParameter>();
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
        changeCountSpinner.editorProperty().addListener(e -> saveToModel());
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
     * Configure table columns for change times.
     * @param nChanges number of change times (number of columns to configure)
     */
    void addChangeTimesColumns(int nChanges) {
        changeTimesTable.getColumns().clear();

        for (int i=0; i<nChanges; i++) {
            int index = i;

            TableColumn<RealParameter, String> thisCol = new TableColumn<>(
                    "Epoch " + (i+1) + "->" + (i+2));
            thisCol.setSortable(false);

            thisCol.setCellValueFactory(p -> new ObservableValueBase<>() {
                @Override
                public String getValue() {
                    return p.getValue().getValue(index).toString();
                }
            });

            thisCol.setCellFactory(TextFieldTableCell.forTableColumn());

            thisCol.setOnEditCommit(event -> event.getRowValue().setValue(
                    index, Double.valueOf(event.getNewValue())));

            changeTimesTable.getColumns().add(thisCol);
        }

        valuesTable.getColumns().clear();

        for (int i=0; i<nChanges+1; i++) {
            int index = i;

            TableColumn<TableColumn<RealParameter,String>, String> thisCol = new TableColumn<>(
                    "Epoch " + (i+1));
            thisCol.setSortable(false);
        }
    }

    /**
     * Populate GUI elements with values/dimensions from current BEASTObject model.
     * Called immediately after init(), and thus after every refreshPanel().
     */
    void loadFromModel() {

        int nChanges = skylineParameter.getChangeCount();
        int nTypes = skylineParameter.getNTypes();

        // Load change times:

        if (nChanges > 0) {
            changeCountSpinner.getValueFactory().setValue(nChanges);
            addChangeTimesColumns(nChanges);

            changeTimesBox.setVisible(true);

            estimateTimesCheckBox.setSelected(((RealParameter)skylineParameter.changeTimesInput.get()).isEstimatedInput.get());

            RealParameter changeTimesParameter = (RealParameter)skylineParameter.changeTimesInput.get();
            changeTimesTable.setItems(FXCollections.observableArrayList(changeTimesParameter));

        } else {
            changeCountSpinner.getValueFactory().setValue(0);
            changeTimesTable.setItems(null);
            changeTimesBox.setVisible(false);
        }

        timesAreAgesCheckBox.setSelected(skylineParameter.timesAreAgesInput.get());

        // Load values

        RealParameter valuesParameter = (RealParameter)skylineParameter.skylineValuesInput.get();

        if (valuesParameter.getDimension()==(nChanges+1)) {
            if (nTypes>1) {
                scalarRatesCheckBox.setSelected(true);
                scalarRatesCheckBox.disableProperty().set(false);
                epochVisualizer.setScalar(true);
            } else {
                scalarRatesCheckBox.setSelected(false);
                scalarRatesCheckBox.disableProperty().set(true);
                epochVisualizer.setScalar(false);
            }
        } else {
            scalarRatesCheckBox.setSelected(false);
            epochVisualizer.setScalar(false);
        }

        estimateValuesCheckBox.setSelected(valuesParameter.isEstimatedInput.get());

        visualizerCheckBox.setSelected(skylineParameter.epochVisualizerDisplayed);
        epochVisualizer.setVisible(skylineParameter.epochVisualizerDisplayed);
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

        // Update table model dimensions

        valuesTableModel.setIntervalCount(nChanges+1);
        valuesTableModel.setScalar(scalarRatesCheckBox.isSelected());

        changeTimesTableModel.setColumnCount(nChanges);
        for (int colIdx=0; colIdx<changeTimesTable.getColumnCount(); colIdx++) {
            if (changeTimesTableModel.getValueAt(0, colIdx) == null)
                changeTimesTableModel.setValueAt(
                        colIdx > 0
                                ? changeTimesTableModel.getValueAt(0, colIdx-1)
                                : 0.0,
                        0, colIdx);
        }

        // Save values

        RealParameter skylineValuesParam = (RealParameter)skylineParameter.skylineValuesInput.get();
        skylineValuesParam.setDimension(valuesTableModel.getRowCount()*(nChanges+1));
        skylineValuesParam.valuesInput.setValue(valuesTableModel.getParameterString(), skylineValuesParam);
        skylineValuesParam.isEstimatedInput.setValue(estimateValuesCheckBox.isSelected(), skylineValuesParam);
        skylineParameter.setInputValue("skylineValues", skylineValuesParam);
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

            StringBuilder changeTimesBuilder = new StringBuilder();
            for (int i=0; i<nChanges; i++)
                changeTimesBuilder.append(" ").append(changeTimesTableModel.getValueAt(0, i));

            changeTimesParam.valuesInput.setValue(changeTimesBuilder.toString(), changeTimesParam);
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

        skylineParameter.epochVisualizerDisplayed = visualizerCheckBox.isSelected();

        skylineParameter.initAndValidate();

        modelSaveInProcess = false;

        refreshPanel();
    }

    abstract String getChangeTimesParameterID();

    private String getPartitionID() {
        return skylineParameter.getID().split("\\.t:")[1];
    }

    protected Tree getTree() {
        return (Tree) doc.pluginmap.get("Tree.t:" + getPartitionID());
    }

    protected TraitSet getTypeTraitSet() {
        BirthDeathMigrationDistribution bdmmPrimeDistrib =
                (BirthDeathMigrationDistribution) doc.pluginmap.get("BDMMPrime.t:" + getPartitionID());

        return bdmmPrimeDistrib.typeTraitSetInput.get();
    }


}
