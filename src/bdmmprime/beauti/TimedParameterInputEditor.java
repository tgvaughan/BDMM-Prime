package bdmmprime.beauti;

import bdmmprime.parameterization.TimedParameter;
import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.inputeditor.InputEditor;
import beastfx.app.util.FXUtils;
import javafx.beans.value.ObservableValueBase;
import javafx.geometry.Insets;
import javafx.scene.control.*;
import javafx.scene.control.cell.TextFieldTableCell;
import javafx.scene.layout.HBox;
import javafx.scene.layout.VBox;

public class TimedParameterInputEditor extends InputEditor.Base {

    TimedParameter timedParameter;

    Spinner<Integer> elementCountSpinner;

    HBox timesEntryRow;
    VBox elementsBox;

    CheckBox scalarValuesCheckBox, timesAreAgesCheckBox;

    TableView<TimedParamValuesEntry> valuesTable;

    CheckBox estimateValuesCheckBox, estimateTimesCheckBox;

    boolean modelSaveInProcess = false;

    public TimedParameterInputEditor(BeautiDoc doc) {
        super(doc);
    }

    @Override
    public Class<?> type() {
        return TimedParameter.class;
    }

    @Override
    public void init(Input<?> input, BEASTInterface beastObject, int itemNr,
                     ExpandOption isExpandOption, boolean addButtons) {

        m_bAddButtons = addButtons;
        m_input = input;
        m_beastObject = beastObject;
        this.itemNr = itemNr;

        timedParameter = (TimedParameter)input.get();

        addInputLabel();

        HBox boxHoriz;
        VBox boxVert;

        // Add elements specific to change times

        boxVert = FXUtils.newVBox();

        boxHoriz = FXUtils.newHBox();
        boxHoriz.getChildren().add(new Label("Number of elements:"));
        elementCountSpinner = new Spinner<>(0, Integer.MAX_VALUE, 0);
        boxHoriz.getChildren().add(elementCountSpinner);

        boxVert.getChildren().add(boxHoriz);

        elementsBox = FXUtils.newVBox();

        VBox timesBox = FXUtils.newVBox();
        timesEntryRow.getChildren().add(new Label("Element times:"));
        HBox timesBoxRow = FXUtils.newHBox();
        timesAreAgesCheckBox = new CheckBox("Times specified as ages");
        timesBoxRow.getChildren().add(timesAreAgesCheckBox);
        estimateTimesCheckBox = new CheckBox("Estimate times");
        timesBoxRow.getChildren().add(estimateTimesCheckBox);
        timesBox.getChildren().add(timesBoxRow);

        elementsBox.getChildren().add(timesBox);

        // Add elements specific to values

        VBox valuesBox = FXUtils.newVBox();
        boxHoriz = FXUtils.newHBox();
        boxHoriz.getChildren().add(new Label("Values:"));
        valuesTable = new TableView();
        VBox valuesTableBoxCol = FXUtils.newVBox();
        valuesTableBoxCol.getChildren().add(valuesTable);
        boxHoriz.getChildren().add(valuesTableBoxCol);

        valuesBox.getChildren().add(boxHoriz);

        boxHoriz = FXUtils.newHBox();
        scalarValuesCheckBox = new CheckBox("Scalar values");
        boxHoriz.getChildren().add(scalarValuesCheckBox);
        estimateValuesCheckBox = new CheckBox("Estimate values");
        boxHoriz.getChildren().add(estimateValuesCheckBox);

        valuesBox.getChildren().add(boxHoriz);
        elementsBox.getChildren().add(valuesBox);

        boxVert.getChildren().add(elementsBox);

        getChildren().add(boxVert);

        loadFromModel();

        // Add event listeners:
        elementCountSpinner.valueProperty().addListener(e -> saveToModel());
        timesAreAgesCheckBox.selectedProperty().addListener(e -> saveToModel());
        estimateTimesCheckBox.selectedProperty().addListener(e -> saveToModel());

        scalarValuesCheckBox.selectedProperty().addListener(e -> saveToModel());
        estimateValuesCheckBox.selectedProperty().addListener(e -> saveToModel());
    }

    /**
     * Ensures that the dimensions of the values RealParameter is consistent
     * with the number of types specified by the TypeSet.
     *
     * This is necessary because the TypeSet can be modified by other input
     * editors, such as the type trait set input editor.
     *
     * Once returning from this method, the number of types reported by the
     * SV itself should match the number reported by the typeset.
     */
    void ensureParamsConsistent() {

        timedParameter.typeSetInput.get().initAndValidate();

        int nTypes = timedParameter.typeSetInput.get().getNTypes();
        int nTimes = timedParameter.getTimeCount();

        if (nTimes==0)
            return;

        RealParameter valuesParam = (RealParameter)timedParameter.valuesInput.get();
        int valuesPerInterval = valuesParam.getDimension() / nTimes;

        if (valuesPerInterval == 1 || valuesPerInterval == nTypes) {
            timedParameter.initAndValidate();
            return;
        }

        StringBuilder valueBuilder = new StringBuilder();

        for (int timeIdx=0; timeIdx<nTimes; timeIdx++) {
            for (int typeIdx = 0; typeIdx < nTypes; typeIdx++) {
                valueBuilder.append(" ");

                if (typeIdx < valuesPerInterval)
                    valueBuilder.append(valuesParam.getValue(timeIdx*valuesPerInterval + typeIdx));
                else
                    valueBuilder.append(valuesParam.getValue(timeIdx*valuesPerInterval + (valuesPerInterval-1)));
            }
        }

        valuesParam.valuesInput.setValue(valueBuilder.toString(), valuesParam);
        valuesParam.initAndValidate();

        timedParameter.initAndValidate();
    }

    /**
     * Configure inputs for element times.
     *
     * @param nTimes number of times
     * @param parameter parameter containing times
     */
    void addTimesColumns(int nTimes, RealParameter parameter) {
        for (int i=0; i<nTimes; i++) {
            timesEntryRow.getChildren().add(new Label("Epoch " + (i+1)));
            TextField textField = new TextField(parameter.getValue(i).toString());

            textField.setPrefWidth(50);
            textField.setPadding(new Insets(0));
            HBox.setMargin(textField, new Insets(0, 10, 0, 0));

            int index = i;
            textField.setOnAction(event -> {
                parameter.setValue(index, Double.valueOf(textField.getText()));
            });

            timesEntryRow.getChildren().add(textField);
        }
    }

    /**
     * Populate GUI elements with values/dimensions from current BEASTObject model.
     * Called immediately after init(), and thus after every refreshPanel().
     */
    void loadFromModel() {

        ensureParamsConsistent();

        int nTimes = timedParameter.getTimeCount();
        int nTypes = timedParameter.getNTypes();

        elementCountSpinner.getValueFactory().setValue(nTimes);

        // Load times:

        if (nTimes > 0) {
            RealParameter timesParameter = (RealParameter)timedParameter.timesInput.get();
            addTimesColumns(nTimes, timesParameter);

            TableColumn<TimedParamValuesEntry, String> typeNameCol = new TableColumn<>("Types");
            typeNameCol.setCellValueFactory(p -> new ObservableValueBase<>() {
                @Override
                public String getValue() {
                    if (p.getValue().type < 0)
                        return "ALL";
                    else
                        return timedParameter.typeSetInput.get().getTypeName(p.getValue().type);
                }
            });
            valuesTable.getColumns().add(typeNameCol);

            for (int i = 0; i < nTimes; i++) {
                TableColumn<TimedParamValuesEntry, String> col =
                        new TableColumn<>("Epoch " + (i+1));
                int index = i;
                col.setCellValueFactory(p -> new ObservableValueBase<>() {
                    @Override
                    public String getValue() {
                        RealParameter param = p.getValue().parameter;
                        int type = p.getValue().type;
                        if (type<0)
                            return String.valueOf(param.getValue(index)); // scalar
                        else
                            return String.valueOf(param.getValue(index*nTypes + type));
                    }
                });

                col.setCellFactory(TextFieldTableCell.forTableColumn());
                col.setOnEditCommit(e -> {
                    RealParameter param = e.getRowValue().parameter;
                    int type = e.getRowValue().type;
                    double newValue = Double.parseDouble(e.getNewValue());
                    if (type<0)
                        param.setValue(index, newValue); // scalar
                    else
                        param.setValue(index*nTypes + type, newValue);
                });

                valuesTable.getColumns().add(col);
            }

            estimateTimesCheckBox.setSelected(timesParameter.isEstimatedInput.get());

            // Load values

            RealParameter valuesParameter = (RealParameter)timedParameter.valuesInput.get();

            if (valuesParameter.getDimension() == nTimes) {
                valuesTable.getItems().add(new TimedParamValuesEntry(valuesParameter, -1));
                if (nTypes > 1) {
                    scalarValuesCheckBox.setSelected(true);
                    scalarValuesCheckBox.disableProperty().set(true);
                } else {
                    scalarValuesCheckBox.setSelected(false);
                    scalarValuesCheckBox.disableProperty().set(false);
                }
            } else {
                for (int i=0; i<nTypes; i++)
                    valuesTable.getItems().add(new TimedParamValuesEntry(valuesParameter, i));

                scalarValuesCheckBox.setSelected(false);
            }

            estimateValuesCheckBox.setSelected(valuesParameter.isEstimatedInput.get());

        } else {
            estimateTimesCheckBox.setSelected(false);
            estimateValuesCheckBox.setSelected(false);

            elementsBox.setVisible(false);
        }

        timesAreAgesCheckBox.setSelected(timedParameter.timesAreAgesInput.get());

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

        int nTypes = timedParameter.getNTypes();
        int nTimes = (int) timeCountSpinnerModel.getValue();

        // Update table model dimensions

        valuesTableModel.setTimeCount(nTimes);
        valuesTableModel.setScalar(scalarValuesCheckBox.isSelected());

        timesTableModel.setColumnCount(nTimes);
        for (int colIdx = 0; colIdx< timesTable.getColumnCount(); colIdx++) {
            if (timesTableModel.getValueAt(0, colIdx) == null)
                timesTableModel.setValueAt(
                        colIdx > 0
                                ? timesTableModel.getValueAt(0, colIdx-1)
                                : 0.0,
                        0, colIdx);
        }

        // Save values and times

        RealParameter valuesParam = (RealParameter)timedParameter.valuesInput.get();
        RealParameter timesParam = (RealParameter)timedParameter.timesInput.get();
        if (nTimes>0) {
            if (valuesParam == null)
                valuesParam = getValuesParam();

            valuesParam.setDimension(valuesTableModel.getRowCount()*nTimes);
            valuesParam.valuesInput.setValue(valuesTableModel.getParameterString(), valuesParam);
            valuesParam.isEstimatedInput.setValue(estimateValuesCheckBox.isSelected(), valuesParam);

            if (timesParam == null)
                timesParam = getTimesParam();


            StringBuilder timesBuilder = new StringBuilder();
            for (int i=0; i<nTimes; i++)
                timesBuilder.append(" ").append(timesTableModel.getValueAt(0, i));

            timesParam.setDimension(nTimes);
            timesParam.valuesInput.setValue(timesBuilder.toString(), timesParam);
            timesParam.isEstimatedInput.setValue(estimateTimesCheckBox.isSelected(), timesParam);

        } else {
            if (timesParam != null)
                timesParam.isEstimatedInput.setValue(false, timesParam);
        }

        if (nTimes>0) {
            timedParameter.setInputValue("values", valuesParam);
            valuesParam.initAndValidate();

            timedParameter.setInputValue("times", timesParam);
            timesParam.initAndValidate();

        } else {
            timedParameter.setInputValue("values", null);
            timedParameter.setInputValue("times", null);
        }

        timedParameter.timesAreAgesInput.setValue(timesAreAgesCheckBox.isSelected(), timedParameter);

        timedParameter.initAndValidate();

        modelSaveInProcess = false;

        refreshPanel();
    }

    public static class TimedParamValuesEntry {

        public TimedParamValuesEntry(RealParameter parameter, int type) {
            this.parameter = parameter;
            this.type = type;
        }

        public RealParameter parameter;
        public int type;
    }

    RealParameter getValuesParam() {
        RealParameter changeTimesParam = (RealParameter)timedParameter.timesInput.get();
        if (changeTimesParam == null) {

            int idx = timedParameter.getID().indexOf("TP");
            String prefix = timedParameter.getID().substring(0, idx);
            String suffix = timedParameter.getID().substring(idx+2);
            String paramID = prefix + suffix;

            changeTimesParam = (RealParameter) doc.pluginmap.get(paramID);
            if (changeTimesParam == null) {
                changeTimesParam = new RealParameter("0.0");
                changeTimesParam.setID(paramID);
            }
        }

        return changeTimesParam;
    }

    RealParameter getTimesParam() {
        RealParameter changeTimesParam = (RealParameter)timedParameter.timesInput.get();
        if (changeTimesParam == null) {

            int idx = timedParameter.getID().indexOf("TP");
            String prefix = timedParameter.getID().substring(0, idx);
            String suffix = timedParameter.getID().substring(idx+2);
            String paramID = prefix + "Times" + suffix;

            changeTimesParam = (RealParameter) doc.pluginmap.get(paramID);
            if (changeTimesParam == null) {
                changeTimesParam = new RealParameter("0.0");
                changeTimesParam.setID(paramID);
            }
        }

        return changeTimesParam;
    }
}
