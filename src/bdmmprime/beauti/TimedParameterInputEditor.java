package bdmmprime.beauti;

import bdmmprime.parameterization.TimedParameter;
import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.inputeditor.InputEditor;
import beastfx.app.util.FXUtils;
import javafx.beans.binding.Bindings;
import javafx.beans.value.ObservableValueBase;
import javafx.geometry.Insets;
import javafx.scene.control.*;
import javafx.scene.control.cell.TextFieldTableCell;
import javafx.scene.layout.*;
import javafx.scene.paint.Color;
import javafx.util.converter.DoubleStringConverter;

import java.util.Arrays;
import java.util.stream.Collectors;

public class TimedParameterInputEditor extends InputEditor.Base {

    TimedParameter timedParameter;

    TableView<TimedParamValuesTableEntry> valuesTable;

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
        pane = FXUtils.newHBox();

        timedParameter = (TimedParameter)input.get();

        addInputLabel();

        ensureValuesConsistency(true);

        // Add elements specific to times

        int timeCount = timedParameter.getTimeCount();

        VBox mainInputBox = FXUtils.newVBox();
        mainInputBox.setBorder(new Border(new BorderStroke(Color.LIGHTGRAY,
                BorderStrokeStyle.SOLID, null, null)));

        HBox boxHoriz = FXUtils.newHBox();
        Label changePointLabel = new Label("Number of times:");
        Spinner<Integer> timeCountSpinner = new Spinner<>(0, Integer.MAX_VALUE, timeCount);
        timeCountSpinner.setEditable(true);
        boxHoriz.getChildren().add(changePointLabel);
        boxHoriz.getChildren().add(timeCountSpinner);

        mainInputBox.getChildren().add(boxHoriz);

        VBox timesAndValuesBox = FXUtils.newVBox();
        HBox timesEntryRow = FXUtils.newHBox();
        timesAndValuesBox.getChildren().add(timesEntryRow);
        HBox timesBoxRow = FXUtils.newHBox();
        CheckBox timesAreAgesCheckBox = new CheckBox("Times specified as ages");
        timesBoxRow.getChildren().add(timesAreAgesCheckBox);
        CheckBox estimateTimesCheckBox = new CheckBox("Estimate change times");
        timesBoxRow.getChildren().add(estimateTimesCheckBox);
        timesAndValuesBox.getChildren().add(timesBoxRow);



        // Add elements specific to values

        boxHoriz = FXUtils.newHBox();
        boxHoriz.getChildren().add(new Label("Values:"));
        valuesTable = new TableView<>();
        valuesTable.getSelectionModel().setCellSelectionEnabled(true);
        valuesTable.setEditable(true);
        valuesTable.setFixedCellSize(25);
        valuesTable.prefHeightProperty().bind(valuesTable.fixedCellSizeProperty()
                .multiply(Bindings.size(valuesTable.getItems()).add(1.1)));
        VBox valuesTableBoxCol = FXUtils.newVBox();
        valuesTableBoxCol.getChildren().add(valuesTable);
        boxHoriz.getChildren().add(valuesTableBoxCol);

        timesAndValuesBox.getChildren().add(boxHoriz);

        boxHoriz = FXUtils.newHBox();
        CheckBox scalarValues = new CheckBox("Scalar values");
        boxHoriz.getChildren().add(scalarValues);
        CheckBox estimateValuesCheckBox = new CheckBox("Estimate values");
        boxHoriz.getChildren().add(estimateValuesCheckBox);

        int nTypes = timedParameter.getNTypes();
        if (nTypes > 1) {
            scalarValues.setSelected(true);
            scalarValues.disableProperty().set(false);
        } else {
            scalarValues.setSelected(false);
            scalarValues.disableProperty().set(true);
        }

        timesAndValuesBox.getChildren().add(boxHoriz);
        mainInputBox.getChildren().add(timesAndValuesBox);

        if (timeCount > 0) {
            updateTimesUI(timesEntryRow);
            updateValuesUI();

            timesAreAgesCheckBox.setSelected(timedParameter.timesAreAgesInput.get());

            estimateTimesCheckBox.setSelected(
                    ((RealParameter) timedParameter.timesInput.get())
                            .isEstimatedInput.get());

            estimateValuesCheckBox.setSelected(
                    ((RealParameter) timedParameter.valuesInput.get())
                            .isEstimatedInput.get());
        } else {
            timesAndValuesBox.setVisible(false);
            timesAndValuesBox.setManaged(false);
        }

        pane.getChildren().add(mainInputBox);
        getChildren().add(pane);

        // Add event listeners:
        timeCountSpinner.valueProperty().addListener((observable, oldValue, newValue) -> {
            System.out.println(oldValue + " -> " + newValue);

            if (newValue > 0) {
                RealParameter times = getTimesParam();
                times.setDimension(newValue);
                sanitiseRealParameter(times);
                timedParameter.timesInput.setValue(times, timedParameter);
            } else {
                if (estimateTimesCheckBox.isSelected())
                    estimateTimesCheckBox.fire();

                if (estimateValuesCheckBox.isSelected())
                    estimateValuesCheckBox.fire();
                timedParameter.timesInput.setValue(null, timedParameter);
            }

            ensureValuesConsistency(scalarValues.isSelected());

            if (newValue > 0) {
                updateTimesUI(timesEntryRow);
                updateValuesUI();


                timesAreAgesCheckBox.setSelected(timedParameter.timesAreAgesInput.get());

                estimateTimesCheckBox.setSelected(getTimesParam().isEstimatedInput.get());

                timesAndValuesBox.setManaged(true);
                timesAndValuesBox.setVisible(true);


            } else {
                timesAndValuesBox.setManaged(false);
                timesAndValuesBox.setVisible(false);
            }

            System.out.println(timedParameter);
            System.out.println(scalarValues.isSelected());
        });

        timesAreAgesCheckBox.selectedProperty().addListener((observable, oldValue, newValue) -> {
            timedParameter.timesAreAgesInput.setValue(newValue, timedParameter);
            timedParameter.initAndValidate();
            System.out.println(timedParameter);
        });

        estimateTimesCheckBox.selectedProperty().addListener((observable, oldValue, newValue) -> {
            RealParameter changeTimes = (RealParameter) timedParameter.timesInput.get();
            changeTimes.isEstimatedInput.setValue(newValue, changeTimes);
            sync();
        });

        scalarValues.selectedProperty().addListener((observable, oldValue, newValue) -> {
            if (newValue)
                getValuesParam().setDimension(timedParameter.getTimeCount());
            else
                getValuesParam().setDimension(timedParameter.getNTypes() * timedParameter.getTimeCount());

            sanitiseRealParameter(getValuesParam());
            ensureValuesConsistency(newValue);
            updateValuesUI();
            System.out.println(timedParameter);
        });

        estimateValuesCheckBox.selectedProperty().addListener((observable, oldValue, newValue) -> {
            RealParameter valuesParam = getValuesParam();
            valuesParam.isEstimatedInput.setValue(newValue, valuesParam);
            hardSync();
        });


    }

    /**
     * Configure inputs for times.
     *
     * @param changeTimesEntryRow HBox containing time inputs
     */
    void updateTimesUI(HBox changeTimesEntryRow) {
        RealParameter parameter = getTimesParam();
        changeTimesEntryRow.getChildren().clear();
        changeTimesEntryRow.getChildren().add(new Label("Times:"));
        for (int i=0; i<parameter.getDimension(); i++) {
            changeTimesEntryRow.getChildren().add(new Label("Time " + (i+1) + ": "));
            TextField textField = new TextField(parameter.getValue(i).toString());

            textField.setPrefWidth(50);
            textField.setPadding(new Insets(0));
            HBox.setMargin(textField, new Insets(0, 10, 0, 0));

            int index = i;
            textField.textProperty().addListener((observable, oldValue, newValue) -> {
                parameter.setValue(index, Double.valueOf(newValue));
                sanitiseRealParameter(parameter);
                timedParameter.initAndValidate();
                System.out.println(timedParameter);
            });

            changeTimesEntryRow.getChildren().add(textField);
        }
    }

    void updateValuesUI() {
        valuesTable.getColumns().clear();
        valuesTable.getItems().clear();

        int nTimes = timedParameter.getTimeCount();
        int nTypes = timedParameter.getNTypes();

        RealParameter valuesParameter = getValuesParam();
        TableColumn<TimedParamValuesTableEntry, String> typeCol = new TableColumn<>("Type");
        typeCol.setCellValueFactory(p -> new ObservableValueBase<>() {
            @Override
            public String getValue() {
                int type = p.getValue().type;
                return type < 0
                        ? "ALL"
                        : timedParameter.typeSetInput.get().getTypeName(type);
            }
        });
        valuesTable.getColumns().add(typeCol);
        for (int i=0; i<nTimes; i++) {
            TableColumn<TimedParamValuesTableEntry, Double> col = new TableColumn<>("Epoch " + (i+1));
            int epochIdx = i;
            col.setCellValueFactory(p -> new ObservableValueBase<>() {
                @Override
                public Double getValue() {
                    int type = p.getValue().type;
                    return type<0
                            ? valuesParameter.getValue(epochIdx)
                            : valuesParameter.getValue(epochIdx*nTypes + type);
                }
            });
            col.setCellFactory(TextFieldTableCell.forTableColumn(new DoubleStringConverter()));
            col.setOnEditCommit(e -> {
                int type = e.getTableView()
                        .getItems().get(e.getTablePosition().getRow()).type;
                if (type < 0) {
                    valuesParameter.setValue(epochIdx, e.getNewValue());
                } else {
                    valuesParameter.setValue(epochIdx*nTypes + type, e.getNewValue());
                }
                sanitiseRealParameter(valuesParameter);
                System.out.println(timedParameter);
            });
            valuesTable.getColumns().add(col);
        }

        if (valuesParameter.getDimension() / nTimes > 1) {
            for (int type=0; type<nTypes; type++)
                valuesTable.getItems().add(new TimedParamValuesTableEntry(type));
        } else {
            valuesTable.getItems().add(new TimedParamValuesTableEntry(-1));
        }
    }

    void ensureValuesConsistency(boolean scalar) {
        int nTypes = timedParameter.typeSetInput.get().getNTypes();
        int nEpochs = timedParameter.timesInput.get() == null
                ? 0
                : timedParameter.timesInput.get().getDimension();

        System.out.println("Number of epochs: " + nEpochs);

        if (nEpochs > 0) {
            RealParameter valuesParam = getValuesParam();

            if (scalar)
                valuesParam.setDimension(nEpochs);
            else
                valuesParam.setDimension(nTypes * nEpochs);

            sanitiseRealParameter(valuesParam);
            timedParameter.valuesInput.setValue(valuesParam, timedParameter);
        } else
            timedParameter.valuesInput.setValue(null, timedParameter);

        timedParameter.initAndValidate();
    }

    void sanitiseRealParameter(RealParameter parameter) {
        parameter.valuesInput.setValue(
                Arrays.stream(parameter.getDoubleValues())
                        .mapToObj(String::valueOf)
                        .collect(Collectors.joining(" ")),
                parameter);
        parameter.initAndValidate();
    }

    RealParameter getTimesParam() {
        RealParameter timesParam = (RealParameter)timedParameter.timesInput.get();
        if (timesParam == null) {

            int idx = timedParameter.getID().indexOf("SP");
            String prefix = timedParameter.getID().substring(0, idx);
            String suffix = timedParameter.getID().substring(idx+2);
            String paramID = prefix + "Times" + suffix;

            timesParam = (RealParameter) doc.pluginmap.get(paramID);
            if (timesParam == null) {
                timesParam = new RealParameter("0.0");
                timesParam.setID(paramID);
            }
        }
        return timesParam;
    }

    RealParameter getValuesParam() {
        RealParameter valuesParam = (RealParameter)timedParameter.valuesInput.get();
        if (valuesParam == null) {

            int idx = timedParameter.getID().indexOf("SP");
            String prefix = timedParameter.getID().substring(0, idx);
            String suffix = timedParameter.getID().substring(idx+2);
            String paramID = prefix + suffix;

            valuesParam = (RealParameter) doc.pluginmap.get(paramID);
            if (valuesParam == null) {
                valuesParam = new RealParameter("0.0");
                valuesParam.isEstimatedInput.setValue(true, valuesParam);
                valuesParam.setID(paramID);
            }
        }
        return valuesParam;
    }

    public static class TimedParamValuesTableEntry {
        int type;

        public TimedParamValuesTableEntry(int type) {
            this.type = type;
        }
    }
}
