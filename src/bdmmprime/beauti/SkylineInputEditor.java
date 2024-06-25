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
import javafx.geometry.Insets;
import javafx.scene.control.*;
import javafx.scene.layout.*;
import javafx.scene.paint.Color;
import javafx.scene.text.Font;
import javafx.util.Duration;

import java.util.Arrays;
import java.util.stream.Collectors;

public abstract class SkylineInputEditor extends InputEditor.Base {

    SkylineParameter skylineParameter;

    TableView<ValuesTableEntry> valuesTable;

    VBox mainInputBox;

    EpochVisualizerPane epochVisualizer;

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
        pane = FXUtils.newHBox();

        skylineParameter = (SkylineParameter) input.get();

        ensureValuesConsistency();

        addInputLabel();

        // Add elements specific to change times

        int nChanges = skylineParameter.getChangeCount();

        mainInputBox = FXUtils.newVBox();
        mainInputBox.setBorder(new Border(new BorderStroke(Color.LIGHTGRAY,
                BorderStrokeStyle.SOLID, null, null)));

        HBox boxHoriz = FXUtils.newHBox();
        Label changePointLabel = new Label("Number of change times:");
        Spinner<Integer> changeCountSpinner = new Spinner<>(0, Integer.MAX_VALUE, nChanges);
        changeCountSpinner.setEditable(true);
        changeCountSpinner.setRepeatDelay(Duration.INDEFINITE); // (Hack around weird race condition I can't solve)
        boxHoriz.getChildren().add(changePointLabel);
        boxHoriz.getChildren().add(changeCountSpinner);

        mainInputBox.getChildren().add(boxHoriz);

        VBox changeTimesBox = FXUtils.newVBox();
        HBox changeTimesEntryRow = FXUtils.newHBox();
        changeTimesBox.getChildren().add(changeTimesEntryRow);
        HBox changeTimesBoxRow = FXUtils.newHBox();
        CheckBox timesAreAgesCheckBox = new CheckBox("Times specified as ages");
        changeTimesBoxRow.getChildren().add(timesAreAgesCheckBox);
        CheckBox estimateTimesCheckBox = new CheckBox("Estimate change times");
        changeTimesBoxRow.getChildren().add(estimateTimesCheckBox);
        changeTimesBox.getChildren().add(changeTimesBoxRow);

        changeTimesBoxRow = FXUtils.newHBox();
        CheckBox timesAreRelativeCheckBox = new CheckBox("Relative to process length");
        changeTimesBoxRow.getChildren().add(timesAreRelativeCheckBox);
        Button distributeChangeTimesButton = new Button("Distribute evenly");
        changeTimesBoxRow.getChildren().add(distributeChangeTimesButton);
        changeTimesBox.getChildren().add(changeTimesBoxRow);

        mainInputBox.getChildren().add(changeTimesBox);

        if (nChanges > 0) {
            updateChangeTimesUI((RealParameter) skylineParameter.changeTimesInput.get(),
                    changeTimesEntryRow);
            timesAreAgesCheckBox.setSelected(skylineParameter.timesAreAgesInput.get());
            timesAreRelativeCheckBox.setSelected(skylineParameter.timesAreRelativeInput.get());

            estimateTimesCheckBox.setSelected(
                    ((RealParameter) skylineParameter.changeTimesInput.get())
                            .isEstimatedInput.get());
        } else {
            changeTimesBox.setVisible(false);
            changeTimesBox.setManaged(false);
        }

        // Add elements specific to values

        RealParameter valuesParameter = (RealParameter) skylineParameter.skylineValuesInput.get();

        boxHoriz = FXUtils.newHBox();
        boxHoriz.getChildren().add(new Label("Values:"));
        valuesTable = new TableView<>();
        valuesTable.getSelectionModel().setCellSelectionEnabled(true);
        valuesTable.setEditable(true);
        VBox valuesTableBoxCol = FXUtils.newVBox();
        valuesTableBoxCol.getChildren().add(valuesTable);
        boxHoriz.getChildren().add(valuesTableBoxCol);

        mainInputBox.getChildren().add(boxHoriz);

        boxHoriz = FXUtils.newHBox();
        CheckBox scalarRatesCheckBox = new CheckBox("Scalar values");
        boxHoriz.getChildren().add(scalarRatesCheckBox);
        CheckBox estimateValuesCheckBox = new CheckBox("Estimate values");
        estimateValuesCheckBox.setSelected(valuesParameter.isEstimatedInput.get());
        boxHoriz.getChildren().add(estimateValuesCheckBox);

        mainInputBox.getChildren().add(boxHoriz);

        boxHoriz = FXUtils.newHBox();
        CheckBox visualizerCheckBox = new CheckBox("Display visualization");
        visualizerCheckBox.setSelected(skylineParameter.epochVisualizerDisplayed);
        boxHoriz.getChildren().add(visualizerCheckBox);
        mainInputBox.getChildren().add(boxHoriz);

        epochVisualizer = new EpochVisualizerPane(getTree(), getTypeTraitSet(), skylineParameter);
        epochVisualizer.widthProperty().bind(mainInputBox.widthProperty().subtract(10));
        epochVisualizer.setVisible(skylineParameter.epochVisualizerDisplayed);
        epochVisualizer.setManaged(skylineParameter.epochVisualizerDisplayed);
        mainInputBox.getChildren().add(epochVisualizer);

        int nTypes = skylineParameter.getNTypes();
        if (valuesParameter.getDimension() == (nChanges + 1)) {
            if (nTypes > 1) {
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

        pane.getChildren().add(mainInputBox);
        getChildren().add(pane);

        // Add event listeners:
        changeCountSpinner.valueProperty().addListener((observable, oldValue, newValue) -> {
            System.out.println(oldValue + " -> " + newValue);

            if (newValue > 0) {
                RealParameter param = (RealParameter) skylineParameter.changeTimesInput.get();
                if (param == null) {
                    if (!doc.pluginmap.containsKey(getChangeTimesParameterID())) {
                        param = new RealParameter("0.0");
                        param.setID(getChangeTimesParameterID());
                    } else {
                        param = (RealParameter) doc.pluginmap.get(getChangeTimesParameterID());
                    }
                    skylineParameter.changeTimesInput.setValue(param, skylineParameter);
                }
                param.setDimension(newValue);
                sanitiseRealParameter(param);
            } else {
                skylineParameter.changeTimesInput.setValue(null, skylineParameter);
            }

            ensureValuesConsistency();

            if (newValue > 0) {
                updateChangeTimesUI((RealParameter) skylineParameter.changeTimesInput.get(),
                        changeTimesEntryRow);
                timesAreAgesCheckBox.setSelected(skylineParameter.timesAreAgesInput.get());

                estimateTimesCheckBox.setSelected(
                        ((RealParameter) skylineParameter.changeTimesInput.get())
                                .isEstimatedInput.get());

                changeTimesBox.setManaged(true);
                changeTimesBox.setVisible(true);
            } else {
                changeTimesBox.setManaged(false);
                changeTimesBox.setVisible(false);
            }

            System.out.println(skylineParameter);
            System.out.println(scalarRatesCheckBox.isSelected());

            updateValuesUI();
        });

        timesAreAgesCheckBox.selectedProperty().addListener((observable, oldValue, newValue) -> {
            skylineParameter.timesAreAgesInput.setValue(newValue, skylineParameter);
            skylineParameter.initAndValidate();
            System.out.println(skylineParameter);
            epochVisualizer.repaintCanvas();
        });

        estimateTimesCheckBox.selectedProperty().addListener((observable, oldValue, newValue) -> {
            RealParameter changeTimes = (RealParameter) skylineParameter.changeTimesInput.get();
            changeTimes.isEstimatedInput.setValue(newValue, changeTimes);
            sync();
        });

        timesAreRelativeCheckBox.selectedProperty().addListener((observable, oldValue, newValue) -> {
            skylineParameter.timesAreRelativeInput.setValue(newValue, skylineParameter);
            skylineParameter.initAndValidate();
            System.out.println(skylineParameter);
            epochVisualizer.repaintCanvas();
        });

        distributeChangeTimesButton.setOnAction(e -> {

            RealParameter changeTimesParam = (RealParameter) skylineParameter.changeTimesInput.get();
            int nTimes = changeTimesParam.getDimension();

            if (skylineParameter.timesAreRelativeInput.get()) {
                for (int i = 0; i < nTimes; i++) {
                    changeTimesParam.setValue(i, ((double) (i + 1)) / (nTimes + 1));
                }
            } else {
                if (nTimes > 1) {
                    for (int i = 0; i < nTimes - 1; i++) {
                        changeTimesParam.setValue(i,
                                (changeTimesParam.getArrayValue(nTimes - 1) * (i + 1)) / (nTimes + 1));
                    }
                }
            }

            sanitiseRealParameter(changeTimesParam);
            sync();
        });

        scalarRatesCheckBox.selectedProperty().addListener((observable, oldValue, newValue) -> {
            skylineParameter.isScalarInput.setValue(newValue, skylineParameter);

            ensureValuesConsistency();
            sanitiseRealParameter(valuesParameter);
            updateValuesUI();
            System.out.println(skylineParameter);
            epochVisualizer.setScalar(newValue);
        });

        estimateValuesCheckBox.selectedProperty().addListener((observable, oldValue, newValue) -> {
                valuesParameter.isEstimatedInput.setValue(newValue, valuesParameter);
                sync();
        });

        visualizerCheckBox.selectedProperty().addListener((observable, oldValue, newValue) -> {
            skylineParameter.epochVisualizerDisplayed = newValue;
            epochVisualizer.setVisible(newValue);
            epochVisualizer.setManaged(newValue);
        });
    }

    /**
     * Configure inputs for change times.
     *
     * @param parameter change times parameter
     * @param changeTimesEntryRow HBox containing time inputs
     */
    void updateChangeTimesUI(RealParameter parameter,
                             HBox changeTimesEntryRow) {
        changeTimesEntryRow.getChildren().clear();
        changeTimesEntryRow.getChildren().add(new Label("Change times:"));
        for (int i=0; i<parameter.getDimension(); i++) {
            changeTimesEntryRow.getChildren().add(new Label("Epoch " + (i+1) + "->" + (i+2) + ": "));
            TextField textField = new TextField(parameter.getValue(i).toString());

            textField.setPrefWidth(50);
            textField.setPadding(new Insets(0));
            HBox.setMargin(textField, new Insets(0, 10, 0, 0));

            int index = i;
            textField.textProperty().addListener((observable, oldValue, newValue) -> {
                parameter.setValue(index, Double.valueOf(newValue));
                sanitiseRealParameter(parameter);
                skylineParameter.initAndValidate();
                System.out.println(skylineParameter);
                epochVisualizer.repaintCanvas();
            });

            changeTimesEntryRow.getChildren().add(textField);
        }
    }

    void sanitiseRealParameter(RealParameter parameter) {
        parameter.valuesInput.setValue(
                Arrays.stream(parameter.getDoubleValues())
                        .mapToObj(String::valueOf)
                        .collect(Collectors.joining(" ")),
                parameter);
        parameter.initAndValidate();
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

    /**
     * Called on initialisation to ensure values parameter is
     * up-to-date with the number of types.  This is necessary because
     * the type count is affected by the TypeTraitSetInputEditor, which
     * calls a sync() when this value changes.
     */
    abstract void ensureValuesConsistency();

    abstract void updateValuesUI();

    public abstract static class ValuesTableEntry { }
}
