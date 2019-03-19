package beast.app.bdmmprime.beauti;

import bdmmprime.parameterization.SkylineParameter;
import beast.app.beauti.BeautiDoc;
import beast.app.draw.InputEditor;
import beast.core.BEASTInterface;
import beast.core.Input;
import beast.core.parameter.RealParameter;

import javax.swing.*;
import javax.swing.border.EtchedBorder;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.DefaultTableModel;
import java.awt.*;

public abstract class SkylineInputEditor extends InputEditor.Base {

    SkylineParameter skylineParameter;

    SpinnerModel changeCountSpinnerModel;
    JCheckBox scalarRatesCheckBox, timesAreAgesCheckBox;

    DefaultTableModel changeTimesTableModel;
    JTable changeTimesTable;
    Box changeTimesBox;

    SkylineValuesTableModel valuesTableModel;
    JTable valuesTable;

    JCheckBox estimateValuesCheckBox, estimateTimesCheckBox;

    Box mainInputBox;

    boolean modelSaveInProcess = false;

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

        Box boxHoriz;

        // Add elements specific to change times

        mainInputBox = Box.createVerticalBox();
        mainInputBox.setBorder(new EtchedBorder());

        boxHoriz = Box.createHorizontalBox();
        JLabel changePointLabel = new JLabel("Number of change times:");
        JSpinner changePointsSpinner = new JSpinner();
        changeCountSpinnerModel =
                new SpinnerNumberModel(0, 0, Integer.MAX_VALUE,1);
        changePointsSpinner.setModel(changeCountSpinnerModel);
        boxHoriz.add(changePointLabel);
        boxHoriz.add(changePointsSpinner);
        boxHoriz.add(makeHorizontalFiller());

        mainInputBox.add(boxHoriz);

        changeTimesBox = Box.createVerticalBox();
        Box changeTimesBoxRow = Box.createHorizontalBox();
        changeTimesBoxRow.add(new JLabel("Change times:"));
        changeTimesTableModel = new DefaultTableModel(new String[] {"Epoch 1->2"}, 1);
        changeTimesTable = new JTable(changeTimesTableModel);
        changeTimesTable.setShowGrid(true);
        changeTimesTable.setGridColor(Color.GRAY);
        changeTimesTable.setCellSelectionEnabled(false);
        Box changeTimesTableBoxCol = Box.createVerticalBox();
        changeTimesTableBoxCol.add(changeTimesTable.getTableHeader());
        changeTimesTableBoxCol.add(changeTimesTable);
        changeTimesBoxRow.add(changeTimesTableBoxCol);
        changeTimesBoxRow.add(makeHorizontalFiller());
        changeTimesBox.add(changeTimesBoxRow);
        changeTimesBoxRow = Box.createHorizontalBox();
        timesAreAgesCheckBox = new JCheckBox("Times specified as ages");
        changeTimesBoxRow.add(timesAreAgesCheckBox);
        estimateTimesCheckBox = new JCheckBox("Estimate change times");
        changeTimesBoxRow.add(estimateTimesCheckBox);
        changeTimesBoxRow.add(makeHorizontalFiller());
        changeTimesBox.add(changeTimesBoxRow);

        mainInputBox.add(changeTimesBox);

        // Add elements specific to values

        boxHoriz = Box.createHorizontalBox();
        boxHoriz.add(new JLabel("Values:"));
        valuesTableModel = getValuesTableModel();
        valuesTable = new JTable(valuesTableModel);
        valuesTable.setShowGrid(true);
        valuesTable.setGridColor(Color.GRAY);
        valuesTable.setCellSelectionEnabled(false);
        Box valuesTableBoxCol = Box.createVerticalBox();
        valuesTableBoxCol.add(valuesTable.getTableHeader());
        valuesTableBoxCol.add(valuesTable);
        boxHoriz.add(valuesTableBoxCol);
        boxHoriz.add(makeHorizontalFiller());

        mainInputBox.add(boxHoriz);

        boxHoriz = Box.createHorizontalBox();
        scalarRatesCheckBox = new JCheckBox("Scalar values");
        boxHoriz.add(scalarRatesCheckBox, Component.LEFT_ALIGNMENT);
        estimateValuesCheckBox = new JCheckBox("Estimate values");
        boxHoriz.add(estimateValuesCheckBox);
        boxHoriz.add(makeHorizontalFiller());

        mainInputBox.add(boxHoriz);

        add(mainInputBox);

        ensureParamsConsistent();
        loadFromModel();

        // Add event listeners:
        changeCountSpinnerModel.addChangeListener(e -> saveToModel());
        changeTimesTableModel.addTableModelListener(e -> saveToModel());
        timesAreAgesCheckBox.addItemListener(e -> saveToModel());
        estimateTimesCheckBox.addItemListener(e -> saveToModel());

        valuesTableModel.addTableModelListener(e -> saveToModel());
        scalarRatesCheckBox.addItemListener(e -> saveToModel());
        estimateValuesCheckBox.addItemListener(e -> saveToModel());
    }

    abstract SkylineValuesTableModel getValuesTableModel();

    /**
     * @return a horizontal filler object.
     */
    public static Box.Filler makeHorizontalFiller() {
        return new Box.Filler(new Dimension(1,1),
                new Dimension(1,1),
                new Dimension(Integer.MAX_VALUE,1));
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
     * Populate GUI elements with values/dimensions from current BEASTObject model.
     * Called immediately after init(), and thus after every refreshPanel().
     */
    void loadFromModel() {

        int nChanges = skylineParameter.getChangeCount();
        int nTypes = skylineParameter.getNTypes();

        // Load change times:

        if (nChanges > 0) {
            changeCountSpinnerModel.setValue(nChanges);
            changeTimesTableModel.setColumnCount(nChanges);

            String[] columnNames = new String[nChanges];
            for (int i=0; i<nChanges; i++) {
                columnNames[i] = "Epoch " + (i+1) + " -> " + (i+2);
            }
            changeTimesTableModel.setColumnIdentifiers(columnNames);

            changeTimesBox.setVisible(true);

            estimateTimesCheckBox.setSelected(skylineParameter.changeTimesInput.get().isEstimatedInput.get());

            RealParameter changeTimesParameter = skylineParameter.changeTimesInput.get();
            for (int i=0; i<changeTimesParameter.getDimension(); i++)
                changeTimesTableModel.setValueAt(changeTimesParameter.getValue(i), 0, i);

        } else {
            changeCountSpinnerModel.setValue(0);
            changeTimesTableModel.setColumnCount(0);
            changeTimesBox.setVisible(false);
        }

        timesAreAgesCheckBox.setSelected(skylineParameter.timesAreAgesInput.get());

        // Load values

        RealParameter valuesParameter = skylineParameter.skylineValuesInput.get();

        valuesTableModel.setIntervalCount(nChanges+1);
        valuesTableModel.loadFromParameter(valuesParameter);

        if (valuesParameter.getDimension()==(nChanges+1)) {
            if (nTypes>1) {
                scalarRatesCheckBox.setSelected(true);
                scalarRatesCheckBox.setEnabled(true);
            } else {
                scalarRatesCheckBox.setSelected(false);
                scalarRatesCheckBox.setEnabled(false);
            }
        } else {
            scalarRatesCheckBox.setSelected(false);
        }

        estimateValuesCheckBox.setSelected(valuesParameter.isEstimatedInput.get());
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

        int nChanges = (int)changeCountSpinnerModel.getValue();

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

        RealParameter skylineValuesParam = skylineParameter.skylineValuesInput.get();
        skylineValuesParam.setDimension(valuesTableModel.getRowCount()*(nChanges+1));
        skylineValuesParam.valuesInput.setValue(valuesTableModel.getParameterString(), skylineValuesParam);
        skylineValuesParam.isEstimatedInput.setValue(estimateValuesCheckBox.isSelected(), skylineValuesParam);
        skylineParameter.setInputValue("skylineValues", skylineValuesParam);
        skylineValuesParam.initAndValidate();

        // Save change times

        RealParameter changeTimesParam = skylineParameter.changeTimesInput.get();
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

        skylineParameter.initAndValidate();

        modelSaveInProcess = false;

        refreshPanel();
    }

    abstract String getChangeTimesParameterID();
}
