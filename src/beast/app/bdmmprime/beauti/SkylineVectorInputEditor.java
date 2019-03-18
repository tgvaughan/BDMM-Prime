package beast.app.bdmmprime.beauti;

import bdmmprime.parameterization.SkylineVectorParameter;
import beast.app.beauti.BeautiDoc;
import beast.app.draw.InputEditor;
import beast.core.BEASTInterface;
import beast.core.Input;
import beast.core.parameter.RealParameter;

import javax.swing.*;
import javax.swing.border.EtchedBorder;
import javax.swing.table.DefaultTableModel;
import java.awt.*;

public class SkylineVectorInputEditor extends InputEditor.Base {

    SkylineVectorParameter skylineVector;

    SpinnerModel changeCountSpinnerModel;
    JCheckBox scalarRatesCheckBox, timesAreAgesCheckBox;

    DefaultTableModel changeTimesTableModel;
    JTable changeTimesTable;
    Box changeTimesBox;

    DefaultTableModel valuesTableModel;
    JTable valuesTable;

    JCheckBox estimateValuesCheckBox, estimateTimesCheckBox;

    boolean modelSaveInProcess = false;

    public SkylineVectorInputEditor(BeautiDoc doc) {
        super(doc);
    }

    @Override
    public Class<?> type() {
        return SkylineVectorParameter.class;
    }

    @Override
    public void init(Input<?> input, BEASTInterface beastObject, int itemNr,
                     ExpandOption isExpandOption, boolean addButtons) {

        m_bAddButtons = addButtons;
        m_input = input;
        m_beastObject = beastObject;
        this.itemNr= itemNr;

        skylineVector = (SkylineVectorParameter)input.get();

        addInputLabel();

        Box boxVert, boxHoriz;

        boxVert = Box.createVerticalBox();
        boxVert.setBorder(new EtchedBorder());

        boxHoriz = Box.createHorizontalBox();
        JLabel changePointLabel = new JLabel("Number of change times:");
        JSpinner changePointsSpinner = new JSpinner();
        changeCountSpinnerModel =
                new SpinnerNumberModel(0, 0, Integer.MAX_VALUE,1);
        changePointsSpinner.setModel(changeCountSpinnerModel);
        boxHoriz.add(changePointLabel);
        boxHoriz.add(changePointsSpinner);

        boxVert.add(boxHoriz);

        changeTimesBox = Box.createVerticalBox();
        Box changeTimesBoxRow = Box.createHorizontalBox();
        changeTimesBoxRow.add(new JLabel("Change times:"));
        changeTimesTableModel = new DefaultTableModel(new String[] {"Epoch 1"}, 1);
        changeTimesTable = new JTable(changeTimesTableModel);
        changeTimesTable.setShowGrid(true);
        changeTimesTable.setGridColor(Color.GRAY);
        changeTimesTable.setCellSelectionEnabled(false);
        Box changeTimesTableBoxCol = Box.createVerticalBox();
        changeTimesTableBoxCol.add(changeTimesTable.getTableHeader());
        changeTimesTableBoxCol.add(changeTimesTable);
        changeTimesBoxRow.add(changeTimesTableBoxCol);
        changeTimesBox.add(changeTimesBoxRow);
        changeTimesBoxRow = Box.createHorizontalBox();
        timesAreAgesCheckBox = new JCheckBox("Times specified as ages");
        changeTimesBoxRow.add(timesAreAgesCheckBox);
        estimateTimesCheckBox = new JCheckBox("Estimate change times");
        changeTimesBoxRow.add(estimateTimesCheckBox);
        changeTimesBox.add(changeTimesBoxRow);

        boxVert.add(changeTimesBox);

        boxHoriz = Box.createHorizontalBox();
        boxHoriz.add(new JLabel("Values:"));
        valuesTableModel = new DefaultTableModel(1,1);
        valuesTable = new JTable(valuesTableModel);
        valuesTable.setShowGrid(true);
        valuesTable.setGridColor(Color.GRAY);
        valuesTable.setCellSelectionEnabled(false);
        boxHoriz.add(valuesTable);

        boxVert.add(boxHoriz);

        boxHoriz = Box.createHorizontalBox();
        scalarRatesCheckBox = new JCheckBox("Scalar values");
        boxHoriz.add(scalarRatesCheckBox, Component.LEFT_ALIGNMENT);
        estimateValuesCheckBox = new JCheckBox("Estimate values");
        boxHoriz.add(estimateValuesCheckBox);

        boxVert.add(boxHoriz);

        add(boxVert);

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
        int nTypes = skylineVector.typeSetInput.get().getNTypes();
        int nIntervals = skylineVector.getChangeCount() + 1;

        RealParameter valuesParam = skylineVector.rateValuesInput.get();
        int valuesPerInterval = valuesParam.getDimension() / nIntervals;

        if (valuesPerInterval == 1 || valuesPerInterval == nTypes) {
            skylineVector.initAndValidate();
            return;
        }

        StringBuilder valueBuilder = new StringBuilder();

        for (int interval=0; interval<nIntervals; interval++) {
            for (int typeIdx = 0; typeIdx < nTypes; typeIdx++) {
                valueBuilder.append(" ");

                if (typeIdx < valuesPerInterval)
                    valueBuilder.append(valuesParam.getValue(interval*valuesPerInterval + typeIdx));
                else
                    valueBuilder.append(valuesParam.getValue(interval*valuesPerInterval + (valuesPerInterval-1)));
            }
        }

        valuesParam.valuesInput.setValue(valueBuilder.toString(), valuesParam);
        valuesParam.initAndValidate();

        skylineVector.initAndValidate();
    }

    /**
     * Populate GUI elements with values/dimensions from current BEASTObject model.
     * Called immediately after init(), and thus after every refreshPanel().
     */
    void loadFromModel() {

        ensureParamsConsistent();

        int nChanges = skylineVector.getChangeCount();
        int nTypes = skylineVector.getNTypes();

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

            estimateTimesCheckBox.setEnabled(true);
            estimateTimesCheckBox.setSelected(skylineVector.changeTimesInput.get().isEstimatedInput.get());

            RealParameter changeTimesParameter = skylineVector.changeTimesInput.get();
            for (int i=0; i<changeTimesParameter.getDimension(); i++)
                changeTimesTableModel.setValueAt(changeTimesParameter.getValue(i), 0, i);

        } else {
            changeCountSpinnerModel.setValue(0);
            changeTimesTableModel.setColumnCount(0);
            changeTimesBox.setVisible(false);

            estimateTimesCheckBox.setEnabled(false);
        }
        timesAreAgesCheckBox.setEnabled(skylineVector.timesAreAgesInput.get());

        // Load values

        valuesTableModel.setColumnCount(nChanges+1);

        RealParameter valuesParameter = skylineVector.rateValuesInput.get();
        if (valuesParameter.getDimension()==(nChanges+1)) {
            if (nTypes>1) {
                scalarRatesCheckBox.setSelected(true);
                scalarRatesCheckBox.setEnabled(true);
            } else {
                scalarRatesCheckBox.setSelected(false);
                scalarRatesCheckBox.setEnabled(false);
            }
            valuesTableModel.setRowCount(1);

            for (int interval=0; interval<nChanges+1; interval++)
                valuesTableModel.setValueAt(valuesParameter.getValue(interval), 0, interval);
        } else {
            scalarRatesCheckBox.setSelected(false);
            valuesTableModel.setRowCount(nTypes);

            int valueIdx = 0;
            for (int interval=0; interval<nChanges+1; interval++) {
                for (int typeIdx = 0; typeIdx<nTypes; typeIdx++) {
                    valuesTableModel.setValueAt(valuesParameter.getValue(valueIdx++), typeIdx, interval);
                }
            }
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

        int nTypes = skylineVector.getNTypes();
        int nChanges = (int)changeCountSpinnerModel.getValue();

        // Update table model dimensions

        valuesTableModel.setColumnCount(nChanges+1);
        if (scalarRatesCheckBox.isSelected()) {
            valuesTableModel.setRowCount(1);
        } else {
            valuesTableModel.setRowCount(nTypes);
        }
        for (int rowIdx=0; rowIdx<valuesTableModel.getRowCount(); rowIdx++) {
            for (int colIdx=0; colIdx<valuesTableModel.getColumnCount(); colIdx++) {
                if (valuesTableModel.getValueAt(rowIdx, colIdx) == null)
                    valuesTableModel.setValueAt(
                            colIdx > 0
                                    ? valuesTableModel.getValueAt(rowIdx, colIdx-1)
                                    : valuesTableModel.getValueAt(rowIdx-1, colIdx),
                            rowIdx, colIdx);
            }
        }

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

        RealParameter rateValuesParam = skylineVector.rateValuesInput.get();
        rateValuesParam.setDimension(
                valuesTableModel.getRowCount()*valuesTableModel.getColumnCount());
        StringBuilder valuesBuilder = new StringBuilder();
        for (int colIdx=0; colIdx<valuesTableModel.getColumnCount(); colIdx++) {
            for (int rowIdx=0; rowIdx<valuesTableModel.getRowCount(); rowIdx++) {
                valuesBuilder.append(" ").append(valuesTableModel.getValueAt(rowIdx, colIdx));
            }
        }
        rateValuesParam.valuesInput.setValue(valuesBuilder.toString(), rateValuesParam);
        rateValuesParam.isEstimatedInput.setValue(estimateValuesCheckBox.isSelected(), rateValuesParam);

        // Save change times

        RealParameter changeTimesParam = skylineVector.changeTimesInput.get();
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

        skylineVector.setInputValue("rateValues", rateValuesParam);
        rateValuesParam.initAndValidate();

        if (nChanges>0) {
            skylineVector.setInputValue("changeTimes", changeTimesParam);
            changeTimesParam.initAndValidate();
        } else {
            skylineVector.setInputValue("changeTimes", null);
        }

        skylineVector.initAndValidate();

        modelSaveInProcess = false;

        refreshPanel();
    }

    String getChangeTimesParameterID() {
        int idx = skylineVector.getID().indexOf("SV");
        String prefix = skylineVector.getID().substring(0, idx);
        String suffix = skylineVector.getID().substring(idx+2);

        return prefix + "ChangeTimes" + suffix;
    }
}
