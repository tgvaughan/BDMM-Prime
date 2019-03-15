package beast.app.bdmmprime.beauti;

import bdmmprime.parameterization.SkylineVectorParameter;
import beast.app.beauti.BeautiDoc;
import beast.app.beauti.PartitionContext;
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
        changeTimesTableModel = new DefaultTableModel(1,1);
        changeTimesTable = new JTable(changeTimesTableModel);
        changeTimesTable.setShowGrid(true);
        changeTimesBoxRow.add(changeTimesTable);
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

    void loadFromModel() {

        int nChanges = skylineVector.getChangeCount();

        // Load changepoints:
        if (nChanges > 0) {
            changeCountSpinnerModel.setValue(nChanges);
            changeTimesTableModel.setColumnCount(nChanges);
            changeTimesBox.setVisible(true);

            estimateTimesCheckBox.setEnabled(true);
            estimateTimesCheckBox.setSelected(skylineVector.changeTimesInput.get().isEstimatedInput.get());
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
            scalarRatesCheckBox.setSelected(true);
            valuesTableModel.setRowCount(1);

            for (int i=0; i<valuesParameter.getDimension(); i++)
                valuesTableModel.setValueAt(valuesParameter.getValue(i), 0, i);
        } else {
            scalarRatesCheckBox.setSelected(false);
            valuesTableModel.setRowCount(skylineVector.getNTypes());

            for (int interval=0; interval<valuesParameter.getDimension(); interval++) {
                for (int typeIdx = 0; typeIdx< skylineVector.getNTypes(); typeIdx++) {
                    valuesTableModel.setValueAt(valuesParameter.getValue(interval), typeIdx, interval);
                }
            }
        }

        estimateValuesCheckBox.setSelected(valuesParameter.isEstimatedInput.get());
    }

    void saveToModel() {

        if (modelSaveInProcess)
            return;

        modelSaveInProcess = true;

        PartitionContext partitionContext = doc.getContextFor(m_beastObject);
        String partitionID = ".t:" + partitionContext.tree;

        int nChanges = (int)changeCountSpinnerModel.getValue();

        // Update table models:

        valuesTableModel.setColumnCount(nChanges+1);
        if (scalarRatesCheckBox.isEnabled()) {
            valuesTableModel.setRowCount(1);
        } else {
            valuesTableModel.setRowCount(skylineVector.getNTypes());
        }
        for (int rowIdx=0; rowIdx<valuesTableModel.getRowCount(); rowIdx++) {
            for (int colIdx=0; colIdx<valuesTableModel.getColumnCount(); colIdx++) {
                if (valuesTableModel.getValueAt(rowIdx,colIdx) == null)
                    valuesTableModel.setValueAt(0.0, rowIdx, colIdx);
            }
        }

        changeTimesTableModel.setColumnCount(nChanges);
        for (int colIdx=0; colIdx<changeTimesTable.getColumnCount(); colIdx++) {
            if (changeTimesTableModel.getValueAt(0, colIdx) == null)
                changeTimesTableModel.setValueAt(0.0, 0, colIdx);
        }

        // Save values

        RealParameter rateValuesParam = skylineVector.rateValuesInput.get();
        rateValuesParam.setDimension(
                changeTimesTableModel.getRowCount()*changeTimesTableModel.getColumnCount());
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
