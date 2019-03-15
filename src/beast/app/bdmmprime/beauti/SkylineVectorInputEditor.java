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

    SkylineVectorParameter skylineVectorParameter;

    SpinnerModel changeCountSpinnerModel;
    JCheckBox scalarRatesCheckBox, timesAreAgesCheckBox;

    DefaultTableModel changeTimesTableModel;
    JTable changeTimesTable;
    Box changeTimesBox;

    DefaultTableModel valuesTableModel;
    JTable valuesTable;

    JCheckBox estimateValuesCheckBox, estimateTimesCheckBox;

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

        skylineVectorParameter = (SkylineVectorParameter)input.get();

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
    }

    void loadFromModel() {

        int nChanges = skylineVectorParameter.getChangeCount();

        // Load changepoints:
        if (nChanges > 0) {
            changeCountSpinnerModel.setValue(nChanges);
            changeTimesTableModel.setColumnCount(nChanges);
            changeTimesBox.setVisible(true);

            estimateTimesCheckBox.setEnabled(true);
            estimateTimesCheckBox.setSelected(skylineVectorParameter.changeTimesInput.get().isEstimatedInput.get());
        } else {
            changeCountSpinnerModel.setValue(0);
            changeTimesTableModel.setColumnCount(0);
            changeTimesBox.setVisible(false);

            estimateTimesCheckBox.setEnabled(false);
        }
        timesAreAgesCheckBox.setEnabled(skylineVectorParameter.timesAreAgesInput.get());

        // Load values

        valuesTableModel.setColumnCount(nChanges+1);

        RealParameter valuesParameter = skylineVectorParameter.rateValuesInput.get();
        if (valuesParameter.getDimension()==(nChanges+1)) {
            scalarRatesCheckBox.setSelected(true);
            valuesTableModel.setRowCount(1);

            for (int i=0; i<valuesParameter.getDimension(); i++)
                valuesTableModel.setValueAt(valuesParameter.getValue(i), 0, i);
        } else {
            scalarRatesCheckBox.setSelected(false);
            valuesTableModel.setRowCount(skylineVectorParameter.getNTypes());

            for (int interval=0; interval<valuesParameter.getDimension(); interval++) {
                for (int typeIdx=0; typeIdx<skylineVectorParameter.getNTypes(); typeIdx++) {
                    valuesTableModel.setValueAt(valuesParameter.getValue(interval), typeIdx, interval);
                }
            }
        }

        estimateValuesCheckBox.setSelected(valuesParameter.isEstimatedInput.get());
    }

    void saveToModel() {

    }
}
