package beast.app.bdmmprime.beauti;

import bdmmprime.parameterization.SkylineVectorParameter;
import bdmmprime.parameterization.TypeSet;
import beast.app.beauti.BeautiDoc;
import beast.app.draw.InputEditor;
import beast.core.BEASTInterface;
import beast.core.Input;
import beast.core.parameter.RealParameter;

import javax.swing.*;
import javax.swing.border.EtchedBorder;
import javax.swing.table.*;
import java.awt.*;

public class SkylineVectorInputEditor extends InputEditor.Base {

    SkylineVectorParameter skylineVector;

    SpinnerModel changeCountSpinnerModel;
    JCheckBox scalarRatesCheckBox, timesAreAgesCheckBox;

    DefaultTableModel changeTimesTableModel;
    JTable changeTimesTable;
    Box changeTimesBox;

    ValuesTableModel valuesTableModel;
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

        // Add elements specific to change times

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
        boxHoriz.add(makeHorizontalFiller());

        boxVert.add(boxHoriz);

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

        boxVert.add(changeTimesBox);

        // Add elements specific to values

        boxHoriz = Box.createHorizontalBox();
        boxHoriz.add(new JLabel("Values:"));
        valuesTableModel = new ValuesTableModel(skylineVector.typeSetInput.get(), true, 1);
        valuesTable = new JTable(valuesTableModel);
        valuesTable.setDefaultRenderer(Object.class, new DefaultTableCellRenderer() {

            Color defaultBG = null;

            @Override
            public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
                JLabel l = (JLabel)super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, column);

                if (defaultBG == null)
                    defaultBG = l.getBackground();

                if (column==0) {
                    l.setBackground(Color.LIGHT_GRAY);
                }  else {
                    l.setBackground(defaultBG);
                }

                return l;
            }
        });
        valuesTable.setShowGrid(true);
        valuesTable.setGridColor(Color.GRAY);
        valuesTable.setCellSelectionEnabled(false);
        Box valuesTableBoxCol = Box.createVerticalBox();
        valuesTableBoxCol.add(valuesTable.getTableHeader());
        valuesTableBoxCol.add(valuesTable);
        boxHoriz.add(valuesTableBoxCol);
        boxHoriz.add(makeHorizontalFiller());

        boxVert.add(boxHoriz);

        boxHoriz = Box.createHorizontalBox();
        scalarRatesCheckBox = new JCheckBox("Scalar values");
        boxHoriz.add(scalarRatesCheckBox, Component.LEFT_ALIGNMENT);
        estimateValuesCheckBox = new JCheckBox("Estimate values");
        boxHoriz.add(estimateValuesCheckBox);
        boxHoriz.add(makeHorizontalFiller());

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

        skylineVector.typeSetInput.get().initAndValidate();

        int nTypes = skylineVector.typeSetInput.get().getNTypes();
        int nIntervals = skylineVector.getChangeCount() + 1;

        RealParameter valuesParam = skylineVector.skylineValuesInput.get();
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
        timesAreAgesCheckBox.setSelected(skylineVector.timesAreAgesInput.get());

        // Load values

        RealParameter valuesParameter = skylineVector.skylineValuesInput.get();

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

        int nTypes = skylineVector.getNTypes();
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

        RealParameter skylineValuesParam = skylineVector.skylineValuesInput.get();
        skylineValuesParam.setDimension(valuesTableModel.getRowCount()*(nChanges+1));
        skylineValuesParam.valuesInput.setValue(valuesTableModel.getParameterString(), skylineValuesParam);
        skylineValuesParam.isEstimatedInput.setValue(estimateValuesCheckBox.isSelected(), skylineValuesParam);
        skylineVector.setInputValue("skylineValues", skylineValuesParam);
        skylineValuesParam.initAndValidate();

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

    /**
     * Table model used to represent SV parameter values.
     */
    class ValuesTableModel extends AbstractTableModel {

        double[] data;
        int nIntervals;
        boolean scalar;
        TypeSet typeSet;

        public ValuesTableModel(TypeSet typeSet, boolean scalar, int nIntervals) {
            this.typeSet = typeSet;
            this.nIntervals = nIntervals;
            this.scalar = scalar;

            int nRows = scalar ? 1 : typeSet.getNTypes();

            this.data = new double[nRows*nIntervals];
        }

        @Override
        public int getRowCount() {
            return scalar ? 1 : typeSet.getNTypes();
        }

        @Override
        public int getColumnCount() {
            return nIntervals+1;
        }

        @Override
        public String getColumnName(int column) {
            return column > 0 ? "Epoch " + column : "Types";
        }

        @Override
        public void setValueAt(Object aValue, int rowIndex, int columnIndex) {
            data[(columnIndex-1)*getRowCount() + rowIndex] = Double.parseDouble((String)aValue);
            fireTableDataChanged();
        }

        @Override
        public Object getValueAt(int rowIndex, int columnIndex) {
            if (columnIndex>0)
                return data[(columnIndex-1)*getRowCount() + rowIndex];
            else {
                if (scalar)
                    return "ALL";
                else
                    return typeSet.getTypeName(rowIndex);
            }
        }

        public void setIntervalCount(int nIntervalsNew) {

            double [] oldData = data;
            int nRows = getRowCount();
            data = new double[nIntervalsNew*nRows];

            for (int interval=0; interval<nIntervalsNew; interval++) {
                for (int row=0; row<nRows; row++) {
                    if (interval< nIntervals)
                        data[interval*nRows + row] = oldData[interval*nRows + row];
                    else {
                        if (interval > 0)
                            data[interval*nRows + row] = oldData[(nIntervals-1)*nRows + row];
                        else {
                            if (row > 0)
                                data[row] = oldData[row-1];
                            else
                                data[0] = 0.0;
                        }

                    }
                }
            }

            nIntervals = nIntervalsNew;

            fireTableDataChanged();
            fireTableStructureChanged();
        }

        public void setScalar(boolean scalar) {

            if (this.scalar == scalar)
                return;

            if (scalar) {

                double [] oldData = data;
                data = new double[nIntervals];

                int nRowsOld = getRowCount();

                for (int interval=0; interval<nIntervals; interval++) {
                    data[interval] = oldData[nRowsOld*interval];
                }

            } else {
                double [] oldData = data;
                data = new double[nIntervals*typeSet.getNTypes()];

                for (int interval=0; interval<nIntervals; interval++) {
                    for (int row=0; row<typeSet.getNTypes(); row++) {
                        data[interval*typeSet.getNTypes()+row] = oldData[interval];
                    }
                }
            }

            this.scalar = scalar;

            fireTableDataChanged();
            fireTableStructureChanged();
        }

        /**
         * @return string appropriate for initializing RealParameter
         * corresponding to table values.
         */
        public String getParameterString() {

            int rowCount = getRowCount();

            StringBuilder valuesBuilder = new StringBuilder();
            for (int interval=0; interval<nIntervals; interval++) {
                for (int row=0; row<rowCount; row++) {
                    valuesBuilder.append(" ").append(data[rowCount*interval + row]);
                }
            }

            return valuesBuilder.toString();
        }

        @Override
        public boolean isCellEditable(int rowIndex, int columnIndex) {
            return columnIndex>0;
        }

        /**
         * @param realParameter parameter from which to import values.
         */
        public void loadFromParameter(RealParameter realParameter) {
            if (realParameter.getDimension()==nIntervals) {
                setScalar(typeSet.getNTypes()>1);

                for (int interval=0; interval<nIntervals; interval++)
                    data[interval] = realParameter.getValue(interval);
            } else {
                setScalar(false);

                int valueIdx = 0;
                for (int interval=0; interval<nIntervals; interval++) {
                    for (int typeIdx = 0; typeIdx<typeSet.getNTypes(); typeIdx++) {
                        data[interval*typeSet.getNTypes() + typeIdx] = realParameter.getValue(valueIdx++);
                    }
                }
            }

            fireTableDataChanged();
        }
    }

    /**
     * @return a horizontal filler object.
     */
    public static Box.Filler makeHorizontalFiller() {
        return new Box.Filler(new Dimension(1,1),
                new Dimension(1,1),
                new Dimension(Integer.MAX_VALUE,1));
    }

}
