package beast.app.bdmmprime.beauti;

import bdmmprime.parameterization.TimedParameter;
import bdmmprime.parameterization.TypeSet;
import beast.app.beauti.BeautiDoc;
import beast.app.draw.InputEditor;
import beast.core.BEASTInterface;
import beast.core.Input;
import beast.core.parameter.RealParameter;

import javax.swing.*;
import javax.swing.border.EtchedBorder;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.DefaultTableModel;
import java.awt.*;

public class TimedParameterInputEditor extends InputEditor.Base {

    TimedParameter timedParameter;

    SpinnerModel timeCountSpinnerModel;
    JCheckBox scalarValuesCheckBox, timesAreAgesCheckBox;

    DefaultTableModel timesTableModel;
    JTable timesTable;
    Box elementsBox;

    ValuesTableModel valuesTableModel;
    JTable valuesTable;

    JCheckBox estimateValuesCheckBox, estimateTimesCheckBox;

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

        Box boxVert, boxHoriz;

        // Add elements specific to change times

        boxVert = Box.createVerticalBox();
        boxVert.setBorder(new EtchedBorder());

        boxHoriz = Box.createHorizontalBox();
        JLabel changePointLabel = new JLabel("Number of elements:");
        JSpinner changePointsSpinner = new JSpinner();
        timeCountSpinnerModel =
                new SpinnerNumberModel(0, 0, Integer.MAX_VALUE,1);
        changePointsSpinner.setModel(timeCountSpinnerModel);
        boxHoriz.add(changePointLabel);
        boxHoriz.add(changePointsSpinner);
        boxHoriz.add(makeHorizontalFiller());

        boxVert.add(boxHoriz);

        elementsBox = Box.createVerticalBox();

        Box timesBox = Box.createVerticalBox();
        Box timesBoxRow = Box.createHorizontalBox();
        timesBoxRow.add(new JLabel("Element times:"));
        timesTableModel = new DefaultTableModel(new String[] {"Time 1"}, 1);
        timesTable = new JTable(timesTableModel);
        timesTable.setShowGrid(true);
        timesTable.setGridColor(Color.GRAY);
        timesTable.setCellSelectionEnabled(false);
        Box timesTableBoxCol = Box.createVerticalBox();
        timesTableBoxCol.add(timesTable.getTableHeader());
        timesTableBoxCol.add(timesTable);
        timesBoxRow.add(timesTableBoxCol);
        timesBoxRow.add(makeHorizontalFiller());
        timesBox.add(timesBoxRow);
        timesBoxRow = Box.createHorizontalBox();
        timesAreAgesCheckBox = new JCheckBox("Times specified as ages");
        timesBoxRow.add(timesAreAgesCheckBox);
        estimateTimesCheckBox = new JCheckBox("Estimate times");
        timesBoxRow.add(estimateTimesCheckBox);
        timesBoxRow.add(makeHorizontalFiller());
        timesBox.add(timesBoxRow);

        elementsBox.add(timesBox);

        // Add elements specific to values

        Box valuesBox = Box.createVerticalBox();
        boxHoriz = Box.createHorizontalBox();
        boxHoriz.add(new JLabel("Values:"));
        valuesTableModel = new ValuesTableModel(timedParameter.typeSetInput.get(), true, 1);
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

        valuesBox.add(boxHoriz);

        boxHoriz = Box.createHorizontalBox();
        scalarValuesCheckBox = new JCheckBox("Scalar values");
        boxHoriz.add(scalarValuesCheckBox, Component.LEFT_ALIGNMENT);
        estimateValuesCheckBox = new JCheckBox("Estimate values");
        boxHoriz.add(estimateValuesCheckBox);
        boxHoriz.add(makeHorizontalFiller());

        valuesBox.add(boxHoriz);
        elementsBox.add(valuesBox);

        boxVert.add(elementsBox);

        add(boxVert);

        loadFromModel();

        // Add event listeners:
        timeCountSpinnerModel.addChangeListener(e -> saveToModel());
        timesTableModel.addTableModelListener(e -> saveToModel());
        timesAreAgesCheckBox.addItemListener(e -> saveToModel());
        estimateTimesCheckBox.addItemListener(e -> saveToModel());

        valuesTableModel.addTableModelListener(e -> saveToModel());
        scalarValuesCheckBox.addItemListener(e -> saveToModel());
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

        timedParameter.typeSetInput.get().initAndValidate();

        int nTypes = timedParameter.typeSetInput.get().getNTypes();
        int nTimes = timedParameter.getTimeCount();

        if (nTimes==0)
            return;

        RealParameter valuesParam = timedParameter.valuesInput.get();
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
     * Populate GUI elements with values/dimensions from current BEASTObject model.
     * Called immediately after init(), and thus after every refreshPanel().
     */
    void loadFromModel() {

        ensureParamsConsistent();

        int nTimes = timedParameter.getTimeCount();
        int nTypes = timedParameter.getNTypes();

        timeCountSpinnerModel.setValue(nTimes);

        if (nTimes > 0) {

            // Load times:

            timesTableModel.setColumnCount(nTimes);

            String[] columnNames = new String[nTimes];
            for (int i = 0; i < nTimes; i++) {
                columnNames[i] = "Time " + (i + 1);
            }
            timesTableModel.setColumnIdentifiers(columnNames);

            estimateTimesCheckBox.setEnabled(true);
            estimateTimesCheckBox.setSelected(timedParameter.timesInput.get().isEstimatedInput.get());

            RealParameter changeTimesParameter = timedParameter.timesInput.get();
            for (int i = 0; i < changeTimesParameter.getDimension(); i++)
                timesTableModel.setValueAt(changeTimesParameter.getValue(i), 0, i);

            timesAreAgesCheckBox.setSelected(timedParameter.timesAreAgesInput.get());

            // Load values

            RealParameter valuesParameter = timedParameter.valuesInput.get();

            valuesTableModel.setTimeCount(nTimes);
            valuesTableModel.loadFromParameter(valuesParameter);

            if (valuesParameter.getDimension() == nTimes) {
                if (nTypes > 1) {
                    scalarValuesCheckBox.setSelected(true);
                    scalarValuesCheckBox.setEnabled(true);
                } else {
                    scalarValuesCheckBox.setSelected(false);
                    scalarValuesCheckBox.setEnabled(false);
                }
            } else {
                scalarValuesCheckBox.setSelected(false);
            }

            estimateValuesCheckBox.setSelected(valuesParameter.isEstimatedInput.get());

        } else {
            timesTableModel.setColumnCount(0);
            valuesTableModel.setTimeCount(0);

            estimateTimesCheckBox.setEnabled(false);

            elementsBox.setVisible(false);
        }
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

        RealParameter valuesParam = timedParameter.valuesInput.get();
        RealParameter timesParam = timedParameter.timesInput.get();
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

        timedParameter.initAndValidate();

        modelSaveInProcess = false;

        refreshPanel();
    }

    RealParameter getValuesParam() {
        RealParameter changeTimesParam = timedParameter.timesInput.get();
        if (changeTimesParam == null) {

            int idx = timedParameter.getID().indexOf("TP");
            String prefix = timedParameter.getID().substring(0, idx);
            String suffix = timedParameter.getID().substring(idx+2);
            String changeTimeId = prefix + "Prob" + suffix;

            changeTimesParam = (RealParameter) doc.pluginmap.get(changeTimeId);
            if (changeTimesParam == null) {
                changeTimesParam = new RealParameter("0.0");
                changeTimesParam.setID(changeTimeId);
            }
        }

        return changeTimesParam;
    }

    RealParameter getTimesParam() {
        RealParameter changeTimesParam = timedParameter.timesInput.get();
        if (changeTimesParam == null) {

            int idx = timedParameter.getID().indexOf("TP");
            String prefix = timedParameter.getID().substring(0, idx);
            String suffix = timedParameter.getID().substring(idx+2);
            String changeTimeId = prefix + "Times" + suffix;

            changeTimesParam = (RealParameter) doc.pluginmap.get(changeTimeId);
            if (changeTimesParam == null) {
                changeTimesParam = new RealParameter("0.0");
                changeTimesParam.setID(changeTimeId);
            }
        }

        return changeTimesParam;
    }

    /**
     * Table model used to represent SV parameter values.
     */
    class ValuesTableModel extends AbstractTableModel {

        double[] data;
        int nTimes;
        boolean scalar;
        TypeSet typeSet;

        public ValuesTableModel(TypeSet typeSet, boolean scalar, int nTimes) {
            this.typeSet = typeSet;
            this.nTimes = nTimes;
            this.scalar = scalar;

            int nRows = scalar ? 1 : typeSet.getNTypes();

            this.data = new double[nRows*nTimes];
        }

        @Override
        public int getRowCount() {
            return scalar ? 1 : typeSet.getNTypes();
        }

        @Override
        public int getColumnCount() {
            return nTimes +1;
        }

        @Override
        public String getColumnName(int column) {
            return column > 0 ? "Time " + column : "Types";
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

        public void setTimeCount(int nTimesNew) {

            double [] oldData = data;
            int nRows = getRowCount();
            data = new double[nTimesNew*nRows];

            for (int timeIdx=0; timeIdx<nTimesNew; timeIdx++) {
                for (int row=0; row<nRows; row++) {
                    if (timeIdx< nTimes)
                        data[timeIdx*nRows + row] = oldData[timeIdx*nRows + row];
                    else {
                        if (timeIdx > 0)
                            data[timeIdx*nRows + row] = oldData[(nTimes -1)*nRows + row];
                        else {
                            if (row > 0)
                                data[row] = oldData[row-1];
                            else
                                data[0] = 0.0;
                        }

                    }
                }
            }

            nTimes = nTimesNew;

            fireTableDataChanged();
            fireTableStructureChanged();
        }

        public void setScalar(boolean scalar) {

            if (this.scalar == scalar)
                return;

            if (scalar) {

                double [] oldData = data;
                data = new double[nTimes];

                int nRowsOld = getRowCount();

                for (int timeIdx = 0; timeIdx< nTimes; timeIdx++) {
                    data[timeIdx] = oldData[nRowsOld*timeIdx];
                }

            } else {
                double [] oldData = data;
                data = new double[nTimes *typeSet.getNTypes()];

                for (int timeIdx = 0; timeIdx< nTimes; timeIdx++) {
                    for (int row=0; row<typeSet.getNTypes(); row++) {
                        data[timeIdx*typeSet.getNTypes()+row] = oldData[timeIdx];
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
            for (int timeIdx = 0; timeIdx< nTimes; timeIdx++) {
                for (int row=0; row<rowCount; row++) {
                    valuesBuilder.append(" ").append(data[rowCount*timeIdx + row]);
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
            if (realParameter.getDimension()== nTimes) {
                setScalar(typeSet.getNTypes()>1);

                for (int timeIdx = 0; timeIdx< nTimes; timeIdx++)
                    data[timeIdx] = realParameter.getValue(timeIdx);
            } else {
                setScalar(false);

                int valueIdx = 0;
                for (int interval = 0; interval< nTimes; interval++) {
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
