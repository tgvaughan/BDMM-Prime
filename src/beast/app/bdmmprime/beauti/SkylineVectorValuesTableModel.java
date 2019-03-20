package beast.app.bdmmprime.beauti;

import bdmmprime.parameterization.TypeSet;
import beast.core.parameter.RealParameter;

import javax.swing.table.AbstractTableModel;

/**
 * Table model used to represent SV parameter values.
 */
class SkylineVectorValuesTableModel extends SkylineValuesTableModel {

    double[][] data;
    int nRows;

    public SkylineVectorValuesTableModel(TypeSet typeSet, boolean scalar, int nIntervals) {
        super(typeSet, scalar, nIntervals);

        nRows = scalar ? 1 : typeSet.getNTypes();
        this.data = new double[nIntervals][nRows];
    }

    @Override
    public int getRowCount() {
        return nRows;
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
        data[columnIndex-1][rowIndex] = Double.parseDouble((String)aValue);
        fireTableDataChanged();
    }

    @Override
    public Object getValueAt(int rowIndex, int columnIndex) {
        if (columnIndex>0)
            return data[columnIndex-1][rowIndex];
        else {
            if (scalar)
                return "ALL";
            else
                return typeSet.getTypeName(rowIndex);
        }
    }

    public void setIntervalCount(int nIntervalsNew) {

        double[][] oldData = data;
        data = new double[nIntervalsNew][nRows];

        for (int interval=0; interval<nIntervalsNew; interval++) {
            for (int row=0; row<nRows; row++) {
                if (interval< nIntervals)
                    data[interval][row] = oldData[interval][row];
                else {
                    if (interval > 0)
                        data[interval][row] = oldData[nIntervals-1][row];
                    else {
                        if (row > 0)
                            data[0][row] = oldData[0][row-1];
                        else
                            data[0][0] = 0.0;
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

            double[][] oldData = data;
            data = new double[nIntervals][1];

            for (int interval=0; interval<nIntervals; interval++) {
                data[interval][0] = oldData[interval][0];
            }

            this.nRows = 1;

        } else {
            double[][] oldData = data;
            data = new double[nIntervals][typeSet.getNTypes()];

            for (int interval=0; interval<nIntervals; interval++) {
                for (int row=0; row<typeSet.getNTypes(); row++) {
                    data[interval][row] = oldData[interval][0];
                }
            }

            this.nRows = typeSet.getNTypes();
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
                valuesBuilder.append(" ").append(data[interval][row]);
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
                data[interval][0] = realParameter.getValue(interval);
        } else {
            setScalar(false);

            int valueIdx = 0;
            for (int interval=0; interval<nIntervals; interval++) {
                for (int typeIdx = 0; typeIdx<typeSet.getNTypes(); typeIdx++) {
                    data[interval][typeIdx] = realParameter.getValue(valueIdx++);
                }
            }
        }

        fireTableDataChanged();
    }
}
