package beast.app.bdmmprime.beauti;

import bdmmprime.parameterization.TypeSet;
import beast.core.parameter.RealParameter;

import javax.swing.table.AbstractTableModel;

/**
 * Table model used to represent SV parameter values.
 */
class SkylineVectorValuesTableModel extends SkylineValuesTableModel {

    double[] data;

    public SkylineVectorValuesTableModel(TypeSet typeSet, boolean scalar, int nIntervals) {
        super(typeSet, scalar, nIntervals);

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
