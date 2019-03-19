package beast.app.bdmmprime.beauti;

import bdmmprime.parameterization.TypeSet;
import beast.core.parameter.RealParameter;

/**
 * Table model used to represent SV parameter values.
 */
class SkylineMatrixValuesTableModel extends SkylineValuesTableModel {

    double[][][] data;
    int nRows;

    public SkylineMatrixValuesTableModel(TypeSet typeSet, boolean scalar, int nIntervals) {
        super(typeSet, scalar, nIntervals);

        nRows = scalar ? 1 : typeSet.getNTypes();
        this.data = new double[nIntervals][nRows][nRows];
    }

    @Override
    public int getRowCount() {
        return nRows;
    }

    @Override
    public int getColumnCount() {
        return nIntervals*nRows + 1;
    }

    @Override
    public String getColumnName(int column) {
        if (column==0)
            return "Types";

        if ((column-1)%nRows == 0)
            return "Epoch " + ((column-1)/nRows + 1);

        else
            return "";
    }

    @Override
    public void setValueAt(Object aValue, int rowIndex, int columnIndex) {
        int interval = (columnIndex-1)/nRows;
        int type1idx = rowIndex;
        int type2idx = (columnIndex-1)%nRows;

        data[interval][type1idx][type2idx] = Double.parseDouble((String)aValue);
        fireTableDataChanged();
    }

    @Override
    public Object getValueAt(int rowIndex, int columnIndex) {
        if (columnIndex==0) {
            if (scalar)
                return "ALL";
            else
                return typeSet.getTypeName(rowIndex);
        }

        int interval = (columnIndex-1)/nRows;
        int type1idx = rowIndex;
        int type2idx = (columnIndex-1)%nRows;

        if (!scalar && type1idx == type2idx)
            return "-";

        return data[interval][type1idx][type2idx];
    }

    public void setIntervalCount(int nIntervalsNew) {

        if (nIntervalsNew == nIntervals)
            return;

        nIntervals = nIntervalsNew;
        data = new double[nIntervals][nRows][nRows];

        fireTableDataChanged();
        fireTableStructureChanged();
    }

    public void setScalar(boolean scalar) {

        if (this.scalar == scalar)
            return;

        this.scalar = scalar;

        nRows = scalar ? 1 : typeSet.getNTypes();

        data = new double[nIntervals][nRows][nRows];

        fireTableDataChanged();
        fireTableStructureChanged();
    }

    /**
     * @return string appropriate for initializing RealParameter
     * corresponding to table values.
     */
    public String getParameterString() {

        StringBuilder valuesBuilder = new StringBuilder();
        for (int interval=0; interval<nIntervals; interval++) {
            if (scalar) {
                valuesBuilder.append(" ").append(data[interval][0][0]);
            } else {
                for (int i = 0; i < nRows; i++) {
                    for (int j = 0; j < nRows; j++) {
                        if (i==j)
                            continue;

                        valuesBuilder.append(" ").append(data[interval][i][j]);
                    }
                }
            }
        }

        return valuesBuilder.toString();
    }

    @Override
    public boolean isCellEditable(int rowIndex, int columnIndex) {
        return columnIndex>0 && (scalar || (columnIndex-1)%nRows != rowIndex);
    }

    /**
     * @param realParameter parameter from which to import values.
     */
    public void loadFromParameter(RealParameter realParameter) {
        if (realParameter.getDimension()==nIntervals) {
            setScalar(true);

            for (int interval=0; interval<nIntervals; interval++)
                data[interval][0][0] = realParameter.getValue(interval);
        } else {
            setScalar(false);

            int valueIdx = 0;
            for (int interval=0; interval<nIntervals; interval++) {
                for (int i = 0; i<typeSet.getNTypes(); i++) {
                    for (int j = 0; j<typeSet.getNTypes(); j++) {
                        if (j==i)
                            continue;

                        data[interval][i][j] = realParameter.getValue(valueIdx++);
                    }
                }
            }
        }

        fireTableDataChanged();
    }
}
