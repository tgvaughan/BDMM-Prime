package beast.app.bdmmprime.beauti;

import bdmmprime.parameterization.TypeSet;
import beast.core.parameter.RealParameter;

import javax.swing.table.AbstractTableModel;

/**
 * Table model used to represent SV parameter values.
 */
abstract class SkylineValuesTableModel extends AbstractTableModel {

    int nIntervals;
    boolean scalar;
    TypeSet typeSet;

    public SkylineValuesTableModel(TypeSet typeSet, boolean scalar, int nIntervals) {
        this.typeSet = typeSet;
        this.nIntervals = nIntervals;
        this.scalar = scalar;
    }

    abstract void setIntervalCount(int nIntervalsNew);

    abstract void setScalar(boolean scalar);

    /**
     * @return string appropriate for initializing RealParameter
     * corresponding to table values.
     */
    abstract String getParameterString();

    /**
     * @param realParameter parameter from which to import values.
     */
    abstract void loadFromParameter(RealParameter realParameter);
}
