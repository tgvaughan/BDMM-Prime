package bdmmprime.parameterization;

import bdmmprime.util.Utils;
import beast.base.core.Description;
import beast.base.core.Function;

import java.io.PrintStream;
import java.util.Arrays;

@Description("Skyline parameter representing a matrix containing.  (Diagonals must be zero.)")
public class SkylineMatrixParameter extends SkylineParameter {

    double[][][] values, storedValues;
    double[][] valuesAtTime;

    boolean inputIsScalar;

    public SkylineMatrixParameter() { }

    public SkylineMatrixParameter(Function changeTimesParam,
                                  Function skylineValuesParam) {
        super(changeTimesParam, skylineValuesParam);
    }

    public SkylineMatrixParameter(Function changeTimesParam,
                                  Function skylineValuesParam,
                                  int nTypes) {
        super(changeTimesParam, skylineValuesParam, nTypes, null);
    }

    public SkylineMatrixParameter(Function changeTimesParam,
                                  Function skylineValuesParam,
                                  int nTypes, Function origin) {
        super(changeTimesParam, skylineValuesParam, nTypes, origin);
    }

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        int totalElementCount = skylineValuesInput.get() != null
                ? skylineValuesInput.get().getDimension()
                : 0;

        if (totalElementCount % nIntervals != 0)
            throw new IllegalArgumentException("Value parameter dimension must " +
                    "be a multiple of the number of intervals.");

        int elementsPerMatrix = totalElementCount/nIntervals;
        inputIsScalar = elementsPerMatrix == 1;

        if (typeSetInput.get() != null) {
            nTypes = typeSetInput.get().getNTypes();

            if (!inputIsScalar && elementsPerMatrix != nTypes*(nTypes-1)) {
                throw new IllegalArgumentException("SkylineMatrix parameter has " +
                        "an incorrect number of elements.");
            }
        } else {
            nTypes = (int) Math.round((1 + Math.sqrt(1 + 4 * elementsPerMatrix)) / 2);
        }

        values = new double[nIntervals][nTypes][nTypes];
        storedValues = new double[nIntervals][nTypes][nTypes];

        valuesAtTime = new double[nTypes][nTypes];
    }

    @Override
    protected void updateValues() {
        if (nTypes<=1)
            return;

        int idx=0;
        for (int interval=0; interval<nIntervals; interval++) {
            for (int i=0; i<nTypes; i++) {
                for (int j=0; j<nTypes; j++) {
                    if (i==j) {
                        values[interval][i][j] = 0.0;
                        continue;
                    }

                    if (inputIsScalar)
                        values[interval][i][j] = skylineValuesInput.get().getArrayValue(interval);
                    else
                        values[interval][i][j] = skylineValuesInput.get().getArrayValue(idx);

                    idx += 1;
                }
            }
        }

        if (timesAreAges)
            Utils.reverseArray(values);
    }

    /**
     * Retrieve value of matrix parameter at particular time (not age).
     *
     * @param time when to evaluate the parameter.
     * @return the matrix value at the chosen time.
     */
    public double[][] getValuesAtTime(double time) {
        update();

        int intervalIdx = getIntervalIdx(time);

        for (int i=0; i<nTypes; i++) {
            System.arraycopy(values[intervalIdx][i], 0,
                    valuesAtTime[i], 0, nTypes);
        }

        return valuesAtTime;
    }

    public int getNTypes() {
        return nTypes;
    }

    @Override
    protected void store() {
        super.store();

        for (int interval=0; interval<nIntervals; interval++) {
            for (int i=0; i<nTypes; i++) {
                System.arraycopy(values[interval][i], 0,
                        storedValues[interval][i], 0, nTypes);
            }
        }
    }

    @Override
    protected void restore() {
        super.restore();

        double[][][] tmp;
        tmp = values;
        values = storedValues;
        storedValues = tmp;
    }

    /*
     * Loggable implementation
     */

    @Override
    public void init(PrintStream out) {
        update();

        if (nTypes == 1)
            return;

        for (int interval=0; interval<nIntervals; interval++) {

            if (interval < nIntervals-1)
                out.print(getID() + "_i" + interval + "_endtime\t");

            if (inputIsScalar) {
                out.print(getID());

                if (nIntervals > 1)
                    out.print("i" + interval);

                out.print("\t");

            } else {
                for (int type1=0; type1<nTypes; type1++) {
                    for (int type2 = 0; type2 < nTypes; type2++) {
                        if (type2 == type1)
                            continue;

                        out.print(getID());

                        if (nIntervals > 1)
                            out.print("i" + interval + "_");

                        if (typeSetInput.get() != null)
                            out.print(typeSetInput.get().getTypeName(type1) + "_to_"
                                    + typeSetInput.get().getTypeName(type2));
                        else
                            out.print(type1 + "_to_" + type2);

                        out.print("\t");
                    }
                }
            }
        }
    }

    @Override
    public void log(long sample, PrintStream out) {

        if (nTypes == 1)
            return;

        for (int interval=0; interval<nIntervals; interval++) {

            if (interval < nIntervals-1)
                out.print(times[interval] + "\t");

            if (inputIsScalar) {
                out.print(values[interval][0][1] + "\t");
            } else {
                for (int type1 = 0; type1 < nTypes; type1++) {
                    for (int type2 = 0; type2 < nTypes; type2++) {
                        if (type2 == type1)
                            continue;

                        out.print(values[interval][type1][type2] + "\t");
                    }
                }
            }
        }
    }

    @Override
    public void close(PrintStream out) {
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append(super.toString());

        sb.append(":");
        for (int i=0; i<getChangeCount()+1; i++) {
            if (i>0)
                sb.append(" (change time ").append(times[i-1]).append(")");
            sb.append(" [");
            for (int t=0; t<getNTypes(); t++)
                sb.append(Arrays.toString(values[i][t]));
            sb.append("]");
        }

        return sb.toString();
    }
}
