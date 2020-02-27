package bdmmprime.parameterization;

import bdmmprime.util.Utils;
import beast.core.parameter.RealParameter;

import java.io.PrintStream;

public class SkylineVectorParameter extends SkylineParameter {

    double[][] values, storedValues;
    double[] valuesAtTime;

    boolean inputIsScalar;

    public SkylineVectorParameter() { }

    public SkylineVectorParameter(RealParameter changeTimesParam,
                                  RealParameter skylineValuesParam) {
        super(changeTimesParam, skylineValuesParam);
    }

    public SkylineVectorParameter(RealParameter changeTimesParam,
                                  RealParameter skylineValuesParam,
                                  int nTypes) {
        super(changeTimesParam, skylineValuesParam, nTypes);
    }


    @Override
    public void initAndValidate() {
        super.initAndValidate();

        if (skylineValuesInput.get().getDimension() % nIntervals != 0)
            throw new IllegalArgumentException("Value parameter dimension must " +
                    "be a multiple of the number of intervals.");

        int valsPerInterval = skylineValuesInput.get().getDimension()/nIntervals;
        inputIsScalar = valsPerInterval==1;

        if (typeSetInput.get() != null) {
            nTypes = typeSetInput.get().getNTypes();

            if (!inputIsScalar && nTypes != valsPerInterval)
                throw new IllegalArgumentException("SkylineVector has an incorrect " +
                        "number of elements.");
        } else {
            nTypes = valsPerInterval;
        }

        values = new double[nIntervals][nTypes];
        storedValues = new double[nIntervals][nTypes];

        valuesAtTime = new double[nTypes];
    }

    @Override
    protected void updateValues() {

        for (int interval=0; interval<nIntervals; interval++) {
            for (int i=0; i<nTypes; i++) {
                if (inputIsScalar)
                    values[interval][i] = skylineValuesInput.get().getValue(interval);
                else
                    values[interval][i] = skylineValuesInput.get().getValue(interval * nTypes + i);
            }
        }

        if (timesAreAges)
            Utils.reverseArray(values);
    }

    /**
     * Retrieve value of vector at a chosen time (not age).
     *
     * @param time when to evaluate the skyline parameter.
     * @return value of the vector at the chosen time.
     */
    protected double[] getValuesAtTime(double time) {
        update();

        int intervalIdx = getIntervalIdx(time);

        System.arraycopy(values[intervalIdx], 0, valuesAtTime, 0, nTypes);

        return valuesAtTime;
    }

    public int getNTypes() {
        return nTypes;
    }

    @Override
    protected void store() {
        super.store();

        for (int interval=0; interval<nIntervals; interval++)
            System.arraycopy(values[interval], 0, storedValues[interval], 0, nTypes);
    }

    @Override
    protected void restore() {
        super.restore();

        double [][] tmp;
        tmp = values;
        values = storedValues;
        storedValues = tmp;
    }

    /*
     * Loggable implementation
     */

    @Override
    public void init(PrintStream out) {

        for (int interval=0; interval<nIntervals; interval++) {

            if (interval < nIntervals-1)
                out.print(getID() + "i" + interval + "_endtime\t");

            if (inputIsScalar) {
                out.print(getID());

                if (nIntervals > 1)
                    out.print("i" + interval);

                out.print("\t");

            } else {
                for (int type = 0; type < nTypes; type++) {
                    out.print(getID());

                    if (nIntervals > 1)
                        out.print("i" + interval + "_");

                    if (typeSetInput.get() != null)
                        out.print(typeSetInput.get().getTypeName(type));
                    else
                        out.print("type" + type);

                    out.print("\t");
                }
            }
        }
    }

    @Override
    public void log(long sample, PrintStream out) {

        for (int interval=0; interval<nIntervals; interval++) {

            if (interval<nIntervals-1)
                out.print(times[interval] + "\t");

            if (inputIsScalar) {
                out.print(values[interval][0] + "\t");
            } else {
                for (int type = 0; type < nTypes; type++)
                    out.print(values[interval][type] + "\t");
            }
        }
    }

    @Override
    public void close(PrintStream out) {
    }
}
