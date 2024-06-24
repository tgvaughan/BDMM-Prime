/*
 * Copyright (C) 2019-2024 Tim Vaughan
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package bdmmprime.parameterization;

import bdmmprime.util.Utils;
import beast.base.inference.parameter.RealParameter;

import java.io.PrintStream;
import java.util.Arrays;

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
        super(changeTimesParam, skylineValuesParam, nTypes, null);
    }

    public SkylineVectorParameter(RealParameter changeTimesParam,
                                  RealParameter skylineValuesParam,
                                  int nTypes, RealParameter origin) {
        super(changeTimesParam, skylineValuesParam, nTypes, origin);
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
                    values[interval][i] = skylineValuesInput.get().getArrayValue(interval);
                else
                    values[interval][i] = skylineValuesInput.get().getArrayValue(interval * nTypes + i);
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
    public double[] getValuesAtTime(double time) {
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

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append(super.toString());

        sb.append(":");
        for (int i=0; i<getChangeCount()+1; i++) {
            if (i>0)
                sb.append(" (change time ").append(times[i-1]).append(")");
            sb.append(" ").append(Arrays.toString(values[i]));
        }

        return sb.toString();
    }
}
