/*
 * Copyright (C) 2012 Tim Vaughan
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package beast.evolution.tree;

import beast.core.Description;
import beast.core.Input;

/**
 * @author Tim Vaughan, edited by Denise Kuehnert for use in bdmm
 */
@Description("Basic plugin describing a simple Markovian migration model.")
public class MigrationModelForward extends SCMigrationModel {


    public Input<Boolean> rateMatrixIsForwardInput = new Input<Boolean>(
            "rateMatrixIsForward",
            "Optional boolean parameter specifying if rate matrix is specified in forward time (default: true)",
            true);

    public MigrationModelForward() { }

    @Override
    public void initAndValidate() {

        super.initAndValidate();
    }


    /**
     * Obtain offset into "rate matrix" and associated flag arrays.
     *
     * @param i
     * @param j
     * @return Offset (or -1 if i==j)
     */
    @Override
    protected int getArrayOffset(int i, int j) {

        if (i == j)
            throw new RuntimeException("Programmer error: requested migration "
                    + "rate array offset for diagonal element of "
                    + "migration rate matrix.");

        if (rateMatrixIsForwardInput.get()) {
            int temp = i;
            i = j;
            j = temp;
        }

        if (rateMatrixIsSquare) {
            return i * nTypes + j;
        } else {
            if (symmetricRateMatrix) {
                if (j < i)
                    return i * (i - 1) / 2 + j;
                else
                    return j * (j - 1) / 2 + i;
            } else {
                if (j > i)
                    j -= 1;
                return i * (nTypes - 1) + j;
            }
        }
    }

}
