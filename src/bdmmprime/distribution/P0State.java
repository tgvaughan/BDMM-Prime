/*
 * Copyright (C) 2016-2025 ETH Zurich
 * Copyright (C) 2013-2018 Denise Kühnert
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

package bdmmprime.distribution;

/**
 * Class containing the values of P0.
 */
public class P0State {

    int dimension;
    public double[] p0;

    public P0State(int nTypes) {
        dimension = nTypes;
        p0 = new double[nTypes];
    }

    public P0State(double[] p0) {
        dimension = p0.length;
        this.p0 = p0;
    }

    @Override
    public String toString() {
	    StringBuilder sb = new StringBuilder();

        for (int type=0; type<dimension; type++) {
            if (type>0)
                sb.append(" ");

            sb.append("p0[").append(type).append("]=").append(p0[type]);
        }

        return sb.toString();
    }
}
