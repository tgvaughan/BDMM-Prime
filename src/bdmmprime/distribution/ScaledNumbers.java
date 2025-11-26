/*
 * Copyright (C) 2016-2025 ETH Zurich
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

import beast.base.core.Description;

/**
 * Created by Jeremie Scire (jscire) on 24.06.16.
 */

@Description("Set of initial conditions for a system of ODEs, all increased " +
        "by the same factor to prevent underflowing (which appears with " +
        "numbers too small to be distinguished from zero, when using the " +
        "type 'double')")
public class ScaledNumbers {

	// scale factor(sampling)
	private int factor;

	// set(sampling) of initial conditions on which the scale factor(sampling) was/were applied
	private double[] equation;

	public ScaledNumbers(int f, double[] e){
		this.factor = f;
		this.equation = e;
	}

	public ScaledNumbers(){
		this.factor = 0;
		this.equation = new double[0];
	}

	public double[] getEquation(){
		return this.equation;
	}
	
	public void setEquation(double[] newEquation){
		this.equation = newEquation;
	}

	public int getScalingFactor(){
		return this.factor;
	}
	
	public void setScalingFactor(int newFactor){
		this.factor = newFactor; 
	}
	
	public void augmentFactor(int increase) {
		factor += increase;
	}
}
