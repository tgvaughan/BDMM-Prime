/*
 * Copyright (C) 2019-2024 ETH Zurich
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

package bdmmprime.trajectories;

/******************************************************************************
 * Author: Keith Schwarz (htiek@cs.stanford.edu)
 *
 * An implementation of the alias method implemented using Vose's algorithm.
 * The alias method allows for efficient sampling of random values from a
 * discrete probability distribution (i.e. rolling a loaded die) in O(1) time
 * each after O(n) preprocessing time.
 *
 * For a complete writeup on the alias method, including the intuition and
 * important proofs, please see the article "Darts, Dice, and Coins: Smpling
 * from a Discrete Distribution" at
 *
 *                 http://www.keithschwarz.com/darts-dice-coins/
 */
import beast.base.util.Randomizer;

import java.util.ArrayDeque;
import java.util.Deque;

public final class ReplacementSampler {

    // The probability and alias tables.
    private final int[] alias;
    private final double[] probability;

    /**
     * Constructs a new AliasMethod to sample from a discrete distribution and
     * hand back outcomes based on the probability distribution.
     * <p>
     * Given as input a list of probabilities corresponding to outcomes 0, 1,
     * ..., n - 1, along with the random number generator that should be used
     * as the underlying generator, this constructor creates the probability
     * and alias tables needed to efficiently sample from this distribution.
     *
     * @param probabilities The list of probabilities.
     */
    public ReplacementSampler(double[] probabilities) {

        // Begin by doing basic structural checks on the inputs.
        if (probabilities == null)
            throw new NullPointerException();

        if (probabilities.length == 0)
            throw new IllegalArgumentException("Probability vector must be nonempty.");

        // Allocate space for the probability and alias tables.
        probability = new double[probabilities.length];
        alias = new int[probabilities.length];


        // Compute the average probability and cache it for later use.
        final double average = 1.0 / probabilities.length;

        // Make a copy of the probabilities list, since we will be making
        // changes to it.
        double [] probsPrime = new double[probabilities.length];
        System.arraycopy(probabilities, 0, probsPrime, 0, probabilities.length);

        // Create two stacks to act as worklists as we populate the tables.
        Deque<Integer> small = new ArrayDeque<>();
        Deque<Integer> large = new ArrayDeque<>();

        // Populate the stacks with the input probabilities.
        for (int i = 0; i < probsPrime.length; ++i) {
            /* If the probability is below the average probability, then we add
             * it to the small list; otherwise we add it to the large list.
             */
            if (probsPrime[i] >= average)
                large.add(i);
            else
                small.add(i);
        }

        /* As a note: in the mathematical specification of the algorithm, we
         * will always exhaust the small list before the big list.  However,
         * due to floating point inaccuracies, this is not necessarily true.
         * Consequently, this inner loop (which tries to pair small and large
         * elements) will have to check that both lists aren't empty.
         */
        while (!small.isEmpty() && !large.isEmpty()) {
            // Get the index of the small and the large probabilities.
            int less = small.removeLast();
            int more = large.removeLast();

            // These probabilities have not yet been scaled up to be such that
            // 1/n is given weight 1.0.  We do this here instead.
            probability[less] = probsPrime[less] * probsPrime.length;
            alias[less] = more;

            // Decrease the probability of the larger one by the appropriate amount.
            probsPrime[more] = (probsPrime[more] + probsPrime[less]) - average;

            /* If the new probability is less than the average, add it into the
             * small list; otherwise add it to the large list.
             */
            if (probsPrime[more] >= 1.0 / probsPrime.length)
                large.add(more);
            else
                small.add(more);
        }

        /* At this point, everything is in one list, which means that the
         * remaining probabilities should all be 1/n.  Based on this, set them
         * appropriately.  Due to numerical issues, we can't be sure which
         * stack will hold the entries, so we empty both.
         */
        while (!small.isEmpty())
            probability[small.removeLast()] = 1.0;
        while (!large.isEmpty())
            probability[large.removeLast()] = 1.0;
    }

    /**
     * Samples a value from the underlying distribution.
     *
     * @return A random value sampled from the underlying distribution.
     */
    public int next() {
        // Generate a fair die roll to determine which column to inspect.
        int column = Randomizer.nextInt(probability.length);

        // Generate a biased coin toss to determine which option to pick.
        boolean coinToss = Randomizer.nextDouble() < probability[column];

        // Based on the outcome, return either the column or its alias.
        return coinToss? column : alias[column];
    }
}
