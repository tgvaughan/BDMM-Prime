/*
 * Copyright (c) 2017-2026 ETH Zürich
 *
 * This file is part of bdmm-prime.
 *
 * bdmm-prime is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at
 * your option) any later version.
 *
 * bdmm-prime is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with bdmm-prime. If not, see <https://www.gnu.org/licenses/>.
 */

package bdmmprime.distribution.flow.utils;

import bdmmprime.distribution.flow.flowSystems.BaseFlow;

import java.util.Arrays;

/**
 * This class caches the likelihood of the subtree below each node of a tree, including the edge
 * leading to the node.
 * <p>
 * Every entry is stored together with a fingerprint of everything the likelihood of the node depends
 * on besides the likelihoods of its children. An entry can be reused as long as its fingerprint
 * still matches, the flow is still the same and none of the children had to be recomputed.
 * <p>
 * Every node has two entries, such that the entries of the stored state survive while the entries of
 * a proposed state are computed.
 */
public class SubtreeLikelihoodCache {

    /**
     * Everything the likelihood of a subtree depends on besides the flow and the likelihoods of the
     * children subtrees.
     *
     * @param timeEdgeStart the time of the parent of the node.
     * @param timeEdgeEnd   the time of the node.
     * @param child1Nr      the number of the first child or -1 if the node is a leaf.
     * @param child2Nr      the number of the second child or -1 if the node is a leaf.
     * @param nodeType      the type index of the sample at the node or -1 if the node is not sampled.
     */
    public record Fingerprint(double timeEdgeStart, double timeEdgeEnd, int child1Nr, int child2Nr, int nodeType) { }

    private final Fingerprint[][] fingerprints;
    private final double[][][] likelihoods;
    private final double[][] logScalingFactors;

    private final int[] currentEntries;
    private final int[] storedEntries;

    private final boolean[] currentEntryIsValid;
    private final boolean[] storedEntryIsValid;

    private BaseFlow currentFlow;
    private BaseFlow storedFlow;

    public SubtreeLikelihoodCache(int numNodes) {
        this.fingerprints = new Fingerprint[2][numNodes];
        this.likelihoods = new double[2][numNodes][];
        this.logScalingFactors = new double[2][numNodes];

        this.currentEntries = new int[numNodes];
        this.storedEntries = new int[numNodes];

        this.currentEntryIsValid = new boolean[numNodes];
        this.storedEntryIsValid = new boolean[numNodes];
    }

    /**
     * Invalidates all current entries if they were computed using a different flow.
     */
    public void invalidateIfFlowChanged(BaseFlow flow) {
        if (flow == this.currentFlow) return;

        this.invalidate();
        this.currentFlow = flow;
    }

    /**
     * Invalidates all current entries. The stored entries are kept.
     */
    public void invalidate() {
        Arrays.fill(this.currentEntryIsValid, false);
    }

    /**
     * Returns true if the current entry of the given node was computed for the given fingerprint.
     */
    public boolean isValid(int nodeNr, Fingerprint fingerprint) {
        return this.currentEntryIsValid[nodeNr]
                && fingerprint.equals(this.fingerprints[this.currentEntries[nodeNr]][nodeNr]);
    }

    /**
     * Stores the likelihood of the subtree of the given node as the current entry. The stored entry
     * of the node is not overwritten.
     *
     * @param likelihood       the scaled per-type likelihood of the subtree. The array is not copied.
     * @param logScalingFactor the log of the factor the likelihood was scaled with.
     */
    public void put(int nodeNr, Fingerprint fingerprint, double[] likelihood, double logScalingFactor) {
        int entry = 1 - this.storedEntries[nodeNr];

        this.fingerprints[entry][nodeNr] = fingerprint;
        this.likelihoods[entry][nodeNr] = likelihood;
        this.logScalingFactors[entry][nodeNr] = logScalingFactor;

        this.currentEntries[nodeNr] = entry;
        this.currentEntryIsValid[nodeNr] = true;
    }

    /**
     * Returns the scaled per-type likelihood of the subtree of the given node. The array must not be modified.
     */
    public double[] getLikelihood(int nodeNr) {
        return this.likelihoods[this.currentEntries[nodeNr]][nodeNr];
    }

    /**
     * Returns the log of the factor the likelihood of the subtree of the given node was scaled with.
     */
    public double getLogScalingFactor(int nodeNr) {
        return this.logScalingFactors[this.currentEntries[nodeNr]][nodeNr];
    }

    /** Caching **/

    public void store() {
        System.arraycopy(this.currentEntries, 0, this.storedEntries, 0, this.currentEntries.length);
        System.arraycopy(this.currentEntryIsValid, 0, this.storedEntryIsValid, 0, this.currentEntryIsValid.length);
        this.storedFlow = this.currentFlow;
    }

    public void restore() {
        System.arraycopy(this.storedEntries, 0, this.currentEntries, 0, this.storedEntries.length);
        System.arraycopy(this.storedEntryIsValid, 0, this.currentEntryIsValid, 0, this.storedEntryIsValid.length);
        this.currentFlow = this.storedFlow;
    }

    /**
     * Invalidates all current and stored entries.
     */
    public void clear() {
        Arrays.fill(this.currentEntryIsValid, false);
        Arrays.fill(this.storedEntryIsValid, false);
        this.currentFlow = null;
        this.storedFlow = null;
    }

}
