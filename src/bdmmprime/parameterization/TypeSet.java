/*
 * Copyright (C) 2019-2025 ETH Zurich
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

import beast.base.core.BEASTObject;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.tree.TraitSet;

import java.util.*;

/**
 * Ordered set of type names.
 *
 * @author ETH Zurich
 */
public class TypeSet extends BEASTObject {

    public Input<String> valueInput = new Input<>("value",
            "Comma-delmited list of types.");

    public Input<TraitSet> typeTraitSetInput = new Input<>("typeTraitSet",
            "Type trait set defining list of types.");

    public Input<Boolean> allowAmbiguousTypesInput = new Input<>("allowAmbiguousTypes",
            "Allow specification of ambiguous types. (Default true.)", true);

    public Input<String> unknownTypeIndicatorInput = new Input<>("unknownTypeIdentifier",
            "String used to identify completely unknown types. (Default is '?'.)", "?");

    public Input<String> ambiguousTypeDelimiterInput = new Input<>("ambiguousTypeDelimiter",
            "Regular expression used to delimit sequences of possible type names. " +
                    "(Default is '\\|'.)", "\\|");

    protected SortedSet<String> typeNameSet;

    protected boolean allowAmbiguousTypes;
    protected String unknownTypeIdentifier, ambiguousTypeDelimiter;

    public TypeSet() {}

    /**
     * Constructor to produce type set containing types provided as arguments.
     * Useful for testing.
     *
     * @param typeNames varargs array of type names to include
     */
    public TypeSet(String ... typeNames) {
        initAndValidate();

        typeNameSet.addAll(Arrays.asList(typeNames));
    }

    /**
     * Constructor to produce type set with n distinct types.
     * Useful for testing.
     *
     * @param n number of unique types:
     */
    public TypeSet(int n) {
        initAndValidate();

        for (int i=0; i<n; i++) {
            typeNameSet.add(String.valueOf(i));
        }
    }

    @Override
    public void initAndValidate() {
        typeNameSet = new TreeSet<>();

        allowAmbiguousTypes = allowAmbiguousTypesInput.get();
        unknownTypeIdentifier = unknownTypeIndicatorInput.get();
        ambiguousTypeDelimiter = ambiguousTypeDelimiterInput.get();

        if (valueInput.get() != null) {
            for (String typeName : valueInput.get().split(","))
                if (!typeName.isEmpty())
                    typeNameSet.add(typeName);
        }

        if (typeTraitSetInput.get() != null)
            addTypesFromTypeTraitSet(typeTraitSetInput.get());

        // Report type<->index map
        Log.info("TypeSet " + getID() + ":");
        for (int i=0; i<getNTypes(); i++)
            Log.info("\t" + getTypeName(i) + " (" + i + ")");
    }

    /**
     * Incorporates all of the traits present in the given trait set into the type set.
     *
     * @param typeTraitSet
     */
    public void addTypesFromTypeTraitSet(TraitSet typeTraitSet) {
        for (String taxon : typeTraitSet.taxaInput.get().getTaxaNames()) {
            String typeName = typeTraitSet.getStringValue(taxon);

            if (allowAmbiguousTypes) {
                if (typeName.equals(unknownTypeIdentifier))
                    continue;

                String[] typeNameSplit = typeName.split(ambiguousTypeDelimiterInput.get());
                if (typeNameSplit.length>1) {
                    Arrays.stream(typeNameSplit).forEach(s -> typeNameSet.add(s.trim()));
                    continue;
                }
            }

            typeNameSet.add(typeName);
        }
    }

    /**
     * @return the number of unique types defined in this type set
     */
    public int getNTypes() {
        return typeNameSet.size();
    }

    /**
     * Returns the numerical index corresponding to the given type name.
     * <p>
     * In the instance that the type name matches the unknown type identifier,
     * or type ambiguities are enabled and the input string includes the
     * ambiguous type list delimitation character, the return value will
     * be -bits, where bits is a strictly positive integer specifying which
     * of the possible types are to be included in the ambiguity set.
     *
     * @param typeName name of type
     * @return numerical index representing type, or a negative value specifying
     * an ambiguity set.
     */
    public int getTypeIndex(String typeName) {

        if (allowAmbiguousTypes) {
            if (typeName.equals(unknownTypeIdentifier))
                return -((1 << getNTypes()) - 1); // Any possible type

            String[] typeNameSplit = typeName.split(ambiguousTypeDelimiter);
            if (typeNameSplit.length>1) {
                // Calculate integer whose binary representation indicates
                // presence/absence of types in ambiguous set
                int bits = Arrays.stream(typeNameSplit)
                        .map(s -> getTypeIndex(s.trim()))
                        .reduce(0, (acc, typeIdx) -> acc + (1 << typeIdx));

                return -bits;
            }
        }

        if (typeNameSet.contains(typeName))
            return (new ArrayList<>(typeNameSet).indexOf(typeName));
        else
            throw new IllegalArgumentException("TypeSet does not contain type with name " + typeName);
    }

    /**
     * Check whether given typeIdx represents an ambiguous type.
     *
     * @param typeIdx type index to test
     * @return true if typeIdx represents an ambiguous type set
     */
    public boolean isAmbiguousTypeIndex(int typeIdx) {
        return typeIdx<0;
    }

    /**
     * Test to see whether the ambiguous type set represented by ambigIdx
     * excludes the type with the given type index.
     *
     * @param ambigIdx ambiguous type "index"
     * @param typeIdx type index to check for presence of
     * @return true if abigIdx excludes type represented by typeIdx
     */
    public boolean ambiguityExcludesType(int ambigIdx, int typeIdx) {
        if (ambigIdx<0)
            return ((-ambigIdx) & (1 << typeIdx)) == 0;
        else
            return ambigIdx != typeIdx;
    }

    /**
     * @param typeIdx numerical index representing type
     * @return name of type
     */
    public String getTypeName(int typeIdx) {
        if (typeIdx<0)
            return unknownTypeIdentifier;

        if (typeIdx<typeNameSet.size())
            return (new ArrayList<>(typeNameSet).get(typeIdx));
        else
            return "type_" + typeIdx;
    }

    /**
     * @param typeName name of type
     * @return true iff this TypeSet contains a type with name typeName.
     */
    public boolean containsTypeWithName(String typeName) {
        return typeNameSet.contains(typeName);
    }

    /**
     * @return list of type names ordered according to type index
     */
    public List<String> getTypesAsList() {
        return new ArrayList<>(typeNameSet);
    }


}
