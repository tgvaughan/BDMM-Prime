package bdmm.parameterization;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.evolution.tree.TraitSet;

import java.util.*;

/**
 * Ordered set of type names.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class TypeSet extends BEASTObject {

    public Input<String> valueInput = new Input<>("value", "Comma-delmited list of types.");
    public Input<TraitSet> typeTraitSetInput = new Input<>("typeTraitSet", "Type trait set defining list of types.");

    public Input<String> unknownTypeIndicatorInput = new Input<>("UnknownTypeIdentifier",
            "String used to identify unknown types. (Default is '?'.)", "?");

    protected SortedSet<String> typeNameSet;
    protected String unknownTypeIdentifier;

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
        unknownTypeIdentifier = unknownTypeIndicatorInput.get().toLowerCase();

        if (valueInput.get() != null) {
            for (String typeName : valueInput.get().split(","))
                if (!typeName.isEmpty())
                    typeNameSet.add(typeName);
        }

        if (typeTraitSetInput.get() != null)
            addTypesFromTypeTraitSet(typeTraitSetInput.get());
    }

    /**
     * Incorporates all of the traits present in the given trait set into the type set.
     *
     * @param typeTraitSet
     */
    public void addTypesFromTypeTraitSet(TraitSet typeTraitSet) {
        for (String taxon : typeTraitSet.taxaInput.get().getTaxaNames()) {
            String typeName = typeTraitSet.getStringValue(taxon);

            if (!typeName.toLowerCase().equals(unknownTypeIdentifier))
                typeNameSet.add(typeTraitSet.getStringValue(taxon));
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
     * In the instance that the type name matches the unknown type identifier,
     * returns -1.
     *
     * @param typeName name of type
     * @return numerical index representing type.
     */
    public int getTypeIndex(String typeName) {
        if (typeName.toLowerCase().equals(unknownTypeIdentifier))
            return -1;

        if (typeNameSet.contains(typeName))
            return (new ArrayList<>(typeNameSet).indexOf(typeName));
        else
            throw new IllegalArgumentException("TypeSet does not contain type with name " + typeName);
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
