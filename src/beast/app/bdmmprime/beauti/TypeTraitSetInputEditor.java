/*
 * Copyright (C) 2015 Tim Vaughan (tgvaughan@gmail.com)
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
package beast.app.bdmmprime.beauti;

import bdmmprime.distributions.BirthDeathMigrationDistribution;
import bdmmprime.parameterization.TypeSet;
import bdmmprime.util.InitializedTraitSet;
import beast.app.beauti.BeautiDoc;
import beast.app.beauti.GuessPatternDialog;
import beast.app.draw.InputEditor;
import beast.core.BEASTInterface;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.TraitSet;

import javax.swing.*;
import javax.swing.table.AbstractTableModel;
import java.awt.event.ActionEvent;
import java.util.stream.Collectors;

/**
 * BEAUti input editor for type traits.
 *
 * @author Tim Vaughan (tgvaughan@gmail.com)
 */
public class TypeTraitSetInputEditor extends InputEditor.Base {

    TypeTraitTableModel tableModel;
    TraitSet traitSet;
    TaxonSet taxonSet;

    public TypeTraitSetInputEditor(BeautiDoc doc) {
        super(doc);

    }

    @Override
    public Class<?> type() {

        // Hack to ensure that this input editor is definitely used for
        // BDMM-Prime analyses, as InitializedTraitSet is used in place of
        // TraitSet in the BDMM-Prime BEAUti template.
        return InitializedTraitSet.class;
    }

    @Override
    public void init(Input<?> input, BEASTInterface plugin, int itemNr, ExpandOption bExpandOption, boolean bAddButtons) {

        traitSet = (TraitSet)input.get();
        taxonSet = traitSet.taxaInput.get();
        tableModel = new TypeTraitTableModel(traitSet);
        JTable table = new JTable(tableModel);

        JButton guessButton = new JButton("Auto-configure");
        guessButton.addActionListener((ActionEvent e) -> {
            GuessPatternDialog dlg = new GuessPatternDialog(null,
                ".*(\\d\\d\\d\\d).*");
            
            String traitString = "";
            switch(dlg.showDialog("Auto-configure locations")) {
                case canceled:
                    return;

                case trait:
                    traitString = dlg.getTraitMap().entrySet().stream()
                            .map(entry -> entry.getKey() + "=" + entry.getValue())
                            .collect(Collectors.joining(","));
                    break;

                case pattern:
                    StringBuilder traitStringBuilder = new StringBuilder();
                    for (String taxonName : taxonSet.asStringList()) {
                        String matchString = dlg.match(taxonName);
                        if (matchString == null || matchString.isEmpty())
                            return;
                        
                        if (traitStringBuilder.length()>0)
                            traitStringBuilder.append(",");
                        
                        traitStringBuilder.append(taxonName)
                            .append("=")
                            .append(matchString);
                    }
                    traitString = traitStringBuilder.toString();
                    break;
            }

            traitSet.traitsInput.setValue(traitString, traitSet);

            try {
                traitSet.initAndValidate();
            } catch (Exception ex) {
                System.err.println("Error setting type trait.");
                ex.printStackTrace();
            }

            updateFrequencies();
            refreshPanel();
        });

        JButton clearButton = new JButton("Clear");
        clearButton.addActionListener((ActionEvent e) -> {
            StringBuilder traitStringBuilder = new StringBuilder();
            for (String taxonName : taxonSet.asStringList()) {
                if (traitStringBuilder.length()>0)
                    traitStringBuilder.append(",");
                traitStringBuilder.append(taxonName).append("=0");
            }
            traitSet.traitsInput.setValue(traitStringBuilder.toString(), traitSet);
            try {
                traitSet.initAndValidate();
            } catch (Exception ex) {
                System.err.println("Error clearing type trait.");
                ex.printStackTrace();
            }

            updateFrequencies();
            refreshPanel();
        });

        Box boxVert = Box.createVerticalBox();

        Box boxHoriz = Box.createHorizontalBox();
        boxHoriz.add(Box.createHorizontalGlue());
        boxHoriz.add(guessButton);
        boxHoriz.add(clearButton);
        boxVert.add(boxHoriz);
        boxVert.add(new JScrollPane(table));

        add(boxVert);
    }

    class TypeTraitTableModel extends AbstractTableModel {

        TraitSet typeTraitSet;

        public TypeTraitTableModel(TraitSet typeTraitSet) {
            this.typeTraitSet = typeTraitSet;
        }

        @Override
        public int getRowCount() {
            return typeTraitSet.taxaInput.get().getTaxonCount();
        }

        @Override
        public int getColumnCount() {
            return 2;
        }

        @Override
        public Object getValueAt(int rowIndex, int columnIndex) {
            if (columnIndex<0 || columnIndex>=getRowCount())
                return null;

            switch(columnIndex) {
                case 0:
                    // Taxon name:
                    return typeTraitSet.taxaInput.get().getTaxonId(rowIndex);
                case 1:
                    // Type:
                    return typeTraitSet.getStringValue(rowIndex);
                default:
                    return null;
            }
        }

        @Override
        public boolean isCellEditable(int rowIndex, int columnIndex) {
            return columnIndex == 1;
        }

        @Override
        public void setValueAt(Object aValue, int rowIndex, int columnIndex) {
            String taxon = taxonSet.getTaxonId(rowIndex);
            String traitString = traitSet.traitsInput.get();
            int startIdx = traitString.indexOf(taxon + "=");
            int endIdx = traitString.indexOf(",", startIdx);

            String newTraitString = traitString.substring(0, startIdx);
            newTraitString += taxon + "=" + (String)aValue;
            if (endIdx>=0)
                newTraitString += traitString.substring(endIdx);

            traitSet.traitsInput.setValue(newTraitString, traitSet);
            try {
                traitSet.initAndValidate();
            } catch (Exception ex) {
                System.err.println("Error setting type trait value.");
                ex.printStackTrace();
            }

            fireTableCellUpdated(rowIndex, columnIndex);

            updateFrequencies();
            refreshPanel();
        }

        @Override
        public String getColumnName(int column) {
            switch(column) {
                case 0:
                    return "Name";
                case 1:
                    return "Location";
                default:
                    return null;
            }
        }
    }

    /**
     * Ugly hack to keep equilibrium type frequency parameter dimension up to date.
     */
    void updateFrequencies() {

        for (BEASTInterface beastInterface : traitSet.getOutputs()) {
            if (!(beastInterface instanceof BirthDeathMigrationDistribution))
                continue;

            BirthDeathMigrationDistribution bdmmDistr =
                    (BirthDeathMigrationDistribution)beastInterface;

            TypeSet typeSet = bdmmDistr.parameterizationInput.get().typeSetInput.get();
            typeSet.initAndValidate();
            int nTypes = typeSet.getNTypes();

            RealParameter frequencies = bdmmDistr.frequenciesInput.get();

            if (frequencies.getDimension() == nTypes)
                continue;

            StringBuilder frequenciesBuilder = new StringBuilder();

            for (int typeIdx=0; typeIdx<nTypes; typeIdx++) {
                frequenciesBuilder.append(" ").append(1.0/nTypes);
            }

            frequencies.setDimension(nTypes);
            frequencies.valuesInput.setValue(frequenciesBuilder.toString(), frequencies);

            try {
                frequencies.initAndValidate();
                bdmmDistr.initAndValidate();
            } catch (Exception ex) {
                System.err.println("Error updating BDMM geo frequencies.");
            }
        }
    }
}
