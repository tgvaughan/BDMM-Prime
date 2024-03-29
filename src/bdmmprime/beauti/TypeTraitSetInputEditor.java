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
package bdmmprime.beauti;

import bdmmprime.distribution.BirthDeathMigrationDistribution;
import bdmmprime.parameterization.TypeSet;
import bdmmprime.util.InitializedTraitSet;
import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.TraitSet;
import beast.base.inference.parameter.RealParameter;
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.inputeditor.BeautiPanel;
import beastfx.app.inputeditor.GuessPatternDialog;
import beastfx.app.inputeditor.InputEditor;
import beastfx.app.util.FXUtils;
import javafx.scene.control.Button;
import javafx.scene.control.TableColumn;
import javafx.scene.control.TableView;
import javafx.scene.control.cell.PropertyValueFactory;
import javafx.scene.control.cell.TextFieldTableCell;
import javafx.scene.layout.*;
import javafx.scene.paint.Color;

import java.util.stream.Collectors;

/**
 * BEAUti input editor for type traits.
 *
 * @author Tim Vaughan (tgvaughan@gmail.com)
 */
public class TypeTraitSetInputEditor extends InputEditor.Base {

    TableView<TaxonEntry> typeTable;
    TraitSet traitSet;
    TaxonSet taxonSet;

    public static class TaxonEntry {
        String taxon;
        TraitSet traitSet;

        public TaxonEntry (String taxon, TraitSet traitSet) {
            this.taxon = taxon;
            this.traitSet = traitSet;
        }

        public String getTaxon() {
            return taxon;
        }

        public String getType() {
            return traitSet.getStringValue(taxon);
        }

        public void setType(String newType) {
            String newInitString =
                    traitSet.taxaInput.get().getTaxaNames().stream()
                    .map(n -> n + "=" + (n.equals(taxon) ?  newType : traitSet.getStringValue(n)))
                    .collect(Collectors.joining(","));

            traitSet.traitsInput.setValue(newInitString, traitSet);
            traitSet.initAndValidate();
        }
    }

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
    public void init(Input<?> input, BEASTInterface beastObject, int itemNr, ExpandOption bExpandOption, boolean bAddButtons) {

        m_bAddButtons = bAddButtons;
        m_input = input;
        m_beastObject =beastObject;
        this.itemNr = itemNr;
        pane = FXUtils.newHBox();

        traitSet = (TraitSet)input.get();
        taxonSet = traitSet.taxaInput.get();

        typeTable = new TableView<>();
        typeTable.setEditable(true);
        typeTable.setPrefWidth(600);

        TableColumn<TaxonEntry,String> taxonNameCol = new TableColumn<>("Sample Name");
        taxonNameCol.setCellValueFactory(new PropertyValueFactory<>("taxon"));
        typeTable.getColumns().add(taxonNameCol);

        TableColumn<TaxonEntry,String> typeCol = new TableColumn<>("Type");
        typeCol.setCellValueFactory(new PropertyValueFactory<>("type"));
        typeCol.setCellFactory(TextFieldTableCell.forTableColumn());
        typeCol.setOnEditCommit(e -> {
            e.getRowValue().setType(e.getNewValue());
            refreshPanel();
        });
        typeTable.getColumns().add(typeCol);

        for (String taxon : taxonSet.asStringList())
            typeTable.getItems().add(new TaxonEntry(taxon, traitSet));

        Button guessButton = new Button("Auto-configure");
        guessButton.setOnAction(e -> {
            GuessPatternDialog dlg = new GuessPatternDialog(null,
                ".*(\\d\\d\\d\\d).*");
            
            String traitString = "";
            switch(dlg.showDialog("Auto-configure types")) {
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
                        
                        if (!traitStringBuilder.isEmpty())
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

            refreshPanel();
        });

        Button clearButton = new Button("Clear");
        clearButton.setOnAction(e -> {
            StringBuilder traitStringBuilder = new StringBuilder();
            for (String taxonName : taxonSet.asStringList()) {
                if (!traitStringBuilder.isEmpty())
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

            refreshPanel();
        });

        addInputLabel();

        VBox boxVert = FXUtils.newVBox();
        boxVert.setBorder(new Border(new BorderStroke(Color.LIGHTGRAY,
                BorderStrokeStyle.SOLID, null, null)));

        HBox boxHoriz = FXUtils.newHBox();
        boxHoriz.getChildren().add(guessButton);
        boxHoriz.getChildren().add(clearButton);
        boxVert.getChildren().add(boxHoriz);
        boxVert.getChildren().add(typeTable);
        pane.getChildren().add(boxVert);

        getChildren().add(pane);
    }


    /**
     * Ugly hack to keep equilibrium type frequency parameter dimension up to date.
     */
    void updateFrequencies() {

        for (BEASTInterface beastInterface : traitSet.getOutputs()) {
            if (!(beastInterface instanceof BirthDeathMigrationDistribution bdmmDistr))
                continue;

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
                System.err.println("Error updating root/origin type frequencies.");
            }
        }
    }

    @Override
    public void refreshPanel() {
        typeTable.refresh();
        updateFrequencies();
        sync();
        super.refreshPanel();
    }
}
