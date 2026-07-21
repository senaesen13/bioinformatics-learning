#!/usr/bin/env python3
"""
GEM Model Parser & Topology Network Generator
Extracts Metabolite-Gene network associations directly from SBML (.xml) or YAML (.yml) Genome-Scale Metabolic Models (GEMs).
"""

import os
import sys
import xml.etree.ElementTree as ET
import pandas as pd

def parse_human_gem_xml(xml_path, out_csv_path):
    print(f"Parsing Genome-Scale Metabolic Model: {xml_path}")
    tree = ET.parse(xml_path)
    root = tree.getroot()

    ns = {
        "sbml": "http://www.sbml.org/sbml/level3/version1/core",
        "fbc": "http://www.sbml.org/sbml/level3/version1/fbc/version2"
    }

    # 1. Map Gene Product ID -> Gene Symbol (from listOfGeneProducts)
    gene_id_to_symbol = {}
    for gp in root.findall(".//fbc:listOfGeneProducts/fbc:geneProduct", ns):
        gp_id = gp.get("{http://www.sbml.org/sbml/level3/version1/fbc/version2}id") or gp.get("id")
        gp_name = gp.get("{http://www.sbml.org/sbml/level3/version1/fbc/version2}label") or gp.get("name") or gp.get("id")
        if gp_id:
            gene_id_to_symbol[gp_id] = gp_name

    print(f"Loaded {len(gene_id_to_symbol)} gene products from GEM.")

    # 2. Map Species ID -> Species Name & Compartment
    species_info = {}
    for species in root.findall(".//sbml:listOfSpecies/sbml:species", ns):
        s_id = species.get("id")
        name = species.get("name")
        comp = species.get("compartment")
        species_info[s_id] = {"name": name, "compartment": comp}

    print(f"Loaded {len(species_info)} species/metabolites from GEM.")

    # 3. Associate Metabolites with Enzyme Genes via Reactions
    records = []
    seen = set()

    for rxn in root.findall(".//sbml:listOfReactions/sbml:reaction", ns):
        rxn_id = rxn.get("id")
        
        # Collect species
        rxn_species = []
        for s_ref in rxn.findall(".//sbml:listOfReactants/sbml:speciesReference", ns) + rxn.findall(".//sbml:listOfProducts/sbml:speciesReference", ns):
            rxn_species.append(s_ref.get("species"))
        
        # Collect geneProductRef IDs
        genes = []
        for g_ref in rxn.findall(".//fbc:geneProductRef", ns):
            g_id = g_ref.get("{http://www.sbml.org/sbml/level3/version1/fbc/version2}geneProduct") or g_ref.get("geneProduct")
            if g_id:
                g_sym = gene_id_to_symbol.get(g_id, g_id)
                genes.append(g_sym)
        
        genes = list(set(genes))
        
        for s_id in rxn_species:
            s_data = species_info.get(s_id, {"name": s_id, "compartment": ""})
            met_name = s_data["name"]
            comp = s_data["compartment"]
            
            for g in genes:
                key = (met_name, g)
                if key not in seen:
                    seen.add(key)
                    records.append({
                        "metabolite_id": s_id,
                        "metabolite_name": met_name,
                        "compartment": comp,
                        "gene_symbol": g,
                        "reaction_id": rxn_id
                    })

    df_out = pd.DataFrame(records)
    df_out.to_csv(out_csv_path, index=False)
    print(f"Extracted {len(df_out)} GEM metabolite-gene network edges across {df_out['metabolite_name'].nunique()} metabolites.")
    print(f"Saved GEM network mapping to: {out_csv_path}")
    return df_out

if __name__ == "__main__":
    base_dir = os.path.dirname(os.path.abspath(__file__))
    gem_file = os.path.join(base_dir, "data", "Human-GEM.xml")
    out_file = os.path.join(base_dir, "data", "human_gem_topology_network.csv")
    
    if not os.path.exists(gem_file):
        gem_file = "/Users/k2254978/Desktop/Work/ToolBoxProgs/Human_GEM_2.0/model/Human-GEM.xml"
        
    parse_human_gem_xml(gem_file, out_file)
