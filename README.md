# lab3_Bio2
# Protein-Protein Interaction (PPI) Network Analyzer

## 🧬 Project Overview
This is a comprehensive web-based tool built with **Streamlit** that allows researchers and students to visualize and analyze Protein-Protein Interaction (PPI) networks.

The application dynamically fetches interaction data from two major biological databases—**BioGRID** and **STRING DB**—constructs a network graph, and calculates advanced topological metrics (centralities) to identify key "hub" proteins.

## ✨ Key Features

* **Dual Database Integration:**
    * **BioGRID:** Supports Access Key authentication and organism selection (Human, Mouse, Yeast).
    * **STRING DB:** Supports limit controls and species selection without requiring an API key.
* **Network Visualization:**
    * Generates 2D network graphs using `NetworkX` and `Matplotlib`.
    * Automatically highlights "Hub" proteins (Top 5) in **orange** and peripheral nodes in **blue**.
* **Topological Analysis:** Calculates five distinct centrality measures:
    * Degree Centrality
    * Betweenness Centrality
    * Closeness Centrality
    * Eigenvector Centrality
    * PageRank
* **Data Export:** View and explore the raw interaction dataframes.
