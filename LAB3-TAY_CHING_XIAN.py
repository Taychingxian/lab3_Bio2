"""
Protein-Protein Interaction Network Analyzer
Author: Tay Ching Xian
Description: A Streamlit application for analyzing and visualizing protein-protein 
             interactions from BioGRID and STRING databases.
"""

import streamlit as st
import requests
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import io

def get_biogrid_interactions(gene_name, access_key, organism=9606):
    try:
        biogrid_url = "https://webservice.thebiogrid.org/interactions"
        params = {
            "accessKey": access_key,
            "format": "json",
            "searchNames": True,
            "geneList": gene_name,
            "organism": organism,
            "searchbiogridids": True,
            "includeInteractors": True
        }
        response = requests.get(biogrid_url, params=params)
        if response.status_code == 200:
            network = response.json()
            return network
        else:
            st.error(f"Error: Status code {response.status_code}")
            return None
    except Exception as e:
        st.error(f"Error retrieving BioGRID data: {e}")
        return None


def get_string_interactions(protein_name, limit=20, species=9606):
    try:
        string_url = "https://string-db.org/api/json/network"
        params = {
            "identifiers": protein_name,
            "species": species,
            "limit": limit
        }
        response = requests.get(string_url, params=params)
        if response.status_code == 200:
            network = response.json()
            return network
        else:
            st.error(f"Error: Status code {response.status_code}")
            return None
    except Exception as e:
        st.error(f"Error retrieving STRING data: {e}")
        return None


def analyze_network(network_graph):
    num_nodes = network_graph.number_of_nodes()
    num_edges = network_graph.number_of_edges()
    degree_centrality = nx.degree_centrality(network_graph)
    top_5 = sorted(degree_centrality.items(), key=lambda x: -x[1])[:5]
    
    return {
        'num_nodes': num_nodes,
        'num_edges': num_edges,
        'degree_centrality': degree_centrality,
        'top_5': top_5
    }


def visualize_network(network_graph, top_proteins=None, seed=123):
    fig, ax = plt.subplots(figsize=(12, 10))
    
    slayout = nx.spring_layout(network_graph, seed=seed)
    
    nx.draw(network_graph, slayout, with_labels=True, node_size=800, 
            node_color='#0ea5e9', font_size=9, font_color='white',
            edge_color='#94a3b8', width=1.5, ax=ax)
    
    if top_proteins:
        nx.draw_networkx_nodes(network_graph, slayout, nodelist=top_proteins, 
                               node_size=1000, node_color='#f59e0b', ax=ax)
    
    ax.set_title('Protein-Protein Interaction Network', fontsize=16, fontweight='bold')
    ax.axis('off')
    
    return fig


def main():
    st.set_page_config(
        page_title="PPI Network Analyzer",
        page_icon="🧬",
        layout="wide",
        initial_sidebar_state="expanded"
    )
    
    st.markdown("""
        <style>
        [data-testid="stSidebar"] {
            background-color: #1e293b;
        }
        [data-testid="stSidebar"] label {
            color: white !important;
            font-weight: 500;
        }
        .stButton button {
            background-color: #0ea5e9;
            color: white;
            border: none;
            border-radius: 8px;
            padding: 12px 24px;
            font-size: 16px;
            font-weight: 600;
            width: 100%;
        }
        .stButton button:hover {
            background-color: #0284c7;
        }
        [data-testid="stMetricValue"] {
            font-size: 24px;
            color: #0ea5e9;
        }
        hr {
            margin: 20px 0;
        }
        </style>
    """, unsafe_allow_html=True)
    
    st.title("🧬 Protein-Protein Interaction Network Analyzer")
    st.markdown("Analyze and visualize protein interactions from BioGRID and STRING databases")
    st.markdown("---")
    
    st.sidebar.header("Data Source Selection")
    
    database = st.sidebar.radio(
        "Choose Database:",
        ["STRING DB", "BioGRID DB"],
        help="STRING: Comprehensive functional protein associations\nBioGRID: Experimental protein interactions"
    )
    
    st.sidebar.markdown("---")
    st.sidebar.header("Input Parameters")
    
    if database == "BioGRID DB":
        access_key = st.sidebar.text_input(
            "BioGRID Access Key:",
            type="password",
            help="Get your access key from https://webservice.thebiogrid.org"
        )
        
        gene_name = st.sidebar.text_input(
            "Gene Symbol:",
            value="TP53",
            help="Enter a gene symbol (e.g., TP53, BRCA1, MB)"
        ).strip().upper()
        
        organism = st.sidebar.selectbox(
            "Organism:",
            [("Human (9606)", 9606), ("Mouse (10090)", 10090), ("Yeast (559292)", 559292)],
            format_func=lambda x: x[0]
        )[1]
        
    else:
        protein_name = st.sidebar.text_input(
            "Protein Name:",
            value="TP53",
            help="Enter a protein name (e.g., TP53, BRCA1, p53)"
        ).strip()
        
        limit = st.sidebar.slider(
            "Max Interactions:",
            min_value=5,
            max_value=50,
            value=20,
            help="Maximum number of interaction partners to retrieve"
        )
        
        species = st.sidebar.selectbox(
            "Species:",
            [("Human (9606)", 9606), ("Mouse (10090)", 10090), ("Yeast (4932)", 4932)],
            format_func=lambda x: x[0]
        )[1]
    
    analyze_button = st.sidebar.button("🔬 Analyze Network", type="primary", use_container_width=True)
    
    st.sidebar.markdown("---")
    st.sidebar.caption("Developed by Tay Ching Xian | Bio2 Lab 3")
    
    if analyze_button:
        if database == "BioGRID DB":
            if not access_key:
                st.warning("⚠️ Please enter your BioGRID access key")
                st.info("📝 Register at https://webservice.thebiogrid.org to get your access key")
                return
            
            if not gene_name:
                st.warning("⚠️ Please enter a gene symbol")
                return
            
            with st.spinner(f"Fetching interactions for {gene_name} from BioGRID..."):
                network_data = get_biogrid_interactions(gene_name, access_key, organism)
                
                if network_data:
                    network_df = pd.DataFrame.from_dict(network_data, orient='index')
                    
                    if network_df.empty:
                        st.warning(f"No interactions found for {gene_name}")
                        return
                    
                    network_df['OFFICIAL_SYMBOL_A'] = network_df['OFFICIAL_SYMBOL_A'].str.upper()
                    network_df['OFFICIAL_SYMBOL_B'] = network_df['OFFICIAL_SYMBOL_B'].str.upper()
                    
                    ppi_data = network_df[['OFFICIAL_SYMBOL_A', 'OFFICIAL_SYMBOL_B']]
                    ppi_data.columns = ['protein1', 'protein2']
                    
                    source_name = "BioGRID"
                    query_name = gene_name
                    
        else:
            if not protein_name:
                st.warning("⚠️ Please enter a protein name")
                return
            
            with st.spinner(f"Fetching interactions for {protein_name} from STRING..."):
                network_data = get_string_interactions(protein_name, limit, species)
                
                if network_data:
                    network_df = pd.json_normalize(network_data)
                    
                    if network_df.empty:
                        st.warning(f"No interactions found for {protein_name}")
                        return
                    
                    ppi_data = network_df[['preferredName_A', 'preferredName_B', 'score']]
                    ppi_data.columns = ['protein1', 'protein2', 'score']
                    
                    source_name = "STRING DB"
                    query_name = protein_name
        
        if network_data:
            st.success(f"✅ Successfully retrieved {len(ppi_data)} interactions from {source_name}")
            
            with st.spinner("Building network graph..."):
                network_graph = nx.from_pandas_edgelist(ppi_data, "protein1", "protein2")
                analysis = analyze_network(network_graph)
            
            col1, col2 = st.columns([1, 1.5], gap="large")
            
            with col1:
                st.subheader("📊 Network Statistics")
                
                st.markdown("#### Network Size")
                metric_col1, metric_col2 = st.columns(2)
                with metric_col1:
                    st.metric("Nodes (Proteins)", analysis['num_nodes'])
                with metric_col2:
                    st.metric("Edges (Interactions)", analysis['num_edges'])
                
                st.markdown("---")
                
                st.markdown("#### Top 5 Hub Proteins")
                st.markdown("*Proteins with highest degree centrality*")
                
                for i, (protein, centrality) in enumerate(analysis['top_5'], 1):
                    st.markdown(f"**{i}. {protein}** - Centrality: {centrality:.3f}")
                
                st.markdown("---")
                
                st.markdown("#### Interaction Data")
                st.dataframe(ppi_data.head(10), use_container_width=True, height=300)
                
                st.info(f"**Query Protein:** {query_name}\n\n**Database:** {source_name}")
            
            with col2:
                st.subheader("🔬 Network Visualization")
                
                with st.spinner("Generating network visualization..."):
                    top_protein_names = [p[0] for p in analysis['top_5']]
                    fig = visualize_network(network_graph, top_protein_names)
                    
                    st.pyplot(fig)
                
                st.markdown("<br>", unsafe_allow_html=True)
                
                info_col1, info_col2 = st.columns(2)
                
                with info_col1:
                    st.markdown("""
                        <div style='background: linear-gradient(135deg, #0ea5e9 0%, #0284c7 100%); 
                                    color: white; padding: 15px; border-radius: 8px;'>
                            <strong>🔵 Blue Nodes</strong>
                            <p style='margin: 10px 0 0 0; font-size: 14px;'>
                                Standard interaction partners
                            </p>
                        </div>
                    """, unsafe_allow_html=True)
                
                with info_col2:
                    st.markdown("""
                        <div style='background: linear-gradient(135deg, #f59e0b 0%, #d97706 100%); 
                                    color: white; padding: 15px; border-radius: 8px;'>
                            <strong>🟠 Orange Nodes</strong>
                            <p style='margin: 10px 0 0 0; font-size: 14px;'>
                                Top 5 hub proteins (high centrality)
                            </p>
                        </div>
                    """, unsafe_allow_html=True)
    
    st.markdown("---")
    st.markdown("""
        <div style='text-align: center; color: #64748b;'>
            <p>Data from <a href='https://thebiogrid.org/' target='_blank'>BioGRID</a> 
            and <a href='https://string-db.org/' target='_blank'>STRING DB</a></p>
        </div>
    """, unsafe_allow_html=True)


if __name__ == "__main__":
    main()
