#TAY CHING XIAN (A23CS0307)
import streamlit as st
import requests
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx

st.set_page_config(page_title="PPI Network Analyzer", page_icon="🧬", layout="wide")

st.markdown("""
    <style>
    .stApp {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
    }
    h1 {
        color: #2d3748;
        font-weight: 700;
    }
    .stButton>button {
        background-color: #0ea5e9;
        color: white;
        border-radius: 5px;
        padding: 0.5rem 2rem;
        font-weight: 600;
    }
    </style>
    """, unsafe_allow_html=True)

def retrieve_ppi_biogrid(target_protein):
    access_key = st.session_state.get('biogrid_key', '')
    organism_id = st.session_state.get('organism_id', '9606')
    
    url = "https://webservice.thebiogrid.org/interactions"
    params = {
        "accessKey": access_key,
        "format": "json",
        "searchNames": "true",
        "geneList": target_protein,
        "organism": organism_id,
        "includeInteractors": "true"
    }
    
    try:
        response = requests.get(url, params=params)
        if response.status_code == 200:
            data = response.json()
            if data:
                interactions = []
                for key, value in data.items():
                    interactions.append({
                        'Protein_A': value['OFFICIAL_SYMBOL_A'],
                        'Protein_B': value['OFFICIAL_SYMBOL_B'],
                        'Experimental_System': value['EXPERIMENTAL_SYSTEM']
                    })
                return pd.DataFrame(interactions)
    except Exception as e:
        st.error(f"Error: {str(e)}")
    
    return pd.DataFrame()

def retrieve_ppi_string(target_protein):
    species = st.session_state.get('species_id', '9606')
    limit = st.session_state.get('max_interactions', 20)
    
    url = "https://string-db.org/api/json/network"
    params = {
        "identifiers": target_protein,
        "species": species,
        "limit": limit
    }
    
    try:
        response = requests.get(url, params=params)
        if response.status_code == 200:
            data = response.json()
            if data:
                interactions = []
                for interaction in data:
                    interactions.append({
                        'Protein_A': interaction['preferredName_A'],
                        'Protein_B': interaction['preferredName_B'],
                        'Score': interaction['score']
                    })
                return pd.DataFrame(interactions)
    except Exception as e:
        st.error(f"Error: {str(e)}")
    
    return pd.DataFrame()

def generate_network(dataframe):
    G = nx.Graph()
    
    if 'Protein_A' in dataframe.columns and 'Protein_B' in dataframe.columns:
        for _, row in dataframe.iterrows():
            G.add_edge(row['Protein_A'], row['Protein_B'])
    
    return G

def get_centralities(network_graph):
    degree_cent = nx.degree_centrality(network_graph)
    betweenness_cent = nx.betweenness_centrality(network_graph)
    closeness_cent = nx.closeness_centrality(network_graph)
    eigenvector_cent = nx.eigenvector_centrality(network_graph, max_iter=1000)
    pagerank_cent = nx.pagerank(network_graph)
    
    return [degree_cent, betweenness_cent, closeness_cent, eigenvector_cent, pagerank_cent]

st.title("🧬 Protein-Protein Interaction Network Analyzer")
st.markdown("Analyze PPI networks from BioGRID and STRING databases")
st.markdown("---")

col_input1, col_input2 = st.columns([2, 1])

with col_input1:
    protein_id = st.text_input("Enter Protein ID", placeholder="e.g., TP53, BRCA1, MB")

with col_input2:
    database = st.selectbox("Select Database", ["BioGRID", "STRING"])

if database == "BioGRID":
    access_key = st.text_input("BioGRID Access Key", type="password", 
                               help="Get your free key at https://webservice.thebiogrid.org")
    organism = st.selectbox("Select Organism", 
                           [("Human", "9606"), ("Mouse", "10090"), ("Yeast", "559292")],
                           format_func=lambda x: x[0])
    st.session_state['biogrid_key'] = access_key
    st.session_state['organism_id'] = organism[1]
else:
    max_interactions = st.slider("Max Interactions", 5, 50, 20)
    species = st.selectbox("Select Species",
                          [("Human", "9606"), ("Mouse", "10090"), ("Yeast", "4932")],
                          format_func=lambda x: x[0])
    st.session_state['species_id'] = species[1]
    st.session_state['max_interactions'] = max_interactions

st.markdown("---")

if st.button("🔍 Analyze Network", use_container_width=True):
    if database == "BioGRID" and not access_key:
        st.error("❌ Please enter your BioGRID access key")
    elif not protein_id:
        st.error("❌ Please enter a protein ID")
    else:
        with st.spinner(f"Fetching data from {database}..."):
            if database == "BioGRID":
                df_ppi = retrieve_ppi_biogrid(protein_id)
            else:
                df_ppi = retrieve_ppi_string(protein_id)
            
            if not df_ppi.empty:
                st.success(f"✅ Found {len(df_ppi)} interactions!")
                
                network = generate_network(df_ppi)
                
                if network.number_of_nodes() > 0:
                    centralities_list = get_centralities(network)
                    
                    col1, col2 = st.columns(2)
                    
                    with col1:
                        st.subheader("PPI data information")
                        
                        st.markdown("**📊 PPI DataFrame**")
                        st.dataframe(df_ppi, use_container_width=True, height=300)
                        
                        st.markdown("**📈 Network Details**")
                        detail_col1, detail_col2 = st.columns(2)
                        with detail_col1:
                            st.metric("Number of Nodes", network.number_of_nodes())
                        with detail_col2:
                            st.metric("Number of Edges", network.number_of_edges())
                        
                        st.markdown("**🌐 Network Visualization**")
                        fig, ax = plt.subplots(figsize=(10, 8))
                        pos = nx.spring_layout(network, seed=42)
                        
                        degree_cent = centralities_list[0]
                        top_5_nodes = sorted(degree_cent.items(), key=lambda x: x[1], reverse=True)[:5]
                        hub_names = [node for node, _ in top_5_nodes]
                        node_colors = ['#f59e0b' if node in hub_names else '#0ea5e9' for node in network.nodes()]
                        
                        nx.draw(network, pos, 
                                node_color=node_colors,
                                node_size=600,
                                with_labels=True,
                                font_size=9,
                                font_weight='bold',
                                edge_color='#94a3b8',
                                width=2,
                                alpha=0.8,
                                ax=ax)
                        
                        ax.set_title("Protein-Protein Interaction Network", 
                                    fontsize=14, fontweight='bold', pad=20)
                        ax.axis('off')
                        plt.tight_layout()
                        st.pyplot(fig)
                        plt.close()
                    
                    with col2:
                        st.subheader("Centrality Measures")
                        
                        centrality_names = [
                            "Degree Centrality",
                            "Betweenness Centrality", 
                            "Closeness Centrality",
                            "Eigenvector Centrality",
                            "PageRank Centrality"
                        ]
                        
                        centrality_descriptions = [
                            "📌 Number of direct connections",
                            "🔗 Importance in connecting different parts",
                            "📍 Average distance to all other nodes",
                            "⭐ Influence based on connected neighbors",
                            "🎯 Importance based on network structure"
                        ]
                        
                        centrality_explanations = [
                            "Proteins with more connections have higher degree centrality.",
                            "Proteins that bridge different groups have higher betweenness.",
                            "Proteins close to all others have higher closeness centrality.",
                            "Proteins connected to important proteins have higher eigenvector centrality.",
                            "Proteins receiving connections from important proteins have higher PageRank."
                        ]
                        
                        for i, (name, desc, explain, cent_dict) in enumerate(zip(
                            centrality_names, centrality_descriptions, 
                            centrality_explanations, centralities_list)):
                            
                            with st.expander(f"{desc} {name}", expanded=(i==0)):
                                st.markdown(f"*{explain}*")
                                st.markdown("---")
                                
                                sorted_cent = sorted(cent_dict.items(), 
                                                   key=lambda x: x[1], reverse=True)[:10]
                                
                                cent_df = pd.DataFrame(sorted_cent, 
                                                      columns=['Protein', 'Centrality Score'])
                                cent_df['Rank'] = range(1, len(cent_df) + 1)
                                cent_df = cent_df[['Rank', 'Protein', 'Centrality Score']]
                                
                                st.dataframe(cent_df, use_container_width=True, 
                                           hide_index=True, height=250)
                                
                                fig_bar, ax_bar = plt.subplots(figsize=(8, 5))
                                proteins = [p for p, _ in sorted_cent]
                                scores = [s for _, s in sorted_cent]
                                
                                bars = ax_bar.barh(proteins, scores, color='#0ea5e9')
                                
                                for j, bar in enumerate(bars):
                                    if j < 5:
                                        bar.set_color('#f59e0b')
                                
                                ax_bar.set_xlabel('Centrality Score', fontweight='bold')
                                ax_bar.set_ylabel('Protein', fontweight='bold')
                                ax_bar.set_title(f'Top 10 Proteins by {name}', 
                                               fontweight='bold', pad=15)
                                ax_bar.invert_yaxis()
                                ax_bar.grid(axis='x', alpha=0.3)
                                plt.tight_layout()
                                st.pyplot(fig_bar)
                                plt.close()
                                
                                if i == 0:
                                    st.markdown("**🏆 Top 5 Hub Proteins:**")
                                    for rank, (prot, score) in enumerate(sorted_cent[:5], 1):
                                        st.markdown(f"{rank}. **{prot}** - Score: {score:.4f}")
                else:
                    st.warning("⚠️ Network is empty. No nodes to analyze.")
            else:
                st.error("❌ No interactions found. Please check your input or try a different protein/database.")

st.markdown("---")
st.markdown("""
    <div style='text-align: center; padding: 20px;'>
        <p style='color: #64748b; font-size: 14px;'>
            <strong>Developed by Tay Ching Xian | Bio2 Lab 3 - November 2025</strong><br>
            Data from <a href='https://thebiogrid.org/' target='_blank'>BioGRID</a> 
            and <a href='https://string-db.org/' target='_blank'>STRING DB</a>
        </p>
    </div>
""", unsafe_allow_html=True)
