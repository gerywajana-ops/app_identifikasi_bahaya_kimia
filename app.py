#!/usr/bin/env python3
"""
Chemical Hazard Identifier - Streamlit Application
Aplikasi identifikasi bahaya kimia menggunakan PubChem API dan GHS Classification
"""

import streamlit as st
import requests
import json
from PIL import Image
from io import BytesIO
import pandas as pd
from dataclasses import dataclass
from typing import List, Optional, Dict

# =============================================================================
# CONFIGURATION
# =============================================================================

st.set_page_config(
    page_title="Chemical Hazard Identifier",
    page_icon="⚗️",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for styling & Animasi
st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap');
    
    * { font-family: 'Inter', sans-serif; }
    
    .main-header {
        text-align: center;
        padding: 2rem 0;
        background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%);
        border-radius: 15px;
        margin-bottom: 2rem;
        color: white;
        transition: transform 0.3s ease;
    }
    
    .main-header:hover {
        transform: translateY(-2px);
    }
    
    .main-header h1 {
        font-size: 2.5rem;
        font-weight: 700;
        margin-bottom: 0.5rem;
        color: white !important;
    }
    
    .main-header p {
        font-size: 1.1rem;
        opacity: 0.9;
        margin-bottom: 0;
        color: #e0e0e0 !important;
    }
    
    /* Animasi Card dan Hover */
    .hazard-card, .info-card {
        transition: all 0.3s ease;
    }
    .hazard-card:hover, .info-card:hover {
        transform: translateY(-4px);
        box-shadow: 0 8px 16px rgba(0,0,0,0.15);
    }
    
    .hazard-card {
        padding: 1.5rem;
        border-radius: 12px;
        margin-bottom: 1rem;
        border-left: 5px solid;
        color: #1a1a1a !important;
    }
    
    .hazard-physical { background-color: #ffe0b2 !important; border-left-color: #ff9800; }
    .hazard-health { background-color: #f8bbd0 !important; border-left-color: #e91e63; }
    .hazard-environmental { background-color: #c8e6c9 !important; border-left-color: #4caf50; }
    
    .info-card {
        background-color: #262730;
        border-radius: 12px;
        padding: 1.5rem;
        margin-bottom: 1rem;
        border: 1px solid #464855;
    }
    
    .search-container {
        padding: 1.5rem;
        border-radius: 15px;
        margin-bottom: 2rem;
        transition: box-shadow 0.3s ease;
    }
    .search-container:focus-within {
        box-shadow: 0 4px 12px rgba(79, 70, 229, 0.2);
    }
    
    .hazard-statement {
        padding: 12px 16px;
        margin: 6px 0;
        border-radius: 8px;
        font-size: 0.95rem;
        color: #1a1a1a !important; 
        font-weight: 500;
        transition: transform 0.2s ease;
    }
    .hazard-statement:hover {
        transform: scale(1.01);
    }
    
    .precautionary-statement {
        padding: 8px 12px;
        margin: 4px 0;
        border-radius: 6px;
        background-color: #1e1e24;
        color: #ffffff;
        font-size: 0.85rem;
        border-left: 3px solid #3f51b5;
        transition: background-color 0.2s;
    }
    .precautionary-statement:hover {
        background-color: #2d2d35;
    }
    
    /* Animasi Tombol Streamlit Bawaan */
    .stButton > button {
        transition: all 0.3s ease !important;
    }
    .stButton > button:hover {
        transform: translateY(-2px) !important;
        box-shadow: 0 4px 10px rgba(0,0,0,0.15) !important;
    }
</style>
""", unsafe_allow_html=True)


# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass
class HazardInfo:
    """Data class untuk menyimpan informasi bahaya kimia"""
    hazard_class: str
    category: str
    statement: str
    pictogram_code: str
    pictogram_name: str
    severity: str  # high, medium, low

@dataclass
class ChemicalCompound:
    """Data class untuk menyimpan data senyawa kimia"""
    cid: int
    name: str
    iupac_name: str
    molecular_formula: str
    molecular_weight: float
    synonyms: List[str]
    description: str
    hazards: List[HazardInfo]
    physical_properties: Dict
    safety_info: Dict
    pictogram_urls: List[str]


# =============================================================================
# PUBCHEM API FUNCTIONS
# =============================================================================

def get_cid_by_name(compound_name: str) -> Optional[int]:
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{compound_name}/cids/JSON"
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            cids = data.get('IdentifierList', {}).get('CID', [])
            return cids[0] if cids else None
        return None
    except Exception as e:
        st.error(f"Error mencari CID: {e}")
        return None

def get_compound_properties(cid: int) -> Optional[Dict]:
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/MolecularFormula,MolecularWeight,IUPACName,Charge,HBondDonorCount,HBondAcceptorCount,RotatableBondCount,ExactMass,MonoisotopicMass,TPSA,XLogP,IsomericSMILES/JSON"
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            properties = data.get('PropertyTable', {}).get('Properties', [])
            return properties[0] if properties else None
        return None
    except Exception as e:
        return None

def get_compound_synonyms(cid: int) -> List[str]:
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/synonyms/JSON"
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            synonyms = data.get('InformationList', {}).get('Information', [{}])[0].get('Synonym', [])
            return synonyms[:15]
        return []
    except:
        return []

def get_ghs_hazards(cid: int) -> List[HazardInfo]:
    hazards = []
    seen_statements = set()
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON/?heading=GHS+Classification"
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            sections = data.get('Record', {}).get('Section', [])
            for section in sections:
                subsections = section.get('Section', [])
                for sub in subsections:
                    if sub.get('TOCHeading') == 'GHS Classification':
                        info = sub.get('Information', [])
                        for item in info:
                            value = item.get('Value', {})
                            if 'StringWithMarkup' in value:
                                for markup in value['StringWithMarkup']:
                                    string = markup.get('String', '')
                                    if 'Hazard Class' in string or 'Category' in string or string.startswith('H'):
                                        if string not in seen_statements:
                                            seen_statements.add(string)
                                            if string.startswith('H'):
                                                hazards.append(parse_hazard_code(string))
                                            else:
                                                hazards.append(parse_hazard_statement(string))
        
        if not hazards:
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON/?heading=Hazard+Statements"
            response = requests.get(url, timeout=10)
            if response.status_code == 200:
                data = response.json()
                sections = data.get('Record', {}).get('Section', [])
                for section in sections:
                    subsections = section.get('Section', [])
                    for sub in subsections:
                        info = sub.get('Information', [])
                        for item in info:
                            value = item.get('Value', {})
                            if 'StringWithMarkup' in value:
                                for markup in value['StringWithMarkup']:
                                    string = markup.get('String', '')
                                    if string.startswith('H') and string not in seen_statements:
                                        seen_statements.add(string)
                                        hazards.append(parse_hazard_code(string))
    except:
        pass
    return hazards

def parse_hazard_statement(statement: str) -> HazardInfo:
    statement_lower = statement.lower()
    
    if 'explosive' in statement_lower or 'explosion' in statement_lower:
        return HazardInfo('Fisika', 'Explosive', statement, 'GHS01', 'Bahaya Ledakan', 'high')
    elif 'flammable' in statement_lower or 'fire' in statement_lower or 'pyrophor' in statement_lower:
        return HazardInfo('Fisika', 'Flammable', statement, 'GHS02', 'Bahaya Api', 'high')
    elif 'oxidiz' in statement_lower:
        return HazardInfo('Fisika', 'Oxidizing', statement, 'GHS03', 'Mengoksidasi', 'high')
    elif 'compressed' in statement_lower or 'gas under press' in statement_lower:
        return HazardInfo('Fisika', 'Compressed Gas', statement, 'GHS04', 'Gas Bertekanan', 'medium')
    elif 'corrosive' in statement_lower or 'corrosivity' in statement_lower or 'eye damag' in statement_lower:
        return HazardInfo('Kesehatan', 'Corrosive', statement, 'GHS05', 'Korosif', 'high')
    elif 'toxic' in statement_lower or 'fatal' in statement_lower or 'poison' in statement_lower:
        return HazardInfo('Kesehatan', 'Acute Toxicity', statement, 'GHS06', 'Toksisitas Akut', 'high')
    elif 'carcinogen' in statement_lower or 'mutagen' in statement_lower or 'target organ' in statement_lower or 'health hazard' in statement_lower:
        return HazardInfo('Kesehatan', 'Health Hazard', statement, 'GHS08', 'Bahaya Kesehatan', 'medium')
    elif 'irritant' in statement_lower or 'harmful' in statement_lower or 'sensitiz' in statement_lower:
        return HazardInfo('Kesehatan', 'Irritant', statement, 'GHS07', 'Pengiritasi', 'medium')
    elif 'environment' in statement_lower or 'aquatic' in statement_lower:
        return HazardInfo('Lingkungan', 'Environmental Hazard', statement, 'GHS09', 'Bahaya Lingkungan', 'medium')
    else:
        return HazardInfo('Umum', 'General Hazard', statement, 'GHS07', 'Bahaya Umum', 'low')

def parse_hazard_code(code: str) -> HazardInfo:
    h_code = code.split(':')[0].strip()
    description = code.split(':')[1].strip() if ':' in code else code
    return parse_hazard_statement(f"{h_code}: {description}")

def get_pictogram_url(pictogram_code: str) -> str:
    pictogram_urls = {
        'GHS01': 'https://upload.wikimedia.org/wikipedia/commons/thumb/9/98/GHS-pictogram-explos.svg/240px-GHS-pictogram-explos.svg.png',
        'GHS02': 'https://upload.wikimedia.org/wikipedia/commons/thumb/6/6f/GHS-pictogram-flamme.svg/240px-GHS-pictogram-flamme.svg.png',
        'GHS03': 'https://upload.wikimedia.org/wikipedia/commons/thumb/e/e5/GHS-pictogram-rondflam.svg/240px-GHS-pictogram-rondflam.svg.png',
        'GHS04': 'https://upload.wikimedia.org/wikipedia/commons/thumb/d/d1/GHS-pictogram-bottle.svg/240px-GHS-pictogram-bottle.svg.png',
        'GHS05': 'https://upload.wikimedia.org/wikipedia/commons/thumb/3/3c/GHS-pictogram-acid.svg/240px-GHS-pictogram-acid.svg.png',
        'GHS06': 'https://upload.wikimedia.org/wikipedia/commons/thumb/9/98/GHS-pictogram-skull.svg/240px-GHS-pictogram-skull.svg.png',
        'GHS07': 'https://upload.wikimedia.org/wikipedia/commons/thumb/2/20/GHS-pictogram-exclam.svg/240px-GHS-pictogram-exclam.svg.png',
        'GHS08': 'https://upload.wikimedia.org/wikipedia/commons/thumb/2/21/GHS-pictogram-silhouette.svg/240px-GHS-pictogram-silhouette.svg.png',
        'GHS09': 'https://upload.wikimedia.org/wikipedia/commons/thumb/b/b3/GHS-pictogram-pollu.svg/240px-GHS-pictogram-pollu.svg.png',
    }
    return pictogram_urls.get(pictogram_code.upper(), '')

def get_compound_2d_structure(cid: int) -> Optional[str]:
    return f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/PNG?image_size=large"

def get_cas_number(cid: int) -> str:
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/xrefs/RegistryID/JSON"
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            registry_ids = data.get('InformationList', {}).get('Information', [{}])[0].get('RegistryID', [])
            for rid in registry_ids:
                if '-' in rid and len(rid.split('-')) == 3:
                    return rid
            return registry_ids[0] if registry_ids else 'N/A'
        return 'N/A'
    except:
        return 'N/A'

# =============================================================================
# UI COMPONENTS
# =============================================================================

def render_header():
    st.markdown("""
    <div class="main-header">
        <h1>⚗️ Chemical Hazard Identifier</h1>
        <p>Aplikasi Identifikasi Bahaya Kimia Berbasis GHS & PubChem Database</p>
    </div>
    """, unsafe_allow_html=True)

def render_quick_search():
    """Render tombol pencarian cepat dengan trigger langsung"""
    st.markdown("### 🔥 Pencarian Cepat")
    common_compounds = [
        ('Methanol', '⚗️'), ('Ethanol', '🍷'), ('Acetone', '💅'),
        ('Sulfuric acid', '🧪'), ('Hydrochloric acid', '🧪'), ('Sodium hydroxide', '🧂'),
        ('Hydrogen peroxide', '💧'), ('Benzene', '⛽'), ('Formaldehyde', '🏠'),
        ('Ammonia', '💨'), ('Toluene', '🎨'), ('Nitric acid', '🧪'),
    ]
    cols = st.columns(4)
    for i, (name, icon) in enumerate(common_compounds):
        with cols[i % 4]:
            if st.button(f"{icon} {name}", key=f"quick_{i}", use_container_width=True):
                # Langsung mengisi session state kolom input teks utama
                st.session_state['search_input'] = name
                st.session_state['trigger_search'] = True
                st.rerun()

def render_hazard_badge(severity: str) -> str:
    colors = {'high': ('🔴', '#f44336', 'Tinggi'), 'medium': ('🟡', '#ff9800', 'Sedang'), 'low': ('🟢', '#4caf50', 'Rendah')}
    emoji, color, label = colors.get(severity, ('⚪', '#9e9e9e', 'Tidak Diketahui'))
    return f'<span style="background-color: {color}; color: white; padding: 2px 10px; border-radius: 12px; font-size: 0.75rem; font-weight: 600;">{emoji} {label}</span>'

def render_compound_overview(compound: ChemicalCompound):
    st.markdown('<div class="result-container">', unsafe_allow_html=True)
    col1, col2 = st.columns([1, 2])
    with col1:
        structure_url = get_compound_2d_structure(compound.cid)
        try:
            response = requests.get(structure_url, timeout=10)
            if response.status_code == 200:
                st.image(Image.open(BytesIO(response.content)), caption=f"Struktur 2D: {compound.name}", use_container_width=True)
        except:
            st.info("Gambar struktur tidak tersedia")
    
    with col2:
        st.markdown(f"""
        <div class="info-card">
            <h2 style="color: #ffffff; margin-bottom: 0.5rem;">{compound.name}</h2>
            <p style="color: #e0e0e0; font-style: italic; font-weight: 500; margin-bottom: 0;">{compound.iupac_name}</p>
        </div>
        """, unsafe_allow_html=True)
        
        mw_display = f"{float(compound.molecular_weight):.3f} g/mol" if compound.molecular_weight else 'N/A'
        properties_data = {
            'Properti': ['CID PubChem', 'Rumus Molekul', 'Massa Molekul', 'Nomor CAS'],
            'Nilai': [str(compound.cid), f"<b>{compound.molecular_formula}</b>", mw_display, get_cas_number(compound.cid)]
        }
        st.markdown(pd.DataFrame(properties_data).to_html(escape=False, index=False), unsafe_allow_html=True)
    st.markdown('</div>', unsafe_allow_html=True)

def render_synonyms(synonyms: List[str]):
    """Render sinonim dengan warna kontras tajam (Putih di atas Indigo)"""
    if not synonyms: return
    st.markdown("### 🏷️ Sinonim")
    cols = st.columns(5)
    for i, syn in enumerate(synonyms):
        with cols[i % 5]:
            st.markdown(f"""
            <div style="
                background-color: #3f51b5; 
                color: #ffffff !important; 
                padding: 8px 12px; 
                border-radius: 20px; 
                font-size: 0.85rem; 
                font-weight: 600; 
                text-align: center; 
                margin: 4px 0;
                box-shadow: 0 2px 5px rgba(0,0,0,0.2);
                white-space: nowrap;
                overflow: hidden;
                text-overflow: ellipsis;
                transition: transform 0.2s ease;
            " title="{syn}" onmouseover="this.style.transform='scale(1.05)'" onmouseout="this.style.transform='scale(1)'">
                {syn}
            </div>
            """, unsafe_allow_html=True)

def render_pictograms(hazards: List[HazardInfo]):
    """Render piktogram GHS anti-gagal menggunakan pemindaian kata kunci"""
    st.markdown("### ⚠️ Pictogram Bahaya GHS")
    
    if not hazards:
        st.info("Tidak ada data piktogram GHS yang tersedia (Daftar bahaya kosong).")
        return
        
    detected_codes = set()
    for h in hazards:
        stmt_text = ""
        if h.statement: stmt_text += " " + h.statement.lower()
        if h.pictogram_code: stmt_text += " " + h.pictogram_code.lower()
        
        # Pemindaian Ekstraksi
        if 'ghs01' in stmt_text or 'explos' in stmt_text: detected_codes.add('GHS01')
        if 'ghs02' in stmt_text or 'flamm' in stmt_text or 'pyrophor' in stmt_text: detected_codes.add('GHS02')
        if 'ghs03' in stmt_text or 'oxidiz' in stmt_text: detected_codes.add('GHS03')
        if 'ghs04' in stmt_text or 'gas under press' in stmt_text: detected_codes.add('GHS04')
        if 'ghs05' in stmt_text or 'corros' in stmt_text or 'eye damag' in stmt_text: detected_codes.add('GHS05')
        if 'ghs06' in stmt_text or 'toxic' in stmt_text or 'fatal' in stmt_text or 'poison' in stmt_text: detected_codes.add('GHS06')
        if 'ghs07' in stmt_text or 'irritat' in stmt_text or 'harmful' in stmt_text: detected_codes.add('GHS07')
        if 'ghs08' in stmt_text or 'carcinogen' in stmt_text or 'mutagen' in stmt_text or 'target organ' in stmt_text: detected_codes.add('GHS08')
        if 'ghs09' in stmt_text or 'aquatic' in stmt_text or 'environment' in stmt_text: detected_codes.add('GHS09')

    if not detected_codes:
        st.info("Senyawa tergolong aman atau tidak memerlukan piktogram bahaya GHS khusus.")
        return
        
    ghs_names = {
        'GHS01': 'Explosive', 'GHS02': 'Flammable', 'GHS03': 'Oxidizing', 'GHS04': 'Compressed Gas',
        'GHS05': 'Corrosive', 'GHS06': 'Acute Toxicity', 'GHS07': 'Harmful/Irritant',
        'GHS08': 'Health Hazard', 'GHS09': 'Environmental'
    }
    
    cols = st.columns(min(len(detected_codes), 4))
    for i, code in enumerate(sorted(list(detected_codes))):
        url = get_pictogram_url(code)
        with cols[i % min(len(detected_codes), 4)]:
            if url:
                # Menggunakan st.image langsung dengan URL PNG agar aman dan tidak blank
                st.image(url, caption=ghs_names.get(code, code), use_container_width=True)


# =============================================================================
# MAIN LOGIC
# =============================================================================

def main():
    render_header()
    render_quick_search()
    
    st.markdown("### 🔍 Cari Senyawa Kimia")
    
    # Menyiapkan session state input agar terhubung dengan tombol quick search
    if 'search_input' not in st.session_state:
        st.session_state['search_input'] = ""
        
    search_query = st.text_input(
        "Nama Senyawa / Rumus / CID",
        key="search_input",
        placeholder="Contoh: methanol, sulfuric acid, NaOH..."
    )
    
    search_button = st.button("🔍 Identifikasi Bahaya", type="primary", use_container_width=True)
    
    # Deteksi apabila pencarian dipicu via tombol utama ATAU tombol pencarian cepat
    should_search = search_button or st.session_state.get('trigger_search', False)
    
    if should_search and search_query.strip():
        # Reset trigger
        if 'trigger_search' in st.session_state:
            st.session_state['trigger_search'] = False
            
        with st.spinner("🔬 Menganalisis senyawa dan mengidentifikasi bahaya..."):
            query = search_query.strip()
            cid = int(query) if query.isdigit() else get_cid_by_name(query)
            
            if cid:
                properties = get_compound_properties(cid)
                if properties:
                    synonyms = get_compound_synonyms(cid)
                    hazards = get_ghs_hazards(cid)
                    
                    compound = ChemicalCompound(
                        cid=cid,
                        name=properties.get('IUPACName', query.capitalize()),
                        iupac_name=properties.get('IUPACName', 'N/A'),
                        molecular_formula=properties.get('MolecularFormula', 'N/A'),
                        molecular_weight=properties.get('MolecularWeight', 0.0),
                        synonyms=synonyms,
                        description='',
                        hazards=hazards,
                        physical_properties=properties,
                        safety_info={},
                        pictogram_urls=[]
                    )
                    
                    st.success(f"✅ Senyawa ditemukan! CID: {cid}")
                    st.divider()
                    
                    render_compound_overview(compound)
                    if compound.synonyms:
                        render_synonyms(compound.synonyms)
                        
                    st.divider()
                    render_pictograms(compound.hazards)
                    
                    st.divider()
                    st.markdown("### 🏷️ Daftar Klasifikasi Bahaya Lengkap")
                    if hazards:
                        for h in hazards:
                            bg_color = '#ffcdd2' if h.severity == 'high' else '#ffe082' if h.severity == 'medium' else '#c8e6c9'
                            border_color = '#b71c1c' if h.severity == 'high' else '#e65100' if h.severity == 'medium' else '#1b5e20'
                            st.markdown(f"""
                            <div class="hazard-statement" style="background: {bg_color}; border-left: 5px solid {border_color}; color: #1a1a1a !important;">
                                {h.statement} {render_hazard_badge(h.severity)}
                            </div>
                            """, unsafe_allow_html=True)
                    else:
                        st.info("Senyawa tergolong aman atau rincian klasifikasi tidak terdaftar.")
                        
                else:
                    st.error(f"❌ Server PubChem tidak merespons untuk properti CID {cid}.")
            else:
                st.error(f"❌ Senyawa '{query}' tidak ditemukan di database.")
                
    elif should_search and not search_query.strip():
        st.warning("⚠️ Masukkan nama senyawa terlebih dahulu!")
        if 'trigger_search' in st.session_state:
            st.session_state['trigger_search'] = False

if __name__ == "__main__":
    main()
