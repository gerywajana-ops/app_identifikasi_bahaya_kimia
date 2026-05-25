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
import base64
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

# Custom CSS for styling
st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght=300;400;500;600;700&display=swap');
    
    * { font-family: 'Inter', sans-serif; }
    
    .main-header {
        text-align: center;
        padding: 2rem 0;
        background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%);
        border-radius: 15px;
        margin-bottom: 2rem;
        color: white;
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
    }
    
    .hazard-statement {
        padding: 12px 16px;
        margin: 6px 0;
        border-radius: 8px;
        font-size: 0.95rem;
        color: #1a1a1a !important;
        font-weight: 500;
    }
    
    .precautionary-statement {
        padding: 8px 12px;
        margin: 4px 0;
        border-radius: 6px;
        background-color: #1e1e24;
        color: #ffffff;
        font-size: 0.85rem;
        border-left: 3px solid #3f51b5;
    }
</style>
""", unsafe_allow_html=True)

# =============================================================================
# DATA MODELS
# =============================================================================

@dataclass
class HazardInfo:
    pictogram_code: str
    pictogram_name: str
    statement: str
    severity: str         # 'high', 'medium', 'low'
    hazard_class: str     # 'Fisika', 'Kesehatan', 'Lingkungan'

@dataclass
class ChemicalCompound:
    cid: int
    name: str
    iupac_name: str
    molecular_formula: str
    molecular_weight: str
    synonyms: List[str]

# =============================================================================
# PUBCHEM API FUNCTIONS
# =============================================================================

def get_cid_by_name(compound_name: str) -> Optional[int]:
    """Mendapatkan CID (Compound ID) dari PubChem berdasarkan nama"""
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
    """Mendapatkan properti senyawa dari PubChem"""
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/MolecularFormula,MolecularWeight,IUPACName,Charge,HBondDonorCount,HBondAcceptorCount,RotatableBondCount,ExactMass,MonoisotopicMass,TPSA,XLogP,IsomericSMILES/JSON"
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            properties = data.get('PropertyTable', {}).get('Properties', [])
            return properties[0] if properties else None
        return None
    except Exception as e:
        st.error(f"Error mendapatkan properti: {e}")
        return None


def get_compound_synonyms(cid: int) -> List[str]:
    """Mendapatkan sinonim senyawa"""
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/synonyms/JSON"
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            synonyms = data.get('InformationList', {}).get('Information', [{}])[0].get('Synonym', [])
            return synonyms[:10]
        return []
    except:
        return []


def get_ghs_hazards(cid: int) -> List[HazardInfo]:
    """Mendapatkan data bahaya GHS dari PubChem secara rekursif menyeluruh"""
    hazards = []
    seen_codes = set()
    
    def parse_json_recursive(node):
        if isinstance(node, dict):
            for k, v in node.items():
                if k == 'String' and isinstance(v, str):
                    if v.startswith('H') and ':' in v:
                        parts = v.split(':', 1)
                        h_code = parts[0].strip()
                        if len(h_code) >= 4 and h_code[1:].isdigit():
                            if h_code not in seen_codes:
                                seen_codes.add(h_code)
                                hazards.append(parse_hazard_code(v))
                    elif any(kwd in v.lower() for kwd in ['flammable', 'toxic', 'corrosive', 'irritant', 'harmful', 'fatal']):
                        parsed = parse_hazard_statement(v)
                        fake_code = f"{parsed.hazard_class}_{parsed.statement[:20]}"
                        if fake_code not in seen_codes:
                            seen_codes.add(fake_code)
                            hazards.append(parsed)
                else:
                    parse_json_recursive(v)
        elif isinstance(node, list):
            for item in node:
                parse_json_recursive(item)

    try:
        # Gunakan timeout yang sedikit lebih longgar (15 detik) untuk antisipasi lag server
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON/?heading=Safety+and+Hazards"
        response = requests.get(url, timeout=15)
        if response.status_code == 200:
            parse_json_recursive(response.json())
            
        if not hazards:
            url_fallback = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON"
            response_fb = requests.get(url_fallback, timeout=15)
            if response_fb.status_code == 200:
                parse_json_recursive(response_fb.json())
    except requests.exceptions.Timeout:
        st.warning("⚠️ Koneksi ke PubChem terlalu lambat (Timeout). Silakan coba klik Identifikasi kembali.")
    except Exception as e:
        st.warning(f"⚠️ Sistem gagal membaca GHS karena masalah teknis: {e}")
        
    return hazards

def parse_hazard_code(hazard_string: str) -> HazardInfo:
    parts = hazard_string.split(':', 1)
    code = parts[0].strip()
    statement = parts[1].strip() if len(parts) > 1 else hazard_string
    
    code_num_str = ''.join(filter(str.isdigit(), code))
    code_num = int(code_num_str) if code_num_str else 0
    
    if 200 <= code_num <= 299:
        hazard_class = 'Fisika'
        severity = 'high' if code_num in [220, 222, 224, 225, 240, 241] else 'medium'
    elif 300 <= code_num <= 399:
        hazard_class = 'Kesehatan'
        severity = 'high' if code_num in [300, 310, 330, 314, 318, 340, 350, 360] else 'medium'
    elif 400 <= code_num <= 499:
        hazard_class = 'Lingkungan'
        severity = 'medium' if code_num in [400, 410] else 'low'
    else:
        hazard_class = 'Kesehatan'
        severity = 'medium'
        
    pic_code = 'GHS07'
    if hazard_class == 'Fisika': pic_code = 'GHS02'
    elif hazard_class == 'Lingkungan': pic_code = 'GHS09'
    elif severity == 'high' and code_num in [300, 310, 330]: pic_code = 'GHS06'
    elif severity == 'high' and code_num in [314, 318]: pic_code = 'GHS05'
    elif code_num in [340, 350, 360, 370, 372]: pic_code = 'GHS08'
    
    return HazardInfo(
        pictogram_code=pic_code,
        pictogram_name=f"Bahaya {hazard_class}",
        statement=f"{code}: {statement}",
        severity=severity,
        hazard_class=hazard_class
    )


def parse_hazard_statement(statement_string: str) -> HazardInfo:
    stmt_lower = statement_string.lower()
    if any(kwd in stmt_lower for kwd in ['flamm', 'explos', 'oxidiz', 'pyrophor', 'reactive', 'gas under press']):
        hazard_class = 'Fisika'
        pic_code = 'GHS02' if 'flamm' in stmt_lower else 'GHS01' if 'explos' in stmt_lower else 'GHS03'
    elif any(kwd in stmt_lower for kwd in ['aquatic', 'environment', 'pollut']):
        hazard_class = 'Lingkungan'
        pic_code = 'GHS09'
    else:
        hazard_class = 'Kesehatan'
        if any(kwd in stmt_lower for kwd in ['corros', 'skin burn', 'eye damag']): pic_code = 'GHS05'
        elif any(kwd in stmt_lower for kwd in ['fatal', 'acute tox', 'poison']): pic_code = 'GHS06'
        elif any(kwd in stmt_lower for kwd in ['carcinogen', 'mutagen']): pic_code = 'GHS08'
        else: pic_code = 'GHS07'

    severity = 'high' if any(kwd in stmt_lower for kwd in ['fatal', 'danger', 'severe']) else 'medium'
    return HazardInfo(
        pictogram_code=pic_code,
        pictogram_name=f"Bahaya {hazard_class}",
        statement=statement_string,
        severity=severity,
        hazard_class=hazard_class
    )


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


def get_precautionary_statements(cid: int) -> List[str]:
    precautions = []
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON/?heading=Precautionary+Statement+Codes"
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            sections = data.get('Record', {}).get('Section', [])
            for section in sections:
                for sub in section.get('Section', []):
                    for item in sub.get('Information', []):
                        value = item.get('Value', {})
                        if 'StringWithMarkup' in value:
                            for markup in value['StringWithMarkup']:
                                string = markup.get('String', '')
                                if string and string.startswith('P'):
                                    precautions.append(string)
    except:
        pass
    return precautions


def get_nfpa_diamond(cid: int) -> Dict:
    nfpa = {'health': '0', 'flammability': '0', 'reactivity': '0', 'special': ''}
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON/?heading=NFPA+704+Diamond"
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            sections = data.get('Record', {}).get('Section', [])
            for section in sections:
                for item in section.get('Information', []):
                    value = item.get('Value', {})
                    if 'StringWithMarkup' in value:
                        for markup in value['StringWithMarkup']:
                            string = markup.get('String', '')
                            if 'Health' in string: nfpa['health'] = string.split(':')[-1].strip()
                            elif 'Flammability' in string: nfpa['flammability'] = string.split(':')[-1].strip()
                            elif 'Stability' in string or 'Reactivity' in string: nfpa['reactivity'] = string.split(':')[-1].strip()
                            elif 'Special' in string: nfpa['special'] = string.split(':')[-1].strip()
    except:
        pass
    return nfpa


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
# UI RENDER COMPONENTS
# =============================================================================

def render_header():
    st.markdown("""
    <div class="main-header">
        <h1>⚗️ Chemical Hazard Identifier</h1>
        <p>Aplikasi Identifikasi Bahaya Kimia Berbasis GHS & PubChem Database</p>
    </div>
    """, unsafe_allow_html=True)


def render_search_section() -> tuple:
    st.markdown('<div class="search-container">', unsafe_allow_html=True)
    col1, col2 = st.columns([4, 1])
    with col1:
        search_query = st.text_input(
            "Masukkan Nama Senyawa / Rumus Molekul / CID PubChem",
            placeholder="Contoh: methanol, H2SO4, ethanol, sodium hydroxide, acetone...",
            label_visibility="collapsed"
        )
    with col2:
        search_button = st.button("🔍 Identifikasi", type="primary", use_container_width=True)
    st.markdown('</div>', unsafe_allow_html=True)
    return search_query, search_button


def render_sidebar():
    with st.sidebar:
        st.markdown("""
        <div style="text-align: center; padding: 1rem 0;">
            <h2 style="color: #4f46e5;">⚗️ ChemHazard ID</h2>
            <p style="font-size: 0.9rem; color: #888;">Aplikasi Identifikasi Bahaya Kimia</p>
        </div>
        """, unsafe_allow_html=True)
        st.divider()
        st.markdown("""
        ### 📖 Panduan Penggunaan
        **1. Masukkan Nama Senyawa**
        - Ketik nama senyawa kimia (bahasa Inggris)
        - Atau masukkan rumus molekul (contoh: H2SO4)
        
        **2. Klik "Identifikasi"**
        - Sistem akan menyisir database PubChem
        
        ### 📋 Contoh Senyawa
        - Methanol, Ethanol, Acetone
        - Sulfuric acid (H2SO4)
        - Sodium hydroxide (NaOH)
        """)
        st.divider()
        st.markdown("""
        <b>Legenda Tingkat Bahaya:</b><br>
        🔴 <b>Tinggi</b> - Bahaya serius<br>
        🟡 <b>Sedang</b> - Bahaya moderat<br>
        🟢 <b>Rendah</b> - Bahaya minimal
        """, unsafe_allow_html=True)


def render_hazard_badge(severity: str) -> str:
    colors = {'high': ('🔴', '#f44336', 'Tinggi'), 'medium': ('🟡', '#ff9800', 'Sedang'), 'low': ('🟢', '#4caf50', 'Rendah')}
    emoji, color, label = colors.get(severity, ('⚪', '#9e9e9e', 'Tidak Diketahui'))
    return f'<span style="background-color: {color}; color: white; padding: 2px 10px; border-radius: 12px; font-size: 0.75rem; font-weight: 600;">{emoji} {label}</span>'


def render_compound_overview(compound: ChemicalCompound):
    col1, col2 = st.columns([1, 2])
    with col1:
        structure_url = get_compound_2d_structure(compound.cid)
        try:
            response = requests.get(structure_url, timeout=10)
            if response.status_code == 200:
                st.image(Image.open(BytesIO(response.content)), caption=f"Struktur 2D: {compound.name}", use_container_width=True)
            else:
                st.info("Gambar struktur tidak tersedia")
        except:
            st.info("Gambar struktur tidak dapat dimuat")
            
    with col2:
        st.markdown(f"""
        <div class="info-card">
            <h2 style="color: #4f46e5; margin-bottom: 0.5rem;">{compound.name}</h2>
            <p style="color: #aaa; font-style: italic;">IUPAC: {compound.iupac_name}</p>
        </div>
        """, unsafe_allow_html=True)
        
        properties_data = {
            'Properti': ['CID PubChem', 'Rumus Molekul', 'Massa Molekul', 'Nomor CAS'],
            'Nilai': [str(compound.cid), f"<b>{compound.molecular_formula}</b>", f"{float(compound.molecular_weight):.3f} g/mol" if compound.molecular_weight else 'N/A', get_cas_number(compound.cid)]
        }
        st.markdown(pd.DataFrame(properties_data).to_html(escape=False, index=False), unsafe_allow_html=True)


def render_physical_properties(properties: Dict):
    if not properties: return
    st.markdown("### 📊 Properti Fisikokimia Tambahan")
    props = {
        'Donor Ikatan H': str(properties.get('HBondDonorCount', '0')),
        'Akseptor Ikatan H': str(properties.get('HBondAcceptorCount', '0')),
        'TPSA': f"{properties.get('TPSA', 'N/A')} Å²",
        'XLogP': str(properties.get('XLogP', 'N/A'))
    }
    cols = st.columns(4)
    for i, (key, value) in enumerate(props.items()):
        with cols[i]:
            st.markdown(f"""
            <div style="background: #262730; padding: 12px; border-radius: 8px; border: 1px solid #464855; text-align: center;">
                <small style="color: #aaa;">{key}</small><br><b style="color: white; font-size: 1.1rem;">{value}</b>
            </div>
            """, unsafe_allow_html=True)


def render_pictograms(hazards: List[HazardInfo]):
    """Render piktogram GHS berdasarkan hasil ekstraksi data hazards"""
    st.markdown("### ⚠️ Pictogram Bahaya GHS")
    
    if not hazards:
        st.info("Tidak ada data piktogram GHS yang tersedia (Daftar bahaya kosong).")
        return
        
    detected_codes = set()
    
    # Perulangan untuk memindai seluruh data bahaya yang ada
    for h in hazards:
        stmt_text = ""
        if h.statement: stmt_text += " " + h.statement.lower()
        if h.pictogram_code: stmt_text += " " + h.pictogram_code.lower()
        if h.pictogram_name: stmt_text += " " + h.pictogram_name.lower()
        
        # Deteksi Kode Langsung internal objek
        if 'ghs01' in stmt_text: detected_codes.add('GHS01')
        if 'ghs02' in stmt_text: detected_codes.add('GHS02')
        if 'ghs03' in stmt_text: detected_codes.add('GHS03')
        if 'ghs04' in stmt_text: detected_codes.add('GHS04')
        if 'ghs05' in stmt_text: detected_codes.add('GHS05')
        if 'GHS06' in stmt_text: detected_codes.add('GHS06')
        if 'ghs07' in stmt_text: detected_codes.add('GHS07')
        if 'ghs08' in stmt_text: detected_codes.add('GHS08')
        if 'ghs09' in stmt_text: detected_codes.add('GHS09')
        
        # Scan kata kunci teks (Fallback mandiri jika kode tidak tersemat)
        if 'flamm' in stmt_text or 'pyrophor' in stmt_text: detected_codes.add('GHS02')
        if 'toxic' in stmt_text or 'fatal' in stmt_text or 'poison' in stmt_text: detected_codes.add('GHS06')
        if 'corros' in stmt_text or 'eye damag' in stmt_text or 'skin burn' in stmt_text: detected_codes.add('GHS05')
        if 'explos' in stmt_text: detected_codes.add('GHS01')
        if 'oxidiz' in stmt_text: detected_codes.add('GHS03')
        if 'gas under press' in stmt_text or 'compressed gas' in stmt_text: detected_codes.add('GHS04')
        if 'irritat' in stmt_text or 'harmful' in stmt_text or 'sensitiz' in stmt_text: detected_codes.add('GHS07')
        if 'carcinogen' in stmt_text or 'mutagen' in stmt_text or 'target organ' in stmt_text: detected_codes.add('GHS08')
        if 'aquatic' in stmt_text or 'environment' in stmt_text: detected_codes.add('GHS09')

    if not detected_codes:
        st.info("Senyawa tergolong aman atau tidak memerlukan piktogram bahaya GHS khusus.")
        return
        
    ghs_names = {
        'GHS01': 'Explosive (Mudah Meledak)', 'GHS02': 'Flammable (Mudah Terbakar)',
        'GHS03': 'Oxidizing (Pengoksidasi)', 'GHS04': 'Gases Under Pressure (Gas Bertekanan)',
        'GHS05': 'Corrosive (Korosif / Merusak)', 'GHS06': 'Acute Toxicity (Beracun)',
        'GHS07': 'Harmful / Irritant (Iritasi)', 'GHS08': 'Health Hazard (Bahaya Kronis)',
        'GHS09': 'Environmental Hazard (Bahaya Lingkungan)'
    }
    
    cols = st.columns(min(len(detected_codes), 4))
    for i, code in enumerate(sorted(list(detected_codes))):
        url = get_pictogram_url(code)
        with cols[i % min(len(detected_codes), 4)]:
            if url:
                st.image(url, caption=ghs_names.get(code, code), use_container_width=True)


def get_safety_recommendations(hazards: List[HazardInfo]) -> Dict:
    rec = {'ppe': set(), 'handling': set(), 'storage': set(), 'emergency': set(), 'disposal': set()}
    for h in hazards:
        code = h.pictogram_code
        if code in ['GHS05', 'GHS06']:
            rec['ppe'].update(['Sarung tangan kimia (Nitril/Neoprena)', 'Kacamata goggle / Face shield', 'Jas lab lengan panjang'])
            rec['handling'].update(['Wajib gunakan di dalam lemari asam (Fume Hood)', 'Hindari hirup uap gas langsung'])
        if code == 'GHS02':
            rec['handling'].update(['Jauhkan dari percikan api, panas, dan permukaan panas', 'Dilarang merokok di sekitar area'])
            rec['storage'].update(['Simpan di lemari khusus bahan mudah terbakar (Flammable Cabinet)', 'Tempatkan di ruangan sejuk'])
        if code == 'GHS09':
            rec['disposal'].update(['Buang sebagai Limbah B3 resmi', 'Jangan buang langsung ke wastafel / saluran air'])
    
    # Standar fallback perlindungan umum
    if not rec['ppe']: rec['ppe'].update(['Kacamata keselamatan standar', 'Sarung tangan pelindung', 'Jas laboratorium'])
    if not rec['handling']: rec['handling'].update(['Tangani dengan prinsip kehati-hatian laboratorium standar'])
    if not rec['storage']: rec['storage'].update(['Simpan di wadah tertutup rapat pada area kering'])
    if not rec['emergency']: rec['emergency'].update(['Bilas dengan air mengalir jika terpapar kontak langsung', 'Gunakan APAR yang sesuai jika terjadi api'])
    if not rec['disposal']: rec['disposal'].update(['Kumpulkan ke dalam wadah jeriken limbah kimia terpisah'])
    
    return {k: list(v) for k, v in rec.items()}


def render_footer():
    st.markdown("<p style='text-align: center; color: #666; font-size: 0.8rem;'>ChemHazard ID v2.5 © 2026 | Powered by PubChem API</p>", unsafe_allow_html=True)

# =============================================================================
# MAIN APPLICATION LOGIC
# =============================================================================

def main():
    render_header()
    render_sidebar()
    
    search_query, search_button = render_search_section()
    
    if search_button and search_query.strip():
        with st.spinner("Menghubungkan ke server PubChem..."):
            query = search_query.strip()
            cid = None
            
            if query.isdigit():
                cid = int(query)
            else:
                cid = get_cid_by_name(query)
                
            if cid:
                props = get_compound_properties(cid)
                if props:
                    synonyms = get_compound_synonyms(cid)
                    compound = ChemicalCompound(
                        cid=cid,
                        name=props.get('IUPACName', query.capitalize()),
                        iupac_name=props.get('IUPACName', 'N/A'),
                        molecular_formula=props.get('MolecularFormula', 'N/A'),
                        molecular_weight=props.get('MolecularWeight', '0'),
                        synonyms=synonyms
                    )
                    
                    # Ambil daftar bahaya GHS
                    hazards = get_ghs_hazards(cid)
                    
                    # 1. Render Informasi Utama & Piktogram
                    render_compound_overview(compound)
                    st.divider()
                    render_pictograms(hazards)
                    st.divider()
                    
                    # 2. Buat Tampilan Tab Data Klasifikasi & Solusi Rekomendasi
                    tabs = st.tabs(["📋 Semua Pernyataan Bahaya", "🛡️ Panduan Keselamatan K3", "📊 Properti Fisikokimia"])
                    
                    with tabs[0]:
                        if hazards:
                            for h in hazards:
                                bg_color = '#ffcdd2' if h.severity == 'high' else '#ffe082' if h.severity == 'medium' else '#c8e6c9'
                                border_color = '#b71c1c' if h.severity == 'high' else '#e65100' if h.severity == 'medium' else '#1b5e20'
                                st.markdown(f"""
                                <div style="background: {bg_color}; border-left: 5px solid {border_color}; padding: 12px 15px; margin: 6px 0; border-radius: 0 8px 8px 0; color: #1a1a1a !important; font-weight: 500;">
                                    {h.statement} {render_hazard_badge(h.severity)}
                                </div>
                                """, unsafe_allow_html=True)
                        else:
                            st.info("Senyawa tergolong aman atau klasifikasi bahaya spesifik tidak terdaftar di PubChem.")
                            
                    with tabs[1]:
                        recs = get_safety_recommendations(hazards)
                        c1, c2 = st.columns(2)
                        with c1:
                            st.markdown("#### 🥽 Alat Pelindung Diri (APD)")
                            for p in recs['ppe']: st.markdown(f"- {p}")
                            st.markdown("#### 🖐️ Penanganan Operasional")
                            for h in recs['handling']: st.markdown(f"- {h}")
                        with c2:
                            st.markdown("#### 📦 Sistem Penyimpanan")
                            for s in recs['storage']: st.markdown(f"- {s}")
                            st.markdown("#### 🚨 Tanggap Darurat & Pembuangan")
                            for e in recs['emergency']: st.markdown(f"- {e}")
                            for d in recs['disposal']: st.markdown(f"- {d}")
                            
                    with tabs[2]:
                        render_physical_properties(props)
                        st.markdown("### 🔷 Data Rating NFPA 704")
                        nfpa = get_nfpa_diamond(cid)
                        nfpa_df = pd.DataFrame({
                            'Parameter': ['Health (Kesehatan)', 'Flammability (Kemudahan Terbakar)', 'Reactivity (Reaktivitas)', 'Special (Khusus)'],
                            'Rating': [nfpa['health'], nfpa['flammability'], nfpa['reactivity'], nfpa['special'] or 'None']
                        })
                        st.dataframe(nfpa_df, use_container_width=True, hide_index=True)
                        
                else:
                    st.error(f"❌ CID {cid} ditemukan, tetapi server properti PubChem tidak merespon.")
            else:
                st.error(f"❌ Nama senyawa atau formula '{query}' tidak ditemukan di database PubChem.")
                
    elif search_button and not search_query.strip():
        st.warning("⚠️ Masukkan nama senyawa kimia terlebih dahulu!")
        
    st.divider()
    render_footer()


if __name__ == "__main__":
    main()
