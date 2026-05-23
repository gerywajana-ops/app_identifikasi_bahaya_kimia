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
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap');
    
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
    }
    
    .main-header p {
        font-size: 1.1rem;
        opacity: 0.9;
        margin-bottom: 0;
    }
    
    .hazard-card {
        padding: 1.5rem;
        border-radius: 12px;
        margin-bottom: 1rem;
        border-left: 5px solid;
    }
    
    .hazard-physical { background-color: #fff3e0; border-left-color: #ff9800; }
    .hazard-health { background-color: #fce4ec; border-left-color: #e91e63; }
    .hazard-environmental { background-color: #e8f5e9; border-left-color: #4caf50; }
    
    .info-card {
        background-color: #f8f9fa;
        border-radius: 12px;
        padding: 1.5rem;
        margin-bottom: 1rem;
        border: 1px solid #e0e0e0;
    }
    
    .ghs-badge {
        display: inline-block;
        padding: 4px 12px;
        border-radius: 20px;
        font-size: 0.8rem;
        font-weight: 600;
        margin: 2px;
    }
    
    .badge-explosive { background: #ff5722; color: white; }
    .badge-flammable { background: #ff9800; color: white; }
    .badge-oxidizing { background: #ffeb3b; color: #333; }
    .badge-toxic { background: #f44336; color: white; }
    .badge-corrosive { background: #9c27b0; color: white; }
    .badge-environmental { background: #4caf50; color: white; }
    .badge-compressed { background: #2196f3; color: white; }
    .badge-health { background: #e91e63; color: white; }
    
    .search-container {
        background: white;
        padding: 2rem;
        border-radius: 15px;
        box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        margin-bottom: 2rem;
    }
    
    .result-container {
        animation: fadeIn 0.5s ease-in;
    }
    
    @keyframes fadeIn {
        from { opacity: 0; transform: translateY(20px); }
        to { opacity: 1; transform: translateY(0); }
    }
    
    .pictogram-img {
        border-radius: 8px;
        margin: 5px;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    }
    
    .safety-section {
        background: linear-gradient(135deg, #e3f2fd 0%, #bbdefb 100%);
        border-radius: 12px;
        padding: 1.5rem;
        margin-top: 1rem;
    }
    
    .pictogram-container {
        display: flex;
        flex-wrap: wrap;
        gap: 10px;
        margin-top: 10px;
    }
    
    .property-table {
        width: 100%;
    }
    
    .sidebar-info {
        font-size: 0.85rem;
    }
    
    .hazard-statement {
        padding: 8px 12px;
        margin: 4px 0;
        border-radius: 8px;
        font-size: 0.9rem;
    }
    
    .precautionary-statement {
        padding: 6px 10px;
        margin: 3px 0;
        border-radius: 6px;
        background-color: #e8eaf6;
        font-size: 0.85rem;
    }
    
    .footer {
        text-align: center;
        padding: 2rem;
        color: #666;
        font-size: 0.85rem;
    }
    
    .stTabs [data-baseweb="tab-list"] {
        gap: 8px;
    }
    
    .stTabs [data-baseweb="tab"] {
        border-radius: 8px 8px 0 0;
        padding: 10px 20px;
        font-weight: 500;
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

BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pulpug/rest"
PROLOG = "https://pubchem.ncbi.nlm.nih.gov"


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


def get_compound_description(cid: int) -> str:
    """Mendapatkan deskripsi senyawa"""
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON/?heading=Depositor-Supplied+Synonyms"
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            sections = data.get('Record', {}).get('Section', [])
            for section in sections:
                if section.get('TOCHeading') == 'Names and Identifiers':
                    return section.get('Description', 'Deskripsi tidak tersedia')
        return "Deskripsi tidak tersedia"
    except:
        return "Deskripsi tidak tersedia"


def get_compound_synonyms(cid: int) -> List[str]:
    """Mendapatkan sinonim senyawa"""
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/synonyms/JSON"
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            synonyms = data.get('InformationList', {}).get('Information', [{}])[0].get('Synonym', [])
            return synonyms[:10]  # Ambil 10 sinonim pertama
        return []
    except:
        return []


def get_ghs_hazards(cid: int) -> List[HazardInfo]:
    """Mendapatkan informasi bahaya GHS dari PubChem"""
    hazards = []
    try:
        # GHS Classification
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
                                    if 'Hazard Class' in string or 'Category' in string:
                                        hazards.append(parse_hazard_statement(string))
        
        # Jika tidak ada data GHS, gunakan hazard statements umum
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
                                    if string.startswith('H'):
                                        hazards.append(parse_hazard_code(string))
                                        
    except Exception as e:
        st.warning(f"Tidak dapat mengambil data GHS: {e}")
    
    return hazards


def parse_hazard_statement(statement: str) -> HazardInfo:
    """Parse pernyataan bahaya menjadi objek HazardInfo"""
    statement = statement.lower()
    
    if 'explosive' in statement or 'explosion' in statement:
        return HazardInfo('Fisika', 'Explosive', statement, 'GHS01', 'Bahaya Ledakan', 'high')
    elif 'flammable' in statement or 'flammability' in statement or 'fire' in statement:
        return HazardInfo('Fisika', 'Flammable', statement, 'GHS02', 'Bahaya Api', 'high')
    elif 'oxidiz' in statement:
        return HazardInfo('Fisika', 'Oxidizing', statement, 'GHS03', 'Mengoksidasi', 'high')
    elif 'compressed' in statement or 'gas' in statement:
        return HazardInfo('Fisika', 'Compressed Gas', statement, 'GHS04', 'Gas Bertekanan', 'medium')
    elif 'corrosive' in statement or 'corrosivity' in statement:
        return HazardInfo('Kesehatan', 'Corrosive', statement, 'GHS05', 'Korosif', 'high')
    elif 'toxic' in statement or 'acute toxicity' in statement:
        return HazardInfo('Kesehatan', 'Acute Toxicity', statement, 'GHS06', 'Toksisitas Akut', 'high')
    elif 'harmful' in statement or 'health hazard' in statement:
        return HazardInfo('Kesehatan', 'Health Hazard', statement, 'GHS08', 'Bahaya Kesehatan', 'medium')
    elif 'irritant' in statement or 'irritation' in statement:
        return HazardInfo('Kesehatan', 'Irritant', statement, 'GHS07', 'Pengiritasi', 'medium')
    elif 'environment' in statement or 'aquatic' in statement:
        return HazardInfo('Lingkungan', 'Environmental Hazard', statement, 'GHS09', 'Bahaya Lingkungan', 'medium')
    else:
        return HazardInfo('Umum', 'General Hazard', statement, 'GHS07', 'Bahaya Umum', 'low')


def parse_hazard_code(code: str) -> HazardInfo:
    """Parse kode bahaya H-code menjadi objek HazardInfo"""
    h_code = code.split(':')[0].strip()
    description = code.split(':')[1].strip() if ':' in code else code
    
    # Mapping H-codes ke hazard classes
    h_mapping = {
        'H200': ('Fisika', 'Bahaya Ledakan', 'GHS01', 'high'),
        'H201': ('Fisika', 'Bahaya Ledakan', 'GHS01', 'high'),
        'H202': ('Fisika', 'Bahaya Ledakan', 'GHS01', 'high'),
        'H203': ('Fisika', 'Bahaya Ledakan', 'GHS01', 'medium'),
        'H204': ('Fisika', 'Bahaya Ledakan', 'GHS01', 'medium'),
        'H205': ('Fisika', 'Bahaya Ledakan', 'GHS01', 'medium'),
        'H220': ('Fisika', 'Gas Sangat Mudah Terbakar', 'GHS02', 'high'),
        'H221': ('Fisika', 'Gas Mudah Terbakar', 'GHS02', 'high'),
        'H222': ('Fisika', 'Aerosol Mudah Terbakar', 'GHS02', 'high'),
        'H223': ('Fisika', 'Aerosol Mudah Terbakar', 'GHS02', 'medium'),
        'H224': ('Fisika', 'Cairan Sangat Mudah Terbakar', 'GHS02', 'high'),
        'H225': ('Fisika', 'Cairan Mudah Terbakar', 'GHS02', 'high'),
        'H226': ('Fisika', 'Cairan Mudah Terbakar', 'GHS02', 'medium'),
        'H228': ('Fisika', 'Padatan Mudah Terbakar', 'GHS02', 'medium'),
        'H240': ('Fisika', 'Mengoksidasi', 'GHS03', 'high'),
        'H241': ('Fisika', 'Mengoksidasi', 'GHS03', 'high'),
        'H242': ('Fisika', 'Mengoksidasi', 'GHS03', 'medium'),
        'H250': ('Kesehatan', 'Pirit spontan di udara', 'GHS02', 'high'),
        'H251': ('Fisika', 'Mudah terbakar; pengoksidasi', 'GHS03', 'high'),
        'H252': ('Fisika', 'Mudah terbakar dalam jumlah besar', 'GHS02', 'medium'),
        'H260': ('Fesehatan', 'Melepaskan gas mudah terbakar', 'GHS02', 'high'),
        'H261': ('Fisika', 'Melepaskan gas mudah terbakar', 'GHS02', 'medium'),
        'H270': ('Fisika', 'Mengoksidasi', 'GHS03', 'high'),
        'H271': ('Fisika', 'Mengoksidasi', 'GHS03', 'high'),
        'H272': ('Fisika', 'Mengoksidasi', 'GHS03', 'medium'),
        'H280': ('Fisika', 'Gas Bertekanan', 'GHS04', 'medium'),
        'H281': ('Fisika', 'Gas Bertekanan', 'GHS04', 'medium'),
        'H290': ('Umum', 'Korosif untuk logam', 'GHS05', 'medium'),
        'H300': ('Kesehatan', 'Toksisitas Akut', 'GHS06', 'high'),
        'H301': ('Kesehatan', 'Toksisitas Akut', 'GHS06', 'high'),
        'H302': ('Kesehatan', 'Toksisitas Akut', 'GHS07', 'medium'),
        'H304': ('Kesehatan', 'Toksisitas Akut', 'GHS08', 'high'),
        'H310': ('Kesehatan', 'Toksisitas Akut', 'GHS06', 'high'),
        'H311': ('Kesehatan', 'Toksisitas Akut', 'GHS06', 'high'),
        'H312': ('Kesehatan', 'Toksisitas Akut', 'GHS07', 'medium'),
        'H314': ('Kesehatan', 'Korosif', 'GHS05', 'high'),
        'H315': ('Kesehatan', 'Pengiritasi', 'GHS07', 'medium'),
        'H317': ('Kesehatan', 'Alergi Kulit', 'GHS07', 'medium'),
        'H318': ('Kesehatan', 'Korosif Serius Mata', 'GHS05', 'high'),
        'H319': ('Kesehatan', 'Iritasi Serius Mata', 'GHS07', 'medium'),
        'H330': ('Kesehatan', 'Toksisitas Akut', 'GHS06', 'high'),
        'H331': ('Kesehatan', 'Toksisitas Akut', 'GHS06', 'high'),
        'H332': ('Kesehatan', 'Toksisitas Akut', 'GHS07', 'medium'),
        'H334': ('Kesehatan', 'Alergi Pernapasan', 'GHS08', 'high'),
        'H335': ('Kesehatan', 'Iritasi Pernapasan', 'GHS07', 'medium'),
        'H336': ('Kesehatan', 'Efek Narotik', 'GHS07', 'medium'),
        'H340': ('Kesehatan', 'Mutagenisitas', 'GHS08', 'high'),
        'H341': ('Kesehatan', 'Mutagenisitas', 'GHS08', 'medium'),
        'H350': ('Kesehatan', 'Karsinogenik', 'GHS08', 'high'),
        'H351': ('Kesehatan', 'Karsinogenik', 'GHS08', 'medium'),
        'H360': ('Kesehatan', 'Reproduktif Toksisitas', 'GHS08', 'high'),
        'H361': ('Kesehatan', 'Reproduktif Toksisitas', 'GHS08', 'medium'),
        'H362': ('Kesehatan', 'Reproduktif Toksisitas', 'GHS08', 'medium'),
        'H370': ('Kesehatan', 'Toksisitas Organ Target', 'GHS08', 'high'),
        'H371': ('Kesehatan', 'Toksisitas Organ Target', 'GHS08', 'medium'),
        'H372': ('Kesehatan', 'Toksisitas Organ Target', 'GHS08', 'high'),
        'H373': ('Kesehatan', 'Toksisitas Organ Target', 'GHS08', 'medium'),
        'H400': ('Lingkungan', 'Toksisitas Akuatik Akut', 'GHS09', 'high'),
        'H410': ('Lingkungan', 'Toksisitas Akuatik Kronis', 'GHS09', 'high'),
        'H411': ('Lingkungan', 'Toksisitas Akuatik Kronis', 'GHS09', 'medium'),
        'H412': ('Lingkungan', 'Toksisitas Akuatik Kronis', 'GHS09', 'medium'),
        'H413': ('Lingkungan', 'Toksisitas Akuatik Kronis', 'GHS09', 'low'),
    }
    
    hazard_class, pictogram_name, pictogram_code, severity = h_mapping.get(
        h_code, ('Umum', 'Bahaya Umum', 'GHS07', 'low')
    )
    
    return HazardInfo(hazard_class, h_code, description, pictogram_code, pictogram_name, severity)


def get_pictogram_url(pictogram_code: str) -> str:
    """Mendapatkan URL pictogram GHS"""
    # Menggunakan Wikimedia Commons untuk GHS pictograms
    pictogram_urls = {
        'GHS01': 'https://upload.wikimedia.org/wikipedia/commons/thumb/9/98/GHS-pictogram-explos.svg/120px-GHS-pictogram-explos.svg.png',
        'GHS02': 'https://upload.wikimedia.org/wikipedia/commons/thumb/6/6f/GHS-pictogram-flamme.svg/120px-GHS-pictogram-flamme.svg.png',
        'GHS03': 'https://upload.wikimedia.org/wikipedia/commons/thumb/e/e5/GHS-pictogram-rondflam.svg/120px-GHS-pictogram-rondflam.svg.png',
        'GHS04': 'https://upload.wikimedia.org/wikipedia/commons/thumb/d/d1/GHS-pictogram-bottle.svg/120px-GHS-pictogram-bottle.svg.png',
        'GHS05': 'https://upload.wikimedia.org/wikipedia/commons/thumb/3/3c/GHS-pictogram-acid.svg/120px-GHS-pictogram-acid.svg.png',
        'GHS06': 'https://upload.wikimedia.org/wikipedia/commons/thumb/9/98/GHS-pictogram-skull.svg/120px-GHS-pictogram-skull.svg.png',
        'GHS07': 'https://upload.wikimedia.org/wikipedia/commons/thumb/2/20/GHS-pictogram-exclam.svg/120px-GHS-pictogram-exclam.svg.png',
        'GHS08': 'https://upload.wikimedia.org/wikipedia/commons/thumb/2/21/GHS-pictogram-silhouette.svg/120px-GHS-pictogram-silhouette.svg.png',
        'GHS09': 'https://upload.wikimedia.org/wikipedia/commons/thumb/b/b3/GHS-pictogram-pollu.svg/120px-GHS-pictogram-pollu.svg.png',
    }
    return pictogram_urls.get(pictogram_code, '')


def get_compound_2d_structure(cid: int) -> Optional[str]:
    """Mendapatkan URL gambar struktur 2D senyawa"""
    return f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/PNG?image_size=large"


def get_safety_recommendations(hazards: List[HazardInfo]) -> Dict:
    """Menghasilkan rekomendasi keselamatan berdasarkan bahaya"""
    recommendations = {
        'ppe': set(),
        'handling': set(),
        'storage': set(),
        'emergency': set(),
        'disposal': set()
    }
    
    for hazard in hazards:
        code = hazard.pictogram_code
        
        # PPE recommendations
        if code in ['GHS05', 'GHS06']:
            recommendations['ppe'].add('Sarung tangan tahan bahan kimia (neoprena/nitril)')
            recommendations['ppe'].add('Kacamata pelindung / face shield')
            recommendations['ppe'].add('Jas laboratorium tahan bahan kimia')
            recommendations['ppe'].add('Respirator (jika terdapat uap/aerosol)')
        elif code in ['GHS08', 'GHS07']:
            recommendations['ppe'].add('Sarung tangan pelindung')
            recommendations['ppe'].add('Kacamata pelindung')
            recommendations['ppe'].add('Jas laboratorium')
        elif code in ['GHS02', 'GHS03']:
            recommendations['ppe'].add('Sarung tangan tahan panas')
            recommendations['ppe'].add('Kacamata pelindung')
            recommendations['ppe'].add('Pakaian anti api')
        
        # Handling recommendations
        if code in ['GHS06', 'GHS05']:
            recommendations['handling'].add('Gunakan di dalam lemari asam/fume hood')
            recommendations['handling'].add('Hindari kontak langsung dengan kulit dan mata')
            recommendations['handling'].add('Jangan makan, minum, atau merokok saat menangani')
            recommendations['handling'].add('Cuci tangan segera setelah menangani')
        elif code == 'GHS02':
            recommendations['handling'].add('Jauhkan dari sumber api dan panas')
            recommendations['handling'].add('Gunakan peralatan anti percikan api')
            recommendations['handling'].add('Pastikan ventilasi yang baik')
        elif code == 'GHS03':
            recommendations['handling'].add('Jauhkan dari bahan mudah terbakar')
            recommendations['handling'].add('Hindari kontak dengan bahan organik')
        
        # Storage recommendations
        if code == 'GHS02':
            recommendations['storage'].add('Simpan di tempat sejuk dan kering')
            recommendations['storage'].add('Jauhkan dari sumber api dan oksidator')
            recommendations['storage'].add('Gunakan lemari penyimpanan bahan mudah terbakar')
        elif code in ['GHS06', 'GHS05']:
            recommendations['storage'].add('Simpan di lemari keamanan bahan kimia')
            recommendations['storage'].add('Kunci dan beri label dengan jelas')
            recommendations['storage'].add('Pisahkan dari bahan yang tidak kompatibel')
        elif code == 'GHS09':
            recommendations['storage'].add('Simpan di area kedap tumpah')
            recommendations['storage'].add('Pisahkan dari saluran pembuangan')
        
        # Emergency recommendations
        if code in ['GHS06', 'GHS05']:
            recommendations['emergency'].add('Tersedia stasiun pencuci mata dan safety shower')
            recommendations['emergency'].add('Tersedia kotak P3K dengan antidot spesifik')
            recommendations['emergency'].add('Ketahui nomor darurat: 118 (Pemadam), 119 (Ambulans)')
        if code == 'GHS05':
            recommendations['emergency'].add('Jika terkena kulit: Bilas dengan air minimal 15 menit')
            recommendations['emergency'].add('Jika terkena mata: Bilas mata segera dengan air bersih')
        
        # Disposal recommendations
        if code in ['GHS06', 'GHS05', 'GHS09']:
            recommendations['disposal'].add('Buang sebagai limbah B3 (Bahan Berbahaya dan Beracun)')
            recommendations['disposal'].add('Ikuti regulasi pembuangan limbah kimia setempat')
            recommendations['disposal'].add('Jangan buang ke saluran pembuangan umum')
        else:
            recommendations['disposal'].add('Buang sesuai regulasi limbah kimia')
    
    return {k: list(v) for k, v in recommendations.items()}


def get_all_hazard_info(cid: int) -> List[Dict]:
    """Mendapatkan semua informasi bahaya dari PubChem"""
    all_hazards = []
    try:
        # Coba ambil GHS Hazard Statements
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON/?heading=Hazard+Statements"
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            sections = data.get('Record', {}).get('Section', [])
            for section in sections:
                for sub in section.get('Section', []):
                    for subsub in sub.get('Section', []):
                        info = subsub.get('Information', [])
                        for item in info:
                            value = item.get('Value', {})
                            if 'StringWithMarkup' in value:
                                for markup in value['StringWithMarkup']:
                                    string = markup.get('String', '')
                                    if string and string.startswith('H'):
                                        parts = string.split(': ', 1)
                                        if len(parts) == 2:
                                            all_hazards.append({
                                                'code': parts[0].strip(),
                                                'statement': parts[1].strip()
                                            })
                                        else:
                                            all_hazards.append({
                                                'code': string[:4] if string.startswith('H') else 'N/A',
                                                'statement': string
                                            })
    except:
        pass
    
    return all_hazards


def get_precautionary_statements(cid: int) -> List[str]:
    """Mendapatkan pernyataan pencegahan (P-codes)"""
    precautions = []
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON/?heading=Precautionary+Statement+Codes"
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            sections = data.get('Record', {}).get('Section', [])
            for section in sections:
                for sub in section.get('Section', []):
                    info = sub.get('Information', [])
                    for item in info:
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
    """Mendapatkan rating NFPA 704 Diamond (jika tersedia)"""
    nfpa = {'health': 'N/A', 'flammability': 'N/A', 'reactivity': 'N/A', 'special': ''}
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON/?heading=NFPA+704+Diamond"
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            sections = data.get('Record', {}).get('Section', [])
            for section in sections:
                info = section.get('Information', [])
                for item in info:
                    value = item.get('Value', {})
                    if 'StringWithMarkup' in value:
                        for markup in value['StringWithMarkup']:
                            string = markup.get('String', '')
                            if 'Health' in string:
                                nfpa['health'] = string.split(':')[-1].strip()
                            elif 'Flammability' in string:
                                nfpa['flammability'] = string.split(':')[-1].strip()
                            elif 'Stability' in string or 'Reactivity' in string:
                                nfpa['reactivity'] = string.split(':')[-1].strip()
                            elif 'Special' in string:
                                nfpa['special'] = string.split(':')[-1].strip()
    except:
        pass
    
    return nfpa


def get_cas_number(cid: int) -> str:
    """Mendapatkan nomor CAS"""
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/xrefs/RegistryID/JSON"
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            registry_ids = data.get('InformationList', {}).get('Information', [{}])[0].get('RegistryID', [])
            for rid in registry_ids:
                if '-' in rid and len(rid.split('-')) == 3:
                    parts = rid.split('-')
                    if all(p.isdigit() for p in parts[:2]) and (parts[2].isdigit() or len(parts[2]) == 1):
                        return rid
            return registry_ids[0] if registry_ids else 'N/A'
        return 'N/A'
    except:
        return 'N/A'


# =============================================================================
# UI COMPONENTS
# =============================================================================

def render_header():
    """Render header aplikasi"""
    st.markdown("""
    <div class="main-header">
        <h1>⚗️ Chemical Hazard Identifier</h1>
        <p>Aplikasi Identifikasi Bahaya Kimia Berbasis GHS & PubChem Database</p>
    </div>
    """, unsafe_allow_html=True)


def render_search_section() -> str:
    """Render bagian pencarian"""
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
    """Render sidebar informasi"""
    with st.sidebar:
        st.markdown("""
        <div style="text-align: center; padding: 1rem 0;">
            <h2 style="color: #0f3460;">⚗️ ChemHazard ID</h2>
            <p style="font-size: 0.9rem; color: #666;">Aplikasi Identifikasi Bahaya Kimia</p>
        </div>
        """, unsafe_allow_html=True)
        
        st.divider()
        
        st.markdown("""
        ### 📖 Panduan Penggunaan
        
        **1. Masukkan Nama Senyawa**
        - Ketik nama senyawa kimia (bahasa Inggris)
        - Atau masukkan rumus molekul (contoh: H2SO4)
        - Atau nomor CID PubChem
        
        **2. Klik "Identifikasi"**
        - Aplikasi akan mencari database PubChem
        - Menampilkan informasi bahaya GHS
        
        **3. Lihat Hasil**
        - Properti fisikokimia
        - Klasifikasi bahaya GHS
        - Pictogram bahaya
        - Rekomendasi keselamatan
        
        ### 📋 Contoh Senyawa
        - Methanol, Ethanol, Acetone
        - Sulfuric acid (H2SO4)
        - Sodium hydroxide (NaOH)
        - Hydrogen peroxide (H2O2)
        - Benzene, Toluene, Xylene
        - Formaldehyde, Acetic acid
        - Ammonia, Chlorine gas
        
        ### ⚠️ Disclaimer
        Aplikasi ini menggunakan data dari **PubChem NIH** dan **GHS Classification**. 
        Selalu rujuk SDS (Safety Data Sheet) resmi untuk informasi keselamatan yang lengkap.
        
        **Sumber Data:**
        - PubChem (NCBI/NIH)
        - Globally Harmonized System (GHS)
        - Wikipedia Commons (Pictograms)
        """)
        
        st.divider()
        
        st.markdown("""
        <div class="sidebar-info">
        <b>Legenda Tingkat Bahaya:</b><br>
        🔴 <b>Tinggi</b> - Bahaya serius, penanganan khusus diperlukan<br>
        🟡 <b>Sedang</span> - Bahaya moderat, perlindungan standar<br>
        🟢 <b>Rendah</b> - Bahaya minimal, tetap berhati-hati
        </div>
        """, unsafe_allow_html=True)


def render_hazard_badge(severity: str) -> str:
    """Render badge tingkat bahaya"""
    colors = {
        'high': ('🔴', '#f44336', 'Tinggi'),
        'medium': ('🟡', '#ff9800', 'Sedang'),
        'low': ('🟢', '#4caf50', 'Rendah')
    }
    emoji, color, label = colors.get(severity, ('⚪', '#9e9e9e', 'Tidak Diketahui'))
    return f'<span style="background-color: {color}; color: white; padding: 2px 10px; border-radius: 12px; font-size: 0.75rem; font-weight: 600;">{emoji} {label}</span>'


def render_compound_overview(compound: ChemicalCompound):
    """Render overview senyawa"""
    st.markdown('<div class="result-container">', unsafe_allow_html=True)
    
    col1, col2 = st.columns([1, 2])
    
    with col1:
        # Tampilkan struktur 2D
        structure_url = get_compound_2d_structure(compound.cid)
        try:
            response = requests.get(structure_url, timeout=10)
            if response.status_code == 200:
                img = Image.open(BytesIO(response.content))
                st.image(img, caption=f"Struktur 2D: {compound.name}", use_container_width=True)
            else:
                st.info("Gambar struktur tidak tersedia")
        except:
            st.info("Gambar struktur tidak dapat dimuat")
    
    with col2:
        st.markdown(f"""
        <div class="info-card">
            <h2 style="color: #0f3460; margin-bottom: 1rem;">{compound.name}</h2>
            <p style="color: #666; font-style: italic; margin-bottom: 1rem;">{compound.iupac_name}</p>
        </div>
        """, unsafe_allow_html=True)
        
        # Tabel properti
        properties_data = {
            'Properti': [
                'CID PubChem',
                'Nama IUPAC',
                'Rumus Molekul',
                'Massa Molekul',
                'Nomor CAS'
            ],
            'Nilai': [
                str(compound.cid),
                compound.iupac_name,
                f"<b>{compound.molecular_formula}</b>",
                f"{compound.molecular_weight:.3f} g/mol" if compound.molecular_weight else 'N/A',
                get_cas_number(compound.cid)
            ]
        }
        
        df = pd.DataFrame(properties_data)
        st.markdown(df.to_html(escape=False, index=False), unsafe_allow_html=True)
    
    st.markdown('</div>', unsafe_allow_html=True)


def render_physical_properties(properties: Dict):
    """Render properti fisikokimia"""
    if not properties:
        st.info("Data properti fisikokimia tidak tersedia")
        return
    
    st.markdown("### 📊 Properti Fisikokimia")
    
    props = {
        'Massa Molekul': f"{properties.get('MolecularWeight', 'N/A')} g/mol",
        'Rumus Molekul': properties.get('MolecularFormula', 'N/A'),
        'Nama IUPAC': properties.get('IUPACName', 'N/A'),
        'Muatan': str(properties.get('Charge', 'N/A')),
        'Donor Ikatan H': str(properties.get('HBondDonorCount', 'N/A')),
        'Akseptor Ikatan H': str(properties.get('HBondAcceptorCount', 'N/A')),
        'Ikatan Rotatable': str(properties.get('RotatableBondCount', 'N/A')),
        'TPSA': f"{properties.get('TPSA', 'N/A')} Å²" if properties.get('TPSA') else 'N/A',
        'XLogP': str(properties.get('XLogP', 'N/A')),
        'SMILES Isomerik': properties.get('IsomericSMILES', 'N/A')[:50] + '...' if properties.get('IsomericSMILES') and len(properties.get('IsomericSMILES')) > 50 else properties.get('IsomericSMILES', 'N/A')
    }
    
    # Buat 3 kolom
    cols = st.columns(3)
    items = list(props.items())
    for i, (key, value) in enumerate(items):
        with cols[i % 3]:
            st.markdown(f"""
            <div style="background: #f5f5f5; padding: 10px; border-radius: 8px; margin-bottom: 8px;">
                <small style="color: #666;">{key}</small><br>
                <b style="color: #333;">{value}</b>
            </div>
            """, unsafe_allow_html=True)


def render_pictograms(hazards: List[HazardInfo]):
    """Render pictogram GHS"""
    st.markdown("### ⚠️ Pictogram Bahaya GHS")
    
    unique_pictograms = {}
    for hazard in hazards:
        if hazard.pictogram_code not in unique_pictograms:
            unique_pictograms[hazard.pictogram_code] = hazard
    
    if not unique_pictograms:
        st.info("Tidak ada data pictogram GHS yang tersedia untuk senyawa ini")
        return
    
    cols = st.columns(min(len(unique_pictograms), 4))
    
    for i, (code, hazard) in enumerate(unique_pictograms.items()):
        url = get_pictogram_url(code)
        with cols[i % 4]:
            if url:
                try:
                    response = requests.get(url, timeout=10)
                    if response.status_code == 200:
                        img = Image.open(BytesIO(response.content))
                        st.image(img, caption=hazard.pictogram_name, use_container_width=True)
                except:
                    st.warning(f"Gagal memuat: {hazard.pictogram_name}")
            
            st.markdown(render_hazard_badge(hazard.severity), unsafe_allow_html=True)
            st.markdown(f"<small>{hazard.hazard_class}</small>", unsafe_allow_html=True)


def render_hazard_classification(hazards: List[HazardInfo], cid: int):
    """Render klasifikasi bahaya lengkap dengan H-codes"""
    st.markdown("### 🏷️ Klasifikasi Bahaya GHS")
    
    # Ambil semua hazard statements
    all_hazards = get_all_hazard_info(cid)
    
    if not all_hazards and not hazards:
        st.info("Data klasifikasi bahaya GHS tidak tersedia untuk senyawa ini")
        return
    
    # Kelompokkan berdasarkan kategori
    physical_hazards = [h for h in hazards if h.hazard_class == 'Fisika']
    health_hazards = [h for h in hazards if h.hazard_class == 'Kesehatan']
    env_hazards = [h for h in hazards if h.hazard_class == 'Lingkungan']
    other_hazards = [h for h in hazards if h.hazard_class not in ['Fisika', 'Kesehatan', 'Lingkungan']]
    
    tabs = st.tabs(['📋 Semua Pernyataan Bahaya', '💨 Bahaya Fisika', '☠️ Bahaya Kesehatan', '🌿 Bahaya Lingkungan'])
    
    with tabs[0]:
        if all_hazards:
            for h in all_hazards:
                # Tentukan severity berdasarkan kode
                code_num = int(h['code'][1:]) if h['code'][1:].isdigit() else 0
                if code_num <= 205 or code_num in [300, 301, 304, 310, 311, 314, 318, 330, 331, 340, 350, 360]:
                    severity = 'high'
                    bg_color = '#ffebee'
                    border_color = '#f44336'
                elif code_num <= 373:
                    severity = 'medium'
                    bg_color = '#fff8e1'
                    border_color = '#ff9800'
                else:
                    severity = 'low'
                    bg_color = '#e8f5e9'
                    border_color = '#4caf50'
                
                st.markdown(f"""
                <div style="background: {bg_color}; border-left: 4px solid {border_color}; padding: 10px 15px; margin: 5px 0; border-radius: 0 8px 8px 0;">
                    <b>{h['code']}</b>: {h['statement']} {render_hazard_badge(severity)}
                </div>
                """, unsafe_allow_html=True)
        else:
            for h in hazards:
                st.markdown(f"""
                <div class="hazard-statement" style="background: {'#ffebee' if h.severity == 'high' else '#fff8e1' if h.severity == 'medium' else '#e8f5e9'};">
                    <b>{h.pictogram_code}</b>: {h.statement} {render_hazard_badge(h.severity)}
                </div>
                """, unsafe_allow_html=True)
    
    with tabs[1]:
        if physical_hazards:
            for h in physical_hazards:
                st.markdown(f"""
                <div class="hazard-card hazard-physical">
                    <b>{h.pictogram_name}</b> {render_hazard_badge(h.severity)}<br>
                    <small>{h.statement}</small>
                </div>
                """, unsafe_allow_html=True)
        else:
            st.info("Tidak ada data bahaya fisika yang terdeteksi")
    
    with tabs[2]:
        if health_hazards:
            for h in health_hazards:
                st.markdown(f"""
                <div class="hazard-card hazard-health">
                    <b>{h.pictogram_name}</b> {render_hazard_badge(h.severity)}<br>
                    <small>{h.statement}</small>
                </div>
                """, unsafe_allow_html=True)
        else:
            st.info("Tidak ada data bahaya kesehatan yang terdeteksi")
    
    with tabs[3]:
        if env_hazards:
            for h in env_hazards:
                st.markdown(f"""
                <div class="hazard-card hazard-environmental">
                    <b>{h.pictogram_name}</b> {render_hazard_badge(h.severity)}<br>
                    <small>{h.statement}</small>
                </div>
                """, unsafe_allow_html=True)
        else:
            st.info("Tidak ada data bahaya lingkungan yang terdeteksi")


def render_nfpa_diamond(cid: int):
    """Render NFPA 704 Diamond"""
    nfpa = get_nfpa_diamond(cid)
    
    st.markdown("### 🔷 NFPA 704 Diamond")
    
    col1, col2, col3 = st.columns([1, 1, 1])
    
    # Render diamond sederhana dengan st.metric
    with col1:
        st.metric("🔴 Health", nfpa['health'])
    with col2:
        st.metric("🔥 Flammability", nfpa['flammability'])
    with col3:
        st.metric("💥 Reactivity", nfpa['reactivity'])
    
    if nfpa['special']:
        st.info(f"**Special Hazard:** {nfpa['special']}")


def render_precautionary_statements(cid: int):
    """Render pernyataan pencegahan"""
    precautions = get_precautionary_statements(cid)
    
    if not precautions:
        return
    
    st.markdown("### 🛡️ Pernyataan Pencegahan (P-Codes)")
    
    # Kelompokkan berdasarkan kategori
    general = [p for p in precautions if p.startswith('P1') or p.startswith('P0')]
    prevention = [p for p in precautions if p.startswith('P2')]
    response_group = [p for p in precautions if p.startswith('P3')]
    storage = [p for p in precautions if p.startswith('P4')]
    disposal = [p for p in precautions if p.startswith('P5')]
    
    tabs = st.tabs(['Umum', 'Pencegahan', 'Respons', 'Penyimpanan', 'Pembuangan'])
    
    with tabs[0]:
        for p in general[:10]:
            st.markdown(f'<div class="precautionary-statement">{p}</div>', unsafe_allow_html=True)
    with tabs[1]:
        for p in prevention[:10]:
            st.markdown(f'<div class="precautionary-statement">{p}</div>', unsafe_allow_html=True)
    with tabs[2]:
        for p in response_group[:10]:
            st.markdown(f'<div class="precautionary-statement">{p}</div>', unsafe_allow_html=True)
    with tabs[3]:
        for p in storage[:10]:
            st.markdown(f'<div class="precautionary-statement">{p}</div>', unsafe_allow_html=True)
    with tabs[4]:
        for p in disposal[:10]:
            st.markdown(f'<div class="precautionary-statement">{p}</div>', unsafe_allow_html=True)


def render_safety_recommendations(hazards: List[HazardInfo]):
    """Render rekomendasi keselamatan"""
    recommendations = get_safety_recommendations(hazards)
    
    st.markdown("""
    <div class="safety-section">
        <h3 style="color: #1565c0; margin-bottom: 1rem;">🛡️ Rekomendasi Keselamatan</h3>
    """, unsafe_allow_html=True)
    
    tabs = st.tabs(['👷 APD / PPE', '⚙️ Penanganan', '📦 Penyimpanan', '🚨 Darurat', '🗑️ Pembuangan'])
    
    with tabs[0]:
        if recommendations['ppe']:
            for item in recommendations['ppe']:
                st.markdown(f"- {item}")
        else:
            st.info("Gunakan APD standar laboratorium (jas lab, kacamata pelindung, sarung tangan)")
    
    with tabs[1]:
        if recommendations['handling']:
            for item in recommendations['handling']:
                st.markdown(f"- {item}")
        else:
            st.info("Ikuti prosedur penanganan bahan kimia standar")
    
    with tabs[2]:
        if recommendations['storage']:
            for item in recommendations['storage']:
                st.markdown(f"- {item}")
        else:
            st.info("Simpan sesuai regulasi bahan kimia umum")
    
    with tabs[3]:
        if recommendations['emergency']:
            for item in recommendations['emergency']:
                st.markdown(f"- {item}")
        else:
            st.info("Tersedia stasiun pencuci mata dan safety shower")
    
    with tabs[4]:
        if recommendations['disposal']:
            for item in recommendations['disposal']:
                st.markdown(f"- {item}")
        else:
            st.info("Buang sesuai regulasi limbah kimia setempat")
    
    st.markdown("</div>", unsafe_allow_html=True)


def render_synonyms(synonyms: List[str]):
    """Render sinonim senyawa"""
    if not synonyms:
        return
    
    st.markdown("### 🏷️ Sinonim")
    synonym_cols = st.columns(5)
    for i, syn in enumerate(synonyms[:15]):
        with synonym_cols[i % 5]:
            st.markdown(f"""<span style="background: #e8eaf6; padding: 3px 10px; border-radius: 15px; font-size: 0.8rem; display: inline-block; margin: 2px;">{syn}</span>""", unsafe_allow_html=True)


def render_footer():
    """Render footer"""
    st.markdown("""
    <div class="footer">
        <p>⚗️ Chemical Hazard Identifier | Powered by <b>PubChem NIH API</b> & <b>GHS Classification</b></p>
        <p style="font-size: 0.8rem;">
        Data disediakan oleh National Center for Biotechnology Information (NCBI).<br>
        Selalu rujuk SDS (Safety Data Sheet) resmi untuk informasi keselamatan yang lengkap dan akurat.
        </p>
        <p style="font-size: 0.75rem; margin-top: 1rem;">
        Dibuat dengan ❤️ menggunakan Streamlit | © 2024
        </p>
    </div>
    """, unsafe_allow_html=True)


def render_quick_search():
    """Render tombol pencarian cepat senyawa umum"""
    st.markdown("### 🔥 Pencarian Cepat")
    
    common_compounds = [
        ('Methanol', '⚗️'),
        ('Ethanol', '🍷'),
        ('Acetone', '💅'),
        ('Sulfuric acid', '🧪'),
        ('Hydrochloric acid', '🧪'),
        ('Sodium hydroxide', '🧂'),
        ('Hydrogen peroxide', '💧'),
        ('Benzene', '⛽'),
        ('Formaldehyde', '🏠'),
        ('Ammonia', '💨'),
        ('Toluene', '🎨'),
        ('Nitric acid', '🧪'),
    ]
    
    cols = st.columns(4)
    for i, (name, icon) in enumerate(common_compounds):
        with cols[i % 4]:
            if st.button(f"{icon} {name}", key=f"quick_{i}", use_container_width=True):
                st.session_state['quick_search'] = name
                st.rerun()


# =============================================================================
# MAIN APPLICATION
# =============================================================================

def main():
    """Fungsi utama aplikasi"""
    
    # Render header
    render_header()
    
    # Render sidebar
    render_sidebar()
    
    # Render pencarian cepat
    render_quick_search()
    
    # Render bagian pencarian
    st.markdown("### 🔍 Cari Senyawa Kimia")
    
    # Check quick search
    quick_search = st.session_state.get('quick_search', '')
    
    search_query = st.text_input(
        "Nama Senyawa / Rumus / CID",
        value=quick_search,
        placeholder="Contoh: methanol, sulfuric acid, NaOH, 702 (CID)...",
        key="search_input"
    )
    
    search_button = st.button("🔍 Identifikasi Bahaya", type="primary", use_container_width=True)
    
    # Clear quick search setelah digunakan
    if quick_search:
        del st.session_state['quick_search']
    
    # Proses pencarian
    if search_button and search_query.strip():
        with st.spinner("🔬 Menganalisis senyawa dan mengidentifikasi bahaya..."):
            
            # Cek apakah input adalah CID (angka saja)
            if search_query.strip().isdigit():
                cid = int(search_query.strip())
            else:
                # Cari CID berdasarkan nama
                cid = get_cid_by_name(search_query.strip())
            
            if cid:
                # Dapatkan properti senyawa
                properties = get_compound_properties(cid)
                
                if properties:
                    # Dapatkan data tambahan
                    synonyms = get_compound_synonyms(cid)
                    hazards = get_ghs_hazards(cid)
                    
                    # Buat objek ChemicalCompound
                    compound = ChemicalCompound(
                        cid=cid,
                        name=search_query.strip().title(),
                        iupac_name=properties.get('IUPACName', 'N/A'),
                        molecular_formula=properties.get('MolecularFormula', 'N/A'),
                        molecular_weight=properties.get('MolecularWeight', 0.0),
                        synonyms=synonyms,
                        description='',
                        hazards=hazards,
                        physical_properties=properties if properties else {},
                        safety_info={},
                        pictogram_urls=[]
                    )
                    
                    # Render hasil
                    st.success(f"✅ Senyawa ditemukan! CID: {cid}")
                    
                    st.divider()
                    
                    # Tabs untuk hasil
                    tab1, tab2, tab3, tab4, tab5 = st.tabs([
                        '📋 Overview',
                        '⚠️ Bahaya & GHS',
                        '📊 Properti',
                        '🛡️ Keselamatan',
                        '🏷️ Lainnya'
                    ])
                    
                    with tab1:
                        render_compound_overview(compound)
                        if compound.synonyms:
                            render_synonyms(compound.synonyms)
                    
                    with tab2:
                        render_pictograms(compound.hazards)
                        st.divider()
                        render_hazard_classification(compound.hazards, cid)
                        st.divider()
                        render_nfpa_diamond(cid)
                        st.divider()
                        render_precautionary_statements(cid)
                    
                    with tab3:
                        render_physical_properties(compound.physical_properties)
                    
                    with tab4:
                        render_safety_recommendations(compound.hazards)
                    
                    with tab5:
                        st.markdown("### 📚 Referensi & Tautan")
                        st.markdown(f"""
                        - [📖 Lihat di PubChem](https://pubchem.ncbi.nlm.nih.gov/compound/{cid})
                        - [📄 Download SD Format](https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d)
                        - [💾 Download JSON Data](https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/JSON)
                        - [🖼️ Struktur 2D (PNG)](https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/PNG?image_size=large)
                        """)
                        
                        st.markdown("### 📋 Export Data")
                        if st.button("📥 Export sebagai JSON", use_container_width=True):
                            export_data = {
                                'cid': compound.cid,
                                'name': compound.name,
                                'iupac_name': compound.iupac_name,
                                'molecular_formula': compound.molecular_formula,
                                'molecular_weight': compound.molecular_weight,
                                'cas_number': get_cas_number(cid),
                                'synonyms': compound.synonyms,
                                'hazards': [
                                    {
                                        'class': h.hazard_class,
                                        'category': h.category,
                                        'statement': h.statement,
                                        'pictogram': h.pictogram_name,
                                        'severity': h.severity
                                    } for h in compound.hazards
                                ],
                                'properties': compound.physical_properties
                            }
                            st.download_button(
                                label="⬇️ Download JSON",
                                data=json.dumps(export_data, indent=2, ensure_ascii=False),
                                file_name=f"{compound.name.replace(' ', '_')}_hazard_data.json",
                                mime="application/json",
                                use_container_width=True
                            )
                        
                        # NFPA Summary
                        st.markdown("### 🔷 Ringkasan NFPA 704")
                        nfpa = get_nfpa_diamond(cid)
                        nfpa_df = pd.DataFrame({
                            'Parameter': ['Health (Kesehatan)', 'Flammability (Kemudahan Terbakar)', 
                                         'Reactivity (Reaktivitas)', 'Special (Khusus)'],
                            'Rating': [nfpa['health'], nfpa['flammability'], nfpa['reactivity'], nfpa['special'] or 'None']
                        })
                        st.dataframe(nfpa_df, use_container_width=True, hide_index=True)
                
                else:
                    st.error(f"❌ CID {cid} ditemukan tetapi data properti tidak tersedia.")
            else:
                st.error(f"❌ Senyawa '{search_query}' tidak ditemukan di database PubChem.\n\nCoba gunakan nama lain atau periksa ejaan.")
    
    elif search_button and not search_query.strip():
        st.warning("⚠️ Silakan masukkan nama senyawa terlebih dahulu!")
    
    # Render footer
    st.divider()
    render_footer()


if __name__ == "__main__":
    main()
