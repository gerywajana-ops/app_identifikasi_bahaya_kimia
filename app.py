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
    
    /* Kontainer Bahaya dengan warna teks gelap agar terbaca jelas di background pastel */
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
        color: #1a1a1a !important; /* Memaksa teks berwarna gelap */
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
    """Mendapatkan data bahaya GHS dari PubChem tanpa filter ketat agar tidak kosong"""
    hazards = []
    seen_codes = set()
    
    def parse_json_recursive(node):
        """Menyisir seluruh isi JSON tanpa peduli struktur heading-nya"""
        if isinstance(node, dict):
            for k, v in node.items():
                if k == 'String' and isinstance(v, str):
                    # Deteksi H-code resmi (Contoh: H225, H314)
                    if v.startswith('H') and ':' in v:
                        parts = v.split(':', 1)
                        h_code = parts[0].strip()
                        # Validasi format kode H standar (Panjang minimal 4 huruf, sisanya angka)
                        if len(h_code) >= 4 and h_code[1:].isdigit():
                            if h_code not in seen_codes:
                                seen_codes.add(h_code)
                                hazards.append(parse_hazard_code(v))
                                
                    # Deteksi jika teks berisi klausa bahaya langsung tanpa kode titik dua
                    elif any(kwd in v.lower() for kwd in ['flammable', 'toxic', 'corrosive', 'irritant', 'harmful', 'fatal']):
                        parsed = parse_hazard_statement(v)
                        fake_code = f"{parsed.hazard_class}_{parsed.category}"
                        if fake_code not in seen_codes:
                            seen_codes.add(fake_code)
                            hazards.append(parsed)
                else:
                    parse_json_recursive(v)
        elif isinstance(node, list):
            for item in node:
                parse_json_recursive(item)

    try:
        # Ambil ringkasan data Safety & Hazards secara utuh
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON/?heading=Safety+and+Hazards"
        response = requests.get(url, timeout=12)
        if response.status_code == 200:
            parse_json_recursive(response.json())
            
        # Jika masih kosong, tembak data mentah record compound sebagai pertahanan terakhir
        if not hazards:
            url_fallback = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON"
            response_fb = requests.get(url_fallback, timeout=12)
            if response_fb.status_code == 200:
                parse_json_recursive(response_fb.json())
                
    except Exception as e:
        st.warning(f"Sistem gagal membaca GHS: {e}")
        
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
        'H250': ('Fisika', 'Pirit spontan di udara', 'GHS02', 'high'),
        'H251': ('Fisika', 'Mudah terbakar; pengoksidasi', 'GHS03', 'high'),
        'H252': ('Fisika', 'Mudah terbakar dalam jumlah besar', 'GHS02', 'medium'),
        'H260': ('Fisika', 'Melepaskan gas mudah terbakar', 'GHS02', 'high'),
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
    """Mendapatkan URL gambar piktogram GHS format PNG yang stabil dan anti-blokir"""
    base_url = "https://raw.githubusercontent.com/ajmendez/ghs-pictograms/master/png/"
    
    pictogram_files = {
        'GHS01': 'ghs01-explos.png',
        'GHS02': 'ghs02-flamme.png',
        'GHS03': 'ghs03-rondflam.png',
        'GHS04': 'ghs04-bouteille.png',
        'GHS05': 'ghs05-acid.png',
        'GHS06': 'ghs06-skull.png',
        'GHS07': 'ghs07-exclam.png',
        'GHS08': 'ghs08-silhouet.png',
        'GHS09': 'ghs09-pollu.png',
    }
    
    filename = pictogram_files.get(pictogram_code.strip().upper())
    if filename:
        return base_url + filename
    return ""
    
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
    """Mendapatkan rating NFPA 704 Diamond dengan pemindaian JSON menyeluruh"""
    nfpa = {'health': 'N/A', 'flammability': 'N/A', 'reactivity': 'N/A', 'special': ''}
    try:
        # Gunakan API utama tanpa filter heading agar struktur datanya utuh
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON"
        response = requests.get(url, timeout=15)
        
        if response.status_code == 200:
            data = response.json()
            
            # Fungsi internal untuk mencari node "NFPA 704 Diamond" di manapun dia berada
            def find_nfpa_section(node):
                if isinstance(node, dict):
                    if node.get('TOCHeading') == 'NFPA 704 Diamond':
                        return node
                    for key, value in node.items():
                        result = find_nfpa_section(value)
                        if result:
                            return result
                elif isinstance(node, list):
                    for item in node:
                        result = find_nfpa_section(item)
                        if result:
                            return result
                return None

            nfpa_node = find_nfpa_section(data)
            
            # Jika bagian NFPA ditemukan, ekstrak nilainya
            if nfpa_node:
                info_list = nfpa_node.get('Information', [])
                for info in info_list:
                    name = info.get('Name', '')
                    value_obj = info.get('Value', {})
                    
                    if 'StringWithMarkup' in value_obj:
                        markup_list = value_obj.get('StringWithMarkup', [])
                        if markup_list:
                            raw_string = markup_list[0].get('String', '').strip()
                            
                            # PubChem memberikan nilai seperti "3 out of 4" atau hanya "2"
                            # Kita ambil karakter angka pertamanya saja
                            clean_value = raw_string[0] if raw_string and raw_string[0].isdigit() else raw_string
                            
                            if 'Health' in name:
                                nfpa['health'] = clean_value
                            elif 'Flammability' in name:
                                nfpa['flammability'] = clean_value
                            elif 'Instability' in name or 'Reactivity' in name:
                                nfpa['reactivity'] = clean_value
                            elif 'Special' in name:
                                nfpa['special'] = raw_string if raw_string != 'none' else ''
                                
    except Exception as e:
        print(f"Error deep parsing NFPA: {e}")
        
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
        🟡 <b>Sedang</b> - Bahaya moderat, perlindungan standar<br>
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
            # Tambahkan headers agar tidak diblokir oleh API PubChem
            headers = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"}
            response = requests.get(structure_url, headers=headers, timeout=10)
        
            if response.status_code == 200:
                img = Image.open(BytesIO(response.content))
                st.image(img, caption=f"Struktur 2D: {compound.name}", use_container_width=True)
            else:
                # Menampilkan status code error untuk memudahkan tracking
                st.info(f"Gambar struktur tidak tersedia (Status: {response.status_code})")
        except Exception as e:
            st.info(f"Gambar struktur tidak dapat dimuat: {e}")
    
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
                f"{float(compound.molecular_weight):.3f} g/mol" if compound.molecular_weight else 'N/A',
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
    """Render piktogram GHS berdasarkan hasil ekstraksi data hazards"""
    st.markdown("### ⚠️ Pictogram Bahaya GHS")
    if not hazards:
        st.info("Tidak ada data piktogram GHS yang tersedia (Daftar bahaya kosong).")
        return

    detected_codes = set()
    for h in hazards:
        if not h or not hasattr(h, 'statement'):
            continue
            
        stmt_text = ""
        if h.statement:
            stmt_text += " " + str(h.statement).lower()
        if h.pictogram_code:
            stmt_text += " " + str(h.pictogram_code).lower()
        if h.pictogram_name:
            stmt_text += " " + str(h.pictogram_name).lower()

        if 'flamm' in stmt_text or 'pyrophor' in stmt_text:
            detected_codes.add('GHS02')
        if 'toxic' in stmt_text or 'fatal' in stmt_text or 'poison' in stmt_text:
            detected_codes.add('GHS06')
        if 'corros' in stmt_text or 'eye damag' in stmt_text or 'skin burn' in stmt_text:
            detected_codes.add('GHS05')
        if 'explos' in stmt_text:
            detected_codes.add('GHS01')
        if 'oxidiz' in stmt_text:
            detected_codes.add('GHS03')
        if 'gas under press' in stmt_text or 'compressed gas' in stmt_text:
            detected_codes.add('GHS04')
        if 'irritat' in stmt_text or 'harmful' in stmt_text or 'sensitiz' in stmt_text:
            detected_codes.add('GHS07')
        if 'carcinogen' in stmt_text or 'mutagen' in stmt_text or 'respiratory' in stmt_text or 'target organ' in stmt_text:
            detected_codes.add('GHS08')
        if 'aquatic' in stmt_text or 'toxic to aqua' in stmt_text or 'environment' in stmt_text:
            detected_codes.add('GHS09')

    # Filter ketat agar hanya kode GHS valid yang diproses
    detected_codes = {code for code in detected_codes if code and str(code).startswith('GHS')}

    if not detected_codes:
        st.info("Senyawa tergolong aman atau tidak memerlukan piktogram bahaya GHS khusus.")
        return

    ghs_names = {
        'GHS01': 'Explosive (Mudah Meledak)',
        'GHS02': 'Flammable (Mudah Terbakar)',
        'GHS03': 'Oxidizing (Pengoksidasi)',
        'GHS04': 'Gases Under Pressure (Gas Bertekanan)',
        'GHS05': 'Corrosive (Korosif / Merusak)',
        'GHS06': 'Acute Toxicity (Beracun)',
        'GHS07': 'Harmful / Irritant (Iritasi / Bahaya Ringan)',
        'GHS08': 'Health Hazard (Bahaya Kesehatan Kronis)',
        'GHS09': 'Environmental Hazard (Bahaya Lingkungan)'
    }

    # Membuat tata letak grid menggunakan st.columns bawaan Streamlit
    valid_codes = sorted(list(detected_codes))
    cols = st.columns(min(len(valid_codes), 4))
    
    for i, code in enumerate(valid_codes):
        url = get_pictogram_url(code)
        with cols[i % 4]:
            if url:
                # Menampilkan gambar PNG langsung via komponen st.image yang aman
                st.image(url, caption=ghs_names.get(code, code), use_container_width=True)
            else:
                st.warning(f"⚠️ {ghs_names.get(code, code)}")
                
def render_hazard_classification(hazards: List[HazardInfo], cid: int):
    """Render klasifikasi bahaya lengkap dengan pemaksaan warna teks gelap agar terbaca"""
    st.markdown("### 🏷️ Klasifikasi Bahaya GHS")
    
    all_hazards = get_all_hazard_info(cid)
    
    if not all_hazards and not hazards:
        st.info("Data klasifikasi bahaya GHS tidak tersedia untuk senyawa ini")
        return
    
    physical_hazards = [h for h in hazards if h.hazard_class == 'Fisika']
    health_hazards = [h for h in hazards if h.hazard_class == 'Kesehatan']
    env_hazards = [h for h in hazards if h.hazard_class == 'Lingkungan']
    
    tabs = st.tabs(['📋 Semua Pernyataan Bahaya', '💨 Bahaya Fisika', '☠️ Bahaya Kesehatan', '🌿 Bahaya Lingkungan'])
    
    with tabs[0]:
        if all_hazards:
            for h in all_hazards:
                code_num = int(h['code'][1:]) if h['code'][1:].isdigit() else 0
                if code_num <= 205 or code_num in [300, 301, 304, 310, 311, 314, 318, 330, 331, 340, 350, 360]:
                    severity = 'high'
                    bg_color = '#ffcdd2' # Merah pastel lebih gelap dikit
                    border_color = '#b71c1c'
                elif code_num <= 373:
                    severity = 'medium'
                    bg_color = '#ffe082' # Kuning/Oren pastel
                    border_color = '#e65100'
                else:
                    severity = 'low'
                    bg_color = '#c8e6c9' # Hijau pastel
                    border_color = '#1b5e20'
                
                # Tambahan style color: #1a1a1a untuk memaksa warna tulisan menjadi gelap
                st.markdown(f"""
                <div style="background: {bg_color}; border-left: 5px solid {border_color}; padding: 12px 15px; margin: 6px 0; border-radius: 0 8px 8px 0; color: #1a1a1a !important; font-weight: 500;">
                    <b style="color: #000000 !important;">{h['code']}</b>: {h['statement']} {render_hazard_badge(severity)}
                </div>
                """, unsafe_allow_html=True)
        else:
            for h in hazards:
                bg = '#ffcdd2' if h.severity == 'high' else '#ffe082' if h.severity == 'medium' else '#c8e6c9'
                st.markdown(f"""
                <div class="hazard-statement" style="background: {bg}; color: #1a1a1a !important;">
                    <b>{h.pictogram_code}</b>: {h.statement} {render_hazard_badge(h.severity)}
                </div>
                """, unsafe_allow_html=True)
                
    # Sisa render sub-tabs (Bahaya Fisika, Kesehatan, Lingkungan)
    with tabs[1]:
        if physical_hazards:
            for h in physical_hazards:
                st.markdown(f"""
                <div class="hazard-card hazard-physical">
                    <b style="color: #1a1a1a;">{h.pictogram_name}</b> {render_hazard_badge(h.severity)}<br>
                    <span style="color: #222222;">{h.statement}</span>
                </div>
                """, unsafe_allow_html=True)
        else:
            st.info("Tidak ada data bahaya fisika yang terdeteksi")
    
    with tabs[2]:
        if health_hazards:
            for h in health_hazards:
                st.markdown(f"""
                <div class="hazard-card hazard-health">
                    <b style="color: #1a1a1a;">{h.pictogram_name}</b> {render_hazard_badge(h.severity)}<br>
                    <span style="color: #222222;">{h.statement}</span>
                </div>
                """, unsafe_allow_html=True)
        else:
            st.info("Tidak ada data bahaya kesehatan yang terdeteksi")
    
    with tabs[3]:
        if env_hazards:
            for h in env_hazards:
                st.markdown(f"""
                <div class="hazard-card hazard-environmental">
                    <b style="color: #1a1a1a;">{h.pictogram_name}</b> {render_hazard_badge(h.severity)}<br>
                    <span style="color: #222222;">{h.statement}</span>
                </div>
                """, unsafe_allow_html=True)
        else:
            st.info("Tidak ada data bahaya lingkungan yang terdeteksi")

def render_nfpa_diamond(cid: int):
    """Render NFPA 704 Diamond"""
    nfpa = get_nfpa_diamond(cid)
    
    st.markdown("### 🔷 NFPA 704 Diamond")
    
    col1, col2, col3 = st.columns([1, 1, 1])
    
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
    """Render sinonim senyawa dengan warna teks tegas dan kontras tinggi di tema gelap/terang"""
    if not synonyms:
        return
    
    st.markdown("### 🏷️ Sinonim")
    synonym_cols = st.columns(5)
    for i, syn in enumerate(synonyms[:15]):
        with synonym_cols[i % 5]:
            # Menggunakan warna background biru indigo gelap dengan teks putih tebal agar terbaca
            st.markdown(f"""
            <div style="
                background-color: #3f51b5; 
                color: #ffffff !important; 
                padding: 6px 12px; 
                border-radius: 20px; 
                font-size: 0.85rem; 
                font-weight: 600; 
                text-align: center; 
                margin: 4px 0;
                box-shadow: 0 2px 4px rgba(0,0,0,0.15);
                white-space: nowrap;
                overflow: hidden;
                text-overflow: ellipsis;
            " title="{syn}">
                {syn}
            </div>
            """, unsafe_allow_html=True)


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
        Dibuat dengan ❤️ oleh kelompok 2 | © 2026
        </p>
    </div>
    """, unsafe_allow_html=True)


def render_quick_search():
    """Render tombol pencarian cepat senyawa umum yang langsung berfungsi"""
    st.markdown("### 🔥 Pencarian Cepat")
    
    common_compounds = [
        ('Methanol', '⚗️'),
        ('Ethanol', '🍷'),
        ('Acetone', '💅'),
        ('Sulfuric acid', '🧪'),
        ('Hydrochloric acid', '🧪'),
        ('Sodium hydroxide',  '🧂'),
        ('Hydrogen peroxide',  '💧'),
        ('Benzene', '⛽'),
        ('Formaldehyde', '🏠'),
        ('Ammonia', '💨'),
        ('Toluene', '🎨'),
        ('Nitric acid',  '🧪'),
    ]
    
    cols = st.columns(4)
    for i, (name, icon) in enumerate(common_compounds):
        with cols[i % 4]:
            # Ketika tombol diklik, langsung set 'search_input' dan trigger rerun
            if st.button(f"{icon} {name}", key=f"quick_{i}", use_container_width=True):
                st.session_state['search_input'] = name
                st.session_state['trigger_search'] = True
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
    
    quick_search = st.session_state.get('quick_search', '')
    
    search_query = st.text_input(
        "Nama Senyawa / Rumus / CID",
        value=quick_search,
        placeholder="Contoh: methanol, sulfuric acid, NaOH, 702 (CID)...",
        key="search_input"
    )
    
    search_button = st.button("🔍 Identifikasi Bahaya", type="primary", use_container_width=True)
    
    if quick_search:
        del st.session_state['quick_search']
    
    if search_button and search_query.strip():
        with st.spinner("🔬 Menganalisis senyawa dan mengidentifikasi bahaya..."):
            
            if search_query.strip().isdigit():
                cid = int(search_query.strip())
            else:
                cid = get_cid_by_name(search_query.strip())
            
            if cid:
                properties = get_compound_properties(cid)
                
                if properties:
                    synonyms = get_compound_synonyms(cid)
                    hazards = get_ghs_hazards(cid)
                    
                    # Logika penentuan nama tampilan agar angka CID berubah menjadi nama asli senyawa kimia
                    compound_name = synonyms[0].title() if (search_query.strip().isdigit() and synonyms) else search_query.strip().title()
                    
                    compound = ChemicalCompound(
                        cid=cid,
                        name=compound_name,
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
                    
                    st.success(f"✅ Senyawa ditemukan! CID: {cid}")
                    st.divider()
                    
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
    
    st.divider()
    render_footer()


if __name__ == "__main__":
    main()
