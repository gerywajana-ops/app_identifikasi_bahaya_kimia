# Chemical Hazard Identifier

Aplikasi identifikasi bahaya kimia berbasis Streamlit yang menggunakan database PubChem NIH dan sistem klasifikasi GHS (Globally Harmonized System).

## Fitur Aplikasi

- Pencarian senyawa kimia berdasarkan nama, rumus molekul, atau CID PubChem
- Tampilan struktur 2D senyawa
- Informasi properti fisikokimia lengkap
- Klasifikasi bahaya GHS dengan pictogram
- Pernyataan bahaya (H-codes) dan pencegahan (P-codes)
- NFPA 704 Diamond
- Rekomendasi keselamatan (APD, penanganan, penyimpanan, darurat, pembuangan)
- Pencarian cepat untuk senyawa umum
- Export data dalam format JSON

## Cara Menjalankan

### 1. Install dependencies

```bash
pip install -r requirements.txt
```

### 2. Jalankan aplikasi

```bash
streamlit run app.py
```

Aplikasi akan berjalan di browser pada alamat `http://localhost:8501`

## Data Sources

- **PubChem** - Database senyawa kimia oleh NCBI/NIH
- **GHS Classification** - Sistem klasifikasi bahaya global
- **Wikimedia Commons** - GHS Pictograms

## Disclaimer

Aplikasi ini menggunakan data dari PubChem NIH dan GHS Classification. Selalu rujuk SDS (Safety Data Sheet) resmi untuk informasi keselamatan yang lengkap dan akurat.
