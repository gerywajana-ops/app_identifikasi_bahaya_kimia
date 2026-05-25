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
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON/?heading=Safety+and+Hazards"
        response = requests.get(url, timeout=15)
        if response.status_code == 200:
            parse_json_recursive(response.json())
            
        if not hazards:
            url_fallback = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON"
            response_fb = requests.get(url_fallback, timeout=15)
            if response_fb.status_code == 200:
                parse_json_recursive(response_fb.json())
    except Exception as e:
        pass
        
    return hazards
