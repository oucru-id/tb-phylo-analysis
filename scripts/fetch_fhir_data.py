import requests
import json
import argparse
import os
import sys
from urllib.parse import urljoin, urlparse
from typing import Optional

def get_headers(api_key: Optional[str] = None, bearer_token: Optional[str] = None):
    headers = {'Accept': 'application/fhir+json'}
    if bearer_token:
        headers['Authorization'] = f'Bearer {bearer_token}'
        return headers
    if api_key:
        headers['X-API-Key'] = api_key
    return headers

def fetch_oauth_token(token_url: str, client_id: str, client_secret: str, scope: Optional[str] = None) -> Optional[str]:
    try:
        data = {
            'grant_type': 'client_credentials',
            'client_id': client_id,
            'client_secret': client_secret
        }
        if scope:
            data['scope'] = scope
        resp = requests.post(token_url, data=data)
        resp.raise_for_status()
        payload = resp.json()
        return payload.get('access_token')
    except Exception:
        return None

def handle_pagination_url(base_url, next_url):
    if not next_url:
        return None
        
    if not next_url.startswith('http'):
        return urljoin(base_url, next_url)
    
    next_parsed = urlparse(next_url)
    base_parsed = urlparse(base_url)
    
    if next_parsed.netloc != base_parsed.netloc:
        return f"{base_url}/Observation?{next_parsed.query}"
    
    return next_url

def fetch_data(base_url, api_key=None, since=None, token_url=None, client_id=None, client_secret=None, scope=None):
    bearer = None
    if token_url and client_id and client_secret:
        bearer = fetch_oauth_token(token_url, client_id, client_secret, scope)
    headers = get_headers(api_key, bearer)
    base_url = base_url.rstrip('/')
    possible_bases = [base_url, f"{base_url}/fhir"]
    active_base = possible_bases[0]
    
    print(f"Connecting to FHIR Server: {active_base}")
    print(f"Searching for Patients with Genetic Data (Code 69548-6)")
    
    patients = set()
    url = f"{active_base}/Observation?code=69548-6&_count=1000"
    
    if since:
        print(f"  > Filtering for data updated after: {since}")
        url += f"&_lastUpdated=gt{since}"
    
    while url:
        try:
            resp = requests.get(url, headers=headers)
            if resp.status_code == 404 and active_base == possible_bases[0]:
                print("  Received 404 for /Observation; retrying with /fhir prefix")
                active_base = possible_bases[1]
                url = f"{active_base}/Observation?code=69548-6&_count=1000"
                if since:
                    url += f"&_lastUpdated=gt{since}"
                continue
            resp.raise_for_status()
            data = resp.json()
            
            if 'entry' in data:
                for entry in data['entry']:
                    res = entry.get('resource', {})
                    subj = res.get('subject', {}).get('reference', '')
                    if subj.startswith('Patient/'):
                        pat_id = subj.split('/')[-1]
                        patients.add(pat_id)
            
            next_link = next((l['url'] for l in data.get('link', []) if l['relation'] == 'next'), None)
            url = handle_pagination_url(active_base, next_link)

        except Exception as e:
            print(f"Error fetching search list: {e}")
            break
            
    print(f"Found {len(patients)} patients with variant data.")
    
    if len(patients) == 0:
        empty_bundle = {"resourceType": "Bundle", "type": "transaction", "entry": []}
        with open("fhir_no_data.json", "w") as f:
            json.dump(empty_bundle, f, indent=2)
        print("  No patients found; wrote fhir_no_data.json")
        return

    for pid in patients:
        print(f"Fetching full bundle for Patient {pid}")
        patient_resources = []
        
        try:
            p_resp = requests.get(f"{active_base}/Patient/{pid}", headers=headers)
            if p_resp.ok:
                patient_resources.append(p_resp.json())
        except Exception as e:
            print(f"  Error fetching Patient/{pid}: {e}")


        obs_url = f"{active_base}/Observation?patient={pid}&_count=1000"
        
        while obs_url:
            try:
                o_resp = requests.get(obs_url, headers=headers)
                if not o_resp.ok:
                    print(f"  Failed to fetch observations: {o_resp.status_code}")
                    break
                    
                o_data = o_resp.json()
                entries = o_data.get('entry', [])
                if entries:
                    print(f"  - Downloaded {len(entries)} observations...")
                    for e in entries:
                        patient_resources.append(e['resource'])
                
                next_ob_link = next((l['url'] for l in o_data.get('link', []) if l['relation'] == 'next'), None)
                obs_url = handle_pagination_url(active_base, next_ob_link)
                
            except Exception as e:
                print(f"  Error fetching variants for {pid}: {e}")
                break
            
        try:
            d_resp = requests.get(f"{active_base}/DiagnosticReport?patient={pid}", headers=headers)
            if d_resp.ok:
                d_data = d_resp.json()
                for e in d_data.get('entry', []):
                    patient_resources.append(e['resource'])
        except Exception as e:
            print(f"  Error fetching reports for {pid}: {e}")

        bundle = {
            "resourceType": "Bundle",
            "type": "transaction", 
            "entry": [{"resource": r} for r in patient_resources]
        }
        
        fname = f"{pid}.fhir.json"
        with open(fname, 'w') as f:
            json.dump(bundle, f, indent=2)
        print(f"  Saved {len(patient_resources)} resources to {fname}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--url', required=True, help="FHIR Server URL")
    parser.add_argument('--auth', help="API Key")
    parser.add_argument('--token-url', help="OAuth2 token endpoint URL")
    parser.add_argument('--client-id', help="OAuth2 client ID")
    parser.add_argument('--client-secret', help="OAuth2 client secret")
    parser.add_argument('--scope', default="openid", help="OAuth2 scope")
    parser.add_argument('--since', help="Filter data updated after YYYY-MM-DD")
    args = parser.parse_args()
    
    fetch_data(
        args.url,
        api_key=args.auth,
        since=args.since,
        token_url=args.token_url,
        client_id=args.client_id,
        client_secret=args.client_secret,
        scope=args.scope
    )
