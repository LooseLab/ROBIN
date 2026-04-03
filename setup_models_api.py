#!/usr/bin/env python3
"""
Updated setup script that uses GitHub API for private repositories.
"""

import os
import sys
import json
import hashlib
import urllib.request
import urllib.error
from pathlib import Path

def get_release_assets(github_token):
    """Get release assets from GitHub API"""
    headers = {"Authorization": f"Bearer {github_token}"}
    
    try:
        request = urllib.request.Request(
            "https://api.github.com/repos/LooseLab/ROBIN/releases/tags/v0.5",
            headers=headers
        )
        with urllib.request.urlopen(request) as response:
            release_data = json.loads(response.read().decode())
            return release_data['assets']
    except urllib.error.HTTPError as e:
        if e.code == 401:
            raise RuntimeError("Authentication failed. Check your GitHub token.")
        elif e.code == 404:
            raise RuntimeError("Release v0.5 not found.")
        else:
            raise RuntimeError(f"HTTP error {e.code}: {e.reason}")

def download_asset(asset_info, target_path, github_token):
    """Download asset using GitHub API URL"""
    headers = {
        "Authorization": f"Bearer {github_token}",
        "Accept": "application/octet-stream"
    }
    
    print(f"Downloading {asset_info['name']}...")
    print(f"API URL: {asset_info['url']}")
    
    try:
        request = urllib.request.Request(asset_info['url'], headers=headers)
        with urllib.request.urlopen(request) as response:
            with open(target_path, 'wb') as f:
                f.write(response.read())
        print(f"SUCCESS: Downloaded {asset_info['name']}")
        return True
    except urllib.error.HTTPError as e:
        print(f"ERROR: HTTP Error {e.code}: {e.reason}")
        if e.code == 401:
            print("Authentication failed. Check your GitHub token permissions.")
        elif e.code == 404:
            print("Asset not found. Check if the release exists.")
        return False
    except Exception as e:
        print(f"ERROR: Error downloading {asset_info['name']}: {e}")
        return False

def calculate_sha256(file_path):
    """Calculate SHA256 hash of a file"""
    sha256_hash = hashlib.sha256()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            sha256_hash.update(chunk)
    return sha256_hash.hexdigest()

def verify_checksum(file_path, expected_sha256):
    """Verify file checksum"""
    calculated_sha256 = calculate_sha256(file_path)
    if calculated_sha256 == expected_sha256:
        print(f"SUCCESS: Checksum verified for {file_path.name}")
        return True
    else:
        print(f"ERROR: Checksum mismatch for {file_path.name}")
        print(f"Expected: {expected_sha256}")
        print(f"Got:      {calculated_sha256}")
        return False

def main():
    """Main setup function"""
    print("Setting up Robin model assets using GitHub API...")
    
    # Get GitHub token
    github_token = os.getenv('GITHUB_TOKEN')
    if not github_token:
        print("ERROR: GITHUB_TOKEN environment variable not set")
        print("Set it with: export GITHUB_TOKEN=your_token")
        sys.exit(1)
    
    print(f"Using GitHub token: {github_token[:8]}...")
    
    # Get release assets
    try:
        assets = get_release_assets(github_token)
        print(f"Found {len(assets)} assets in release")
    except Exception as e:
        print(f"ERROR: Error getting release assets: {e}")
        sys.exit(1)
    
    # Create models directory
    models_dir = Path("src/robin/models")
    models_dir.mkdir(parents=True, exist_ok=True)
    print(f"Models will be downloaded to: {models_dir}")
    
    # Expected checksums from assets.json
    expected_checksums = {
        "general.zip": "706e043f3b9f248e51e9e57b0a01085efb93ae01d8a52fcff53418c7b39bbe26",
        "Capper_et_al_NN_v2.pkl": "00b05dc5ad2c6f617d95517422ef47181f891d60a0c87fbe91af542d574cbbd5",
        "pancan_devel_v5i_NN_v2.pkl": "4d4193ff98c5bc066a5eaf480924aeecf91d61dd85dbe6ad649419b596236806"
    }
    
    # Download each asset
    success_count = 0
    for asset in assets:
        asset_name = asset['name']
        
        # Skip SHA256SUMS file
        if asset_name == "SHA256SUMS":
            continue
            
        target_path = models_dir / asset_name
        
        # Check if file already exists
        if target_path.exists():
            print(f"SKIP: Skipping {asset_name} - already exists")
            # Verify existing file
            if asset_name in expected_checksums:
                if verify_checksum(target_path, expected_checksums[asset_name]):
                    success_count += 1
            continue
        
        # Download the asset
        if download_asset(asset, target_path, github_token):
            # Verify checksum
            if asset_name in expected_checksums:
                if verify_checksum(target_path, expected_checksums[asset_name]):
                    success_count += 1
                else:
                    # Remove corrupted file
                    target_path.unlink()
                    print(f"REMOVED: Removed corrupted file: {asset_name}")
            else:
                success_count += 1
    
    print(f"\nDownload Results: {success_count}/{len([a for a in assets if a['name'] != 'SHA256SUMS'])} successful")
    
    if success_count == len([a for a in assets if a['name'] != 'SHA256SUMS']):
        print("SUCCESS: All model assets downloaded successfully!")
    else:
        print("WARNING: Some downloads failed. Check the output above.")
        sys.exit(1)

if __name__ == "__main__":
    main()
