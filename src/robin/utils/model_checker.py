#!/usr/bin/env python3
"""
Standalone model validation utility for ROBIN.

This module provides functions to check for required model files and provide
helpful guidance to users on how to download them.

This module is designed to be importable without requiring the full ROBIN
package dependencies.
"""

import os
import sys
from pathlib import Path
from typing import List, Optional, Tuple

CLINVAR_VCF_GZ_NAME = "clinvar.vcf.gz"
CLINVAR_TBI_NAME = "clinvar.vcf.gz.tbi"


def get_models_directory(project_root: Optional[Path] = None) -> Path:
    """
    Get the path to the models directory.

    Uses the models module's DIR constant, which works regardless of
    where robin is run from since it's based on the installed package location.
    """
    if project_root is None:
        try:
            # Use the models module's DIR constant - this is the most reliable method
            # as it works regardless of current working directory
            from robin import models

            return models.DIR
        except ImportError:
            # Fallback: try to find models directory relative to this file
            # This should only happen in unusual circumstances
            candidate_roots = [
                Path(__file__).parent.parent.parent.parent,  # project root
                Path(__file__).parent.parent.parent,  # src directory
                Path.cwd(),  # current working directory (fallback)
            ]

            for candidate_root in candidate_roots:
                # Try project_root/src/robin/models first
                models_dir = candidate_root / "src" / "robin" / "models"
                if models_dir.exists():
                    return models_dir

                # If candidate_root is actually the src directory, try src/robin/models directly
                models_dir = candidate_root / "robin" / "models"
                if models_dir.exists():
                    return models_dir

            # Final fallback: use current working directory
            project_root = Path.cwd()

    models_dir = project_root / "src" / "robin" / "models"
    return models_dir


def get_required_models() -> List[Tuple[str, str]]:
    """
    Get list of required model files.

    Returns:
        List of tuples containing (model_name, filename)
    """
    return [
        ("general_model", "general.zip"),
        ("capper_model", "Capper_et_al_NN_v2.pkl"),
        ("pancan_model", "pancan_devel_v5i_NN_v2.pkl"),
    ]


def get_resources_directory(project_root: Optional[Path] = None) -> Path:
    """
    Get the path to the robin resources directory.

    Uses a file-relative strategy to avoid importing `robin` (which may import
    optional heavy dependencies).
    """

    # This file is src/robin/utils/model_checker.py
    # -> project root is src/robin/.. (or use `project_root` if provided)
    if project_root is None:
        this_file = Path(__file__).resolve()
        # .../src/robin/utils/model_checker.py -> .../src/robin
        robin_dir = this_file.parent.parent
        resources_dir = robin_dir / "resources"
        return resources_dir

    # Provided project_root strategy (dev / alternate layouts)
    return project_root / "src" / "robin" / "resources"


def check_clinvar_files(
    resources_dir: Optional[Path] = None,
) -> Tuple[bool, List[str], List[str]]:
    """
    Check if required ClinVar files exist.

    Returns:
        (all_present, missing_files, present_files)
    """

    resources_dir = resources_dir or get_resources_directory()
    required = [CLINVAR_VCF_GZ_NAME, CLINVAR_TBI_NAME]

    missing_files: List[str] = []
    present_files: List[str] = []

    for filename in required:
        path = resources_dir / filename
        if path.exists() and path.is_file() and path.stat().st_size > 0:
            present_files.append(filename)
        else:
            missing_files.append(filename)

    return (len(missing_files) == 0, missing_files, present_files)


def _ensure_clinvar_or_exit(resources_dir: Optional[Path] = None) -> None:
    """
    Ensure ClinVar files exist; if missing, warn and attempt download.
    """

    resources_dir = resources_dir or get_resources_directory()
    all_present, missing_files, present_files = check_clinvar_files(resources_dir)

    if all_present:
        return

    print("\n" + "=" * 60)
    print("ROBIN CLINVAR STATUS CHECK")
    print("=" * 60)
    print("❌ Missing required ClinVar VCF files:")
    for filename in missing_files:
        print(f"   ✗ {filename}")

    if present_files:
        print("\n✅ Present ClinVar files:")
        for filename in present_files:
            print(f"   ✓ {filename}")

    print(f"\n📁 ClinVar resources directory: {resources_dir}")

    print("\n⚠️  Attempting to download/generate missing ClinVar files...")
    try:
        from robin.utils.clinvar_manager import ensure_clinvar_files

        # download_if_missing=True will either download or convert formats
        ensure_clinvar_files(resources_dir=resources_dir, download_if_missing=True)
    except Exception as e:
        print(f"\n❌ Failed to download/generate ClinVar files: {e}")
        sys.exit(1)

    all_present_after, missing_after, _present_after = check_clinvar_files(
        resources_dir
    )
    if not all_present_after:
        print("\n❌ ClinVar files are still missing after download attempt:")
        for filename in missing_after:
            print(f"   ✗ {filename}")
        print("Please ensure network access is available and retry.")
        sys.exit(1)

    print("\n🎉 ClinVar files are now available.")


def check_model_files(
    project_root: Optional[Path] = None,
) -> Tuple[bool, List[str], List[str]]:
    """
    Check if all required model files are present.

    Args:
        project_root: Path to the project root directory

    Returns:
        Tuple of (all_present, missing_files, present_files)
    """
    models_dir = get_models_directory(project_root)
    required_models = get_required_models()

    missing_files = []
    present_files = []

    for model_name, filename in required_models:
        model_path = models_dir / filename
        if model_path.exists() and model_path.stat().st_size > 0:
            present_files.append(filename)
        else:
            missing_files.append(filename)

    all_present = len(missing_files) == 0
    return all_present, missing_files, present_files


def _download_missing_models(missing_files, models_dir, project_root):
    """Download missing model files (same asset manifest logic as ``robin utils update-models``)."""
    import hashlib
    import json
    import os
    import urllib.error
    import urllib.request

    print("\n🔄 Attempting to download missing models...")

    # Load assets manifest
    try:
        assets_file = project_root / "assets.json"
        if not assets_file.exists():
            print("❌ assets.json not found. Cannot download models automatically.")
            return False

        with open(assets_file, "r") as f:
            manifest = json.load(f)
    except Exception as e:
        print(f"❌ Failed to load assets manifest: {e}")
        return False

    # Asset name mapping
    asset_mapping = {
        "general.zip": "general_model",
        "Capper_et_al_NN_v2.pkl": "capper_model",
        "pancan_devel_v5i_NN_v2.pkl": "pancan_model",
    }

    github_token = os.getenv("GITHUB_TOKEN")
    if not github_token:
        print("ℹ️  No GITHUB_TOKEN found. Trying public download...")

    success_count = 0
    for filename in missing_files:
        if filename not in asset_mapping:
            print(f"⚠️  Unknown model file: {filename}")
            continue

        asset_name = asset_mapping[filename]
        if asset_name not in manifest["assets"]:
            print(f"❌ Asset '{asset_name}' not found in manifest")
            continue

        asset_info = manifest["assets"][asset_name]
        asset_url = asset_info["url"]
        expected_sha256 = asset_info["sha256"]

        target_path = models_dir / filename

        try:
            print(f"\n📥 Downloading {filename}...")

            # Download the file
            headers = {}
            if github_token:
                headers["Authorization"] = f"Bearer {github_token}"

            request = urllib.request.Request(asset_url, headers=headers)

            with urllib.request.urlopen(request) as response:
                with open(target_path, "wb") as f:
                    f.write(response.read())

            # Verify checksum
            print("🔍 Verifying checksum...")
            sha256_hash = hashlib.sha256()
            with open(target_path, "rb") as f:
                for chunk in iter(lambda: f.read(4096), b""):
                    sha256_hash.update(chunk)
            calculated_sha256 = sha256_hash.hexdigest()

            if calculated_sha256 != expected_sha256:
                print(f"❌ Checksum mismatch for {filename}")
                print(f"Expected: {expected_sha256}")
                print(f"Got:      {calculated_sha256}")
                target_path.unlink()
                continue

            print(f"✅ Successfully downloaded {filename}")
            success_count += 1

        except urllib.error.HTTPError as e:
            if e.code == 401:
                print(f"❌ Authentication failed for {filename}. Need GitHub token.")
            elif e.code == 404:
                print(f"❌ Asset not found: {filename}")
            else:
                print(f"❌ HTTP error {e.code} downloading {filename}: {e.reason}")
        except Exception as e:
            print(f"❌ Failed to download {filename}: {e}")

    return success_count == len(missing_files)


def print_model_status(project_root: Optional[Path] = None) -> bool:
    """
    Print the status of model files and return whether all are present.

    Args:
        project_root: Path to the project root directory

    Returns:
        True if all models are present, False otherwise
    """
    all_present, missing_files, present_files = check_model_files(project_root)

    print("\n" + "=" * 60)
    print("ROBIN MODEL STATUS CHECK")
    print("=" * 60)

    if all_present:
        print("✅ All required model files are present:")
        for filename in present_files:
            print(f"   ✓ {filename}")
        print("\n🎉 ROBIN is ready to run!")
        return True
    else:
        print("❌ Missing required model files:")
        for filename in missing_files:
            print(f"   ✗ {filename}")

        if present_files:
            print("\n✅ Present model files:")
            for filename in present_files:
                print(f"   ✓ {filename}")

        models_dir = get_models_directory(project_root)
        print(f"\n📁 Models directory: {models_dir}")
        print(f"📁 Current working directory: {Path.cwd()}")

        # Ask user if they want to download
        print("\n" + "=" * 60)
        print("AUTOMATIC DOWNLOAD OPTION")
        print("=" * 60)
        print("Would you like to automatically download the missing model files?")
        print("This uses the same logic as: robin utils update-models")

        try:
            response = input("\nDownload missing models? [Y/n]: ").strip().lower()
            if response in ["", "y", "yes"]:
                if _download_missing_models(
                    missing_files, models_dir, project_root or Path.cwd()
                ):
                    print("\n🎉 All models downloaded successfully!")
                    print("ROBIN is now ready to run.")
                    return True
                else:
                    print("\n⚠️  Some downloads failed. Trying alternative method...")
                    print("\n" + "=" * 60)
                    print("FALLBACK TO API METHOD")
                    print("=" * 60)
                    print("The automatic download failed. You can try the API method:")
                    print()
                    print("1. Set a GitHub token:")
                    print("   export GITHUB_TOKEN=your_github_token")
                    print()
                    print("2. Run the API download script:")
                    print("   python setup_models_api.py")
                    print()
                    print("3. Or from a shell (after pip install -e .):")
                    print("   robin utils update-models")
                    print("   robin utils update-clinvar")
                    print()
                    print("After downloading, run ROBIN again.")
                    print("=" * 60)
                    return False
            else:
                print("\n" + "=" * 60)
                print("MANUAL DOWNLOAD INSTRUCTIONS")
                print("=" * 60)
                print("To download the missing model files manually:")
                print()
                print("  robin utils update-models")
                print("  robin utils update-clinvar")
                print()
                print(
                    "For private GitHub assets, set GITHUB_TOKEN before update-models:"
                )
                print("   export GITHUB_TOKEN=your_github_token")
                print("   robin utils update-models")
                print()
                print("Legacy API script (optional):")
                print("   python setup_models_api.py")
                print()
                print("Create a token at: https://github.com/settings/tokens")
                print()
                print("After downloading, you can run ROBIN normally.")
                print("=" * 60)
                return False
        except KeyboardInterrupt:
            print("\n\n⚠️  Download cancelled by user.")
            return False


def validate_models_or_exit(project_root: Optional[Path] = None) -> None:
    """
    Check for required models and exit with helpful message if any are missing.

    This function is designed to be called at application startup.

    Args:
        project_root: Path to the project root directory
    """
    # ClinVar is required for target/variant annotation steps.
    _ensure_clinvar_or_exit()

    if not print_model_status(project_root):
        print("\n⚠️  ROBIN cannot run without the required model files.")
        print("Please download them using one of the methods above and try again.")
        sys.exit(1)


def get_model_path(
    filename: str, project_root: Optional[Path] = None
) -> Optional[Path]:
    """
    Get the full path to a model file if it exists.

    Args:
        filename: Name of the model file
        project_root: Path to the project root directory

    Returns:
        Path to the model file if it exists, None otherwise
    """
    models_dir = get_models_directory(project_root)
    model_path = models_dir / filename

    if model_path.exists() and model_path.stat().st_size > 0:
        return model_path

    return None


if __name__ == "__main__":
    # Allow running this script directly for testing
    validate_models_or_exit()
