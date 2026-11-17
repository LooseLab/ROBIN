"""Sample ID generation and identifier-manifest helpers."""

from __future__ import annotations

import base64
import hashlib
import json
import re
from dataclasses import dataclass
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

SAMPLE_IDENTIFIER_MANIFEST_FILENAME = "sample_identifier_manifest.json"

# Salt for deriving encryption key from DOB (fixed so the same DOB always produces the same key).
_IDENTIFIER_MANIFEST_KEY_SALT = b"robin_sample_manifest_v1"

# Conservative filesystem-safe sample ID (MinKNOW / folder name).
_CUSTOM_SAMPLE_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._+-]{0,127}$")
_DOB_RE = re.compile(r"^\d{4}-\d{2}-\d{2}$")

_ENCRYPTED_MANIFEST_KEYS = (
    "first_name",
    "last_name",
    "dob",
    "nhs_number",
    "notes",
)


def generate_sample_id_md5(
    test_id: str,
    *,
    first_name: str = "",
    last_name: str = "",
    date_of_birth: str = "",
) -> str:
    """Return the ROBIN MD5 sample ID used by the Sample ID generator page."""
    parts = [
        test_id.strip(),
        first_name.strip(),
        last_name.strip(),
        date_of_birth.strip(),
    ]
    if not parts[0]:
        raise ValueError("Test ID is required")
    payload = "|".join(parts)
    return hashlib.md5(payload.encode("utf-8")).hexdigest()


def md5_payload(
    test_id: str,
    *,
    first_name: str = "",
    last_name: str = "",
    date_of_birth: str = "",
) -> str:
    """Return the pipe-separated string hashed for the MD5 sample ID."""
    return "|".join(
        [
            test_id.strip(),
            first_name.strip(),
            last_name.strip(),
            date_of_birth.strip(),
        ]
    )


def normalize_dob(value: Any) -> str:
    """Normalize a date input to YYYY-MM-DD, or '' if empty."""
    if value is None:
        return ""
    if hasattr(value, "strftime"):
        return value.strftime("%Y-%m-%d")
    text = str(value).strip()
    return text


def validate_dob(dob: str) -> str:
    """Return normalized DOB or raise ValueError if invalid."""
    dob = normalize_dob(dob)
    if not dob:
        raise ValueError("Date of birth is required")
    if not _DOB_RE.match(dob):
        raise ValueError("Date of birth must be YYYY-MM-DD")
    try:
        date.fromisoformat(dob)
    except ValueError as exc:
        raise ValueError("Date of birth must be a valid calendar date") from exc
    return dob


def validate_custom_sample_id(sample_id: str) -> str:
    """Return a validated custom sample ID or raise ValueError."""
    sample_id = (sample_id or "").strip()
    if not sample_id:
        raise ValueError("MinKNOW RUN ID is required")
    if sample_id in {".", ".."} or "/" in sample_id or "\\" in sample_id:
        raise ValueError("MinKNOW RUN ID cannot contain path separators")
    if not _CUSTOM_SAMPLE_ID_RE.match(sample_id):
        raise ValueError(
            "MinKNOW RUN ID must start with a letter or digit and use only "
            "letters, digits, '.', '_', '+', or '-' (max 128 characters)"
        )
    return sample_id


def has_encrypted_identifier_fields(
    *,
    first_name: str = "",
    last_name: str = "",
    nhs_number: str = "",
    notes: str = "",
) -> bool:
    """True when the user supplied fields that must be stored encrypted."""
    return bool(
        (first_name or "").strip()
        or (last_name or "").strip()
        or (nhs_number or "").strip()
        or (notes or "").strip()
    )


@dataclass(frozen=True)
class SampleRegistrationRequest:
    """Resolved public sample ID plus optional identifier fields for the manifest."""

    sample_id: str
    id_source: str
    test_id: str = ""
    first_name: str = ""
    last_name: str = ""
    dob: str = ""
    nhs_number: str = ""
    notes: str = ""
    derived_md5: str = ""
    preview: str = ""


def build_sample_registration(
    *,
    mode: str,
    custom_sample_id: str = "",
    test_id: str = "",
    first_name: str = "",
    last_name: str = "",
    dob: Any = "",
    nhs_number: str = "",
    notes: str = "",
) -> SampleRegistrationRequest:
    """
    Validate inputs and resolve the public sample ID for registration / MinKNOW start.

    ``mode`` is ``custom`` (MinKNOW RUN ID) or ``md5``. Encrypted PII/notes require DOB.
    """
    mode_norm = (mode or "custom").strip().lower()
    if mode_norm not in {"custom", "md5"}:
        mode_norm = "custom"

    first_name_val = (first_name or "").strip()
    last_name_val = (last_name or "").strip()
    nhs_val = (nhs_number or "").strip()
    notes_val = (notes or "").strip()
    dob_str = normalize_dob(dob)

    if has_encrypted_identifier_fields(
        first_name=first_name_val,
        last_name=last_name_val,
        nhs_number=nhs_val,
        notes=notes_val,
    ):
        dob_str = validate_dob(dob_str)

    if mode_norm == "custom":
        sample_id = validate_custom_sample_id(custom_sample_id)
        test_id_val = (test_id or "").strip()
        return SampleRegistrationRequest(
            sample_id=sample_id,
            id_source="custom",
            test_id=test_id_val,
            first_name=first_name_val,
            last_name=last_name_val,
            dob=dob_str,
            nhs_number=nhs_val,
            notes=notes_val,
            preview=f"MinKNOW RUN ID: {sample_id}",
        )

    test_id_val = (test_id or "").strip()
    if not test_id_val:
        raise ValueError("Test ID is required")
    payload = md5_payload(
        test_id_val,
        first_name=first_name_val,
        last_name=last_name_val,
        date_of_birth=dob_str,
    )
    sample_id = generate_sample_id_md5(
        test_id_val,
        first_name=first_name_val,
        last_name=last_name_val,
        date_of_birth=dob_str,
    )
    return SampleRegistrationRequest(
        sample_id=sample_id,
        id_source="md5",
        test_id=test_id_val,
        first_name=first_name_val,
        last_name=last_name_val,
        dob=dob_str,
        nhs_number=nhs_val,
        notes=notes_val,
        derived_md5=sample_id,
        preview=f"MD5 of: {payload!r}",
    )


def _derive_key_from_dob(dob: str) -> bytes:
    """Derive a Fernet key from the date of birth for encrypting manifest fields."""
    from cryptography.hazmat.primitives import hashes
    from cryptography.hazmat.primitives.kdf.pbkdf2 import PBKDF2HMAC

    kdf = PBKDF2HMAC(
        algorithm=hashes.SHA256(),
        length=32,
        salt=_IDENTIFIER_MANIFEST_KEY_SALT,
        iterations=100_000,
    )
    key_bytes = kdf.derive(dob.encode("utf-8"))
    return base64.urlsafe_b64encode(key_bytes)


def encrypt_identifier_manifest_field(plaintext: str, dob: str) -> str:
    """Encrypt a string using a key derived from DOB; returns base64-encoded ciphertext."""
    from cryptography.fernet import Fernet

    key = _derive_key_from_dob(dob)
    f = Fernet(key)
    ciphertext = f.encrypt(plaintext.encode("utf-8"))
    return ciphertext.decode("ascii")


def decrypt_identifier_manifest_field(ciphertext_b64: str, dob: str) -> str:
    """Decrypt a base64-encoded ciphertext using a key derived from DOB."""
    from cryptography.fernet import Fernet

    key = _derive_key_from_dob(dob)
    f = Fernet(key)
    plaintext = f.decrypt(ciphertext_b64.encode("ascii"))
    return plaintext.decode("utf-8")


def get_test_id_from_manifest(sample_dir: Path) -> str:
    """Read test_id from sample_identifier_manifest.json in the sample directory if present."""
    manifest_path = sample_dir / SAMPLE_IDENTIFIER_MANIFEST_FILENAME
    if not manifest_path.exists():
        return ""
    try:
        with open(manifest_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        return str(data.get("test_id", "") or "").strip()
    except Exception:
        return ""


def load_manifest_encrypted_fields(
    sample_dir: Optional[Path],
) -> Optional[Dict[str, str]]:
    """Load encrypted identifier fields from sample_identifier_manifest.json if present."""
    if not sample_dir or not sample_dir.exists():
        return None
    manifest_path = sample_dir / SAMPLE_IDENTIFIER_MANIFEST_FILENAME
    if not manifest_path.exists():
        return None
    try:
        with open(manifest_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        out: Dict[str, str] = {}
        for key in _ENCRYPTED_MANIFEST_KEYS:
            value = data.get(key)
            if value:
                out[key] = str(value)
        return out or None
    except Exception:
        return None


def save_sample_identifier_manifest(
    work_directory: str | Path,
    sample_id: str,
    *,
    test_id: str = "",
    first_name: str = "",
    last_name: str = "",
    dob: str = "",
    nhs_number: str = "",
    notes: str = "",
    id_source: str = "md5",
    derived_md5: str = "",
) -> Tuple[bool, str]:
    """
    Create the sample output folder and write the identifier manifest.

    Always registers ``sample_id`` (and optional plaintext ``test_id``) so a later
    run with the same MinKNOW/sample folder name links to this registration.

    If first name, last name, hospital number, or notes are provided, ``dob`` is
    required and those fields (plus DOB) are stored encrypted under a DOB-derived key.

    Returns (success, message).
    """
    work = str(work_directory or "").strip()
    if not work:
        return False, "No output directory configured (work directory not set)."

    sample_id = (sample_id or "").strip()
    if not sample_id:
        return False, "Sample ID is required."

    test_id = (test_id or "").strip()
    first_name = (first_name or "").strip()
    last_name = (last_name or "").strip()
    dob = normalize_dob(dob)
    nhs_number = (nhs_number or "").strip()
    notes = (notes or "").strip()
    id_source = (id_source or "md5").strip().lower()
    derived_md5 = (derived_md5 or "").strip()

    store_encrypted = has_encrypted_identifier_fields(
        first_name=first_name,
        last_name=last_name,
        nhs_number=nhs_number,
        notes=notes,
    )
    if store_encrypted:
        try:
            dob = validate_dob(dob)
        except ValueError as exc:
            return False, str(exc)

    base = Path(work)
    try:
        sample_dir = base / sample_id
        sample_dir.mkdir(parents=True, exist_ok=True)
    except OSError as e:
        return False, f"Cannot create sample folder: {e}"

    manifest: Dict[str, Any] = {
        "sample_id": sample_id,
        "id_source": id_source if id_source in {"md5", "custom"} else "custom",
        "created_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%fZ"),
    }
    if test_id:
        manifest["test_id"] = test_id
    if derived_md5:
        manifest["derived_md5"] = derived_md5

    if store_encrypted:
        try:
            if first_name:
                manifest["first_name"] = encrypt_identifier_manifest_field(
                    first_name, dob
                )
            if last_name:
                manifest["last_name"] = encrypt_identifier_manifest_field(
                    last_name, dob
                )
            manifest["dob"] = encrypt_identifier_manifest_field(dob, dob)
            if nhs_number:
                manifest["nhs_number"] = encrypt_identifier_manifest_field(
                    nhs_number, dob
                )
            if notes:
                manifest["notes"] = encrypt_identifier_manifest_field(notes, dob)
        except Exception as e:
            return False, f"Encryption failed: {e}"

    manifest_path = sample_dir / SAMPLE_IDENTIFIER_MANIFEST_FILENAME
    try:
        with open(manifest_path, "w", encoding="utf-8") as f:
            json.dump(manifest, f, indent=2)
    except OSError as e:
        return False, f"Cannot write manifest: {e}"

    if store_encrypted:
        return True, f"Registered sample ID with encrypted identifiers at {sample_dir}"
    return True, f"Registered sample ID at {sample_dir}"


def save_sample_registration(
    work_directory: str | Path,
    registration: SampleRegistrationRequest,
) -> Tuple[bool, str]:
    """Persist a resolved :class:`SampleRegistrationRequest` to the work directory."""
    return save_sample_identifier_manifest(
        work_directory,
        registration.sample_id,
        test_id=registration.test_id,
        first_name=registration.first_name,
        last_name=registration.last_name,
        dob=registration.dob,
        nhs_number=registration.nhs_number,
        notes=registration.notes,
        id_source=registration.id_source,
        derived_md5=registration.derived_md5,
    )
