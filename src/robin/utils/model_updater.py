from __future__ import annotations

import hashlib
import json
import os
import re
import urllib.error
import urllib.request
from pathlib import Path
from typing import Optional, Tuple, List, Dict, Any

# https://github.com/{owner}/{repo}/releases/download/{tag}/{filename}
_GITHUB_RELEASE_DOWNLOAD_RE = re.compile(
    r"^https://github\.com/(?P<owner>[^/]+)/(?P<repo>[^/]+)/releases/download/(?P<tag>[^/]+)/(?P<name>[^/?#]+)$"
)


def _default_repo_root_guess() -> Optional[Path]:
    """
    Best-effort guess for repo root when running from source checkout.
    When installed as a package, this may not exist.
    """
    try:
        here = Path(__file__).resolve()
        # .../src/robin/utils/model_updater.py -> .../ (repo root)
        return here.parents[3]
    except Exception:
        return None


def _resolve_assets_manifest_path(manifest_path: Optional[str] = None) -> Optional[Path]:
    """
    Resolve assets manifest path.

    Priority:
    1) explicit argument
    2) env ROBIN_ASSETS_MANIFEST
    3) repo-root assets.json (dev checkout)
    4) packaged robin.resources/assets.json (installed package)
    """
    if manifest_path:
        return Path(manifest_path).expanduser().resolve()

    env = os.getenv("ROBIN_ASSETS_MANIFEST", "").strip()
    if env:
        return Path(env).expanduser().resolve()

    root = _default_repo_root_guess()
    if root:
        candidate = root / "assets.json"
        if candidate.exists():
            return candidate.resolve()

    # Installed-package fallback: look for robin.resources/assets.json
    try:
        import importlib.resources as importlib_resources

        p = importlib_resources.files("robin.resources").joinpath("assets.json")
        # `p` may be a Traversable; try converting to filesystem path.
        try:
            fs = Path(p)  # type: ignore[arg-type]
            if fs.exists():
                return fs.resolve()
        except TypeError:
            # Not a real filesystem path (e.g. zipped); handle below by reading content.
            return None
    except Exception:
        pass

    return None


def _load_assets_manifest(manifest_path: Optional[str]) -> Dict[str, Any]:
    """
    Load assets manifest JSON.
    If manifest_path is None, attempt to load packaged robin.resources/assets.json.
    """
    if manifest_path:
        mp = Path(manifest_path).expanduser().resolve()
        with mp.open("r", encoding="utf-8") as f:
            return json.load(f)

    # Packaged fallback
    import importlib.resources as importlib_resources

    txt = importlib_resources.files("robin.resources").joinpath("assets.json").read_text(encoding="utf-8")
    return json.loads(txt)


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _parse_github_release_download_url(url: str) -> Optional[Tuple[str, str, str, str]]:
    m = _GITHUB_RELEASE_DOWNLOAD_RE.match(url.strip())
    if not m:
        return None
    g = m.groupdict()
    return (g["owner"], g["repo"], g["tag"], g["name"])


def _stream_http_response(resp, target_path: Path, *, label: str) -> None:
    target_path.parent.mkdir(parents=True, exist_ok=True)
    total: Optional[int] = None
    try:
        cl = resp.headers.get("Content-Length")
        if cl:
            total = int(cl)
    except Exception:
        total = None

    try:
        import click  # type: ignore
    except Exception:
        click = None  # type: ignore

    chunk_size = 1024 * 1024  # 1 MiB
    downloaded = 0

    if click is not None and total and total > 0:
        with click.progressbar(length=total, label=label, show_eta=True, show_percent=True) as bar:
            with target_path.open("wb") as f:
                while True:
                    chunk = resp.read(chunk_size)
                    if not chunk:
                        break
                    f.write(chunk)
                    downloaded += len(chunk)
                    bar.update(len(chunk))
    else:
        last_reported_mb = -1
        with target_path.open("wb") as f:
            while True:
                chunk = resp.read(chunk_size)
                if not chunk:
                    break
                f.write(chunk)
                downloaded += len(chunk)
                if click is not None:
                    mb = downloaded // (1024 * 1024)
                    if mb // 16 != last_reported_mb // 16:
                        last_reported_mb = mb
                        click.echo(f"{label}: {mb} MiB downloaded...")


def _download_github_release_via_api(
    owner: str,
    repo: str,
    tag: str,
    filename: str,
    target_path: Path,
    token: str,
    *,
    label: str,
) -> None:
    """Download a release asset using the GitHub API (required for private repositories)."""
    from urllib.parse import quote

    # Keep dots etc. in tags like v0.5 (quote(..., safe="") would produce v0%2E5).
    safe_tag = quote(tag, safe="/:@+-.~_|")
    meta_url = f"https://api.github.com/repos/{owner}/{repo}/releases/tags/{safe_tag}"
    meta_req = urllib.request.Request(
        meta_url,
        headers={
            "Authorization": f"Bearer {token}",
            "Accept": "application/vnd.github+json",
            "X-GitHub-Api-Version": "2022-11-28",
            "User-Agent": "ROBIN-model-updater",
        },
    )
    with urllib.request.urlopen(meta_req) as meta_resp:
        release = json.loads(meta_resp.read().decode())

    asset = None
    for a in release.get("assets") or []:
        if a.get("name") == filename:
            asset = a
            break
    if not asset:
        raise FileNotFoundError(
            f"No release asset named {filename!r} on {owner}/{repo} tag {tag!r}."
        )

    api_asset_url = str(asset["url"])
    dl_req = urllib.request.Request(
        api_asset_url,
        headers={
            "Authorization": f"Bearer {token}",
            "Accept": "application/octet-stream",
            "User-Agent": "ROBIN-model-updater",
        },
    )
    with urllib.request.urlopen(dl_req) as resp:
        _stream_http_response(resp, target_path, label=label)


def _download_manifest_asset(
    url: str,
    target_path: Path,
    github_token: Optional[str],
    *,
    label: str,
) -> None:
    """
    Download using the manifest URL. For github.com release assets with a token,
    prefer the GitHub API so private repositories work (anonymous release URLs 404).
    """
    parsed = _parse_github_release_download_url(url)
    if github_token and parsed:
        owner, repo, tag, name = parsed
        try:
            _download_github_release_via_api(
                owner, repo, tag, name, target_path, github_token, label=label
            )
            return
        except urllib.error.HTTPError as e:
            if e.code != 404:
                raise
        except FileNotFoundError:
            raise

    headers: Dict[str, str] = {}
    if github_token:
        headers["Authorization"] = f"Bearer {github_token}"
    req = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(req) as resp:
        _stream_http_response(resp, target_path, label=label)


def _download(url: str, target_path: Path, github_token: Optional[str], *, label: str = "Downloading") -> None:
    _download_manifest_asset(url, target_path, github_token, label=label)


def update_models(
    models_dir: Path,
    github_token: Optional[str] = None,
    manifest_path: Optional[str] = None,
    overwrite: bool = False,
) -> Tuple[bool, List[str]]:
    """
    Download/verify required model assets into models_dir.

    Returns:
        (ok, messages)
    """
    messages: List[str] = []
    token = github_token if github_token is not None else os.getenv("GITHUB_TOKEN")

    # Lazy import so importing `robin` doesn't require heavy deps.
    try:
        from robin.utils.model_checker import get_required_models
    except Exception as e:
        return False, [f"Could not load required models list: {e}"]

    # Resolve manifest path (or fall back to packaged resources)
    mp = _resolve_assets_manifest_path(manifest_path)
    try:
        manifest = _load_assets_manifest(str(mp) if mp else None)
    except Exception:
        return False, [
            "Could not locate assets manifest (assets.json).",
            "Pass --manifest /path/to/assets.json or set ROBIN_ASSETS_MANIFEST=/path/to/assets.json.",
        ]

    if mp:
        messages.append(f"Using assets manifest: {mp}")
    else:
        messages.append("Using bundled package assets manifest (robin.resources/assets.json).")

    models_dir = Path(models_dir).expanduser().resolve()
    models_dir.mkdir(parents=True, exist_ok=True)

    required = get_required_models()
    for asset_key, filename in required:
        if asset_key not in manifest.get("assets", {}):
            messages.append(f"Asset key not found in manifest: {asset_key}")
            return False, messages

        asset_info = manifest["assets"][asset_key]
        url = str(asset_info.get("url") or "")
        expected_sha256 = str(asset_info.get("sha256") or "")
        if not url or not expected_sha256:
            messages.append(f"Manifest entry incomplete for {asset_key} (missing url/sha256).")
            return False, messages

        target_path = models_dir / filename
        if target_path.exists() and not overwrite:
            messages.append(f"Skipping {filename} (already exists).")
            continue
        try:
            _download(url, target_path, token, label=f"Downloading {filename}")
            got = _sha256_file(target_path)
            if got != expected_sha256:
                try:
                    target_path.unlink()
                except Exception:
                    pass
                messages.append(f"Checksum mismatch for {filename} (expected {expected_sha256}, got {got}).")
                return False, messages
            messages.append(f"Downloaded {filename}.")
        except Exception as e:
            messages.append(f"Failed to download {filename} from {url!r}: {e}")
            err_s = str(e).lower()
            if "404" in err_s:
                messages.append(
                    "HTTP 404 usually means the URL in your assets.json does not match a published release "
                    "(e.g. an older ROBIN install still points at v0.0.1 or another repo). "
                    "Reinstall from current source (`pip install -e .`) or run "
                    "`robin utils update-models --manifest /path/to/ROBIN/assets.json`. "
                    "For private repos, also set GITHUB_TOKEN."
                )
            return False, messages

    return True, messages

