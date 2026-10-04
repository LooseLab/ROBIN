#!/usr/bin/env bash
# Activate micromamba via the image ENTRYPOINT, then start ROBIN.
set -euo pipefail

_robin_prepare_config_dir() {
  local config_home="${XDG_CONFIG_HOME:-${HOME:-/home/mambauser}/.config}"
  local config_dir="${config_home}/robin"
  if mkdir -p "${config_dir}" 2>/dev/null && [[ -w "${config_dir}" ]]; then
    return 0
  fi
  # Named volumes and Docker-created mount parents are often root-owned.
  local fallback="${HOME:-/tmp}/.local"
  export XDG_CONFIG_HOME="${fallback}"
  mkdir -p "${fallback}/robin"
  echo "Warning: ${config_dir} is not writable; using ${fallback}/robin for GUI/auth state." >&2
}

_robin_prepare_config_dir

if [[ "${ROBIN_SKIP_MODEL_CHECK:-0}" != "1" ]]; then
  python - <<'PY' || true
from robin.utils.model_checker import check_model_files

ok, missing, _present = check_model_files()
if not ok:
    print(
        "ROBIN models are missing: "
        + ", ".join(missing)
        + "\nDownload with: docker compose --profile setup run --rm models"
        + "\nor: docker compose exec robin robin utils update-models",
        flush=True,
    )
PY
fi

exec "$@"
