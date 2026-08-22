import json
import logging
import os
import time
from dataclasses import dataclass
from typing import Optional

import requests


@dataclass
class AuthToken:
    access_token: str
    token_type: str = "bearer"
    raw: Optional[dict] = None


class MNPFlexClient:
    """
    Standalone client for MNP-Flex API:
      - authenticate (password grant)
      - upload sample
      - wait for results (bundle summary + plots)
      - save results to disk
      - delete uploaded sample
    """

    def __init__(
        self,
        *,
        base_url: str = "https://app.epignostix.com",
        username: str,
        password: str,
        verify_ssl: bool = False,
        timeout: int = 120,
        client_id: str = "",
        client_secret: str = "",
        scope: str = "",
    ):
        self.base_url = base_url.rstrip("/")
        self.verify_ssl = verify_ssl
        self.timeout = timeout

        self.session = requests.Session()
        self.session.headers.update(
            {
                "accept": "application/json",
                "user-agent": "mnpflex-standalone-client/0.1",
            }
        )
        self._token: Optional[AuthToken] = None
        self.authenticate(
            username=username,
            password=password,
            client_id=client_id,
            client_secret=client_secret,
            scope=scope,
        )

    # ----------------------------
    # Auth + helpers
    # ----------------------------
    def _url(self, path: str) -> str:
        return f"{self.base_url}{path}"

    def _auth_headers(self) -> dict:
        if not self._token or not self._token.access_token:
            raise RuntimeError("Not authenticated. Call authenticate() first.")
        return {"Authorization": f"Bearer {self._token.access_token}"}

    def authenticate(
        self,
        *,
        username: str,
        password: str,
        client_id: str = "",
        client_secret: str = "",
        scope: str = "",
    ) -> AuthToken:
        token_url = self._url("/api/v1/auth/token")
        data = {
            "grant_type": "password",
            "username": username,
            "password": password,
            "client_id": client_id,
            "client_secret": client_secret,
            "scope": scope,
        }
        resp = self.session.post(
            token_url,
            data=data,
            verify=self.verify_ssl,
            timeout=self.timeout,
        )
        resp.raise_for_status()
        payload = resp.json()
        token = AuthToken(
            access_token=payload.get("access_token", ""),
            token_type=payload.get("token_type", "bearer"),
            raw=payload,
        )
        if not token.access_token:
            raise RuntimeError(f"Authentication failed: {payload}")
        self._token = token
        return token

    # ----------------------------
    # API endpoints
    # ----------------------------
    def upload_sample(
        self,
        bed_file_path: str,
        *,
        sample_identifier: str,
        workflow_id: int,
        used_technology: str = "NA",
        target_coverage: str = "NA",
        extraction_type: str = "NA",
        sex: str = "NA",
        keep_filename: bool = False,
        localisation: str = "",
        diagnosis: str = "",
        age: str = "",
        comment: str = "",
        group_id: str = "",
    ) -> dict:
        if not os.path.exists(bed_file_path):
            raise FileNotFoundError(f"BED file not found: {bed_file_path}")

        url = self._url("/api/v1/mnpflex_sample")
        params = {
            "sample_identifier": sample_identifier,
            "used_technology": used_technology,
            "target_coverage": target_coverage,
            "extraction_type": extraction_type,
            "sex": sex,
            "workflow_id": str(workflow_id),
            "keep_filename": str(keep_filename).lower(),
        }
        data = {
            "localisation": localisation,
            "diagnosis": diagnosis,
            "age": age,
            "comment": comment,
            "group_id": group_id,
        }

        filename = os.path.basename(bed_file_path)
        is_gz = filename.endswith(".gz")
        mime = "application/x-gzip" if is_gz else "text/plain"

        with open(bed_file_path, "rb") as f:
            files = {"bed_file": (filename, f, mime)}
            resp = self.session.put(
                url,
                params=params,
                data=data,
                files=files,
                headers=self._auth_headers(),
                verify=self.verify_ssl,
                timeout=max(self.timeout, 300),
                allow_redirects=True,
            )
        resp.raise_for_status()
        return resp.json()

    def list_workflow_runs_by_entity(self, entity_id: int) -> list[dict]:
        url = self._url(f"/api/v1/workflow_runs/by_entity/{entity_id}")
        resp = self.session.get(
            url,
            headers=self._auth_headers(),
            verify=self.verify_ssl,
            timeout=self.timeout,
        )
        resp.raise_for_status()
        out = resp.json()
        if isinstance(out, list):
            return out
        if isinstance(out, dict) and "items" in out:
            return out["items"]
        raise TypeError(f"Unexpected response shape: {type(out)}")

    def get_workflow_run_details(self, workflow_run_id: int) -> dict:
        url = self._url(f"/api/v1/workflow_runs/{workflow_run_id}")
        resp = self.session.get(
            url,
            headers=self._auth_headers(),
            verify=self.verify_ssl,
            timeout=self.timeout,
        )
        resp.raise_for_status()
        return resp.json()

    def get_bundle_summary(self, workflow_run_id: int, task_result_id: int) -> dict:
        url = self._url(
            f"/api/v1/mnpflex_sample/analysis/bundle_summary/{workflow_run_id}/{task_result_id}"
        )
        try:
            resp = self.session.get(
                url,
                headers=self._auth_headers(),
                verify=self.verify_ssl,
                timeout=self.timeout,
            )
            resp.raise_for_status()
            return resp.json()
        except requests.HTTPError:
            status = getattr(resp, "status_code", None) if "resp" in locals() else None
            # Response text is usually safe; still trim to a small window.
            text_snippet = ""
            try:
                text_snippet = (resp.text or "")[:500] if "resp" in locals() else ""
            except Exception:
                text_snippet = ""

            log_fn = (
                logging.warning
                if status is not None and status >= 500
                else logging.error
            )
            # Avoid exc_info on every transient 5xx; the bulk runner will capture the final traceback.
            log_kwargs = {}
            if not (status is not None and status >= 500):
                log_kwargs["exc_info"] = True

            log_fn(
                "MNPFlex get_bundle_summary failed: status=%s url=%s body_snippet=%r",
                status,
                url,
                text_snippet,
                **log_kwargs,
            )
            raise

    def get_plot(
        self,
        *,
        plot_type: str,
        workflow_run_id: int,
        task_result_id: int,
        save_path: Optional[str] = None,
        optional: bool = False,
    ) -> Optional[bytes]:
        plot_endpoints = {
            "qc_coverage": "qc_coverage_plot",
            "qc_methylation_density": "qc_methylation_density_plot",
            "mgmt_region": "mgmt_region_plot",
        }
        if plot_type not in plot_endpoints:
            raise ValueError(
                f"plot_type must be one of {list(plot_endpoints.keys())}, got {plot_type}"
            )
        endpoint = plot_endpoints[plot_type]
        url = self._url(
            f"/api/v1/mnpflex_sample/analysis/{endpoint}/{workflow_run_id}/{task_result_id}"
        )
        resp = self.session.get(
            url,
            headers=self._auth_headers(),
            verify=self.verify_ssl,
            timeout=self.timeout,
        )
        if optional and resp.status_code == 404:
            logging.warning(
                "MNPFlex get_plot(%s) returned 404 (no artifact); skipping. "
                "workflow_run_id=%s task_result_id=%s url=%s",
                plot_type,
                workflow_run_id,
                task_result_id,
                url,
            )
            return None
        resp.raise_for_status()
        if save_path:
            with open(save_path, "wb") as f:
                f.write(resp.content)
        return resp.content

    def delete_sample(self, sample_id: int) -> dict:
        url = self._url(f"/api/v1/mnpflex_sample/{sample_id}")
        resp = self.session.delete(
            url,
            headers=self._auth_headers(),
            verify=self.verify_ssl,
            timeout=self.timeout,
        )
        resp.raise_for_status()
        return resp.json()

    # ----------------------------
    # High-level workflow
    # ----------------------------
    def wait_for_results(
        self,
        *,
        sample_id: int,
        workflow_id: Optional[int] = None,
        poll_seconds: int = 30,
        max_wait_seconds: int = 60 * 60,
    ) -> tuple[int, int]:
        """
        Wait for workflow to finish and return (workflow_run_id, task_result_id).
        """
        deadline = time.time() + max_wait_seconds

        workflow_run_id = None
        while time.time() < deadline:
            runs = self.list_workflow_runs_by_entity(sample_id)
            if runs:
                if workflow_id is None:
                    workflow_run_id = runs[0].get("id")
                else:
                    workflow_run_id = next(
                        (
                            r.get("id")
                            for r in runs
                            if r.get("workflow_id") == workflow_id
                        ),
                        None,
                    )
            if workflow_run_id:
                break
            time.sleep(poll_seconds)

        if not workflow_run_id:
            raise TimeoutError("No workflow run found before timeout.")

        while time.time() < deadline:
            details = self.get_workflow_run_details(workflow_run_id)
            task_runs = details.get("task_runs", [])
            bundle = next(
                (
                    tr
                    for tr in task_runs
                    if tr.get("task", {}).get("task_name") == "mnpflex_bundle"
                ),
                None,
            )
            if bundle:
                status = bundle.get("status")
                task_result_id = bundle.get("task_result_id")
                if task_result_id and status in {"complete", "completed", "success"}:
                    return workflow_run_id, task_result_id
                if status in {"failed", "error", "canceled", "cancelled"}:
                    error_detail = (
                        bundle.get("error_message")
                        or bundle.get("message")
                        or bundle.get("details")
                    )
                    raise RuntimeError(
                        "Workflow task failed. "
                        f"status={status}, task_result_id={task_result_id}, "
                        f"detail={error_detail}, bundle={bundle}"
                    )
            time.sleep(poll_seconds)

        raise TimeoutError("Results not ready before timeout.")

    def retrieve_and_save_results(
        self,
        *,
        sample_id: int,
        workflow_id: Optional[int],
        output_dir: str,
        poll_seconds: int = 30,
        max_wait_seconds: int = 60 * 60,
    ) -> dict:
        os.makedirs(output_dir, exist_ok=True)
        workflow_run_id, task_result_id = self.wait_for_results(
            sample_id=sample_id,
            workflow_id=workflow_id,
            poll_seconds=poll_seconds,
            max_wait_seconds=max_wait_seconds,
        )

        # MNP-Flex occasionally returns transient 5xx for result materialization.
        # Keep a small retry window so a single backend hiccup doesn't fail the whole sample.
        def _retry(
            callable_fn,
            *,
            what: str,
            attempts: int = 3,
            base_sleep_s: float = 2.0,
        ):
            last_exc: Optional[Exception] = None
            for attempt in range(1, attempts + 1):
                try:
                    return callable_fn()
                except requests.HTTPError as e:
                    last_exc = e
                    status = getattr(getattr(e, "response", None), "status_code", None)
                    if status is not None and status >= 500 and attempt < attempts:
                        sleep_s = base_sleep_s * attempt
                        logging.warning(
                            "MNPFlex transient 5xx during %s for sample_id=%s (attempt %d/%d, status=%s). Sleeping %.1fs",
                            what,
                            sample_id,
                            attempt,
                            attempts,
                            status,
                            sleep_s,
                        )
                        time.sleep(sleep_s)
                        continue
                    # Either non-5xx, or we ran out of retries.
                    logging.debug(
                        "MNPFlex not retrying %s for sample_id=%s (attempt %d/%d, status=%s)",
                        what,
                        sample_id,
                        attempt,
                        attempts,
                        status,
                    )
                    raise
                except Exception as e:
                    last_exc = e
                    raise
            if last_exc:
                raise last_exc

        summary = _retry(
            lambda: self.get_bundle_summary(workflow_run_id, task_result_id),
            what="get_bundle_summary",
            attempts=6,
            base_sleep_s=5.0,
        )
        summary_path = os.path.join(output_dir, "bundle_summary.json")
        with open(summary_path, "w") as f:
            json.dump(summary, f, indent=2)

        _retry(
            lambda: self.get_plot(
                plot_type="qc_coverage",
                workflow_run_id=workflow_run_id,
                task_result_id=task_result_id,
                save_path=os.path.join(output_dir, "qc_coverage_plot.png"),
            ),
            what="get_plot(qc_coverage)",
        )
        _retry(
            lambda: self.get_plot(
                plot_type="qc_methylation_density",
                workflow_run_id=workflow_run_id,
                task_result_id=task_result_id,
                save_path=os.path.join(output_dir, "qc_methylation_density_plot.png"),
            ),
            what="get_plot(qc_methylation_density)",
        )
        # MGMT region PNG is not always exposed at this analysis URL (workflow/version/sample).
        # Treat 404 as missing optional artifact so QC + bundle still succeed.
        _retry(
            lambda: self.get_plot(
                plot_type="mgmt_region",
                workflow_run_id=workflow_run_id,
                task_result_id=task_result_id,
                save_path=os.path.join(output_dir, "mgmt_region_plot.png"),
                optional=True,
            ),
            what="get_plot(mgmt_region)",
        )

        return {
            "workflow_run_id": workflow_run_id,
            "task_result_id": task_result_id,
            "bundle_summary_path": summary_path,
        }

    def upload_retrieve_cleanup(
        self,
        *,
        bed_file_path: str,
        sample_identifier: str,
        workflow_id: int,
        output_dir: str,
        poll_seconds: int = 30,
        max_wait_seconds: int = 60 * 60,
    ) -> dict:
        upload_result = self.upload_sample(
            bed_file_path,
            sample_identifier=sample_identifier,
            workflow_id=workflow_id,
        )
        sample_id = upload_result.get("id") or upload_result.get("sample_id")
        if not sample_id:
            raise RuntimeError(f"Upload response missing sample id: {upload_result}")

        result = self.retrieve_and_save_results(
            sample_id=sample_id,
            workflow_id=workflow_id,
            output_dir=output_dir,
            poll_seconds=poll_seconds,
            max_wait_seconds=max_wait_seconds,
        )
        self.delete_sample(sample_id)
        return {"sample_id": sample_id, **result}


if __name__ == "__main__":
    # Example usage (use env vars to avoid hardcoding secrets)
    USERNAME = os.environ.get("EPIGNOSTIX_USERNAME", "")
    PASSWORD = os.environ.get("EPIGNOSTIX_PASSWORD", "")

    client = MNPFlexClient(
        username=USERNAME,
        password=PASSWORD,
        verify_ssl=True,
    )

    result = client.upload_retrieve_cleanup(
        bed_file_path="PATH_TO_BEDFILE/file.bed.gz",
        sample_identifier="YOUR_SAMPLE_IDENTIFIER",
        workflow_id=18,
        output_dir="mnpflex_results_output",
    )
    print("Done:", result)
