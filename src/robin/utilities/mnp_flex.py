import requests
from typing import Any, Dict, Tuple, List
from collections import deque
import logging
import os


class APIClient:
    """
    A client for interacting with the API endpoints.

    Endpoints include:
    - Authentication
    - Sample management
    - User management
    """

    def __init__(self, base_url: str, verify_ssl: bool = True):
        """
        Initializes the API client.

        :param base_url: The base URL for the API (e.g., "https://api.example.com")
        :param verify_ssl: Whether to verify SSL certificates (default is True)
        """
        self.base_url = base_url
        self.token = None
        self.verify_ssl = verify_ssl  # Set the SSL verification globally

    def authenticate(
        self,
        username: str,
        password: str,
        client_id: str,
        client_secret: str,
        grant_type: str = "",
        scope: str = "",
    ) -> None:
        """
        Authenticates a user and stores the access token.

        :param username: The username of the user
        :param password: The password of the user
        :param client_id: The client ID for the OAuth2 authentication
        :param client_secret: The client secret for the OAuth2 authentication
        :param grant_type: The grant type for the OAuth2 authentication (default is "password")
        :param scope: The scope of the access request (optional)
        """
        url = f"{self.base_url}/api/v1/auth/token"

        payload = {
            "grant_type": grant_type,
            "username": username,
            "password": password,
            "client_id": client_id,
            "client_secret": client_secret,
            "scope": scope,
        }

        response = requests.post(
            url,
            data=payload,
            verify=self.verify_ssl,  # Use global SSL verification setting
        )

        if response.status_code == 200:
            self.token = response.json().get("access_token")
        else:
            raise Exception(f"Failed to authenticate: {response.json()}")

    def get_headers(self) -> Dict[str, str]:
        """
        Constructs headers for authenticated requests.

        :return: A dictionary with authorization headers
        """
        if not self.token:
            raise Exception("Authentication required. Call authenticate() first.")
        return {"Authorization": f"Bearer {self.token}"}

    def list_samples(self) -> Any:
        """
        Lists all available samples.

        :return: The response from the samples listing API
        """
        url = f"{self.base_url}/api/v1/samples/"
        response = requests.get(url, headers=self.get_headers(), verify=self.verify_ssl)
        return response.json()

    def upload_sample(
        self, file_path: str, sample_name: str, disclaimer_confirmed: bool
    ) -> Any:
        """
        Uploads a sample BED file to the server.

        Args:
            file_path (str): The path to the local BED file to upload.
            sample_name (str): The name of the sample.
            disclaimer_confirmed (bool): Whether the disclaimer is confirmed.

        Returns:
            The response from the server.
        """
        url = f"{self.base_url}/api/v1/samples/"

        # Log file details
        logging.info(f"Preparing to upload file: {file_path}")
        try:
            file_size = os.path.getsize(file_path)
            logging.info(f"File size: {file_size} bytes")

            # Check first few lines of file
            with open(file_path, "r") as f:
                first_lines = [next(f) for _ in range(5)]
                logging.info("First 5 lines of file:")
                for line in first_lines:
                    logging.info(line.strip())
        except Exception as e:
            logging.error(f"Error checking file: {str(e)}")

        # Prepare the query parameters
        params = {
            "sample_name": sample_name,
            "disclaimer_confirmed": str(disclaimer_confirmed).lower(),
        }

        # Open the file in binary mode
        with open(file_path, "rb") as file:
            files = {"file": (file_path, file)}

            # Send the POST request with the file and query parameters
            logging.info(f"Sending PUT request to {url}")
            logging.info(f"Parameters: {params}")

            response = requests.put(
                url,
                headers=self.get_headers(),
                files=files,
                params=params,
                verify=self.verify_ssl,
            )

            logging.info(f"Response status code: {response.status_code}")
            logging.info(f"Response text: {response.text}")

            if response.status_code == 200:
                return response.json()
            else:
                raise Exception(
                    f"Failed to upload file: {response.status_code}, {response.text}"
                )

    def get_sample(self, sample_id: int) -> Any:
        """
        Fetches details for a specific sample.

        :param sample_id: The ID of the sample
        :return: The response from the sample retrieval API
        """
        url = f"{self.base_url}/api/v1/samples/{sample_id}"
        response = requests.get(url, headers=self.get_headers(), verify=self.verify_ssl)
        return response.json()

    def get_sample_report(self, sample_id: int) -> Any:
        """
        Fetches a sample report for a specific sample if it exists.

        :param sample_id: The ID of the sample
        :return: The response content (could be a file, JSON, or plain text)
        """
        url = f"{self.base_url}/api/v1/samples/download_sample_result/{sample_id}"

        # Perform the GET request
        response = requests.get(url, headers=self.get_headers(), verify=self.verify_ssl)

        # Check the response status
        if response.status_code == 200:
            # Determine the content type
            content_type = response.headers.get("Content-Type")

            if "application/json" in content_type:
                # If the response is JSON, return it
                return response.json()
            elif (
                "application/pdf" in content_type
                or "application/octet-stream" in content_type
            ):
                # If the response is a file (e.g., PDF or binary data), return the raw content
                return response.content
            else:
                # Fallback for text responses
                return response.text
        else:
            raise Exception(
                f"Failed to fetch the report: {response.status_code}, {response.text}"
            )

    def delete_sample(self, sample_id: int) -> Any:
        """
        Deletes a specific sample.

        :param sample_id: The ID of the sample to delete
        :return: The response from the sample deletion API
        """
        url = f"{self.base_url}/api/v1/samples/{sample_id}"
        response = requests.delete(
            url, headers=self.get_headers(), verify=self.verify_ssl
        )
        return response.json()

    def process_streaming(self, file1: str, file2: str, out_file: str) -> None:
        """
        One‐pass merge‐join of sorted inputs.  If the inputs aren't already
        sorted, we sort them in‐memory first (O(N log N) time, O(N) memory).
        Windows with no matching CpGs in file2 are skipped entirely.
        """

        # --- helper to sort chromosome positions lexicographically ---
        def chr_pos_key(x):
            # x[0] is chr (str), x[1] is start (int)
            return (x[0], x[1])

        # --- load and sort file1 windows ---
        windows: List[Tuple[str, int, int, str]] = []
        with open(file1) as f1:
            next(f1)  # skip header
            for line in f1:
                chr_, s, e, _, _, id_ = line.rstrip("\n").split()
                windows.append((chr_, int(s), int(e), id_))
        windows.sort(key=chr_pos_key)

        # --- load and sort file2 CpGs ---
        cpgs: List[Tuple[str, int, int, float, float]] = []
        with open(file2) as f2:
            for line in f2:
                cols = line.rstrip("\n").split("\t")
                cpgs.append(
                    (
                        cols[0],  # chr
                        int(cols[1]),  # start
                        int(cols[2]),  # end
                        float(cols[9]),  # coverage
                        float(cols[10]),  # methylation fraction
                    )
                )
        cpgs.sort(key=chr_pos_key)

        # --- now do the one-pass merge ---
        buffer = deque()
        idx2 = 0
        n2 = len(cpgs)

        with open(out_file, "w") as out:
            out.write("chr start end methylation_percentage coverage IlmnID\n")

            for chr1, w0, w1, id1 in windows:
                # advance idx2 until cpgs[idx2].start ≥ window start
                while idx2 < n2 and (cpgs[idx2][0], cpgs[idx2][1]) < (chr1, w0):
                    idx2 += 1

                # buffer all CpGs with start ≤ window end
                j = idx2
                while j < n2 and (cpgs[j][0], cpgs[j][1]) <= (chr1, w1):
                    buffer.append(cpgs[j])
                    j += 1
                idx2 = j

                # drop any CpGs not fully inside [w0, w1] or wrong chr
                while buffer and (
                    buffer[0][0] != chr1  # different chromosome
                    or buffer[0][1] < w0  # start before window
                    or buffer[0][2] > w1  # end after window
                ):
                    buffer.popleft()

                # compute coverage & methylation
                tot_cov = 0.0
                tot_mr = 0.0
                for _, _, _, cov2, pct2 in buffer:
                    tot_cov += cov2
                    tot_mr += cov2 * pct2 * 0.01

                # skip windows with no coverage
                if tot_cov == 0.0:
                    continue

                meth_pct = tot_mr / tot_cov * 100.0
                out.write(f"{chr1} {w0} {w1} {meth_pct:.2f} {int(tot_cov)} {id1}\n")
