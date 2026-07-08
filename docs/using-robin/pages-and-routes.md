## Sample ID generator {#sample-id-generator}

**What it is:** A form to **register a sample identifier** that ROBIN will later link to a sequencing run. You can:

- **Generate an MD5 run name** from **Test ID** (required) plus optional first name, last name, and date of birth, **or**  
- **Use your own MinKNOW run ID** (custom identifier)

Optionally store **first name**, **last name**, **hospital number**, and free-text **notes** encrypted alongside the sample. If you enter any of those fields, **date of birth is required** (it is the unlock key for the encrypted data).

**Use the same value in MinKNOW:** After you **Register sample ID**, copy the resulting **Sample ID** and enter it as the **sample ID** for your sequencing run in **MinKNOW** (or use the same registration controls on the Sequencer page). ROBIN creates a folder named with that ID under the work directory and writes a small manifest there. When BAMs and run folders appear under the same name, ROBIN **matches** the run to the registration **automatically**. Do not invent a different ID in MinKNOW after registering—use the registered ID verbatim (whether MD5 or custom).

The browser title may read **Generate Sample Identifier**.

![Sample ID generator: fields, generated MD5 identifier, Generate and Copy](../images/GenerateSampleID.png)

**What you’ll see:**

- Toggle: **Use my sample ID** (default) or **Generate MD5 ID**. Switching modes clears the registered Sample ID field so you don’t copy the wrong value.  
- For custom: **MinKNOW RUN ID** (required) and optional **Test ID**. For MD5: **Test ID** (required).  
- Optional encrypted fields: first name, last name, date of birth, hospital number, notes.  
- **Register sample ID** — creates/uses the public sample ID and writes the manifest when a work directory is configured.  
- **Copy to clipboard** — copies the registered ID when one exists.

**How it works:**

1. **Public sample ID / MinKNOW run ID** — Either the **MD5** digest of `Test ID|First name|Last name|DOB` (pipe-separated; blank optionals become empty strings; DOB as **YYYY-MM-DD**), or the **MinKNOW RUN ID** you enter. This string is used for the Nanopore run itself and is the folder name ROBIN links identifier data to.  
2. **Optional encryption** — If you provide first name, last name, hospital number, or notes, DOB is mandatory. Those fields (and DOB) are stored **encrypted** with a key derived from DOB. Test ID remains plaintext in the manifest when provided.  
3. **ID-only registration** — You can register just a sample ID (MD5 or custom) with no encrypted PII; ROBIN still creates the folder/manifest so the run can be linked when detected.  
4. **After registering** — If a **work directory** is configured, ROBIN writes `sample_identifier_manifest.json` under `{work directory}/{sample ID}/`. If no work directory is set, you still get the ID in the UI, but the persist step cannot run—watch the on-screen notification.

The same registration controls (including encrypted notes) are available on the **Sequencer (MinKNOW)** page under **Register sample identifiers**, so you do not need to visit this page before starting a run.

Use this when your lab needs either a deterministic MD5 scheme or its own naming convention before or during a run.

*Bookmark path: `/sample_id_generator`.*

---

## Next

- [Reading your results](sample-results.md) — understand each block on the sample page.  
- [Troubleshooting](troubleshooting.md) — if something doesn’t look right.
