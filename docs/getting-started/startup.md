# What happens when you start ROBIN

!!! abstract "What this page covers"
    The usual order of checks and prompts when you run **`robin workflow`**: models, reference, **research disclaimer** (`I agree`), optional **large-BAM** warning, **Ray vs threads**, **NiceGUI** URL, **GUI password**, then **watching BAMs**.  
    Install: [Installation](installation.md). Flags: [`robin workflow`](../cli/workflow.md). Always confirm behaviour with **`robin workflow --help`** on your install.

---

## At a glance (startup order)

| Step | What happens |
|------|----------------|
| 1 | **Model assets** checked — exit with instructions if missing (`robin utils update-models`). |
| 2 | **Reference** validated if you pass `--reference` / `-r`. |
| 3 | **Disclaimer** — you must type **`I agree`** exactly (case-sensitive). |
| 4 | **Large BAMs** — optional warning if `ROBIN_PROCESS_LARGE_BAMS` is set (not for live runs). |
| 5 | **Configuration summary** printed (paths, `--center`, steps, logging, Ray, etc.). |
| 6 | **Execution engine** — Ray (default) or `--no-use-ray` (threaded). |
| 7 | **NiceGUI** — if `--with-gui` (default); on **Ray**, GUI usually needs **`--work-dir`** — see [NiceGUI](#nicegui-workflow-monitor-default-on). |
| 8 | **GUI password** — set or verify (terminal prompts). |
| 9 | **Watch** input directory for `*.bam` and schedule jobs. |

---

## 1. Model assets

ROBIN checks **required model files** (same manifest as `robin utils update-models`). If anything is missing, the process exits and tells you to run **`robin utils update-models`** (and **`GITHUB_TOKEN`** if assets are on private GitHub).

---

## 2. Reference genome (`--reference` / `-r`)

If you pass a reference **FASTA**, ROBIN validates it and ensures an index (e.g. `.fai`) exists or can be created. On failure, fix the path or omit `--reference` only if your workflow truly does not need it.

---

## 3. Research disclaimer (`I agree`)

ROBIN prints the **research-use** text and waits for **`I agree`** exactly (case-sensitive) before the workflow starts.

---

## 4. Optional warning: large BAMs

If **`ROBIN_PROCESS_LARGE_BAMS`** is enabled, ROBIN warns that this mode must **not** be used with **live** sequencing runs.

---

## 5. Configuration summary

The terminal prints paths, **`--center`**, workflow steps (including automatic **`preprocessing`** / **`bed_conversion`** insertion), logging, Ray/threading mode, and related options.

---

## 6. Execution engine (Ray vs threading)

| Mode | Flag | Notes |
|------|------|--------|
| **Ray** (default) | `--use-ray` | Ray Core driver; optional dashboard and CPU limits. |
| **Threading** | `--no-use-ray` | Threaded workers instead of Ray. |

---

## 7. NiceGUI workflow monitor (default: on) {#nicegui-workflow-monitor-default-on}

With **`--with-gui`** (default), ROBIN starts a **browser** workflow monitor.

### When the GUI starts

| Mode | Behaviour |
|------|------------|
| **Ray (default)** | GUI starts **inside** the Ray driver **only if** you pass **`--work-dir`** (`-d`). Without `--work-dir`, the driver may skip the GUI and report that **`--work-dir` was not provided**. |
| **`--no-use-ray`** | GUI starts when **`--with-gui`**, using **`--work-dir`** if set, else the **watched BAM directory** as context. |

**Practical takeaway:** on the default **Ray** path, pass **both** the data directory (positional `PATH`) **and** **`--work-dir`** if you want the browser UI.

### URLs

The terminal prints a base URL such as **`http://<gui-host>:<gui-port>`** (defaults: **`0.0.0.0`**, **`8081`**). Common routes include:

- Welcome: `/`  
- Workflow monitor: `/robin`  
- Sample views: under `/live_data` (exact routes are printed at startup)

Use **`--gui-host`** / **`--gui-port`** to change bind address and port; **`--no-gui`** disables the web UI.

```bash
robin workflow ... --no-gui
```

---

## 8. GUI password (terminal)

Access to the web UI is password-protected (**`argon2-cffi`**). Passwords are **not** echoed.

### First run (no password yet)

If stdin is a **TTY**, you are prompted:

```text
Set GUI password:
Confirm GUI password:
```

Enter the same password twice. If **no TTY** (some automation), startup **fails** with a message to run from an interactive terminal first.

### Later runs

If stdin is a **TTY**:

```text
GUI password:
```

Wrong password → **Invalid password.** and the GUI does not start.

### Change password without a full workflow

```bash
robin password set
```

See [GUI password](../cli/password.md).

---

## 9. Workflow hooks (optional)

If the GUI starts, ROBIN may install **workflow hooks** for live updates. If hook installation fails, you may see a message that the GUI will show **static** information only.

---

## 10. Watching BAMs and shutdown

The runner **watches** the input directory for **`*.bam`** files (subject to `--no-process-existing`, `--no-watch`, etc.) and schedules jobs. **Ctrl+C** stops the run; ROBIN attempts **graceful shutdown** (workflow manager, Ray, GUI), though exit may take a moment.

---

## Quick reference — prompts

1. Model check  
2. (If `-r`) Reference validation  
3. Type **`I agree`**  
4. (If GUI) **GUI password** — set twice first time, or single verify later  
5. Open the printed URL in a browser  

---

## See also

- [Quickstart](quickstart.md)  
- [`robin workflow`](../cli/workflow.md)  
- [GUI password](../cli/password.md)  
- [CLI overview](../cli/index.md)  
