# Panel commands

Target **panels** define which genes/regions drive **target**, **CNV**, and **fusion** analyses and which packaged BED files ROBIN uses. Built-in panels ship with the install; you can add **custom** panels from BED files.

Commands:

- [`robin list-panels`](#robin-list-panels)
- [`robin add-panel`](#robin-add-panel)
- [`robin remove-panel`](#robin-remove-panel)

All require the interactive **disclaimer** (`I agree`).

---

## `robin list-panels`

Lists **built-in** panels (`rCNS2`, `AML`, …) and **custom** panels registered in `robin.resources`.

```bash
robin list-panels
```

Use the printed names with **`robin workflow --target-panel`** and **`robin utils sequencing-files --panel`**.

---

## `robin add-panel`

Register a new panel from a BED file.

```bash
robin add-panel /path/to/panel.bed MyPanelName
```

### BED format

- At least **4 columns**: `chrom`, `start`, `end`, `gene` (or comma-separated genes).
- **6-column** BED is allowed (`score`, `strand`).
- Regions spanning multiple genes can list **comma-separated** gene symbols in the name column.

### Options

| Flag | Meaning |
|------|---------|
| `--validate-only` | Check the BED only; do not copy files into resources. |

### Reserved names

You **cannot** use `rCNS2` or `AML` as custom names.

### What ROBIN stores

- **Processed** unique-gene BED: `{PanelName}_panel_name_uniq.bed` under `robin.resources`.
- **Original upload** copy: `{PanelName}_panel_source.bed` (for traceability and [`sequencing-files`](utils.md#robin-utils-sequencing-files)).

If a panel with the same name already exists, the command fails until you [`remove-panel`](#robin-remove-panel).

---

## `robin remove-panel`

Remove a **custom** panel (not built-ins).

```bash
robin remove-panel MyPanelName
robin remove-panel MyPanelName --force
```

| Flag | Meaning |
|------|---------|
| `-f` / `--force` | Skip the `yes` confirmation (destructive). |

Built-in panels **cannot** be removed.

---

## Using panels in a workflow

```bash
robin workflow /data/bams \
  -w target,cnv,fusion \
  --target-panel MyPanelName \
  --center Sherwood \
  -d ~/results \
  --reference ~/references/hg38.fa
```

---

## Related

- [Utilities — `sequencing-files`](utils.md#robin-utils-sequencing-files) (copy panel BED + reference for MinKNOW)  
- [`robin workflow`](workflow.md)  
