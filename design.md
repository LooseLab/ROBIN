Design System Specification: Editorial Bioinformatics

1. Creative North Star: The Digital Curator

The bioinformatics landscape is often cluttered with dense, overwhelming data. This system acts as a "Digital Curator," transforming complex genomic output into high-fidelity, actionable insights. We prioritize clear visual hierarchy, data density without clutter, and semantic color usage.



2. Core Visual Principles

2.1 Information Architecture (The "Hybrid" Layout)





The Anchor (Header): Large, bold typography for the Sample ID (e.g., "Sample_104"). Primary actions like "Generate Report" must be high-contrast and easily accessible.



The Summary Strip: Technical metadata (Run Time, Device, Flow Cell) is presented in a compact, horizontal strip using secondary typography and supporting iconography.



Modular Analytical Cards: Distinct sections for Classification, Coverage, and Analysis are housed in cards with subtle borders or depth to maintain separation.

2.2 Typography





Headlines: Manrope (Bold/ExtraBold) for primary anchors.



Data Labels: Inter (Medium/SemiBold) for high legibility at small sizes.



Monospace: JetBrains Mono for technical strings, DNA sequences, or precise numeric data.



3. Thematic Depth & Gradients (The "Editorial" Look)

To achieve a professional, high-fidelity aesthetic, surfaces should not be flat.





Subtle Gradients: Use light vertical gradients (e.g., from-slate-50 to-white in Light Mode, or from-slate-900 to-zinc-950 in Dark Mode) to create structure.



Accent Glows: Use low-opacity radial gradients behind primary actions or critical data points to indicate focus and importance.



Borders: Keep borders minimal (1px) and semi-transparent to define boundaries without adding visual weight.



4. Shell & Navigation Theming (Header & Footer)

The **app shell** (header and footer) uses a **continuous emerald brand fill** so navigation and system status stay recognizable and aligned with the primary Emerald 600 / 400 semantic palette. Main content areas remain on white / slate gradients (§6); the shell is intentionally distinct.

4.1 The Header (The "Anchor")

The header serves as the primary navigation and status center.

**Light Mode:** Horizontal gradient **left → right**: darker **Emerald 700** (#15803D) on the left, lighter **Emerald 600** (#16A34A) on the right, with a soft **1px** light edge (semi-transparent white) along the bottom. **Text, icons, and meters** use **white** (#FFFFFF) for contrast; optional subtle text shadow for legibility.

**Dark Mode:** Same geometry: **Emerald 800** (#065f46) → **Emerald 700** (#047857) **left → right**. Text and icons remain **white**; primary accents on controls may use **Emerald 400** where a luminous cue is needed.

4.2 The Footer (Metadata & System Status)

The footer uses the same two emerald stops, with a **horizontal gradient right → left** (darker on the right, lighter on the left) so it mirrors the header. It stays minimal: secondary system information in **Inter** at small/medium sizes, **white** on the green background. Links and secondary actions use light-on-emerald contrast (e.g. white / near-white, or translucent white fills on hover).





5. Data Visualization Style Guide (CRITICAL)

5.1 Color Semantics & Meaning

Visualizations must use color to convey status and biological significance consistently.





Positive Signals / Genomic Gains / Success:





Light Mode: Emerald 600 (#16A34A)



Dark Mode: Emerald 400 (#34D399) with a subtle "neon" glow.



Negative Signals / Genomic Losses / Deviations:





Light Mode: Red 600 (#DC2626)



Dark Mode: Coral Red / Rose 400 (#FB7185)



Warning / Medium Confidence:





Light Mode: Amber 500 (#F59E0B)



Dark Mode: Amber 300 (#FCD34D)

5.2 Genomic Mapping (CNV Plots)





Segmented Chromosomes: All genomic plots (1-22, X, Y) must be explicitly segmented.





Vertical Grids:* Use subtle vertical lines to delineate chromosome boundaries.



Alternating Stripes: Use very light/dark alternating background shades for odd/even chromosomes to aid horizontal tracking across the genome.



Baseline Alignment: A distinct horizontal baseline (typically Y=0) must be visible to anchor gains (above) and losses (below).



6. Theme Implementation Summary

6.1 Light Mode (Professional Lab)





Background: White (#FFFFFF) / Slate 50 (#F8FAFC).



Cards: White with a 1px Slate 200 border.



Text: Slate 900 for headlines, Slate 600 for secondary text.

6.2 Dark Mode (Midnight Analysis)





Background: Zinc 950 (#09090B) / Slate 950 (#020617).



Cards: Semi-transparent Slate 900 (rgba(15, 23, 42, 0.5)) with a 1px Slate 800 border.



Luminance: Use high-saturation versions of semantic colors to "pop" against the dark background.



7. Run Summary Section (Web & Mobile)

To ensure consistent implementation across web and mobile, use the following styling for **Run Summary** (and similar high-density metadata blocks). Layout is **not** restricted to a 2×2 grid: choose **2×2**, a **single horizontal strip** (e.g. 4×1), **3×2**, or other arrangements that fit the viewport while preserving legibility and touch targets.

7.1 Structural Layout (Responsive Logic)

**Desktop / web:** Prefer a **high-density** layout—e.g. a 2×2 grid **or** a 4×1 horizontal strip—to **maximize vertical space** for analysis plots below. Other grid shapes are allowed when they improve scanability or fit the data (e.g. long technical strings on a full-width row).

**Mobile:** Default to a **single-column vertical stack**. A **2×2** (or similar) grid is acceptable **when space permits**, provided touch targets stay large enough and labels are not truncated.

**Padding:** Use at least **16px** (`p-4`) padding **inside each cell** (or the equivalent inner padding for the container when using a unified card). For **dense** layouts, **12px–16px** horizontal padding per cell (`p-3`–`p-4` or `px-4`) is acceptable if touch targets and scanability remain good.

**Outer container vs grid (full width):** Apply **horizontal padding on the outer Run Summary shell** (e.g. ~12–16px left/right). The **inner CSS grid should be `width: 100%`** within that padded area. **Do not** combine `width: 100%` on the grid with **additional horizontal margins** on the same element—this overflows the parent and clips the right-hand gutter when the parent uses `overflow: hidden`. Use **vertical margin only** on the grid if needed (e.g. margin-top below the section title).

**Section title:** The title bar may use a **full-bleed bottom border** edge-to-edge on the shell; restore text alignment with the same horizontal padding as the grid (e.g. negative horizontal margin on the title row equal to the shell padding, with matching padding on the title text).

**Long values:** Fields such as **basecall model** should use a **full-width grid row** (span all columns) so long strings wrap cleanly. Shorter follow-up rows (e.g. sample ID) may span a subset of columns; define column spans explicitly in CSS so layout does not depend on utility classes that might be missing from a build.

**Flex wrapper:** The outer shell should behave as a **flex column** with **stretch** alignment so the grid always fills the available width (avoid flex children that shrink-wrap to content width).

7.2 Light Mode Styling (Clinical & Precise)

**Container:** `bg-white`, a very subtle border `border-slate-200`, and soft elevation `shadow-sm`.

**Dividers:** When using a multi-cell grid, use thin internal rules: `divide-slate-100` (or equivalent) between cells—e.g. a **1px gap** between cells with a neutral track color, or hairline borders.

**Iconography:** **Solid Emerald** `#10b981` for icons. Icon size **20px–24px** (use the upper end when values are set to **1rem** for balance).

**Typography:**

- **Label:** Small caps line: **`text-[10px]`** (minimum) or **`~11px`** for improved legibility; **`uppercase`** **`tracking-wider`** **`text-slate-500`** **`font-semibold`**
- **Value:** **`text-sm`** **`font-medium`** **`text-slate-900`** as the baseline; **`text-base` (1rem)** is allowed for **higher scanability** on dense dashboards. Use **JetBrains Mono** (or the app monospace) for technical values.

7.3 Dark Mode Styling (High-Tech & Reduced Strain)

**Container:** Semi-transparent surface, e.g. `bg-slate-800/40`, with a defined border `border-slate-700` for depth.

**Dividers:** `divide-slate-700/50` between grid cells (or the same 1px-gap “track” approach as light mode with a darker track color).

**Iconography:** **Neon green** `#22c55e` for icons. Optional subtle glow: `drop-shadow-[0_0_5px_rgba(34,197,94,0.3)]`.

**Typography:**

- **Label:** Same treatment as §7.2 with **`text-slate-400`** in dark mode.
- **Value:** **`text-sm`** **`font-medium`** **`text-white`** as the baseline; **`text-base` (1rem)** is allowed to match light-mode density. Use monospace for technical strings.

7.4 Component Hierarchy (Internal Layout)

Within **each** grid cell, use this structure:

1. **Icon:** **Left of** the text stack (recommended) to save vertical space, or **above** the text (centered) when a hero treatment is desired.
2. **Text stack:** A vertical column: **Label** (top) → **Value** (bottom). Optional third lines (hints) should be avoided in dense layouts; use **tooltips** or the **title** attribute when extra context is required without increasing row height.
3. **Spacing:** **12px** (`gap-3`) between the icon and the text stack.

7.5 Implementation Note (Tailwind / NiceGUI / Quasar)

**Dark mode:** This app uses **Quasar**’s **`body--dark`** class on `<body>`. **Tailwind `dark:` variants do not apply** unless the project also sets a `dark` class on a parent (e.g. `html`). For Run Summary, **light and dark tokens must be applied with explicit `.body--dark` rules in CSS** (or equivalent), not only with `dark:` utilities.

**Shell pattern:** Prefer a **block or flex column** wrapper with semantic classes (e.g. `run-summary-shell`) so the inner grid is not sized incorrectly by framework column helpers. Pair **§7.2** surfaces with default/light selectors and **§7.3** surfaces with `.body--dark` selectors.

**Example (dark surface only—mirror tokens for light per §7.2):**

`ui.element('div').classes('run-summary-shell')` with styles such as `border border-slate-700 rounded-xl bg-slate-800/40` under `.body--dark`, and light counterparts without `.body--dark`.

7.6 Reference Implementation (ROBIN)

The **Run Summary** on the sample detail page uses **dedicated CSS** in `styles.css` (classes `run-summary-shell`, `run-summary-grid`, `run-summary-cell`, etc.) so **layout and colors stay correct in dark mode** without relying on Tailwind `dark:`. The grid is **responsive** (1 → 2 → 4 columns by breakpoint), with **explicit `grid-column: span …`** rules for full-width rows. **Hints** are implemented as **tooltips / `title` on values** to keep row height compact.



8. Classification Insight Cards (Diagnostic Pulse)

These cards summarize each classifier on the sample detail page and should stay visually aligned with **Run Summary** (§7): same **`.body--dark`** approach in CSS, not Tailwind `dark:` alone.

8.1 Structural Layout & Hierarchy

**Grid:** Responsive **1 column (mobile) → 2 columns (tablet, from ~640px) → 4 columns (desktop, from ~1280px)**.

**Internal padding:** **`p-4` (16px) on mobile**, **`p-6` (24px) on `md` and up**, applied inside each card.

**Alignment:** **Vertically stack** content: model header (icon + name), main call, confidence bar, numeric meta, confidence label line, footer description.

**Interactivity:** The **entire card** is the hit target (`cursor-pointer`); **click** scrolls to the matching **detailed classification expansion** (anchor id `classification-detail-{tool_step}`).

8.2 Visual Encoding (Confidence Bar)

Pair the **numeric confidence** with a **thin horizontal bar** (height **~6–8px**, **fully rounded** ends) for fast scanning.

**Fill color by score (same band logic as legacy badges):**

- **High (≥80%):** Emerald light **`#10b981`** / neon green dark **`#22c55e`**
- **Medium (50–79%):** Amber **`#f59e0b`**
- **Low (<50%):** Rose **`#f43f5e`**

**Track (unfilled):** Light **`bg-slate-100`**; dark **`bg-slate-700`**.

8.3 Light Mode (Clean & Academic)

**Card:** **`bg-white`**, border **`border-slate-200`**, **`shadow-sm`**; on hover **`shadow-md`** (and optional slight elevation change).

**Model name (header):** **`text-[10px]`** **`uppercase`** **`tracking-widest`** **`text-slate-500`** **`font-bold`**, with bottom spacing (`mb-1` equivalent).

**Main result:** **`text-base`** **`font-semibold`** **`text-slate-900`**.

8.4 Dark Mode (Focused & Modern)

**Card:** **`bg-slate-800/60`**, border **`border-slate-700/50`**; on hover **`border-green-500/30`** (or similar luminous edge).

**Model name:** **`text-slate-400`** (same uppercase / tracking as light).

**Main result:** **`text-base`** **`font-semibold`** **`text-white`**.

8.5 Mobile

**Truncation:** On **narrow viewports**, the **main result** line should use **ellipsis** (`truncate` / `text-overflow: ellipsis` + `nowrap`) so long biological labels do not break card height; provide full text via **`title`** on hover/tap. From **`md` up**, allow **normal wrapping** again.

8.6 Implementation Note (NiceGUI)

Use a **wrapper** (`ui.element('div')` or `ui.card`) with classes such as **`hover:shadow-md transition-all cursor-pointer`**, and **`on('click', …)`** to scroll. Prefer **dedicated CSS classes** in `styles.css` with **`.body--dark`** for dark tokens, consistent with §7.5.

8.7 Reference Implementation (ROBIN)

The sample **Classification details** strip uses **`classification-insight-*`** classes in `styles.css` and **`_create_classification_dashboard_card_with_data`** in `summary.py**. Expansion panels in `classification.py` expose **`id="classification-detail-{tool_step}"`** for scroll targets.


8.8 Classification Charts (ECharts — Sturgeon, NanoDX, PanNanoDX, Random Forest)

The **Classification** section shell (**`classification-insight-shell`**, **`id=classification-section`**) wraps per-classifier **expansion** panels. Inside each panel, **summary** rows use **`classification-insight-card`** + **`classification-insight-icon`** / **`classification-insight-result`** / **`classification-insight-meta`**, aligned with §8.1–§8.4. **Charts** are **Apache ECharts** via NiceGUI **`ui.echart`** in **`classification.py`**.

**Horizontal bar (Top classes):**

- **Signal-first colour:** The **highest-scoring** class (leading bar) uses **brand green** (**emerald `#10b981`**); all other bars use **muted slate** (**`#94a3b8`**). No rainbow per-class colours.
- **Thresholds:** Two **vertical** reference lines on the **value (X) axis** at the classifier’s **`medium`** and **`high`** bands from **`CLASSIFIER_CONFIDENCE_THRESHOLDS`** in **`src/robin/classification_config.py`** (e.g. Sturgeon **80% / 95%**, NanoDX **25% / 50%**, Random Forest **65% / 85%**). Labels read **Medium (…%)** and **High (…%)** — diagnostic zones, not a dense tick grid. If a key is missing, **`DEFAULT_CONFIDENCE_THRESHOLDS`** is used.
- **Bar geometry:** **Thicker** bars (**fixed pixel width**), **rounded** ends (**pill** aesthetic, **`borderRadius`** on bar items) to align visually with Run Summary / Insight Cards.

**Time series (Confidence over time):**

- **Leader vs secondary:** The **current top class** (first of the top five series) is drawn **last** with a **thick** **neon green** stroke (**`#22c55e`**) and a **soft glow** (**`lineStyle.shadowBlur` / `shadowColor`**). Other candidate classes use **thin dashed** **slate grey** lines (**`#94a3b8`**).
- **Confidence zones:** **Shaded bands** on the **Y axis**, aligned with the same **`medium`** / **`high`** values: **medium–high** uses a **subtle slate** fill; **high–100%** uses a **subtle emerald** tint so the top band matches **`get_confidence_level`** semantics in **`classification_config.py`**.
- **Threshold lines:** **Horizontal** dashed lines at **`medium`** and **`high`** with **Medium (…%)** / **High (…%)** labelling, matching the bar chart.
- **X-axis:** **Time** axis uses **reduced tick density** (**`splitNumber`**, **`axisLabel.hideOverlap`**) and, when data exist, **min/max** from the series span so the axis is not an overly dense timestamp grid. Prefer a **small number of readable** ticks over every raw timestamp.

**Legend:** **Scrollable** legend. **Desktop / wide:** **vertical**, **top-right** (**`grid.right`** reserves ~20% so the plot stays readable). **Mobile / narrow:** **`media`** queries move the legend to a **horizontal strip under the plot** (full chart width, scrollable) so the **grid** can use almost the full width — see §8.8.1. Each line series sets **`color`** and **`itemStyle.color`** to the same value as **`lineStyle.color`** so **legend swatches match the stroke** (leader = neon green, others = muted slate); ECharts otherwise applies the default palette to the legend only.

**8.8.1 Mobile and narrow viewports**

- **Container height:** Charts use **taller min-heights on small screens** (Tailwind **`min-h-*` / `h-*`** on the **`ui.echart`** element) so bars and lines are not squeezed into an ultra-wide aspect ratio on phones.
- **ECharts `media`:** **`_echart_media_responsive_bar()`** and **`_echart_media_responsive_ts()`** in **`classification.py`** adjust **`grid`**, **title** font size, **category label** font size / **truncate** width, and (for the time series) **legend** position. Breakpoints use **`maxWidth: 640`** and **`maxWidth: 400`** (container pixels) so layout updates when the user rotates the device or resizes the window.
- **Time series:** Legend **below** the plot on narrow widths avoids compressing the **x** (time) axis; keep **`type: 'scroll'`** on the legend so many classes remain usable.

**Tooltips:** **Insight-style** chrome: **dark** background (**`rgba(15,23,42,0.96)`**), **light** text, **border**. **Bar** chart tooltip uses an explicit **“Confidence: …%”** line; **time-series** tooltips rely on the **default axis** content (series name + numeric value with **Y** in **%**) under the same dark styling. **Axis** pointer line stays subtle (**slate**).

**Implementation note:** Plot **reference lines** and **shaded bands** are driven by **`_confidence_thresholds_for_classifier()`** in **`classification.py`**, which reads **`CLASSIFIER_CONFIDENCE_THRESHOLDS`** in **`src/robin/classification_config.py`** (with **`DEFAULT_CONFIDENCE_THRESHOLDS`** as fallback) so charts align with **`get_confidence_level`**, reporting, and other consumers. The **Classification details** insight cards (§8.2) use their own **confidence bar** tier logic for the strip; **ECharts** use the config’s **`medium`** and **`high`** per classifier key (**`sturgeon`**, **`nanodx`**, **`pannanodx`**, **`random_forest`**).


9. Analysis Details Strip (Coverage, CNV, MGMT, Fusion)

This strip sits below **Classification details** on the sample detail page. It should stay visually aligned with **Run Summary** (§7) and **Classification Insight Cards** (§8): same **responsive grid** (`classification-insight-grid`), same **card shell** (`classification-insight-card`), **emerald** header icons (`classification-insight-icon`), and **explicit `.body--dark` rules** in CSS—not Tailwind `dark:` alone.

9.1 Layout & Hierarchy

**Shell:** Reuse **`classification-insight-shell`** and **`classification-insight-heading`** for the outer block and section title (**“Analysis details”**).

**Grid:** Same breakpoints as §8: **1 → 2 → 4** columns (mobile → tablet → desktop).

**Cards:** Same hierarchy as §8 for **header** (icon + short label) and **primary line** (`classification-insight-result`), then **meta** (`classification-insight-meta`), optional **tier-coloured labels** (`classification-insight-level--*`), **footer** (`classification-insight-foot`). **Unlike §8 classification cards, analysis cards do not use progress / confidence bars**—encode depth, CNV counts, and methylation with **text, pills, and coloured level lines** only.

**Interactivity:** The **entire card** is clickable (`cursor-pointer`); **click** runs **`_scroll_to_analysis_detail`** in `summary.py`, which scrolls to **`id="analysis-detail-{key}"`**. Keys used: **`coverage`**, **`cnv`**, **`mgmt`**, **`fusion`**.

9.2 Key metrics (no bars)

- **Coverage:** **Target depth** on the analysis strip uses **`classification-insight-level--*`** from bands: **≥30×** high, **≥20×** medium, **below 20×** low, with **threshold pills** (§9.3). The **expandable Coverage** block on the sample page adds a separate **Sufficient / Moderate / Low** badge on target estimated depth — see **§9.5** (≥30× / ≥15× / &lt;15×).
- **CNV:** **Gained vs total** and bin/variance as **meta**; **gained / lost** as **pills** (`analysis-insight-pill--*`).
- **MGMT:** **Methylation %** line uses **`classification-insight-level--*`** with the same **>10% / >5% / else** band logic as the legacy badges.

9.3 Threshold Pills (Coverage)

The **≥30× / ≥20× / ≥10× / below 10×** legend uses **`analysis-insight-pill`** and **`analysis-insight-pill--emerald|sky|amber|rose`** in `styles.css` so light and **`.body--dark`** backgrounds stay readable.

9.4 Reference implementation (ROBIN)

- **`summary.py`:** **`_analysis_section`** (shell + grid), **`_scroll_to_analysis_detail`**, **`_create_*_dashboard_card_with_data`** for each analysis card.
- **`styles.css`:** Reuses **`classification-insight-*`**; adds **`analysis-insight-pill*`** for the coverage legend; scoped **`#analysis-detail-cnv`** / **`#analysis-detail-mgmt`** / **`#analysis-detail-fusion`** for dark-mode copy on detail sections.
- **Scroll targets:** Outer section cards set **`id=analysis-detail-coverage`**, **`id=analysis-detail-cnv`**, **`id=analysis-detail-mgmt`**, **`id=analysis-detail-fusion`** in **`coverage.py`**, **`cnv.py`**, **`mgmt.py`**, **`fusion.py`** respectively.

**CNV detail UI & ECharts** — see **§9.6**. **Fusion detail** — see **§9.7**.

9.5 Coverage Plot Strategy

**Design intent**

1. **Signal-to-noise:** Avoid a spaghetti of per-target lines. Emphasise the **population mean** and **outliers**; dim targets that stay within a **normal band** (±2 SD of the mean coverage time series).
2. **On-target vs off-target:** **On-target** = brand green family; **off-target** = muted slate/grey so enrichment is obvious at a glance.
3. **Quality badge:** A single interpretive label (**Sufficient**, **Moderate**, **Low**) beside target estimated coverage, before the user reads charts. Thresholds in the UI: **≥30×** Sufficient, **≥15×** Moderate, **&lt;15×** Low (targets depth).

### A. Per Chromosome Target Coverage (grouped bar chart)

- **Layout:** Grouped bars — **On Target** vs **Off Target** per chromosome (`chr1`…`chr22`, `X`, `Y`). **X-axis** uses **short labels** (`1`, `2`, …, `X`, `Y`) to reduce overlap on narrow viewports.
- **Colour — dark mode:** On-target **`#22c55e`** (neon green); off-target **`#475569`** (slate-600).
- **Colour — light mode:** On-target **`#10b981`** (emerald); off-target **`#cbd5e1`** (slate-300).
- **Implementation:** **`echart_target_cov`** and **`_update_target_cov()`** in **`coverage.py`**; helpers **`_cov_on_off_colors()`**, **`_short_chr_label()`**, **`_cov_chrome_palette()`**.

### B. Coverage over time (accumulation)

- **Style:** Smooth **line + area**: cumulative **estimated coverage (×)** vs time; **semi-transparent fill** under the primary line (green tint).
- **Interaction:** **Crosshair** (`axisPointer: cross`) and axis tooltip. Data source: **`coverage_time_chart.npy`** (`[timestamp_ms, coverage]`); tooltip shows **estimated coverage (×)** at the crosshair (total reads are not stored in this file).
- **Implementation:** **`echart_time`**, **`_update_time()`** in **`coverage.py`**; **`_cov_line_primary()`**, **`_cov_area_fill_primary()`**.

### C. Individual target coverage (outlier analysis) — signal-first trend

Full visual and interaction spec: **Outlier Coverage & High-Density Trend Analysis** (see table at end of this subsection).

#### Design Specification: Outlier Coverage & High-Density Trend Analysis

##### 1. Plot archetype: the “signal-first” trend

In bioinformatics, outlier analysis often involves tracking multiple targets over time. To prevent visual “spaghetti,” the plot prioritises the **population mean** and **statistically significant deviations** while minimising background noise.

##### 2. Layering & signal-to-noise ratio (SNR)

**2.1 Background architecture**

- **The normal range (the “envelope”):** A shaded background band represents normal variance (**±2 standard deviations** from the **mean coverage time series** — mean of per-timepoint means).
  - *Light mode:* Slate 100 (`#F1F5F9`) at **50%** opacity.
  - *Dark mode:* Slate 800 (`#1E293B`) at **30%** opacity.
- **Grid lines:** Horizontal and vertical grids stay **minimal** and **semi-transparent** so they give scale without competing with data lines.

**2.2 The reference mean**

- **Style:** Thick **dashed** line for the population / run average.
- **Light mode:** Slate 800 (`#1E293B`) dashed stroke.
- **Dark mode:** Slate 200 (`#E2E8F0`) dashed stroke with a subtle **~1px glow** (`shadowBlur` on the line).

**2.3 Data line hierarchy**

- **Inactive / normal targets:** Non-outlier lines at **~15–20%** opacity (“ghost” population trend).
- **Active / outlier targets:** Genes crossing the threshold use **high-saturation** brand colours (**Emerald, Cyan, Rose, Amber** family — theme-tuned for light vs luminous dark). **Solid 2px** stroke.

##### 3. Visual encoding & semantic markers

**3.1 Point markers (anomalies)**

- **Gains / positive deviations:** Triangle-up (or circle) — **Emerald 400** (dark) / **Emerald 600** (light).
- **Losses / negative deviations:** Triangle-down — **Red 400** (dark) / **Red 600** (light).
- **Interaction:** Marker **emphasis** uses **shadow / glow** on hover (ECharts `emphasis.itemStyle`) so points read as interactive.

**3.2 Time-series smoothing**

- **Spline-style curves:** `smooth` + moderate tension (Bezier-style) on line series so raw sequencing jitter is softened without hiding breakpoints.

##### 4. Legend & interactivity

**4.1 Dynamic legend**

- **Placement:** **Right**, vertical **scroll** on wide viewports; **narrow viewports** fall back to **horizontal scroll** under the grid (media `option` in `coverage.py`, same idea as classification charts).
- **Focus behaviour:** **Legend hover** + ECharts **`emphasis.focus: 'series'`** / **`blur`** so one series is emphasised and others **dim** (~10% opacity). Native legend **click** still toggles series visibility (JSON-safe options; no custom JS formatters).

**4.2 Crosshair & scrubbing**

- **Vertical scrubber:** Axis pointer **`type: 'line'`** on **`x`** — thin full-height line following the pointer.
- **Tooltip:** **Axis** trigger (default) — shows values at the scrubbed time for visible series, including **Mean (population)** for context. Do **not** use JS **`formatter` functions** in options (NiceGUI serialises to JSON; use defaults / string templates only).

##### 5. Thematic implementation summary

| Element | Light mode (clinical) | Dark mode (midnight) |
| :--- | :--- | :--- |
| **Normal band** | Slate-100 envelope | Slate-800 envelope |
| **Mean line** | Dark slate, dashed | Slate-200 dashed + glow |
| **Outlier lines** | Solid high-contrast brand palette | Luminous emerald / cyan / rose / amber |
| **Anomalies** | Semantic red / green markers | Same hues, brighter + emphasis glow |

*Design spec version 1.0 | October 2026*

---

- **Implementation:** **`_plot_target_coverage_over_time()`** in **`coverage.py`** (`target_coverage_time.csv`); helpers **`_cov_band_fill_rgba()`**, **`_cov_mean_line_color()`**, **`_cov_outlier_line_color()`**, **`_cov_anomaly_*_color()`**, **`_cov_target_time_split()`**, **`_echart_target_cov_time_media()`**, **`_COV_MAX_DIM_GENES`**, **`_cov_legend_label`**. In-range suffix **`· in-range`**.
- **Tooltip:** Axis tooltip only — **no** JS `formatter` functions in options.

**Reference implementation:** **`coverage.py`** — module-level **`_cov_*`** theme helpers; charts updated when **`_refresh_coverage_sync`** runs (and bar colours follow storage **`dark_mode`** on refresh).


9.6 CNV detail section (Copy number)

The **expandable CNV block** (`#analysis-detail-cnv` in **`cnv.py`**) follows the same **Digital Curator** shell as **MGMT** and **Coverage**: **`classification-insight-shell`**, **`classification-insight-heading`** (title **without emoji** — e.g. **“Copy number (CNV)”**), a **`classification-insight-card`** summary row (**`person`** icon, status / genetic sex / bin width / variance using **`classification-insight-*`**), **toolbar** controls with **`classification-insight-meta`** labels, and **ECharts** in **`target-coverage-panel__plot-wrap`** (absolute and difference scatter plots). Section breaks use **`mgmt-detail-separator`** and **`target-coverage-panel__meta-label`** for **CNV events** / **Cytoband summary** headings.

**Plots (light + dark):** Theme is **not** inferred from Tailwind `dark:`. **`_is_dark_mode()`** reads **`app.storage.user['dark_mode']`** (Quasar **`body--dark`**). **`_apply_cnv_echart_chrome()`** sets **`backgroundColor`**, **`textStyle`**, **title** / **axis** / **tooltip** / **legend** colours from **`_cnv_echart_palette()`** (slate surfaces §6; insight-style tooltip background **§8**). Charts are re-chromed after every data refresh and **`_sync_cnv_echarts_theme_if_needed()`** runs on a short timer (and on **`0.05s` / `0.45s`** once) so toggling **Dark Mode** updates axes without waiting for a file refresh.

**Semantics:** Gains/losses and breakpoint overlays follow §5 where applicable; **breakpoint** **`markLine`** styling remains a high-visibility red on both themes.

**CSS:** **`#analysis-detail-cnv`** scoped **`.body--dark`** text fixes in **`styles.css`** (table/label greys).

**Reference implementation:** **`cnv.py`** — **`_apply_cnv_echart_chrome`**, **`_sync_cnv_echarts_theme_if_needed`**, state key **`cnv_plot_theme_dark`**.


9.7 Fusion detail section

The **expandable Fusion block** (`#analysis-detail-fusion` in **`fusion.py`**) uses the same **insight shell** as §9.6: **`classification-insight-shell`**, **`classification-insight-heading`** (**“Fusion analysis”** — **no emoji** in the section title), a top **`classification-insight-card`** with **`merge`** icon (**`classification-insight-icon`**), **Candidate summary** / target and genome-wide pair–group counts (**`classification-insight-meta`**), and a short **`classification-insight-foot`**. **Target panel** and **Genome-wide** subsections use **`target-coverage-panel__meta-label`** and **`mgmt-detail-separator`** between them. Tables and plots fill the shell below (matplotlib / DNA Features Viewer as implemented).

**CSS:** **`#analysis-detail-fusion`** scoped **`.body--dark`** text fixes in **`styles.css`** for legacy Tailwind greys on tables and status copy.

**Reference implementation:** **`fusion.py`** — **`add_fusion_section`**.


10. MNP-Flex Results (Epignostix)

The **MNP-Flex results** block on the sample detail page follows the same **insight strip** vocabulary as §8–§9: **`classification-insight-shell`**, section title via **`classification-insight-heading`**, **`classification-insight-card`** for sub-panels, **`classification-insight-icon`** / **`classification-insight-model`** / **`classification-insight-meta`** / **`classification-insight-foot`**, and **explicit `.body--dark` rules** in CSS—not Tailwind `dark:` alone.

10.1 Layout

- **Shell:** Full-width **`classification-insight-shell`** with **`id=mnpflex-results`** for deep-linking / scroll targets.
- **Toolbar:** **`classification-insight-meta`** for **Last updated**; **toolbar state** uses **`mnpflex-toolbar-badge--idle|running|busy`** (idle / running / preparing subset). Primary actions use Quasar **`color="primary"`** and **`outline`** for secondary.
- **Classifier:** One **insight card** with icon + label **Classifier**, then three **meta** lines (name, version, type).
- **Hierarchical summary:** Second card with **account_tree** icon; hierarchy table, **Top path** row with score badge, optional **unconfirmed classification** notice (**`mnpflex-notice`**, **amber** light / warm dark).
- **Hierarchy aggregates:** Four **insight cards** in **`classification-insight-grid`** (1 → 2 → 4 columns): **Subclass / Class / Family / Superfamily** with name + **score** badge. Duplicate set inside **Top 10 classifier scores** expansion when hierarchy is empty.
- **QC & MGMT:** **`classification-insight-grid--2`** (1 column mobile, **2 columns from ~1024px**): **Quality control** and **MGMT methylation** cards with **verified** / **science** icons, status badges, meta lines, **QC plots** / **MGMT plot** expansions.

10.2 Badges (no Tailwind-only colours)

- **Score badges** (hierarchy, aggregates, top path): **`mnpflex-score-badge`** tiers **`--neutral`** and **`--5`** … **`--1`** by numeric score bands (see `styles.css`).
- **QC / MGMT status:** **`mnpflex-status-badge--neutral|pass|warn|fail|meth|unmeth`**.
- **Errors:** **`mnpflex-error-text`** for inline error lines.

10.3 Reference implementation (ROBIN)

- **`mnpflex.py`:** **`add_mnpflex_section`** builds the UI; **`_score_badge_classes`** / **`_status_badge_classes`** return semantic **`mnpflex-*`** classes.
- **`styles.css`:** **`mnpflex-*`** and **`classification-insight-grid--2`** (two-column grid for QC/MGMT).


Version 1.16 | April 2026