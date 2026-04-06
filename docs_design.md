# R.O.B.I.N. Documentation Design System: The Digital Curator

This file is the **source of truth** for documentation styling. It is not published on the docs site; implementation lives in `docs/stylesheets/extra.css` and `mkdocs.yml`.

## 1. Overview and creative north star

Documentation in the R.O.B.I.N. ecosystem is more than just a list of features; it is an "Editorial Guide" for clinicians and researchers. It must feel clinical, authoritative, and high-fidelity, matching the premium aesthetic of the main analysis tools.

## 2. Information architecture: the modular approach

Documentation content should be broken down into high-impact cards rather than long, unbroken columns of text. This reduces cognitive load and allows for better scanning.

### 2.1 Hero and navigation

- **The breadcrumb:** Always provide clear context (for example, Docs > Genomic Pipelines > Sturgeon).
- **The table of contents:** On desktop, use a persistent sidebar for navigation. On mobile, use a clean "In this article" accordion.

### 2.2 Content blocks (card-based)

- **The feature card:** Use for high-level benefits (for example, "Rapid Results"). Include a modern line icon, a bold Manrope headline, and concise Inter body text.
- **The "case study" or "impact" card:** Use for real-world evidence. These can feature large background images with high-contrast text overlays to add visual variety.

## 3. Thematic depth and gradients

Documentation surfaces must follow the "Editorial" look defined in the core design system.

- **Subtle gradients:** Use light vertical gradients to prevent surfaces from feeling flat and digital.
- **Accent glows:** Use low-opacity radial gradients behind primary headers or key CTA buttons to drive focus.
- **Borders:** Keep borders minimal (1px) and semi-transparent. Avoid heavy shadows; prefer tonal shifts in background color to define depth.

## 4. Typography and color

### 4.1 Typeface pairings

- **Headlines:** Manrope (Bold/ExtraBold).
- **Body and labels:** Inter (Medium).
- **Technical snippets:** JetBrains Mono for API calls or code examples.

### 4.2 Semantic callouts (badges and alerts)

Use color to guide the user's attention without overwhelming them.

- **Information:** Slate 100 (light) / Slate 800 (dark).
- **Success/Action:** Emerald 600 (light) / Emerald 400 with glow (dark).
- **Warning/Critical:** Amber 500 (light) / Rose 400 (dark). Use these sparingly so they retain their alert value.

## 5. Theming implementation

### 5.1 Light mode (clinical clarity)

- **Background:** White (`#FFFFFF`) / Slate 50 (`#F8FAFC`).
- **Text:** Slate 900 for headlines, Slate 600 for body.
- **Cards:** Pure white with soft-shadowed edges and 1px Slate 200 borders.

### 5.2 Dark mode (midnight authority)

- **Background:** Zinc 950 (`#09090B`) / Slate 950 (`#020617`).
- **Text:** Slate 50 for headlines, Slate 300 for body.
- **Luminance:** Use high-saturation neon versions of brand colors (`#22c55e`) for interactive elements like links and buttons.

## 6. MkDocs implementation notes

Site-wide styling is applied through `docs/stylesheets/extra.css` (Material for MkDocs). Theme options are in `mkdocs.yml`.

- **Feature cards:** use a `div` with class `robin-feature-card` (for example inside a `.grid` block). Optional emoji or an element with class `robin-card-icon` can sit above the headline.
- **Impact / case-study cards:** use class `robin-impact-card`. For a photographic background, set a CSS variable on the element, for example `style="--robin-impact-bg: url('../images/example.jpg');"`.

---

Version 1.0 | October 2026
