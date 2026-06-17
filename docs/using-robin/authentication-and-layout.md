# First steps and navigation

!!! abstract "What this page covers"
    Open the **URL**, **sign in** if asked, and use the **menu** (☰), **Dark mode**, **Links**, and **Log out**. Explains the **top bar** (CPU/RAM) and **Quit** vs closing the tab.

---

## Opening ROBIN

Paste the address into your browser. Common forms:

- `http://127.0.0.1:8081` or `http://localhost:8081`  
- On another machine: use the **host** and **port** your administrator gave you  

If nothing loads: [Troubleshooting — I can’t open the page](troubleshooting.md#i-cant-open-the-page).

---

## Signing in

If the GUI is enabled, you’ll see **Sign in** with **username** and **password** fields.

- Enter the credentials your administrator provided (or the `admin` account created at setup).  
- If you are signing in with a **temporary** password, you must **choose a new password** before continuing.  
- On first login, accept the **research-use agreement** — acceptance is recorded per user in the audit database.  
- Press **Log in** or **Enter**.  

After sign-in, ROBIN opens the page you wanted or **Welcome**.

For account setup and auditing, see [Users, consent, and auditing](audit-and-users.md).

---

## Top bar (every screen)

| Element | Purpose |
|---------|---------|
| **Menu (☰)** | Opens main navigation (see below). |
| **Title** | Full ROBIN name on large screens; **R.O.B.I.N** on small screens. |
| **Viewing:** … | Which machine the browser is talking to. |
| **CPU** / **RAM** | Quick load indicators — **not** a full system monitor. |
| **Logo** | ROBIN branding. |

---

## Main menu (☰)

| Goal | Menu item |
|------|-----------|
| Home | **Home** |
| Table of all samples | **View Samples** |
| Create a library ID from test ID + optional fields | **Generate Sample ID** |
| Add/remove BAM watch folders | **Watched Folders** *(may need a specific pipeline mode)* |
| Pipeline status | **Activity Monitor** |
| Published docs (new tab) | **Documentation** |
| Dark theme | **Dark Mode** |
| Remote access (advanced) | **Allow Remote Access** |
| Update your password | **Change password** |
| User accounts and audit log *(admin only)* | **Administration** |
| End web session | **LOG OUT** |
| Stop ROBIN on this machine *(some setups)* | **Quit** — read the warning; closing the **tab** often leaves analysis running. |

!!! tip "“Workflow” in the menu"
    Some builds list **Workflow** separately. In the usual **`robin workflow`** setup, use **Activity Monitor** for pipeline status. If **Workflow** errors or does nothing, prefer **Activity Monitor**.

**Close** dismisses the menu.

---

## Footer

**Links** — shortcuts to GitHub, papers, lab protocol, Oxford Nanopore, and related sites. **Info** may show extra help on some builds.

---

## Quit vs closing the browser

If you see **Quit R.O.B.I.N?**:

- **Cancel** — dialog closes; analysis may continue.  
- **Really quit** — can stop analysis on that machine.  
- Often you can **close the tab** and leave ROBIN running — ask your operator if unsure.  

---

## Next

[Tour of the screens](pages-and-routes.md) — each page and what to click.
