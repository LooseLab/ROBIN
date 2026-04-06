# Troubleshooting

!!! abstract "What this page covers"
    Common **browser** issues: can’t load the page, sign-in, **Unknown sample**, watched folders errors, dark mode, **Workflow** menu, failed reports. For server install issues, see the [README — Common issues](https://github.com/LooseLab/ROBIN/blob/main/README.md#common-issues) and [CLI reference](../cli/index.md).

---

## I can’t open the page {#i-cant-open-the-page}

| Check | Action |
|-------|--------|
| **URL** | Copy exactly from the terminal or your admin (`http://…`, correct **port**, often **8081**). |
| **Network** | Same LAN/VPN if ROBIN runs on another machine. |
| **Firewall** | IT may need to allow the port. |
| **Process running** | The web UI is only up while **`robin workflow`** is running. Restart or ask your admin. |
| **Work directory** | Some setups need **`--work-dir`** for the main dashboard — [What happens at startup](../getting-started/startup.md). |

---

## I can’t sign in {#i-cant-sign-in}

- Use the **same password** as for the web UI (often set at first **GUI password** prompt in the terminal).  
- **Caps Lock** off; retry after a short wait.  
- If no password was ever set, run ROBIN **from a terminal** once — [At startup](../getting-started/startup.md).  
- Ask an admin to reset: **`robin password set`**.  

---

## It says “Unknown sample” {#unknown-sample}

ROBIN only knows samples it has **seen** (or that are already tracked). If the ID is wrong or BAMs have not arrived:

- Wait for **BAMs** in the watched folder.  
- Open **All samples** and **View** the correct **row**.  
- Confirm the **library ID** matches MinKNOW.  

---

## “Watched folders” mentions Ray or the workflow {#watched-folders-ray}

That screen **adds or removes input folders**. It only works in configurations that support it. Most clinical users do **not** change this mid-run — ask bioinformatics if you see an error.

---

## Dark mode looks wrong {#dark-mode}

Toggle **Dark Mode** in the menu off and on. If a **chart** is wrong, **refresh** or leave the sample and reopen — plots sometimes update after a moment.

---

## The “Workflow” menu item doesn’t work {#workflow-menu}

Use **Activity Monitor** — that is the **workflow monitor** in the standard setup.

---

## Report or download failed {#report-failed}

- Wait for spinners / notifications to finish.  
- Retry; large PDFs can take time.  
- Ask your admin to check **disk space** and **write permissions** on the output folder.  

---

## Still stuck {#still-stuck}

- [README — Common issues](https://github.com/LooseLab/ROBIN/blob/main/README.md#common-issues)  
- [Command-line reference](../cli/index.md) (staff who run ROBIN)  
