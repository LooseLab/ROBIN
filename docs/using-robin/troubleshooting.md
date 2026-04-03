# Troubleshooting

## I can’t open the page

- **Check the address** — copy it exactly from the terminal or the note your administrator gave you (`http://…` and the correct **port**, often `8081`).  
- **Same network** — if ROBIN runs on another computer, you may need to be on the hospital network or VPN.  
- **Firewall** — IT may need to allow the port through.  
- **Workflow still running** — ROBIN only serves the web UI while the `robin workflow` process is active. If someone stopped it, restart the workflow or ask your admin.  
- **Work directory** — in some setups the main browser dashboard only appears if the workflow was started with an output **work directory** set; see [What happens at startup](../getting-started/startup.md).

---

## I can’t sign in

- Use the **same password** your team agreed for the web browser (often set the first time ROBIN asked for a **GUI password** in the terminal).  
- **Caps Lock** off; try again after a short delay.  
- If you never set a password, run ROBIN **from a terminal** once so it can prompt you—see [At startup](../getting-started/startup.md).  
- Ask your administrator to **reset** or **re-set** the password if needed (`robin password set`).

---

## It says “Unknown sample”

ROBIN only knows about samples it has **seen in this session** (or that are already in its tracking list). If you opened a link with the wrong ID, or data hasn’t arrived yet:

- Wait for **BAM files** to appear in the watched folder.  
- Go back to **All samples** and click **View** on the correct row in the **table**.  
- Confirm the **library ID** matches MinKNOW.

---

## “Watched folders” says something about Ray or the workflow

That page is for **adding or removing input folders**. It only works when ROBIN was started in a configuration that supports it. **Clinical users** rarely need this during a run—ask your bioinformatics team if you see an error here.

---

## Dark mode looks wrong

Toggle **Dark Mode** in the menu off and on. If a **chart** still looks wrong, **refresh the page** or leave the sample and open it again—plots sometimes catch up after a moment.

---

## The “Workflow” menu item doesn’t work

Use **Activity Monitor** instead. That’s the **workflow monitor** screen in the standard setup.

---

## Report or download failed

- Wait for any **spinner** or **notification** to finish.  
- Try again; large PDFs can take a minute.  
- Ask your admin to check that the **output folder** has space and that you have **permission** to write there.

---

## Still stuck

- [README — Common issues](https://github.com/LooseLab/ROBIN/blob/main/README.md#common-issues)  
- [Command-line reference](../cli/index.md) (for staff who run ROBIN)  
