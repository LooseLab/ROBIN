# Using ROBIN

Welcome. This guide is for **anyone using the ROBIN web app in a browser**—whether you are new to the tool or only need a reminder of where things live.

ROBIN’s interface shows how your sequencing run is progressing and surfaces **classification**, **coverage**, **CNV**, **MGMT**, **fusion**, and other results as they become available. You do **not** need to read code or command-line manuals to follow this guide; we link to technical docs only where it helps.

## Before you start

1. **ROBIN must already be running** on your computer or server (someone usually starts it with `robin workflow` from a terminal). When it starts, it prints a **web address**—for example `http://127.0.0.1:8081`. Open that address in **Chrome, Edge, Firefox, or Safari**.
2. You may be asked for a **password** in the browser or in the terminal first—see [First steps and navigation](authentication-and-layout.md).
3. The window title will show **ROBIN** so you can find the tab easily.

## What you’ll learn here

| Guide | You’ll learn how to… |
|-------|----------------------|
| [First steps and navigation](authentication-and-layout.md) | Sign in, use the menu, turn on **Dark mode**, find **Links** and **Log out**, and understand the top bar (CPU/RAM). |
| [Tour of the screens](pages-and-routes.md) | Move from the **Welcome** page to **samples**, the **workflow monitor**, **watched folders**, and the **sample ID** helper. |
| [Reading your results](sample-results.md) | Read the **Run summary**, **Classification**, **Analysis**, and **reports** on a sample page. |
| [MNP-Flex setup](../getting-started/mnpflex.md) | Configure Epignostix credentials, submit samples, and troubleshoot the optional MNP-Flex integration. |
| [Troubleshooting](troubleshooting.md) | Fix common problems (can’t connect, wrong sample, dark mode, reports). |

## A simple workflow

1. Open the link you were given → land on **Welcome** (or sign in first).  
2. Open **View Samples** (or **View All Samples**) → click **View** on the row for your **library / sample** to open it.  
3. Scroll through **Run summary**, then **Classification details**, then **Analysis details** as numbers appear.  
4. Use **Generate report** when you need a **PDF** (and optional data export if your team enabled it).  
5. To see overall pipeline health, open **Activity Monitor** from the menu.

When you’re ready for installation, startup prompts, or command-line options, see [Quickstart](../getting-started/quickstart.md) and [What happens at startup](../getting-started/startup.md).
