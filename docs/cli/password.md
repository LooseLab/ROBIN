# `robin password`

Configure the **default admin account** used for GUI sign-in.

For multi-user account management, see **[`robin users`](users.md)**.

For the full startup sequence, see **[What happens at startup](../getting-started/startup.md)**.

## `robin password set`

```bash
robin password set
```

- Prompts twice for the new password (no echo).
- If no GUI users exist yet, creates the default **`admin`** account.
- If **`admin`** already exists, asks for confirmation before replacing that password.

Optional:

```bash
robin password set --username admin
```

If the security module fails to import (minimal install / broken env), the command prints an error and exits non-zero.

## When it applies

- **`robin password set`** updates the user database (`security.db`), not a standalone password file.
- On first GUI launch with no users, ROBIN can also prompt in the terminal for the same default admin password.
- Prefer **`robin users bootstrap-admin`** when you want an explicit bootstrap command with audit logging only from that path.

The GUI login screen uses **username + password** accounts stored in `security.db`.

## Related

- [`robin users`](users.md)  
- [`robin audit`](audit.md)  
- [`robin workflow` GUI options](workflow.md#gui-nicegui)  
- [CLI overview](index.md)  
