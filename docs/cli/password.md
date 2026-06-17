# `robin password`

Configure the legacy **GUI password hash file** used when bootstrapping the first admin user.

For multi-user account management, see **[`robin users`](users.md)**.

For the full startup sequence, see **[What happens at startup](../getting-started/startup.md)**.

## `robin password set`

```bash
robin password set
```

- Prompts twice for the new password (no echo).
- Replaces any previously stored password.

If the GUI password module fails to import (minimal install / broken env), the command prints an error and exits non-zero.

## When it applies

- **`robin password set`** writes `gui_password_hash` under your config directory.
- On first GUI launch with no users, ROBIN can **bootstrap** user `admin` from this hash.
- Prefer **`robin users bootstrap-admin`** for explicit first-time setup.

The GUI login screen uses **username + password** accounts stored in `security.db`.

## Related

- [`robin users`](users.md)  
- [`robin audit`](audit.md)  
- [`robin workflow` GUI options](workflow.md#gui-nicegui)  
- [CLI overview](index.md)  
