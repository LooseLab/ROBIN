# `robin users`

Manage GUI user accounts, roles, and research-use consent status.

User and audit data are stored in SQLite at `~/.config/robin/security.db` (or `%APPDATA%\robin\security.db` on Windows).

## First-time setup

Create the first admin account explicitly:

```bash
robin users bootstrap-admin
```

Or import the existing single-password hash from a prior ROBIN install:

```bash
robin password set          # if not already set
robin users bootstrap-admin --from-legacy-hash
```

This creates user `admin` with the same password as before.

Then sign in to the GUI with **username + password** (not password alone).

## User commands

| Command | Purpose |
|---------|---------|
| `robin users bootstrap-admin` | Create the first admin account |
| `robin users create <username> [--role admin\|user]` | Add a user |
| `robin users set-password <username>` | Reset a user's password |
| `robin users list` | List accounts |
| `robin users activate <username>` | Re-enable a deactivated account |
| `robin users deactivate <username>` | Disable login (cannot remove last active admin) |
| `robin users grant-role <username> <admin\|user>` | Assign a role |
| `robin users revoke-role <username> <admin\|user>` | Remove a role |
| `robin users consent-status [--version v1]` | Show who has accepted the active consent version |

## Consent versioning

Research-use agreement acceptance is recorded **per user** at login. The active version defaults to `v1` and can be overridden:

```bash
export ROBIN_CONSENT_VERSION=v2
```

When the version changes, each user must accept again on next login. Acceptance is stored in the `consents` table and logged as `consent.accepted` in the audit trail.

## Related

- [`robin audit`](audit.md) — query and export audit events  
- [`robin password set`](password.md) — legacy GUI password file (used for bootstrap import)  
- [First steps and navigation](../using-robin/authentication-and-layout.md)
