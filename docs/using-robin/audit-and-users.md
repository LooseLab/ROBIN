# Users, consent, and auditing

ROBIN supports **named user accounts** with role-based access and a full **audit trail** of GUI activity.

All users can view all samples in the monitored directory; access control is focused on **accountability** (who did what, when, from which IP) rather than dataset isolation.

---

## Accounts and roles

| Role | Typical use |
|------|-------------|
| `admin` | Initial bootstrap account; can manage users via CLI |
| `user` | Standard operator account |

Roles are assigned when creating accounts (`robin users create --role user`) or adjusted with `grant-role` / `revoke-role`.

---

## Signing in

1. Open the ROBIN GUI URL (e.g. `http://localhost:8081`).
2. Enter **username** and **password**.
3. On first login (or after a consent version change), accept the **research-use agreement**.
4. Use the app normally; actions are logged to the audit database.

Log out via **LOG OUT** in the menu. This ends your session and records `auth.logout`.

---

## Research-use consent

The startup “I agree” dialog has been replaced by **per-user consent at login**:

- Each acceptance is stored with user ID, consent version, timestamp, IP, and session metadata.
- Default consent version is `v1`.
- Set `ROBIN_CONSENT_VERSION` (e.g. `v2`) to require re-acceptance after policy updates.

Check consent status:

```bash
robin users consent-status
```

---

## What is audited

Examples of recorded events:

- Login success and failure
- Logout
- Consent acceptance
- Viewing sample list or individual sample pages
- Generating reports (PDF, CSV ZIP, sample-tracking TSV) from the GUI
- Downloading existing report files via the API
- Starting runs/jobs from the GUI (SNP, MNP-Flex, finalize, etc.)
- Admin CLI changes (user create, password reset, role changes)

Query recent events:

```bash
robin audit list --limit 50
```

Export for compliance review:

```bash
robin audit export --out robin_audit.csv
```

---

## Administration (CLI)

See the CLI references:

- [`robin users`](../cli/users.md) — account management  
- [`robin audit`](../cli/audit.md) — audit query/export  

Admins can also use the GUI **Administration** page (`/admin`, menu item visible when signed in as `admin`):

- View users, roles, consent status, and last login
- Create users, reset passwords, activate/deactivate accounts, grant/revoke admin
- Filter and export audit events to CSV

---

## Data location

SQLite database: `~/.config/robin/security.db`

Back up this file to preserve users, consent records, and audit history.
