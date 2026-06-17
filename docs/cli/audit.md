# `robin audit`

Query and export GUI audit events (login, sample views, report exports, admin actions).

Audit events are append-only records in `~/.config/robin/security.db`.

## `robin audit list`

```bash
robin audit list
robin audit list --user admin --event auth.login.success --limit 20
robin audit list --from-ts 2026-06-01T00:00:00+00:00 --to-ts 2026-06-30T23:59:59+00:00
```

## `robin audit export`

```bash
robin audit export --out audit_june.csv
robin audit export --user alice --event report.exported --out alice_exports.csv
```

Exported CSV columns include timestamp, user, event type, target, result, IP, session/request IDs, and a JSON details field.

## Common event types

| Event | When recorded |
|-------|----------------|
| `auth.login.success` / `auth.login.failure` | GUI sign-in |
| `auth.logout` | User logs out |
| `auth.password.changed` | User sets a new password (forced first-login or voluntary) |
| `consent.accepted` | User accepts research-use agreement |
| `sample.viewed` / `sample.list.viewed` | Sample pages opened |
| `sample.audit.viewed` / `sample.audit.exported` | Per-sample audit history opened or exported |
| `report.generated` | PDF/CSV/TSV report built from the GUI |
| `report.exported` | Existing report file downloaded via the API |
| `run.started` | Manual run/job submission from GUI (SNP, MNP-Flex, finalize, etc.) |
| `admin.user.*` | User management CLI actions |

## Related

- [`robin users`](users.md)  
- [Users, consent, and auditing](../using-robin/audit-and-users.md)
