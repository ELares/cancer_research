# Security Policy

## Supported Versions

| Component | Version | Supported |
|-----------|---------|-----------|
| Python pipeline (scripts/) | latest main | ✅ |
| ferroptosis-core (Rust) | 0.2.0 | ✅ |
| ferroptosis-python bindings | 0.1.0 | ✅ |
| Manuscript & analysis | latest main | ✅ |

## Reporting a Vulnerability

If you discover a security vulnerability in this project, please report it responsibly:

1. **Do NOT open a public GitHub issue.** Security vulnerabilities should be reported privately.
2. **Email:** Send details to the repository owner via [GitHub private vulnerability reporting](https://github.com/ELares/cancer_research/security/advisories/new).
3. **Include:** A description of the vulnerability, steps to reproduce, and the potential impact.

## What to Expect

- **Acknowledgment** within 72 hours of your report.
- **Assessment** of severity and impact within 1 week.
- **Fix or mitigation** as soon as practical, prioritized by severity.
- **Credit** in the fix commit and release notes (unless you prefer anonymity).

## Scope

This project contains:

- **Python scripts** that fetch data from external APIs (PubMed, OpenAlex, Semantic Scholar) using API keys stored in `.env` (gitignored). Vulnerabilities in API key handling, data parsing, or injection are in scope.
- **Rust simulation binaries** that process numerical data. These do not handle untrusted user input in production but memory safety issues are still relevant.
- **GitHub Actions workflows** that run CI/CD. Workflow permission escalation or secret exposure is in scope.
- **A corpus of research articles** stored as markdown files. These are not executable but could contain injection vectors if processed unsafely.

## Out of Scope

- Vulnerabilities in third-party dependencies (report these to the upstream project, though we appreciate a heads-up so we can update).
- The scientific correctness of simulation results or manuscript claims (use GitHub issues for these).
- The bundled OpenStax textbooks in `books/` (these are CC BY 4.0 reference materials, not project code).

## Security Practices

- API keys are stored in `.env` (gitignored) and never logged or committed.
- GitHub Actions workflows use explicit `permissions` blocks with least privilege.
- Dependabot monitors Python (pip), Rust (Cargo), and GitHub Actions dependencies weekly.
- The Rust simulation engine uses safe Rust with no `unsafe` blocks in ferroptosis-core.
