"""Custom OpenAPI documentation theme, metadata, and styling for VNtyper Online API.

Provides structured, accessible descriptions, tags metadata, and custom
styling for Swagger UI adhering to WCAG AA contrast standards.
"""

from typing import Any

API_SUMMARY = "High-speed MUC1 VNTR genotyping and cohort analysis pipeline for ADTKD-MUC1"

API_DESCRIPTION = """The **VNtyper Online API** provides programmatic access to **VNtyper 2**, an advanced pipeline engineered for high-speed genotyping of *MUC1* Variable Number Tandem Repeats (VNTR) in Autosomal Dominant Tubulointerstitial Kidney Disease (**ADTKD-MUC1**) using Short-Read Sequencing (SRS) alignments.

---

## API Overview

The API supports automated end-to-end processing of genomic alignment files (BAM/BAI and CRAM), rapid variant calling, asynchronous background execution, and passphrase-protected multi-sample cohort management.

---

## Key Capabilities

- **High-Speed Variant Calling**: Genotyping via k-mer analysis using Kestrel for rapid local haplotype assembly and frameshift detection.
- **Repeat-Unit Validation**: Optional code-adVNTR profile-HMM repeat-unit validation for in-depth motif architecture analysis.
- **Asynchronous Processing**: Scalable Celery task queue architecture with real-time job status polling and automatic artifact retention.
- **Passphrase-Protected Cohorts**: Group multiple related samples into secure, collaborative cohorts with credential protection, joint multi-sample genotyping, and aggregated call tables.
- **Runtime Configuration Discovery**: Query server-side defaults and administrator policy enforcement flags via `GET /options-config/`.
- **Anonymized Research Donations**: Opt-in framework for contributing anonymous run summaries and variant calls to support algorithm benchmarking.
- **Platform Usage Metrics**: Transparent system throughput tracking including cumulative jobs processed, unique users, and recent 24-hour activity.

---

## Analysis Modes

| Mode | Engine | Best For | Typical Processing Time |
| :--- | :--- | :--- | :--- |
| **Normal Mode** | Kestrel local haplotype assembly | Rapid screening & routine ADTKD-MUC1 diagnostics | ~1–3 minutes |
| **adVNTR Mode** | Profile-HMM repeat-unit counting | Deep motif validation & repeat-unit expansion analysis | ~15–45 minutes |

*Tip: Query `GET /options-config/` to inspect whether your server enforces or defaults adVNTR or normal analysis modes.*

---

## Workflow Guide

### 1. Single-Sample Genotyping

1. **Submit Alignment**: Upload a BAM file with BAI index to `POST /run-job/`. Specify optional notification email, cohort alias, or analysis mode.
2. **Poll Status**: Periodically query `GET /job-status/{job_id}/` to monitor job progress and queue position.
3. **Download Results**: Once the status is `completed`, retrieve the complete results bundle from `GET /download/{job_id}/`.

### 2. Passphrase-Protected Cohort Analysis

1. **Create Cohort**: Register a new cohort via `POST /create-cohort/` with an alias and passphrase.
2. **Join Samples**: Submit multiple samples referencing the cohort alias and passphrase.
3. **Trigger Analysis**: Request collective analysis with `POST /cohort-analysis/` passing the cohort authorization header.
4. **Download Summary**: Retrieve aggregated cohort call summaries via `GET /cohort-download/`.

---

## Security & Privacy

- **Cohort Authentication**: All operations modifying or querying cohorts require the cohort passphrase. We recommend using the **`X-Cohort-Passphrase`** HTTP header to protect credentials from appearing in URL logs.
- **Rate Limiting**: Tiered token-bucket rate limiting applies to all public endpoints. If exceeded, the API returns `429 Too Many Requests`.
- **Data Retention**: Uploaded alignments and result archives are retained for **3 days** (`MAX_RESULT_AGE_DAYS`) before automatic purging.

---

## Citations & Research

If you use VNtyper Online in your research, please cite:

> **Popp B, Saei H, et al.** *VNtyper 2 enables open-access short-read genotyping of MUC1 VNTR variants in ADTKD at high-speed.* **medRxiv** (2026). [doi:10.64898/2026.05.27.26352937](https://doi.org/10.64898/2026.05.27.26352937)
"""

API_TAGS_METADATA: list[dict[str, Any]] = [
    {
        "name": "General",
        "description": "System health checks, service discovery, and API & bioinformatics tool version discovery.",
    },
    {
        "name": "Job Management",
        "description": (
            "Submit genomic alignment files (BAM), configure analysis options, query asynchronous "
            "job execution status, and download final result packages."
        ),
    },
    {
        "name": "Cohort Management",
        "description": (
            "Group multiple related samples into secure, passphrase-protected cohorts for collective genotyping, "
            "cross-sample allele frequency metrics, and aggregated analysis."
        ),
    },
    {
        "name": "Usage Statistics",
        "description": "Access anonymized platform usage metrics, including cumulative jobs, unique users, and recent 24-hour activity.",
    },
    {
        "name": "Research Data Donations",
        "description": "Opt-in anonymous research donation endpoints to help benchmark and improve VNTR genotyping accuracy.",
    },
]

SWAGGER_CUSTOM_CSS = """
/* Custom VNtyper Swagger UI Theme - WCAG AA Compliant & Clean Aesthetic */
body {
  margin: 0;
  padding: 0;
  background-color: #f8fafc;
  color: #0f172a;
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
  -webkit-font-smoothing: antialiased;
}

.swagger-ui {
  color: #0f172a;
  font-family: inherit;
}

.swagger-ui .topbar {
  display: none !important;
}

.swagger-ui .wrapper {
  max-width: 1360px;
  padding: 0 24px 48px;
}

/* Info Header Card */
.swagger-ui .info {
  margin: 32px 0 24px;
  background: #ffffff;
  border: 1px solid #e2e8f0;
  border-radius: 12px;
  padding: 32px;
  box-shadow: 0 1px 3px rgba(0, 0, 0, 0.05);
}

.swagger-ui .info .title {
  font-size: 30px;
  font-weight: 800;
  color: #0f172a;
  letter-spacing: -0.025em;
  display: flex;
  align-items: center;
  flex-wrap: wrap;
  gap: 10px;
}

.swagger-ui .info .title small {
  background: #e2e8f0;
  color: #1e293b;
  font-size: 13px;
  font-weight: 700;
  padding: 6px 12px !important;
  border-radius: 6px !important;
  display: inline-block !important;
  line-height: 1.2 !important;
  top: auto;
}

.swagger-ui .info .title small.version-stamp {
  background: #0369a1;
  color: #ffffff;
  border: none;
  font-size: 12px;
  font-weight: 700;
  padding: 6px 12px !important;
  border-radius: 6px !important;
}

.swagger-ui .info .base-url {
  font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
  font-size: 13px;
  color: #475569;
  margin-top: 8px;
}

.swagger-ui .info .description {
  font-size: 15px;
  line-height: 1.7;
  color: #334155;
  margin-top: 20px;
}

.swagger-ui .info .description p,
.swagger-ui .info .description li,
.swagger-ui .info .description blockquote,
.swagger-ui .opblock-tag small,
.swagger-ui .opblock-tag p,
.swagger-ui .opblock-tag-section p {
  max-width: 72ch !important;
  margin-left: 0 !important;
  margin-right: auto !important;
}

.swagger-ui .info .description h2 {
  font-size: 20px;
  font-weight: 700;
  color: #0f172a;
  margin: 28px 0 12px;
  padding-bottom: 8px;
  border-bottom: 2px solid #e2e8f0;
}

.swagger-ui .info .description h3 {
  font-size: 16px;
  font-weight: 600;
  color: #1e293b;
  margin: 20px 0 8px;
}

.swagger-ui .info .description hr {
  border: 0;
  height: 1px;
  background: #e2e8f0;
  margin: 24px 0;
}

.swagger-ui .info .description table {
  width: 100%;
  max-width: 900px;
  border-collapse: collapse;
  margin: 16px 0;
  border-radius: 8px;
  overflow: hidden;
  border: 1px solid #cbd5e1;
}

.swagger-ui .info .description th {
  background: #f1f5f9;
  color: #0f172a;
  font-weight: 700;
  padding: 10px 14px;
  border: 1px solid #cbd5e1;
}

.swagger-ui .info .description td {
  padding: 10px 14px;
  border: 1px solid #cbd5e1;
  color: #334155;
}

.swagger-ui .info .description tr:nth-child(even) {
  background: #f8fafc;
}

.swagger-ui .info .description pre {
  background: #0f172a;
  color: #f8fafc;
  padding: 14px 18px;
  border-radius: 8px;
  overflow-x: auto;
  font-size: 13px;
  line-height: 1.5;
  border: 1px solid #334155;
  max-width: 900px;
}

.swagger-ui .info .description code {
  font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
  font-size: 13px;
  background: #f1f5f9;
  color: #0f172a;
  padding: 2px 6px;
  border-radius: 4px;
  border: 1px solid #e2e8f0;
}

.swagger-ui .info .description pre code {
  background: transparent;
  border: none;
  padding: 0;
  color: inherit;
}

.swagger-ui .info .description a,
.swagger-ui .info a {
  color: #0369a1;
  text-decoration: underline;
  text-underline-offset: 2px;
  font-weight: 600;
}

.swagger-ui .info .description a:hover,
.swagger-ui .info a:hover {
  color: #075985;
}

/* Tag Sections & Filter */
.swagger-ui .opblock-tag-section {
  margin-top: 36px !important;
  margin-bottom: 20px !important;
}

.swagger-ui .opblock-tag {
  font-size: 19px;
  font-weight: 700;
  color: #0f172a;
  border-bottom: 2px solid #cbd5e1;
  padding: 12px 0 8px;
  margin-top: 0 !important;
  margin-bottom: 10px !important;
}

.swagger-ui .opblock-tag small {
  font-size: 14px;
  font-weight: 400;
  color: #475569;
  padding: 6px 12px !important;
  display: inline-block !important;
}

.swagger-ui .filter {
  margin-bottom: 24px;
}

.swagger-ui .filter .operation-filter-input {
  border: 1px solid #cbd5e1;
  border-radius: 8px;
  padding: 10px 16px;
  font-size: 14px;
  width: 100%;
  max-width: 440px;
  box-shadow: 0 1px 2px rgba(0, 0, 0, 0.05);
}

.swagger-ui .filter .operation-filter-input:focus {
  border-color: #0369a1;
  outline: 2px solid rgba(3, 105, 161, 0.25);
}

/* Method Badges - WCAG AA Contrast Compliant */
.swagger-ui .opblock.opblock-get .opblock-summary-method {
  background: #0369a1 !important; /* Sky 700: 5.9:1 contrast with white */
  color: #ffffff !important;
  font-weight: 700;
  border-radius: 6px;
  min-width: 80px;
}

.swagger-ui .opblock.opblock-post .opblock-summary-method {
  background: #15803d !important; /* Green 700: 5.0:1 contrast with white */
  color: #ffffff !important;
  font-weight: 700;
  border-radius: 6px;
  min-width: 80px;
}

.swagger-ui .opblock.opblock-delete .opblock-summary-method {
  background: #b91c1c !important; /* Red 700: 5.6:1 contrast with white */
  color: #ffffff !important;
  font-weight: 700;
  border-radius: 6px;
  min-width: 80px;
}

.swagger-ui .opblock.opblock-put .opblock-summary-method {
  background: #b45309 !important; /* Amber 700: 4.6:1 contrast with white */
  color: #ffffff !important;
  font-weight: 700;
  border-radius: 6px;
  min-width: 80px;
}

.swagger-ui .opblock.opblock-patch .opblock-summary-method {
  background: #0f766e !important; /* Teal 700: 5.1:1 contrast with white */
  color: #ffffff !important;
  font-weight: 700;
  border-radius: 6px;
  min-width: 80px;
}

/* Operation Cards */
.swagger-ui .opblock {
  border-radius: 8px;
  margin: 0 0 12px;
  box-shadow: 0 1px 2px rgba(0, 0, 0, 0.04);
  background: #ffffff;
  transition: all 0.15s ease;
}

.swagger-ui .opblock:hover {
  box-shadow: 0 4px 12px rgba(0, 0, 0, 0.08);
}

.swagger-ui .opblock.opblock-get {
  border: 1px solid #bae6fd;
  background: rgba(3, 105, 161, 0.03);
}

.swagger-ui .opblock.opblock-post {
  border: 1px solid #bbf7d0;
  background: rgba(21, 128, 61, 0.03);
}

.swagger-ui .opblock .opblock-summary {
  padding: 10px 16px;
}

.swagger-ui .opblock .opblock-summary-path {
  font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
  font-size: 14px;
  font-weight: 600;
  color: #0f172a;
}

.swagger-ui .opblock .opblock-summary-description {
  font-size: 13px;
  color: #334155;
  font-weight: 400;
}

/* Models & JSON Schema Contrast Fixes */
.swagger-ui .json-schema-2020-12-expand-deep-button {
  color: #0f172a !important;
  background-color: #cbd5e1 !important;
  font-weight: 600 !important;
  border-radius: 4px;
  padding: 4px 10px !important;
  border: 1px solid #94a3b8 !important;
}

.swagger-ui .json-schema-2020-12-expand-deep-button:hover {
  background-color: #94a3b8 !important;
}

.swagger-ui .json-schema-2020-12-accordion {
  background-color: #f1f5f9 !important;
  color: #0f172a !important;
  font-weight: 600 !important;
}

.swagger-ui .parameter__name {
  color: #0f172a;
  font-weight: 600;
}

.swagger-ui .parameter__type,
.swagger-ui .prop-format,
.swagger-ui .parameter__in {
  color: #334155 !important;
}

.swagger-ui .tabli .tablinks {
  color: #334155 !important;
  font-weight: 600;
}

.swagger-ui .tabli.active .tablinks {
  color: #0369a1 !important;
  border-bottom-color: #0369a1 !important;
}

.swagger-ui .btn.execute {
  background-color: #0369a1 !important;
  border-color: #0369a1 !important;
  color: #ffffff !important;
  font-weight: 700;
  border-radius: 6px;
  padding: 8px 24px;
}

.swagger-ui .btn.execute:hover {
  background-color: #075985 !important;
  border-color: #075985 !important;
}

.swagger-ui .btn.try-out__btn {
  border-radius: 6px;
  font-weight: 600;
  color: #1e293b;
  border-color: #cbd5e1;
}

.swagger-ui select {
  border-radius: 6px;
  border: 1px solid #cbd5e1;
  padding: 6px 12px;
  color: #0f172a;
}
"""


def render_custom_swagger_ui_html(
    openapi_url: str,
    title: str,
    oauth2_redirect_url: str | None = None,
) -> str:
    """Generate custom Swagger UI HTML embedding accessible CSS and configuration."""
    redirect_script = ""
    if oauth2_redirect_url:
        redirect_script = f"oauth2RedirectUrl: window.location.origin + '{oauth2_redirect_url}',"

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{title}</title>
  <link type="text/css" rel="stylesheet" href="https://cdn.jsdelivr.net/npm/swagger-ui-dist@5/swagger-ui.css">
  <link rel="icon" type="image/svg+xml" href="data:image/svg+xml,<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 100 100'><text y='.9em' font-size='90'>🧬</text></svg>">
  <style>{SWAGGER_CUSTOM_CSS}</style>
</head>
<body>
  <div id="swagger-ui"></div>
  <script src="https://cdn.jsdelivr.net/npm/swagger-ui-dist@5/swagger-ui-bundle.js"></script>
  <script>
    window.addEventListener('DOMContentLoaded', () => {{
      window.ui = SwaggerUIBundle({{
        url: '{openapi_url}',
        dom_id: '#swagger-ui',
        docExpansion: 'list',
        defaultModelsExpandDepth: -1,
        filter: true,
        displayRequestDuration: true,
        tryItOutEnabled: true,
        syntaxHighlight: {{ theme: 'monokai' }},
        {redirect_script}
        presets: [
          SwaggerUIBundle.presets.apis,
          SwaggerUIBundle.SwaggerUIStandalonePreset
        ]
      }});
    }});
  </script>
</body>
</html>
"""
