"""Custom OpenAPI documentation theme, metadata, and styling for VNtyper Online API.

Provides structured, accessible descriptions, tags metadata, and custom
styling for Swagger UI adhering to WCAG AA contrast standards.
"""

from typing import Any

API_SUMMARY = "Genotyping and cohort analysis for MUC1 VNTR in ADTKD-MUC1"

API_DESCRIPTION = """Genotyping pipeline for *MUC1* Variable Number Tandem Repeats (VNTR) in Autosomal Dominant Tubulointerstitial Kidney Disease (**ADTKD-MUC1**) from short-read sequencing data (BAM/CRAM).

## Overview

| Feature | Details | Action / Reference |
| :--- | :--- | :--- |
| **Pipeline** | K-mer haplotype assembly via Kestrel and optional profile-HMM repeat counting | [VNtyper 2 Repository](https://github.com/hassansaei/vntyper) |
| **Workflow** | Submit alignment files, monitor task execution, and retrieve results | `POST /run-job/` → `GET /download/` |
| **Cohorts** | Passphrase-protected sample grouping for joint multi-sample genotyping | `POST /create-cohort/` |
| **Data Retention** | Automatic purging of uploaded files and result archives after 3 days | 3-day retention policy |
| **Validation** | Genotyping methodology and benchmarking in ADTKD-MUC1 cohorts | [Popp & Saei et al., medRxiv (2026)](https://doi.org/10.64898/2026.05.27.26352937) |
"""

API_TAGS_METADATA: list[dict[str, Any]] = [
    {
        "name": "General",
        "description": "Health checks and version endpoints.",
    },
    {
        "name": "Job Management",
        "description": "Submit alignments (BAM/CRAM), check job status, and download results.",
    },
    {
        "name": "Cohort Management",
        "description": "Passphrase-protected sample cohorts and joint analysis.",
    },
    {
        "name": "Usage Statistics",
        "description": "Platform throughput and run statistics.",
    },
    {
        "name": "Research Data Donations",
        "description": "Optional anonymous variant data contributions for benchmark validation.",
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
  width: 100% !important;
  max-width: 100% !important;
  padding: 0 32px 48px !important;
  box-sizing: border-box !important;
}

/* Info Header Card */
.swagger-ui .info {
  margin: 24px 0 20px !important;
  background: #ffffff;
  border: 1px solid #e2e8f0;
  border-radius: 12px;
  padding: 28px 32px !important;
  box-shadow: 0 1px 3px rgba(0, 0, 0, 0.05);
  box-sizing: border-box !important;
  width: 100% !important;
}

.swagger-ui .info .title {
  font-size: 28px;
  font-weight: 800;
  color: #0f172a;
  letter-spacing: -0.025em;
  display: flex;
  align-items: center;
  flex-wrap: wrap;
  gap: 10px;
  margin-bottom: 8px;
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
  margin-top: 4px;
}

.swagger-ui .info .description {
  font-size: 15px;
  line-height: 1.6;
  color: #334155;
  margin-top: 16px;
  width: 100%;
}

.swagger-ui .info .description p {
  margin: 10px 0;
  max-width: 75ch;
  line-height: 1.6;
}

.swagger-ui .info .description h2 {
  font-size: 18px;
  font-weight: 700;
  color: #0f172a;
  margin: 24px 0 12px;
  padding-bottom: 6px;
  border-bottom: 1px solid #e2e8f0;
}

.swagger-ui .info .description table {
  width: 100% !important;
  border-collapse: collapse !important;
  margin: 16px 0 !important;
  border-radius: 8px !important;
  overflow: hidden !important;
  border: 1px solid #e2e8f0 !important;
  box-sizing: border-box !important;
}

.swagger-ui .info .description th {
  background: #f1f5f9 !important;
  color: #0f172a !important;
  font-weight: 700 !important;
  font-size: 13px !important;
  text-transform: uppercase !important;
  letter-spacing: 0.05em !important;
  padding: 10px 16px !important;
  border: 1px solid #e2e8f0 !important;
  text-align: left !important;
}

.swagger-ui .info .description td {
  padding: 12px 16px !important;
  border: 1px solid #e2e8f0 !important;
  color: #334155 !important;
  font-size: 14px !important;
  line-height: 1.5 !important;
}

.swagger-ui .info .description tr:nth-child(even) td {
  background: #f8fafc !important;
}

.swagger-ui .info .description th:first-child,
.swagger-ui .info .description td:first-child {
  white-space: nowrap !important;
  width: 130px !important;
  font-weight: 600 !important;
  color: #0f172a !important;
}

.swagger-ui .info .description th:nth-child(2),
.swagger-ui .info .description td:nth-child(2) {
  width: 55% !important;
}

.swagger-ui .info .description th:nth-child(3),
.swagger-ui .info .description td:nth-child(3) {
  width: 35% !important;
}

.swagger-ui .info .description code {
  font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
  font-size: 12.5px;
  background: #e2e8f0;
  color: #0f172a;
  padding: 2px 6px;
  border-radius: 4px;
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

/* Scheme / Server Selector Card */
.swagger-ui .scheme-container {
  background: #ffffff !important;
  border: 1px solid #e2e8f0 !important;
  border-radius: 12px !important;
  padding: 16px 28px !important;
  box-shadow: 0 1px 3px rgba(0, 0, 0, 0.05) !important;
  margin: 0 0 20px !important;
  box-sizing: border-box !important;
  width: 100% !important;
}

.swagger-ui .scheme-container .schemes {
  margin: 0 !important;
  padding: 0 !important;
  display: flex !important;
  align-items: center !important;
  gap: 12px !important;
}

.swagger-ui .scheme-container .schemes > label {
  font-weight: 600 !important;
  color: #0f172a !important;
}

/* Tag Sections & Filter */
.swagger-ui .filter {
  margin: 0 0 20px !important;
  padding: 0 !important;
  width: 100% !important;
}

.swagger-ui .filter .operation-filter-input {
  width: 100% !important;
  max-width: 100% !important;
  border: 1px solid #cbd5e1 !important;
  border-radius: 8px !important;
  padding: 10px 16px !important;
  font-size: 14px !important;
  box-sizing: border-box !important;
  box-shadow: 0 1px 2px rgba(0, 0, 0, 0.05) !important;
}

.swagger-ui .filter .operation-filter-input:focus {
  border-color: #0369a1 !important;
  outline: 2px solid rgba(3, 105, 161, 0.25) !important;
}

.swagger-ui .opblock-tag-section {
  margin-top: 28px !important;
  margin-bottom: 20px !important;
  width: 100% !important;
}

.swagger-ui .opblock-tag {
  font-size: 18px;
  font-weight: 700;
  color: #0f172a;
  border-bottom: 2px solid #cbd5e1;
  padding: 12px 0 8px;
  margin-top: 0 !important;
  margin-bottom: 10px !important;
}

.swagger-ui .opblock-tag small,
.swagger-ui .opblock-tag p,
.swagger-ui .opblock-tag-section .renderedMarkdown p {
  font-size: 14px;
  font-weight: 400;
  color: #475569;
  padding: 6px 12px !important;
  display: inline-block !important;
  max-width: 75ch !important;
  margin: 0 !important;
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
    favicon_url: str = "/api/favicon.svg",
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
  <link rel="icon" type="image/svg+xml" href="{favicon_url}">
  <link rel="shortcut icon" href="/api/favicon.ico">
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
