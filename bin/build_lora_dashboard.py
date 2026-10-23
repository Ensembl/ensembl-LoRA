#!/usr/bin/env python
import argparse
import html
import json
import os


def main():
    args = get_args()
    os.makedirs(args.output_dir, exist_ok=True)
    with open(args.dashboard_data) as handle:
        data = json.load(handle)

    with open(os.path.join(args.output_dir, "dashboard_data.json"), "w") as handle:
        json.dump(data, handle, indent=2)

    with open(os.path.join(args.output_dir, "index.html"), "w") as handle:
        handle.write(render_html(data))


def get_args():
    parser = argparse.ArgumentParser(description="Render the static LoRA dashboard")
    parser.add_argument("--dashboard_data", required=True)
    parser.add_argument("--output_dir", required=True)
    return parser.parse_args()


def render_html(data):
    run = data["run_summary"]
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>LoRA Dashboard: {esc(run['run_id'])} {esc(run['genome'])}</title>
<style>
:root {{ --ink:#172033; --muted:#607085; --line:#d8e1ea; --soft:#f6f8fb; --blue:#2563eb; --green:#16865a; --amber:#c76a21; --purple:#6d5bb8; }}
body {{ margin:0; font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif; color:var(--ink); background:var(--soft); }}
header {{ background:var(--ink); color:white; padding:22px 30px; }}
h1 {{ margin:0 0 5px; font-size:28px; letter-spacing:0; }}
main {{ max-width:1500px; margin:0 auto; padding:22px 30px 44px; }}
h2 {{ margin:0 0 10px; font-size:20px; }}
h3 {{ margin:15px 0 8px; font-size:15px; }}
.muted {{ color:var(--muted); font-size:13px; line-height:1.4; }}
header .muted {{ color:#cbd5e1; }}
.tabs {{ display:flex; gap:8px; flex-wrap:wrap; margin:0 0 18px; }}
button {{ border:1px solid #cbd5df; background:white; color:var(--ink); border-radius:7px; padding:8px 11px; cursor:pointer; font-weight:650; }}
button.active {{ background:var(--ink); color:white; border-color:var(--ink); }}
.view {{ display:none; }}
.view.active {{ display:block; }}
.grid {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(210px,1fr)); gap:12px; }}
.two {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(390px,1fr)); gap:14px; }}
.card {{ background:white; border:1px solid var(--line); border-radius:8px; padding:15px; }}
.metric {{ border-left:4px solid var(--blue); min-height:104px; }}
.metric.good {{ border-left-color:var(--green); }}
.metric.warn {{ border-left-color:var(--amber); }}
.metric.purple {{ border-left-color:var(--purple); }}
.value {{ font-size:30px; line-height:1; font-weight:820; font-variant-numeric:tabular-nums; margin:7px 0; }}
.barrow {{ display:grid; grid-template-columns:minmax(150px,250px) minmax(0,1fr) 82px; gap:10px; align-items:center; margin:8px 0; font-size:13px; }}
.track {{ height:12px; background:#e6edf4; border-radius:999px; overflow:hidden; }}
.fill {{ height:100%; background:var(--blue); }}
.barval {{ text-align:right; font-variant-numeric:tabular-nums; }}
table {{ width:100%; border-collapse:collapse; background:white; border:1px solid var(--line); font-size:12.5px; }}
th,td {{ padding:8px 9px; border-bottom:1px solid #e7edf3; text-align:left; vertical-align:top; }}
th {{ background:#eef3f8; font-weight:760; white-space:nowrap; }}
.callout {{ border-left:4px solid var(--amber); background:#fff8ef; border-radius:6px; padding:11px 12px; font-size:13px; line-height:1.42; }}
.mono {{ font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace; font-size:12px; color:#2d3d52; }}
@media (max-width:620px) {{ main {{ padding:16px; }} .barrow {{ grid-template-columns:1fr; gap:4px; }} .barval {{ text-align:left; }} }}
</style>
</head>
<body>
<header>
<h1>LoRA Dashboard</h1>
<div class="muted">{esc(run['run_id'])} · {esc(run['genome'])} · generated {esc(data['generated_at'])}</div>
</header>
<main>
<div class="tabs">
<button class="active" data-view="overview">Overview</button>
<button data-view="checks">Data Checks</button>
<button data-view="models">Models</button>
<button data-view="sqanti">SQANTI3</button>
<button data-view="files">Files</button>
<button data-view="help">Help</button>
</div>
{render_overview(data)}
{render_data_checks(data)}
{render_models(data)}
{render_sqanti(data)}
{render_files(data)}
{render_help()}
<script>
document.querySelectorAll('button[data-view]').forEach(button => {{
  button.addEventListener('click', () => {{
    document.querySelectorAll('button[data-view]').forEach(item => item.classList.toggle('active', item === button));
    document.querySelectorAll('.view').forEach(panel => panel.classList.toggle('active', panel.dataset.panel === button.dataset.view));
  }});
}});
</script>
</main>
</body>
</html>
"""


def render_overview(data):
    run = data["run_summary"]
    pre = data["models"]["prefiltered"]
    post = data["models"]["postfiltered"]
    sqanti_available = data["sqanti3"].get("postfiltered", {}).get("available", False)
    return f"""
<section class="view active" data-panel="overview">
<div class="grid">
{metric("Prefiltered models", pre["transcript_count"], "post-tmerge non-redundant set")}
{metric("Filtered models", post["transcript_count"], f"full-length threshold {esc(run['full_length_threshold'])}", "good")}
{metric("Input reads", run["input_reads"], esc(run.get("input_reads_source", "not available")), "purple")}
{metric("Canonical reads", run["canonical_reads"], "unique reads represented in LoRA models", "warn")}
</div>
<div class="two" style="margin-top:14px">
<div class="card"><h2>Run Summary</h2>{key_table([
("Genome", run["genome"]),
("Run ID", run["run_id"]),
("Full-length threshold", run["full_length_threshold"]),
("GenomeSize parameter", run["genomeSize"]),
("Input files", run.get("input_file_count", "not_available")),
("Input samples", run["input_sample_count"]),
("Input tissues", run.get("input_tissue_count", "not_available")),
("Input cell lines", run.get("input_cell_line_count", "not_available")),
("Samples missing tissue/cell line", run.get("input_samples_missing_tissue_cell_line", "not_available")),
("Grouped input files", run.get("grouped_input_files", "not_available")),
("Collapsed sample groups", run.get("collapsed_sample_groups", "not_available")),
("SQANTI3", "available" if sqanti_available else "not available"),
])}</div>
<div class="card"><h2>Input Reads By Method</h2>{bars(run.get("read_counts_by_method", {}))}
<h3>Reads By Tissue</h3>{bars(run.get("read_counts_by_tissue", {}))}
<h3>Reads By Cell Line</h3>{bars(run.get("read_counts_by_cell_line", {}))}</div>
</div>
<div class="two" style="margin-top:14px">
<div class="card"><h2>Samples Per Method</h2>{bars(run.get("samples_per_method", {}))}
<h3>Tissues</h3>{bars(run.get("tissue_counts", {}))}
<h3>Cell Lines</h3>{bars(run.get("cell_line_counts", {}))}</div>
<div class="card"><h2>Read Metadata Completeness</h2>{bars(run.get("read_counts_missing_metadata", {}))}
<h3>Tissue / Cell Line Read Counts</h3>{bars(run.get("read_counts_by_tissue_cell_line", {}))}</div>
</div>
</section>
"""


def render_models(data):
    pre = data["models"]["prefiltered"]
    post = data["models"]["postfiltered"]
    return f"""
<section class="view" data-panel="models">
<div class="two">
<div class="card"><h2>Prefiltered Model Set</h2>{model_table(pre)}</div>
<div class="card"><h2>Filtered Model Set</h2>{model_table(post)}</div>
</div>
<div class="two" style="margin-top:14px">
<div class="card"><h2>Prefiltered Support</h2>{bars(pre.get("support_bins", {}))}</div>
<div class="card"><h2>Filtered Support</h2>{bars(post.get("support_bins", {}))}</div>
</div>
</section>
"""


def render_data_checks(data):
    checks = data.get("data_checks", {})
    status_counts = checks.get("status_counts", {})
    stages = checks.get("stages", {})
    cards = "".join(
        metric(label.title(), status_counts.get(label, 0), "data checks", cls)
        for label, cls in [
            ("passed", "good"),
            ("warning", "warn"),
            ("failed", "warn"),
            ("skipped", "purple"),
        ]
    )
    stage_sections = []
    for stage, rows in stages.items():
        stage_sections.append(f'<div class="card"><h2>{esc(stage)}</h2>{table(rows, ["check_name", "status", "observed", "expected", "message", "source_file"])}</div>')
    if not stage_sections:
        stage_sections.append('<div class="card"><h2>Data Checks</h2><div class="muted">No data-check rows were available.</div></div>')
    return f"""
<section class="view" data-panel="checks">
<div class="grid">{cards}</div>
<div style="margin-top:14px; display:grid; gap:14px">
{''.join(stage_sections)}
</div>
</section>
"""


def render_sqanti(data):
    sqanti = data["sqanti3"]
    if not sqanti.get("postfiltered", {}).get("available", False):
        return """
<section class="view" data-panel="sqanti">
<div class="card"><h2>SQANTI3 Not Available</h2><div class="callout">No post-filtered SQANTI3 classification file was found. This usually means no reference GTF was supplied for this genome or SQANTI3 was disabled.</div></div>
</section>
"""
    return f"""
<section class="view" data-panel="sqanti">
<div class="two">
{sqanti_stage_card("Prefiltered", sqanti.get("prefiltered", {}))}
{sqanti_stage_card("Post Filtered", sqanti.get("postfiltered", {}))}
</div>
<div class="two" style="margin-top:14px">
<div class="card"><h2>Structural Category Change</h2>{table(sqanti.get("delta_structural_category", []), ["label", "prefiltered", "postfiltered"])}</div>
<div class="card"><h2>Coding Status Change</h2>{table(sqanti.get("delta_coding", []), ["label", "prefiltered", "postfiltered"])}</div>
</div>
</section>
"""


def sqanti_stage_card(title, stage):
    if not stage.get("available"):
        return f'<div class="card"><h2>{esc(title)}</h2><div class="muted">Not available</div></div>'
    return f"""
<div class="card"><h2>{esc(title)}</h2>
<div class="grid">
{metric("Classified", stage.get("classification_rows", 0), "classification rows")}
{metric("Non-canonical", stage.get("non_canonical", 0), "from all_canonical", "warn")}
</div>
<h3>Structural Category</h3>{bars(list_to_dict(stage.get("structural_category", [])))}
<h3>Subcategory Grouped By Structural Category</h3>{table(stage.get("subcategory_by_structural_category", [])[:16], ["structural_category", "subcategory", "count"])}
<h3>Canonical, Coding, and CDS Type</h3>
{small_summary_tables(stage)}
<h3>Transcript-End Evidence</h3>
{transcript_end_evidence_tables(stage)}
</div>
"""


def small_summary_tables(stage):
    rows = []
    for label, key in [("all_canonical", "all_canonical"), ("coding", "coding"), ("CDS_type", "CDS_type")]:
        for item in stage.get(key, []):
            rows.append({"field": label, "label": item["label"], "count": item["count"]})
    return table(rows, ["field", "label", "count"])


def transcript_end_evidence_tables(stage):
    rows = []
    for label, key in [("within_CAGE_peak", "within_CAGE_peak"), ("within_polyA_site", "within_polyA_site")]:
        for item in stage.get(key, []):
            rows.append({"field": label, "label": item["label"], "count": item["count"]})
    if not rows:
        return '<div class="muted">Not available</div>'
    return table(rows, ["field", "label", "count"])


def render_files(data):
    rows = [{"file": key, "path": value} for key, value in data.get("files", {}).items()]
    rows.extend({"file": key, "path": value} for key, value in data.get("summary_files", {}).items())
    return f'<section class="view" data-panel="files"><div class="card"><h2>Output File Map</h2>{table(rows, ["file", "path"])}</div></section>'


def render_help():
    return """
<section class="view" data-panel="help">
<div class="two">
<div class="card"><h2>LoRA Terms</h2>
<table><tbody>
<tr><th>Prefiltered models</th><td>Non-redundant transcript models after tmerge and before full-length filtering.</td></tr>
<tr><th>Filtered models</th><td>Transcript models retained after LoRA full-length support filtering.</td></tr>
<tr><th>Contains file</th><td>LoRA table linking transcript models to supporting reads.</td></tr>
<tr><th>Canonical reads</th><td>Reads retained after LoRA splice-site orientation/canonical filtering and represented in model support.</td></tr>
</tbody></table></div>
<div class="card"><h2>SQANTI3 Terms</h2>
<table><tbody>
<tr><th>Full-Splice Match (FSM)</th><td>All splice junctions match a known reference transcript.</td></tr>
<tr><th>Incomplete-Splice Match (ISM)</th><td>A transcript matches consecutive junctions from a reference but is incomplete.</td></tr>
<tr><th>Novel In Catalog (NIC)</th><td>Uses known splice donors/acceptors in a new combination.</td></tr>
<tr><th>Novel Not in Catalog (NNC)</th><td>Contains at least one novel splice donor or acceptor.</td></tr>
<tr><th>Open Reading Frame (ORF), Coding Sequence (CDS), Nonsense-Mediated Decay (NMD)</th><td>Functional triage fields from SQANTI3 when ORF prediction is enabled.</td></tr>
</tbody></table></div>
</div>
</section>
"""


def model_table(summary):
    return key_table([
        ("Transcript count", summary["transcript_count"]),
        ("Single-exon", summary["single_exon"]),
        ("Multi-exon", summary["multi_exon"]),
        ("Median length", summary["median_length"]),
        ("Min length", summary["min_length"]),
        ("Max length", summary["max_length"]),
    ])


def metric(label, value, hint, cls=""):
    return f'<div class="card metric {cls}"><div class="muted">{esc(label)}</div><div class="value">{esc(value)}</div><div class="muted">{esc(hint)}</div></div>'


def key_table(rows):
    return "<table><tbody>" + "".join(f"<tr><th>{esc(k)}</th><td>{esc(v)}</td></tr>" for k, v in rows) + "</tbody></table>"


def bars(values):
    items = sorted(values.items(), key=lambda item: item[1], reverse=True) if isinstance(values, dict) else values
    if not items:
        return '<div class="muted">Not available</div>'
    max_value = max(int(value) for _label, value in items) or 1
    output = []
    for label, value in items:
        width = max(2, int(value) / max_value * 100)
        output.append(f'<div class="barrow"><div>{esc(label)}</div><div class="track"><div class="fill" style="width:{width:.1f}%"></div></div><div class="barval">{esc(value)}</div></div>')
    return "".join(output)


def table(rows, cols):
    if not rows:
        return '<div class="muted">Not available</div>'
    return "<table><thead><tr>" + "".join(f"<th>{esc(col)}</th>" for col in cols) + "</tr></thead><tbody>" + "".join(
        "<tr>" + "".join(f"<td>{esc(row.get(col, ''))}</td>" for col in cols) + "</tr>" for row in rows
    ) + "</tbody></table>"


def list_to_dict(rows):
    return {row["label"]: row["count"] for row in rows}


def esc(value):
    return html.escape(str(value))


if __name__ == "__main__":
    main()
