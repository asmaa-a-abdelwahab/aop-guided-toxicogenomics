#!/usr/bin/env python3
"""Render an AOP-Wiki network from AOP fingerprint enrichment results."""

import argparse
import math
import os
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--enrichment",
    type=Path,
    default=PROJECT_ROOT / "outputs" / "AOP" / "AOP_fingerprint_Literature_enriched.csv",
    help="AOP fingerprint enrichment CSV",
)
parser.add_argument(
    "--ke-metadata",
    type=Path,
    default=PROJECT_ROOT / "data" / "aop_wiki" / "aop_ke_ec.tsv",
    help="AOP-Wiki KE metadata TSV",
)
parser.add_argument(
    "--relationships",
    type=Path,
    default=PROJECT_ROOT / "data" / "aop_wiki" / "aop_ke_ker.tsv",
    help="AOP-Wiki KE relationship TSV",
)
parser.add_argument(
    "--roles",
    type=Path,
    default=PROJECT_ROOT / "data" / "aop_wiki" / "aop_ke_mie_ao.tsv",
    help="AOP-Wiki KE role TSV",
)
parser.add_argument(
    "-o",
    "--outdir",
    type=Path,
    default=PROJECT_ROOT / "outputs" / "AOP_network",
    help="Output directory",
)
args = parser.parse_args()

try:
    import matplotlib.pyplot as plt
    import networkx as nx
    import numpy as np
    import pandas as pd
except ImportError as exc:
    parser.error(
        f"Missing Python dependency '{exc.name}'. "
        "Run: python -m pip install -r requirements.txt"
    )

try:
    from pyvis.network import Network
    HAS_PYVIS = True
except ImportError:
    HAS_PYVIS = False

in_enrich = args.enrichment
in_ke_ec = args.ke_metadata
in_ker = args.relationships
in_mieao = args.roles
out_dir = args.outdir

missing_inputs = [
    path for path in (in_enrich, in_ke_ec, in_ker, in_mieao) if not path.is_file()
]
if missing_inputs:
    formatted = "\n- ".join(str(path) for path in missing_inputs)
    parser.error(f"Required input files not found:\n- {formatted}")

out_dir.mkdir(parents=True, exist_ok=True)

SAFE_MODE = True
MAX_NODES = 500
MAX_EDGES = 2000

alpha = 0.05
min_neglog10 = -math.log10(alpha)

PNG_FILE = out_dir / "AOP_enrichment_network.png"
PNG_DPI = 300
FIGSIZE = (14, 10)

HTML_FILE = out_dir / "AOP_enrichment_network.html"
NODES_CSV = out_dir / "subgraph_nodes.csv"
EDGES_CSV = out_dir / "subgraph_edges.csv"

# ------------------------------- HELPERS --------------------------------------
def coalesce(df, candidates, new_name, as_numeric=False):
    for c in candidates:
        if c in df.columns:
            s = df[c]
            if as_numeric:
                s = pd.to_numeric(s, errors="coerce")
            df[new_name] = s
            return df
    df[new_name] = np.nan if as_numeric else pd.NA
    return df

def trimlower(s):
    return s.astype(str).str.strip().str.lower()

def html_escape(x):
    return (x.replace("&","&amp;")
             .replace("<","&lt;")
             .replace(">","&gt;")
             .replace('"',"&quot;")
             .replace("'","&#39;"))

def largest_weakly_connected_component(G):
    if G.number_of_nodes() == 0:
        return G
    comps = list(nx.weakly_connected_components(G))
    comps.sort(key=len, reverse=True)
    return G.subgraph(comps[0]).copy()

# ------------------------------ LOAD INPUTS -----------------------------------
enrich_raw = pd.read_csv(in_enrich)
ke_ec      = pd.read_csv(in_ke_ec, sep="\t")
ke_ker     = pd.read_csv(in_ker,   sep="\t")
mieao      = pd.read_csv(in_mieao, sep="\t")

# Normalize key columns
ke_ec = (coalesce(ke_ec, ["KE_id","KE.ID","KEID","KE_wiki_id","Event_ID"], "KE_id")
         .pipe(coalesce, ["KE_name","Event","Event_Name","KE Title","KE_title"], "KE_name"))
ke_ec["KE_id"] = ke_ec["KE_id"].astype(str)
ke_ec["KE_name"] = ke_ec["KE_name"].fillna(ke_ec["KE_id"])
ke_ec["KE_name_norm"] = trimlower(ke_ec["KE_name"])

ke_ker = (coalesce(ke_ker, ["AOP_id","AOP.ID","AOP","AOP_wiki_id"], "AOP_id")
          .pipe(coalesce, ["AOP_title","AOP_Title","Pathway","AOP_label"], "AOP_label")
          .pipe(coalesce, ["KE_up_id","KE_upstream_id","from","source","KE_source"], "from_id")
          .pipe(coalesce, ["KE_down_id","KE_downstream_id","to","target","KE_target"], "to_id")
          .pipe(coalesce, ["KER_id","KER.ID","Relationship_ID","Rel_ID"], "KER_id"))
for col in ["AOP_id","from_id","to_id"]:
    ke_ker[col] = ke_ker[col].astype(str)

mieao = (coalesce(mieao, ["AOP_id","AOP.ID","AOP_wiki_id"], "AOP_id")
         .pipe(coalesce, ["AOP_title","AOP_Title","AOP_label"], "AOP_label")
         .pipe(coalesce, ["KE_id","KE.ID","Event_ID"], "KE_id")
         .pipe(coalesce, ["KE_name","Event","Event_Name","KE Title","KE_title"], "KE_name")
         .pipe(coalesce, ["Role","KE_role","Node_Type"], "Role"))
mieao["AOP_id"] = mieao["AOP_id"].astype(str)
mieao["KE_id"]  = mieao["KE_id"].astype(str)
mieao["AOP_label"] = mieao["AOP_label"].fillna(mieao["AOP_id"])
mieao["AOP_label_norm"] = trimlower(mieao["AOP_label"])

# ---------------------- PROCESS ENRICHMENT TABLE ------------------------------
enrich = (enrich_raw
          .pipe(coalesce, ["KE","KeyEvent","Key_Event","KE_Name","KE_name","key_event"], "KE_name")
          .pipe(coalesce, ["KE_ID","KEID","KeyEventID","KE_Id","KE.id","KE_Wiki_ID","KE_AOPWiki_ID"], "KE_id")
          .pipe(coalesce, ["AOP","AOP_ID","AOP.id","AOPID","AOP_Title","Pathway","Pathway_Name"], "AOP_label")
          .pipe(coalesce, ["AOP_Wiki_ID","AOP_ID_num","AOP_Number","AOP_No","AOP.ID"], "AOP_id")
          .pipe(coalesce, ["padj","p.adjust","adj.p","FDR"], "padj", as_numeric=True)
          .pipe(coalesce, ["neglog10_padj","score","NES","-log10(padj)","minuslog10padj"], "neglog10_padj", as_numeric=True))

enrich["neglog10_padj"] = np.where(enrich["neglog10_padj"].isna() & ~enrich["padj"].isna(),
                                   -np.log10(enrich["padj"]), enrich["neglog10_padj"])
enrich["is_signif"] = np.where(~enrich["neglog10_padj"].isna(),
                               enrich["neglog10_padj"] >= min_neglog10,
                               np.where(~enrich["padj"].isna(), enrich["padj"] <= alpha, False))

enrich["KE_name_norm"]   = np.where(enrich["KE_name"].isna(), pd.NA, trimlower(enrich["KE_name"]))
enrich["AOP_label_norm"] = np.where(enrich["AOP_label"].isna(), pd.NA, trimlower(enrich["AOP_label"]))

# Map KE_name -> KE_id when missing
if enrich["KE_id"].isna().all() and enrich["KE_name_norm"].notna().any():
    tmp = ke_ec[["KE_id","KE_name_norm"]]
    enrich = enrich.merge(tmp, on="KE_name_norm", how="left", suffixes=("",".map"))
    if "KE_id.map" in enrich.columns:
        enrich["KE_id"] = enrich["KE_id"].fillna(enrich["KE_id.map"])
        enrich = enrich.drop(columns=["KE_id.map"])

# Map AOP_label -> AOP_id when missing
if enrich["AOP_id"].isna().all() and enrich["AOP_label_norm"].notna().any():
    tmp = mieao[["AOP_id","AOP_label_norm"]].drop_duplicates()
    enrich = enrich.merge(tmp, on="AOP_label_norm", how="left", suffixes=("",".map"))
    if "AOP_id.map" in enrich.columns:
        enrich["AOP_id"] = enrich["AOP_id"].fillna(enrich["AOP_id.map"])
        enrich = enrich.drop(columns=["AOP_id.map"])

# Fallback if no rows pass threshold
if not enrich["is_signif"].any():
    score = enrich["neglog10_padj"].copy()
    score = score.fillna(-np.log10(enrich["padj"].fillna(1.0)))
    enrich["fallback_rank"] = score.rank(ascending=False, method="dense")
    enrich["is_signif"] = enrich["fallback_rank"] <= 25

enriched_KEs = (enrich.loc[enrich["is_signif"], ["KE_id","KE_name","neglog10_padj","padj","AOP_id","AOP_label"]]
                .drop_duplicates())
enriched_AOPs = (enrich.loc[enrich["is_signif"], ["AOP_id","AOP_label"]]
                 .drop_duplicates())

print(">>> Enriched KEs n = {0}".format(len(enriched_KEs)))
print(">>> Enriched AOPs n = {0}".format(len(enriched_AOPs)))

# Second-pass name->ID if any enriched KE lacks ID
if enriched_KEs["KE_id"].isna().any() and enriched_KEs["KE_name"].notna().any():
    e2 = enriched_KEs.copy()
    e2["KE_name_norm"] = trimlower(e2["KE_name"])
    e2 = e2.merge(ke_ec[["KE_id","KE_name_norm"]], on="KE_name_norm", how="left", suffixes=("",".map"))
    if "KE_id.map" in e2.columns:
        e2["KE_id"] = e2["KE_id"].fillna(e2["KE_id.map"])
        e2 = e2.drop(columns=["KE_id.map"])
    e2 = e2.drop(columns=["KE_name_norm"])
    enriched_KEs = e2.drop_duplicates()

print(">>> Non-NA enriched KE_ids: {0}".format((~enriched_KEs["KE_id"].isna()).sum()))
print(">>> Non-NA enriched AOP_ids: {0}".format((~enriched_AOPs["AOP_id"].isna()).sum()))

# -------------------------- SUBGRAPH SELECTION --------------------------------
touched_aops = pd.concat([
    enriched_AOPs.loc[enriched_AOPs["AOP_id"].notna(), ["AOP_id"]],
    mieao.loc[mieao["KE_id"].isin(enriched_KEs["KE_id"].dropna()), ["AOP_id"]]
]).drop_duplicates()["AOP_id"].tolist()

if len(touched_aops) == 0:
    print("!!! No AOP IDs detected; will fall back to KERs that touch enriched KE IDs.")

if len(touched_aops) > 0:
    ker_sub = ke_ker[ke_ker["AOP_id"].isin(touched_aops)].copy()
else:
    valid_ke_ids = enriched_KEs["KE_id"].dropna().astype(str).unique().tolist()
    ker_sub = ke_ker[(ke_ker["from_id"].isin(valid_ke_ids)) | (ke_ker["to_id"].isin(valid_ke_ids))].copy()

if ker_sub.empty:
    print(">>> Still empty; attempting KE name-based matching.")
    ker_named = ke_ker.merge(ke_ec[["KE_id","KE_name"]], left_on="from_id", right_on="KE_id", how="left")
    ker_named = ker_named.rename(columns={"KE_name":"from_name"}).drop(columns=["KE_id"])
    ker_named = ker_named.merge(ke_ec[["KE_id","KE_name"]], left_on="to_id", right_on="KE_id", how="left")
    ker_named = ker_named.rename(columns={"KE_name":"to_name"}).drop(columns=["KE_id"])
    ker_sub = ker_named[ker_named["from_name"].isin(enriched_KEs["KE_name"]) |
                        ker_named["to_name"].isin(enriched_KEs["KE_name"])].copy()

ker_sub = ker_sub[ker_sub["from_id"].notna() & ker_sub["to_id"].notna()].copy()
print(">>> ker_sub rows: {0}".format(len(ker_sub)))

if ker_sub.empty:
    print("!!! No matching KERs found. Check that IDs/titles align between enrichment and AOP-Wiki extracts.")
    pd.DataFrame().to_csv(NODES_CSV, index=False)
    pd.DataFrame().to_csv(EDGES_CSV, index=False)
    raise SystemExit(0)

# ----------------------------- BUILD TABLES -----------------------------------
node_ids = pd.unique(pd.concat([ker_sub["from_id"], ker_sub["to_id"]], axis=0).astype(str))
nodes = pd.DataFrame({"KE_id": node_ids})
nodes = nodes.merge(ke_ec[["KE_id","KE_name"]], on="KE_id", how="left")

role_df = (mieao.groupby("KE_id", as_index=False)
           .agg(Role=("Role", lambda s: "; ".join(sorted(set([x for x in s.dropna().astype(str)]))) if s.notna().any() else pd.NA),
                AOP_ids=("AOP_id", lambda s: "; ".join(sorted(set([str(x) for x in s.dropna()])))),
                AOP_labels=("AOP_label", lambda s: "; ".join(sorted(set([str(x) for x in s.dropna()]))))
                ))
nodes = nodes.merge(role_df, on="KE_id", how="left")
nodes["Role"] = nodes["Role"].fillna("KE")

en_k = enriched_KEs[["KE_id","neglog10_padj","padj"]].drop_duplicates()
nodes = nodes.merge(en_k, on="KE_id", how="left")
nodes["is_enriched"]  = nodes["neglog10_padj"].fillna(-np.inf) >= min_neglog10
nodes["enrich_score"] = nodes["neglog10_padj"].fillna(0)
nodes["node_label"]   = nodes["KE_name"].fillna(nodes["KE_id"])

edges = ker_sub.copy()
if "KER_id" not in edges.columns:
    edges["KER_id"] = edges["from_id"].astype(str) + "->" + edges["to_id"].astype(str)
edges = edges.rename(columns={"from_id":"from","to_id":"to"})
edges["is_in_enriched_aop"] = edges["AOP_id"].astype(str).isin(enriched_AOPs["AOP_id"].astype(str).dropna())

# ------------------------------ BUILD GRAPH -----------------------------------
G = nx.DiGraph()
for _, r in nodes.iterrows():
    G.add_node(r["KE_id"],
               label=str(r["node_label"]),
               Role=str(r["Role"]) if pd.notna(r["Role"]) else "KE",
               is_enriched=bool(r["is_enriched"]),
               padj=None if pd.isna(r["padj"]) else float(r["padj"]),
               nlog10=None if pd.isna(r["neglog10_padj"]) else float(r["neglog10_padj"]),
               AOP_labels=None if pd.isna(r["AOP_labels"]) else str(r["AOP_labels"]))

for _, r in edges.iterrows():
    G.add_edge(str(r["from"]), str(r["to"]),
               AOP_id=None if pd.isna(r["AOP_id"]) else str(r["AOP_id"]),
               AOP_label=None if pd.isna(r["AOP_label"]) else str(r["AOP_label"]),
               is_in_enriched_aop=bool(r["is_in_enriched_aop"]))

G = largest_weakly_connected_component(G)

deg = dict(G.degree())
nx.set_node_attributes(G, deg, "deg")
btw = {}
if G.number_of_nodes() <= 5000:
    btw = nx.betweenness_centrality(G, normalized=True)
    nx.set_node_attributes(G, btw, "btw")

# ------------------------------- SAVE CSVs ------------------------------------
nodes_out = []
for n, d in G.nodes(data=True):
    nodes_out.append({
        "KE_id": n,
        "label": d.get("label"),
        "Role": d.get("Role"),
        "is_enriched": d.get("is_enriched"),
        "padj": d.get("padj"),
        "-log10(padj)": d.get("nlog10"),
        "AOP_labels": d.get("AOP_labels"),
        "deg": d.get("deg"),
        "btw": d.get("btw")
    })
pd.DataFrame(nodes_out).to_csv(NODES_CSV, index=False)

edges_out = []
for u, v, d in G.edges(data=True):
    edges_out.append({
        "from": u,
        "to": v,
        "AOP_id": d.get("AOP_id"),
        "AOP_label": d.get("AOP_label"),
        "is_in_enriched_aop": d.get("is_in_enriched_aop")
    })
pd.DataFrame(edges_out).to_csv(EDGES_CSV, index=False)

print(">>> Saved CSVs: subgraph_nodes.csv, subgraph_edges.csv")

# ------------------------------ STATIC PNG ------------------------------------
def node_color(d):
    return "#d62728" if d.get("is_enriched") else "#9ecae1"

def node_shape(d):
    role = (d.get("Role") or "").upper()
    if "AO" in role:
        return "^"
    if "MIE" in role:
        return "s"
    return "o"

# Try Graphviz layout first; fallback to spring_layout
try:
    from networkx.drawing.nx_agraph import graphviz_layout
    pos = graphviz_layout(G, prog="dot")
except Exception:
    pos = nx.spring_layout(G, seed=42, k=1/max(1, math.sqrt(max(1, G.number_of_nodes()))))

plt.figure(figsize=FIGSIZE)
edge_colors = ["#9E9E9E" if d.get("is_in_enriched_aop") else "#D0D0D0" for _,_,d in G.edges(data=True)]
nx.draw_networkx_edges(G, pos, arrows=True, arrowstyle="-|>", arrowsize=8, edge_color=edge_colors, width=1.2, alpha=0.6)

for shp in ["o","s","^"]:
    nlist = [n for n,d in G.nodes(data=True) if node_shape(d)==shp]
    if not nlist:
        continue
    colors = [node_color(G.nodes[n]) for n in nlist]
    sizes  = [80 + 8*math.sqrt(max(0, G.nodes[n].get("deg",0))) for n in nlist]
    nx.draw_networkx_nodes(G, pos, nodelist=nlist, node_color=colors, node_size=sizes,
                           node_shape=shp, linewidths=0.4, edgecolors="#333333", alpha=0.95)

# labels: enriched + top-degree 20
deg_sorted = sorted(G.nodes(), key=lambda x: G.nodes[x].get("deg",0), reverse=True)[:20]
labels_nodes = set([n for n,d in G.nodes(data=True) if d.get("is_enriched")]) | set(deg_sorted)
labels = dict((n, G.nodes[n].get("label","")) for n in labels_nodes)
nx.draw_networkx_labels(G, pos, labels=labels, font_size=8, font_weight="bold")

plt.title("AOP-Wiki Network with Enriched KEs / Pathways\nRed=enriched | Square=MIE | Triangle=AO | Circle=KE", fontsize=11)
plt.axis("off")
plt.tight_layout()
plt.savefig(PNG_FILE, dpi=PNG_DPI)
plt.close()
print(">>> Saved PNG:", os.path.basename(PNG_FILE))

# --------------------------- INTERACTIVE HTML ---------------------------------
if HAS_PYVIS:
    nt = Network(height="800px", width="100%", directed=True, notebook=False, bgcolor="#ffffff", font_color="#222")
    nt.barnes_hut()

    nodes_for_interactive = list(G.nodes())
    edges_for_interactive = list(G.edges())

    if SAFE_MODE:
        enriched_nodes = [n for n,d in G.nodes(data=True) if d.get("is_enriched")]
        top_deg_nodes  = [n for n,_ in sorted(G.degree(), key=lambda kv: kv[1], reverse=True)[:MAX_NODES]]
        keep = set(enriched_nodes + top_deg_nodes)
        if len(keep) < MAX_NODES:
            for n in nodes_for_interactive:
                if len(keep) >= MAX_NODES:
                    break
                keep.add(n)
        nodes_for_interactive = [n for n in nodes_for_interactive if n in keep]
        edges_for_interactive = [(u,v) for u,v in edges_for_interactive if u in keep and v in keep]
        if len(edges_for_interactive) > MAX_EDGES:
            edges_for_interactive = edges_for_interactive[:MAX_EDGES]

    for n in nodes_for_interactive:
        d = G.nodes[n]
        title = "<b>{0}</b>".format(html_escape(str(d.get("label",""))))
        if d.get("Role"):
            title += "<br/>Role: {0}".format(html_escape(str(d["Role"])))
        if d.get("AOP_labels"):
            title += "<br/>AOPs: {0}".format(html_escape(str(d["AOP_labels"])))
        if d.get("padj") is not None:
            title += "<br/>padj: {0}".format("{0:.3g}".format(d["padj"]))
        if d.get("nlog10") is not None:
            title += "<br/>-log10(padj): {0}".format("{0:.3f}".format(d["nlog10"]))
        color = "#d62728" if d.get("is_enriched") else "#6aaed6"
        shape = "dot"
        role  = (d.get("Role") or "").upper()
        if "AO" in role:
            shape = "triangle"
        elif "MIE" in role:
            shape = "square"
        size = max(8, min(30, 8*math.sqrt(max(1, d.get("deg",1)))))
        nt.add_node(n, label=str(d.get("label","")), title=title, color=color, shape=shape, value=size)

    for (u,v) in edges_for_interactive:
        d = G.edges[u, v]
        title = "AOP: {0}".format(d.get("AOP_label") or d.get("AOP_id") or "")
        color = "#9E9E9E" if d.get("is_in_enriched_aop") else "#D0D0D0"
        nt.add_edge(u, v, title=title, arrows="to", color=color)

    nt.set_options("""
    var options = {
      physics: { stabilization: true, solver: 'forceAtlas2Based' },
      interaction: { hover: true, tooltipDelay: 120 },
      nodes: { borderWidth: 1 }
    }
    """)
    nt.save_graph(str(HTML_FILE))
    print(">>> Saved HTML:", os.path.basename(HTML_FILE))
else:
    print(">>> PyVis not installed; skipping interactive HTML. Install with: pip install pyvis")

print(">>> Done.")
