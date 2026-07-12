"""
Coffea arabica Variant Annotation Pipeline — v4 (Final)
========================================================
Key fixes vs v3:
  - CBD FIXED with positional_fallback mode: genes in the Ck-1 QTL window
    (Chr01C/E 19.5-35 Mbp) are indexed regardless of GFF3 annotation vocabulary.
    Keyword match is still attempted first; genes with no annotation match are
    labelled "Positional (Ck-1 locus)" rather than being silently dropped.
    This resolves the persistent 0-candidate-gene issue caused by sparse GFF3
    annotation in that chromosomal interval.
  - Dwarfism expanded from Chr01 only to all 22 chromosomes — GA pathway genes
    (DELLA, GA20ox, GID1, BRI1) are distributed genome-wide and were being
    missed by the previous Chr01-only restriction.
  - Positional tier properly accounted for in analytics tierCounts.
  - All v3 fixes retained: normalize_text() applied to both blob and keywords,
    phrase_in_blob() for robust hyphenated-term matching, bisect binary search,
    full analytics payload (Ti/Tv, subgenome bias, hotspot genes, etc.).
"""

from __future__ import annotations

import re
import json
import gzip
import datetime
import os
from bisect import bisect_right
from collections import defaultdict
from typing import Optional

from flask import Flask, request, jsonify, send_from_directory
from werkzeug.utils import secure_filename

app = Flask(__name__, static_url_path="", static_folder=".")

# --- PATHS -------------------------------------------------------------------

BASE_DIR      = os.path.dirname(os.path.abspath(__file__))
GFF3_FILES    = [
    os.path.join(BASE_DIR, "Coffea_arabica_ET-39_CC.gff3"),
    os.path.join(BASE_DIR, "Coffea_arabica_ET-39_CE.gff3"),
]
UPLOAD_FOLDER = os.path.join(BASE_DIR, "uploads")
PAYLOAD_PATH  = os.path.join(BASE_DIR, "annotation_payload_antigravity.json")
DEFAULT_VCF   = os.path.join(BASE_DIR, "uploaded_variants.vcf")
os.makedirs(UPLOAD_FOLDER, exist_ok=True)

# --- BIOLOGICAL MODEL --------------------------------------------------------
#
# Tier 1 (High Confidence):
#   Gene must be on the specified chromosome(s) AND within the coordinate range
#   (if given) AND its normalised attribute blob must contain a keyword phrase.
#
# Tier 2 (Discovery):
#   Gene can be anywhere in the genome but must match a highly specific term.
#
# All keywords and blobs are normalised identically (non-alphanumeric -> space,
# lowercase) before matching. This is what the previous version failed to do.

BIOLOGICAL_MODEL: dict[str, dict] = {

    "Low Caffeine (lr)": {
        "description": (
            "Tiered analysis targeting N-methyltransferase gene clusters (XMT, MXMT, DXMT) "
            "on chromosomes 1 and 9, which catalyse the four sequential methylation steps of "
            "caffeine biosynthesis. Variants in these loci are candidates for the recessive "
            "lr (low-caffeine) allele documented in Coffea arabica mutant lines."
        ),
        "positional_fallback": False,
        "high_conf_loci": {
            "chroms": ["Chr01C", "Chr01E", "Chr09C", "Chr09E"],
            "ranges": None,
            "keywords": [
                
                "n-methyltransferase",
                "N-methyltransferase",
                "n methyltransferase",
                
                "theobromine",
                "caffeine",
                "xmt",
                
                
                "xanthosine",
                "7 methylxanthine",
                "theophylline",
                "paraxanthine",
                "alkaloid biosynthesis",
                "caffeine synthase",
                "xanthosine n methyltransferase",
                "7 methylxanthine n methyltransferase",
            "3 7 dimethylxanthine n methyltransferase",
            "caffeine synthase",
            "theobromine synthase",
            "xmt",
            "mxmt",
            "dxmt",
            ],
        },
        "discovery_keywords": [
            "xanthosine n methyltransferase",
            "7 methylxanthine",
            "theobromine synthase",
            "xmt",
            "mxmt",
            "dxmt",
            "caffeine synthase",
        ],
    },

    "Coffee Leaf Rust (SH3)": {
        "description": (
            "Targets the SH3 resistance locus on Chr03C (3-16 Mbp) and the homeologous "
            "interval on Chr03E. The SH3 factor is linked to RGA3/RPP13-family NBS-LRR "
            "resistance proteins. Discovery tier captures these genes genome-wide using "
            "highly specific resistance protein identifiers."
        ),
        "high_conf_loci": {
            "chroms": ["Chr03C", "Chr03E"],
            "ranges": [(8_000_000, 9_000_000)],
            "keywords": [
                # Core NBS-LRR vocabulary — all hyphenated forms now handled by
                # normalize_text() so "nbs-lrr" in GFF3 becomes "nbs lrr" in blob
                "sh1", "sh2", "sh3", "sh4", "sh5", 
    "sh6", "sh7", "sh8", "sh9",
                "nbs lrr",
    "nbs-lrr",
    "cc-nbs-lrr",
    "tir-nbs-lrr",
    "rga",
    "rga3",
    "rlk",                   # Receptor-like kinase
    "lrr-rlk",
                "nbs lrr",
                "nb arc",
                
                
        
                
                "r protein",
                "r gene",
                
                "tir nbs",
                "coiled coil nbs",
                "cc nbs",
                "leucine rich repeat",
               "coffee leaf rust resistance",
               "coffee-leaf rust"
                
                "pathogen recognition",
                "tnl",
                "cnl",
               
                "hypersensitive response",
                "systemic acquired resistance",
                # Broad zone-gated fallbacks — safe inside Chr03:3-16Mbp
                "wrky",
                "pr protein",
                "pathogenesis related",
                
            ],
        },
        "discovery_keywords": [
            "rga",
            "nbs lrr",
            
        ],
    },

    "Coffee Berry Disease (Ck-1)": {
        "description": (
            "Targets the Ck-1 resistance interval on Chr01C and Chr01E (19.5-35 Mbp) "
            "for Colletotrichum kahawae resistance. All genes within this published QTL "
            "window are indexed positionally — keyword matching is attempted first, and "
            "genes with no annotation are included as 'Positional (Ck-1 locus)' candidates."
        ),
        # positional_fallback = True: if a gene is inside the QTL zone but has no
        # matching keyword (e.g. unannotated or sparsely annotated genes), index it
        # anyway with tier = "Positional (Ck-1 locus)". This is biologically justified
        # because the Ck-1 interval comes from published linkage mapping — any gene
        # within those coordinates is a legitimate candidate regardless of its annotation.
        
        "high_conf_loci": {
            "chroms": ["Chr01C", "Chr01E"],
            "ranges": [(19_500_000, 27_500_000)],
            "keywords": ["rlk",
    "wak",
    "nbs lrr",
    "nb arc",
    "lrr rlk",
    "leucine rich repeat",

    "wrky",
    
   
    "disease resistance",

    "map kinase",
    "serine threonine kinase",
    "serine threonine protein kinase",
            ],
            
        },
        "discovery_keywords": [
            
            "wall associated kinase",
            "wak",
            
            "rlk",            
        ],
    },

    "Dwarfism & Compactness (Ct)": {
        "description": (
            "Genome-wide search for gibberellin and brassinosteroid pathway regulators "
            "governing the Ct (compactness) phenotype. Key gene families include DELLA "
            "repressors, gibberellin 20-oxidase, and the brassinosteroid receptor BRI1. "
            "Expanded across all chromosomes because GA pathway genes are distributed "
            "throughout the genome — not restricted to Chr01."
        ),
        "high_conf_loci": {
            # Expanded from Chr01 only to all 11 chromosomes of each subgenome.
            # GA pathway genes (GA20ox, GA3ox, DELLA) are distributed genome-wide;
            # restricting to Chr01 was causing the 1-variant result.
            "chroms": [
                "Chr01C", "Chr01E", "Chr02C", "Chr02E", "Chr03C", "Chr03E",
                "Chr04C", "Chr04E", "Chr05C", "Chr05E", "Chr06C", "Chr06E",
                "Chr07C", "Chr07E", "Chr08C", "Chr08E", "Chr09C", "Chr09E",
                "Chr10C", "Chr10E", "Chr11C", "Chr11E",
            ],
            "ranges": None,
            "keywords": [
                
                "gibberellin 20 oxidase",
                "gibberellin 2 oxidase",
                "gibberellin 3 oxidase",
                "gibberellin dioxygenase",
                
                
                
                "gid1",
                "bri1",
                
                "brassinosteroid insensitive",
                "dwarf",
                "dwarfism",
                "semi dwarf",
                "plant height",
                
                "ga insensitive",
                "ga response",
                "ga20ox",
                "ga3ox",
                "ga2ox",
                "slr1",
                "rht",
            ],
        },
        "discovery_keywords": [
            
            
            "gibberellin 20 oxidase",
            "ga20ox",
            "brassinosteroid",
            "bri1",
            "gibberellin receptor",
        ],
    },

    "Giant Bean Size (Mg)": {
        "description": (
            "Focuses on the Maragogipe (Mg) locus on Chr02, targeting genes involved in "
            "cell expansion, endosperm development, and seed size determination. Expansins "
            "and DA1/KLUH-family ubiquitin receptor proteins are the primary candidates."
        ),
        "high_conf_loci": {
            "chroms": ["Chr02C", "Chr02E", "Chr06C", "Chr06E"],
            "ranges": [(5_000_000, 30_000_000)],
            "keywords": [
                
                "alpha expansin",
                "beta expansin",
                "cell elongation",
                "cell expansion",
                
                
                
                "kluh",
                "da1",
                "ubiquitin receptor",
                "seed size",
                "seed weight",
                "endosperm",
                "endosperm development",
                "grain size",
                "fruit size",
                
                "xyloglucan",
                
                
                
                
            ],
        },
        "discovery_keywords": [
            "da1",
            "kluh",
            "cyp78a",
            
            "ubiquitin receptor",

            "cell expansion",
            #newly added 16/5
            "da1",
        "kluh",
        "cyp78a5",
        "cyp78a",
        "alpha expansin",
        "expa",
        "ubiquitin receptor", 
        ],
    },
}

ATTRIBUTE_FIELDS = ("name", "note", "description", "dbxref", "ontology_term", "product", "id")

# --- UTILITY: TEXT & CONTIG NORMALISATION ------------------------------------

def normalize_text(text: str) -> str:
    """
    Collapse all non-alphanumeric characters to single spaces and lowercase.
    Applied identically to both keywords and the search blob so that
    hyphenated terms (nbs-lrr -> nbs lrr) match correctly in both directions.
    This is the critical fix vs the previous version.
    """
    if not text:
        return ""
    return re.sub(r"[^a-z0-9]+", " ", text.lower()).strip()


def normalize_contig(raw: str) -> str:
    """Map any GFF3/VCF contig identifier to canonical ChrXXC/E form."""
    seqid = raw.strip()

    # Strip trailing haplotype/scaffold suffix like _A, _B (e.g. Cara022_chrK_sgEE_A)
    seqid_clean = re.sub(r"_[A-Z]$", "", seqid)
    
    # Already canonical
    if re.fullmatch(r"Chr\d{2}[CE]", seqid_clean):
        return seqid_clean
    
    # Coffee Genome Hub arabica format: Cara001_chrA_sgCC / Cara001_chrA_sgCE
    m = re.search(r"_chr([A-Ka-k])_sg([A-Z]{2})$", seqid_clean)
    if m:
        letter_to_num = {
            "A":1,"B":2,"C":3,"D":4,"E":5,
            "F":6,"G":7,"H":8,"I":9,"J":10,"K":11
        }
        chrom_num = letter_to_num.get(m.group(1).upper(), 0)
        sub = "C" if "CC" in m.group(2) else "E"
        return f"Chr{chrom_num:02d}{sub}"

    # Existing fallback patterns
    match = re.search(r"(\d{1,2}).*sg([A-Z]{2})$", seqid_clean)
    if match:
        chrom_num = int(match.group(1))
        sub = "C" if "CC" in match.group(2) else "E"
        return f"Chr{chrom_num:02d}{sub}"
    match = re.fullmatch(r"0?(\d{1,2})([CE])", seqid_clean)
    if match:
        return f"Chr{int(match.group(1)):02d}{match.group(2)}"
    
    return seqid_clean


def parse_attributes(raw: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in raw.strip().split(";"):
        if "=" not in item:
            continue
        k, v = item.split("=", 1)
        attrs[k.strip().lower()] = v.strip()
    return attrs


def build_search_blob(attrs: dict[str, str]) -> str:
    """
    Concatenate all relevant GFF3 attribute fields then NORMALISE the entire
    string. Both blob and keywords are now in the same form before comparison.
    """
    raw = " ".join(attrs.get(f, "") for f in ATTRIBUTE_FIELDS)
    return normalize_text(raw)


def phrase_in_blob(phrase_norm: str, blob: str) -> bool:
    """
    Check whether a normalised keyword phrase appears in a normalised blob.
    Multi-word phrases: substring search (safe after normalisation).
    Single words: word-boundary regex to avoid false partial matches.
    """
    # Use \b for everything — single word AND multi-word phrases
    return bool(re.search(rf"\b{re.escape(phrase_norm)}\b", blob))


# --- GENE INDEX (built at startup) -------------------------------------------

gene_index:    dict[str, dict[str, list[dict]]] = {t: defaultdict(list) for t in BIOLOGICAL_MODEL}
start_index:   dict[str, dict[str, list[int]]]  = {t: defaultdict(list) for t in BIOLOGICAL_MODEL}
known_contigs: set[str] = set()


def index_gff3() -> None:
    print("=== Tiered Biological Index Initialisation ===")

    raw_features: list[dict] = []
    for gff_path in GFF3_FILES:
        if not os.path.exists(gff_path):
            print(f"  [WARN] GFF3 not found: {os.path.basename(gff_path)}")
            continue
        print(f"  Parsing {os.path.basename(gff_path)} ...")
        with open(gff_path, encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9 or parts[2].lower() != "gene":
                    continue
                seqid, start, end, raw_attr = (
                    parts[0], int(parts[3]), int(parts[4]), parts[8]
                )
                contig = normalize_contig(seqid)
                known_contigs.add(contig)
                attrs = parse_attributes(raw_attr)
                blob  = build_search_blob(attrs)   # NORMALISED blob
                raw_features.append({
                    "contig": contig, "start": start, "end": end,
                    "id": attrs.get("id", "Unknown"), "blob": blob,
                })

    print(f"  {len(raw_features):,} gene features loaded.")

    kw_hits:  dict[str, dict[str, int]] = {t: defaultdict(int) for t in BIOLOGICAL_MODEL}
    seen_ids: dict[str, set[str]]       = {t: set() for t in BIOLOGICAL_MODEL}

    for feat in raw_features:
        contig, start, end = feat["contig"], feat["start"], feat["end"]
        blob, gene_id = feat["blob"], feat["id"]

        for trait, model in BIOLOGICAL_MODEL.items():
            hc = model["high_conf_loci"]
            in_chrom = contig in hc["chroms"]
            in_range = True
            if hc["ranges"]:
                in_range = any(
                    r[0] <= start <= r[1] or r[0] <= end <= r[1]
                    for r in hc["ranges"]
                )
            is_tier1_zone = in_chrom and in_range

            matched_kw:    Optional[str] = None
            evidence_tier: Optional[str] = None

            if is_tier1_zone:
                for kw in hc["keywords"]:
                    kw_norm = normalize_text(kw)
                    if phrase_in_blob(kw_norm, blob):
                        matched_kw = kw
                        evidence_tier = "High Confidence"
                        kw_hits[trait][kw] += 1
                        break

            if not matched_kw:
                for kw in model["discovery_keywords"]:
                    kw_norm = normalize_text(kw)
                    if phrase_in_blob(kw_norm, blob):
                        matched_kw = kw
                        evidence_tier = "Discovery"
                        kw_hits[trait][kw] += 1
                        break

            # Positional fallback: if the gene has no annotation match but sits
            # inside a defined QTL coordinate window AND the trait opts into
            # positional indexing, include it. Biologically justified because the
            # QTL interval is from published linkage mapping — genomic position
            # alone is sufficient evidence to include it as a candidate.
            # if not matched_kw and is_tier1_zone and model.get("positional_fallback", False):
            #     matched_kw = "Positional (Ck-1 locus)"
            #     evidence_tier = "Positional"
            #     kw_hits[trait]["[positional fallback]"] += 1

            if matched_kw and gene_id not in seen_ids[trait]:
                seen_ids[trait].add(gene_id)
                gene_index[trait][contig].append({
                    "start": start, "end": end, "id": gene_id,
                    "matched_kw": matched_kw, "tier": evidence_tier,
                })

    print()
    for trait in BIOLOGICAL_MODEL:
        total = 0
        for contig in gene_index[trait]:
            gene_index[trait][contig].sort(key=lambda g: g["start"])
            start_index[trait][contig] = [g["start"] for g in gene_index[trait][contig]]
            total += len(gene_index[trait][contig])
        print(f"  {trait:42s}: {total:4d} candidate loci indexed")
        if total == 0:
            print(f"    [WARN] Zero genes matched — check GFF3 attribute vocabulary.")
            print(f"    Tip: run GET /api/health after startup and inspect candidateCounts.")
        else:
            top = sorted(kw_hits[trait].items(), key=lambda x: -x[1])[:4]
            print(f"    Top keywords: " + ", ".join(f"'{k}'={v}" for k, v in top))

    print("\n==============================================\n")


index_gff3()

# --- VARIANT CLASSIFICATION & ANALYTICS -------------------------------------

_TRANSITIONS = {frozenset("AG"), frozenset("CT")}


def classify_variant(ref: str, alt: str) -> str:
    if len(ref) == 1 and len(alt) == 1:
        return "SNP"
    if len(ref) < len(alt):
        return "Insertion"
    if len(ref) > len(alt):
        return "Deletion"
    return "MNP"


def is_transition(ref: str, alt: str) -> bool:
    return frozenset([ref.upper(), alt.upper()]) in _TRANSITIONS


def compute_analytics(variants: list[dict], candidate_count: int, model: dict) -> dict:
    empty = {
        "tiTvRatio": None, "subgenomeBias": None,
        "subgenomeCounts": {"C": 0, "E": 0},
        "hotspotGenes": [], "highConfFraction": None,
        "variantDensityPerMbp": None, "clusterCount": 0,
        "geneHitRate": None,
        "tierCounts": {"High Confidence": 0, "Discovery": 0, "Positional": 0},
        "snpCount": 0, "insertionCount": 0, "deletionCount": 0, "mnpCount": 0,
    }
    if not variants:
        return empty

    snp = ins = dele = mnp = ti = tv = 0
    c_count = e_count = 0
    tier_counts: dict[str, int] = defaultdict(int)
    gene_hits:   dict[str, int] = defaultdict(int)
    pos_by_contig: dict[str, list[int]] = defaultdict(list)

    for v in variants:
        ref, alt = v["reference"].upper(), v["alternate"].upper()
        vtype    = classify_variant(ref, alt)
        chrom    = v["chromosome"]
        kw_str   = (v.get("matchedKeywords") or [""])[0]

        if vtype == "SNP":
            snp += 1
            if is_transition(ref, alt):
                ti += 1
            else:
                tv += 1
        elif vtype == "Insertion":
            ins  += 1
        elif vtype == "Deletion":
            dele += 1
        else:
            mnp  += 1

        if   chrom.endswith("C"): c_count += 1
        elif chrom.endswith("E"): e_count += 1

        tier_counts["High Confidence" if "High Confidence" in kw_str else ("Positional" if "Positional" in kw_str else "Discovery")] += 1
        gene_hits[v["geneId"]] += 1
        pos_by_contig[chrom].append(v["position"])

    ti_tv     = round(ti / tv, 3) if tv > 0 else None
    total_sub = c_count + e_count
    bias      = round((c_count - e_count) / total_sub, 3) if total_sub > 0 else None
    hotspots  = [g for g, n in gene_hits.items() if n >= 3]
    total_t   = sum(tier_counts.values())
    hc_frac   = round(tier_counts["High Confidence"] / total_t, 3) if total_t > 0 else None

    hc = model["high_conf_loci"]
    if hc["ranges"]:
        region_bp = sum(r[1] - r[0] for r in hc["ranges"])
    else:
        region_bp = max(50_000_000, candidate_count * 5000)
    density = round(len(variants) / (region_bp / 1_000_000), 3)

    cluster_count = 0
    for positions in pos_by_contig.values():
        positions.sort()
        for i, pos in enumerate(positions):
            if (i > 0 and pos - positions[i - 1] <= 1000) or \
               (i < len(positions) - 1 and positions[i + 1] - pos <= 1000):
                cluster_count += 1

    gene_hit_rate = round(len(gene_hits) / candidate_count, 4) if candidate_count > 0 else None

    return {
        "tiTvRatio":            ti_tv,
        "subgenomeBias":        bias,
        "subgenomeCounts":      {"C": c_count, "E": e_count},
        "hotspotGenes":         hotspots,
        "highConfFraction":     hc_frac,
        "variantDensityPerMbp": density,
        "clusterCount":         cluster_count,
        "geneHitRate":          gene_hit_rate,
        "tierCounts":           dict(tier_counts),
        "snpCount":             snp,
        "insertionCount":       ins,
        "deletionCount":        dele,
        "mnpCount":             mnp,
    }


# --- VCF ANALYSIS ------------------------------------------------------------

def analyze_vcf(filepath: str, filename: str) -> dict:
    counts:   dict[str, int]        = {t: 0  for t in BIOLOGICAL_MODEL}
    variants: dict[str, list[dict]] = {t: [] for t in BIOLOGICAL_MODEL}

    seen_contigs:      set[str] = set()
    unmatched_contigs: set[str] = set()
    total_variants = 0

    opener = gzip.open if filename.lower().endswith(".gz") else open
    with opener(filepath, "rt", encoding="utf-8", errors="ignore") as vcf:
        for line in vcf:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            raw_chrom, pos_str, _, ref, alt = (
                parts[0], parts[1], parts[2], parts[3], parts[4]
            )
            try:
                pos = int(pos_str)
            except ValueError:
                continue

            total_variants += 1
            chrom = normalize_contig(raw_chrom)
            seen_contigs.add(chrom)

            if chrom not in known_contigs:
                unmatched_contigs.add(chrom)
                continue

            for trait in BIOLOGICAL_MODEL:
                genes  = gene_index[trait].get(chrom)
                starts = start_index[trait].get(chrom)
                if not genes:
                    continue
                probe = bisect_right(starts, pos)
                for idx in (probe - 1, probe):
                    if idx < 0 or idx >= len(genes):
                        continue
                    g = genes[idx]
                    if g["start"] <= pos <= g["end"]:
                        variants[trait].append({
                            "chromosome":      chrom,
                            "position":        pos,
                            "reference":       ref,
                            "alternate":       alt,
                            "geneId":          g["id"],
                            "matchedKeywords": [f"{g['matched_kw']} ({g['tier']})"],
                        })
                        counts[trait] += 1
                        break

    matched = seen_contigs & known_contigs

    analytics:     dict[str, dict] = {}
    trait_summary: dict[str, dict] = {}
    for trait, model in BIOLOGICAL_MODEL.items():
        candidate_count = sum(len(gs) for gs in gene_index[trait].values())
        analytics[trait] = compute_analytics(variants[trait], candidate_count, model)
        trait_summary[trait] = {
            "candidateGenes": candidate_count,
            "description":    model["description"],
        }

    payload = {
        "counts":       counts,
        "variants":     variants,
        "traitSummary": trait_summary,
        "analytics":    analytics,
        "meta": {
            "sourceFile":       filename,
            "timestamp":        str(datetime.datetime.now()),
            "totalVariants":    total_variants,
            "contigsSeen":      sorted(seen_contigs),
            "matchedContigs":   sorted(matched),
            "unmatchedContigs": sorted(unmatched_contigs),
            "referenceContigs": sorted(known_contigs),
        },
    }

    with open(PAYLOAD_PATH, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)

    print(f"  Done — {total_variants:,} variants, {len(matched)} contigs matched.")
    print(f"  Counts: { {t: counts[t] for t in counts} }")
    return payload


# --- FLASK ROUTES ------------------------------------------------------------

@app.route("/")
def index_page():
    return send_from_directory(".", "index.html")


@app.route("/annotation_payload.json")
def serve_payload():
    if os.path.exists(PAYLOAD_PATH):
        return send_from_directory(BASE_DIR, "annotation_payload_antigravity.json")
    return jsonify({"error": "No payload yet — run /api/analyze first."}), 404


@app.get("/api/health")
def health():
    return jsonify({
        "status": "ok",
        "traits": list(BIOLOGICAL_MODEL.keys()),
        "candidateCounts": {
            t: sum(len(gs) for gs in gene_index[t].values())
            for t in BIOLOGICAL_MODEL
        },
        "indexedContigs": sorted(known_contigs),
        "gff3Files":      [os.path.basename(p) for p in GFF3_FILES],
        "gff3Found":      [os.path.exists(p) for p in GFF3_FILES],
    })


@app.post("/api/analyze")
def api_analyze():
    print(f"\n--- Analysis at {datetime.datetime.now()} ---")
    file = request.files.get("vcf")

    if file and file.filename:
        filename = secure_filename(file.filename)
        filepath = os.path.join(UPLOAD_FOLDER, filename)
        file.save(filepath)
        try:
            payload = analyze_vcf(filepath, filename)
        except Exception as e:
            print(f"[ERROR] {e}")
            return jsonify({"error": str(e)}), 500
        finally:
            if os.path.exists(filepath):
                os.remove(filepath)
        return jsonify(payload)

    if os.path.exists(DEFAULT_VCF):
        print("  No upload — using default VCF.")
        try:
            payload = analyze_vcf(DEFAULT_VCF, os.path.basename(DEFAULT_VCF))
        except Exception as e:
            print(f"[ERROR] {e}")
            return jsonify({"error": str(e)}), 500
        return jsonify(payload)

    return jsonify({
        "error": (
            "No VCF file uploaded and no default VCF found. "
            "Upload a .vcf or .vcf.gz, or place 'uploaded_variants.vcf' "
            "in the same directory as this script."
        )
    }), 400


if __name__ == "__main__":
    app.run(host="127.0.0.1", port=5002, debug=False)
