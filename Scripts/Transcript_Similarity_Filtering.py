#!/usr/bin/env python3
import requests
import subprocess
import tempfile
import os
import gzip
import re
import pandas as pd
from collections import defaultdict
from Bio import SeqIO
import random
import time
import concurrent.futures
import requests
import threading


# cDNA + ncRNA CONFIG
files = [
    "Homo_sapiens.GRCh38.cdna.all.fa", 
    "Homo_sapiens.GRCh38.ncrna.fa"
]
records = []

for f in files:
    for record in SeqIO.parse(f, "fasta"):
        records.append(record)

SeqIO.write(records, "sequence_data.fa", "fasta")


# CONFIG
GTF_FILE = "Homo_sapiens.GRCh38.115.gtf.gz"
CDNA_FASTA = "sequence_data.fa"
PROTEIN_FASTA = "Homo_sapiens.GRCh38.pep.all.fa"
UNIPROT_FASTA = "Homo_sapiens.GRCh38.115.uniprot.tsv"

APPRIS_FILE = "appris_data.appris.txt"
GENCODE_GTF = "gencode_transcripts_only.gtf"

OUTPUT_RETAINED = "retained_isoforms_global.csv"
OUTPUT_EXCLUDED = "excluded_isoforms_global.csv"

LENGTH_THRESHOLD = 50
IDENTITY_THRESHOLD = 96.0  # percent
HEADERS = {"Content-Type": "application/json"}

def parse_gtf(gtf_file):
    records = []
    with gzip.open(gtf_file, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs = fields
            attrs_dict = dict(re.findall(r'(\S+) "([^"]+)"', attrs))
            
            gene_id = attrs_dict.get("gene_id", "NA").split('.')[0]  # remove version
            transcript_id = attrs_dict.get("transcript_id", "NA").split('.')[0]  # remove version
            
            tags = attrs_dict.get("tag", "")
            records.append({
                "gene_id": gene_id,
                "gene_name": attrs_dict.get("gene_name"),
                "gene_biotype": attrs_dict.get("gene_biotype"),
                "transcript_id": transcript_id,
                "transcript_name": attrs_dict.get("transcript_name"),
                "transcript_biotype": attrs_dict.get("transcript_biotype"),
                "tags": tags.split(",") if tags else [],
                "tsl": attrs_dict.get("transcript_support_level", "NA"),
                "ccds_id": attrs_dict.get("ccds_id", "NA"),
                "appris": attrs_dict.get("appris", "NA"),
                "chrom": chrom,
                "strand": strand,
            })
    return pd.DataFrame(records)

_gtf_df = parse_gtf(GTF_FILE)

def clean_tid(tid):
    return tid.split(".")[0] if tid and tid != "NA" else "NA"


# cDNA FASTA (version-stripped)
_cdna_dict = {rec.id.split('.')[0]: str(rec.seq) for rec in SeqIO.parse(CDNA_FASTA, "fasta")}


# Protein FASTA (version-stripped)
_pep_dict = {rec.id.split()[0].split('.')[0]: len(rec.seq) for rec in SeqIO.parse(PROTEIN_FASTA, "fasta")}


print(f"Loaded {_gtf_df.shape[0]} transcripts from GTF.")

# FETCH FUNCTIONS
def fetch_transcripts(gene_id, appris_data, gencode_data):
    """Local version of fetch_transcripts — reads from preloaded GTF/FASTA."""
    subset = _gtf_df[_gtf_df["gene_id"] == gene_id]
    transcripts = []

    for _, row in subset.iterrows():
        tid = row["transcript_id"].strip()
        seq = _cdna_dict.get(tid) or _cdna_dict.get(tid.split('.')[0], "")
        length = len(seq)
        pep_len = _pep_dict.get(tid, "NA")
        tags = row["tags"]

        transcripts.append({
            "Transcript_ID": row["transcript_id"],
            "Transcript_Name": row["transcript_name"],
            "Transcript_Biotype": row["transcript_biotype"],
            "Nucleotide_Length": length,
            "Protein_Length": pep_len,
            "APPRIS": row["appris"],
            "APPRIS_Combined": row["appris"],
            "MANE": "MANE_Select" in tags,
            "GENCODE": ";".join(tags),
            "CCDS": row["ccds_id"],
            "UniProt_ID": "NA",
            "TRIFID": 0.0,
            "TSL": row["tsl"],
            "seq": seq
        })
    return transcripts


def fetch_ccds(transcript_id):
    """Local CCDS fetch from GTF (already included in transcript record)"""
    entry = _gtf_df[_gtf_df["transcript_id"] == transcript_id]
    return entry["ccds_id"].iloc[0] if not entry.empty else "NA"


def fetch_protein_length(protein_id):
    """Return protein length from preloaded FASTA dict"""
    pid = protein_id.split('.')[0]
    return _pep_dict.get(pid, "NA")


def fetch_uniprot_from_protein(protein_id):
    """No direct UniProt mapping offline — placeholder"""
    return "NA"

# FUNCTIONS
def fetch_gene_info(gene_id):
    """Fetch gene info locally from GTF"""
    rows = _gtf_df[_gtf_df["gene_id"] == gene_id]
    if rows.empty:
        return {"Gene_ID": gene_id, "Gene_Name": "NA", "Gene_Biotype": "NA"}
    gene_name = rows.iloc[0]["gene_name"] if "gene_name" in rows.columns else "NA"
    gene_biotype = rows.iloc[0]["gene_biotype"] if "gene_biotype" in rows.columns else "NA"
    return {"Gene_ID": gene_id, "Gene_Name": gene_name, "Gene_Biotype": gene_biotype}

    
# APPRIS metadata
def load_appris(appris_file):
    df = pd.read_csv(appris_file, sep="\t", dtype=str)
    appris = {}
    for _, row in df.iterrows():
        tid = str(row["Transcript ID"]).split('.')[0]  # remove version
        ann = str(row.get("APPRIS Annotation", "NA")).strip().lower()
        
        if ann.startswith("appris_principal_"):
            num = ann.split("_")[-1]
            appris_ann = f"PRINCIPAL:{num}"
        elif ann.startswith("appris_alternative_"):
            num = ann.split("_")[-1]
            appris_ann = f"ALTERNATIVE:{num}"
        else:
            appris_ann = ann.upper() if ann != "na" else "NA"
        
        appris[tid] = {
            "APPRIS": appris_ann,
            "MANE": str(row.get("MANE", "")).strip().upper() == "MANE_SELECT",
            "TRIFID": float(row.get("Trifid Score", 0.0)) if row.get("Trifid Score") not in (None, "") else 0.0,
            "TSL": row.get("Transcript support level", "NA")
        }
    return appris

appris_data = load_appris(APPRIS_FILE)

# GENCODE metadata
def load_gencode(gencode_gtf):
    keep_tags = {
        "appris_principal_1","appris_principal_2","appris_principal_3",
        "appris_principal_4","appris_principal_5","appris_alternative_1","appris_alternative_2",
        "MANE_Select","GENCODE_Primary","Ensembl_canonical","basic",
        "CCDS","CAGE_supported_TSS","454_RNA_Seq_supported",
        "RNA_Seq_supported_only","RNA_Seq_supported_partial",
        "nested_454_RNA_Seq_supported","RP_supported_TIS","exp_conf",
        "dotter_confirmed","alternative_3_UTR","alternative_5_UTR"
    }
    gencode = defaultdict(lambda: {"tags": [], "TSL": "NA"})
    
    with open(gencode_gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[2] != "transcript":
                continue
            tid = None
            tags = []
            tsl_value = "NA"
            for kv in parts[8].split(";"):
                kv = kv.strip()
                if not kv: continue
                pieces = kv.split()
                if len(pieces) >= 2:
                    key = pieces[0]
                    value = pieces[1].strip('"')
                    if key == "transcript_id":
                        tid = value.split('.')[0]  # remove version
                    elif key == "tag" and value in keep_tags:
                        tags.append(value)
                    elif key == "transcript_support_level":
                        tsl_value = value
            if tid:
                if tags:
                    gencode[tid]["tags"].extend(tags)
                if tsl_value != "NA":
                    gencode[tid]["TSL"] = tsl_value
    return gencode

gencode_data = load_gencode(GENCODE_GTF)


def fetch_transcripts_from_gtf(gene_id, appris_data, gencode_data, transcript_seq_dict=_cdna_dict):
    """
    Fetch transcript metadata + sequence from local GTF (_gtf_df) and merge with APPRIS + GENCODE
    
    transcript_seq_dict: optional dict {Transcript_ID: seq}, preloaded from FASTA
    """
    # Filter GTF for this gene
    gene_df = _gtf_df[_gtf_df["gene_id"] == gene_id]
    if gene_df.empty:
        return []

    transcripts = []
    # group by transcript_id
    for tid, tgroup in gene_df.groupby("transcript_id"):
        tid_clean = tid.strip()
        tver = tid_clean.split(".")[0]

        # Metadata from APPRIS + GENCODE
        meta_appris = appris_data.get(tver, {})
        meta_gencode = gencode_data.get(tver, {})

        # Merge TSL
        tsl_appris = meta_appris.get("TSL", "NA")
        tsl_gencode = meta_gencode.get("TSL", "NA")
        tsl_final = tsl_gencode if tsl_gencode not in (None, "", "NA") else tsl_appris

        # CCDS
        ccds_id = tgroup["ccds"].iloc[0] if "ccds" in tgroup.columns else "NA"
        if not ccds_id or ccds_id in ("", "NA"):
            ccds_id = "NA"

        # Protein length
        protein_len = tgroup["protein_length"].iloc[0] if "protein_length" in tgroup.columns else 0

        # UniProt ID
        uniprot_id = tgroup["uniprot_id"].iloc[0] if "uniprot_id" in tgroup.columns else "NA"

        # APPRIS combined
        appris_ann = meta_appris.get("APPRIS", "NA")
        gencode_tags = meta_gencode.get("tags", [])
        gencode_appris = [tag for tag in gencode_tags if tag.startswith("appris_")]

        combined_appris = set()
        if appris_ann != "NA":
            combined_appris.add(appris_ann)
        combined_appris.update(gencode_appris)
        appris_combined = ";".join(sorted(combined_appris)) if combined_appris else "NA"

        # Sequence
        seq_rec = transcript_seq_dict.get(tid_clean)
        seq = str(seq_rec) if seq_rec else ""
        if transcript_seq_dict and tid_clean in transcript_seq_dict:
            seq = transcript_seq_dict[tid_clean]
        length = len(seq)

        # Build transcript dictionary
        transcripts.append({
            "Transcript_ID": tver,
            "Transcript_Name": tgroup["transcript_name"].iloc[0] if "transcript_name" in tgroup.columns else tid_clean,
            "Transcript_Biotype": tgroup["transcript_biotype"].iloc[0] if "transcript_biotype" in tgroup.columns else "NA",
            "Nucleotide_Length": length,
            "Protein_Length": protein_len,
            "APPRIS": appris_ann,
            "APPRIS_Combined": appris_combined,
            "MANE": meta_appris.get("MANE", False),
            "GENCODE": ";".join(gencode_tags),
            "CCDS": ccds_id,
            "UniProt_ID": uniprot_id,
            "TRIFID": meta_appris.get("TRIFID", 0.0),
            "TSL": tsl_final,
            "TSL_APPRIS": tsl_appris,
            "TSL_GENCODE": tsl_gencode,
            "seq": seq
        })

    return transcripts

def group_by_length(transcripts, tolerance=50):
    """Group transcripts within +/- tolerance bp length"""
    transcripts = sorted(transcripts, key=lambda x: x["Nucleotide_Length"])
    unique, groups, used = [], [], set()
    for i, t in enumerate(transcripts):
        if t["Transcript_ID"] in used:
            continue
        group = [t]
        for j in range(i+1, len(transcripts)):
            if abs(transcripts[j]["Nucleotide_Length"] - t["Nucleotide_Length"]) <= tolerance:
                group.append(transcripts[j])
                used.add(transcripts[j]["Transcript_ID"])
        if len(group) == 1:
            unique.append(t)
        else:
            groups.append(group)
        for g in group:
            used.add(g["Transcript_ID"])
    return unique, groups

def run_mafft(sequences):
    """Run MAFFT and compute pairwise identity ignoring gaps. Returns pairwise percent id."""
    import subprocess
    from Bio import SeqIO
    import tempfile, os
    from collections import defaultdict

    # write sequences to temp fasta
    with tempfile.NamedTemporaryFile("w+", delete=False, suffix=".fa") as fasta:
        for s in sequences:
            fasta.write(f">{s['Transcript_ID']}\n{s['seq']}\n")
        fasta_name = fasta.name

    aln_name = fasta_name + ".aln"
    try:
        # run MAFFT and capture stdout to aln_name
        with open(aln_name, "w") as outfh:
            subprocess.run(["mafft", "--auto", fasta_name], stdout=outfh, stderr=subprocess.DEVNULL, check=True)
        alignment = list(SeqIO.parse(aln_name, "fasta"))
    except subprocess.CalledProcessError as e:
        print(f"Warning: MAFFT failed for cluster ({e}); skipping cluster.")
        alignment = []
    finally:
        # cleanup input fasta and aln if present
        try:
            os.remove(fasta_name)
        except OSError:
            pass
        try:
            if os.path.exists(aln_name):
                os.remove(aln_name)
        except OSError:
            pass

    scores = defaultdict(dict)
    if not alignment:
        return scores

    for i in range(len(alignment)):
        for j in range(i+1, len(alignment)):
            a, b = str(alignment[i].seq), str(alignment[j].seq)
            matches = compared = 0
            for x, y in zip(a, b):
                if x == "-" or y == "-":
                    continue
                compared += 1
                if x == y:
                    matches += 1
            pid = 100.0 * matches / compared if compared else 0.0
            scores[alignment[i].id][alignment[j].id] = pid
            scores[alignment[j].id][alignment[i].id] = pid
    return scores

# GENCODE priority dictionary & helper
gencode_priority = {
    "MANE_Select": 1,
    "GENCODE_Primary": 1,
    "Ensembl_canonical": 1,
    "454_RNA_Seq_supported": 2,
    "nested_454_RNA_Seq_supported": 2,
    "RNA_Seq_supported_only": 2,
    "RNA_Seq_supported_partial": 2,
    "RP_supported_TIS": 2,
    "exp_conf": 2,
    "CCDS": 3,
    "CAGE_supported_TSS": 3,
    "basic": 3,
    "appris_principal_1": 3,
    "appris_principal_2": 3,
    "appris_principal_3": 3,
    "appris_principal_4": 3,
    "appris_principal_5": 3,
    "appris_alternative_1": 4,
    "appris_alternative_2": 4,
    "dotter_confirmed": 4,
    "alternative_3_UTR": 4,
    "alternative_5_UTR": 4,
}

def best_gencode_score(tags: str):
    """Return best GENCODE tag rank and how many times it occurs"""
    if not tags or pd.isna(tags):
        return (999, 0)  # no tags
    scores = []
    for tag in str(tags).replace(";", ",").split(","):
        tag = tag.strip()
        if tag in gencode_priority:
            scores.append(gencode_priority[tag])
    if not scores:
        return (999, 0)
    best = min(scores)
    return (best, sum(1 for s in scores if s == best))


# Small helper ranking functions
def appris_priority(appris_value: str) -> int:
    """Assign numeric rank to APPRIS annotations. Lower = better"""
    if not isinstance(appris_value, str):
        return 99
    appris_value = appris_value.strip().upper()
    if appris_value.startswith("PRINCIPAL:"):
        try:
            return int(appris_value.split(":")[1])
        except:
            return 99
    if appris_value.startswith("ALTERNATIVE:"):
        try:
            return 10 + int(appris_value.split(":")[1])
        except:
            return 20
    return 99

def tsl_priority(tsl_value: str) -> int:
    """Assign numeric rank to TSL. Lower = better"""
    tsl_order = {"1": 1, "2": 2, "3": 3, "4": 4, "5": 5}
    if tsl_value is None:
        return 99
    return tsl_order.get(str(tsl_value).strip(), 99)


# Filtering function (canonical tie-breaker) — currently unused; filter_strong_metadata is called directly
def filter_metadata_with_gencode(cluster):
    """
    Input: list of dicts or pandas DataFrame for a cluster of similar transcripts
    Output: (retained_df, excluded_df)

    Logic:
    - Strong metadata (CCDS, UniProt, MANE, APPRIS principal, TRIFID >0) is prioritized
    - Weak transcripts (no metadata at all) are excluded if tied
    - Tie-breaking order for strong metadata:
        1. MANE
        2. CCDS uniqueness
        3. Unique UniProt ID
        4. APPRIS priority
        5. TRIFID score
        6. TSL
        7. GENCODE tags
        8. Sequence length (only if no metadata)
    """

    if isinstance(cluster, list):
        cluster = pd.DataFrame(cluster)
    elif not isinstance(cluster, pd.DataFrame):
        raise TypeError(f"Expected DataFrame or list, got {type(cluster)}")

    retained, excluded = [], []
    cluster = cluster.fillna("NA").replace({pd.NA: "NA", None: "NA"})

    def drop_and_mark_winner(df, winner_idx, winner_reason, loser_reason_prefix="Excluded"):
        df.loc[winner_idx, "Reason"] = winner_reason
        winner_row = df.loc[winner_idx]
        losers = df.drop(winner_idx)
        for idx, _ in losers.iterrows():
            prev_reason = df.loc[idx]["Reason"] if "Reason" in df.columns and pd.notna(df.loc[idx]["Reason"]) else "lost tie-break"
            df.loc[idx, "Reason"] = f"{loser_reason_prefix}: {prev_reason}"
            excluded.append(df.loc[idx].to_dict())
        retained.append(winner_row.to_dict())
        return pd.DataFrame(retained), pd.DataFrame(excluded)


# Module-level function: determine if a transcript has strong metadata
def is_strong_metadata(t):
    """Return True if transcript has strong metadata (MANE, CCDS, UniProt, APPRIS principal, TRIFID >0)"""
    try:
        trifid_val = float(t.get("TRIFID", 0.0))
    except (ValueError, TypeError):
        trifid_val = 0.0

    # Correct UniProt check: allow IDs with periods
    up_field = str(t.get("UniProt_ID", "")).strip()
    has_uniprot = any(uid.strip() not in ("", "NA", "nan") for uid in up_field.replace(",", ";").split(";"))

    appris_field = str(t.get("APPRIS", "")).upper()
    ccds_field = t.get("CCDS", "")

    return (
        t.get("MANE") is True
        or (ccds_field not in ("NA", "", None))
        or has_uniprot
        or appris_field.startswith("PRINCIPAL")
        or trifid_val > 0
    )


# Module-level function: metadata-based filtering within a cluster
def filter_strong_metadata(cluster):
    """
    Input: transcript cluster 
    Output: (retained_list, excluded_list)

    Logic:
    - For strong transcripts:
        * Identify and retain the best one in each CCDS or UniProt group
        * All other strong transcripts that were compared are NOT excluded;
          they are re-ranked together with weak transcripts using the same metadata hierarchy
    - If no strong transcripts exist:
        * Rank all by metadata hierarchy + sequence length; remove only the weakest, shortest transcript
    """

    # Normalize input
    if isinstance(cluster, list):
        df = pd.DataFrame(cluster)
    else:
        df = cluster.copy()

    df = df.fillna("NA").replace({pd.NA: "NA", None: "NA"})

    # helper ranking columns
    df["MANE_rank"] = df["MANE"].apply(lambda x: 0 if x is True else 1)
    df["APPRIS_rank"] = df["APPRIS"].astype(str).apply(appris_priority)
    df["TRIFID_num"] = pd.to_numeric(df["TRIFID"], errors="coerce").fillna(0.0)
    df["TSL_rank"] = df["TSL"].astype(str).apply(tsl_priority)
    df["GENCODE_rank"] = df["GENCODE"].astype(str).apply(lambda x: best_gencode_score(x)[0])
    df["Length_rank"] = -df["Nucleotide_Length"].astype(float)

    # function to get best row by metadata key
    def best_index_by_metadata(subdf):
        keys = subdf.apply(
            lambda r: (
                int(r["MANE_rank"]),
                int(r["APPRIS_rank"]),
                -float(r["TRIFID_num"]),
                int(r["TSL_rank"]),
                int(r["GENCODE_rank"]),
                float(r["Length_rank"]),
            ),
            axis=1,
        )
        return keys.idxmin()

    # identify strong and weak transcripts
    strong_mask = df.apply(is_strong_metadata, axis=1)
    strong_df = df[strong_mask].copy()
    weak_df = df[~strong_mask].copy()

    retained_rows, excluded_rows = [], []
    retained_ids = set()
    compared_strong_ids = set()

    # Step 1: Handle strong transcripts (retain best per CCDS/UniProt)
    if not strong_df.empty:
        processed = set()

        # (a) CCDS grouping
        for ccds_val, group in strong_df[~strong_df["CCDS"].isin(["", "NA"])].groupby("CCDS"):
            best_idx = best_index_by_metadata(group)
            best_row = group.loc[best_idx].to_dict()
            best_row["Reason"] = "Retained: best strong transcript within CCDS group"
            retained_rows.append(best_row)
            retained_ids.add(best_row["Transcript_ID"])
            compared_strong_ids.update(group["Transcript_ID"].tolist())
            processed.update(group["Transcript_ID"].tolist())

        # (b) UniProt overlap grouping for remaining strong transcripts
        remaining = strong_df[~strong_df["Transcript_ID"].isin(processed)].copy()

        def normalize_uniprot_set(s):
            if not isinstance(s, str) or s in ("", "NA"):
                return set()
            parts = [x.strip() for x in s.replace(",", ";").split(";") if x.strip() and x.strip().upper() != "NA"]
            return set(parts)

        tid_to_up = {row["Transcript_ID"]: normalize_uniprot_set(row["UniProt_ID"]) for _, row in remaining.iterrows()}
        tids = list(tid_to_up.keys())

        seen = set()
        for tid in tids:
            if tid in seen:
                continue
            comp = {tid}
            stack = [tid]
            while stack:
                cur = stack.pop()
                seen.add(cur)
                for other in tids:
                    if other not in comp and tid_to_up[cur] & tid_to_up[other]:
                        comp.add(other)
                        stack.append(other)
            comp_rows = remaining[remaining["Transcript_ID"].isin(comp)]
            if len(comp_rows) == 1:
                r = comp_rows.iloc[0].to_dict()
                if r["Transcript_ID"] not in retained_ids:
                    r["Reason"] = "Retained: strong metadata (unique UniProt)"
                    retained_rows.append(r)
                    retained_ids.add(r["Transcript_ID"])
                continue
            best_idx = best_index_by_metadata(comp_rows)
            best_row = comp_rows.loc[best_idx].to_dict()
            best_row["Reason"] = "Retained: best strong transcript within UniProt group"
            retained_rows.append(best_row)
            retained_ids.add(best_row["Transcript_ID"])
            compared_strong_ids.update(comp_rows["Transcript_ID"].tolist())

        # (c) Strong transcripts not in any CCDS/UniProt group
        remaining = strong_df[~strong_df["Transcript_ID"].isin(processed | set(tids))].copy()
        for _, row in remaining.iterrows():
            r = row.to_dict()
            if r["Transcript_ID"] not in retained_ids:
                r["Reason"] = "Retained: strong metadata (unique)"
                retained_rows.append(r)
                retained_ids.add(r["Transcript_ID"])
            compared_strong_ids.add(r["Transcript_ID"])

        # (d) Any strong transcripts not retained (but compared) -> move to re-ranking pool
        re_rank_df = df[df["Transcript_ID"].isin(compared_strong_ids - retained_ids)].copy()
        weak_df = pd.concat([weak_df, re_rank_df], ignore_index=True)

    # Step 2: Rank remaining transcripts
    if not weak_df.empty:
        rank_df = weak_df.copy()

        rank_df["CCDS_rank"] = rank_df["CCDS"].apply(lambda x: 0 if (isinstance(x, str) and x not in ("", "NA")) else 1)
        rank_df["UniProt_rank"] = rank_df["UniProt_ID"].apply(
            lambda x: 0 if (isinstance(x, str) and any(u.strip() not in ("", "NA", "nan")
                                                       for u in str(x).replace(",", ";").split(";"))) else 1
        )

        rank_df["sort_key"] = rank_df.apply(
            lambda r: (
                int(r["MANE_rank"]),
                int(r["CCDS_rank"]),
                int(r["UniProt_rank"]),
                int(r["APPRIS_rank"]),
                -float(r["TRIFID_num"]),
                int(r["TSL_rank"]),
                int(r["GENCODE_rank"]),
                float(r["Length_rank"]),
            ), axis=1
        )

        rank_df = rank_df.sort_values(by="sort_key", ascending=True).reset_index(drop=True)
        max_key = max(rank_df["sort_key"].tolist())
        weakest_df = rank_df[rank_df["sort_key"].apply(lambda k: k == max_key)].copy()
        min_len = weakest_df["Nucleotide_Length"].astype(float).min()
        to_exclude = weakest_df[weakest_df["Nucleotide_Length"].astype(float) == min_len].iloc[[0]]
        exclude_ids = set(to_exclude["Transcript_ID"].tolist())

        for _, row in rank_df.iterrows():
            r = row.to_dict()
            if row["Transcript_ID"] in exclude_ids:
                r["Reason"] = "Excluded: weakest-ranking, shortest transcript"
                excluded_rows.append(r)
            else:
                r["Reason"] = "Retained: ranked by metadata hierarchy"
                retained_rows.append(r)
                retained_ids.add(r["Transcript_ID"])

    # Final cleanup
    retained_ids = {r["Transcript_ID"] for r in retained_rows}
    excluded_rows = [r for r in excluded_rows if r["Transcript_ID"] not in retained_ids]

    return retained_rows, excluded_rows


def process_gene(gene_id, appris_data, gencode_data):
    """Process transcripts: retain unique-length transcripts, then filter by biotype + sequence similarity"""
    transcripts = fetch_transcripts_from_gtf(gene_id, appris_data, gencode_data)
    if not transcripts:
        print(f"Warning: Gene {gene_id} has no transcripts")
        return [], []

    retained, excluded = [], []

    # Step 1: Cluster transcripts by length ± LENGTH_THRESHOLD
    unique, groups = group_by_length(transcripts, tolerance=LENGTH_THRESHOLD)

    # Step 2: Retain all unique-length transcripts (regardless of metadata)
    for t in unique:
        if is_strong_metadata(t):
            t["Reason"] = "Retained: unique length + strong metadata"
        else:
            t["Reason"] = "Retained: unique length (no strong metadata)"
        t["Similarity_Score"] = 0.0
        retained.append(t)

    # Step 3: Process grouped transcripts
    for group in groups:
        retained_ids = {r["Transcript_ID"] for r in retained}
        group = [t for t in group if t["Transcript_ID"] not in retained_ids]
        if not group:
            continue

        # Step 3a: Sub-cluster by biotype
        biotype_groups = defaultdict(list)
        for t in group:
            biotype = t.get("Transcript_Biotype", "NA")
            biotype_groups[biotype].append(t)

        # Step 3b: Process each biotype group separately
        for biotype, seq_group in biotype_groups.items():

            # Remove missing-sequence transcripts first
            seq_group = [t for t in seq_group if t.get("seq")]
            for t in biotype_groups[biotype]:
                if not t.get("seq"):
                    t["Reason"] = f"Excluded: no sequence ({biotype})"
                    t["Similarity_Score"] = 0.0
                    excluded.append(t)
            if len(seq_group) == 0:
                continue

            # Single transcript: retain directly
            if len(seq_group) == 1:
                t = seq_group[0]
                if is_strong_metadata(t):
                    t["Reason"] = f"Retained: unique within biotype ({biotype}) + strong metadata"
                else:
                    t["Reason"] = f"Retained: unique within biotype ({biotype})"
                t["Similarity_Score"] = 0.0
                retained.append(t)
                continue

            # Step 3c: Align and compute pairwise identity
            scores = run_mafft(seq_group)

            # Compute each transcript's highest similarity to any other
            for t in seq_group:
                tid = t["Transcript_ID"]
                highest = max(scores.get(tid, {}).values(), default=0.0)
                t["Similarity_Score"] = round(highest, 5)

            # Step 3d: Separate low- and high-similarity transcripts
            low_sim, high_sim = [], []
            for t in seq_group:
                if t["Similarity_Score"] < IDENTITY_THRESHOLD:
                    low_sim.append(t)
                else:
                    high_sim.append(t)

            # Automatically retain below-threshold transcripts
            for t in low_sim:
                t["Reason"] = f"Retained: low similarity (<{IDENTITY_THRESHOLD}%) ({biotype})"
            retained.extend(low_sim)

            # Step 3e: Apply metadata-based filtering to >= threshold transcripts
            if high_sim:
                r, e = filter_strong_metadata(high_sim)
                for x in r:
                    x["Reason"] += f" ({biotype})"
                for x in e:
                    x["Reason"] += f" ({biotype})"
                retained.extend(r)
                excluded.extend(e)

    # Step 4: Any transcript not yet classified is excluded
    retained_ids = {t["Transcript_ID"] for t in retained}
    excluded_ids = {t["Transcript_ID"] for t in excluded}
    for t in transcripts:
        tid = t["Transcript_ID"]
        if tid not in retained_ids and tid not in excluded_ids:
            if is_strong_metadata(t):
                t["Reason"] = "Retained: strong metadata (caught late)"
                t["Similarity_Score"] = 0.0
                retained.append(t)
            else:
                t["Reason"] = "Excluded: not retained by filtering"
                t["Similarity_Score"] = 0.0
                excluded.append(t)

    return retained, excluded

FAILED_GENES = []  # global tracker for timeouts or failed requests
FAILED_LOCK = threading.Lock()  # lock for thread-safe updates

def process_one_gene(gene_id):
    """Process a single gene and return retained + excluded transcripts"""
    print(f"Processing gene {gene_id}...")
    try:
        gene_info = fetch_gene_info(gene_id)
        retained, excluded = process_gene(gene_id, appris_data, gencode_data)

        # Add gene-level info
        for t in retained + excluded:
            t["Gene_ID"] = gene_info["Gene_ID"]
            t["Gene_Name"] = gene_info["Gene_Name"]
            t["Gene_Biotype"] = gene_info["Gene_Biotype"]
            t.setdefault("Avg_Similarity", "NA")
            t.setdefault("Num_Comparisons", "NA")

        return retained, excluded

    except requests.exceptions.Timeout:
        print(f"Timeout fetching {gene_id}, skipping for now.")
        with FAILED_LOCK:
            FAILED_GENES.append(gene_id)
        return [], []

    except requests.exceptions.RequestException as e:
        print(f"Request failed for {gene_id}: {e}")
        with FAILED_LOCK:
            FAILED_GENES.append(gene_id)
        return [], []

    except Exception as e:
        print(f"Unexpected error for {gene_id}: {e}")
        with FAILED_LOCK:
            FAILED_GENES.append(gene_id)
        return [], []

# MAIN SCRIPT
if __name__ == "__main__":

    # Load metadata
    appris_data = load_appris(APPRIS_FILE)
    gencode_data = load_gencode(GENCODE_GTF)

    # Read gene IDs
    GENE_LIST_FILE = "HS_Gene_IDs115.txt"
    with open(GENE_LIST_FILE) as f:
        gene_list = [line.strip() for line in f if line.strip()]

    all_retained = []
    all_excluded = []

    headers = [
        "Gene_ID","Gene_Name","Gene_Biotype","Protein_Length","Transcript_ID",
        "Transcript_Biotype","Nucleotide_Length","APPRIS","APPRIS_Combined","MANE",
        "GENCODE","CCDS","UniProt_ID","TRIFID","TSL","Similarity_Score","Num_Comparisons","Reason"
    ]

    # Pass 1: main parallel run
    max_workers = 4
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        results = list(executor.map(process_one_gene, gene_list))

    # Aggregate results from first pass
    for retained, excluded in results:
        all_retained.extend(retained)
        all_excluded.extend(excluded)

    # Pass 2: retry failed genes 
    if FAILED_GENES:
        print(f"\nRetrying {len(FAILED_GENES)} failed genes...\n")
        time.sleep(5)  # small pause before retry

        retry_failed = FAILED_GENES.copy()
        FAILED_GENES.clear()

        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            retry_results = list(executor.map(process_one_gene, retry_failed))

        for retained, excluded in retry_results:
            all_retained.extend(retained)
            all_excluded.extend(excluded)

    # Save final combined outputs
    pd.DataFrame(all_retained).reindex(columns=headers).to_csv(OUTPUT_RETAINED, index=False)
    pd.DataFrame(all_excluded).reindex(columns=headers).to_csv(OUTPUT_EXCLUDED, index=False)

    print(f"Done. Total Retained: {len(all_retained)} | Total Excluded: {len(all_excluded)}")
    if FAILED_GENES:
        print(f"Still failed after retry: {len(FAILED_GENES)} genes")
