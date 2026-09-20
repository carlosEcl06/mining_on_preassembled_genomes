#!/usr/bin/env python3
# Pulls NCBI BioSample metadata for the near-duplicate pairs (near_duplicates.tsv) and writes near_duplicates_biosample.tsv.
#!/usr/bin/env python3
import csv
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
import json

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

def fetch(url):
    with urllib.request.urlopen(url, timeout=60) as r:
        return r.read()

pairs = list(csv.DictReader(
    open("/home/clima/lab/projects/bgcs_in_cancer/mt/mining_on_preassembled_genomes/analysis/phylogenetics/results/grafgen/near_duplicates.tsv"),
    delimiter="\t"
))
accessions = sorted(set(p["sample_a"] for p in pairs) | set(p["sample_b"] for p in pairs))
print(f"{len(accessions)} unique accessions across {len(pairs)} pairs")

# Batch esearch to resolve accessions -> assembly UIDs
term = " OR ".join(f"{a}[Assembly Accession]" for a in accessions)
url = f"{BASE}/esearch.fcgi?db=assembly&term={urllib.parse.quote(term)}&retmax=200&retmode=json"
uids = json.loads(fetch(url))["esearchresult"]["idlist"]
print(f"{len(uids)} assembly UIDs resolved")
time.sleep(0.4)

# Batch esummary to get biosample accessions per assembly UID
url = f"{BASE}/esummary.fcgi?db=assembly&id={','.join(uids)}&retmode=json"
summ = json.loads(fetch(url))["result"]
acc_to_biosample = {}
biosample_ids = []
for uid in uids:
    doc = summ[uid]
    acc_to_biosample[doc["assemblyaccession"]] = doc["biosampleid"]
    biosample_ids.append(doc["biosampleid"])
time.sleep(0.4)

# Batch efetch full BioSample XML records
url = f"{BASE}/efetch.fcgi?db=biosample&id={','.join(biosample_ids)}&retmode=xml"
xml_data = fetch(url)
root = ET.fromstring(xml_data)

biosample_info = {}
for bs in root.findall("BioSample"):
    bs_id = bs.get("id")
    attrs = {a.get("attribute_name"): (a.text or "").strip() for a in bs.findall(".//Attribute")}
    strain = ""
    for idnode in bs.findall(".//Id"):
        if idnode.get("db_label") == "Sample name":
            strain = idnode.text or ""
    bioproject = ""
    link = bs.find(".//Link[@target='bioproject']")
    if link is not None:
        bioproject = link.get("label", "")
    biosample_info[bs_id] = {
        "sample_name": strain,
        "bioproject": bioproject,
        "collection_date": attrs.get("collection_date", ""),
        "geo_loc_name": attrs.get("geo_loc_name", ""),
        "host_age": attrs.get("host_age", ""),
        "host_sex": attrs.get("host_sex", ""),
        "host_disease": attrs.get("host_disease", ""),
        "host_tissue_sampled": attrs.get("host_tissue_sampled", ""),
        "isolation_source": attrs.get("isolation_source", ""),
    }

acc_info = {acc: biosample_info.get(bsid, {}) for acc, bsid in acc_to_biosample.items()}

out_rows = []
for p in pairs:
    a, b, d = p["sample_a"], p["sample_b"], p["snp_distance"]
    ia, ib = acc_info.get(a, {}), acc_info.get(b, {})
    same_bioproject = bool(ia.get("bioproject")) and ia.get("bioproject") == ib.get("bioproject")
    same_date = bool(ia.get("collection_date")) and ia.get("collection_date") == ib.get("collection_date")
    same_age_sex = (ia.get("host_age"), ia.get("host_sex")) == (ib.get("host_age"), ib.get("host_sex")) and ia.get("host_age")

    import re
    def norm(name):
        return re.sub(r"^[A-Za-z]{1,4}-", "", name)  # strip a leading country/lab code like "MX-"
    name_a, name_b = ia.get("sample_name", ""), ib.get("sample_name", "")
    norm_a, norm_b = norm(name_a), norm(name_b)
    identical_sample_id = bool(norm_a) and norm_a == norm_b  # same isolate, redeposited under a different accession

    strain_a_root = name_a.rsplit("-", 1)[0] if "-" in name_a else ""
    strain_b_root = name_b.rsplit("-", 1)[0] if "-" in name_b else ""
    same_strain_root = bool(strain_a_root) and strain_a_root == strain_b_root

    verdict_signals = sum([same_bioproject, same_date, bool(same_age_sex), same_strain_root])
    if identical_sample_id:
        verdict = "likely duplicate deposit (identical sample ID)"
    elif verdict_signals >= 3:
        verdict = "likely same patient"
    elif verdict_signals >= 1:
        verdict = "possible same patient"
    else:
        verdict = "likely independent"
    out_rows.append({
        "sample_a": a, "sample_b": b, "snp_distance": d, "verdict": verdict,
        "sample_name_a": ia.get("sample_name", ""), "sample_name_b": ib.get("sample_name", ""),
        "tissue_a": ia.get("host_tissue_sampled", ""), "tissue_b": ib.get("host_tissue_sampled", ""),
        "age_a": ia.get("host_age", ""), "sex_a": ia.get("host_sex", ""),
        "age_b": ia.get("host_age", ""), "sex_b": ib.get("host_sex", ""),
        "bioproject_a": ia.get("bioproject", ""), "bioproject_b": ib.get("bioproject", ""),
        "collection_date_a": ia.get("collection_date", ""), "collection_date_b": ib.get("collection_date", ""),
        "isolation_source_a": ia.get("isolation_source", ""), "isolation_source_b": ib.get("isolation_source", ""),
    })

out_path = "/home/clima/lab/projects/bgcs_in_cancer/mt/mining_on_preassembled_genomes/analysis/phylogenetics/results/grafgen/near_duplicates_biosample.tsv"
with open(out_path, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(out_rows)

verdict_counts = {}
for r in out_rows:
    verdict_counts[r["verdict"]] = verdict_counts.get(r["verdict"], 0) + 1
print("Verdict counts:", verdict_counts)
print(f"Wrote {out_path}")
