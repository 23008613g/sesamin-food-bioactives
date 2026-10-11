"""Audit metadata and map the existing count header; no expression analysis."""
from pathlib import Path
import csv
import datetime
import gzip
import hashlib
import json
import re
import sys
import xml.etree.ElementTree as ET

OUT = Path(__file__).resolve().parent
PROJECT = OUT.parents[2]
COUNT = PROJECT / "downloads/GSE247407_gene_count.txt.gz"
SOFT = OUT / "GSE247407_family.soft.gz"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    with gzip.open(SOFT, "rt", encoding="utf-8") as handle:
        soft = handle.read()
    samples = {}
    current = None
    for line in soft.splitlines():
        if line.startswith("^SAMPLE = "):
            current = line.split(" = ", 1)[1]
            samples[current] = {}
        elif line.startswith("^"):
            current = None
        elif current and line.startswith("!") and " = " in line:
            key, value = line.split(" = ", 1)
            samples[current].setdefault(key, []).append(value)

    # Deliberately read only the count-file header, not gene-expression rows.
    with gzip.open(COUNT, "rt", encoding="utf-8") as handle:
        header_line = handle.readline()
    header = header_line.rstrip("\r\n").split("\t")
    expr_columns = header[1:31]
    assert len(header) == 40 and header[31] == "gene_name"
    assert len(expr_columns) == len(set(expr_columns)) == 30
    by_column = {}
    for gsm, fields in samples.items():
        chars = dict(x.split(": ", 1) for x in fields["!Sample_characteristics_ch1"])
        columns = [x.split(": ", 1)[1] for x in fields["!Sample_description"]
                   if x.startswith("Column in matrix: ")]
        assert len(columns) == 1
        column = columns[0]
        assert column not in by_column
        rep = re.search(r"rep([123])$", fields["!Sample_title"][0]).group(1)
        donor = chars["individual"]
        condition = chars["treatment"]
        expected_column = donor + ("C" if condition == "Control" else "F") + rep
        assert column == expected_column
        relations = fields.get("!Sample_relation", [])
        by_column[column] = {
            "matrix_column": column,
            "gsm": gsm,
            "title": fields["!Sample_title"][0],
            "donor_id": donor,
            "condition": condition,
            "replicate_label": rep,
            "donor_condition_id": f"{donor}_{condition}",
            "population_inference_unit": "donor_reported_identity_not_independently_verified",
            "donor_identity_status": "author_reports_five_donors_label_relationship_not_verified",
            "replicate_nature": "unresolved_not_additional_donors",
            "exposure_days": 7,
            "fatty_acid_composition": "not_verified" if condition == "FFA" else "no_added_FFA_reported",
            "fatty_acid_concentration": "not_verified" if condition == "FFA" else "not_applicable",
            "vehicle_composition": "not_verified",
            "donor_sex": "not_verified",
            "donor_age": "not_verified",
            "donor_ethnicity": "not_verified",
            "culture_batch": "not_reported_in_GEO",
            "sequencing_batch": "not_reported_in_GEO",
            "biosample": next((x.rsplit("/", 1)[-1] for x in relations if x.startswith("BioSample:")), ""),
            "sra_experiment": next((x.rsplit("=", 1)[-1] for x in relations if x.startswith("SRA:")), ""),
            "metadata_source": SOFT.name,
        }
    assert set(by_column) == set(expr_columns)
    manifest = []
    for idx, col in enumerate(expr_columns, 1):
        manifest.append({"count_column_1_based": idx, **by_column[col]})
    donor_groups = {}
    for row in manifest:
        donor_groups.setdefault(row["donor_id"], {}).setdefault(row["condition"], []).append(row["matrix_column"])
    assert len(donor_groups) == 5
    assert all(set(x) == {"Control", "FFA"} and all(len(y) == 3 for y in x.values())
               for x in donor_groups.values())
    with (OUT / "sample_manifest.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest[0]))
        writer.writeheader()
        writer.writerows(manifest)

    pubmed = ET.parse(OUT / "pubmed_records.xml")
    publications = []
    for article in pubmed.findall(".//PubmedArticle"):
        ids = {x.get("IdType"): x.text for x in article.findall("./PubmedData/ArticleIdList/ArticleId")}
        publications.append({
            "pmid": article.findtext("./MedlineCitation/PMID"),
            "title": "".join(article.find("./MedlineCitation/Article/ArticleTitle").itertext()),
            "doi": ids.get("doi"), "pmcid": ids.get("pmc"),
            "publication_types": [x.text for x in article.findall(".//PublicationTypeList/PublicationType")],
            "correction_links": [dict(x.attrib) for x in article.findall("./MedlineCitation/CommentsCorrectionsList/CommentsCorrections")],
        })
    xml_root = ET.parse(OUT / "PMC13453441.xml").getroot()
    body_text = " ".join(xml_root.find("body").itertext())
    terms = ["CXCL8", "IL-8", "interleukin-8", "VNN1", "vanin-1"]
    keyword_presence = {term: bool(re.search(re.escape(term), body_text, flags=re.I)) for term in terms}
    metadata = {
        "created_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "dataset": "GSE247407",
        "scope": "Metadata/design audit and true count-column mapping only; no expression results generated.",
        "skill": "nature-academic-search; multi-source-search with official stdlib API fallback",
        "url_requests": {"web_attempts": 3, "successful_small_source_downloads": 5, "total_attempted_urls": 8,
                         "additional_authorized_focused_reads": 2},
        "source_files": [
            {"file": name, "sha256": digest(OUT / name)}
            for name in ["GSE247407_family.soft.gz", "pubmed_records.xml", "PMC13453441.xml",
                         "europepmc_original.json", "biosample_F95_F95P.xml"]
        ],
        "local_count_file_read_scope": "header only",
        "count_header_sha256_utf8": hashlib.sha256(header_line.encode("utf-8")).hexdigest(),
        "matrix_total_columns": len(header),
        "expression_columns": expr_columns,
        "annotation_columns": [header[0], *header[31:]],
        "manifest_rows": len(manifest),
        "unique_GSM": len({x["gsm"] for x in manifest}),
        "unique_biosamples": len({x["biosample"] for x in manifest}),
        "unique_SRA_experiments": len({x["sra_experiment"] for x in manifest}),
        "donor_condition_groups": donor_groups,
        "mapping_checks": {
            "exact_GEO_description_to_count_header": True,
            "title_donor_condition_rep_suffix_agree": True,
            "same_donor_control_FFA_pairing": True,
            "cross_condition_rep1_to_rep1_pairing_proven": False,
            "all_submitted_libraries_present": True,
        },
        "donor_independence": {
            "source_claim": "GEO: 5 patients; primary paper abstract: five donors varying in sex and ethnicity.",
            "verified_scope": "Five distinct individual labels with balanced within-label treatment; five donors is an author report, not independently established identity.",
            "not_verified": "Genetic identity/kinship, cell-lot identity documents, and donor clinical table.",
            "F95_vs_F95P": "GEO and separately fetched BioSample records label F95 and F95P, without explaining whether they are unrelated donors versus related cell lots/subbatches. Separate genetic identity was not tested; neither automatic merging nor confirmed independent n is justified.",
            "donor_labels_n": 5,
            "author_reported_donors_n": 5,
            "independently_confirmed_donors_n": None,
            "libraries_n": 30,
        },
        "protocol": {
            "model": "Cryopreserved primary human hepatocytes; 3D collagen sandwich",
            "exposure_days": 7,
            "comparison": "FFA versus no added FFA control within each donor",
            "FA_species": None,
            "FA_dose": None,
            "FA_molar_ratio": None,
            "vehicle": None,
            "medium_changes": None,
            "replicate_nature": "GEO says 3 reps only; technical split versus separate wells/culture batches not established.",
            "sex_age_ethnicity_per_donor": "Not present in retrieved GEO metadata; variation in sex/ethnicity reported only in abstract.",
            "batch_and_randomization": "Not documented in retrieved metadata.",
            "sesamin_intervention": False,
        },
        "processing": {
            "assembly": "hg38",
            "alignment": "Hisat2 v2.0.5",
            "counting": "FeatureCounts v1.5.0-p3",
            "instrument": "Illumina NovaSeq 6000",
            "library": "poly-T selected mRNA; paired-end alignment described",
            "exact_gene_annotation_release": "not_reported_in_retrieved_metadata",
        },
        "sample_exclusion": {
            "GEO_to_existing_count_header_missing": [],
            "original_paper_exclusion_rules": "unresolved; original full text/supplement not retrieved",
            "this_audit_excluded_samples": [],
        },
        "publications": publications,
        "novelty_checks": {
            "2024_original_full_text_status": "not retrieved; PubMed record has no PMCID; EuropePMC core record lists subscription required and inEPMC=N",
            "2024_CXCL8_VNN1_prior_reporting": "unresolved, cannot infer absence from abstract",
            "2024_already_reported_general_effects": "Steatosis, insulin resistance, mitochondrial dysfunction, inflammation and MASLD-related gene-signature changes.",
            "2026_role": "Secondary multispecies reanalysis of the same GSE247407 human data; not an independent human cohort or correction.",
            "2026_body_keyword_presence": keyword_presence,
            "2026_keyword_limit": "Main XML body search only; not every figure/supplement. Negative keyword search is not novelty evidence.",
        },
        "recommendations": [
            "Estimate FFA-Control within each donor label. There are five labels and author-reported donors, not fifteen or thirty independently confirmed people; label identity relationships remain unresolved.",
            "Do not pair replicate suffixes across conditions without culture-block evidence.",
            "For an exploratory display, normalize each library using all eligible genes jointly, summarize the three libraries per donor-condition, then report five paired donor effects and all library points.",
            "Do not normalize only the selected candidate genes or compare raw counts across libraries.",
            "If later verified to be technical sequencing replicates of the same RNA preparation, sum raw counts before normalization; if separate culture replicates, preserve nested structure or summarize by donor-condition and label the loss of within-donor variance.",
            "Avoid treating donor-condition sums as ten unrelated samples: donor blocking is still necessary. Unequal library depth requires normalization even after aggregation.",
            "Donor blocking accounts for fixed donor differences, but treatment-confounded culture or sequencing batches remain unresolved.",
            "Because F95/F95P identity relation is not documented, avoid asserting independently confirmed n=5; keep each label separate descriptively and show leave-one-label-out or omit-both-F95-label sensitivity if needed, clearly without deciding identity from expression similarity.",
            "Use direction, donor heterogeneity, and sensitivity analyses as descriptive background; the data cannot independently validate sesamin reversal or establish a new mechanism.",
        ],
        "analysis_not_performed": ["differential expression", "normalization", "gene effect calculation", "statistical testing", "raw-read donor identity test", "AI model scoring"],
        "verification": {
            "manifest_sha256": digest(OUT / "sample_manifest.csv"),
            "script_sha256": digest(Path(__file__)),
            "python": sys.version,
            "interpreter": sys.executable,
            "assertions_passed": True,
        },
    }
    (OUT / "metadata_finding.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps({"manifest_rows": len(manifest), "donors": len(donor_groups),
                      "mapped_columns": len(by_column), "keyword_presence_2026": keyword_presence,
                      "assertions_passed": True}))


if __name__ == "__main__":
    main()
