"""Independent bounded arithmetic review; never imports/runs production analysis."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, sys, math
import numpy as np
import pandas as pd

R = Path(__file__).resolve().parent
P = R.parents[3]
S = P/"next_steps/analysis_strengthening_15/single_cell"
I = P/"next_steps/publication_package_14/data/single_cell"
OLD = P/"next_steps/single_cell_extension_14/analysis"
checks=[]
def sha(f):
    h=hashlib.sha256()
    with Path(f).open("rb") as stream:
        for b in iter(lambda:stream.read(1048576),b""):h.update(b)
    return h.hexdigest()
def check(name,a,b,atol=1e-10):
    aa=np.asarray(a,dtype=float);bb=np.asarray(b,dtype=float)
    assert aa.shape==bb.shape,(name,aa.shape,bb.shape)
    assert np.array_equal(np.isnan(aa),np.isnan(bb)),name+" NA mismatch"
    ok=np.isfinite(aa)&np.isfinite(bb)
    error=float(np.max(np.abs(aa[ok]-bb[ok]))) if np.any(ok) else 0.
    assert np.allclose(aa,bb,rtol=0,atol=atol,equal_nan=True),(name,error)
    checks.append({"name":name,"elements":aa.size,"NA":int(np.isnan(aa).sum()),"max_abs_error":error,"absolute_tolerance":atol,"passed":True})
def boolean(name,condition):
    assert bool(condition),name
    checks.append({"name":name,"passed":True})
def compare_table(name,own,keys,cols,atol=1e-10,location=S):
    expected=pd.read_csv(location/name).set_index(keys).sort_index()
    actual=own.set_index(keys).sort_index()
    boolean(name+" exact key set",actual.index.equals(expected.index))
    for c in cols:check(name+" "+c,actual[c],expected[c],atol)
def save(name,value):
    (R/name).write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf8")

summary=json.loads((S/"analysis_summary.json").read_text(encoding="utf8"))
oldreview=json.loads((P/"next_steps/single_cell_extension_14/review/review.json").read_text(encoding="utf8"))
oldhash={x["file"]:x["sha256"] for x in oldreview["source_stability"]}
manifest=json.loads((S/"artifact_manifest.json").read_text(encoding="utf8"))
input_files=[Path(x["path"]) for x in summary["inputs"]]
protected=input_files+[S/x["file"] for x in manifest]+[S/"artifact_manifest.json"]
start_hash={str(f):sha(f) for f in protected}
source_checks=[]
for x in summary["inputs"]:
    f=Path(x["path"])
    original=OLD/x["file"]
    observed=sha(f)
    boolean("source matches current production provenance "+x["file"],observed==x["sha256"])
    boolean("source matches original stage14 review "+x["file"],observed==oldhash[x["file"]])
    boolean("stage14 original copy unchanged "+x["file"],sha(original)==oldhash[x["file"]])
    source_checks.append({"file":x["file"],"current_sha256":observed,"stage14_frozen_sha256":oldhash[x["file"]],"original_copy_sha256":sha(original),"passed":True})
for x in manifest:
    boolean("production artifact frozen "+x["file"],sha(S/x["file"])==x["sha256"])
boolean("frozen protocol",sha(S/"protocol.md")==summary["protocol_sha256"]=="4ab6ea1a388b1c5fc66e1e3c172acbb6baeb7d6391d630357a25995e68034edd")

data=pd.read_csv(I/"cell_level_umap_qc_expression.csv.gz")
coverage=pd.read_csv(I/"candidate_gene_coverage.csv")
genes=coverage.symbol.tolist()
available=coverage.covered.to_numpy(dtype=bool)
boolean("25 fixed/24 measured/AOC4P absent",len(genes)==25 and available.sum()==24 and coverage.loc[~available,"symbol"].tolist()==["AOC4P"])
boolean("library barcode key unique",not data.duplicated(["gsm","barcode"]).any() and data.cell_key.is_unique)
boolean("key built from library and barcode",data.cell_key.eq(data.gsm+"|"+data.barcode).all())
boolean("24034 cells",len(data)==24034)
liblist=sorted(data.gsm.unique())
boolean("fixed source libraries",liblist==["GSM6322577","GSM6322578","GSM6322579","GSM6322580"])
check("per-library counts",data.groupby("gsm").size(),[3093,5510,7399,8032],0)
boolean("control then PA source labels",[data.loc[data.gsm==l,"condition"].unique().tolist() for l in liblist]==[["CTRL"],["CTRL"],["PA"],["PA"]])
boolean("positive raw all-gene denominator",(data.raw_total_counts>0).all())
counts=data[[g+"__counts" for g in genes]].copy()
counts.columns=genes
y=np.log1p(counts.divide(data.raw_total_counts,axis=0)*10000)
stored=data[[g+"__log1p_cp10k" for g in genes]].to_numpy()
check("per-cell all25 normalized values",y.to_numpy(),stored,1e-10)
boolean("absent AOC4P source NA",counts.AOC4P.isna().all() and y.AOC4P.isna().all())
boolean("measured counts nonnegative integer",np.all(counts.loc[:,available]>=0) and np.all(np.mod(counts.loc[:,available],1)==0))

def groups(fields):
    rows=[]
    for key,idx in data.groupby(fields,sort=True).groups.items():
        key=key if isinstance(key,tuple) else (key,)
        labels=dict(zip(fields,key));xx=counts.loc[idx];yy=y.loc[idx]
        for j,g in enumerate(genes):
            row={**labels,"symbol":g,"n_cells":len(idx)}
            if available[j]:
                det=int((xx[g]>0).sum())
                row.update(n_detected_cells=det,detected_pct=100*det/len(idx),mean_log1p_cp10k=float(yy[g].mean()),sum_raw_counts=float(xx[g].sum()))
            else:
                row.update(n_detected_cells=np.nan,detected_pct=np.nan,mean_log1p_cp10k=np.nan,sum_raw_counts=np.nan)
            rows.append(row)
    return pd.DataFrame(rows)
metrics=["n_cells","n_detected_cells","detected_pct","mean_log1p_cp10k","sum_raw_counts"]
cat=groups(["gsm","celltype_category"])
sub=groups(["gsm","celltype_scTypedb"])
lib=groups(["gsm"])
compare_table("library_category_candidates_all25.csv",cat,["gsm","celltype_category","symbol"],metrics)
compare_table("library_subtype_candidates_all25.csv",sub,["gsm","celltype_scTypedb","symbol"],metrics)
compare_table("library_candidates_all25.csv",lib,["gsm","symbol"],metrics)
compare_table("candidate25_library_subtype_all10.csv",sub,["gsm","celltype_scTypedb","symbol"],metrics,atol=1e-9,location=I)

categories=["Hepatocytes","Cholangiocytes","Fibroblasts","Hepatic stellate cells"]
sampling=data.groupby(["gsm","celltype_category"]).size().unstack().reindex(index=liblist,columns=categories)
boolean("all four categories in all four libraries",(sampling>0).all().all())
fractions=sampling.divide(sampling.sum(axis=1),axis=0)
check("four fractions sum1",fractions.sum(axis=1),np.ones(4),1e-14)
catidx=cat.set_index(["gsm","celltype_category","symbol"])
locrows=[]
for l in liblist:
    for g in genes:
        ch=catidx.loc[(l,"Cholangiocytes",g)];he=catidx.loc[(l,"Hepatocytes",g)]
        locrows.append({"gsm":l,"symbol":g,"cholangiocyte_cells":ch.n_cells,"hepatocyte_cells":he.n_cells,
        "cholangiocyte_detection_pct":ch.detected_pct,"hepatocyte_detection_pct":he.detected_pct,
        "detection_difference_percentage_points":ch.detected_pct-he.detected_pct,
        "cholangiocyte_mean_log1p_cp10k":ch.mean_log1p_cp10k,"hepatocyte_mean_log1p_cp10k":he.mean_log1p_cp10k,
        "mean_expression_difference":ch.mean_log1p_cp10k-he.mean_log1p_cp10k})
localization=pd.DataFrame(locrows)
compare_table("cholangiocyte_vs_hepatocyte_per_library_all25.csv",localization,["gsm","symbol"],
["cholangiocyte_cells","hepatocyte_cells","cholangiocyte_detection_pct","hepatocyte_detection_pct","detection_difference_percentage_points","cholangiocyte_mean_log1p_cp10k","hepatocyte_mean_log1p_cp10k","mean_expression_difference"])
areg=localization.loc[localization.symbol=="AREG"].copy()
boolean("AREG detection positive in each technical library",(areg.detection_difference_percentage_points>0).all())
boolean("AREG mean expression positive in each technical library",(areg.mean_expression_difference>0).all())
areg.to_csv(R/"AREG_independent_per_library.csv",index=False,float_format="%.17g")

# Independently derive profiles from explicit per-cell weights, not production aggregates.
# EL weight per retained cell = 1/(two libraries * cells in that cell's library).
# Pooled weight per cell = 1/(all retained cells in that condition).
profiles=[];condition_values={};covariances={}
for condition in ["CTRL","PA"]:
    idx=data.index[data.condition==condition]
    subset=data.loc[idx]
    selected=sorted(subset.gsm.unique())
    boolean("two libraries "+condition,len(selected)==2)
    n=subset.groupby("gsm").size()
    p=fractions.loc[selected].to_numpy()
    mu=np.stack([catidx.loc[(l,c)].reindex(genes).mean_log1p_cp10k.to_numpy() for l in selected for c in categories]).reshape(2,4,25)
    # For exactly two equally weighted libraries, population covariance = dp*dmu/4.
    covariance=np.sum((p[0]-p[1])[:,None]*(mu[0]-mu[1])/4,axis=0)
    covariances[condition]=covariance
    direct_el=np.vstack([y.loc[data.index[data.gsm==l]].mean().to_numpy() for l in selected]).mean(axis=0)
    direct_pool=y.loc[idx].mean().to_numpy()
    for scheme in ["equal_library_actual","pooled_cell_actual","separate_average_reference_diagnostic"]:
        if scheme=="equal_library_actual":
            w=1/(2*subset.gsm.map(n).to_numpy())
        elif scheme=="pooled_cell_actual":
            w=np.full(len(idx),1/len(idx))
        else:
            w=None
        pp=[];mm=[]
        for k,c in enumerate(categories):
            belongs=subset.celltype_category.eq(c).to_numpy()
            if w is not None:
                mass=w[belongs].sum()
                yy=y.loc[idx[belongs]].to_numpy()
                mean=(w[belongs,None]*yy).sum(axis=0)/mass
            else:
                mass=float(p[:,k].mean())
                mean=mu[:,k,:].mean(axis=0)
            pp.append(mass);mm.append(mean)
            for j,g in enumerate(genes):
                profiles.append({"symbol":g,"condition":condition,"scheme":scheme,"celltype_category":c,
                   "category_fraction":mass,"category_mean_log1p_cp10k":mean[j],"category_contribution_to_profile_mean":mass*mean[j]})
        pp=np.array(pp);mm=np.array(mm)
        reconstructed=(pp[:,None]*mm).sum(axis=0)
        direct=direct_el if scheme=="equal_library_actual" else direct_pool if scheme=="pooled_cell_actual" else reconstructed
        check("profile direct condition mean "+condition+" "+scheme,reconstructed,direct,1e-12)
        condition_values[(condition,scheme)]=(pp,mm,direct)
    ref=condition_values[(condition,"separate_average_reference_diagnostic")][2]
    check("population covariance reconstructs EL mean "+condition,ref+covariance,direct_el,1e-12)
profiledf=pd.DataFrame(profiles)
compare_table("condition_category_profiles_all25.csv",profiledf,["symbol","condition","scheme","celltype_category"],
["category_fraction","category_mean_log1p_cp10k","category_contribution_to_profile_mean"],atol=1e-12)
terms=[];totals=[];covrows=[]
for scheme in ["equal_library_actual","pooled_cell_actual","separate_average_reference_diagnostic"]:
    p0,m0,d0=condition_values[("CTRL",scheme)]
    p1,m1,d1=condition_values[("PA",scheme)]
    c=(p1-p0)[:,None]*(m1+m0)/2
    w=(m1-m0)*(p1+p0)[:,None]/2
    delta=d1-d0
    check("category-level product identity "+scheme,c+w,p1[:,None]*m1-p0[:,None]*m0,1e-12)
    check("all24 total identity "+scheme,c.sum(axis=0)+w.sum(axis=0),delta,1e-12)
    reverse_c=(p0-p1)[:,None]*(m0+m1)/2
    reverse_w=(m0-m1)*(p0+p1)[:,None]/2
    check("swap reverses C "+scheme,reverse_c,-c,1e-14)
    check("swap reverses W "+scheme,reverse_w,-w,1e-14)
    for j,g in enumerate(genes):
        totals.append({"symbol":g,"scheme":scheme,"control_profile_mean":d0[j],"PA_profile_mean":d1[j],
          "total_PA_minus_control":delta[j],"composition_total":c[:,j].sum(),"within_category_total":w[:,j].sum()})
        for k,name in enumerate(categories):
            terms.append({"symbol":g,"scheme":scheme,"celltype_category":name,"control_fraction":p0[k],"PA_fraction":p1[k],
            "control_mean_log1p_cp10k":m0[k,j],"PA_mean_log1p_cp10k":m1[k,j],"composition_term":c[k,j],
            "within_category_term":w[k,j],"category_total_difference":c[k,j]+w[k,j]})
tdf=pd.DataFrame(terms);tot=pd.DataFrame(totals)
compare_table("kitagawa_terms_all25.csv",tdf,["symbol","scheme","celltype_category"],
["control_fraction","PA_fraction","control_mean_log1p_cp10k","PA_mean_log1p_cp10k","composition_term","within_category_term","category_total_difference"],atol=1e-12)
compare_table("kitagawa_summary_all25.csv",tot,["symbol","scheme"],
["control_profile_mean","PA_profile_mean","total_PA_minus_control","composition_total","within_category_total"],atol=1e-12)
pivot=tot.pivot(index="symbol",columns="scheme",values="total_PA_minus_control").reindex(genes)
resid=covariances["PA"]-covariances["CTRL"]
check("diagnostic plus covariance equals actual EL contrast",
pivot.separate_average_reference_diagnostic.to_numpy()+resid,pivot.equal_library_actual.to_numpy(),1e-12)
prod_total=pd.read_csv(S/"kitagawa_summary_all25.csv")
ref=prod_total.loc[prod_total.scheme=="separate_average_reference_diagnostic"].set_index("symbol").reindex(genes)
check("reported diagnostic covariance delta",resid,ref.diagnostic_covariance_delta,1e-12)
covdf=pd.DataFrame({"symbol":genes,"CTRL_population_covariance":covariances["CTRL"],"PA_population_covariance":covariances["PA"],
"PA_minus_CTRL_covariance":resid,"actual_EL_delta":pivot.equal_library_actual.to_numpy(),"reference_delta":pivot.separate_average_reference_diagnostic.to_numpy()})
covdf.to_csv(R/"covariance_independent_check.csv",index=False,na_rep="NA",float_format="%.17g")
agree=np.sign(pivot.equal_library_actual.to_numpy())==np.sign(pivot.pooled_cell_actual.to_numpy())
disagree=[g for g,a,measured in zip(genes,agree,available) if measured and not a]
boolean("21 of24 weighting signs agree",int(agree[available].sum())==summary["weighted_total_sign_agreement"]==21)
boolean("three discordant candidates retained",disagree==summary["weighted_total_sign_disagreements"]==["PSTPIP2","BCAS1","C2orf49-DT"])
boolean("AOC4P all decomposition metrics remain NA",tot.loc[tot.symbol=="AOC4P",["control_profile_mean","PA_profile_mean","total_PA_minus_control","composition_total","within_category_total"]].isna().all().all())
zero_count=int((cat.loc[cat.symbol!="AOC4P","n_detected_cells"]==0).sum())
boolean("measured zero detections remain real zeros",zero_count>0)
boolean("no P-value outputs",not any(c.lower() in ["p","pvalue","p_value","q_value","fdr"] for fn in S.glob("*.csv") for c in pd.read_csv(fn,nrows=0).columns))
boolean("all protected inputs unchanged",all(sha(f)==h for f,h in start_hash.items()))

# Quantitative scope counters count compared array entries, not independent tests.
numeric=[c for c in checks if "elements" in c]
review={
"status":"PASS_BOUNDED_TECHNICAL_REVIEW_WITH_CAVEATS","confidence":"Share with caveats",
"reviewer":"Codex internal sub-agent /root/figures_support12","executed_utc":datetime.now(timezone.utc).isoformat(),
"scope":"Independent arithmetic from stage14 per-cell table; source stability, full25 NA completeness, all-library AREG localization, 24 measured candidates x4 categories under both actual-mean schemes and the diagnostic reference scheme.",
"not_performed":["New raw MTX/h5ad processing","Statistical inference from cells/libraries","New biological or causal validation","New literature search","Final manuscript/figure visual review","Claude or human expert review"],
"environment":{"executable":sys.executable,"python":sys.version,"numpy":np.__version__,"pandas":pd.__version__},
"production_script_sha256":sha(S/"analyze_single_cell.py"),"protocol_sha256":sha(S/"protocol.md"),
"review_script_sha256":sha(__file__),"production_manifest_sha256":sha(S/"artifact_manifest.json"),
"stage14_source_stability":source_checks,"protected_source_and_artifact_sha256":start_hash,
"all_checks_passed":len(checks),"numeric_comparison_batches":len(numeric),"numeric_elements_compared":sum(c["elements"] for c in numeric),
"checks":checks,"AREG_per_library":areg.to_dict("records"),
"weighting_sign_agreement":{"numerator":int(agree[available].sum()),"denominator":24,"discordant":disagree},
"completeness":{"cells":len(data),"unique_library_barcode_keys":data.cell_key.nunique(),"barcode_labels_repeated_across_libraries":int((data.groupby("barcode").gsm.nunique()>1).sum()),"candidate_count":25,"measured_candidates":24,"missing":"AOC4P","category_zero_detection_rows_retained":zero_count,"HB1_cells_per_library":data.loc[data.celltype_scTypedb=="HB1"].groupby("gsm").size().to_dict()},
"mathematical_review":{"EL_mu":"Direct cell weights 1/(2*N_library) produce category mass p and weighted category mean mu; sum_k p*mu reconstructs the mean of the two library-wide means.","pooled_mu":"Equal weights per cell within condition reconstruct the direct pooled mean.","covariance":"For two equally weighted libraries, sum_k (p1-p2)*(mu1-mu2)/4 equals actual EL mean minus product of separate average profiles; PA-minus-CTRL residual restores the actual contrast.","algebra":"Per category, C+W=p_PA*mu_PA-p_CTRL*mu_CTRL. This is an exact algebraic identity, with floating-point verification and no causal assumption."},
"blocking_findings":[],"required_caveats":["Four technical libraries from one differentiation, not four donors or independent biological replicates.","Control/PA only; no sesamin.","Detection differences can reflect sequencing depth, recovery and technical composition.","Within-category contribution can contain finer subtype/state/measurement changes; it is not cell-autonomous or causal.","The separate-average diagnostic total is a reference mixture, not an observed equal-library mean.","Missing AOC4P remains NA and measured zeros remain zero.","PSTPIP2, BCAS1 and C2orf49-DT change total direction under different weighting and must not be hidden."],
"source_unchanged_after_review":True
}
save("review.json",review)
lines=[
"# Bounded independent review: stage15 single-cell analysis","",
"**Assessment: Share with caveats. No numerical or source-integrity blocker found.**",
"",
"This review independently recomputed the requested quantities from the frozen stage14 per-cell table. It did not run or import the production analysis and did not reopen raw sequencing matrices.",
"",
"## Verified findings","",
f"- All {len(checks)} check batches passed; arithmetic counts are not biological replication.",
"- All five source inputs match both the original stage14 review hashes and the read-only original analysis copies. Production artifacts and sources remained unchanged.",
"- All 24,034 cells, 25 candidates, four libraries and four broad categories were retained. AOC4P remains unavailable; measured zero counts remain zeros.",
"- All 1,000 library/subtype rows and all 400 library/category rows agree with independent group aggregation. Normalization uses raw all-gene UMI denominators.",
"",
"| Library | Cholangiocyte AREG detection (%) | Hepatocyte AREG detection (%) | Difference (percentage points) |",
"|---|---:|---:|---:|"]
for row in areg.itertuples():
    lines.append(f"| {row.gsm} | {row.cholangiocyte_detection_pct:.6f} | {row.hepatocyte_detection_pct:.6f} | {row.detection_difference_percentage_points:+.6f} |")
lines += [
"",
"AREG has greater detection and mean expression in the author Cholangiocytes category than the Hepatocytes category in each library. These are descriptive within-library comparisons, not four independent replications.",
"",
"## Mixture estimands and identities","",
"The independent implementation uses explicit cell weights. Giving each cell weight 1/(2*N_library) makes every library contribute one half of its condition, and the resulting category-weighted means exactly reconstruct the average of the two library-wide means. Equal per-cell weights reconstruct the pooled-cell mean. Therefore, the primary mu weights are appropriate for the stated actual equal-library mixture.",
"",
"All measured candidates and categories satisfy the symmetric composition plus within-category identity under both actual-mean schemes. The diagnostic separate-average profile was checked separately; its finite-library covariance residual restores the actual equal-library contrast. Swapping conditions negates both terms.",
"",
"Total direction agrees under equal-library and pooled-cell weights for 21/24 candidates. PSTPIP2, BCAS1 and C2orf49-DT remain explicit discordances.",
"",
"## Required interpretation limits","",
"- Four technical libraries arise from one differentiation, with no sesamin arm. No cell/library P values or biological independence claims are supported.",
"- The within-category term includes possible finer-subtype, cell-state, recovery and depth effects. It is not a cell-autonomous or causal contribution.",
"- The diagnostic reference-mixture total is not the observed equal-library total; the covariance residual is necessary.",
"- This is an internal Codex technical review of code, arithmetic and sources, not Claude/human expert review, new experimental validation or a publication guarantee.",
"",
"## Reproducible review artifacts","",
"- check_single_cell.py: independent implementation.",
"- review.json: checks, source/artifact hashes, scope and caveats.",
"- AREG_independent_per_library.csv: actual four-library detection and mean comparisons.",
"- covariance_independent_check.csv: all25 diagnostic reconstruction inputs.",
"",
"Manuscript-specific wording and figures were not included in this numerical pass; those may be checked once their final added text is provided.",
]
(R/"review.md").write_text("\n".join(lines)+"\n",encoding="utf8")
print(json.dumps({"status":review["status"],"checks":len(checks),"numeric_batches":len(numeric),"AREG":review["AREG_per_library"],"weighting":review["weighting_sign_agreement"],"HB1":review["completeness"]["HB1_cells_per_library"],"review_sha256":sha(R/"review.json")},ensure_ascii=True))

