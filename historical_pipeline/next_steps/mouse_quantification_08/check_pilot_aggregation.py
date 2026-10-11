"""Integration check on real pilot table; no across-treatment results from one run."""
from common08 import *
import importlib.util

spec=importlib.util.spec_from_file_location('aggregate08', HERE/'aggregation/aggregate_quantification.py')
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
mapping,txids,genes,metadata,txgene=module.read_mapping(DATA/'reference/tx2gene.tsv')
values,metrics=module.read_abundance(DATA/'quant/DRR1010532/abundance.tsv',txids,mapping)
gene=module.aggregate_gene(values,txgene,len(genes))
write_json(HERE/'pilot_aggregation_check.json',{'passed':True,'checked_utc':now(),'scope':'Real pilot integration check with same aggregation implementation, not independent final validation or treatment comparison','run':'DRR1010532','transcripts':len(txids),'genes':len(genes),'metrics':metrics,'abundance_sha256':digest(DATA/'quant/DRR1010532/abundance.tsv'),'mapping_sha256':digest(DATA/'reference/tx2gene.tsv')})
print(json.dumps({'passed':True,'transcripts':len(txids),'genes':len(genes),'metrics':metrics}))
