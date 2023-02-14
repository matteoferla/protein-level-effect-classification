## Analysis

Submitting to Venus, in bulk.

Venus caches results for a few hours, so one can do a test:

```python
from typing import Dict, Tuple

class Analyze:
    """
    Expects columns uniprot, uniprot_mutation
    """
    results = {}
    
    def __init__(self):
        self.errors: Dict[Tuple[str, str], Exception] = {}
        
    def __call__(self, row):
        if isinstance(row.uniprot, float):
            return {}
        if isinstance(row.uniprot_mutation, float):
            return {}
        if 'missense' not in row['Consequence (most severe)']:
            return {}
        if (row.uniprot, row.uniprot_mutation) in results:
            return results[row.uniprot, row.uniprot_mutation]
        try:
            analysis = venus.analyse(gene=row.uniprot, mutation=row.uniprot_mutation)
            self.results[row.uniprot, row.uniprot_mutation] = analysis
            return analysis
        except Exception as error:
            self.errors[row.uniprot, row.uniprot_mutation] = error
            self.results[row.uniprot, row.uniprot_mutation] = {}
            return {}
        
analyze = Analyze()
# subset test:
ref.sample(frac=1).apply(analyze, axis=1)
```

What went wrong?

```python
import operator
errors = pd.concat([
                    pd.DataFrame(analyze.errors.keys(), columns=['uniprot', 'mutation']),
                    pd.Series(analyze.errors.values(), name='error'),
                    ],axis=1)
errors['error_type'] = errors.error.apply(type).apply(operator.attrgetter('__name__'))
errors
```

The above may need tweaking if getting reused due the stored data
```python
analyze.results = {k: v for k, v in analyze.results.items() if len(v)}
```
Store? This API step can take a while.
```python
import pickle

with open('results.pkl', 'wb') as fh:
    pickle.dump(results, fh)
```

Note that data in results is not straightforward (sorry).

```python
# data is results looks like:
analysis = results[('Q06330', 'p.Leu179Val')]
analysis['protein']['gene_name'], \
analysis['protein']['uniprot'], \
analysis['mutation']['clean_mutation'], \
analysis['ddG']['ddG']
```

As a result some polishing is needed.
Say `ref.apply(analyze, axis=1)` was run in full, then:

```python
def summarize(data) -> dict:
    nan = float('nan')
    if isinstance(data, float):
        return {}
    if len(data) == 0:
        return {}
    raw_distance = data['structural']['distance_to_closest_ligand']
    distance_to_closest_ligand = raw_distance if raw_distance and raw_distance != 'None' and raw_distance is not None else nan
    neighbours=data['structural']['neighbours']
    neighbor_ptm_details=[n for n in neighbours if len(n['ptms'])]
    neighbor_other_chain_details=[n for n in neighbours if n['other_chain'] and len(n['resn']) == 1]
    neighbor_gnomAD_details=[n for n in neighbours for g in n['gnomads']]
    neighbor_clinvar_details=[n for n in neighbours for g in n['clinvars']]
    general = dict(gene_name=data['protein']['gene_name'],
                uniprot=data['protein']['uniprot'],
                mutation=data['mutation']['clean_mutation'],
                ddG=data['ddG']['ddG'],
                dsol=data['ddG']['dsol'],
                rsa=data['structural']['RSA'],
                distance_to_closest_ligand=distance_to_closest_ligand,
                closest_ligand=data['structural']['closest_ligand'] if distance_to_closest_ligand and data['structural']['closest_ligand'] != 'None' else '',
                neighbor_ptm_details=neighbor_ptm_details,
                neighbor_other_chain_details=neighbor_other_chain_details, # filter for peptide
                neighbor_gnomAD_details=neighbor_gnomAD_details,
                neighbor_clinvar_details=neighbor_clinvar_details,
#                 N_neighbor_other_chain=sum([1 for n in protein.structural.neighbours if n['other_chain'] and len(n['resn']) == 1]),
#                 N_neighbor_ptms=sum([1 for n in protein.structural.neighbours if len(n['ptms'])]),
               )
    return {**general}

keys: pdt.Series[Tuple[str, str]] = pd.Series(zip(ref.uniprot, ref.uniprot_mutation))
data_series = keys.apply(lambda k: results[k] if k in results else float('nan'))
extra = pd.DataFrame( data_series.apply(summarize).to_list() )

ref = pd.concat([
                ref.drop(columns=set(extra.columns).intersection(ref.columns)), 
                extra.drop(columns=[c for c in extra.columns if '_details' in c])
                ], 
                axis=1)
```

## Data from Venus's data
An alternative for gene to uniprot is:

```python
import json, os

from michelanglo_protein import ProteinAnalyser, Structure, Mutation, global_settings
global_settings.startup('👾👾👾/protein-data')
human = json.load(open(os.path.join(global_settings.dictionary_folder, 'taxid9606-names2uniprot.json')))
ref['uniprot_fom_gene'] = ref.gene.map(human)
```
which is what Venus will use anyway. But can give extra details, for example:

```python
import operator

protein_series: pdt.Series[ProteinAnalyser] = ref.uniprot_fom_gene.apply(lambda uniprot: ProteinAnalyser(uniprot=str(uniprot), taxid=9606))
protein_series.apply(lambda p: p.load() if p.exists() else None)
ref['protein_length'] = protein_series.apply(len)
ref['pLI']: pdt.Series[float] = protein_series.apply(lambda p: p.pLI).astype(float)
ref['pRec']: pdt.Series[float] = protein_series.apply(lambda p: p.pRec).astype(float)
ref['pNull']: pdt.Series[float] = protein_series.apply(lambda p: p.pNull).astype(float)

def in_clinvar(a):
    (resi, p) = a
    if not resi:
        return ''
    return ','.join([c.mutation for c in p.clinvar if c.x == resi])
    
def in_gnomAD(a):
    (resi, p) = a
    if not resi:
        return ''
    return ','.join([c.mutation for c in p.gnomAD if c.x == resi])

def N_del(p):
    #set(map(operator.attrgetter('consequence'), p.gnomAD))
    return sum([g.N for g in p.gnomAD if g.consequence != 'missense_variant' and g.N])

is_valid = lambda m: isinstance(m, str) and m.find('+') == -1 and m.find('-') == -1
resi: pdt.Series[int] = ref.uniprot_mutation.apply(lambda m: Mutation(m).residue_index if is_valid(m) else 0)
ref['in_clinvar']: pdt.Series[str] = pd.Series(zip(resi, protein_series)).apply(in_clinvar)
ref['in_gnomAD']: pdt.Series[str] = pd.Series(zip(resi, protein_series)).apply(in_gnomAD)
ref['N_del']: pdt.Series[int] = protein_series.apply(N_del).astype(int)
```

Regarding dominance as predicted from the dataset, one can get the mismatches
```python
pLI_dominance_mismatch: pdt.Series[bool] = ref.dominant & (ref.pLI < 0.7)
```

