## Clean-up

A table may not contain everything in the right format, say we have:
```python
import pandas as pd

ref = pd.read_csv('👾👾👾.csv')
```

Which contains a column `HGVS` (Human Genome Variation Society, a nomenclature standard),
but we want a Uniprot and a Uniprot variant, [Venus](https://venus.sgc.ox.ac.uk/) can be queried to fix it.
The following adds the columns 'uniprot', 'uniprot_mutation', 'link' to `ref`

```python
import operator
from michelanglo_api import VenusAPI
import pandera.typing as pdt
from typing import Dict, Any

venus = VenusAPI()

def target_split(target, empty='') -> Dict[str, str]:
    """
    keys: uniprot, uniprot_mutation, link
    """
    blank = {'uniprot': empty, 'uniprot_mutation': empty, 'link': empty}
    if isinstance(target, float):  # nan
        return blank
    if 'ENST' not in target:
        return blank
    try:
        target = target.split(',')[-1]
        #print(target)
        enst, mutation = target.split(':')
        temp = venus.post_json('venus_transcript',data=dict(enst=enst, mutation=mutation) )
        uniprot = temp['uniprot']
        uniprot_mutation = temp['mutation']
        link = f"https://venus.sgc.ox.ac.uk/venus?uniprot={uniprot}"+\
              f"&species=9606&mutation={uniprot_mutation}"
        return dict(uniprot = uniprot, 
                     uniprot_mutation=uniprot_mutation,
                     link=link)
    except Exception as error:
        print(f'{error.__class__.__name__}: {error}')
        return blank

extras: pdt.Series[Dict[str, str] = ref['HGVS'].apply(target_split)
for k in ('uniprot', 'uniprot_mutation', 'link'):
    ref[k] = extra.apply(operator.itemgetter(k).astype(str)
```

Without fail, some entries will be wrong. It's easier to fix manually:

```python
i = 20
ref.at[i, 'uniprot'] = 'Q9NZC7'
ref.at[i, 'uniprot_mutation'] = 'p.His173Met'
ref.at[i, 'link'] = f"https://venus.sgc.ox.ac.uk/venus?uniprot={ref.at[i, 'uniprot']}&species=9606&mutation={ref.at[i, 'uniprot_mutation']}"
```

Many of the clean-up steps will be database dependent, e.g. splitting up compound-hets 
or `ref = ref.loc[:,~ref.columns.duplicated()].copy()` to remove duplicate columns etc etc.

Regarding cmp-htz it may be worth storing that detail though:
```python
ref['is_dominant'] = ~((ref.zygosity in ('hom', 'homo', 'hmz')) | (ref['Case Id'].duplicated(keep=False)))
```
(parenthetically htz/hmz are the recommended abbreviations, not homo/het)