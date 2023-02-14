> Not finished moving

## Other



```python
import os
from michelanglo_protein import ProteinCore
ProteinCore.settings.startup(data_folder=os.environ['MICHELANGLO_PROTEIN_DATA'])

import json

with open(os.path.join(ProteinCore.settings.dictionary_folder, 'taxid9606-names2uniprot.json'), 'r') as fh:
    human = json.load(fh)
table['uniprot'] = pd.Series(zip(table.uniprot, table.gene.str.split(' ',expand=True)[0].map(human))).apply(lambda ab: ab[0] if str(ab[0]) != 'nan' else ab[1])
```

```python
*(map(operator.itemgetter('value'), warehouse.proteins)),

# make dictionaries of uni, mut -> value
import functools, operator
ref = table
get_sub: Callable[[pd.Series,], Dict[Tuple[str, str], Any]] = lambda series: dict(zip(zip(ref.uniprot.values, ref.uniprot_mutation.values), series))
to_pLI = get_sub(ref.pLI)
#to_disease = get_sub(ref.Disease)
to_uniprot = get_sub(ref.uniprot)
to_note = get_sub(ref.Notes.fillna(''))
to_case = get_sub(ref.pid)
to_previous = get_sub(ref.previous_page.fillna(''))
onto = ref.classification\
           .apply(lambda v: v if isinstance(v, set) else {Ontology.UNKNOWN})\
           .apply(functools.partial(map, operator.attrgetter('name')))\
           .apply(functools.partial(map, str.lower))\
           .apply(','.join)
to_csq = get_sub(onto)
is_dominant = dict(zip(zip(ref.uniprot.values, ref.uniprot_mutation.values), ref.Dominant))
```

```python
class Quick:
    
    def __init__(self, df: pd.DataFrame,
                 gene_column='gene',
                 mutation_column='uniprot_mutation',
                 entry_column='classification'):
        self.df = df
        self.gene_column = gene_column
        self.mutation_column = mutation_column
        self.entry_column = entry_column
        
    def get_row(self, gene, mutation):
        mask: pd.Series = (self.df[self.gene_column] == gene) & (self.df[self.mutation_column] == mutation)
        matched: pd.DataFrame = self.df.loc[mask]
        assert len(matched)
        return matched.reset_index().iloc[0]

    def get_index(self, gene, mutation):
        first = self.get_row(gene, mutation)
        return first['index']

    def set(self, gene, mutation, value):
        index = self.get_index(gene, mutation)
        self.df.at[index, self.entry_column] = value

    def add(self, gene, mutation, value):
        index = self.get_index(gene, mutation)
        self.df.at[index, self.entry_column].add(value)
        
    def remove(self, gene, mutation, value):
        index = self.get_index(gene, mutation)
        self.df.at[index, self.entry_column].remove(value)
        
q = Quick(table, entry_column='classification')
q.remove('RTTN','p.Val891fs', Ontology.SURFACE)
q.add('RTTN','p.Val891fs', Ontology.NOT_MISSENSE)


```