```python
import pandas as pd
from enum import Enum
from typing import Set

class Ontology(Enum):
    """
    Version 2. copied from Mike_clean
    """
    NOT_MISSENSE = -1
    UNKNOWN = 0
    DESTABILIZATION = 1
    DESTABILIZATION_CAVITY = 1.1
    DESTABILIZATION_CLASH = 1.2
    DESTABILIZATION_BBTORSION = 1.3
    ALTERED_ENERGETICS = 1.4
    LIGAND = 2
    MISLOCALISATION = 3
    PTM = 4
    PORE = 5
    INTERFACE = 6
    INTERFACE_DNA = 6.1
    INTERFACE_INTERDOMAIN = 6.2
    # qualitative...
    BURIED = 10
    SURFACE = 11
    CHARGE_CHANGE = 12
    GNOMAD = 13
    
    @classmethod
    def ontologize(cls, row:pd.Series):
        """
        Reads multiple cells to make a verdict...
        e.g. ``distance_to_closest_ligand``, ``closest_ligand``, ``ddG`` and ``rsa``
        """
        properties: Set[cls] = set()
        if not row.is_missense:
            return {cls.NOT_MISSENSE}
        if row.distance_to_closest_ligand < 10 and \
            not any([v in row.closest_ligand for v in ('UNK', 'UNL')]):
            print(row.closest_ligand)
            properties.add(cls.LIGAND)
        properties.add(cls.BURIED if row.rsa < 0.3 else cls.SURFACE)
        if row.ddG > 1.5:
            properties.add(cls.DESTABILIZATION)
        return properties if properties else {cls.UNKNOWN}
    
    @classmethod
    def reontologize(cls, value: str, sep='; '):
        """
        Given a string from a cell, return a set of Ontologies
        """
        if not isinstance(value, str):
            return set()
        return {cls[v.upper()] for v in value.split(sep)}
```
Go!
```python
ref['ontology'] = ref.apply(ontologize, axis=1)
```
But stuff may be wrong...
```python
def add_onto(gene, mutation, value):
    index = get_index(gene, mutation)
    ref.at[index, 'ontology'].add(value)
    
add_onto('ATP8B1', 'p.Arg934Gly', Ontology.DESTABILIZATION_CAVITY)
add_onto('CFI', 'p.Pro64Leu', Ontology.DESTABILIZATION_BBTORSION)
add_onto('KDM2B','p.Trp1102Ser', Ontology.DESTABILIZATION_CAVITY)
add_onto('PLPPR5','p.Leu13Phe', Ontology.MISLOCALISATION)
add_onto('COG7', 'p.Tyr500Cys', Ontology.UNKNOWN)  # substantial changes: pocket filled.
add_onto('PSTPIP1','p.Glu250Lys', Ontology.CHARGE_CHANGE)
add_onto('VHL', 'p.Arg200Trp', Ontology.CHARGE_CHANGE)  # surface hydrophobic
```