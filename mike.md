## Upload to Michelanglo

Aim: have a warehouse page and a viewing page

### Warehouse
Make warehouse
```python
from michelanglo_api import MikeAPI, MikePage, Prolink

mike = MikeAPI('👾👾👾', '🤖🤖🤖')
warehouse = mike.convert_pdb(code='1ubq', filename='')
warehouse.retrieve()
```
Subsequent runs:
```python
warehouse = mike.get_page('👾👾-👾👾-👾👾-👾👾-👾👾👾')
```
If you lost the link to `warehouse`...

```python
for page in mike.owned_pages:
    page.retrieve()
    if 'warehouse' in page.title.lower():
        break
warehouse = page
warehouse.retrieve()
print(warehouse.title)
```
Edit warehouse

```python
warehouse.title = 'Structure warehouse'
warehouse.description = 'This will likely crash your browser. Not meant for general use.'
warehouse.commit()
warehouse.show_link()
```
# get for first time
```python
import re
warehouse.retrieve()
warehouse.protein = []
warehouse.pdbs = {}
broken = []
for uniport, mutation in results:
    result = results[(uniport, mutation)]
    if 'ddG' not in result:
        broken.append((uniport, mutation))
        continue
    forename = f'{uniport}_{mutation.replace(".", "")}'
    assert forename.replace("_", "").isalnum(), f'{forename} has some random chars... how??'
    result['warehouse_name'] = forename
    for aftname in ('native', 'mutant'):
        name = f'{forename}_{aftname}'
        warehouse.append_pdbblock(varname=name,
                                  pdbblock=result['ddG'][aftname],
                                  chain_definitions=result['structural']['chain_definitions'])
print(broken)
warehouse.commit()
```
Add power to the warehouse by monkeypatching it:
```python
def get_i(self, name):
    for i, entry in enumerate(self.proteins):
        if entry['value'] == name:
            return i
    else:
        raise Exception
        
def get_pdb_url(self, name):
    i = self.get_i(name)
    return f'/save_pdb?uuid={self.page}&key=None&index={i}'
        
warehouse.get_i = types.MethodType(get_i, warehouse)
warehouse.get_pdb_url = types.MethodType(get_pdb_url, warehouse)
```

### Gallery
```python
gallery = mike.convert_pdb(code='1ubq')
gallery.retrieve()
gallery.show_link()
```

monkeypatcha:
```python
import types

def append_url(self, name:str, url:str, ext:str='pdb', verify:bool=True):
    if verify: # test url
        response = requests.get(url)
        response.raise_for_status()
    self.proteins.append(dict(name=name,
                          type='url',
                          ext=ext,
                          value=url)
                        )
    
gallery.append_url = types.MethodType(append_url, gallery)

# blank
gallery.loadfun = ''
gallery.description = ''
gallery.proteins = []
gallery.pdbs = {}
gallery.image = False  # url
gallery.columns_text = 6
gallery.title='HICF2 Variants'
gallery.append_url = types.MethodType(append_url, gallery)
gallery.commit()
```

Now, for testing I suggest making a dummy when needed:

```python
import copy

fauxllery = mike.convert_pdb(code='1ubq')
gallery = mike.get_page('👾👾-👾👾-👾👾-👾👾-👾👾👾')
fauxllery.__dict__.update(**gallery.__dict__)
fauxllery.proteins = copy.deepcopy(gallery.proteins)
fauxllery.title='Debug'
fauxllery.commit()
```
Something
```python
from michelanglo_api import MikeAPI, MikePage, Prolink

mike = MikeAPI('👾👾👾', '👾👾👾')

is_dominant = dict(zip(zip(ref.uniprot.values, ref.uniprot_mutation.values), ref.Dominant))
to_pLI = dict(zip(zip(ref.uniprot.values, ref.uniprot_mutation.values), ref.pLI))

# make_fragment_table is gen'ed dynamically... this isnt

class Dormant(Prolink):
    """
    'Dormant' because it is not awake on load
    """
    def __str__(self):
        text = super().__str__()
        return text.replace('class="prolink" data-toggle="protein"', 
                            'type="button" class="btn btn-info dormant_prolink"')\
                   .replace('<span', '<button')\
                   .replace('</span', '</button')

#gallery.title='HICF2 Variants'
gallery.proteins = []
gallery.pdbs = {}

gallery.loadfun = '''

$('.dormant_prolink').click(event => {
    // Convert the data and then click.
    NGL.getStage().removeAllComponents();
   let jq = $(event.currentTarget);
   let data = jq.data();
   let wanted = [data.load];
   if (data.focus.indexOf('overlay') !== -1) {
       wanted.push(data.focus.split(' ')[1]);
   }
    // get only the named ones
    // get only the urls
    // the item returned by $.get is a deferred (faux promise)
    // done does the same
   let thenables = myData.proteins
                     .filter(e => wanted.includes(e.name) )
                     .filter(e => e.type === 'url')
                     .map(entry => $.get(entry.value).done(text => {entry.value = text; entry.type='data'})  );
   Promise.all(thenables).then(() => NGL.specialOps.prolink(event.currentTarget));
});

'''

gallery.description = '## HICF2 Structures\n'
gallery.description += 'Clicking on the icon for wild-type/mutant/overlay will load '+\
                      'and show those structures. '+\
                      'A PDB file is generally a few megabytes, hence the delay.\n'
gallery.description += '''
<div class="alert alert-primary" role="alert">
  Button legend: W=wild type, M=mutant, O=overlay.
</div>
'''

gallery.description += '<table class="table" width="100%">\n' 
gallery.description += '<thead class="thead-light"><tr>'
gallery.description += '<th scope="col">Gene</th>'
gallery.description += '<th scope="col">Mutation</th>'
gallery.description += '<th scope="col">Dominant genotype</th>'
gallery.description += '<th scope="col">pLI</th>'
gallery.description += '<th scope="col">∆∆G</th>'
gallery.description += '<th scope="col">Show</th>'
gallery.description += '<th scope="col">Venus</th>'
gallery.description += '</tr></thead>'
gallery.description += '<tbody>'
for uniport, mutation in results:
    result: Dict[str, Any] = results[(uniport, mutation)]
    gallery.description += '<tr>'
    gene = result['protein']['gene_name']
    
    print(gene)
    
    gallery.description += f'<td>{gene}</td>'
    gallery.description += f'<td>{mutation}</td>'
    dominance = {True: '<i class="far fa-check"></i>', 
                 False: '<i class="far fa-times"></i>', 
                 None: '<i class="far fa-question"></i>', }\
                [is_dominant.get((uniport, mutation), None)]
    pLI = to_pLI.get((uniport, mutation), None)
    gallery.description += f'<td>{dominance}</td>'
    gallery.description += f'<td>{pLI:.2f}</td>'
    if 'ddG' not in result:
        broken.append((uniport, mutation))
        gallery.description += f'<td><i class="far fa-times"></i></td>' * 4
        gallery.description += '</tr>'
        continue
    sele = f"{result['mutation']['residue_index']}:A"
    gallery.description += f"<td>{result['ddG']['ddG']:.1f}</td>"
    gallery.description += '<div class="btn-group" role="group"><td>'
    for aftname in ('native', 'mutant'):
        forename = result['warehouse_name']
        name = f'{forename}_{aftname}'
        gallery.append_url(name=name, url=warehouse.get_pdb_url(name), verify=False)
        #'<i class="fas fa-crosshairs"></i>'
        pro = Dormant(text={'native': 'W', 'mutant': 'M'}[aftname],
                       load=name,
                       title=aftname,                      
                       focus='residue',
                       selection=sele,
                       hetero=True,
                       )
        gallery.description += str(pro)
    pro = Dormant(text='O',
                       load=f'{forename}_native',
                       title='overlay',                      
                       focus=f'overlay {forename}_mutant',
                       selection=sele,
                       hetero=True,
                       )
    gallery.description += str(pro)
    gallery.description += '</div></td>'
    
    gallery.description += f'''<td>
    <a href="https://michelanglo.sgc.ox.ac.uk/venus?uniprot={uniport}&species=9606&mutation={mutation}" 
    target="_blank" type="button" class="btn btn-success">
    <i class="far fa-calculator"></i>
    </a></td>'''
    gallery.description += '</tr>\n'
gallery.description += '</tbody>'
gallery.description += '</table>'
gallery.commit()
print('done')
```