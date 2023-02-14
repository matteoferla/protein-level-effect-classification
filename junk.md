## Junk

The notes herein may be been incorrectly polished and may call some function that are decorative. For example:

```python
from typing import Optional, Dict, Any
from IPython.display import display, HTML

def entitle(title, 
            caption:Optional[str]=None, 
            level:int=1,
            styles:Optional[Dict[str, Any]]=None) -> None:
    style: str = ''
    if styles:
        style: str = 'style="' + ' '.join([f'{k}:{v};' for k, v in styles.items()]) +'"'
    html: str = f'<h{level} {styles}>{title}</h{level}>'
    if caption:
        html += f'<p>{caption}<p>'
    display(HTML(html))
```