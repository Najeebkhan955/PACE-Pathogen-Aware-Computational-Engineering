#!/usr/bin/env python3
"""Extract original embedded PACE figures from a DOCX without resizing.
Usage: python extract_pace_figures.py Manuscript.docx --out figures
Requires: python-docx
"""
import argparse,csv,re,hashlib
from pathlib import Path
from docx import Document
from docx.oxml.ns import qn
ap=argparse.ArgumentParser()
ap.add_argument('manuscript')
ap.add_argument('--out',default='figures')
args=ap.parse_args()
out=Path(args.out)
out.mkdir(parents=True,exist_ok=True)
doc=Document(args.manuscript)
rows=[]
seen=set()
for i,para in enumerate(doc.paragraphs):
    blips=para._p.xpath('.//a:blip')
    if not blips:
        continue
    matched=None
    for next_para in doc.paragraphs[i:i+7]:
        m=re.match(r'\s*Figure\s+(\d+)\s*[.:]',next_para.text,re.I)
        if m:
            matched=(int(m.group(1)),next_para.text.strip())
            break
    if not matched:
        continue
    num,caption=matched
    for k,blip in enumerate(blips,1):
        rid=blip.get(qn('r:embed'))
        if not rid or (num,rid) in seen:
            continue
        seen.add((num,rid))
        part=doc.part.related_parts[rid]
        extension=Path(part.partname).suffix or '.bin'
        filename=f'PACE_Figure_{num:02d}'+(f'_part{k}' if len(blips)>1 else '')+extension
        data=part.blob
        (out/filename).write_bytes(data)
        rows.append(dict(figure=num,file=filename,bytes=len(data),
                         sha256=hashlib.sha256(data).hexdigest(),caption=caption))
with (out/'FIGURE_MANIFEST.csv').open('w',newline='',encoding='utf-8') as file:
    writer=csv.DictWriter(file,fieldnames=['figure','file','bytes','sha256','caption'])
    writer.writeheader()
    writer.writerows(rows)
print(f'Extracted {len(rows)} original figures to {out}')
