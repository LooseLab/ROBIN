"""Pure construction of temporary SV boundary overlays for CNV charts."""

SERIES_NAME = 'SV selected boundaries'

def boundary_series(row, cnv_map, selected):
    offsets={};offset=0
    for chrom,values in sorted(cnv_map.items(),key=lambda item: (int(item[0][3:]) if item[0].startswith('chr') and item[0][3:].isdigit() else 23 if item[0]=='chrX' else 24 if item[0]=='chrY' else 999, item[0])):
        # Match the CNV renderer's chromosome acceptance and concatenated bins.
        if chrom not in {*(f'chr{i}' for i in range(0,23)), 'chrX','chrY'}:
            continue
        offsets[chrom]=offset;offset+=len(values)
    points=[]
    for chrom,pos in [(row['c1'],row['p1']),(row['c2'],row['p2'])]:
        if selected!='All' and chrom!=selected:continue
        if selected=='All' and chrom not in offsets:continue
        x=pos+(offsets[chrom]*row['bin_width'] if selected=='All' else 0)
        points.append(x)
    midpoint_line=row['c1']==row['c2'] and abs(row['p2']-row['p1'])<10000
    width=max(1,row['bin_width']*.05)
    if row['c1']==row['c2'] and len(points)==2:
        start,end=sorted(points)
        ranges=[((start+end)/2,(start+end)/2)] if midpoint_line else [(start,end)]
    else:
        ranges=[(max(0,x-width),x+width) for x in points]
    return dict(id='robin-sv-highlight',midpoint_line=midpoint_line,name=SERIES_NAME,type='line',data=[],silent=True,
                markArea=dict(silent=True,animation=False,label=dict(show=False),
                              itemStyle=dict(color='#8b5cf6',opacity=row.get('opacity',.22),borderWidth=0),
                              data=[[dict(xAxis=start),dict(xAxis=end)] for start,end in ranges]))
