"""Versioned experimental evidence ranking, not a calibrated probability."""
import math

def confidence(row):
    n=max(0,row['support']);hq=row.get('high_quality_support')
    measured=min(n,row.get('mapping_measured_ids') or 0)
    strong=min(measured,row.get('mapq30_support') or 0)
    errors=min(2*n,row.get('substitution_measured_segments') or 0)
    rate=row.get('substitution_rate_mean')
    support=30*min(1,math.log1p(n)/math.log(11))
    mapping=25*strong/n if n else 0
    anchors=10*min(n,hq or 0)/n if n else 0
    error=20*max(0,1-rate/.08)*errors/(2*n) if n and rate is not None else 0
    spread=[row.get('breakpoint1_sd'),row.get('breakpoint2_sd')]
    consistency=0
    if measured>=3 and all(v is not None for v in spread):
        variability=max(spread)
        size_sd=row.get('event_size_sd');mean_size=row.get('event_size_mean')
        if size_sd is not None and mean_size is not None:
            variability=max(variability,100*size_sd/max(50,.1*mean_size))
        consistency=15*max(0,1-variability/100)*measured/n
    cap=100;notes=[]
    if n<3:cap=39 if n==2 else 25;notes.append('fewer than 3 supporting IDs')
    if measured<n or errors<2*n:cap=min(cap,69);notes.append('mapping/error measurements incomplete')
    if strong<3 or (measured and strong/measured<.5):
        cap=min(cap,69);notes.append('insufficient MAPQ>=30 support')
    if rate is not None and rate>.05:
        cap=min(cap,69);notes.append('elevated whole-segment substitution burden')
    if measured>=3 and any(v is not None and v>50 for v in spread):
        cap=min(cap,69);notes.append('breakpoint dispersion above 50 bp SD')
    if measured>=3 and row.get('event_size_sd') is not None and row.get('event_size_mean') is not None and row['event_size_sd']>max(50,.1*row['event_size_mean']):
        cap=min(cap,69);notes.append('event-size dispersion above tolerance')
    flags=set(filter(None,row.get('mask_flags','').split(',')))
    if flags-{'tandem_repeat','centromere_satellite'}:cap=min(cap,69);notes.append('other region flags')
    if row.get('status')=='reference_boundary':cap=5;notes.append('reference-boundary flag')
    score=round(min(cap,support+mapping+anchors+error+consistency))
    reason=(f'Support {support:.1f}/30; MAPQ>=30 on both sides {mapping:.1f}/25; '
            f'local anchors {anchors:.1f}/10; substitution burden {error:.1f}/20; '
            f'breakpoint/size consistency {consistency:.1f}/15; cap {cap}. ')
    if hq is None:notes.append('local-anchor quality not measured')
    if rate is None:notes.append('substitution burden unknown')
    if flags & {'tandem_repeat','centromere_satellite'}:notes.append('repeat/centromere annotations retained without score penalty; unique placement not established')
    return dict(confidence_score=score,confidence='high' if score>=70 else 'medium' if score>=40 else 'low',
                confidence_reason=reason+'; '.join(notes),score_version='mapping_v2')
