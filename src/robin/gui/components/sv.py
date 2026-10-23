"""Structural-variant sample viewer, backed by a read-only SQLite snapshot."""
import csv
import io
import logging
from pathlib import Path

from nicegui import run, ui
from robin.gui.theme import styled_table, client_timer
from .sv_data import load_page, load_evidence

logger=logging.getLogger(__name__)


def add_sv_section(launcher, sample_dir):
    sample_dir=Path(sample_dir)
    with ui.element('div').classes('classification-insight-shell w-full min-w-0'):
        ui.label('Structural variants').classes('classification-insight-heading text-headline-small')
        ui.label('Experimental SV evidence. Local-anchor HQ requires 100 consecutive aligned CIGAR bases on each side; it does not prove unique placement or override repeat masks. Provisional means sufficient unmasked read support; '
                 'it does not establish somatic status. Masked calls are retained for review.').classes('classification-insight-meta')
        ui.label('Breakpoints use zero-based interbase coordinates. Inversion-like and duplication-like describe single junctions; complete events and copy counts are not established. Select a call to highlight its boundaries on the CNV profile for 3 seconds; use Supporting reads for read IDs.').classes('text-xs text-grey-7')
        summary=ui.label('Loading SV results…').classes('classification-insight-meta')
        with ui.row().classes('w-full items-center gap-3'):
            status=ui.select({'all':'All statuses','provisional':'Provisional','masked_candidate':'Masked candidate',
                'candidate':'Below-threshold candidate','reference_boundary':'Reference boundary'},value='all',label='Status').classes('w-56')
            mask=ui.select({'all':'Any mask','unmasked':'Unmasked','tandem_repeat':'Tandem repeat',
                'centromere_satellite':'Centromere / satellite','chromosome_end_margin':'Chromosome end',
                'custom_mask':'Custom mask','mitochondrial_reference_boundary':'Mitochondrial origin'},value='all',label='Region').classes('w-56')
            minimum=ui.number('Minimum supporting IDs',value=3,min=0,step=1).classes('w-48')
            search=ui.input('Chromosome or call ID',placeholder='e.g. chr9').props('debounce=500').classes('w-48')
        show_low_mapq=ui.checkbox('Show low mapping quality calls (MAPQ<55)',value=False)
        ui.label('Uses mean minimum MAPQ across supporting reads. Calls awaiting MAPQ measurements remain visible.').classes('text-xs text-grey-7')
        columns=[dict(name=k,label=label,field=k,align='left') for k,label in [
            ('id','Call ID'),('kind','Signature'),('topology','Junction topology'),('interval_bp','Interval (bp)'),('breakend1','Breakend 1'),('breakend2','Breakend 2'),
            ('size','Size (bp; 0 = unassigned)'),('support','Supporting IDs'),('high_quality_support','Local-anchor HQ IDs'),('mean_min_mapq','Mean minimum MAPQ'),('mapping_measured_ids','Mapping measured IDs'),('substitution_measured_segments','Error measured segments'),('confidence_score','Evidence score'),('confidence_reason','Info'),('status','Evidence status'),('mask_flags','Mask reasons')]]
        numeric={'id','size','interval_bp','support','high_quality_support','confidence_score','mapq30_support','mean_min_mapq','substitution_rate_mean','mapping_measured_ids','substitution_measured_segments','breakpoint1_sd','breakpoint2_sd','event_size_sd'}
        for column in columns:
            if column['name'] in numeric:column['align']='right'
        _,table=styled_table(columns=columns,rows=[],row_key='id',pagination=0,selection='single')
        table.add_slot('body-cell-confidence_score', r'''<q-td :props="props">
          <div style="min-width:150px" role="meter" aria-label="Experimental SV evidence score"
            aria-valuemin="0" aria-valuemax="100" :aria-valuenow="props.row.confidence_score">
            <div class="text-caption text-right">{{ props.row.confidence_score }} / 100</div>
            <q-linear-progress :value="props.row.confidence_score / 100" color="green-6"
              track-color="grey-3" rounded size="12px" />
            <q-tooltip>{{ props.row.confidence_reason }} Experimental ranking, not a probability.</q-tooltip>
          </div>
        </q-td>''')
        table.add_slot('body-cell-confidence_reason', r'''<q-td :props="props" @click.stop>
          <q-btn flat round dense icon="info_outline" aria-label="Show confidence reasons">
            <q-tooltip style="max-width:480px; white-space:normal">{{ props.row.confidence_reason }}</q-tooltip>
            <q-menu><div class="q-pa-md" style="max-width:480px; white-space:normal">
              <div class="text-subtitle2 q-mb-sm">Confidence reasons</div>
              <div>{{ props.row.confidence_reason }}</div>
              <div class="q-mt-sm">MAPQ ≥30 IDs: {{ props.row.mapq30_support ?? 'Unknown' }}</div>
              <div>Substitution rate: {{ props.row.substitution_rate_mean == null ? 'Unknown' : (100 * props.row.substitution_rate_mean).toFixed(2) + '%' }}</div>
              <div>Breakpoint 1 SD: {{ props.row.breakpoint1_sd == null ? 'Unknown' : props.row.breakpoint1_sd.toFixed(1) + ' bp' }}</div>
              <div>Breakpoint 2 SD: {{ props.row.breakpoint2_sd == null ? 'Unknown' : props.row.breakpoint2_sd.toFixed(1) + ' bp' }}</div>
              <div>Event size SD: {{ props.row.event_size_sd == null ? 'Unknown' : props.row.event_size_sd.toFixed(1) + ' bp' }}</div>
              <div class="text-caption q-mt-sm">Experimental evidence ranking, not a probability.</div>
            </div></q-menu>
          </q-btn>
        </q-td>''')
        with ui.row().classes('items-center gap-3'):
            previous=ui.button('Previous');page_label=ui.label('Page 1');following=ui.button('Next')
            reload_button=ui.button('Refresh',icon='refresh')
            download=ui.button('Download this page',icon='download')
            inspect_reads=ui.button('Supporting reads',icon='list').props('disable')
        message=ui.label('').classes('text-xs text-grey-7')
        with ui.dialog() as dialog, ui.card().classes('w-full max-w-4xl'):
            title=ui.label('Supporting reads').classes('text-h6')
            evidence_note=ui.label('')
            _,evidence=styled_table(columns=[dict(name=k,label=label,field=k,align='left') for k,label in
                [('read_id','Read ID'),('p1','Observed breakend 1'),('p2','Observed breakend 2'),('anchor1','Anchor 1 (bp)'),('anchor2','Anchor 2 (bp)'),('mapq1','MAPQ 1'),('mapq2','MAPQ 2'),('error1','Substitution rate 1'),('error2','Substitution rate 2'),('repeat_context','Repeat intervals / edge distances')]],
                rows=[],row_key='read_id',pagination=20)
            ui.button('Close',on_click=dialog.close)
        table.add_slot('header', r'''<q-tr :props="props">
          <q-th auto-width></q-th>
          <q-th v-for="col in props.cols" :key="col.name" :props="props"
            :style="['id','size','interval_bp','support','high_quality_support','confidence_score','mapq30_support','mean_min_mapq','substitution_rate_mean','mapping_measured_ids','substitution_measured_segments','breakpoint1_sd','breakpoint2_sd','event_size_sd'].includes(col.name) ? 'cursor:pointer' : ''"
            @click="['id','size','interval_bp','support','high_quality_support','confidence_score','mapq30_support','mean_min_mapq','substitution_rate_mean','mapping_measured_ids','substitution_measured_segments','breakpoint1_sd','breakpoint2_sd','event_size_sd'].includes(col.name) && $parent.$emit('numeric_sort', col.name)">
            {{ col.label }} <span v-if="['id','size','interval_bp','support','high_quality_support','confidence_score','mapq30_support','mean_min_mapq','substitution_rate_mean','mapping_measured_ids','substitution_measured_segments','breakpoint1_sd','breakpoint2_sd','event_size_sd'].includes(col.name)">↕</span>
          </q-th>
        </q-tr>''')
        ui.label('Evidence score v2: support /30, MAPQ≥30 on both sides /25, local anchors /10, whole-segment substitution burden /20, breakpoint and size consistency /15. Missing metrics cap at 69; one read caps at 25, two at 39. Repeat/centromere tags have no penalty. Experimental ranking, not a probability; partial backfill is labelled by measured counts.').classes('text-xs text-grey-7')
        state=dict(page=1,busy=False,rows=[],sort_by='support',descending=True)

        async def refresh():
            if state['busy']:return
            state['busy']=True
            try:
                result=await run.io_bound(load_page,sample_dir,status=status.value,mask=mask.value,
                    search=search.value or '',min_support=max(0,int(minimum.value or 0)),page=state['page'],sort_by=state['sort_by'],descending=state['descending'],show_low_mapq=show_low_mapq.value)
                if result is None:return
                state['page']=result['page'];state['rows']=result['rows']
                table.rows=result['rows'];table.selected=[];table.update()
                state['selected']=None;inspect_reads.set_enabled(False)
                counts=result['counts']
                summary.set_text(' · '.join(f'{value:,} {key.replace("_"," ")}' for key,value in counts.items())
                    if result['available'] else 'No SV checkpoint yet. This panel will refresh automatically.')
                page_label.set_text(f"Page {result['page']} of {result['pages']} · {result['total']:,} matching calls")
                previous.set_enabled(result['page']>1);following.set_enabled(result['page']<result['pages'])
                download.set_enabled(bool(result['rows']))
                message.set_text('Showing up to 50 calls per page. Click numeric headers to sort all matching calls; click again to reverse. Set minimum support to 0 for low-support candidates. Refreshes every 10 seconds.')
            except Exception:
                logger.exception('Could not refresh SV viewer')
                message.set_text('SV data is temporarily unavailable; keeping the previous view. Try Refresh.')
            finally:state['busy']=False

        async def sort_numeric(event):
            key=event.args
            if key not in numeric:return
            state['descending']=not state['descending'] if state['sort_by']==key else True
            state['sort_by']=key;state['page']=1
            await refresh()
        table.on('numeric_sort',sort_numeric)

        async def change_filter():
            state['page']=1;await refresh()
        for control in [status,mask,minimum,search,show_low_mapq]:control.on_value_change(change_filter)
        async def move(delta):
            state['page']=max(1,state['page']+delta);await refresh()
        previous.on_click(lambda:move(-1));following.on_click(lambda:move(1));reload_button.on_click(refresh)

        def select_row(row):
            state['selected']=row
            inspect_reads.set_enabled(row is not None)
            if row is None:return
            callback=getattr(ui.context.client,'_sv_cnv_highlighters',{}).get(str(sample_dir))
            if callback:callback(row)
            else:ui.notify('CNV profile is not available on this page',type='info')

        async def inspect_selected():
            row=state.get('selected')
            if row is None:return
            try:
                data=await run.io_bound(load_evidence,sample_dir,row['id'])
                if data is None:return
                title.set_text(f"SV {row['id']} · supporting read IDs")
                evidence_note.set_text(f"Showing {len(data['rows']):,} of {data['total']:,} supporting IDs. "
                                      'Related parent/duplex molecules are not yet reconciled.')
                evidence.rows=data['rows'];evidence.update();dialog.open()
            except Exception:
                logger.exception('Could not load SV evidence');ui.notify('Could not load SV evidence',type='warning')
        def selected(event):
            select_row(event.selection[0] if event.selection else None)
        def row_clicked(event):
            args=event.args
            clicked=args[1] if isinstance(args,(list,tuple)) and len(args)>1 else args
            if not isinstance(clicked,dict):return
            row=next((r for r in state['rows'] if r['id']==clicked.get('id')),None)
            if row is None:return
            table.selected=[row];table.update();select_row(row)
        table.on_select(selected)
        table.on('rowClick',row_clicked,[[],['id'],None])
        inspect_reads.on_click(inspect_selected)

        def export_page():
            if not state['rows']:return
            authorize=getattr(launcher,'_require_export_or_notify',None)
            if not callable(authorize) or not authorize():return
            output=io.StringIO();writer=csv.DictWriter(output,fieldnames=list(state['rows'][0]),delimiter='\t')
            writer.writeheader();writer.writerows(state['rows'])
            ui.download(output.getvalue().encode(),filename=f'{sample_dir.name}_sv_page_{state["page"]}.tsv')
        download.on_click(export_page)
        client_timer(10,refresh,immediate=True)
