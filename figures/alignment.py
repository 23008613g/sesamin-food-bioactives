"""Project-local plot-rectangle check; no third-party audit tool is vendored.

This small check reports physical panel geometry and same-row alignment. It is
not the fuller historical figure-tool collision/glyph audit. Those earlier
checks are described separately in the manuscript execution history.
"""
from pathlib import Path
import json

def require_matplotlib_panel_alignment(fig, axes, panel_ids=None, json_out=None,
        overlay_svg=None, tolerance_pt=1.5, gutter_tolerance_pt=1.5,
        require_panel_labels=True, strict=True, exemptions=(), **unused):
    fig.canvas.draw()
    ids=panel_ids or list('abcdefghijklmnopqrstuvwxyz')[:len(axes)]
    width,height=fig.get_size_inches()*72
    panels=[]
    for name,ax in zip(ids,axes):
        x,y,w,h=ax.get_position().bounds
        panels.append({'id':name,'left_pt':float(x*width),'bottom_pt':float(y*height),
                       'width_pt':float(w*width),'height_pt':float(h*height),
                       'top_pt':float((y+h)*height),
                       'panel_label_present':any(t.get_text()==name for t in ax.texts)})
    failures=[]
    comparisons=[]
    for i,a in enumerate(panels):
        if require_panel_labels and not a['panel_label_present']:
            failures.append('Missing panel label '+a['id'])
        for b in panels[i+1:]:
            if abs(a['top_pt']-b['top_pt'])<=tolerance_pt:
                error=abs(a['bottom_pt']-b['bottom_pt'])
                comparisons.append({'panels':[a['id'],b['id']],'bottom_difference_pt':error})
                if error>tolerance_pt:failures.append('Unequal row bottoms: '+a['id']+b['id'])
                exempt=any('panel-width' in e.get('checks',[]) and
                           {a['id'],b['id']}.issubset(e.get('panels',[])) for e in exemptions)
                if not exempt and abs(a['width_pt']-b['width_pt'])>tolerance_pt:
                    failures.append('Unequal row widths without exemption: '+a['id']+b['id'])
    result={'status':'PASS' if not failures else 'FAIL','audit':'project-local rectangle/label check',
            'width_mm':float(width/72*25.4),'height_mm':float(height/72*25.4),
            'panels':panels,'comparisons':comparisons,'exemptions':list(exemptions),'failures':failures}
    if json_out:
        Path(json_out).parent.mkdir(parents=True,exist_ok=True)
        Path(json_out).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    if overlay_svg:
        # Exporting a fresh plot is useful for manual inspection. This is not a
        # claimed collision overlay from the excluded third-party audit tool.
        fig.savefig(overlay_svg)
    if strict and failures:raise AssertionError('; '.join(failures))
    return result
