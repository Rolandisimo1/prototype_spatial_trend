#!/usr/bin/env python
"""EXECUTE the viewer's JavaScript the way a browser would, in QuickJS.

Headless browsers do not run in this sandbox, so this is the substitute for the check that would
have caught the v4 blank page: JSON.parse the payload under real JS strictness, eval the page
script, call render(), then drive every species x metric combination and every interaction, and
assert the DOM actually received content.

The DOM shim below is deliberately small but real enough for the page's own selectors:
class-based querySelectorAll, attribute selectors on data-cid, classList, dataset, innerHTML
reparsing (so elements created by innerHTML are queryable, which is how the map cells appear).
"""
import json, re, sys
import quickjs

SHIM = r"""
function Elem(tag, attrs){
  this.tagName=(tag||'div').toUpperCase();
  this.attrs=attrs||{};
  this._html='';
  this._text='';
  this._kids=[];        // element children created via appendChild
  this._parsed=[];      // element descendants created via innerHTML
  this._cls={};
  var c=this.attrs['class']||'';
  var self=this;
  c.split(/\s+/).forEach(function(k){ if(k) self._cls[k]=1; });
  this.dataset={};
  for(var k in this.attrs){
    if(k.indexOf('data-')===0){
      var camel=k.slice(5).replace(/-([a-z])/g,function(m,p){return p.toUpperCase();});
      this.dataset[camel]=this.attrs[k];
    }
  }
  this.style={};
  this.classList={
    add:function(k){ self._cls[k]=1; },
    remove:function(k){ delete self._cls[k]; },
    contains:function(k){ return !!self._cls[k]; }
  };
  this.onclick=null; this.onmousemove=null; this.onmouseleave=null; this.onmouseenter=null;
}
Elem.prototype.hasClass=function(k){ return !!this._cls[k]; };
Elem.prototype.setAttribute=function(k,v){ this.attrs[k]=String(v); };
Elem.prototype.getAttribute=function(k){ return this.attrs[k]; };
Elem.prototype.appendChild=function(el){ this._kids.push(el); return el; };
Elem.prototype.getBoundingClientRect=function(){ return {left:0,top:0,bottom:0,right:0}; };
Elem.prototype.descendants=function(){
  var out=[];
  this._parsed.forEach(function(e){ out.push(e); });
  this._kids.forEach(function(e){ out.push(e); e.descendants().forEach(function(x){out.push(x);}); });
  return out;
};
// minimal attribute scanner: our own generated markup, so tags are well formed
function parseHTML(html){
  var out=[], re=/<([a-zA-Z][\w-]*)((?:\s+[\w:-]+\s*=\s*"[^"]*")*)\s*\/?>/g, m;
  while((m=re.exec(html))!==null){
    var attrs={}, ar=/([\w:-]+)\s*=\s*"([^"]*)"/g, a;
    while((a=ar.exec(m[2]))!==null){ attrs[a[1]]=a[2]; }
    out.push(new Elem(m[1], attrs));
  }
  return out;
}
// className assignment. The page sets b.className='tab unsup' on tabs it creates via
// createElement, and the original shim only read the class attribute at construction time, so
// every such element looked class-less and the tab probes reported zero tabs. Setting className
// must update the class set the selector engine consults.
Object.defineProperty(Elem.prototype,'className',{
  get:function(){ return Object.keys(this._cls).join(' '); },
  set:function(v){
    var self=this; this._cls={};
    String(v).split(/\s+/).forEach(function(k){ if(k) self._cls[k]=1; });
    this.attrs['class']=String(v);
  }
});
Object.defineProperty(Elem.prototype,'innerHTML',{
  get:function(){ return this._html; },
  set:function(v){ this._html=String(v); this._parsed=parseHTML(this._html); this._kids=[]; }
});
Object.defineProperty(Elem.prototype,'textContent',{
  get:function(){ return this._text; },
  set:function(v){ this._text=String(v); }
});
function matchSel(el, sel){
  // supports  .a  .a.b  tag  .a[data-x="y"]
  var attr=null;
  var am=sel.match(/\[([\w-]+)="([^"]*)"\]$/);
  if(am){ attr=[am[1],am[2]]; sel=sel.slice(0, am.index); }
  var parts=sel.split('.').filter(function(s){return s.length;});
  for(var i=0;i<parts.length;i++){ if(!el.hasClass(parts[i])) return false; }
  if(attr && el.attrs[attr[0]]!==attr[1]) return false;
  return true;
}
Elem.prototype.querySelectorAll=function(sel){
  return this.descendants().filter(function(e){ return matchSel(e,sel); });
};
Elem.prototype.querySelector=function(sel){
  var r=this.querySelectorAll(sel); return r.length?r[0]:null;
};

var __ids={};
var document={
  _root:new Elem('body',{}),
  getElementById:function(id){
    if(!__ids[id]){ __ids[id]=new Elem('div',{id:id}); document._root.appendChild(__ids[id]); }
    return __ids[id];
  },
  createElement:function(t){ return new Elem(t,{}); },
  querySelectorAll:function(sel){
    var out=[];
    for(var k in __ids){ __ids[k].querySelectorAll(sel).forEach(function(e){out.push(e);}); }
    return out;
  },
  addEventListener:function(){}
};
var window={innerWidth:1400, innerHeight:900};
function setTimeout(f,ms){ return 0; }
var console={log:function(){}};
"""

PROBE = r"""
function __map(){ return document.getElementById('map'); }
function __n(sel){ return __map().querySelectorAll(sel).length; }
function __strip(s){ return String(s).replace(/<[^>]*>/g,''); }

/* Every species x measure x map-type combination. Asserts cells are drawn, the badge states the
   verdict this combination actually has in covariate_model_cv_fixed.csv (beats_fix), and that an
   unsupported surface is marked in a way a supported one is not. */
function __probeCombos(){
  var res=[], empty=0, bad_verdict=0, bad_mark=0;
  var types=['cov','position'];
  for(var t=0;t<types.length;t++){
    for(var i=0;i<D.sp_order.length;i++){
      for(var j=0;j<D.metrics.length;j++){
        MAPTYPE=types[t]; SP=D.sp_order[i]; MET=D.metrics[j].key;
        PICKED=null; LAYER='value'; SCALE='pooled'; CURVEMODE='deviation';
        render();
        var n = (MAPTYPE==='cov') ? __n('.covcell') : __n('.cell');
        if(n===0) empty++;
        var badge=__strip(document.getElementById('verdictbadge').innerHTML);
        var stamp=__strip(document.getElementById('mapstamp').innerHTML);
        var mh=__map().innerHTML;
        var rec=(D.cov.cv[SP]||{})[MET];
        var st=supportOf(SP,MET).state;
        // the badge word must match what beats_fix says for THIS combination
        var want={supported:'Supported', unsupported:'Not supported', weak:'Weakly supported',
                  untested:'Not tested'}[st];
        var vok = badge.indexOf(want)>=0;
        if(!vok) bad_verdict++;
        // an unsupported/weak surface must carry BOTH a stamp and a whole-surface hatch;
        // a supported one must carry neither
        var hasHatch = mh.indexOf('url(#hatchunsup)')>=0;
        var hasStamp = stamp.length>0;
        var markOK = (st==='supported') ? (!hasHatch && !hasStamp) : (hasHatch && hasStamp);
        if(!markOK) bad_mark++;
        // measure tabs must strike through an unsupported measure BEFORE it is clicked
        var mtabs=document.getElementById('mtabs');
        var nUnsup=mtabs.querySelectorAll('.tab.unsup').length;
        res.push({maptype:MAPTYPE, sp:SP, met:MET, cells:n, state:st,
                  badge_ok:vok, mark_ok:markOK, hatched:hasHatch, stamped:hasStamp,
                  beats_fix: rec? rec.beats : null,
                  unsup_tabs:nUnsup,
                  curve_drawn: document.getElementById('curves').innerHTML.indexOf('<svg')>=0,
                  dec_readout: document.getElementById('decreadout').innerHTML.indexOf('decbar')>=0,
                  cv_rows:(document.getElementById('cvtable').innerHTML.match(/<tr/g)||[]).length,
                  sup_rows:(document.getElementById('suptable').innerHTML.match(/suprow/g)||[]).length});
      }
    }
  }
  return {results:res, empty:empty, bad_verdict:bad_verdict, bad_mark:bad_mark};
}

/* Cell clicks on BOTH map types in BOTH curve modes, plus every toggle. A click must select a cell
   and redraw the curve panel; in deviation mode the panel must contain the two-panel deviation
   chart and its key, in absolute mode it must not. */
function __probeClicks(){
  var out=[];
  var types=['cov','position'], modes=['deviation','absolute'];
  for(var t=0;t<types.length;t++){
    for(var mo=0;mo<modes.length;mo++){
      for(var i=0;i<D.sp_order.length;i++){
        MAPTYPE=types[t]; SP=D.sp_order[i]; MET=D.metrics[0].key;
        PICKED=null; CURVEMODE=modes[mo]; LAYER='value'; SCALE='pooled';
        render();
        var sel = (MAPTYPE==='cov')?'.covcell':'.cell';
        var cells=__map().querySelectorAll(sel);
        if(!cells.length){ out.push({maptype:MAPTYPE, mode:modes[mo], sp:SP, err:'no cells'});
          continue; }
        var pick=cells[Math.floor(cells.length/2)];
        pick.onclick();
        var box=document.getElementById('curves').innerHTML;
        var nsel=__map().querySelectorAll(sel+'.sel').length;
        out.push({maptype:MAPTYPE, mode:modes[mo], sp:SP,
                  picked:String(PICKED).slice(0,14),
                  selected_on_map:nsel,
                  has_svg:box.indexOf('<svg')>=0,
                  has_devkey:box.indexOf('devkey')>=0,
                  has_two_panels:(box.match(/polygon/g)||[]).length>=2,
                  labelled:box.indexOf('nationwide rhythm')>=0,
                  clear_link:box.indexOf('clearpick')>=0});
      }
    }
  }
  return out;
}

/* Every control: map-type switch, curve-mode switch, value/uncertainty, pooled/species scale,
   species tabs, measure tabs, panel modal, clear-selection. */
function __probeToggles(){
  var log=[];
  function snap(tag){
    log.push({step:tag, maptype:MAPTYPE, mode:CURVEMODE, layer:LAYER, scale:SCALE,
              map_cells:__n('.cell')+__n('.covcell'),
              curve_svg:document.getElementById('curves').innerHTML.indexOf('<svg')>=0,
              legend:document.getElementById('legsw').innerHTML.length>0});
  }
  MAPTYPE='cov'; SP=D.sp_order[0]; MET=D.metrics[0].key; PICKED=null; CURVEMODE='deviation';
  LAYER='value'; SCALE='pooled'; render(); snap('start');
  document.getElementById('btnpos').onclick(); snap('btnpos');
  document.getElementById('btncov').onclick(); snap('btncov');
  document.getElementById('btnci').onclick();  snap('btnci');
  document.getElementById('btnval').onclick(); snap('btnval');
  document.getElementById('btnsp').onclick();  snap('btnsp');
  document.getElementById('btnpool').onclick();snap('btnpool');
  document.getElementById('btnabs').onclick(); snap('btnabs');
  document.getElementById('btndev').onclick(); snap('btndev');
  // species and measure tabs
  var sp_t=document.getElementById('sptabs').querySelectorAll('.tab');
  var m_t=document.getElementById('mtabs').querySelectorAll('.tab');
  var ntabs={species:sp_t.length, metrics:m_t.length};
  for(var i=0;i<sp_t.length;i++){ sp_t[i].onclick(); }
  snap('all_species_tabs');
  var m2=document.getElementById('mtabs').querySelectorAll('.tab');
  for(var j=0;j<m2.length;j++){ m2[j].onclick(); }
  snap('all_metric_tabs');
  // click a cell then clear it
  var c=__map().querySelectorAll('.covcell'); var cleared=null;
  if(c.length){ c[0].onclick();
    var cp=document.getElementById('clearpick');
    if(cp){ cp.onclick({preventDefault:function(){}}); cleared=(PICKED===null); } }
  snap('clear_selection');
  // the covariate-relationship modal
  var minis=document.getElementById('panels').querySelectorAll('.mini');
  var modal=null;
  if(minis.length){ minis[0].onclick();
    modal={open:document.getElementById('modal').classList.contains('on'),
           body:document.getElementById('modalbody').innerHTML.indexOf('<svg')>=0};
    document.getElementById('modalx').onclick();
    modal.closed=!document.getElementById('modal').classList.contains('on'); }
  return {log:log, ntabs:ntabs, cleared:cleared, modal:modal};
}

/* The deviation chart must actually show the species reference and a non-trivial departure, and the
   readout percentages must match the payload decomposition. */
function __probeDeviation(){
  var out=[];
  for(var i=0;i<D.sp_order.length;i++){
    SP=D.sp_order[i]; MAPTYPE='cov'; MET='pct_nocturnal'; CURVEMODE='deviation'; PICKED=null;
    render();
    var cells=__map().querySelectorAll('.covcell');
    var ci=+cells[Math.floor(cells.length/2)].dataset.ci;
    var cell=covCurve(SP,ci), mean=D.cov_mean[SP];
    var maxdev=0, sumabs=0;
    for(var k=0;k<cell.length;k++){ var d=Math.abs(cell[k]-mean[k]);
      if(d>maxdev)maxdev=d; sumabs+=d; }
    cells[Math.floor(cells.length/2)].onclick();
    var box=document.getElementById('curves').innerHTML;
    var ro=__strip(document.getElementById('decreadout').innerHTML);
    out.push({sp:SP, mean_len:mean.length, cell_len:cell.length,
              max_dev_pp:+(maxdev*100).toFixed(4), sum_abs_dev:+(sumabs*100).toFixed(3),
              readout_has_species_pct: ro.indexOf(D.dec[SP].pct_species.toFixed(0)+'%')>=0,
              readout_has_cell_pct: ro.indexOf(D.dec[SP].pct_cell.toFixed(0)+'%')>=0,
              chart_has_dashed: box.indexOf('stroke-dasharray="7,4"')>=0,
              chart_has_shaded: box.indexOf('#7b52ab')>=0});
  }
  return out;
}
"""


def main(path="diel_activity_viewer.html"):
    h = open(path).read()
    m = re.search(r'<script id="payload" type="application/json">(.*?)</script>', h, re.S)
    payload = m.group(1)
    js = h.split("<script>")[-1].split("</script>")[0]

    ctx = quickjs.Context()
    ctx.set_memory_limit(1 << 32)
    ctx.eval(SHIM)
    # the page reads its payload out of the DOM, exactly as the browser does: this is the step that
    # would have caught the blank page, because QuickJS's JSON.parse rejects bare NaN
    ctx.eval("document.getElementById('payload')._text = " + json.dumps(payload) + ";")
    ctx.eval(js)          # includes the trailing render() call
    ctx.eval(PROBE)
    combos = json.loads(ctx.eval("JSON.stringify(__probeCombos())"))
    clicks = json.loads(ctx.eval("JSON.stringify(__probeClicks())"))
    toggles = json.loads(ctx.eval("JSON.stringify(__probeToggles())"))
    dev = json.loads(ctx.eval("JSON.stringify(__probeDeviation())"))
    return combos, clicks, toggles, dev


if __name__ == "__main__":
    combos, clicks, toggles, dev = main(
        sys.argv[1] if len(sys.argv) > 1 else "diel_activity_viewer.html")
    json.dump({"combos": combos, "clicks": clicks, "toggles": toggles, "deviation": dev},
              open("js_exec_result.json", "w"), indent=1)
    r = combos["results"]
    print(f"COMBINATIONS: {len(r)} rendered (species x measure x map type),"
          f" {combos['empty']} with zero cells")
    print(f"  verdict-word mismatches vs beats_fix: {combos['bad_verdict']}")
    print(f"  support-marking failures: {combos['bad_mark']}")
    print(f"  cells drawn min/max: {min(x['cells'] for x in r)} / {max(x['cells'] for x in r)}")
    print(f"  curve panel drawn on every combo: {all(x['curve_drawn'] for x in r)}")
    print(f"  decomposition readout on every combo: {all(x['dec_readout'] for x in r)}")
    st = {}
    for x in r:
        st[x["state"]] = st.get(x["state"], 0) + 1
    print("  support states:", st)
    print(f"CLICKS: {len(clicks)} (2 map types x 2 curve modes x 5 species)")
    ok_sel = sum(1 for x in clicks if x.get("selected_on_map") == 1)
    print(f"  cell selected on map: {ok_sel}/{len(clicks)}")
    devm = [x for x in clicks if x["mode"] == "deviation"]
    absm = [x for x in clicks if x["mode"] == "absolute"]
    print(f"  deviation mode: key present {sum(1 for x in devm if x['has_devkey'])}/{len(devm)},"
          f" two shaded panels {sum(1 for x in devm if x['has_two_panels'])}/{len(devm)},"
          f" labelled {sum(1 for x in devm if x['labelled'])}/{len(devm)}")
    print(f"  absolute mode: key absent {sum(1 for x in absm if not x['has_devkey'])}/{len(absm)}")
    print(f"TOGGLES: {len(toggles['log'])} steps, tabs {toggles['ntabs']},"
          f" clear-selection {toggles['cleared']}, modal {toggles['modal']}")
    bad = [x for x in toggles["log"] if x["map_cells"] == 0 or not x["curve_svg"]]
    print(f"  steps leaving an empty map or no curve: {len(bad)}")
    print("DEVIATION CHART:")
    for x in dev:
        print(f"  {x['sp']}: largest departure {x['max_dev_pp']:.3f} pp/bin,"
              f" total {x['sum_abs_dev']:.2f} pp, reference line {x['chart_has_dashed']},"
              f" shading {x['chart_has_shaded']},"
              f" readout {x['readout_has_species_pct'] and x['readout_has_cell_pct']}")
