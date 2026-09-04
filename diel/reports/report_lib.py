
import base64, json, html as _h, pandas as pd, numpy as np

def _img(vid, host):
    """Embed a figure as base64 so the file works offline and travels anywhere."""
    p = host.artifact_path(vid)
    with open(p, "rb") as f:
        return "data:image/png;base64," + base64.b64encode(f.read()).decode()

def _table(df, cols=None, fmt=None, highlight=None):
    d = df if cols is None else df[cols]
    out = ["<table><thead><tr>"]
    for c in d.columns:
        out.append(f"<th>{_h.escape(str(c))}</th>")
    out.append("</tr></thead><tbody>")
    # iterate the FULL frame so `highlight` can read columns not being displayed
    for (_, r), (_, rf) in zip(d.iterrows(), df.iterrows()):
        cls = ' class="hl"' if highlight is not None and highlight(rf) else ""
        out.append(f"<tr{cls}>")
        for c in d.columns:
            v = r[c]
            if fmt and c in fmt:
                s = fmt[c](v)
            elif isinstance(v, (int, np.integer)):
                s = f"{v:,}"
            elif isinstance(v, (float, np.floating)):
                s = "&mdash;" if pd.isna(v) else f"{v:.3f}"
            else:
                s = _h.escape(str(v))
            out.append(f"<td>{s}</td>")
        out.append("</tr>")
    out.append("</tbody></table>")
    return "".join(out)

CSS = """
:root{--ink:#1e2a30;--muted:#5b6b73;--line:#d8dee2;--accent:#1f6f8b;--warn:#b06010;--bad:#c0392b;--good:#2e7d4f;}
*{box-sizing:border-box}
body{margin:0;font:16px/1.65 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;
     color:#1e2a30;background:#fbfbfa;}
header{background:#1f4a5a;color:#fff;padding:34px 40px 28px;}
header h1{margin:0 0 6px;font-size:26px;font-weight:650;letter-spacing:-.2px}
header p{margin:0;opacity:.9;font-size:15px}
nav{position:sticky;top:0;background:#fff;border-bottom:1px solid #d8dee2;padding:10px 40px;z-index:20;
    display:flex;gap:18px;flex-wrap:wrap;font-size:13.5px}
nav a{color:#1f6f8b;text-decoration:none}
nav a:hover{text-decoration:underline}
main{max-width:1080px;margin:0 auto;padding:0 40px 80px}
section{padding:34px 0 10px;border-bottom:1px solid #e6ebee}
h2{font-size:22px;margin:0 0 4px;color:#1f4a5a}
h2 .num{color:#9aa8b0;font-weight:400;margin-right:10px}
h3{font-size:17px;margin:28px 0 6px;color:#28414c}
.lead{font-size:17px;color:#3d4d55;margin:6px 0 18px}
p{margin:12px 0}
figure{margin:22px 0 8px;background:#fff;border:1px solid #d8dee2;border-radius:6px;padding:14px}
figure img{width:100%;height:auto;display:block;border-radius:3px;cursor:zoom-in}
figcaption{font-size:14.5px;color:#4a5a63;margin-top:12px;padding-top:11px;border-top:1px solid #eef2f4}
figcaption b{color:#1f4a5a}
.plain{background:#f2f7f9;border-left:3px solid #1f6f8b;padding:12px 16px;margin:14px 0;font-size:15.5px}
.plain b{color:#1f4a5a}
.callout{border-left:3px solid #b06010;background:#fdf6ee;padding:12px 16px;margin:16px 0;font-size:15.5px}
.callout.bad{border-color:#c0392b;background:#fdf0ee}
.callout.good{border-color:#2e7d4f;background:#eef7f1}
.callout h4{margin:0 0 6px;font-size:15px;color:#1f4a5a}
table{border-collapse:collapse;width:100%;margin:16px 0;font-size:14px;background:#fff}
th,td{border:1px solid #dde3e7;padding:6px 9px;text-align:right}
th{background:#eef2f4;font-weight:600;text-align:right;color:#28414c}
td:first-child,th:first-child{text-align:left}
tr.hl td{background:#eef7f1;font-weight:600}
.note{font-size:14px;color:#5b6b73;font-style:italic}
ul,ol{margin:12px 0;padding-left:26px}
li{margin:6px 0}
.big{font-size:30px;font-weight:650;color:#1f6f8b;line-height:1.1}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:14px;margin:18px 0}
.card{background:#fff;border:1px solid #d8dee2;border-radius:6px;padding:14px 16px}
.card .lbl{font-size:13px;color:#5b6b73;margin-top:4px}
#lb{display:none;position:fixed;inset:0;background:rgba(12,20,24,.93);z-index:100;cursor:zoom-out;
    align-items:center;justify-content:center;padding:20px}
#lb.on{display:flex}
#lb img{max-width:100%;max-height:100%;object-fit:contain}
footer{padding:26px 40px;color:#6b7a82;font-size:13.5px;background:#f2f5f6;border-top:1px solid #d8dee2}
"""

JS = """
document.querySelectorAll('figure img').forEach(function(im){
  im.addEventListener('click', function(){
    var lb=document.getElementById('lb');
    document.getElementById('lbimg').src=im.src; lb.classList.add('on');
  });
});
document.getElementById('lb').addEventListener('click', function(){ this.classList.remove('on'); });
document.addEventListener('keydown', function(e){
  if(e.key==='Escape'){ document.getElementById('lb').classList.remove('on'); }
});
"""
