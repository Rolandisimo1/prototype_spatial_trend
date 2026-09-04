
import base64, html as _h, pandas as pd, numpy as np

def _img64(vid, host):
    with open(host.artifact_path(vid), "rb") as f:
        return "data:image/png;base64," + base64.b64encode(f.read()).decode()

def _tab(df, cols=None, hl=None, fmt=None):
    d = df if cols is None else df[cols]
    o = ["<table><thead><tr>"] + [f"<th>{_h.escape(str(c))}</th>" for c in d.columns] + ["</tr></thead><tbody>"]
    for (_, r), (_, rf) in zip(d.iterrows(), df.iterrows()):
        cls = ' class="hl"' if hl is not None and hl(rf) else ""
        o.append(f"<tr{cls}>")
        for c in d.columns:
            v = r[c]
            if fmt and c in fmt: s = fmt[c](v)
            elif isinstance(v, (bool, np.bool_)): s = "yes" if v else "no"
            elif isinstance(v, (int, np.integer)): s = f"{v:,}"
            elif isinstance(v, (float, np.floating)): s = "&ndash;" if pd.isna(v) else f"{v:.3f}"
            else: s = _h.escape(str(v))
            o.append(f"<td>{s}</td>")
        o.append("</tr>")
    return "".join(o + ["</tbody></table>"])

MEMO_CSS = """
*{box-sizing:border-box}
body{margin:0;font:16px/1.62 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;
     color:#1e2a30;background:#fbfbfa}
header{background:#1f4a5a;color:#fff;padding:32px 44px 26px}
header h1{margin:0 0 6px;font-size:25px;font-weight:650}
header p{margin:0;opacity:.92;font-size:14.5px;max-width:820px}
nav{position:sticky;top:0;background:#fff;border-bottom:1px solid #d8dee2;padding:9px 44px;z-index:20;
    display:flex;gap:16px;flex-wrap:wrap;font-size:13px}
nav a{color:#1f6f8b;text-decoration:none}
main{max-width:1080px;margin:0 auto;padding:0 44px 90px}
section{padding:32px 0 8px;border-bottom:1px solid #e6ebee}
h2{font-size:21px;margin:0 0 6px;color:#1f4a5a}
h2 .num{color:#9aa8b0;font-weight:400;margin-right:9px}
h3{font-size:16.5px;margin:26px 0 6px;color:#28414c}
h4{font-size:14.5px;margin:18px 0 4px;color:#28414c}
.lead{font-size:16.5px;color:#3d4d55;margin:6px 0 16px}
figure{margin:20px 0 8px;background:#fff;border:1px solid #d8dee2;border-radius:6px;padding:13px}
figure img{width:100%;height:auto;display:block;border-radius:3px;cursor:zoom-in}
figcaption{font-size:14px;color:#4a5a63;margin-top:11px;padding-top:10px;border-top:1px solid #eef2f4}
table{border-collapse:collapse;width:100%;margin:14px 0;font-size:13.4px;background:#fff}
th,td{border:1px solid #dde3e7;padding:5px 8px;text-align:right}
th{background:#eef2f4;font-weight:600;color:#28414c}
td:first-child,th:first-child{text-align:left}
tr.hl td{background:#eef7f1;font-weight:600}
.box{border-left:3px solid #1f6f8b;background:#f2f7f9;padding:11px 15px;margin:14px 0;font-size:15px}
.box.warn{border-color:#b06010;background:#fdf6ee}
.box.bad{border-color:#c0392b;background:#fdf0ee}
.box.good{border-color:#2e7d4f;background:#eef7f1}
.box h4{margin:0 0 5px;font-size:14.5px;color:#1f4a5a}
.note{font-size:13.4px;color:#5b6b73;font-style:italic}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(185px,1fr));gap:12px;margin:16px 0}
.card{background:#fff;border:1px solid #d8dee2;border-radius:6px;padding:12px 14px}
.card .big{font-size:26px;font-weight:650;color:#1f6f8b;line-height:1.1}
.card .lbl{font-size:12.6px;color:#5b6b73;margin-top:3px}
code{background:#eef2f4;padding:1px 4px;border-radius:3px;font-size:13px}
ul,ol{margin:11px 0;padding-left:25px}
li{margin:5px 0}
.verdict{display:inline-block;padding:1px 7px;border-radius:9px;font-size:12.4px;font-weight:600}
.v-yes{background:#dff0e5;color:#1e6b40}
.v-no{background:#f6dedb;color:#96291c}
.v-part{background:#fbeedd;color:#8a4f10}
#lb{display:none;position:fixed;inset:0;background:rgba(12,20,24,.94);z-index:100;cursor:zoom-out;
    align-items:center;justify-content:center;padding:18px}
#lb.on{display:flex}
#lb img{max-width:100%;max-height:100%;object-fit:contain}
footer{padding:24px 44px;color:#6b7a82;font-size:13px;background:#f2f5f6;border-top:1px solid #d8dee2}
"""
MEMO_JS = """
document.querySelectorAll('figure img').forEach(function(im){
  im.addEventListener('click', function(){
    document.getElementById('lbimg').src=im.src;
    document.getElementById('lb').classList.add('on');});});
document.getElementById('lb').addEventListener('click',function(){this.classList.remove('on');});
document.addEventListener('keydown',function(e){
  if(e.key==='Escape'){document.getElementById('lb').classList.remove('on');}});
"""
