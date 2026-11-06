"""Build the 30-second, beat-cut HyperFrames launch film from real project data."""
from pathlib import Path
import argparse
import csv
import html
import json
import math

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def esc(value):
    return html.escape(str(value), quote=True)


def cube():
    return '<div class="cube-scaffold" data-layout-ignore aria-hidden="true"><div class="cube">' + ''.join(f'<div class="face face-{i}"><i></i><i></i><b>{"P" if i == 0 else "+"}</b></div>' for i in range(6)) + '</div><div class="cube-shadow"></div></div>'


def wave():
    paths = []
    for row in range(14):
        points = ' '.join(f'{x},{155 + row * 18 + math.sin(x / 74 + row * .19) * (22 + row * 3):.1f}' for x in range(0, 850, 12))
        paths.append(f'<polyline points="{points}" fill="none" stroke="currentColor" stroke-width="1.6"/>')
    return '<svg class="mesh-wave" viewBox="0 0 850 470" data-layout-ignore aria-hidden="true">' + ''.join(paths) + '</svg>'


def domains():
    graphics = [
        '<svg viewBox="0 0 320 260"><circle cx="160" cy="130" r="100"/><circle cx="160" cy="130" r="60"/><path d="M35 130H285M160 8V250"/><path class="radar-arm" d="M160 130L240 70L255 125Z" fill="currentColor" fill-opacity=".22"/><circle class="hot" cx="216" cy="95" r="8"/></svg>',
        '<svg viewBox="0 0 320 260"><path d="M50 232H280M117 230V204H175V230"/><g class="robot-arm"><path class="draw" d="M146 205L109 125L206 62L264 108" stroke-width="14"/><circle cx="146" cy="205" r="15"/><circle cx="109" cy="125" r="15"/><circle cx="206" cy="62" r="15"/><path class="hot" d="M263 108L287 105L299 128M263 108L273 135L295 137"/></g></svg>',
        '<svg viewBox="0 0 320 260"><rect x="75" y="55" width="170" height="150" rx="13"/><rect x="113" y="90" width="94" height="80" rx="7"/><path class="draw" d="M108 55V20M147 55V20M187 55V20M220 55V20M108 205V239M147 205V239M187 205V239M220 205V239M75 90H35M75 130H35M75 170H35M245 90H285M245 130H285M245 170H285"/></svg>',
        '<svg viewBox="0 0 320 260"><path class="draw" d="M39 202L142 58L291 112L150 204Z"/><path d="M142 58L150 204M181 73L188 184M219 85L225 163M255 100L259 135M61 173L274 120M84 142L233 106M110 104L190 89"/></svg>'
    ]
    labels = [('通信雷达', 'SIGNALS'), ('机器人', 'ROBOTICS'), ('嵌入式', 'EMBEDDED'), ('流体计算', 'CFD / 2 CASES')]
    return ''.join(f'<div class="domain part">{graphic}<h3>{label}</h3><p>{english}</p></div>' for graphic, (label, english) in zip(graphics, labels))


def m6_chart():
    with (ROOT / 'projects/m6-wing/results/raw/history.csv').open(encoding='utf-8-sig', newline='') as stream:
        values = [float(row['fy']) for row in csv.DictReader(stream)]
    x = lambda step: 82 + (step - 1) / 399 * 1000
    y = lambda value: 414 - (value - 1750) / (4700 - 1750) * 330
    points = ' '.join(f'{x(i + 1):.2f},{y(v):.2f}' for i, v in enumerate(values))
    ticks = ''.join(f'<text x="{x(step):.1f}" y="463" text-anchor="middle">{step}</text>' for step in (1, 207, 306, 400))
    return f'''<svg class="force-chart" viewBox="0 0 1140 490" role="img" aria-label="真实400步Fy历史与207至306百步稳定窗">
<text x="12" y="36">Fy / N</text><text x="1080" y="36" text-anchor="end">真实 400 步记录</text>
<path d="M82 105H1082M82 218H1082M82 331H1082M82 414H1082" class="chart-grid"/>
<rect class="window" x="{x(207):.1f}" y="74" width="{x(306)-x(207):.1f}" height="340" rx="7"/>
<g class="trace"><polyline points="{points}" fill="none" stroke="currentColor" stroke-width="5" stroke-linejoin="round"/></g>
<path class="marker" d="M{x(306):.1f} 74V414" stroke="#ca713d" stroke-width="3" stroke-dasharray="7 7"/>
<circle class="marker-dot" cx="{x(306):.1f}" cy="{y(values[305]):.1f}" r="9" fill="#ca713d"/>{ticks}</svg>'''


STYLE = r'''
@font-face{font-family:"ERP Sans";src:url("fonts/erp-sans.woff2") format("woff2");font-weight:100 900;font-style:normal;font-display:block}
:root{--paper:#fbfcfb;--ink:#172d30;--teal:#0b746b;--orange:#ca713d;--line:#dce5e2;--muted:#637477}
*{box-sizing:border-box;margin:0;padding:0}html,body,#root{width:1920px;height:1080px;overflow:hidden}body{font-family:"ERP Sans",sans-serif;-webkit-font-smoothing:antialiased}#root{position:relative;background:var(--ink)}
.scene{position:absolute;inset:0;width:1920px;height:1080px;background:var(--paper);color:var(--ink);overflow:hidden}.scene.dark{background:var(--ink);color:var(--paper)}
.scene-content{width:100%;height:100%;padding:54px 96px 120px;display:flex;flex-direction:column;gap:20px;position:relative}
.mast{display:flex;align-items:center;gap:30px;height:44px;flex:none;font-size:25px;letter-spacing:2px}.mast .wordmark{font-size:48px;font-weight:850;letter-spacing:-5px}.wordmark i{font-style:normal;color:var(--orange)}.mast .chapter{opacity:.8}.mast .count{margin-left:auto;font-size:25px;letter-spacing:1px}
.stage{flex:1;min-height:0;padding:20px 0;position:relative;display:flex;flex-direction:column;justify-content:center;gap:28px}.source-note{height:32px;flex:none;line-height:32px;font-size:24px;color:var(--muted);white-space:nowrap}.dark .source-note{color:var(--line)}
.caption-rail{height:76px;flex:none;border-top:1px solid var(--line);display:flex;align-items:center;gap:24px;font-size:34px;letter-spacing:.4px;line-height:1.3}.caption-rail:before{content:"";display:block;width:36px;height:4px;background:var(--orange);flex:none}.dark .caption-rail{border-color:#dce5e24d}
h1{font-size:86px;line-height:1.24;font-weight:850;letter-spacing:-3px}h2{font-size:70px;line-height:1.3;letter-spacing:-2px}h3{font-size:40px;line-height:1.4}p{font-size:32px;line-height:1.6}.accent{color:var(--teal)}.dark .accent{color:var(--orange)}.eyebrow{font-size:26px;letter-spacing:4px;font-weight:650}.small{font-size:29px}.mono{font-family:Consolas,monospace}
.back-grid{position:absolute;inset:-130px;background-image:linear-gradient(#63747720 1px,transparent 1px),linear-gradient(90deg,#63747720 1px,transparent 1px);background-size:96px 96px;transform:perspective(1100px) rotateX(32deg);transform-origin:center;opacity:.6}.dark .back-grid{opacity:.7}.halo{position:absolute;right:-160px;top:-300px;width:1100px;height:1100px;border:1px solid #63747733;border-radius:50%;box-shadow:0 0 0 60px #6374770c,0 0 0 150px #63747709}
.split{display:grid;grid-template-columns:1.1fr .9fr;align-items:center;gap:70px}.hook-copy{display:flex;flex-direction:column;gap:22px}.hook-copy h1{font-size:161px;line-height:1.17;letter-spacing:-5px}.hook-copy h1 span{display:block}.hook-copy .accent{font-size:137px}.hook-copy>p{font-size:30px;letter-spacing:5px}.brand-title h1{font-size:112px}.brand-title{display:flex;flex-direction:column;gap:30px}.brand-title p{font-size:32px}.brand-pills{display:flex;gap:18px}.brand-pills span{font-size:29px;padding:12px 22px;border:1px solid var(--line);border-radius:99px}
.cube-scaffold{width:620px;height:570px;perspective:1100px;display:grid;place-items:center;position:relative}.cube{width:300px;height:300px;display:grid;transform-style:preserve-3d;transform:rotateX(-23deg) rotateY(-32deg) rotateZ(-6deg);color:var(--teal)}.dark .cube{color:var(--line)}.face{grid-area:1/1;width:300px;height:300px;border:3px solid currentColor;background:#0b746b1a;display:grid;place-items:center;position:relative}.face i{position:absolute;inset:35px;border:1px solid currentColor;opacity:.25}.face i:nth-child(2){inset:75px}.face b{font-size:144px;font-weight:800;color:var(--orange)}.face-0{transform:translateZ(150px)}.face-1{transform:rotateY(180deg) translateZ(150px)}.face-2{transform:rotateY(90deg) translateZ(150px)}.face-3{transform:rotateY(-90deg) translateZ(150px)}.face-4{transform:rotateX(90deg) translateZ(150px)}.face-5{transform:rotateX(-90deg) translateZ(150px)}.cube-shadow{position:absolute;bottom:10px;width:430px;height:32px;background:#0b746b19;border-radius:50%;filter:blur(15px)}
.pipe-wrap{display:flex;flex-direction:column;gap:70px}.flow{display:grid;grid-template-columns:repeat(4,1fr);gap:65px;position:relative}.flow-item{display:flex;flex-direction:column;align-items:center;gap:25px}.pipe-node{height:132px;width:132px;border:2px solid #dce5e28a;border-radius:50%;display:grid;place-items:center;font-family:Consolas,monospace;font-size:53px;position:relative;background:var(--ink)}.pipe-node:after{content:"";position:absolute;inset:13px;border:1px solid #dce5e233;border-radius:50%}.flow-item h3{font-size:43px}.flow-item p{font-size:24px;letter-spacing:4px;color:var(--line)}.flow-line{position:absolute;inset:58px 8% auto;width:84%;height:20px;stroke:var(--orange);fill:none;stroke-width:4}
.architecture{display:grid;grid-template-columns:1fr 1fr;align-items:center;gap:50px}.arch-copy{display:flex;flex-direction:column;gap:29px}.arch-copy h1{font-size:85px;line-height:1.35}.directory{font-family:Consolas,monospace;font-size:36px;line-height:1.8;color:var(--teal)}.deck{height:540px;display:grid;place-items:center;perspective:1300px;transform-style:preserve-3d}.sheets{display:grid;transform-style:preserve-3d;transform:rotateX(12deg) rotateY(-16deg)}.sheet{grid-area:1/1;width:560px;height:335px;background:var(--paper);border:2px solid var(--teal);border-radius:18px;padding:32px;display:flex;flex-direction:column;gap:24px;box-shadow:14px 18px 0 #172d3008}.sheet:nth-child(1){transform:translate3d(-65px,62px,-80px) rotateZ(-9deg)}.sheet:nth-child(2){transform:translate3d(4px,9px,0px) rotateZ(3deg)}.sheet:nth-child(3){transform:translate3d(70px,-50px,100px) rotateZ(11deg)}.sheet h3{font-size:45px}.sheet p{font-size:30px;color:var(--teal)}.sheet .lines{height:80px;background:repeating-linear-gradient(0deg,transparent 0 23px,#dce5e2 24px 26px)}
.domain-wrap{display:flex;flex-direction:column;gap:42px}.domains{display:grid;grid-template-columns:repeat(4,1fr);gap:24px}.domain{display:flex;flex-direction:column;align-items:center;gap:16px}.domain svg{height:260px;width:320px;color:var(--line);fill:none;stroke:currentColor;stroke-width:3;stroke-linecap:round;stroke-linejoin:round}.domain svg .hot{stroke:var(--orange);fill:var(--orange)}.domain h3{font-size:39px}.domain p{font-size:22px;letter-spacing:2px;color:var(--line)}
.case-stage{display:grid;grid-template-columns:minmax(0,1fr) 450px;gap:56px;align-items:center}.case-main{display:flex;flex-direction:column;gap:20px;min-width:0}.case-head{display:flex;align-items:baseline;gap:24px}.case-head h1{font-size:82px}.case-head span{font-size:29px;letter-spacing:2px;color:var(--teal)}.dark .case-head span{color:var(--line)}.field-wrap{border:1px solid var(--line);border-radius:15px;padding:18px;background:var(--paper);overflow:hidden}.field-wrap img{display:block;width:100%;height:470px;object-fit:contain}.case-stats{display:flex;flex-direction:column;gap:30px}.stat{display:flex;flex-direction:column;gap:6px}.stat span{font-size:31px}.stat strong{font-family:"ERP Sans",sans-serif;font-size:134px;line-height:1.1;letter-spacing:-6px;color:var(--teal);font-weight:780;white-space:nowrap}.stat strong small{font-size:47px;letter-spacing:0;margin-left:10px;font-weight:550}.dark .stat strong{color:var(--paper)}.stat .compact{font-size:94px;letter-spacing:-3px}.stat p{font-size:27px}.case-stats .tiny-stat{font-size:27px;color:var(--muted)}.dark .tiny-stat{color:var(--line)}
.force-chart{width:100%;height:510px;color:var(--paper);display:block}.force-chart text{font-family:"ERP Sans",sans-serif;fill:var(--line);font-size:28px;stroke:none}.chart-grid{stroke:#dce5e226;stroke-width:1;fill:none}.force-chart .window{fill:var(--orange);fill-opacity:.16}.m6-stats .stat strong{font-size:123px}.m6-stats .stat strong.compact{font-size:91px}.window-copy{font-size:30px;color:var(--line)}
.repro{display:grid;grid-template-columns:510px minmax(0,1fr);gap:55px;align-items:center}.repro-copy{display:flex;flex-direction:column;gap:32px}.repro-copy h1{font-size:82px;line-height:1.35}.command{font-family:Consolas,monospace;color:var(--teal);font-size:31px;line-height:1.8}.browser{border:1px solid var(--line);border-radius:19px;overflow:hidden;background:var(--paper);transform:perspective(1500px) rotateY(-6deg);box-shadow:18px 23px 0 #172d3008}.browser-top{height:44px;border-bottom:1px solid var(--line);display:flex;align-items:center;padding:0 18px;font-size:22px;gap:14px}.browser-top i{height:10px;width:10px;background:var(--orange);border-radius:50%}.browser img{display:block;width:100%;height:510px;object-fit:cover;object-position:top}.browser-surface{overflow:hidden}.mesh-wave{width:760px;height:470px;color:var(--line)}
.invite-layout{display:grid;grid-template-columns:1.2fr .8fr;gap:40px;align-items:center}.invite-layout h1{font-size:117px;line-height:1.28}.invite-layout h1 span{display:block}.closing-layout{display:grid;grid-template-columns:minmax(0,1fr) 470px;gap:60px;align-items:center}.closing-copy{display:flex;flex-direction:column;gap:26px}.closing-copy h1{font-size:101px}.closing-copy .closing-brand{font-size:45px;letter-spacing:1px}.closing-link{font-family:Consolas,monospace;font-size:31px;line-height:1.5;color:var(--line);white-space:nowrap}.cta{display:flex;gap:18px}.cta span{font-size:34px;padding:16px 28px;border:1px solid #dce5e26a;border-radius:9px}.cta span:first-child{background:var(--paper);color:var(--ink)}.closing-logo{font-size:190px;font-weight:850;letter-spacing:-14px;text-align:center}.closing-orbit{border:1px solid #dce5e24a;border-radius:50%;width:460px;height:460px;display:grid;place-items:center;box-shadow:0 0 0 28px #dce5e208,0 0 0 58px #dce5e204}
'''


def content():
    return {
        'hook': ('dark', f'<div class="split"><div class="hook-copy"><p class="eyebrow part">MAKE IDEAS MOVE.</p><h1 class="part"><span>想法。</span><span class="accent">动起来。</span></h1></div><div class="part">{cube()}</div></div>', 'ENGINEERING RESEARCH PROTOTYPES'),
        'brand': ('', f'<div class="split"><div class="brand-title"><p class="eyebrow part">ENGINEERING RESEARCH PROTOTYPES</p><h1 class="part">工程研究原型<span class="accent">。</span></h1><div class="brand-pills part"><span>独立实现</span><span>真实结果</span><span>完整复现</span></div></div><div class="part">{cube()}</div></div>', 'SOURCE · RUN · VERIFY'),
        'pipeline': ('dark', '<div class="pipe-wrap"><h1 class="part">从问题，到<span class="accent">可复现。</span></h1><div class="flow"><svg class="flow-line" viewBox="0 0 1200 20" preserveAspectRatio="none" data-layout-ignore><path class="draw" d="M0 10H1200"/></svg>' + ''.join(f'<div class="flow-item part"><div class="pipe-node">{i:02}</div><h3>{a}</h3><p>{b}</p></div>' for i, (a, b) in enumerate([('问题','QUESTION'),('实现','BUILD'),('验证','VERIFY'),('复现','REPRODUCE')], 1)) + '</div></div>', '问题 → 独立实现 → 验证证据 → 复现入口'),
        'architecture': ('', '<div class="architecture"><div class="arch-copy"><p class="eyebrow part">INDEPENDENT BY DESIGN</p><h1 class="part">一个成果，<br><span class="accent">一个独立目录。</span></h1><p class="directory part">projects/<br>　cylinder-flow/<br>　m6-wing/</p></div><div class="deck part" data-layout-allow-overlap><div class="sheets" aria-label="独立目录包含源码环境及结果材料">' + ''.join(f'<div class="sheet"><h3>{a}</h3><p>{b}</p><div class="lines"></div></div>' for a, b in [('源码','src/ · 独立运行入口'),('环境','README · REPRODUCE'),('结果','results/ · 图表与日志')]) + '</div></div></div>', '依赖 · 输入 · 运行入口'),
        'scenarios': ('dark', '<div class="domain-wrap"><h1 class="part">下一项，<span class="accent">什么方向？</span></h1><div class="domains">' + domains() + '</div></div>', '预留方向 / 已有 CFD × 2'),
        'cylinder': ('', '<div class="case-stage"><div class="case-main"><div class="case-head part"><h1>圆柱绕流</h1><span>CYLINDER FLOW</span></div><div class="field-wrap part"><img src="assets/cylinder-field.png" width="1800" height="792" alt="真实圆柱压力与马赫数采样"></div></div><div class="case-stats"><div class="stat part"><span>单核收敛</span><strong>66<small>步</small></strong></div><div class="stat part"><span>九次正式中位求解时间</span><strong class="compact">0.223<small>s</small></strong></div><p class="tiny-stat part">Cp RMSE 0.0003429</p></div></div>', '历史 Linux 单核 · 9 次中位数'),
        'm6': ('dark', f'<div class="case-stage"><div class="case-main"><div class="case-head part"><h1>M6 机翼</h1><span>148.6 万单元</span></div><div class="part">{m6_chart()}</div></div><div class="case-stats m6-stats"><div class="stat part"><span>并行规模</span><strong>32<small>MPI</small></strong></div><div class="stat part"><span>百步稳定窗确认</span><strong class="compact">265.477<small>s</small></strong></div><p class="window-copy part">207 起点 → 306 确认</p></div></div>', '历史 Linux · 稳定窗 ≠ 精度认证'),
        'reproduce': ('', '<div class="repro"><div class="repro-copy"><p class="eyebrow part">READY TO REPRODUCE</p><h1 class="part">看成果，<br><span class="accent">也能照着跑。</span></h1><p class="command part">&gt; doctor<br>&gt; build / test<br>&gt; run / check-run</p></div><div class="browser part"><div class="browser-top"><i></i>工程研究原型 / 成果总览</div><div class="browser-surface" data-layout-allow-overflow><img src="assets/website-overview.png" width="1440" height="980" alt="真实成果网站"></div></div></div>', '流程示意 · 输入条件见 DATA.md'),
        'invite': ('dark', f'<div class="invite-layout"><div><p class="eyebrow part">YOUR NEXT ENGINEERING IDEA</p><h1 class="part"><span>下一个，</span><span class="accent">做什么？</span></h1></div><div class="part">{wave()}</div></div>', '问题 / 数据 / 设备'),
        'closing': ('dark', '<div class="closing-layout"><div class="closing-copy"><p class="closing-brand part">Engineering Research Prototypes</p><h1 class="part">一起把想法<span class="accent">做出来。</span></h1><div class="cta part"><span>提一个需求 ↗</span><span>浏览成果</span></div><p class="closing-link part">github.com/yuguo-Jack/<br>engineering-research-prototypes</p></div><div class="closing-orbit part"><div class="closing-logo wordmark">erp<i>.</i></div></div></div>', '工程研究原型 · 从问题出发，持续积累。')
    }


def motion(timing):
    data = json.dumps([{k: s[k] for k in ('id', 'start', 'duration')} for s in timing['scenes']])
    return r'''
window.__timelines=window.__timelines||{};
const tl=gsap.timeline({paused:true});
const scenes=SCENES;
const eases=["power4.out","expo.out","back.out(1.12)","circ.out"];
scenes.forEach((scene,i)=>{
 const selector="#scene-"+scene.id,s=scene.start;
 if(i){
  tl.set(selector,{opacity:1},s);
  if(i%3===1)tl.fromTo(selector,{clipPath:"polygon(0% 0%,0% 0%,0% 100%,0% 100%)"},{clipPath:"polygon(0% 0%,100% 0%,100% 100%,0% 100%)",duration:.22,ease:"power3.inOut"},s);
  else if(i%3===2)tl.fromTo(selector,{x:1920},{x:0,duration:.22,ease:"power4.out"},s);
  else tl.fromTo(selector,{clipPath:"inset(48% 0% 48% 0%)"},{clipPath:"inset(0% 0% 0% 0%)",duration:.2,ease:"expo.out"},s);
  tl.set("#scene-"+scenes[i-1].id,{opacity:0},s+.24);
 }
 const parts=document.querySelectorAll(selector+" .part");
 parts.forEach((element,j)=>{
  const options={opacity:0,duration:.3+(j%2)*.08,ease:eases[j%4]};
  if(j%3===0)options.y=24;else if(j%3===1){options.x=-48;options.scale=.97;}else{options.y=44;options.rotation=1.6;}
  if(element.classList.contains("closing-orbit"))delete options.scale;
  tl.from(element,options,s+.04+j*.055);
 });
 tl.from(selector+" .caption-rail",{y:10,opacity:0,duration:.2,ease:"power2.out"},s+.15);
 tl.from(selector+" .source-note",{opacity:0,duration:.2,ease:"sine.out"},s+.2);
 tl.from(selector+" .stage",{scale:1.035,duration:.65,ease:"power4.out"},s);
 tl.to(selector+" .back-grid",{x:45,y:-20,duration:scene.duration,ease:"none"},s);
 tl.to(selector+" .halo",{rotation:24,scale:1.08,duration:scene.duration,ease:"sine.inOut"},s);
 const wire=document.querySelector(selector+" .cube");
 if(wire)tl.to(wire,{rotationY:i===0?66:145,rotationX:i===0?-8:-38,rotationZ:7,duration:scene.duration,ease:"none"},s);
});
document.querySelectorAll(".draw").forEach(path=>{
 const parent=path.closest(".scene"),scene=scenes.find(s=>parent.id==="scene-"+s.id),length=path.getTotalLength();
 tl.fromTo(path,{strokeDasharray:length,strokeDashoffset:length},{strokeDashoffset:0,duration:Math.min(1.4,scene.duration-.45),ease:"power2.out"},scene.start+.35);
});
tl.from("#scene-architecture .sheet",{z:-160,y:110,rotationZ:-8,duration:.55,stagger:.12,ease:"back.out(1.12)"},scenes[3].start+.15);
tl.to("#scene-architecture .sheets",{rotationY:12,rotationX:5,duration:scenes[3].duration,ease:"sine.inOut"},scenes[3].start);
tl.to("#scene-pipeline .pipe-node",{scale:1.12,borderColor:"#ca713d",duration:.2,stagger:.42,yoyo:true,repeat:1,ease:"power2.inOut"},scenes[2].start+.28);
tl.to("#scene-scenarios .radar-arm",{rotation:240,svgOrigin:"160 130",duration:2.1,ease:"none"},scenes[4].start+.1);
tl.from("#scene-scenarios .robot-arm",{rotation:28,svgOrigin:"146 205",duration:1.15,ease:"power3.out"},scenes[4].start+.2);
tl.from("#scene-cylinder .field-wrap img",{scale:1.1,duration:3.8,ease:"power1.out"},scenes[5].start+.12);
tl.from("#scene-m6 .trace",{clipPath:"inset(0% 100% 0% 0%)",duration:3.65,ease:"none"},scenes[6].start+.35);
tl.from("#scene-m6 .window",{scaleX:0,transformOrigin:"left center",duration:.45,ease:"power2.out"},scenes[6].start+2.25);
tl.from("#scene-m6 .marker",{scaleY:0,transformOrigin:"center bottom",duration:.3,ease:"expo.out"},scenes[6].start+3.35);
tl.from("#scene-m6 .marker-dot",{scale:0,transformOrigin:"50% 50%",duration:.35,ease:"back.out(1.5)"},scenes[6].start+3.5);
tl.to("#scene-reproduce .browser",{rotationY:2,rotationX:0,duration:3.3,ease:"power2.out"},scenes[7].start);
tl.to("#scene-reproduce .browser img",{y:-95,duration:2.8,ease:"power2.inOut"},scenes[7].start+.4);
tl.from("#scene-invite .mesh-wave",{rotationY:45,scale:.7,duration:.8,ease:"expo.out"},scenes[8].start+.1);
tl.to("#scene-invite .mesh-wave",{rotationX:18,rotationZ:-7,duration:2.4,ease:"sine.inOut"},scenes[8].start+.1);
tl.from("#scene-closing .closing-orbit",{rotationY:65,scale:.75,duration:.65,ease:"back.out(1.08)"},scenes[9].start+.1);
tl.to("#scene-closing .closing-orbit",{rotationZ:9,duration:2.6,ease:"sine.inOut"},scenes[9].start+.5);
window.__timelines.main=tl;tl.seek(0);
'''.replace('SCENES', data)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--layout-only', action='store_true')
    args = parser.parse_args()
    timing = json.loads((HERE / 'timing.json').read_text(encoding='utf-8'))
    pieces = content()
    scenes = []
    for i, row in enumerate(timing['scenes']):
        theme, body, note = pieces[row['id']]
        scenes.append(f'''<section id="scene-{row['id']}" class="scene {theme}" style="z-index:{i+1};{'opacity:0;' if i else ''}" aria-label="{esc(row['chapter'])}"><div class="back-grid" data-layout-ignore></div><div class="halo" data-layout-ignore></div><div class="scene-content"><header class="mast" data-layout-allow-overlap><span class="wordmark">erp<i>.</i></span><span class="chapter">{esc(row['chapter'])}</span><span class="count">{i+1:02} / 10</span></header><div class="stage">{body}</div><p class="source-note" data-layout-allow-overlap>{esc(note)}</p><div class="caption-rail" data-layout-allow-overlap>{esc(row['caption'])}</div></div></section>''')
    script = 'window.__timelines={main:gsap.timeline({paused:true})};' if args.layout_only else motion(timing)
    page = f'''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=1920,height=1080"><title>工程研究原型 · 30秒</title><link rel="icon" href="data:,"><script src="https://cdn.jsdelivr.net/npm/gsap@3.15.0/dist/gsap.min.js"></script><style>{STYLE}</style></head><body><div id="root" data-composition-id="main" data-start="0" data-duration="30" data-width="1920" data-height="1080" data-fps="30">{''.join(scenes)}<audio id="soundtrack" class="clip" data-start="0" data-duration="30" data-track-index="5" src="assets/master.m4a"></audio></div><script>{script}</script></body></html>\n'''
    (HERE / 'index.html').write_text(page, encoding='utf-8', newline='\n')
    (HERE / 'meta.json').write_text(json.dumps({'title':'工程研究原型 · 30秒快节奏版','duration':30,'width':1920,'height':1080,'fps':30},ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'seconds':30,'scenes':len(scenes),'layout_only':args.layout_only}))


if __name__ == '__main__':
    main()
