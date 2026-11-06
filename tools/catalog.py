"""Generate a deterministic, offline-capable overview from independent project metadata."""
from pathlib import Path
import argparse
import html
import json
import re

ROOT = Path(__file__).resolve().parents[1]
STATUSES = {"draft": "草稿", "prototype": "原型", "reproduced": "已复现", "archived": "已归档"}


def local_file(project, relative):
    if not isinstance(relative, str) or not relative or Path(relative).is_absolute():
        raise ValueError("Expected a project-relative file path")
    path = (project / relative).resolve()
    if not path.is_relative_to(project.resolve()) or not path.is_file():
        raise ValueError(f"Missing or escaping file: {relative}")
    return path


def load_projects(root=ROOT):
    config = json.loads((root / "catalog-config.json").read_text(encoding="utf-8"))
    projects = []
    for path in sorted((root / "projects").glob("*/project.json")):
        project = json.loads(path.read_text(encoding="utf-8"))
        for key in ("schema_version", "id", "title", "summary", "domain", "status", "version", "license", "readme", "reproduce", "results", "data"):
            if key not in project:
                raise ValueError(f"{path.parent.name}: missing {key}")
        if project["schema_version"] != 1 or not re.fullmatch(r"[a-z0-9]+(?:-[a-z0-9]+)*", project["id"]):
            raise ValueError("Unsupported schema or invalid project id")
        if project["id"] != path.parent.name or project["status"] not in STATUSES or project["domain"] not in config["domains"]:
            raise ValueError(f"{path.parent.name}: invalid id, status or domain")
        for key in ("readme", "reproduce", "results"):
            local_file(path.parent, project[key])
        for key in ("hero", "metrics", "validation"):
            if project.get(key):
                local_file(path.parent, project[key])
        if project["status"] in ("reproduced", "archived") and not project.get("validation"):
            raise ValueError("Verified/archived projects require a validation record")
        if project["data"].get("availability") not in ("not-required", "included", "external-permission-required"):
            raise ValueError("Invalid data availability")
        if project["data"].get("availability") == "external-permission-required":
            local_file(path.parent, project["data"]["manifest"])
        for figure in project.get("figures", []):
            local_file(path.parent, figure["path"])
        projects.append(project)
    if len({p["id"] for p in projects}) != len(projects):
        raise ValueError("Duplicate project id")
    return config, projects


def e(value):
    return html.escape(str(value), quote=True)


def header(prefix=""):
    return f'''<header class="topbar"><a class="brand" href="{prefix}index.html"><span class="brand-mark">eri<span>.</span></span><span>Engineering Research Ideas</span></a><nav aria-label="主导航"><a href="{prefix}index.html#projects">成果目录</a><a href="{prefix}index.html#reproduce">如何复现</a><a href="{prefix}CONTRIBUTING.md">组织规则</a></nav></header>'''


def shell(title, body, prefix=""):
    return f'''<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><meta name="description" content="工程原型、可复现代码与实际结果"><title>{e(title)}</title><link rel="stylesheet" href="{prefix}site/styles.css"></head>
<body><a class="skip-link" href="#main">跳到内容</a>{header(prefix)}<main id="main">{body}</main><footer><span>Engineering Research Ideas</span><span>源码 · 环境 · 结果 · 复现</span><a href="{prefix}docs/MATERIALS.md">资料与授权</a></footer><script src="{prefix}site/app.js" defer></script></body></html>\n'''


def metric_tiles(project):
    return ''.join(f'<div class="metric"><span class="metric-label">{e(m["label"])}</span><strong>{e(m["value"])}<small>{e(m.get("unit", ""))}</small></strong><span class="metric-note">{e(m.get("note", ""))}</span></div>' for m in project.get("highlights", []))


def card(project, config):
    slug = project["id"]
    image = f'<img src="projects/{slug}/{e(project["hero"])}" alt="{e(project["title"])}关键结果" width="700" height="420" loading="lazy">' if project.get("hero") else '<div class="text-cover">原型整理中</div>'
    search = ' '.join([project["title"], project["summary"], *project.get("tags", []), project["domain"]])
    tags = ''.join(f'<span>{e(t)}</span>' for t in project.get("tags", []))
    first = next(iter(project.get("highlights", [])), None)
    number = f'<div class="card-number"><strong>{e(first["value"])}</strong><span>{e(first.get("unit", ""))}</span><small>{e(first["label"])}</small></div>' if first else ''
    return f'''<article class="project-card" data-domain="{e(project["domain"])}" data-search="{e(search)}"><a class="card-visual" href="site/projects/{slug}.html">{image}<span class="status">{e(STATUSES[project["status"]])}</span></a><div class="card-body"><div class="card-domain">{e(config["domains"][project["domain"]])}</div><h3><a href="site/projects/{slug}.html">{e(project["title"])}</a></h3><p>{e(project["summary"])}</p>{number}<div class="tags">{tags}</div><div class="card-links"><a href="site/projects/{slug}.html">查看成果 <span aria-hidden="true">↗</span></a><a href="projects/{slug}/{e(project["reproduce"])}">复现说明</a></div></div></article>'''


def overview(config, projects):
    domains = [key for key in config["domains"] if any(p["domain"] == key for p in projects)]
    filters = '<button type="button" data-filter="all" aria-pressed="true">全部成果</button>' + ''.join(f'<button type="button" data-filter="{e(d)}" aria-pressed="false">{e(config["domains"][d])}</button>' for d in domains)
    cards = ''.join(card(p, config) for p in projects)
    directions = ''.join(f'<span>{e(x)}</span>' for x in config["future_directions"])
    return shell(config["title"] + " · Engineering Research Ideas", f'''
<section class="hero"><div><p class="eyebrow">ENGINEERING RESEARCH IDEAS</p><h1>让灵感成为<br><em>可复现的工程原型。</em></h1><p class="hero-copy">记录做出来的东西，也留下让别人重新做出来的方法。<br>从计算模型到运行程序，从关键结果到验证证据。</p><div class="hero-actions"><a class="button primary" href="#projects">浏览研究成果 <span aria-hidden="true">↓</span></a><a class="quiet-link" href="CONTRIBUTING.md">了解组织方式 ↗</a></div><div class="hero-counts"><div><strong>{len(projects):02}</strong><span>独立成果</span></div><div><strong>{sum(p["status"] in ("reproduced", "archived") for p in projects):02}</strong><span>附验证记录</span></div><div><strong>01</strong><span>统一组织规则</span></div></div></div><aside class="hero-aside"><span class="aside-label">从想法到可以复现</span><ol class="pipeline"><li><span>01</span><div><strong>一个明确的问题</strong><p>目标、输入、范围与约束</p></div></li><li><span>02</span><div><strong>一个独立的原型</strong><p>源码、依赖与运行入口</p></div></li><li><span>03</span><div><strong>一份可核对的结果</strong><p>指标、图表与实际证据</p></div></li><li><span>04</span><div><strong>一条完整的复现路径</strong><p>从环境到运行，再到验证</p></div></li></ol><div class="aside-footer"><span class="dot"></span>持续积累 · 当前为本地预览</div></aside></section>
<section id="projects" class="projects-section"><div class="section-heading"><div><p class="eyebrow">PROJECT COLLECTION</p><h2>已有成果</h2></div><p>每个项目独立运行，结果有据可查。</p></div><div class="collection-toolbar"><div class="filters" role="group" aria-label="按领域筛选">{filters}</div><label class="search"><span>搜索成果</span><input id="project-search" type="search" placeholder="名称、领域或技术关键词" autocomplete="off"></label></div><p id="result-count" class="result-count" aria-live="polite">共 {len(projects)} 项成果</p><div id="project-grid" class="project-grid">{cards}</div><div id="empty-state" class="empty-state" hidden><h3>暂时没有匹配成果</h3><p>换一个关键词，或返回全部成果。</p><button type="button" id="clear-search">清除筛选</button></div></section>
<section id="reproduce" class="reproduce-section"><div><p class="eyebrow">REPRODUCIBLE BY DESIGN</p><h2>能看到结果，<br>也能找到如何运行。</h2><p>每项成果都说明已验证到哪一步。历史测量、当前检查和外部参考各有清晰口径。</p><a class="quiet-link" href="docs/REPRODUCIBILITY.md">阅读复现约定 ↗</a></div><div class="reproduce-grid"><article><span>01 / ENVIRONMENT</span><h3>环境与依赖</h3><p>系统、工具链、版本和安装条件；可选依赖按需准备。</p></article><article><span>02 / EXECUTION</span><h3>构建与运行</h3><p>完整命令、明确输入和新的输出目录，单项成果可独立使用。</p></article><article><span>03 / EVIDENCE</span><h3>结果与证据</h3><p>原始数值、图表和验证记录，说明指标与比较边界。</p></article><article><span>04 / DATA</span><h3>输入与来源</h3><p>公开数据与外部资料分别声明，获取条件和校验路径明确。</p></article></div></section>
<section class="future-section"><div><p class="eyebrow">ROOM FOR THE NEXT IDEA</p><h2>为下一项成果留好位置。</h2><p>从模板开始，沿用同一套组织方式。以下是后续可收录的方向。</p></div><div class="directions">{directions}</div><a class="button secondary" href="_templates/project/README.md">查看项目模板 ↗</a></section>
''')


def detail(project, config):
    slug = project["id"]
    base = f'../../projects/{slug}/'
    figures = ''.join(f'<figure><a href="{base}{e(f["path"])}"><img src="{base}{e(f["path"])}" alt="{e(f["caption"])}" loading="lazy"></a><figcaption>{e(f["caption"])}</figcaption></figure>' for f in project.get("figures", []))
    commands = '\n'.join(project.get("commands", []))
    limits = ''.join(f'<li>{e(x)}</li>' for x in project.get("limits", []))
    data = "需要有授权的原始输入。可先执行源码测试和已保存结果核对，再按 DATA.md 准备完整工况。" if project["data"]["availability"] == "external-permission-required" else "数据条件以项目 DATA.md 为准。"
    return shell(project["title"], f'''<section class="detail-heading"><a class="back-link" href="../../index.html#projects">← 返回成果目录</a><p class="eyebrow">{e(config["domains"][project["domain"]])} / {e(STATUSES[project["status"]])}</p><h1>{e(project["title"])}</h1><p class="detail-summary">{e(project["summary"])}</p><div class="detail-meta"><span>版本 {e(project["version"])}</span><span>{e(' / '.join(project.get('language', [])))}</span><span>源码 {e(project.get('source_revision', '')[:12])}</span></div></section><section class="metrics-row">{metric_tiles(project)}</section><section class="detail-content"><div class="detail-gallery"><div class="section-heading"><h2>关键结果</h2><a href="{base}{e(project["results"])}">完整结果说明 ↗</a></div>{figures}</div><aside class="detail-guide"><h2>开始复现</h2><p>在该项目目录内执行：</p><pre><code>{e(commands)}</code></pre><button class="copy-button" type="button" data-copy="commands">复制命令</button><p class="copy-status" aria-live="polite"></p><a class="button primary" href="{base}{e(project["reproduce"])}">环境、构建与运行说明 ↗</a><div class="guide-block"><h3>输入条件</h3><p>{e(data)}</p><a href="{base}DATA.md">获取与校验方法 ↗</a></div><div class="guide-block"><h3>结果边界</h3><ul>{limits}</ul></div><div class="guide-block"><h3>源码与验证</h3><a href="{base}README.md">项目 README ↗</a><a href="{base}results/VALIDATION.json">本次验证记录 ↗</a><a href="{base}SOURCE_MANIFEST.json">计算源码身份 ↗</a><a href="{base}LICENSE.md">授权状态 ↗</a></div></aside></section>''', "../../")


def generated_files(root=ROOT):
    config, projects = load_projects(root)
    files = {"index.html": overview(config, projects), "site/catalog.json": json.dumps({"schema_version": 1, "projects": projects}, ensure_ascii=False, indent=2) + "\n"}
    for project in projects:
        files[f'site/projects/{project["id"]}.html'] = detail(project, config)
    readme = (root / "README.md").read_text(encoding="utf-8")
    table = "| 成果 | 范围 | 状态 | 入口 |\n|---|---|---|---|\n" + ''.join(f'| [{p["title"]}](projects/{p["id"]}/README.md) | {config["domains"][p["domain"]]} | {STATUSES[p["status"]]} | [复现](projects/{p["id"]}/REPRODUCE.md) · [结果](projects/{p["id"]}/RESULTS.md) |\n' for p in projects)
    replacement = "<!-- PROJECTS:START -->\n" + table + "<!-- PROJECTS:END -->"
    if readme.count("<!-- PROJECTS:START -->") != 1 or readme.count("<!-- PROJECTS:END -->") != 1:
        raise ValueError("README catalog markers missing or repeated")
    files["README.md"] = re.sub(r"<!-- PROJECTS:START -->.*?<!-- PROJECTS:END -->", lambda _: replacement, readme, flags=re.S)
    return files


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    stale = []
    files = generated_files()
    for name, content in files.items():
        path = ROOT / name
        if args.check:
            if not path.is_file() or path.read_text(encoding="utf-8") != content:
                stale.append(name)
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(content, encoding="utf-8", newline="\n")
    if stale:
        raise SystemExit("Stale generated files; run python tools/catalog.py: " + ", ".join(stale))
    print(json.dumps({"status": "current" if args.check else "generated", "files": len(files)}))


if __name__ == "__main__":
    main()
