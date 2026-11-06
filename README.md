# Engineering Research Ideas · 工程研究灵感

把灵感做成可运行、可核对、可复现的工程原型。

An independent collection of engineering prototypes, with source, reproducible
entry points, measured results, and explicit environment/data requirements.

**[打开成果总览](index.html)** · [组织规则](CONTRIBUTING.md) ·
[复现等级](docs/REPRODUCIBILITY.md) · [新增成果模板](_templates/project/README.md)

当前是本地审阅版本。每个 `projects/<id>/` 是一个独立成果，拥有自己的源码、
依赖、构建入口、结果和说明；复制单个目录即可进行该成果的工作。

## 已收录成果

<!-- PROJECTS:START -->
| 成果 | 范围 | 状态 | 入口 |
|---|---|---|---|
| [圆柱绕流 · 单核有限体积求解](projects/cylinder-flow/README.md) | 计算流体 | 已归档 | [复现](projects/cylinder-flow/REPRODUCE.md) · [结果](projects/cylinder-flow/RESULTS.md) |
| [M6 机翼 · MPI 并行流动求解](projects/m6-wing/README.md) | 计算流体 | 已归档 | [复现](projects/m6-wing/REPRODUCE.md) · [结果](projects/m6-wing/RESULTS.md) |
<!-- PROJECTS:END -->

时间、精度和平台信息以各项目 `RESULTS.md` 为准。历史正式测量、本次本地检查和
外部给定参考分别标注；给定参考与本机结果不自动构成同机性能比较。

## 从这里开始

1. 直接打开 `index.html` 查看总览。页面不依赖联网字体、前端框架或构建服务。
2. 阅读所选成果的 `README.md`，再按 `REPRODUCE.md` 检查环境、编译和运行。
3. 仅想核对已保存结果时，执行项目内的 `python reproduce.py verify-results`。

```powershell
python tools/check_repo.py
python tools/catalog.py --check
python tools/preview.py --port 8765
```

预览地址为 `http://127.0.0.1:8765`；退出时按 Ctrl+C。静态总览也可直接从文件打开。

```text
projects/<id>/       独立成果：代码、环境、复现方法、结果、来源
_templates/project/  新成果模板
site/                由项目元数据生成的展示页面及静态资源
tools/               目录生成、新成果创建和仓库检查工具
schemas/             项目元数据约定
docs/                复现、资料与发布约定
.github/workflows/   仓库检查和不依赖外部算例的数值测试
```

## 新增一个成果

```powershell
python tools/new_project.py my-prototype --title "我的工程原型" --domain robotics
# 填写项目说明、依赖、数据要求、结果和 project.json
python tools/catalog.py
python tools/check_repo.py
```

可使用 `cfd`、`signals`、`robotics`、`embedded`、`controls`、`other` 等领域标签。
新项目默认为草稿，验证通过后再更新状态。不要把未来方向伪装成已完成成果。

## 资料与授权

见 [资料与授权说明](docs/MATERIALS.md) 和各成果的 `LICENSE.md`。
原始竞赛输入与参考文件按其来源授权管理，通过项目的数据清单从本地授权归档导入，
不随 Git 仓库分发。此预览未为继承的源码新增开源授权。
