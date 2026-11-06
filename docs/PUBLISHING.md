# 新成果收录与发布

每次收录围绕一个已经做出的原型展开；领域标签不决定工具或编程语言。
通信雷达、机器人、嵌入式、控制和其他原型都使用同样的项目边界。

1. 从仓库根目录安装维护工具依赖：`python -m pip install -r requirements.txt`。
2. 用 `python tools/new_project.py <id> --title "成果名称" --domain signals` 创建独立草稿。
3. 填写 README、REPRODUCE、RESULTS、DATA、LICENSE 和 project.json；补齐实际代码、
   环境、输入清单、运行入口、关键图与小型原始结果。
4. 按声明范围实际构建或运行，保存带日期、环境和范围的验证记录。未运行完整工况时，
   不以单元测试、旧日志核对或绘图代替完整运行。
5. 检查材料来源、权利状态与 Git 内容。大型输入、原生临时文件、缓存和无效候选
   留在忽略目录或授权的外部附件中。填写审阅状态后再准备发布。
6. 运行下列检查，检查实际 diff，再提交。远端创建、push 和页面发布需有用户授权。

```powershell
python tools/catalog.py
python tools/catalog.py --check
python tools/check_repo.py --release
python -m unittest discover -s tests -v
git diff --check
git status --short
```

`catalog.py` 同步 README 目录、总览、项目详情及可浏览的 HTML 文档。
生成页面由源 Markdown 和项目元数据维护，不直接编辑生成的 HTML。

项目运行仅使用该项目自己的依赖；仓库根目录的 Markdown/元数据工具不参与数值求解。
不同项目的验证命令不同，按各自复现说明执行，不给所有新原型强加 CFD 的入口或环境。

现有结果保持原样。更新算法或工况时，在新输出目录产生新测量，记录差异和新版本，
不要用后续候选覆盖已经归档的正式结果。
