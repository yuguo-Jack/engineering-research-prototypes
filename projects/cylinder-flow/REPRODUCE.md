# 环境、构建、运行与验证

所有命令从本项目目录执行。最低工具为 Python 3.10+、CMake 3.20+、C++17 编译器；
本地包装检查使用 Python 3.12.10、MSVC 19.44、Ninja。不需要 MPI。编译只需要标准 C++ 库。

## 1. 检查与准备

```powershell
python reproduce.py doctor
```

Windows 使用已有 Visual Studio 2022 Build Tools 的 C++ 桌面工具链与 Windows SDK；
如未进入 Developer PowerShell，入口会尝试通过已安装的 vswhere 初始化进程内编译环境。
Linux 需要 C++ 编译器、CMake 和构建工具。MPI 项目还需安装匹配的运行时及开发包。
`单核 CPU、C++17、无需 MPI`。

Ubuntu/Debian 可安装如下依赖，再确认 CMake 不低于 3.20；其他发行版使用对应包管理器。

```bash
sudo apt-get update
sudo apt-get install -y python3 build-essential cmake ninja-build
# 若系统只提供 python3，后续 python 命令均换为 python3
python3 --version
cmake --version
```

Windows 从 [Microsoft C++ 安装说明](https://learn.microsoft.com/en-us/cpp/build/vscpp-step-0-installation?view=msvc-170)
安装或补齐“使用 C++ 的桌面开发”、x64 MSVC 工具和 Windows SDK。
从 [CMake 官方入口](https://cmake.org/download/) 取得 CMake，
从 [Ninja 官方入口](https://ninja-build.org/) 取得 Ninja，将其可执行目录加入 PATH。
本次验证使用 VS 2022/MSVC 19.44；其他版本需自行重新验证。已有工具满足要求时直接复用。
Python 使用已有 3.10+ 安装并确认 `python --version`；源码项目无需额外 pip 包。

`requirements.txt` 说明最小 Python 依赖；核心流程只使用标准库。
重绘或处理原始数据时再安装可选环境：

```powershell
python -m pip install -r requirements-analysis.txt
```

可选依赖固定到本次包装机版本，未证明所有旧 Python 或平台组合兼容。

## 2. 从源码构建与测试

```powershell
python reproduce.py build --jobs 2
python reproduce.py test
```

默认输出到 `build/`，保留 configure.log、compile.log、ctest.log。默认不启用 LTO、
fast-math 或主机专属指令。GCC/Clang 在确认目标后可用 `build --native` 单独实验；
目标变更后必须重建。详细开关见 CMakeLists.txt 和 profiles/。
从 `--native` 切回默认 build 会显式清除该入口缓存的主机指令选项；跨系统或更换编译器时，
使用新的 build 目录（或移走旧 build）后再构建。

## 3. 核对已保存结果

```powershell
python reproduce.py verify-results
```

先核验归档结果 SHA256，再从原始数值记录重新计算正式样本中位数或稳定窗里程碑。
这验证已有数据与说明一致，不会重新执行完整工况。

## 4. 准备原始初场并运行

```powershell
python reproduce.py import-data --archive "PATH_TO_AUTHORIZED_SUBMISSION.zip"
python reproduce.py run --input data/external/cylinder.bin --output runs/my-first-run --dry-run
python reproduce.py run --input data/external/cylinder.bin --output runs/my-first-run
python reproduce.py check-run --output runs/my-first-run
```

导入先验证原归档及输入 SHA256；输入或输出冲突时拒绝覆盖。
运行固定使用 profiles/benchmark.json，保存 stdout.log、stderr.log、run-record.json 和
程序产生的数值结果。完整工况从原始初场开始，不用参考终场初始化。

资源：一个物理核。正式输入约 0.22 MB；一个新运行目录通常较小。

## 5. 结果与绘图

依据 RESULTS.md 比较迭代、目标量、精度和时间；新平台的时间记录为新测量。
可用 `python reproduce.py plot` 从保存的原始数值绘制示例汇总图，输出到
`data/generated/summary.png`，保留归档图片及其 SHA256。

`check-run` 读取新输出而不覆盖它：要求原初场起算的完整 66 步历史、残差满足 1e-9，
并把最终 Cd 与归档原型比较。数值容差为 `max(1e-5 × |参考值|, 1e-7)`；这是本仓库的
回归检查容差，不是赛题官方精度阈值。不检查耗时，也不替代完整压力场误差评估。
检查失败时返回非零状态；保留该次输出定位差异，不要改写归档数据使其通过。

## 历史环境与常见问题

历史正式环境：Linux x86-64 / Hygon 7581B、GCC 12.2、CMake 3.26.5，一个物理求解核。
本次可移植构建不继承原二进制身份和计时。

历史正式环境：GCC 12.2、CMake 3.26.5，Linux x86-64 / Hygon 7581B，单物理核。

- 编译失败：先看 build/configure.log 或 build/compile.log，核实编译器和 SDK。
- 找不到 MPI：确认开发头文件和库，与 launcher 配套；在 doctor 中确认路径。
- 缺少数据：仅编译/测试无需输入；完整工况按 DATA.md 准备，不能随意替换文件。
- 输出目录已存在：换一个新 runs/ 子目录，避免覆盖上次日志。
- 新平台结果或时间不同：先核输入、版本、配置与采样定义，再评估浮点及资源差异。

macOS、其他 CPU/编译器和新的机器人/硬件目标没有由这些记录自动验证。
