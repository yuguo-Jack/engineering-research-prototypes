# 工程研究原型 · 30 秒 HyperFrames 发布片

V2：**整整 30 秒，1920 × 1080、30 fps、140 BPM**。保留网站的深绿、暖白和
橙色，用短句、旋转线框、立体目录、图形切入与真实结果曲线介绍工程。

观看：[网页播放器](https://yuguo-jack.github.io/engineering-research-prototypes/#intro-video)。
成片与校验清单在 `assets/video/`。本目录为可编辑、可重新渲染的源工程。

## 分镜与数据口径

| 时间 | 画面内容 |
|---|---|
| 00:00 | 想法动起来：线框立方体与短句开场 |
| 00:01.7 | 工程研究原型：独立实现、真实结果、完整复现 |
| 00:03.4 | 问题 → 实现 → 验证 → 复现 |
| 00:06.0 | 一个成果一个目录；源码、环境、结果分层展示 |
| 00:08.6 | 通信雷达、机器人、嵌入式与流体计算的应用方向 |
| 00:11.1 | 圆柱绕流：真实采样图、单核 66 步、0.223 s |
| 00:15.4 | M6：32 MPI、真实 400 步历史、稳定窗确认 265.477 s |
| 00:20.6 | 真实网站入口与运行方法 |
| 00:24.0 | 邀请带问题、数据或设备条件来提需求 |
| 00:26.6 | 品牌、仓库地址与行动入口 |

主图在上方，中文字幕位于独立底栏，不覆盖图表。播放器可选逐字字幕默认关闭，
避免与画面底栏重复。通信雷达、机器人与嵌入式标注为面向方向；实际已收录案例为 CFD。

圆柱指标来自 [RESULTS.md](../../projects/cylinder-flow/RESULTS.md) 与正式运行记录：
9 次中位数原值 **0.223238375 s**，66 步，Cp RMSE 0.0003429041。
M6 指标来自 [RESULTS.md](../../projects/m6-wing/RESULTS.md)：1,485,672 单元、32 MPI，
保留完整 400 步。**265.4771408 s 是第 306 步确认 207—306 百步稳定窗的耗时**，
不表示相对参考的精度认证。历史平台耗时不等于观众机器上的保证值或同机加速比。

## 源工程

- `DESIGN.md`：配色、排版、运动、字幕安全区。
- `storyboard.json`：帧数固定的 10 个章节、简短旁白与底栏字幕。
- `timing.json`：实际旁白时长、节奏校正因子与章节位置。
- `build_composition.py` → `index.html`：真实曲线、CSS 3D 和确定性 GSAP 动画。
- `assets/master.m4a`：当前旁白与原创合成配乐母带；重新渲染不必再次生成语音。
- `build_audio.py`：联网合成旁白，按 140 BPM 合成鼓点、低音、和弦与转场声。
- `render.py`：渲染、网页封装、封面、可选字幕轨与 SHA256 清单。
- `fonts/erp-sans.woff2` / `prepare_font.py`：自带的 OFL 字体子集与生成工具。

## 重新渲染

只渲染需要 Python 3.10+、Node.js 22+、FFmpeg/ffprobe、Chromium/Edge，以及
锁定的 HyperFrames 0.9.19。动画使用固定版本 GSAP 3.15.0；首次编译需要网络。
音频和字体再生成使用 Python 3.12+ 与本目录 `requirements.txt`。

在仓库根目录执行：

```powershell
# 已有共享 HyperFrames 可以使用后面的 --cli，免于重复安装。
npm --prefix media/product-intro ci --cache D:/codex/tools/npm-cache
python media/product-intro/build_composition.py
cd media/product-intro
npm run check
cd ../..
python media/product-intro/render.py
python tools/catalog.py
python tools/check_repo.py --release
python tools/preview.py --port 8876
```

使用本机已有共享工具：

```powershell
python media/product-intro/render.py --cli D:/codex/tools/hyperframes/node_modules/hyperframes/bin/hyperframes.mjs --browser "C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe"
```

网站本地地址 `http://127.0.0.1:8876/#intro-video`。在本目录执行
`npm run dev -- --port 3017` 可打开 HyperFrames Studio，选择 `product-intro` 逐帧检查。
不同浏览器、字体和编码器可能产生不同的 MP4 字节；附带成片的身份见
`assets/video/manifest.json`。

## 修改文案与音频

```powershell
python -m pip install -r media/product-intro/requirements.txt
# 修改 storyboard.json，旁白缓存按文案、音色和语速的哈希自动失效。
python media/product-intro/build_audio.py
python media/product-intro/build_composition.py
# 增加当前字体子集以外的字符时：
python media/product-intro/prepare_font.py --source "<Noto Sans SC variable TTF>" --license "<对应 OFL.txt>"
python media/product-intro/render.py
```

旁白使用 edge-tts 的 `zh-CN-YunxiNeural`、语速 +20%，必要时进行有限的保调节奏校正，
超过 1.4 倍会拒绝并要求缩短文案。章节固定到 900 帧，总长始终 30 秒。
`--direct` 仅在当前进程绕过失效代理，不更改系统网络设置。
改变章节时，同步更新 `catalog-config.json` 的跳转点与视频版本标识。

## 灵感和素材来源

参考 [Linear 官方产品发布内容](https://linear.app/changelog/2024-09-19-introducing-linear-mobile)
的产品画面优先与简洁视觉表达，以及
[短片镜头拆解](https://cutmap.co/v/linear-introducing-linear-mobile) 提供的节拍组织思路。
本片图形、镜头、动画与音效独立制作，没有复制参考片的音乐或画面。

字体为重命名的 Noto Sans SC 子集，保留版权与 [SIL OFL 1.1](fonts/OFL.txt)；
完整字体来源 [Google Fonts](https://github.com/google/fonts/tree/main/ofl/notosanssc)。
不分发专有系统字体。HyperFrames 遵循 Apache-2.0；GSAP、edge-tts 与所连接语音服务
各自遵循自身许可和服务条款，不受本仓库材料声明替代。

## 临时内容

`.local/`、`.hyperframes/`、`node_modules/`、渲染帧、检查截图、WAV 和语音缓存不进 Git。
每批结束盘点并清理确认废弃的大文件。最终 MP4、字幕、封面、音频母带、动画源文件
和字体许可保留；不存在反复生成的大型历史归档。
