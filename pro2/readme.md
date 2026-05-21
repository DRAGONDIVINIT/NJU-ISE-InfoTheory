# Project 2：图像有损信源编码 + 信道编码

**课程**：南京大学 ISE · 信息论基础（2026 春）  
**任务**：对图像集进行有损信源编码，经 BSC/BEC 信道（可仿真）传输，经信道编码后恢复图像，并从**准确度、算法复杂度、PSNR** 三方面综合分析结果。

本目录为仓库中的 **Project 2** 实现，主程序为 `pro2.py`。**默认以交互菜单运行**，无需记忆命令行参数。测试图像需事先放入 `images/` 目录。

---

## 目录结构

```
pro2/
├── readme.md      # 本文件
├── pro2.py        # 有损编解码、信道仿真、指标统计、交互菜单
├── pro2.bat       # Windows 双击启动（进入交互菜单）
├── images/        # 图像集（.png / .jpg / .bmp / .gif，需自行准备）
└── output/        # 实验 JSON、还原 PNG（运行后生成）
```

运行完整实验或单张还原后，`output/` 中会出现 `results.json` 及 `*_q<质量>_restored.png` 等文件。

---

## 快速开始（推荐）

### 方式一：双击 `pro2.bat`（Windows）

自动切换到 UTF-8 并启动交互菜单，按提示用数字选择即可。

### 方式二：只运行主程序

```powershell
cd "E:\NJU-ISE-Info Theory\pro2"
pip install numpy Pillow
python pro2.py
```

**首次使用前**请将至少一张测试图放入 `images/`。

### 交互菜单

```
  [1] 运行完整实验
  [2] 单张图像编解码
  [3] 查看最近一次实验 JSON
  [0] 退出
```

| 功能 | 操作说明 |
|------|----------|
| **1 完整实验** | 对 `images/` 中**全部**图像，在 **IDEAL / BSC / BEC** 三种信道、质量 **40 / 75 / 90** 下批量仿真；逐条打印中文指标，输出汇总表，并保存 `output/results.json`。默认 BSC 交叉概率 \(p=0.02\)，BEC 擦除概率 \(\varepsilon=0.03\)，随机种子 `42`。 |
| **2 单张还原** | 列出 `images/` 内图像，输入编号选择；再输入质量因子 1–100（**回车默认 75**）。仅经有损编解码（不经信道误码），保存 `output/<原名>_q<质量>_restored.png`。 |
| **3 查看 JSON** | 显示 `output/results.json` 前 2000 字符；需先执行过 [1]。 |

每项结束后按回车返回菜单。`images/` 为空时会提示将图像放入该目录。

---

## 方案设计

| 项目 | 实现说明 |
|------|----------|
| 有损信源 | 8×8 DCT + JPEG 亮度量化表（质量因子 \(Q\)）+ Zigzag + DC-DPCM + AC 游程 + zlib 压缩 JSON |
| 信道编码 | **Hamming (7,4)** 纠错 + **(5,1) 重复码**（每比特重复 5 次，BSC 多数表决 / BEC 擦除表决） |
| 信道模型 | **IDEAL**（无误码）、**BSC**（以 \(p\) 独立翻转）、**BEC**（以 \(\varepsilon\) 独立擦除为 `E`） |
| 准确度 | 原图与恢复图像像素差 \(\le 1\) 的占比 |
| 复杂度 | 理论 \(O(N)\) 信源（\(N\) 为像素数）+ \(O(B)\) 信道码（\(B\) 为码长）；报告给出各阶段毫秒耗时 |
| 质量评价 | **PSNR**（dB）、像素准确度、码流 **BER**、Hamming 纠正次数 |

### 端到端流程

```
灰度图 → DCT/量化/游程（有损）→ zlib 字节流 → 比特流
      → Hamming(7,4) → 五倍重复 → IDEAL / BSC / BEC
      → 表决 / 擦除填充译码 → Hamming 译码 → 字节流
      → 逆 DCT/反量化 → 恢复图像
```

**IDEAL** 信道不引入误码，恢复图失真仅来自有损信源。**BSC/BEC** 在相同有损码流上叠加信道错误，检验级联信道码能否把比特流恢复回来。

### 信道参数（菜单 [1] 默认）

| 信道 | 参数 | 默认值 |
|------|------|--------|
| BSC | 交叉概率 \(p\) | `0.02` |
| BEC | 擦除概率 \(\varepsilon\) | `0.03` |
| 公共 | 随机种子 | `42` |

### 码流封装（内部 `IMG2`）

```
[4 B] 魔数 "IMG2"
[8 B] 大端 uint32 ×2：信源载荷长度、原始比特长度
[变长] Hamming 编码后的比特流（打包为字节）
```

有损载荷本体为 zlib 压缩的块系数 JSON（含宽高、质量、各 8×8 块 DPCM/RLE 数据）。

### 读图约定

- 统一转为**灰度**（`L`）  
- 长边缩至不超过 **256** 像素（保持比例，LANCZOS 缩放）  
- 支持扩展名：`.png`、`.jpg`、`.jpeg`、`.bmp`、`.gif`

---

## 代码结构（`pro2.py`）

| 符号 | 作用 |
|------|------|
| `encode_grayscale` / `decode_grayscale` | 有损信源编解码（DCT、量化、RLE） |
| `hamming74_encode_nibble` / `hamming74_decode_codeword` | Hamming (7,4) 编解码 |
| `repetition_encode` / `repetition_decode_*` | (5,1) 重复码（硬判决 / 擦除） |
| `bsc_channel` / `bec_channel` | BSC、BEC 仿真 |
| `hamming_decode_erasure_bits` | BEC 下带擦除的 Hamming 译码 |
| `run_single_image` | 单图端到端流水线，返回 `PipelineResult` |
| `run_experiment` / `print_summary` | 批量实验与汇总表 |
| `save_restored` | 单张有损编解码并保存 PNG |
| `interactive_loop` | 交互菜单主循环 |
| `main` | 无参数→交互；有参数→命令行 `--run` |

常量：`MAGIC = b"IMG2"`，`REP_FACTOR = 5`，目录 `IMAGES_DIR`、`OUTPUT_DIR`。

---

## 运行环境

- **Python**：3.8+
- **依赖**：`numpy`、`Pillow`（仅此两项）

安装依赖：

```powershell
pip install numpy Pillow
```

- **编码**：终端 UTF-8；Windows 建议用 `pro2.bat` 或 `chcp 65001` 避免中文乱码

---

## 命令行模式（可选）

需要脚本化时传入参数（**有任意参数时不再进入菜单**）：

```powershell
python pro2.py --run --quality 75 --p-bsc 0.02 --eps-bec 0.03 --seed 42
```

| 参数 | 说明 |
|------|------|
| `--run` | 对 `images/` 全部图像在 IDEAL/BSC/BEC 下实验，写出 `output/results.json` 与各图还原 PNG |
| `--quality` | 有损编码质量 1–100（默认 `75`）；命令行模式仅测**单一**质量档 |
| `--p-bsc` | BSC 交叉概率（默认 `0.02`） |
| `--eps-bec` | BEC 擦除概率（默认 `0.05`，与菜单 [1] 的 `0.03` 不同，可按需修改） |
| `--seed` | 随机种子（默认 `42`） |

---

## 输出报告说明

### 单条结果（每张图、每个信道、每个质量）

| 输出项 | 含义 |
|--------|------|
| 图像 | 文件名（不含扩展名） |
| 信道 | `IDEAL` / `BSC` / `BEC` |
| 质量因子 | 有损编码 \(Q\) |
| PSNR | 峰值信噪比（dB），越大越好 |
| 像素准确度 | \(\|\Delta\| \le 1\) 的像素占比（有损后常 \< 100%） |
| BER（码流） | 重复译码后相对发送 Hamming 码流的比特误码率 |
| 信源字节数 | zlib 压缩后的有损载荷大小 |
| 编码比特数 | 经重复码后的信道比特数 |
| 信道错误数 | BSC 翻转数或 BEC 擦除数 |
| Hamming 纠正 | 单比特纠错触发次数 |
| 耗时 (ms) | 编码、信道、译码及合计毫秒数 |

### 汇总表

菜单 [1] 或 `--run` 结束后打印**汇总报告**，含每张配置的 PSNR、准确度、BER、总耗时及全体平均值。

**解读提示**：若 BSC/BEC 下 PSNR、准确度与 IDEAL 接近，说明在当前 \(p,\varepsilon\) 下级联码足以保护码流；PSNR 主要由有损信源（\(Q\)）决定，IDEAL 与有信道但纠错成功时数值应一致。

---

## 图像集 `images/`

请自行准备测试图并放入 `images/`（如课程常用的 Lenna、Peppers，或自备照片）。仓库中可含示例图，也可仅保留目录说明后自行添加。

交互 [1] 会遍历该目录下全部支持格式；[2] 通过编号选择单张。长边超过 256 像素时会自动缩小以便实验可复现、耗时可控。

---

## 实验报告可写要点

1. 简述 DCT+量化有损编码与率失真关系（\(Q\) 与 PSNR）。  
2. 给出 IDEAL / BSC / BEC 下同一图像的 PSNR、准确度、BER 对比表（可直接引用 `results.json`）。  
3. 说明 BSC、BEC 模型假设及 \(p\)、\(\varepsilon\) 选取；解释重复码 + Hamming 级联的作用。  
4. 结合耗时与 \(O(N)\)、\(O(B)\) 讨论算法复杂度；可对比提高 \(Q\) 或增大 \(p\) 时的变化趋势。
