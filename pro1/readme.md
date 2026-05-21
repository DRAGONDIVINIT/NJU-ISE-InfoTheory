# Project 1：文本无损信源编码（Huffman）

**课程**：南京大学 ISE · 信息论基础（2026 春）  
**任务**：设计一种面向文本的无损信源编码，编程实现，并验证编码效率（平均码长与 Shannon 熵的关系、无损还原）。

本目录为仓库中的 **Project 1** 实现，主程序为 `pro1.py`（Python 标准库，无第三方依赖）。

---

## 目录结构

```
pro1/
├── readme.md      # 本文件
├── pro1.py        # Huffman 编解码、效率分析、命令行入口
└── sample.txt     # 默认信源样例（无参数或找不到输入文件时使用）
```

使用 `-o` 编码后会额外生成 `.huf` 压缩文件（路径由参数指定，不在仓库内）。

---

## 方案设计

| 项目 | 实现说明 |
|------|----------|
| 信源模型 | 将文本视为**字符**符号序列，用 `Counter` 统计经验频率，按无记忆信源处理 |
| 编码算法 | 最小堆（`heapq`）合并结点构建 Huffman 树；同频率用 `index` 打破平局，保证结果稳定 |
| 码字生成 | 自根向下：左子树追加 `0`、右子树追加 `1`；单符号时码字为 `"0"` |
| 比特打包 | 有效比特流末尾补 `0` 至整字节；头中记录有效位数 `bit_count` |
| 无损验证 | 编码后立即 `decode`，断言还原文本与原文一致 |
| 效率验证 | 计算 \(H\)、\(L\)、冗余 \(L-H\)、效率 \(\eta=H/L\)（理论上有 \(H \le L < H+1\)） |
| 参考对比 | 同文本 UTF-8 字节的 `zlib.compress(..., level=9)` 体积（仅参考，非实验必达） |

### 压缩文件格式（`.huf`）

由 `HuffmanCoder.encode` 写出，布局如下：

```
[4 B] 魔数 MAGIC = "HUF1"
[4 B] 大端 uint32：JSON 头长度
[变长] UTF-8 JSON：{ "bit_count": int, "codes": { "字符": "01串", ... } }
[变长] 载荷字节（比特串经 bits_to_bytes 打包）
```

解码时读取头与载荷，按码本逐比特匹配前缀码还原字符。码本以 JSON 存入文件头，短文本时**头开销较大**，压缩比可能小于 1，但不影响用 \(\eta\) 评价码本是否接近熵界。

---

## 代码结构（`pro1.py`）

| 符号 | 作用 |
|------|------|
| `shannon_entropy` | 由频率表计算 \(H(X)\)（bit/符号，以 2 为底） |
| `build_huffman_tree` / `build_codes` | 建树与生成字符→码字映射 |
| `average_code_length` | 平均码长 \(L=\sum p(x)\,l(x)\) |
| `bits_to_bytes` / `bytes_to_bits` | 比特流与字节互转（含尾部填充） |
| `HuffmanCoder` | `fit` 统计频率并建树；`encode` / `decode` 编解码 |
| `utf8_byte_entropy` | 按 UTF-8 **字节**统计的熵（报告参考项） |
| `verify_lossless` | 编码结果往返解码，与原文比对 |
| `analyze` / `print_report` | 汇总指标并打印中文报告 |
| `load_text` / `load_default_sample` | 读取 UTF-8 文本；默认路径为同目录 `sample.txt` |
| `main` | 命令行：`input`、`-o/--output`、`--decode` |

常量：`MAGIC = b"HUF1"`，`SAMPLE_FILE` 指向 `sample.txt`。

---

## 运行环境

- **Python**：3.8+（建议 3.10+）
- **依赖**：仅标准库（`heapq`、`json`、`struct`、`zlib`、`argparse` 等），**无需** `pip install`
- **编码**：读写文本均为 **UTF-8**；Windows 终端若中文乱码，可先执行 `chcp 65001` 或设置控制台为 UTF-8

---

## 如何运行

进入本目录后执行：

```powershell
cd "E:\NJU-ISE-Info Theory\pro1"
```

### 命令行参数

| 用法 | 说明 |
|------|------|
| `python pro1.py` | 读取同目录 `sample.txt`，编码并打印效率报告（**不写**文件） |
| `python pro1.py <文件>` | 对指定 UTF-8 文本编码、打印报告（**不写** `.huf`） |
| `python pro1.py <文件> -o out.huf` | 编码、打印报告，并将压缩结果写入 `out.huf` |
| `python pro1.py --decode out.huf` | 解码已生成的 `out.huf`，原文输出到 **stdout** |

说明：

- `<文件>` 路径不存在时：向 stderr 打印警告，并**回退为 `sample.txt`**。
- 未加 `-o` 时，压缩数据仅在内存中用于报告与无损验证；**解码前必须先**用 `-o` 写出 `.huf`。
- `--decode` 时若文件不存在，程序退出码为 `1` 并提示先编码。

### 示例

```powershell
# 默认样例 + 报告
python pro1.py

# 指定信源
python pro1.py sample.txt

# 自备 UTF-8 文本
python pro1.py test.txt

# 写出压缩文件（解码依赖此步骤）
python pro1.py test.txt -o test.huf

# 解码到控制台，可重定向保存
python pro1.py --decode test.huf | Out-File -Encoding utf8 restored.txt
```

---

## 输出报告说明

运行编码流程后，终端输出中文报告（`print_report`），主要字段如下：

| 输出项 | 含义 |
|--------|------|
| 符号数（字符） | 字符总数 `len(text)` |
| 不同符号种类数 | 频率表中不同字符数 |
| Shannon 熵 H（每符号） | 按**字符**统计的 \(H\)（bit/符号） |
| Shannon 熵 H（每字节） | 按 UTF-8 字节统计的熵 |
| Huffman 平均码长 L | 由码长加权的理论平均码长 |
| 冗余度 (L - H) | \(L - H\) |
| 编码效率 η = H/L | 越接近 1 越接近熵界（Kraft 意义下的最优性） |
| 原文体积 (UTF-8) | 原文 UTF-8 字节数 |
| 压缩后体积 (Huffman) | 内存中压缩 `bytes` 的长度（含 JSON 码本头） |
| 压缩比（原/压） | 原文字节数 / 压缩字节数 |
| 实际比特/符号（含开销） | `len(compressed)×8 / 符号数` |
| zlib 参考 | 标准库 zlib 压缩体积与压缩比 |
| 无损往返验证：通过 | `verify_lossless` 通过 |

**说明**：\(\eta \approx 1\) 表示**码字设计**接近最优；**压缩比**受 JSON 码本头影响，短文本常常小于 1 属正常现象。使用更长、重复度更高的文本可观察压缩比上升。

---

## 实验报告可写要点

1. 简述 Huffman 构造步骤与前缀码性质（Kraft 不等式、无歧义解码）。
2. 给出一次运行的 \(H\)、\(L\)、\(L-H\)、\(\eta\) 及「无损往返验证：通过」截图或表格。
3. 对比 `sample.txt` 与自备长文本：\(\eta\) 与压缩比的变化，解释码本开销。
4. 可选：说明 `zlib` 在字节级体积更小的原因（信源模型为字节、DEFLATE 字典等）。

---

## 默认样例 `sample.txt`

无有效输入文件时，程序读取与 `pro1.py` 同目录的 **`sample.txt`**（UTF-8）。内容为信息论相关中英文段落，整体重复 3 次，便于观察频率与 \(\eta\)。可直接编辑该文件，或复制为其他文件名做对比实验。

---

## 常见问题（Windows）

**现象**：执行 `python` 即报错  
`UnicodeDecodeError: 'gbk' codec can't decode ...`（`init_import_site`）。

**原因**：用户目录 `site-packages` 中 `.pth` 含 UTF-8 中文路径，与系统 GBK locale 冲突。

**处理**：临时禁用用户 site-packages 后再运行：

```powershell
$env:PYTHONNOUSERSITE = "1"
python pro1.py
```

或将问题 `.pth` 重命名为 `.pth.bak`。请使用 `python pro1.py`，勿依赖 Unix shebang。
