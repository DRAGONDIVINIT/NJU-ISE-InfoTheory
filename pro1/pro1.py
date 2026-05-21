"""
Project 1: Lossless source coding for text (Huffman coding).

Design: variable-length prefix codes from symbol frequencies (Huffman algorithm).
Verify: compare average code length L with Shannon entropy H; report efficiency η = H/L.
"""

from __future__ import annotations
import argparse
import heapq
import json
import math
import struct
import sys
import zlib
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

MAGIC = b"HUF1"
SAMPLE_FILE = Path(__file__).resolve().parent / "sample.txt"


@dataclass(order=True)
class _HeapItem:
    freq: int
    index: int
    node: "HuffmanNode" = field(compare=False, default=None)


@dataclass
class HuffmanNode:
    freq: int
    char: Optional[str] = None
    left: Optional["HuffmanNode"] = None
    right: Optional["HuffmanNode"] = None


def shannon_entropy(freq: Counter) -> float:
    """H(X) in bits per symbol (natural log base 2)."""
    total = sum(freq.values())
    if total == 0:
        return 0.0
    h = 0.0
    for count in freq.values():
        if count > 0:
            p = count / total
            h -= p * math.log2(p)
    return h


def build_huffman_tree(freq: Counter) -> HuffmanNode:
    if not freq:
        raise ValueError("empty frequency table")
    if len(freq) == 1:
        (char,) = freq.keys()
        return HuffmanNode(freq=freq[char], char=char)

    heap: List[_HeapItem] = []
    tie = 0
    for char, count in freq.items():
        heapq.heappush(heap, _HeapItem(count, tie, HuffmanNode(freq=count, char=char)))
        tie += 1

    while len(heap) > 1:
        a = heapq.heappop(heap).node
        b = heapq.heappop(heap).node
        merged = HuffmanNode(freq=a.freq + b.freq, left=a, right=b)
        heapq.heappush(heap, _HeapItem(merged.freq, tie, merged))
        tie += 1

    return heap[0].node


def build_codes(root: HuffmanNode) -> Dict[str, str]:
    codes: Dict[str, str] = {}

    def walk(node: HuffmanNode, prefix: str) -> None:
        if node.char is not None:
            codes[node.char] = prefix or "0"
            return
        if node.left:
            walk(node.left, prefix + "0")
        if node.right:
            walk(node.right, prefix + "1")

    walk(root, "")
    return codes


def average_code_length(freq: Counter, codes: Dict[str, str]) -> float:
    total = sum(freq.values())
    if total == 0:
        return 0.0
    return sum((count / total) * len(codes[char]) for char, count in freq.items())


def bits_to_bytes(bit_string: str) -> Tuple[bytes, int]:
    """Pack bits; return (padded_bytes, valid_bit_count)."""
    if not bit_string:
        return b"", 0
    pad = (8 - len(bit_string) % 8) % 8
    padded = bit_string + "0" * pad
    out = bytearray()
    for i in range(0, len(padded), 8):
        out.append(int(padded[i : i + 8], 2))
    return bytes(out), len(bit_string)


def bytes_to_bits(data: bytes, bit_count: int) -> str:
    bits = "".join(f"{b:08b}" for b in data)
    return bits[:bit_count]


class HuffmanCoder:
    def __init__(self, freq: Optional[Counter] = None) -> None:
        self.freq: Counter = freq or Counter()
        self.codes: Dict[str, str] = {}
        self.root: Optional[HuffmanNode] = None

    def fit(self, text: str) -> None:
        self.freq = Counter(text)
        self.root = build_huffman_tree(self.freq)
        self.codes = build_codes(self.root)

    def encode(self, text: str) -> bytes:
        if not self.codes:
            self.fit(text)
        bit_string = "".join(self.codes[c] for c in text)
        payload, bit_count = bits_to_bytes(bit_string)
        header = {
            "bit_count": bit_count,
            "codes": self.codes,
        }
        header_bytes = json.dumps(header, ensure_ascii=False).encode("utf-8")
        return MAGIC + struct.pack(">I", len(header_bytes)) + header_bytes + payload

    def decode(self, data: bytes) -> str:
        if not data.startswith(MAGIC):
            raise ValueError("invalid Huffman bitstream header")
        offset = len(MAGIC)
        (header_len,) = struct.unpack_from(">I", data, offset)
        offset += 4
        header = json.loads(data[offset : offset + header_len].decode("utf-8"))
        offset += header_len
        payload = data[offset:]
        self.codes = header["codes"]
        bit_count = header["bit_count"]
        bits = bytes_to_bits(payload, bit_count)
        reverse = {code: char for char, code in self.codes.items()}
        symbols: List[str] = []
        buf = ""
        for b in bits:
            buf += b
            if buf in reverse:
                symbols.append(reverse[buf])
                buf = ""
        if buf:
            raise ValueError("truncated or corrupt bitstream")
        return "".join(symbols)


def utf8_byte_entropy(text: str) -> float:
    data = text.encode("utf-8")
    freq = Counter(data)
    return shannon_entropy(freq)


def verify_lossless(original: str, compressed: bytes) -> str:
    coder = HuffmanCoder()
    restored = coder.decode(compressed)
    if restored != original:
        raise AssertionError("lossless check failed: decoded text differs from original")
    return restored


def analyze(text: str, compressed: bytes) -> dict:
    freq = Counter(text)
    coder = HuffmanCoder(freq)
    coder.root = build_huffman_tree(freq)
    coder.codes = build_codes(coder.root)

    h_symbol = shannon_entropy(freq)
    h_byte = utf8_byte_entropy(text)
    L = average_code_length(freq, coder.codes)
    eta = h_symbol / L if L > 0 else 0.0

    orig_bytes = len(text.encode("utf-8"))
    zlib_bytes = len(zlib.compress(text.encode("utf-8"), level=9))

    return {
        "symbols": len(text),
        "unique_symbols": len(freq),
        "entropy_per_symbol_bits": h_symbol,
        "entropy_per_byte_bits": h_byte,
        "avg_code_length_bits": L,
        "efficiency_eta": eta,
        "redundancy_bits": L - h_symbol,
        "original_utf8_bytes": orig_bytes,
        "compressed_bytes": len(compressed),
        "compression_ratio": orig_bytes / len(compressed) if compressed else 0.0,
        "bits_per_symbol_actual": (len(compressed) * 8) / len(text) if text else 0.0,
        "zlib_bytes": zlib_bytes,
        "zlib_ratio": orig_bytes / zlib_bytes if zlib_bytes else 0.0,
    }


def print_report(text: str, compressed: bytes, path: Optional[Path] = None) -> None:
    stats = analyze(text, compressed)
    title = f"信源编码效率报告 — {path}" if path else "信源编码效率报告"
    print("=" * 60)
    print(title)
    print("=" * 60)
    print(f"符号数（字符）           : {stats['symbols']}")
    print(f"不同符号种类数           : {stats['unique_symbols']}")
    print(f"Shannon 熵 H（每符号）   : {stats['entropy_per_symbol_bits']:.4f} bit/符号")
    print(f"Shannon 熵 H（每字节）   : {stats['entropy_per_byte_bits']:.4f} bit/字节 (UTF-8)")
    print(f"Huffman 平均码长 L       : {stats['avg_code_length_bits']:.4f} bit/符号")
    print(f"冗余度 (L - H)           : {stats['redundancy_bits']:.4f} bit/符号")
    print(f"编码效率 η = H/L         : {stats['efficiency_eta']:.4f}  (≤ 1，接近 Kraft 界)")
    print("-" * 60)
    print(f"原文体积 (UTF-8)         : {stats['original_utf8_bytes']} 字节")
    print(f"压缩后体积 (Huffman)     : {stats['compressed_bytes']} 字节")
    print(f"压缩比（原/压）          : {stats['compression_ratio']:.4f}")
    print(f"实际比特/符号（含开销）  : {stats['bits_per_symbol_actual']:.4f}")
    print(f"zlib 参考                : {stats['zlib_bytes']} 字节 (压缩比 {stats['zlib_ratio']:.4f})")
    print("=" * 60)
    print("无损往返验证：通过")
    print("说明：η 接近 1 表示码长接近熵界；压缩文件含码本开销，长文本压缩比通常更好。")


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def load_default_sample() -> str:
    if not SAMPLE_FILE.is_file():
        raise FileNotFoundError(f"默认样例文件不存在: {SAMPLE_FILE}")
    return load_text(SAMPLE_FILE)


BASE_DIR = Path(__file__).resolve().parent


def _list_files(directory: Path, suffix: str) -> List[Path]:
    return sorted(
        p for p in directory.iterdir() if p.is_file() and p.suffix.lower() == suffix
    )


def _prompt_line(prompt: str) -> str:
    try:
        return input(prompt).strip()
    except (EOFError, KeyboardInterrupt):
        print()
        raise KeyboardInterrupt


def _pick_path(
    title: str,
    candidates: List[Path],
    default: Optional[Path] = None,
) -> Optional[Path]:
    """从编号列表或路径字符串选择文件；空输入返回 default。"""
    print(title)
    for i, path in enumerate(candidates, 1):
        print(f"  [{i}] {path.name}")
    if default:
        print(f"  [回车] 使用默认：{default.name}")
    print("  或直接输入/拖入文件路径")
    choice = _prompt_line("请选择: ")
    if not choice:
        return default
    if choice.isdigit():
        idx = int(choice)
        if 1 <= idx <= len(candidates):
            return candidates[idx - 1]
        print("编号无效。")
        return None
    path = Path(choice.strip('"').strip("'"))
    if path.is_file():
        return path
    print(f"找不到文件: {path}")
    return None


def run_encode(src_path: Path, out_path: Optional[Path] = None) -> int:
    text = load_text(src_path)
    coder = HuffmanCoder()
    coder.fit(text)
    compressed = coder.encode(text)
    verify_lossless(text, compressed)
    if out_path:
        out_path.write_bytes(compressed)
        print(f"\n已写入 {out_path}（{len(compressed)} 字节）")
    print()
    print_report(text, compressed, src_path)
    return 0


def run_decode(huf_path: Path, restored_path: Optional[Path] = None) -> int:
    data = huf_path.read_bytes()
    text = HuffmanCoder().decode(data)
    if restored_path:
        restored_path.write_text(text, encoding="utf-8")
        print(f"已还原并保存到 {restored_path}（{len(text)} 字符）")
    else:
        print(text)
    return 0


def interactive_loop() -> int:
    print("=" * 60)
    print("  文本 Huffman 无损信源编码")
    print(f"  工作目录: {BASE_DIR}")
    print("=" * 60)

    while True:
        print()
        print("  [1] 编码并分析（打印效率报告）")
        print("  [2] 编码并保存为 .huf")
        print("  [3] 解码 .huf 并保存还原文本")
        print("  [0] 退出")
        try:
            action = _prompt_line("请选择功能 (0-3): ")
        except KeyboardInterrupt:
            print("\n再见。")
            return 0

        if action in ("0", "q", "Q", "exit"):
            print("再见。")
            return 0
        if action not in ("1", "2", "3"):
            print("请输入 0、1、2 或 3。")
            continue

        try:
            if action in ("1", "2"):
                txt_files = _list_files(BASE_DIR, ".txt")
                src = _pick_path("\n可选 UTF-8 文本:", txt_files, default=SAMPLE_FILE)
                if not src:
                    continue
                out: Optional[Path] = None
                if action == "2":
                    default_out = src.with_suffix(".huf")
                    hint = _prompt_line(
                        f"输出 .huf 路径（回车默认 {default_out.name}）: "
                    )
                    out = Path(hint) if hint else default_out
                    if not out.is_absolute():
                        out = BASE_DIR / out
                run_encode(src, out)
            else:
                huf_files = _list_files(BASE_DIR, ".huf")
                if not huf_files:
                    print("\n当前目录下没有 .huf 文件，请先执行 [2] 编码保存。")
                    continue
                huf = _pick_path("\n可选 .huf 文件:", huf_files)
                if not huf:
                    continue
                default_restored = huf.with_name(huf.stem + "_restored.txt")
                hint = _prompt_line(
                    f"还原文本保存路径（回车默认 {default_restored.name}）: "
                )
                restored = Path(hint) if hint else default_restored
                if not restored.is_absolute():
                    restored = BASE_DIR / restored
                run_decode(huf, restored)
        except KeyboardInterrupt:
            print("\n已取消当前操作。")
        except (FileNotFoundError, ValueError, AssertionError) as exc:
            print(f"\n错误: {exc}")

        _prompt_line("\n按回车继续...")


def main(argv: Optional[Iterable[str]] = None) -> int:
    raw = list(argv) if argv is not None else sys.argv[1:]
    if not raw:
        return interactive_loop()

    parser = argparse.ArgumentParser(description="文本 Huffman 无损信源编码")
    parser.add_argument(
        "input",
        nargs="?",
        type=Path,
        help="待压缩的 UTF-8 文本文件（默认 sample.txt）",
    )
    parser.add_argument("-o", "--output", type=Path, help="将压缩结果写入 .huf")
    parser.add_argument("--decode", type=Path, help="解码 .huf 并输出到 stdout")
    args = parser.parse_args(raw)

    if args.decode:
        if not args.decode.is_file():
            print(f"错误：找不到压缩文件 {args.decode}", file=sys.stderr)
            return 1
        return run_decode(args.decode)

    if args.input and args.input.is_file():
        src_path = args.input
    else:
        src_path = SAMPLE_FILE
        if args.input:
            print(f"警告：未找到 {args.input}，改用 {SAMPLE_FILE}。", file=sys.stderr)

    return run_encode(src_path, args.output)


if __name__ == "__main__":
    raise SystemExit(main())
