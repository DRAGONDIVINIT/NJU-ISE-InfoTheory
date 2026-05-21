"""
Project 2: Lossy image source coding + channel coding over BSC and BEC.

Pipeline: image -> DCT/quantize (lossy) -> bitstream -> Hamming (7,4) -> channel
          -> decode -> inverse DCT -> recovered image.
Metrics: PSNR, pixel accuracy, BER, wall-clock time (complexity proxy).
"""

from __future__ import annotations

import argparse
import json
import math
import random
import struct
import sys
import time
import zlib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union

import numpy as np
from PIL import Image

MAGIC = b"IMG2"
BASE_DIR = Path(__file__).resolve().parent
IMAGES_DIR = BASE_DIR / "images"
OUTPUT_DIR = BASE_DIR / "output"

# 信道标识（内部小写，报告/菜单统一为大写缩写，与 BSC、BEC 格式一致）
CHANNEL_IDEAL = "ideal"
CHANNEL_BSC = "bsc"
CHANNEL_BEC = "bec"
DEFAULT_CHANNELS = (CHANNEL_IDEAL, CHANNEL_BSC, CHANNEL_BEC)

# Standard JPEG luminance quantization (quality scales this table).
JPEG_LUMA = np.array(
    [
        [16, 11, 10, 16, 24, 40, 51, 61],
        [12, 12, 14, 19, 26, 58, 60, 55],
        [14, 13, 16, 24, 40, 57, 69, 56],
        [14, 17, 22, 29, 51, 87, 80, 62],
        [18, 22, 37, 56, 68, 109, 103, 77],
        [24, 35, 55, 64, 81, 104, 113, 92],
        [49, 64, 78, 87, 103, 121, 120, 101],
        [72, 92, 95, 98, 112, 100, 103, 99],
    ],
    dtype=np.float64,
)

_DCT8 = None
_IDCT8 = None


def _dct8_matrix() -> Tuple[np.ndarray, np.ndarray]:
    global _DCT8, _IDCT8
    if _DCT8 is not None:
        return _DCT8, _IDCT8
    n = 8
    c = np.zeros((n, n), dtype=np.float64)
    for k in range(n):
        for i in range(n):
            if k == 0:
                c[k, i] = 1.0 / math.sqrt(n)
            else:
                c[k, i] = math.sqrt(2.0 / n) * math.cos(
                    math.pi * k * (2 * i + 1) / (2 * n)
                )
    _DCT8, _IDCT8 = c, c.T
    return _DCT8, _IDCT8


ZIGZAG = [
    (0, 0), (0, 1), (1, 0), (2, 0), (1, 1), (0, 2), (0, 3), (1, 2),
    (2, 1), (3, 0), (4, 0), (3, 1), (2, 2), (1, 3), (0, 4), (0, 5),
    (1, 4), (2, 3), (3, 2), (4, 1), (5, 0), (6, 0), (5, 1), (4, 2),
    (3, 3), (2, 4), (1, 5), (0, 6), (0, 7), (1, 6), (2, 5), (3, 4),
    (4, 3), (5, 2), (6, 1), (7, 0), (7, 1), (6, 2), (5, 3), (4, 4),
    (3, 5), (2, 6), (1, 7), (2, 7), (3, 6), (4, 5), (5, 4), (6, 3),
    (7, 2), (7, 3), (6, 4), (5, 5), (4, 6), (3, 7), (4, 7), (5, 6),
    (6, 5), (7, 4), (7, 5), (6, 6), (7, 7),
]

# ---------------------------------------------------------------------------
# DCT / quantization (lossy source coding)
# ---------------------------------------------------------------------------

def scale_quant_table(quality: int) -> np.ndarray:
    q = max(1, min(100, quality))
    if q < 50:
        scale = 5000 / q
    else:
        scale = 200 - 2 * q
    table = np.floor((JPEG_LUMA * scale + 50) / 100)
    return np.maximum(table, 1)


def dct2(block: np.ndarray) -> np.ndarray:
    c, _ = _dct8_matrix()
    return c @ block @ c.T


def idct2(block: np.ndarray) -> np.ndarray:
    _, ct = _dct8_matrix()
    return ct @ block @ ct.T


def block_to_zigzag(qblock: np.ndarray) -> List[int]:
    return [int(qblock[r, c]) for r, c in ZIGZAG]


def zigzag_to_block(coeffs: Sequence[int]) -> np.ndarray:
    block = np.zeros((8, 8), dtype=np.int32)
    for val, (r, c) in zip(coeffs, ZIGZAG):
        block[r, c] = val
    return block


def rle_encode(coeffs: Sequence[int]) -> List[Tuple[int, int]]:
    """Run-length on AC coefficients (63 values after DC)."""
    ac = list(coeffs[1:])
    runs: List[Tuple[int, int]] = []
    zero_run = 0
    for v in ac:
        if v == 0:
            zero_run += 1
            if zero_run == 16:
                runs.append((15, 0))
                zero_run = 0
        else:
            runs.append((zero_run, int(v)))
            zero_run = 0
    if zero_run > 0:
        runs.append((zero_run, 0))
    runs.append((0, 0))  # EOB
    return runs


def rle_decode(runs: Sequence[Tuple[int, int]], ac_len: int = 63) -> List[int]:
    ac: List[int] = []
    for zeros, val in runs:
        if zeros == 0 and val == 0:
            break
        ac.extend([0] * zeros)
        if val != 0:
            ac.append(val)
    ac.extend([0] * max(0, ac_len - len(ac)))
    return ac[:ac_len]


def encode_grayscale(arr: np.ndarray, quality: int) -> bytes:
    h, w = arr.shape
    pad_h = (8 - h % 8) % 8
    pad_w = (8 - w % 8) % 8
    padded = np.pad(arr.astype(np.float64), ((0, pad_h), (0, pad_w)), mode="edge")
    qtable = scale_quant_table(quality)
    blocks: List[dict] = []
    prev_dc = 0
    ph, pw = padded.shape
    for by in range(0, ph, 8):
        for bx in range(0, pw, 8):
            block = padded[by : by + 8, bx : bx + 8] - 128.0
            dct_b = dct2(block)
            qblock = np.round(dct_b / qtable).astype(np.int32)
            zig = block_to_zigzag(qblock)
            dc = int(zig[0])
            dpcm = dc - prev_dc
            prev_dc = dc
            runs = rle_encode(zig)
            blocks.append({"dpcm": dpcm, "runs": runs})
    payload = {"h": h, "w": w, "quality": quality, "blocks": blocks}
    return zlib.compress(json.dumps(payload).encode("utf-8"), level=9)


def decode_grayscale(data: bytes) -> np.ndarray:
    payload = json.loads(zlib.decompress(data).decode("utf-8"))
    h, w, quality = payload["h"], payload["w"], payload["quality"]
    qtable = scale_quant_table(quality)
    ph = h + (8 - h % 8) % 8
    pw = w + (8 - w % 8) % 8
    out = np.zeros((ph, pw), dtype=np.float64)
    prev_dc = 0
    idx = 0
    for by in range(0, ph, 8):
        for bx in range(0, pw, 8):
            blk = payload["blocks"][idx]
            idx += 1
            dc = prev_dc + int(blk["dpcm"])
            prev_dc = dc
            ac = rle_decode([tuple(x) for x in blk["runs"]])
            zig = [dc] + ac
            qblock = zigzag_to_block(zig) * qtable
            block = idct2(qblock.astype(np.float64)) + 128.0
            out[by : by + 8, bx : bx + 8] = block
    return np.clip(out[:h, :w], 0, 255).astype(np.uint8)


def bytes_to_bits(data: bytes) -> str:
    return "".join(f"{b:08b}" for b in data)


def bits_to_bytes(bits: str) -> bytes:
    pad = (8 - len(bits) % 8) % 8
    padded = bits + "0" * pad
    return bytes(int(padded[i : i + 8], 2) for i in range(0, len(padded), 8))


# ---------------------------------------------------------------------------
# Hamming (7, 4) channel coding
# Layout (1-indexed): [p1, p2, d1, p3, d2, d3, d4]
# ---------------------------------------------------------------------------


def _parity(*bits: int) -> int:
    return sum(bits) % 2


def hamming74_encode_nibble(d: Sequence[int]) -> List[int]:
    d0, d1, d2, d3 = (int(x) for x in d)
    p1 = _parity(d0, d1, d3)
    p2 = _parity(d0, d2, d3)
    p3 = _parity(d1, d2, d3)
    return [p1, p2, d0, p3, d1, d2, d3]


def hamming74_decode_codeword(r: List[int]) -> Tuple[List[int], int]:
    """Return 4 data bits [d0..d3] and correction count."""
    p1, p2, d0, p3, d1, d2, d3 = (int(x) for x in r)
    s1 = _parity(p1, d0, d1, d3)
    s2 = _parity(p2, d0, d2, d3)
    s3 = _parity(p3, d1, d2, d3)
    corrected = 0
    if s1 or s2 or s3:
        err_pos = s1 + 2 * s2 + 4 * s3  # 1..7
        if 1 <= err_pos <= 7:
            r[err_pos - 1] ^= 1
            corrected = 1
            p1, p2, d0, p3, d1, d2, d3 = r
    return [d0, d1, d2, d3], corrected


def hamming_encode_bits(bits: str) -> str:
    out: List[str] = []
    for i in range(0, len(bits), 4):
        chunk = bits[i : i + 4].ljust(4, "0")
        cw = hamming74_encode_nibble(int(b) for b in chunk)
        out.extend(str(b) for b in cw)
    return "".join(out)


REP_FACTOR = 5  # repetition before Hamming; majority vote tolerates 2 errors per 5 bits at p≈0.02


def repetition_encode(bits: str, n: int = REP_FACTOR) -> str:
    return "".join(b * n for b in bits)


def repetition_decode_hard(bits: str, n: int = REP_FACTOR) -> str:
    out: List[str] = []
    for i in range(0, len(bits), n):
        chunk = bits[i : i + n]
        ones = chunk.count("1")
        zeros = chunk.count("0")
        out.append("1" if ones > zeros else "0")
    return "".join(out)


def repetition_decode_erasure(bits: str, n: int = REP_FACTOR) -> str:
    out: List[str] = []
    for i in range(0, len(bits), n):
        chunk = bits[i : i + n]
        known = [c for c in chunk if c in ("0", "1")]
        if not known:
            out.append("0")
            continue
        ones = known.count("1")
        zeros = known.count("0")
        out.append("1" if ones > zeros else "0")
    return "".join(out)


def hamming_decode_bits(bits: str) -> Tuple[str, int]:
    corrected = 0
    out: List[str] = []
    for i in range(0, len(bits), 7):
        chunk = bits[i : i + 7].ljust(7, "0")
        data, c = hamming74_decode_codeword([int(b) for b in chunk])
        corrected += c
        out.extend(str(b) for b in data)
    return "".join(out), corrected


# ---------------------------------------------------------------------------
# Channel models
# ---------------------------------------------------------------------------

Bit = Union[int, str]  # 0/1 or 'E' erasure


def bsc_channel(bits: str, p: float, rng: random.Random) -> Tuple[str, int]:
    """Binary Symmetric Channel: flip each bit independently with probability p."""
    out: List[str] = []
    errors = 0
    for b in bits:
        if rng.random() < p:
            out.append("1" if b == "0" else "0")
            errors += 1
        else:
            out.append(b)
    return "".join(out), errors


def bec_channel(bits: str, epsilon: float, rng: random.Random) -> Tuple[str, int]:
    """Binary Erasure Channel: erase with probability epsilon."""
    out: List[str] = []
    erasures = 0
    for b in bits:
        if rng.random() < epsilon:
            out.append("E")
            erasures += 1
        else:
            out.append(b)
    return "".join(out), erasures


def hamming_decode_erasure_bits(bits: str) -> Tuple[str, int]:
    """Decode Hamming (7,4) with erasures by testing consistent completions."""
    out: List[str] = []
    resolved = 0
    for i in range(0, len(bits), 7):
        chunk = list(bits[i : i + 7].ljust(7, "0"))
        positions = [j for j, c in enumerate(chunk) if c == "E"]
        if not positions:
            data, _ = hamming74_decode_codeword([int(b) for b in chunk])
            out.extend(str(b) for b in data)
            continue

        best: Optional[List[int]] = None

        def search(pos: int, current: List[str]) -> None:
            nonlocal best
            if pos == len(positions):
                trial = [int(x) for x in current]
                s1 = _parity(trial[0], trial[2], trial[4], trial[6])
                s2 = _parity(trial[1], trial[2], trial[5], trial[6])
                s3 = _parity(trial[3], trial[4], trial[5], trial[6])
                if s1 == 0 and s2 == 0 and s3 == 0:
                    best = trial
                return
            idx = positions[pos]
            for bit in ("0", "1"):
                current[idx] = bit
                search(pos + 1, current)

        search(0, chunk)
        if best is None:
            best = [0 if c in ("0", "1") else 0 for c in chunk]
        data, _ = hamming74_decode_codeword(best)
        resolved += len(positions)
        out.extend(str(b) for b in data)
    return "".join(out), resolved


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

def psnr(original: np.ndarray, restored: np.ndarray) -> float:
    orig = original.astype(np.float64)
    rec = restored.astype(np.float64)
    mse = np.mean((orig - rec) ** 2)
    if mse <= 0:
        return float("inf")
    return 10.0 * math.log10(255.0 ** 2 / mse)


def pixel_accuracy(original: np.ndarray, restored: np.ndarray, tol: int = 0) -> float:
    diff = np.abs(original.astype(np.int16) - restored.astype(np.int16))
    match = np.sum(diff <= tol)
    return float(match) / diff.size


def bit_error_rate(sent: str, received: str) -> float:
    n = min(len(sent), len(received))
    if n == 0:
        return 0.0
    errs = sum(1 for i in range(n) if sent[i] != received[i] and received[i] != "E")
    return errs / n


@dataclass
class PipelineResult:
    name: str
    channel: str
    quality: int
    psnr_db: float
    pixel_acc: float
    ber: float
    source_bytes: int
    coded_bits: int
    channel_errors: int
    hamming_corrections: int
    encode_ms: float
    channel_ms: float
    decode_ms: float
    total_ms: float


def load_gray_image(path: Path, max_side: int = 256) -> np.ndarray:
    img = Image.open(path).convert("L")
    w, h = img.size
    if max(w, h) > max_side:
        scale = max_side / max(w, h)
        img = img.resize((int(w * scale), int(h * scale)), Image.Resampling.LANCZOS)
    return np.array(img, dtype=np.uint8)


def pack_bitstream(source_payload: bytes) -> bytes:
    bits = bytes_to_bits(source_payload)
    coded = hamming_encode_bits(bits)
    header = struct.pack(">II", len(source_payload), len(bits))
    return MAGIC + header + bits_to_bytes(coded)


def unpack_bitstream(
    data: bytes,
    channel_bits: str,
    channel: str,
) -> Tuple[bytes, int]:
    if not data.startswith(MAGIC):
        raise ValueError("无效的码流文件头")
    src_len, raw_len = struct.unpack_from(">II", data, len(MAGIC))
    if channel == "bec":
        decoded_bits, corrections = hamming_decode_erasure_bits(channel_bits)
    else:
        decoded_bits, corrections = hamming_decode_bits(channel_bits)
    decoded_bits = decoded_bits[:raw_len]
    payload = bits_to_bytes(decoded_bits)[:src_len]
    return payload, corrections


def run_single_image(
    path: Path,
    channel: str,
    p_bsc: float,
    eps_bec: float,
    quality: int,
    seed: int,
) -> PipelineResult:
    rng = random.Random(seed)
    gray = load_gray_image(path)
    name = path.stem

    t0 = time.perf_counter()
    source_payload = encode_grayscale(gray, quality)
    packet = pack_bitstream(source_payload)
    raw_bits = bytes_to_bits(source_payload)
    hamming_bits = hamming_encode_bits(raw_bits)
    coded_bits = repetition_encode(hamming_bits)
    t1 = time.perf_counter()

    if channel == "bsc":
        rx_rep, ch_err = bsc_channel(coded_bits, p_bsc, rng)
        rx_hamming = repetition_decode_hard(rx_rep)
    elif channel == "bec":
        rx_rep, ch_err = bec_channel(coded_bits, eps_bec, rng)
        rx_hamming = repetition_decode_erasure(rx_rep)
    else:
        rx_rep, ch_err = coded_bits, 0
        rx_hamming = repetition_decode_hard(rx_rep)
    t2 = time.perf_counter()

    payload, corrections = unpack_bitstream(packet, rx_hamming, channel)
    restored = decode_grayscale(payload)
    t3 = time.perf_counter()

    ber = bit_error_rate(hamming_bits, rx_hamming)
    return PipelineResult(
        name=name,
        channel=channel.upper(),
        quality=quality,
        psnr_db=psnr(gray, restored),
        pixel_acc=pixel_accuracy(gray, restored, tol=1),
        ber=ber,
        source_bytes=len(source_payload),
        coded_bits=len(coded_bits),
        channel_errors=ch_err,
        hamming_corrections=corrections,
        encode_ms=(t1 - t0) * 1000,
        channel_ms=(t2 - t1) * 1000,
        decode_ms=(t3 - t2) * 1000,
        total_ms=(t3 - t0) * 1000,
    )


def list_images() -> List[Path]:
    if not IMAGES_DIR.is_dir():
        return []
    return sorted(
        p for p in IMAGES_DIR.iterdir()
        if p.is_file() and p.suffix.lower() in {".png", ".jpg", ".jpeg", ".bmp", ".gif"}
    )


def require_images() -> List[Path]:
    """仅使用 images/ 目录中的图像；无图像时抛出 FileNotFoundError。"""
    imgs = list_images()
    if not imgs:
        raise FileNotFoundError(
            f"未在 {IMAGES_DIR} 中找到图像，请将 .png / .jpg / .bmp / .gif 放入该目录后再运行。"
        )
    return imgs


def save_restored(path: Path, quality: int, out_dir: Path) -> Path:
    gray = load_gray_image(path)
    payload = encode_grayscale(gray, quality)
    restored = decode_grayscale(payload)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{path.stem}_q{quality}_restored.png"
    Image.fromarray(restored).save(out_path)
    return out_path


def _display_width(text: str) -> int:
    """字符串在等宽终端中的显示宽度（中文计 2，ASCII 计 1）。"""
    return sum(2 if ord(c) > 127 else 1 for c in text)


def _pad_display(text: str, width: int, align: str = "left") -> str:
    gap = width - _display_width(text)
    if gap <= 0:
        return text
    if align == "right":
        return " " * gap + text
    if align == "center":
        left = gap // 2
        return " " * left + text + " " * (gap - left)
    return text + " " * gap


def channel_display(ch: str) -> str:
    """统一信道名为 IDEAL / BSC / BEC。"""
    key = ch.strip().lower()
    return {"ideal": "IDEAL", "bsc": "BSC", "bec": "BEC"}.get(key, ch.upper())


def print_result(r: PipelineResult) -> None:
    ch = channel_display(r.channel)
    label_w = 14
    rows = [
        ("图像", r.name),
        ("信道", ch),
        ("质量因子", str(r.quality)),
        ("PSNR", f"{r.psnr_db:.2f} dB"),
        ("像素准确度", f"{r.pixel_acc * 100:.2f}%（±1 灰度）"),
        ("BER（码流）", f"{r.ber * 100:.4f}%"),
        ("信源字节数", str(r.source_bytes)),
        ("编码比特数", str(r.coded_bits)),
        ("信道错误数", str(r.channel_errors)),
        ("Hamming 纠正", str(r.hamming_corrections)),
        (
            "耗时 (ms)",
            f"编码 {r.encode_ms:.1f}，信道 {r.channel_ms:.1f}，"
            f"译码 {r.decode_ms:.1f}，合计 {r.total_ms:.1f}",
        ),
    ]
    for label, value in rows:
        print(f"  {_pad_display(label, label_w)} : {value}")


def run_experiment(
    images: List[Path],
    channels: Sequence[str],
    qualities: Sequence[int],
    p_bsc: float,
    eps_bec: float,
    seed: int,
) -> List[PipelineResult]:
    results: List[PipelineResult] = []
    for path in images:
        for ch in channels:
            for q in qualities:
                r = run_single_image(path, ch, p_bsc, eps_bec, q, seed)
                results.append(r)
    return results


def print_summary(results: List[PipelineResult]) -> None:
    """汇总表：信道列与 IDEAL/BSC/BEC 等宽对齐，数值列右对齐。"""
    col_image = 12
    col_channel = 5  # IDEAL / BSC / BEC
    col_quality = 4
    col_psnr = 8
    col_acc = 8
    col_ber = 9
    col_ms = 8
    line_w = col_image + col_channel + col_quality + col_psnr + col_acc + col_ber + col_ms + 6

    def row(
        img: str,
        ch: str,
        q: str,
        psnr: str,
        acc: str,
        ber: str,
        ms: str,
    ) -> str:
        return (
            f"{_pad_display(img, col_image)} "
            f"{_pad_display(ch, col_channel, 'center')} "
            f"{_pad_display(q, col_quality, 'right')} "
            f"{_pad_display(psnr, col_psnr, 'right')} "
            f"{_pad_display(acc, col_acc, 'right')} "
            f"{_pad_display(ber, col_ber, 'right')} "
            f"{_pad_display(ms, col_ms, 'right')}"
        )

    print("\n" + "=" * line_w)
    print("汇总报告（Project 2 — 有损信源编码 + 信道编码）")
    print("=" * line_w)
    print(
        row("图像", "信道", "质量", "PSNR", "准确度", "BER", "ms")
    )
    print("-" * line_w)
    for r in results:
        print(
            row(
                r.name,
                channel_display(r.channel),
                str(r.quality),
                f"{r.psnr_db:.2f}",
                f"{r.pixel_acc * 100:.2f}",
                f"{r.ber * 100:.4f}",
                f"{r.total_ms:.1f}",
            )
        )
    avg_psnr = sum(r.psnr_db for r in results) / len(results)
    avg_acc = sum(r.pixel_acc for r in results) / len(results)
    avg_time = sum(r.total_ms for r in results) / len(results)
    print("-" * line_w)
    print(
        f"平均 PSNR: {avg_psnr:.2f} dB | 平均像素准确度: {avg_acc * 100:.2f}% | "
        f"平均耗时: {avg_time:.1f} ms"
    )
    print("\n复杂度（理论）:")
    print("  信源编/译码: O(N)，N 为像素数；每块 8×8 DCT 为常数工作量。")
    print(f"  重复码 ({REP_FACTOR},1) + Hamming (7,4): 对信源比特长度 B 为 O(B)。")
    print("=" * line_w)


def export_results_json(results: List[PipelineResult], path: Path) -> None:
    data = [r.__dict__ for r in results]
    path.write_text(json.dumps(data, indent=2), encoding="utf-8")


def _prompt_line(prompt: str) -> str:
    try:
        return input(prompt).strip()
    except (EOFError, KeyboardInterrupt):
        print()
        raise KeyboardInterrupt


def interactive_loop() -> int:
    print("=" * 60)
    print("  Project 2 — 图像有损信源编码 + IDEAL/BSC/BEC 信道")
    print(f"  工作目录: {BASE_DIR}")
    print("=" * 60)

    while True:
        print()
        print("  [1] 运行完整实验（IDEAL + BSC + BEC，打印报告）")
        print("  [2] 单张图像编解码（保存还原 PNG）")
        print("  [3] 查看最近一次实验 JSON（output/results.json）")
        print("  [0] 退出")
        try:
            action = _prompt_line("请选择功能 (0-3): ")
        except KeyboardInterrupt:
            print("\n再见。")
            return 0

        if action in ("0", "q", "Q", "exit", "退出"):
            print("再见。")
            return 0

        try:
            if action == "1":
                imgs = require_images()
                print(f"\n使用 {IMAGES_DIR} 下共 {len(imgs)} 张图像")
                results = run_experiment(
                    imgs,
                    channels=DEFAULT_CHANNELS,
                    qualities=(40, 75, 90),
                    p_bsc=0.02,
                    eps_bec=0.03,
                    seed=42,
                )
                for r in results:
                    print()
                    print_result(r)
                print_summary(results)
                export_results_json(results, OUTPUT_DIR / "results.json")
                print(f"\n已保存 {OUTPUT_DIR / 'results.json'}")
            elif action == "2":
                imgs = require_images()
                print("\n可选图像:")
                for i, p in enumerate(imgs, 1):
                    print(f"  [{i}] {p.name}")
                choice = _prompt_line("请输入图像编号: ")
                if not choice.isdigit() or not (1 <= int(choice) <= len(imgs)):
                    print("编号无效。")
                    continue
                q = int(_prompt_line("质量因子 1-100（回车默认 75）: ") or "75")
                out = save_restored(imgs[int(choice) - 1], q, OUTPUT_DIR)
                print(f"已保存 {out}")
            elif action == "3":
                p = OUTPUT_DIR / "results.json"
                if p.is_file():
                    print(p.read_text(encoding="utf-8")[:2000])
                else:
                    print("请先执行 [1] 生成 results.json")
            else:
                print("请输入 0、1、2 或 3。")
        except KeyboardInterrupt:
            print("\n已取消当前操作。")
        except FileNotFoundError as exc:
            print(f"提示: {exc}")
        except (OSError, ValueError) as exc:
            print(f"错误: {exc}")

        _prompt_line("\n按回车继续...")


def main(argv: Optional[Iterable[str]] = None) -> int:
    raw = list(argv) if argv is not None else sys.argv[1:]
    parser = argparse.ArgumentParser(description="Project 2 图像有损编码 + IDEAL/BSC/BEC 信道")
    parser.add_argument("--run", action="store_true", help="使用 images/ 中图像运行实验后退出")
    parser.add_argument("--p-bsc", type=float, default=0.02, help="BSC 交叉概率 p")
    parser.add_argument("--eps-bec", type=float, default=0.05, help="BEC 擦除概率 ε")
    parser.add_argument("--quality", type=int, default=75, help="有损编码质量因子 1-100")
    parser.add_argument("--seed", type=int, default=42, help="随机种子")
    args = parser.parse_args(raw)

    if args.run or raw:
        imgs = require_images()
        results = run_experiment(
            imgs,
            channels=DEFAULT_CHANNELS,
            qualities=(args.quality,),
            p_bsc=args.p_bsc,
            eps_bec=args.eps_bec,
            seed=args.seed,
        )
        for r in results:
            print()
            print_result(r)
        print_summary(results)
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        export_results_json(results, OUTPUT_DIR / "results.json")
        for path in imgs:
            save_restored(path, args.quality, OUTPUT_DIR)
        return 0

    return interactive_loop()


if __name__ == "__main__":
    raise SystemExit(main())
