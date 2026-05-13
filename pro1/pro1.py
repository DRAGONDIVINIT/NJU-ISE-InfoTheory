import heapq
from collections import Counter, defaultdict

class HuffmanCoding:
    def __init__(self):
        self.heap = []
        self.codes = {}
        self.reverse_codes = {}
    
    class Node:
        def __init__(self, char, freq, left=None, right=None):
            self.char = char
            self.freq = freq
            self.left = left
            self.right = right
        
        def __lt__(self, other):
            return self.freq < other.freq
    
    def build_huffman_tree(self, text):
        # 统计频率
        freq = Counter(text)
        
        # 构建最小堆
        for char, f in freq.items():
            heapq.heappush(self.heap, self.Node(char, f))
        
        # 合并节点
        while len(self.heap) > 1:
            left = heapq.heappop(self.heap)
            right = heapq.heappop(self.heap)
            merged = self.Node(None, left.freq + right.freq, left, right)
            heapq.heappush(self.heap, merged)
        
        return self.heap[0] if self.heap else None
    
    def generate_codes(self, node, code=""):
        if node is None:
            return
        if node.char is not None:
            self.codes[node.char] = code
            self.reverse_codes[code] = node.char
            return
        self.generate_codes(node.left, code + "0")
        self.generate_codes(node.right, code + "1")
    
    def compress(self, text):
        if not text:
            return "", {}
        root = self.build_huffman_tree(text)
        self.generate_codes(root)
        
        # 压缩
        compressed = "".join(self.codes[ch] for ch in text)
        return compressed, self.codes
    
    def decompress(self, compressed_bitstring):
        if not compressed_bitstring:
            return ""
        decoded = []
        current_code = ""
        for bit in compressed_bitstring:
            current_code += bit
            if current_code in self.reverse_codes:
                decoded.append(self.reverse_codes[current_code])
                current_code = ""
        return "".join(decoded)
    
    def calculate_efficiency(self, text, compressed_bitstring):
        original_bits = len(text) * 8  # ASCII 每字符8位
        compressed_bits = len(compressed_bitstring)
        compression_ratio = compressed_bits / original_bits
        avg_code_length = compressed_bits / len(text)
        
        # 计算熵 H(X)
        freq = Counter(text)
        total = len(text)
        entropy = -sum((count/total) * (count/total).bit_length()? 
                       # 更精确的熵计算
                       for count in freq.values())
        # 改用 math.log2
        import math
        entropy = -sum((count/total) * math.log2(count/total) for count in freq.values())
        
        coding_efficiency = entropy / avg_code_length if avg_code_length > 0 else 0
        
        return {
            "original_bits": original_bits,
            "compressed_bits": compressed_bits,
            "compression_ratio": compression_ratio,
            "avg_code_length": avg_code_length,
            "entropy": entropy,
            "coding_efficiency": coding_efficiency
        }


# 示例运行
if __name__ == "__main__":
    # 测试文本
    test_text = """
    this is a test text for lossless source coding using Huffman coding.
    The goal is to compress text and verify coding efficiency.
    """
    
    hc = HuffmanCoding()
    
    # 压缩
    compressed, codes = hc.compress(test_text)
    
    # 解压
    decompressed = hc.decompress(compressed)
    
    # 计算效率指标
    metrics = hc.calculate_efficiency(test_text, compressed)
    
    # 输出结果
    print("=" * 60)
    print("Project 1: Lossless Source Coding - Huffman Coding")
    print("=" * 60)
    print(f"原始文本长度: {len(test_text)} 字符")
    print(f"压缩后比特数: {metrics['compressed_bits']} bits")
    print(f"原始比特数 (ASCII): {metrics['original_bits']} bits")
    print(f"压缩率: {metrics['compression_ratio']:.2%}")
    print(f"平均码长: {metrics['avg_code_length']:.2f} bits/char")
    print(f"信源熵 H(X): {metrics['entropy']:.4f} bits/char")
    print(f"编码效率: {metrics['coding_efficiency']:.2%}")
    print(f"\n解码验证成功: {test_text == decompressed}")
    print("\n码表 (前10个):")
    for i, (ch, code) in enumerate(codes.items()):
        if i >= 10:
            break
        print(f"  '{ch}': {code}")