class Solution:
    def addBinary(self, a: str, b: str) -> str:
        # 利用 python 内置转换为int类型加和再表示为二进制
        return '{0:b}'.format( int(a, 2) + int(b, 2))