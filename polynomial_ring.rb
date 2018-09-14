require "prime"

module Fp

  # モジュラべき乗
  def pow(base, exp, mod)
    ans = 1
    while exp > 0
      ans = (ans * base) % mod if (exp & 1) == 1
      exp >>= 1
      base = base**2 % mod
    end
    ans
  end

  # 剰余環Z/(2**28-2**16+1)Zの原始根を求める
  def primitive_root
    q = 2**28-2**16+1
    pr = 1
    while true
      pr += 2
      next if Prime.prime?(pr) == false
      # 相互法則で得られる値は常に1 = (-1)^(2**27-2**15)*(-1)^((pr-1)/2)
      # (2**28-2**16+1) mod pr が平方非剰余であればよい
      next if ((0..pr-1).map { |i| 1 if (i**2 % pr) == (q % pr)}).count(1) != 0
      # 2**28 - 2**16 = 2**16*3**2*5*7*13より、9,5,7,13乗根でないことを確かめる
      next if pow(pr, (q-1)/9, q) == 1
      next if pow(pr, (q-1)/5, q) == 1
      next if pow(pr, (q-1)/7, q) == 1
      next if pow(pr, (q-1)/13, q) == 1
      break
    end
    # まあ答えは23なのですが
    pr
  end

  def mod_mult_tesla256(a, b)
    q = 2**28 - 2**16 + 1
    a0 = a & 0x3fff
    a1 = a >> 14
    b0 = b & 0x3fff;
    b1 = b >> 14;
    a0b0 = a0*b0;
    a1b1 = a1*b1;
    middle = (a0 + a1)*(b0 + b1) % q
    middle = (middle - a0b0 - a1b1) % q
    a0 = (a0b0 + ((middle & 0x3fff) << 14) + (middle >> 14)*0xffff) % q
    a0 = (a0 - a1b1) % q
    a1 = (((a1b1 >> 12)*0xffff % q) + ((a1b1 & 0xfff) << 16) % q) % q
    (a0 + a1) % q
  end

  module_function :pow, :primitive_root, :mod_mult_tesla256

end

# 多項式環
class PolynomialRing

  include Fp

  def initialize(q:, n:)
    @q = q
    @n = n
    # 正しくは @prim_root = primitive_root()です
    # 実際には11を使っていました。11は11^((q-1)/5)=1となってしまうので原始根ではありません
    # しかし、数論変換が2バタフライなので5で巡回していても理論上問題ないです。
    # 数論変換なしで多項式乗算した結果と照らし合わせて合っていたので問題ないはずです。
    # 本番で動いていたほうを残しておきます。
    @prim_root = 11
    @n_inv = pow(@n, @q-2, @q)
    @stage = Math.log2(@n).to_i
    @omega = [1, pow(@prim_root, (@q-1)/2/@n, @q)]
    (2..2*@n).each do |i|
      @omega += [@omega[i-1]*@omega[1] % @q]
    end
  end

  # 数論変換&負巡回畳み込み
  def dwt!(poly)
    # 負巡回畳み込みなので先にomega^{i/2}をかけておく
    poly.map!.with_index { |e, i| e*@omega[i] % @q}
    # 以下FFTとおなじ
    base = @n
    (0..@stage-1).each do |stage|
      (0..@n/base-1).each do |i|
        offset = i*base;
        (0..base/2-1).each do |j|
          radix2_a = poly[offset+j]
          radix2_b = poly[offset+j+base/2]
          poly[offset+j] = (radix2_a + radix2_b) % @q
          poly[offset+j+base/2] = (radix2_a - radix2_b) * @omega[j*2*@n/base] % @q
        end
      end
      base /= 2
    end
  end

  # 数論逆変換
  def idwt!(poly)
    base = 2
    (0..@stage-1).each do |stage|
      (0..@n/base-1).each do |i|
        offset = i*base;
        (0..base/2-1).each do |j|
          radix2_a = poly[offset+j]
          radix2_b = poly[offset+j+base/2] * @omega[2*@n-j*2*@n/base] % @q
          poly[offset+j] = (radix2_a + radix2_b) % @q
          poly[offset+j+base/2] = (radix2_a - radix2_b) % @q
        end
      end
      base *= 2
    end
    # 負巡回畳み込み
    poly.map!.with_index { |e, i| e*@n_inv*@omega[2*@n-i] % @q}
  end

  # 変数を破壊しないバージョン
  def dwt(poly)
    mem = poly.map.with_index { |e, i| e*@omega[i] % @q}
    base = @n
    (0..@stage-1).each do |stage|
      (0..@n/base-1).each do |i|
        offset = i*base;
        (0..base/2-1).each do |j|
          radix2_a = mem[offset+j]
          radix2_b = mem[offset+j+base/2]
          mem[offset+j] = (radix2_a + radix2_b) % @q
          mem[offset+j+base/2] = (radix2_a - radix2_b) * @omega[j*2*@n/base] % @q
        end
      end
      base /= 2
    end
    mem
  end

  def idwt(poly)
    mem = poly.map { |e| e }
    base = 2
    (0..@stage-1).each do |stage|
      (0..@n/base-1).each do |i|
        offset = i*base;
        (0..base/2-1).each do |j|
          radix2_a = mem[offset+j]
          radix2_b = mem[offset+j+base/2] * @omega[2*@n-j*2*@n/base] % @q
          mem[offset+j] = (radix2_a + radix2_b) % @q
          mem[offset+j+base/2] = (radix2_a - radix2_b) % @q
        end
      end
      base *= 2
    end
    mem.map.with_index { |e, i| e*@n_inv*@omega[2*@n-i] % @q}
  end

  # 多項式乗算の検算用
  def conv(a, b)
    c = Array.new(@n, 0)
    (0..@n-1).to_a.repeated_permutation(2) do |i, j|
      if i + j >= @n
        c[i+j-@n] = (c[i+j-@n] - a[i]*b[j]) % @q
      else
        c[i+j] = (c[i+j] + a[i]*b[j]) % @q
      end
    end
    c
  end

  def add(a, b)
    a.map.with_index { |e, i| (e + b[i]) % @q}
  end

  def sub(a, b)
    a.map.with_index { |e, i| (e - b[i]) % @q}
  end

  # 要素積
  def inner_prod(a, b)
    a.map.with_index { |e, i| (e * b[i]) % @q}
  end

end
