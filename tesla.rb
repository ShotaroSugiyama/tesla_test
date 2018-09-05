require "./const"
require "./polynomial_ring"
require "digest/sha2"

class ChaCha

  def initialize(seed)
    @x = [0x61707876, 0x3320646e, 0x79622d32, 0x6b206574]
    @x += seed
    @x += [0, 0, 0, 0]
    @y
  end

  def quarter_round(a, b, c, d)
    @x[a] = (@x[a] + @x[b]) % 2**32
    @x[d] ^= @x[a]
    @x[d] = ((@x[d] << 16) | (@x[d] >> (32-16)))  % 2**32
    @x[c] = (@x[c] + @x[d]) % 2**32
    @x[b] ^= @x[c]
    @x[b] = ((@x[b] << 12) | (@x[b] >> (32-12)))  % 2**32
    @x[a] = (@x[a] + @x[b]) % 2**32
    @x[d] ^= @x[a]
    @x[d] = ((@x[d] << 8) | (@x[d] >> (32-8)))  % 2**32
    @x[c] = (@x[c] + @x[d]) % 2**32
    @x[b] ^= @x[c]
    @x[b] = ((@x[b] << 7) | (@x[b] >> (32-7)))  % 2**32
  end

  def update
    @y = @x.map { |e| e }
    (0..9).each do |i|
      # Column Round
      quarter_round(0, 4, 8, 12)
      quarter_round(5, 9, 13, 1)
      quarter_round(10, 14, 2, 6)
      quarter_round(15, 3, 7, 11)
      # Diagonal Round
      quarter_round(0, 5, 10, 15)
      quarter_round(1, 6, 11, 12)
      quarter_round(2, 7, 8, 13)
      quarter_round(3, 4, 9, 14)
    end
    @x.map!.with_index { |e, i| (e + @y[i]) % 2**32 }
  end

end

module TESLAFunctions

  Generator = Random.new(0)

  def uniform_sampling
    (0..Const::N-1).map { |i| Generator.rand(Const::Q) }
  end

  def bounded_sampling
    (0..Const::N-1).map { |i|  Generator.rand(-Const::B..Const::B) % Const::Q}
  end

  def gaussian_sampling
    (0..Const::N-1).map do |i|
      x = Generator.rand
      y = Generator.rand
      z = Math.sqrt(-2.0*Math.log(x))*Math.sin(2*Math::PI*y)
      (z*Const::Sigma).to_i % Const::Q
    end
  end

  def d_rounding(v)
    v.map do |e|
      tmp = e % 2**(32-Const::D)
      tmp -= 2**Const::D if tmp > 2**(Const::D-1)
      (e - tmp) >> Const::D
    end
  end

  def hash512(v1, v2, message)
    d = Digest::SHA512.hexdigest(v1.pack("c*")+v2.pack("c*")+message)
    d = d.each_char.each_slice(8).map(&:join)
    d.map { |e| e.hex }
  end

  def hash_f(sig_c)
    shift = Array.new(Const::OmegaPr, 0)
    (0..15).each do |i|
      shift[3*i] = sig_c[i] & 0x3ff
      shift[3*i+1] = (sig_c[i] >> 10) & 0x3ff
      shift[3*i+2] = (sig_c[i] >> 20) & 0x3ff
    end
    (0..4).each do |i|
      shift[48] += (sig_c[i] >> (30 - 2*i)) & (0x3 << 2*i)
      shift[49] += (sig_c[i+5] >> (30 - 2*i)) & (0x3 << 2*i)
      shift[50] += (sig_c[i+10] >> (30 - 2*i)) & (0x3 << 2*i)
    end

    count = 0
    f_digest = Array.new(Const::N, 0)
    (0..Const::OmegaPr-1).each do |i|
      if f_digest[shift[i]] == 0
        f_digest[shift[i]] = 1
        count += 1
      end
      break if count == Const::Omega
    end
    f_digest
  end

  module_function :uniform_sampling, :bounded_sampling, \
  :gaussian_sampling, :d_rounding, :hash512, :hash_f

end

class TESLA256

  class << self
    include TESLAFunctions

    def key_gen
      key = {}
      key[:s] = gaussian_sampling
      key[:e1] = gaussian_sampling
      key[:e2] = gaussian_sampling

      key[:seed] = (0..7).map { |e| Generator.rand(2**32) }
      a = []
      chacha = ChaCha.new(key[:seed])
      # Concat
      (0..2*Const::N/16-1).each do |i|
        a += (chacha.update).map { |e| e % Const::Q }
      end

      p_ring = PolynomialRing.new(q: Const::Q, n: Const::N)
      s = p_ring.dwt(key[:s])
      a1 = p_ring.dwt(a[0..Const::N-1])
      a2 = p_ring.dwt(a[Const::N..-1])
      key[:t1] = p_ring.add(p_ring.idwt(p_ring.inner_prod(s, a1)), key[:e1])
      key[:t2] = p_ring.add(p_ring.idwt(p_ring.inner_prod(s, a2)), key[:e2])
      key
    end

    def sign(key, message)
      signature = {c:[], z:[]}

      a = []
      chacha = ChaCha.new(key[:seed])
      # Concat
      (0..2*Const::N/16-1).each do |i|
        a += (chacha.update).map { |e| e % Const::Q }
      end

      p_ring = PolynomialRing.new(q: Const::Q, n: Const::N)
      a1 = p_ring.dwt(a[0..Const::N-1])
      a2 = p_ring.dwt(a[Const::N..-1])

      s = p_ring.dwt(key[:s])
      e1 = p_ring.dwt(key[:e1])
      e2 = p_ring.dwt(key[:e2])

      loop do
        r = p_ring.dwt(bounded_sampling)

        v1 = p_ring.inner_prod(a1, r)
        v2 = p_ring.inner_prod(a2, r)
        v1_d = d_rounding(p_ring.idwt(v1))
        v2_d = d_rounding(p_ring.idwt(v2))

        sig_c = hash512(v1_d, v2_d, message)
        c = hash_f(sig_c)
        next if c.count(1) != Const::Omega

        c = p_ring.dwt(c)
        sig_z = p_ring.add(r, p_ring.inner_prod(s, c))
        sig_z = p_ring.idwt(sig_z)
        next if sig_z.max_by { |n| [n, Const::Q - n].min } > (Const::B - Const::U)

        w1 = p_ring.sub(v1, p_ring.inner_prod(e1, c))
        next if d_rounding(p_ring.idwt(w1)) != v1_d

        w2 = p_ring.sub(v2, p_ring.inner_prod(e2, c))
        next if d_rounding(p_ring.idwt(w2)) != v2_d

        signature[:c] = sig_c
        signature[:z] = sig_z

        break
      end
      signature
    end

    def verify(key, signature, message)
      output = 1

      a = []
      chacha = ChaCha.new(key[:seed])
      # Concat
      (0..2*Const::N/16-1).each do |i|
        a += (chacha.update).map { |e| e % Const::Q }
      end

      p_ring = PolynomialRing.new(q: Const::Q, n: Const::N)
      a1 = p_ring.dwt(a[0..Const::N-1])
      a2 = p_ring.dwt(a[Const::N..-1])

      c = hash_f(signature[:c])

      c = p_ring.dwt(c)
      z = p_ring.dwt(signature[:z])
      t1 = p_ring.dwt(key[:t1])
      t2 = p_ring.dwt(key[:t2])

      w1 = p_ring.sub(p_ring.inner_prod(a1, z), p_ring.inner_prod(t1, c))
      w2 = p_ring.sub(p_ring.inner_prod(a2, z), p_ring.inner_prod(t2, c))
      w1 = p_ring.idwt(w1)
      w2 = p_ring.idwt(w2)

      c_pr = hash512(d_rounding(w1), d_rounding(w2), message)

      output = 0 if signature[:c] != c_pr
      output = 0 if signature[:z].max_by { |n| [n, Const::Q - n].min } > (Const::B - Const::U)
      output
    end

  end

end
